#!/bin/bash
# Batch runner for hipsta (hippocampal shape/thickness analysis,
# https://github.com/Deep-MI/Hipsta) over a directory of longitudinal
# FreeSurfer subfield segmentations. Companion to run_subfield_segmentation.sh
# (which produces the volume-based hippoSfVolumes.long.txt files consumed by
# extract_hippo_subfields.py); this script's output feeds
# extract_hipsta_thickness.py instead -- both now live in the
# MRI-Lab-Graz/flex-analysis repo, since this repo only produces
# preprocessing output, not analysis of it.
#
# Usage: ./run_hipsta_segmentation.sh [options]
#
# Flag-based rather than positional (unlike run_subfield_segmentation.sh)
# because there are several required paths (FreeSurfer subfield derivatives,
# hipsta container, FreeSurfer license, output dir) with no natural single
# positional argument among them; defaults below match study 129's current
# layout, override any of them for a different study.

set -uo pipefail

usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  --fs-deriv PATH       FreeSurfer derivatives dir with sub-*.long.* session
                         dirs under mri/<hemi>.hippoAmygLabels.long.mgz
                         [default: $FS_DERIV]
  --output-dir PATH      hipsta output directory [default: $OUTPUT_BASE]
  --container PATH       hipsta apptainer/singularity .sif image [default: $CONTAINER]
  --fs-license PATH      FreeSurfer license.txt (host path) [default: $FS_LICENSE_HOST]
  --seg-suffix NAME       Segmentation filename suffix, appended to "<hemi>."
                         [default: $SEG_SUFFIX]
                         IMPORTANT: use the NATIVE-resolution segmentation
                         (hippoAmygLabels.long.mgz, ~0.33mm), NOT the
                         ".FSvoxelSpace." resample (~0.9mm in a typical FS
                         conformed space) -- hipsta's default mask filters are
                         parameterised in voxels (e.g. a 1-voxel-sigma
                         gaussian), tuned for ~0.33mm. At ~0.9mm the same
                         filter erodes through the hippocampal ribbon (only
                         ~2.3 voxels thick there) and most hemispheres fail
                         checkSurface with "contains holes" -- a resolution
                         artifact, not a real segmentation defect. See
                         Deep-MI/hipsta#25.
  --log-dir PATH          Per-run stdout logs and summary.tsv [default: $LOG_DIR]
  --jobs N                Parallel sessions (both hemispheres of a session run
                         sequentially within one job -- see note below)
                         [default: $JOBS]
  --qc                   Enable hipsta's QC plots (off by default: kaleido in
                         the published 0.10.1 image needs Google Chrome, which
                         isn't installed there, and every run dies at
                         qcPlots(); see Deep-MI/hipsta#24. Only pass --qc with
                         a container image where this is fixed.)
  --limit N               Only process the first N sessions (sorted), for a
                         pilot run
  -h, --help              Show this help and exit
EOF
}

HOST_BASE="/data/local/129_PK01"
FS_DERIV="$HOST_BASE/derivatives/freesurfer_hpc"
OUTPUT_BASE="$HOST_BASE/derivatives/hipsta"
CONTAINER="/data/local/software/apptainer_images/all/hipsta_0.10.1.sif"
FS_LICENSE_HOST="/data/local/freesurfer/license.txt"
SEG_SUFFIX="hippoAmygLabels.long.mgz"
LOG_DIR="/data/local/tmp_big/hipsta_batch_logs"
JOBS=8
QC_ARGS="--no-qc"
LIMIT=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fs-deriv) FS_DERIV="$2"; shift 2 ;;
    --output-dir) OUTPUT_BASE="$2"; shift 2 ;;
    --container) CONTAINER="$2"; shift 2 ;;
    --fs-license) FS_LICENSE_HOST="$2"; shift 2 ;;
    --seg-suffix) SEG_SUFFIX="$2"; shift 2 ;;
    --log-dir) LOG_DIR="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --qc) QC_ARGS=""; shift ;;
    --limit) LIMIT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

mkdir -p "$LOG_DIR"
SUMMARY="$LOG_DIR/summary.tsv"

# ------------------------------------------------------------------------------
# preflight

for f in "$CONTAINER" "$FS_LICENSE_HOST"; do
    [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done
[ -d "$FS_DERIV" ] || { echo "ERROR: missing $FS_DERIV" >&2; exit 1; }

sessions=$(cd "$FS_DERIV" && ls -d *.long.* 2>/dev/null | sort)
if [ "$LIMIT" -gt 0 ]; then
    sessions=$(echo "$sessions" | head -n "$LIMIT")
fi
n_sessions=$(echo "$sessions" | grep -c .)
n_hemis=$((n_sessions * 2))
echo "Sessions to process: $n_sessions  ($n_hemis hemispheres)"

# disk check: hipsta output at native (~0.33mm) resolution runs roughly 40MB
# per hemisphere
need_mb=$((n_hemis * 40))
avail_mb=$(df -Pm "$(dirname "$OUTPUT_BASE")" | awk 'NR==2 {print $4}')
echo "Disk: need ~${need_mb} MB, have ${avail_mb} MB free"
if [ "$avail_mb" -lt "$need_mb" ]; then
    echo "ERROR: not enough free space on $(dirname "$OUTPUT_BASE")." >&2
    echo "       Free up space or lower --limit, then re-run." >&2
    exit 1
fi

# The native-resolution segmentations are typically git-annex symlinks in a
# datalad clone and not fetched by default -- retrieve them up front, in one
# pass, rather than racing parallel git-annex calls against the same dataset.
if command -v datalad >/dev/null 2>&1 && (cd "$FS_DERIV" && git rev-parse --git-dir >/dev/null 2>&1); then
    echo "Fetching $SEG_SUFFIX from the annex (if needed) ..."
    get_paths=()
    for s in $sessions; do
        for hemi in lh rh; do
            f="$s/mri/${hemi}.${SEG_SUFFIX}"
            [ -f "$FS_DERIV/$f" ] || get_paths+=("$f")
        done
    done
    if [ "${#get_paths[@]}" -gt 0 ]; then
        echo "  ${#get_paths[@]} files to fetch"
        (cd "$FS_DERIV" && datalad get -J 4 "${get_paths[@]}") || \
            echo "WARNING: some inputs could not be fetched; those hemispheres will be skipped" >&2
    else
        echo "  all inputs already present"
    fi
fi

# Results produced from a different input segmentation (e.g. an earlier run
# against the wrong resolution -- see --seg-suffix above) must not be mixed
# into the same tree. Archive rather than delete, so nothing is lost if this
# turns out to be the wrong call.
MARKER="$OUTPUT_BASE/.hipsta-batch-input"
if [ -d "$OUTPUT_BASE" ] && [ -n "$(ls -A "$OUTPUT_BASE" 2>/dev/null)" ]; then
    if [ ! -f "$MARKER" ] || [ "$(cat "$MARKER")" != "$SEG_SUFFIX" ]; then
        archive="${OUTPUT_BASE}.superseded-$(date +%Y%m%dT%H%M%S)"
        echo "Existing results were produced from a different input segmentation."
        echo "Archiving $OUTPUT_BASE -> $archive"
        mv "$OUTPUT_BASE" "$archive" || exit 1
    else
        echo "Resuming into $OUTPUT_BASE (hemispheres already marked OK are skipped)"
    fi
fi
mkdir -p "$OUTPUT_BASE"
echo "$SEG_SUFFIX" > "$MARKER"

[ -f "$SUMMARY" ] || printf "session\themi\tstatus\tmethod\n" > "$SUMMARY"

# ------------------------------------------------------------------------------

write_summary() {
    flock "$SUMMARY" printf "%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$4" >> "$SUMMARY"
}

run_one() {
    local session="$1" hemi="$2"
    local infile="$FS_DERIV/$session/mri/${hemi}.${SEG_SUFFIX}"
    local outdir="$OUTPUT_BASE/$session"
    local stamp="$outdir/hipsta-status_${hemi}.txt"

    mkdir -p "$outdir"

    # resumable: don't redo a hemisphere that already completed
    if [ -f "$stamp" ] && grep -q "^status: OK" "$stamp"; then
        echo "  skip (already OK): $session $hemi"
        return
    fi

    if [ ! -f "$infile" ]; then
        printf "status: FAIL\nmethod: input-missing\ninput: %s\ndate: %s\n" \
            "$infile" "$(date -Is)" > "$stamp"
        write_summary "$session" "$hemi" "FAIL" "input-missing"
        return
    fi

    # Tier 1: defaults. Tier 2: --long-filter. Tier 3: --long-filter plus a
    # wider gaussian filter. With a native-resolution segmentation tier 1
    # should carry the large majority; the later tiers are a fallback. Tier 3
    # trades "holes" failures for boundary-loop failures, because widening the
    # gaussian pushes the filtered mask the mesh is built from away from the
    # unfiltered labels the head/tail cut planes are derived from.
    local tiers=(
        "default:"
        "long-filter:--long-filter"
        "long-filter+gauss:--long-filter --gauss-filter-size 2 40"
    )

    for tier in "${tiers[@]}"; do
        local method="${tier%%:*}"
        local extra_args="${tier#*:}"
        local runlog="$LOG_DIR/${session}_${hemi}_${method}.log"

        apptainer run \
            -B "$HOST_BASE:$HOST_BASE" \
            -B "$FS_DERIV:$FS_DERIV" \
            -B "$OUTPUT_BASE:$OUTPUT_BASE" \
            -B "$FS_LICENSE_HOST:/opt/freesurfer/.license" \
            "$CONTAINER" \
            --filename "$infile" \
            --hemi "$hemi" --lut freesurfer \
            --outputdir "$outdir" $QC_ARGS $extra_args \
            > "$runlog" 2>&1

        # hipsta always writes <outdir>/logfile.txt with mode="w" and no
        # hemisphere in the name, so keep a per-hemisphere copy before the next
        # run (next tier, or the other hemisphere) clobbers it.
        [ -f "$outdir/logfile.txt" ] && mv "$outdir/logfile.txt" "$outdir/logfile_${hemi}.txt"

        if grep -q "Hipsta finished without errors" "$runlog"; then
            printf "status: OK\nmethod: %s\ninput: %s\ndate: %s\n" \
                "$method" "$infile" "$(date -Is)" > "$stamp"
            write_summary "$session" "$hemi" "OK" "$method"
            return
        fi
    done

    # A failed run still leaves surf/tetra/cut files behind, so record the
    # verdict in the output directory itself rather than leaving it ambiguous.
    printf "status: FAIL\nmethod: all-tiers-exhausted\ninput: %s\ndate: %s\n" \
        "$infile" "$(date -Is)" > "$stamp"
    write_summary "$session" "$hemi" "FAIL" "all-tiers-exhausted"
}

# Both hemispheres share one output directory, so they must not run at the same
# time: they would race on logfile.txt and on the intermediate subdirectories.
# Parallelise across sessions, and do lh then rh within a session.
run_session() {
    local session="$1"
    for hemi in lh rh; do
        run_one "$session" "$hemi"
    done
}

export -f run_one run_session write_summary
export FS_DERIV OUTPUT_BASE CONTAINER FS_LICENSE_HOST HOST_BASE LOG_DIR SUMMARY \
       SEG_SUFFIX QC_ARGS

echo "$sessions" | xargs -P "$JOBS" -n 1 bash -c 'run_session "$0"'

echo
echo "Done. Summary: $SUMMARY"
awk -F'\t' 'NR>1 {c[$3]++} END {for (k in c) printf "  %-6s %d\n", k, c[k]}' "$SUMMARY"
