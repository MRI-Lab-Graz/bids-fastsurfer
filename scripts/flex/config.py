"""
Shared study-definition config loader for the flex pipeline.

Loads a study.json, validates it against configs/flex/study.schema.json
(raising with the JSON-pointer path of the offending field), expands
built-in measure profiles, and resolves the handful of relative paths
(region_groups, moderators.file when given relative to the repo root)
against the config file's own directory or the repo root.

Shared by ingest.py, extract_freesurfer.py (--config mode) and
run_study.py so the three stages agree on one definition of a valid
config -- there is no separate validation logic in the R engine, which
treats a config that reached it as already validated (see
scripts/flex/R/core_config.R).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = REPO_ROOT / "configs" / "flex" / "study.schema.json"

# Built-in measure profiles: {transform, etiv_covariate, hemisphere, aggregate}.
# See docs/FLEX_PIPELINE.md for the rationale behind each one -- carried
# over verbatim from 01_primary_lmm.R (volume) and 29_hippo_thickness_lmm.R
# / 25_brainstem_lmm.R (thickness / volume_midline).
#
# `aggregate` is the FUN used when a measure's ROI has multiple raw rows
# per subject/session/hemisphere (e.g. hippo subfields' head/body split):
# volumes are physically additive (sum), thickness is not (mean) -- this
# is a real, consistent convention across every original script (01/08
# aggregate volume with FUN=sum; 27/29/33 aggregate thickness with
# FUN=mean), not a stylistic choice.
BUILTIN_PROFILES: Dict[str, Dict[str, Any]] = {
    "volume": {"transform": "log", "etiv_covariate": True, "hemisphere": "pooled", "aggregate": "sum"},
    "thickness": {"transform": "identity", "etiv_covariate": False, "hemisphere": "pooled", "aggregate": "mean"},
    "volume_midline": {"transform": "log", "etiv_covariate": True, "hemisphere": "none", "aggregate": "sum"},
}


class ConfigError(ValueError):
    pass


def _validate_schema(cfg: Dict[str, Any]) -> None:
    try:
        import jsonschema
    except ImportError as e:
        raise ConfigError(
            "the 'jsonschema' package is required to validate flex study configs "
            "(pip install jsonschema)"
        ) from e

    with open(SCHEMA_PATH) as f:
        schema = json.load(f)

    validator_cls = jsonschema.validators.validator_for(schema)
    validator_cls.check_schema(schema)
    validator = validator_cls(schema)

    errors = sorted(validator.iter_errors(cfg), key=lambda e: list(e.absolute_path))
    if errors:
        lines = []
        for e in errors:
            pointer = "/" + "/".join(str(p) for p in e.absolute_path)
            lines.append(f"  {pointer or '/'}: {e.message}")
        raise ConfigError(
            f"study config failed schema validation ({len(errors)} error(s)):\n" + "\n".join(lines)
        )


def _validate_semantics(cfg: Dict[str, Any]) -> None:
    """Cross-field checks the JSON Schema itself can't express."""
    dataset_names = set(cfg["datasets"].keys())
    for i, m in enumerate(cfg["measures"]):
        if m["dataset"] not in dataset_names:
            raise ConfigError(f"/measures/{i}/dataset: unknown dataset '{m['dataset']}' "
                               f"(known: {sorted(dataset_names)})")
        profile_names = set(BUILTIN_PROFILES) | set(cfg.get("profiles", {}))
        if m["profile"] not in profile_names:
            raise ConfigError(f"/measures/{i}/profile: unknown profile '{m['profile']}' "
                               f"(known: {sorted(profile_names)})")

    measure_names = {m["name"] for m in cfg["measures"]}
    for i, a in enumerate(cfg.get("analyses", [])):
        for mname in a["measures"]:
            if mname not in measure_names:
                raise ConfigError(f"/analyses/{i}/measures: unknown measure '{mname}' "
                                   f"(known: {sorted(measure_names)})")

    design = cfg["design"]
    contrasts = design.get("contrasts", {})
    if "control" in contrasts and "intervention" in contrasts:
        overlap = set(contrasts["intervention"]) & set(contrasts["control"])
        if overlap:
            raise ConfigError(f"/design/contrasts: group level(s) {sorted(overlap)} appear in "
                               "both 'intervention' and 'control'")

    factorial = design.get("factorial", {})
    for axis, levels in factorial.items():
        seen: Dict[str, str] = {}
        for level_name, group_values in levels.items():
            for g in group_values:
                if g in seen and seen[g] != level_name:
                    raise ConfigError(f"/design/factorial/{axis}: group '{g}' assigned to both "
                                       f"'{seen[g]}' and '{level_name}'")
                seen[g] = level_name


def load_config(path: str | Path, validate: bool = True) -> Dict[str, Any]:
    path = Path(path)
    if not path.is_file():
        raise ConfigError(f"study config not found: {path}")
    with open(path) as f:
        cfg = json.load(f)

    if validate:
        _validate_schema(cfg)
        _validate_semantics(cfg)

    cfg["_config_path"] = str(path.resolve())
    cfg["_config_dir"] = str(path.resolve().parent)
    return cfg


def resolve_profile(cfg: Dict[str, Any], profile_name: str) -> Dict[str, Any]:
    custom = cfg.get("profiles", {})
    if profile_name in custom:
        return custom[profile_name]
    return BUILTIN_PROFILES[profile_name]


def resolve_path(cfg: Dict[str, Any], p: str) -> Path:
    """Relative paths in a config are resolved against the repo root (configs/
    and results paths are conventionally written relative to it, matching
    every other config in this repo -- see configs/meta_pipeline.TEMPLATE.json)."""
    pp = Path(p)
    return pp if pp.is_absolute() else (REPO_ROOT / pp)


def dataset_path(cfg: Dict[str, Any], dataset_name: str) -> Path:
    return resolve_path(cfg, cfg["datasets"][dataset_name]["path"])


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Validate a flex study config")
    p.add_argument("config")
    args = p.parse_args()
    try:
        cfg = load_config(args.config)
    except ConfigError as e:
        print(f"INVALID: {e}")
        raise SystemExit(1)
    print(f"OK: {args.config} ({len(cfg['measures'])} measures, {len(cfg.get('analyses', []))} analyses)")
