#!/usr/bin/env Rscript
# Quick debug script to check ID matching between qdec and aseg files

# Read the files
qdec_file <- "/data/local/129_PK01/derivatives/fastsurfer/stats/qdec.table.dat"
aseg_file <- "/data/local/129_PK01/derivatives/fastsurfer/stasts.aseg.long.table"

cat("Reading qdec file...\n")
qdec <- read.delim(qdec_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

cat("Reading aseg file...\n")
aseg <- read.table(aseg_file, header=TRUE, stringsAsFactors=FALSE)

cat("QDEC file structure:\n")
cat("Columns:", paste(names(qdec), collapse=", "), "\n")
cat("First 3 rows:\n")
print(head(qdec, 3))

cat("\nASEG file structure:\n")
cat("Columns:", paste(names(aseg)[1:5], collapse=", "), "... (truncated)\n")
cat("First column name:", names(aseg)[1], "\n")
cat("First 3 values in first column:\n")
print(head(aseg[,1], 3))

# Extract IDs from aseg (simulate the script logic)
id_col <- names(aseg)[1]  # First column should be the ID column
cat("\nExtracting IDs from aseg using column:", id_col, "\n")

aseg_fsid <- sub("\\.long\\..*", "", aseg[[id_col]])
aseg_fsid_base <- sub(".*\\.long\\.", "", aseg[[id_col]])

cat("Original aseg IDs:\n")
print(head(aseg[[id_col]], 3))
cat("Extracted fsid:\n")
print(head(aseg_fsid, 3))
cat("Extracted fsid.base:\n")
print(head(aseg_fsid_base, 3))

# Check what we have in qdec
cat("\nQDEC IDs:\n")
cat("qdec$fsid:\n")
print(head(qdec$fsid, 3))
if ("fsid-base" %in% names(qdec)) {
  cat("qdec$fsid-base:\n")
  print(head(qdec$`fsid-base`, 3))
}

# Check overlaps
cat("\nID overlap analysis:\n")
cat("Common fsid values:", length(intersect(qdec$fsid, aseg_fsid)), "out of", length(unique(qdec$fsid)), "qdec entries\n")
if ("fsid-base" %in% names(qdec)) {
  cat("Common fsid-base values:", length(intersect(qdec$`fsid-base`, aseg_fsid_base)), "\n")
}

# Show some examples
cat("\nFirst 5 unique fsid values from each:\n")
cat("QDEC fsid:\n")
print(head(unique(qdec$fsid), 5))
cat("ASEG extracted fsid:\n")
print(head(unique(aseg_fsid), 5))

if ("fsid-base" %in% names(qdec)) {
  cat("\nFirst 5 unique fsid-base values from each:\n")
  cat("QDEC fsid-base:\n")
  print(head(unique(qdec$`fsid-base`), 5))
  cat("ASEG extracted fsid.base:\n")
  print(head(unique(aseg_fsid_base), 5))
}