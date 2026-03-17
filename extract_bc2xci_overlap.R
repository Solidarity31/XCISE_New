#!/usr/bin/env Rscript
# extract_bc2xci_overlap.R
#
# Extract cells (barcodes) present in both bc2xci.txt files and show
# which category each cell was assigned to in each pipeline.
#
# Usage:
#   Rscript extract_bc2xci_overlap.R file1_bc2xci.txt file2_bc2xci.txt [output.tsv]
#
# Output:
#   - Console summary of overlap counts per category pair
#   - TSV file (optional 3rd argument, default: bc2xci_overlap.tsv)
#     Columns: Barcode | Assign_pipe1 | Assign_pipe2

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript extract_bc2xci_overlap.R file1_bc2xci.txt file2_bc2xci.txt [output.tsv]")
}

file1   <- args[1]
file2   <- args[2]
outfile <- if (length(args) >= 3) args[3] else "bc2xci_overlap.tsv"

## --- Read bc2xci files -------------------------------------------------------
read_bc2xci <- function(path) {
  df <- read.table(
    path,
    header          = FALSE,
    sep             = "\t",
    quote           = "",
    comment.char    = "",
    stringsAsFactors = FALSE
  )
  if (ncol(df) < 4) {
    stop(sprintf("File '%s' has fewer than 4 columns (expected: Barcode, Val1, Val2, Assign).", path))
  }
  df <- df[, c(1, 4)]
  colnames(df) <- c("Barcode", "Assign")
  df$Barcode <- trimws(df$Barcode)
  # Keep first occurrence of each barcode
  df[!duplicated(df$Barcode), ]
}

cat("Reading:", file1, "\n")
df1 <- read_bc2xci(file1)
cat("Reading:", file2, "\n")
df2 <- read_bc2xci(file2)

pipe1_name <- basename(file1)
pipe2_name <- basename(file2)

## --- Extract overlap (inner join on Barcode) ---------------------------------
overlap <- merge(df1, df2, by = "Barcode", suffixes = c("_pipe1", "_pipe2"))
colnames(overlap) <- c("Barcode", "Assign_pipe1", "Assign_pipe2")

n1       <- nrow(df1)
n2       <- nrow(df2)
n_over   <- nrow(overlap)
only_p1  <- n1 - n_over
only_p2  <- n2 - n_over

## --- Console summary ---------------------------------------------------------
cat("\n=== Overlap summary ===\n")
cat(sprintf("  %-40s %d barcodes\n", pipe1_name, n1))
cat(sprintf("  %-40s %d barcodes\n", pipe2_name, n2))
cat(sprintf("  Overlap (in both):%*d barcodes\n", 40 - nchar("Overlap (in both):") + 2, n_over))
cat(sprintf("  Only in %s: %d\n", pipe1_name, only_p1))
cat(sprintf("  Only in %s: %d\n", pipe2_name, only_p2))

## --- Category breakdown for overlapping cells --------------------------------
tab <- table(
  Pipe1 = overlap$Assign_pipe1,
  Pipe2 = overlap$Assign_pipe2
)

agree    <- sum(diag(tab))
disagree <- n_over - agree

cat(sprintf("\n  Agree (same category in both): %d (%.1f%%)\n",
            agree, 100 * agree / max(n_over, 1)))
cat(sprintf("  Disagree (different category): %d (%.1f%%)\n",
            disagree, 100 * disagree / max(n_over, 1)))

cat("\n=== Category cross-table (rows = pipe1, cols = pipe2) ===\n")
print(tab)

## --- Per-category counts in the overlap --------------------------------------
cat("\n=== Assign_pipe1 distribution in overlap ===\n")
print(sort(table(overlap$Assign_pipe1), decreasing = TRUE))

cat("\n=== Assign_pipe2 distribution in overlap ===\n")
print(sort(table(overlap$Assign_pipe2), decreasing = TRUE))

## --- Write overlap table to TSV ----------------------------------------------
write.table(
  overlap,
  file      = outfile,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
cat(sprintf("\nOverlap table written to: %s\n", outfile))
