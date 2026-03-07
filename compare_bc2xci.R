#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript compare_bc2xci.R file1_bc2xci.txt file2_bc2xci.txt")
}

file1 <- args[1]
file2 <- args[2]

## --- Read & prep -------------------------------------------------------------
read_bc2xci <- function(path) {
  df <- read.table(
    path,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )
  
  if (ncol(df) < 4) {
    stop(sprintf("File %s has fewer than 4 columns.", path))
  }
  
  df <- df[, 1:4]
  colnames(df) <- c("Barcode", "Val1", "Val2", "Assign")
  df$Barcode <- trimws(df$Barcode)
  df[!duplicated(df$Barcode), c("Barcode", "Assign")]
}

df1u <- read_bc2xci(file1)
df2u <- read_bc2xci(file2)

## Pick the shorter as template
if (nrow(df1u) <= nrow(df2u)) {
  df_short <- df1u
  df_long  <- df2u
  short_name <- basename(file1)
  long_name  <- basename(file2)
} else {
  df_short <- df2u
  df_long  <- df1u
  short_name <- basename(file2)
  long_name  <- basename(file1)
}

names(df_short)[2] <- "Assign_short"
names(df_long)[2]  <- "Assign_long"

## Left-join to keep all rows from the shorter file
cmp <- merge(df_short, df_long, by = "Barcode", all.x = TRUE)

## --- Coverage / agreement counts --------------------------------------------
n_short           <- nrow(cmp)
n_overlap         <- sum(!is.na(cmp$Assign_long))
missing_in_longer <- sum(is.na(cmp$Assign_long))
agree             <- sum(!is.na(cmp$Assign_long) & cmp$Assign_short == cmp$Assign_long)
disagree          <- sum(!is.na(cmp$Assign_long) & cmp$Assign_short != cmp$Assign_long)
both_is_both      <- sum(cmp$Assign_short == "Both" & cmp$Assign_long == "Both", na.rm = TRUE)
either_is_both    <- sum(cmp$Assign_short == "Both" | cmp$Assign_long == "Both", na.rm = TRUE)
extra_in_longer   <- sum(!df_long$Barcode %in% df_short$Barcode)

cat("Shorter/template file:", short_name, "\n")
cat("Longer/reference file:", long_name, "\n\n")

cat("Shorter/template rows:", n_short, "\n")
cat("Overlap rows:", n_overlap, "\n")
cat("Missing in longer:", missing_in_longer, "\n")
cat("Agree:", agree, "\n")
cat("Disagree:", disagree, "\n")
cat("'Both' in both:", both_is_both, "\n")
cat("'Both' in either:", either_is_both, "\n")
cat("Extra only in longer:", extra_in_longer, "\n\n")

## --- Build overlap-only table and drop LC/Unknown ----------------------------
overlap <- cmp[!is.na(cmp$Assign_long), ]
names(overlap)[2:3] <- c("Short", "Long")

dropcats <- c("Low_coverage", "LC", "Unknown")
overlap_f <- subset(overlap, !(Short %in% dropcats | Long %in% dropcats))

tab <- table(overlap_f$Short, overlap_f$Long)

cat("Contingency table after dropping Low_coverage / LC / Unknown:\n")
print(tab)

## --- Agreement effect size: Cohen's kappa (unweighted) -----------------------
cohen_kappa_from_table <- function(tab) {
  N <- sum(tab)
  if (N == 0) return(NA_real_)
  po <- sum(diag(tab)) / N
  pe <- sum(rowSums(tab) * colSums(tab)) / (N^2)
  (po - pe) / (1 - pe)
}

kappa <- cohen_kappa_from_table(tab)
cat(sprintf("\nCohen's kappa (unweighted): %.3f\n", kappa))

## --- Tests: McNemar (2x2) or Bowker (k>2) -----------------------------------
bowker_test <- function(tab) {
  if (nrow(tab) != ncol(tab)) stop("Bowker's test requires a square table.")
  
  k <- nrow(tab)
  stat <- 0
  df <- 0
  
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      nij <- tab[i, j]
      nji <- tab[j, i]
      s <- nij + nji
      if (s > 0) {
        stat <- stat + (nij - nji)^2 / s
        df <- df + 1
      }
    }
  }
  
  p <- pchisq(stat, df = df, lower.tail = FALSE)
  list(
    statistic = stat,
    parameter = df,
    p.value = p,
    method = "Bowker's test of symmetry"
  )
}

cat("\n== Agreement test ==\n")
if (nrow(tab) == 2 && ncol(tab) == 2) {
  mc <- mcnemar.test(tab, correct = FALSE)
  print(mc)
  
  b <- tab[1, 2]
  c <- tab[2, 1]
  if ((b + c) < 25) {
    cat("Note: sparse discordant counts; consider exact McNemar (exact2x2::mcnemar.exact).\n")
  }
} else {
  bw <- bowker_test(tab)
  print(bw)
  cat("Tip: For marginal homogeneity you can also run Stuart-Maxwell (DescTools::StuartMaxwellTest).\n")
}