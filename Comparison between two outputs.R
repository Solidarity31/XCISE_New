## --- Read & prep -------------------------------------------------------------
df1 <- read.table("AD1_chrX_XCISE_bc2xci.txt", header = FALSE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
df2 <- read.table("Z1_chrX_XCISE_bc2xci.txt", header = FALSE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
colnames(df1) <- c("Barcode","Val1","Val2","Assign")
colnames(df2) <- c("Barcode","Val1","Val2","Assign")

df1$Barcode <- trimws(df1$Barcode); df2$Barcode <- trimws(df2$Barcode)
df1u <- df1[!duplicated(df1$Barcode), c("Barcode","Assign")]
df2u <- df2[!duplicated(df2$Barcode), c("Barcode","Assign")]

## Pick the shorter as template (for coverage stats)
df_short <- if (nrow(df1u) <= nrow(df2u)) df1u else df2u
df_long  <- if (identical(df_short, df1u)) df2u else df1u
names(df_short)[2] <- "Assign_short"
names(df_long)[2]  <- "Assign_long"

## Left-join to keep all rows from the shorter file
cmp <- merge(df_short, df_long, by = "Barcode", all.x = TRUE)

## Coverage / agreement counts relative to the shorter set
n_short           <- nrow(cmp)
n_overlap         <- sum(!is.na(cmp$Assign_long))
missing_in_longer <- sum(is.na(cmp$Assign_long))
agree             <- sum(!is.na(cmp$Assign_long) & cmp$Assign_short == cmp$Assign_long)
disagree          <- sum(!is.na(cmp$Assign_long) & cmp$Assign_short != cmp$Assign_long)
both_is_both      <- sum(cmp$Assign_short == "Both" & cmp$Assign_long == "Both", na.rm = TRUE)
either_is_both    <- sum(cmp$Assign_short == "Both" |  cmp$Assign_long == "Both", na.rm = TRUE)
extra_in_longer   <- sum(!df_long$Barcode %in% df_short$Barcode)

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
names(overlap)[2:3] <- c("Short","Long")

dropcats <- c("Low_coverage","LC","Unknown")
overlap_f <- subset(overlap, !(Short %in% dropcats | Long %in% dropcats))

tab <- table(overlap_f$Short, overlap_f$Long)
print(tab)

## --- Agreement effect size: Cohen's kappa (unweighted) -----------------------
cohen_kappa_from_table <- function(tab) {
  N  <- sum(tab)
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
  stat <- 0; df <- 0
  for (i in 1:(k-1)) {
    for (j in (i+1):k) {
      nij <- tab[i, j]
      nji <- tab[j, i]
      s <- nij + nji
      if (s > 0) {
        stat <- stat + ( (nij - nji)^2 / s )
        df <- df + 1
      }
    }
  }
  p <- pchisq(stat, df = df, lower.tail = FALSE)
  list(statistic = stat, parameter = df, p.value = p,
       method = "Bowker's test of symmetry")
}

cat("\n== Agreement test ==\n")
if (nrow(tab) == 2 && ncol(tab) == 2) {
  ## McNemar (paired binary)
  mc <- mcnemar.test(tab, correct = FALSE)
  print(mc)
  b <- tab[1,2]; c <- tab[2,1]
  if ((b + c) < 25) cat("Note: sparse discordant counts; consider exact McNemar (exact2x2::mcnemar.exact).\n")
} else {
  ## Multi-class nominal: Bowker
  bw <- bowker_test(tab)
  print(bw)
  cat("Tip: For marginal homogeneity you can also run Stuart–Maxwell (DescTools::StuartMaxwellTest).\n")
}
