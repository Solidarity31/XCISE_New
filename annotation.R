#!/usr/bin/env Rscript
# =============================================================================
#  annotation.R  -  XCISE/EXCISE VCF annotator + chrX haplotype UMI plots
#
#  Required arguments:
#    --vcf=<path>      XCISE/EXCISE output VCF (INFO has AD=ref,alt and
#                      X1A=Ref|Alt|Unk). Works with the original XCISE VCF
#                      (X1A + AD only) and with the excise1.5 VCF (adds X1, X2,
#                      DP, MONO).
#    --fasta=<path>    Reference genome FASTA (REF allele validation)
#    --gtf=<path>      Gene annotation GTF (gene/biotype interval overlap)
#
#  Optional arguments:
#    --xci_ref=<path>  XCI reference table (.xlsx, .tsv or .csv).
#                      Columns are read by position (header row required,
#                      header names ignored):
#                        Col 1 : gene symbol
#                        Col 2 : Ensembl gene ID
#                        Col 3 : chromosome ("X" or "chrX")
#                        Col 4 : start (1-based integer, same build as GTF)
#                        Col 5 : end   (1-based integer)
#                        Col 6 : transcript type / biotype
#                        Col 7 : XCI status (escape/variable/inactive/par)
#                      Two-pass annotation:
#                        Pass 1 - gene-name join for GTF-annotated SNVs.
#                        Pass 2 - coordinate overlap (cols 3-5) for SNVs that
#                                 are still "unknown" after pass 1.
#                      If omitted, the xci_status column is absent.
#    --hap=<path>      excise1.5 haplotype TSV (*_XCISE_haplotypes.tsv:
#                      CHROM POS ID REF ALT X1_allele X2_allele mono).
#                      It is merged with the VCF by CHROM+POS: the haplotype
#                      TSV gives the X1/X2 allele, the VCF gives the AD UMI
#                      counts and the X1A tag. The two phasings are
#                      cross-checked and the merged table is written to
#                      <out>.haplotype_merged.tsv. Without --hap, the phase
#                      comes from the per-SNV VCF X1A tag alone.
#    --plot=TRUE       Draw chrX haplotype UMI plots (default: FALSE)
#    --out=<path>      Output TSV path (default: <vcf_basename>.annot.tsv)
#    --nearest=FALSE   Annotate intergenic SNVs to nearest gene (default: TRUE)
#    --summary=TRUE    Write <out>.summary.txt + <out>.gene_summary.tsv
#                      (default: FALSE)
#
#  Example:
#    Rscript annotation.R \
#      --vcf=Z1_chrX_XCISE.vcf \
#      --hap=Z1_chrX_XCISE_haplotypes.tsv \
#      --fasta=/path/to/genome.fa \
#      --gtf=/path/to/genes.gtf \
#      --xci_ref=excise_xci_reference.tsv \
#      --summary=TRUE --plot=TRUE \
#      --out=Z1_output.tsv
#
#  Outputs (<base> = --out without .tsv):
#    <out>                          annotated SNV table (columns below)
#    <base>.haplotype_merged.tsv    VCF x haplotype-TSV merge   (--hap)
#    <base>.summary.txt             human-readable summary      (--summary)
#    <base>.gene_summary.tsv        per-gene table              (--summary)
#    <base>_chrX_haplotype.{pdf,png}           UMI along chrX   (--plot)
#    <base>_chrX_haplotype_umi_dist.{pdf,png}  X1 vs X2 UMIs    (--plot)
#    <base>_chrX_haplotype_x1frac.{pdf,png}    per-SNV X1 frac  (--plot)
#
#  Annotated table columns:
#    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO,
#    [x1_allele, x2_allele]          (--hap)
#    gene_merged, biotype_merged,
#    [xci_status]                    (--xci_ref)
#    [nearest_gene, nearest_strand, nearest_side, nearest_dist_bp] (--nearest)
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── 0. Helpers ───────────────────────────────────────────────────────────────
is_true <- function(x) toupper(x) %in% c("1", "TRUE", "T", "YES", "Y")
warn    <- function(...) message("WARNING: ", ...)

# Blank / "NA" / "nan" strings (e.g. from Excel or pandas exports) -> NA
clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "na", "nan", "NaN", "N/A")] <- NA_character_
  x
}

xci_levels   <- c("escape", "variable", "inactive", "par")
xci_priority <- c(escape = 3L, variable = 2L, inactive = 1L, par = 0L)
prio_of <- function(s) {
  p <- unname(xci_priority[s])
  p[is.na(p)] <- -1L
  p
}

# Value of key=value tag in a VCF INFO vector (NA if absent)
info_tag <- function(info, tag) {
  out <- rep(NA_character_, length(info))
  hit <- grepl(paste0("(^|;)", tag, "="), info)
  out[hit] <- sub(paste0("^(.*;)?", tag, "=([^;]*).*$"), "\\2", info[hit])
  out
}
info_flag <- function(info, tag) grepl(paste0("(^|;)", tag, "(;|$)"), info)

# Merge GTF and XCI-reference biotypes for one gene (vectorised)
merge_bio <- function(gtf_b, xci_b) {
  out <- fifelse(gtf_b == xci_b, gtf_b, paste0(gtf_b, "/", xci_b))
  out <- fifelse(is.na(xci_b), gtf_b, out)
  out <- fifelse(is.na(gtf_b), xci_b, out)
  fifelse(is.na(gtf_b) & is.na(xci_b), "unknown", out)
}

# ── 1. Arguments ─────────────────────────────────────────────────────────────
known <- c("vcf", "fasta", "gtf", "xci_ref", "hap", "out",
           "nearest", "summary", "plot")
opt <- list(nearest = "TRUE", summary = "FALSE", plot = "FALSE")
for (a in commandArgs(trailingOnly = TRUE)) {
  if (!grepl("^--[^=]+=", a)) { warn("ignoring malformed argument: ", a); next }
  key <- sub("^--([^=]+)=.*$", "\\1", a)
  if (!key %in% known) { warn("ignoring unknown argument: --", key); next }
  opt[[key]] <- sub("^--[^=]+=", "", a)   # split on the first '=' only
}

for (k in c("vcf", "fasta", "gtf"))
  if (is.null(opt[[k]])) stop("--", k, " is required")

opt$nearest <- is_true(opt$nearest)
opt$summary <- is_true(opt$summary)
opt$plot    <- is_true(opt$plot)

vcf_path   <- path.expand(opt$vcf)
fasta_path <- path.expand(opt$fasta)
gtf_path   <- path.expand(opt$gtf)
for (p in c(vcf_path, fasta_path, gtf_path))
  if (!file.exists(p)) stop("File not found: ", p)

use_xci_ref <- !is.null(opt$xci_ref)
xci_path    <- if (use_xci_ref) path.expand(opt$xci_ref) else NA_character_
if (use_xci_ref && !file.exists(xci_path)) stop("xci_ref not found: ", xci_path)

use_hap  <- !is.null(opt$hap)
hap_path <- if (use_hap) path.expand(opt$hap) else NA_character_
if (use_hap && !file.exists(hap_path)) stop("hap file not found: ", hap_path)

if (is.null(opt$out)) {
  base    <- tools::file_path_sans_ext(basename(vcf_path))
  opt$out <- file.path(dirname(vcf_path), paste0(base, ".annot.tsv"))
}
# All side outputs derive from out_base, so they never overwrite --out
out_base <- sub("\\.tsv$", "", opt$out)

# ── 2. XCI reference (optional) ──────────────────────────────────────────────
xci_name  <- NULL   # one row per gene name          (pass 1)
xci_coord <- NULL   # interval table                 (pass 2)

if (use_xci_ref) {
  ext <- tolower(tools::file_ext(xci_path))
  raw <- if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("readxl required for Excel input: install.packages('readxl')")
    as.data.table(readxl::read_excel(xci_path, sheet = 1L, col_names = TRUE))
  } else {
    fread(xci_path, sep = if (ext == "csv") "," else "\t",
          header = TRUE, fill = TRUE, colClasses = "character")
  }
  if (ncol(raw) < 7)
    stop(sprintf("xci_ref must have at least 7 columns (found %d)", ncol(raw)))

  xci_raw <- data.table(
    gene       = clean_chr(raw[[1]]),
    chrom      = sub("^chr", "", clean_chr(raw[[3]])),
    start      = suppressWarnings(as.integer(raw[[4]])),
    end        = suppressWarnings(as.integer(raw[[5]])),
    biotype    = clean_chr(raw[[6]]),
    xci_status = tolower(clean_chr(raw[[7]]))
  )[!is.na(gene) & !is.na(xci_status)]
  xci_raw[, prio := prio_of(xci_status)]

  cat(sprintf("XCI reference: %d genes (%s)\n", nrow(xci_raw),
              paste(sprintf("%s=%d", c(xci_levels, "other"),
                            c(sapply(xci_levels, function(s) sum(xci_raw$xci_status == s)),
                              sum(!xci_raw$xci_status %in% xci_levels))),
                    collapse = ", ")))
  other <- xci_raw[!xci_status %in% xci_levels]
  if (nrow(other)) {
    warn(nrow(other), " row(s) with unrecognised xci_status:")
    print(other[, .(gene, xci_status)])
  }

  # A gene listed twice keeps its highest-priority status
  n_dup <- sum(duplicated(xci_raw$gene))
  if (n_dup) warn(n_dup, " duplicated gene name(s) in xci_ref; keeping highest-priority status")
  xci_name <- xci_raw[order(gene, -prio)][!duplicated(gene),
                                          .(gene, xci_bio = biotype, xci_status)]

  xci_coord <- xci_raw[!is.na(chrom) & !is.na(start) & !is.na(end) &
                         start > 0 & end >= start,
                       .(chrom, start, end, gene, biotype, xci_status, prio)]
  setkey(xci_coord, chrom, start, end)
  cat(sprintf("  %d intervals indexed for pass-2 coordinate fallback\n", nrow(xci_coord)))
} else {
  message("NOTE: --xci_ref not provided; XCI status will not be annotated.")
}

# ── 3. Read VCF (vectorised) ─────────────────────────────────────────────────
# Only the first 8 fixed VCF columns are used. Meta lines are dropped wherever
# they occur (excise1.5 writes its ##INFO lines after #CHROM).
lines <- readLines(vcf_path, warn = FALSE)
if (!any(startsWith(lines, "#CHROM"))) stop("VCF missing #CHROM header line")
body  <- lines[!startsWith(lines, "#") & nzchar(lines)]
rm(lines)

vcf_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
if (length(body) == 0L) {
  vcf <- setnames(as.data.table(rep(list(character(0)), 8)), vcf_cols)
} else {
  n_tab <- nchar(body) - nchar(gsub("\t", "", body, fixed = TRUE))
  if (any(n_tab < 6L))
    stop("VCF data line with < 7 fields, e.g.: ", body[which(n_tab < 6L)[1]])
  # 7 fields = FILTER missing (older XCISE output): insert PASS before INFO
  if (any(n_tab == 6L)) {
    message("NOTE: ", sum(n_tab == 6L), " VCF row(s) lack FILTER; inserting PASS.")
    body[n_tab == 6L] <- sub("\t([^\t]*)$", "\tPASS\t\\1", body[n_tab == 6L])
  }
  vcf <- fread(text = body, sep = "\t", header = FALSE, quote = "",
               colClasses = "character", fill = TRUE, na.strings = NULL)
  if (ncol(vcf) > 8L) {
    message("NOTE: dropping ", ncol(vcf) - 8L, " FORMAT/sample column(s) from VCF.")
    vcf <- vcf[, 1:8]
  }
  setnames(vcf, vcf_cols)
}
rm(body)

vcf[, POS   := suppressWarnings(as.integer(POS))]
vcf[, CHROM := sub("^chr", "", CHROM)]
n_before <- nrow(vcf)
vcf <- vcf[!is.na(POS) & !is.na(CHROM) & CHROM != ""]
if (nrow(vcf) < n_before)
  message(sprintf("Dropped %d VCF row(s) with invalid CHROM/POS", n_before - nrow(vcf)))
if (nrow(vcf) == 0L) stop("No variants in VCF")
cat(sprintf("VCF: %d variants\n", nrow(vcf)))

# ── 4. Parse INFO: X1A, AD, MONO ─────────────────────────────────────────────
ad <- info_tag(vcf$INFO, "AD")
ad[!grepl(",", ad, fixed = TRUE)] <- NA_character_
vcf[, `:=`(
  X1A     = info_tag(INFO, "X1A"),
  ref_umi = suppressWarnings(as.integer(sub(",.*$", "", ad))),
  alt_umi = suppressWarnings(as.integer(sub("^[^,]*,", "", ad))),
  is_mono = info_flag(INFO, "MONO")
)]
rm(ad)

# ── 5. FASTA: REF allele validation ──────────────────────────────────────────
if (requireNamespace("Rsamtools",     quietly = TRUE) &&
    requireNamespace("GenomicRanges", quietly = TRUE) &&
    requireNamespace("IRanges",       quietly = TRUE)) {
  tryCatch({
    if (!file.exists(paste0(fasta_path, ".fai"))) {
      message("FASTA index not found; creating ...")
      Rsamtools::indexFa(fasta_path)
    }
    fa <- Rsamtools::FaFile(fasta_path)
    open(fa)
    fa_seqs <- as.character(GenomicRanges::seqnames(Rsamtools::scanFaIndex(fa)))
    chr_pfx <- if (any(grepl("^chr", fa_seqs))) "chr" else ""
    gr <- GenomicRanges::GRanges(
      seqnames = paste0(chr_pfx, vcf$CHROM),
      ranges   = IRanges::IRanges(start = vcf$POS, width = nchar(vcf$REF)))
    valid   <- as.character(GenomicRanges::seqnames(gr)) %in% fa_seqs
    ref_seq <- rep(NA_character_, nrow(vcf))
    if (any(valid))
      ref_seq[valid] <- as.character(Rsamtools::scanFa(fa, gr[valid]))
    close(fa)

    # toupper: soft-masked (lower-case) FASTA bases are not mismatches
    mm <- which(!is.na(ref_seq) & toupper(vcf$REF) != toupper(ref_seq))
    cat(sprintf("FASTA REF check: %d checked, %d mismatches, %d contig-absent\n",
                sum(valid), length(mm), sum(!valid)))
    if (length(mm)) {
      warn("REF mismatches (first 5 shown) - check genome build:")
      print(head(data.table(CHROM = vcf$CHROM[mm], POS = vcf$POS[mm],
                            REF_vcf = vcf$REF[mm], REF_fasta = ref_seq[mm]), 5))
    }
  }, error = function(e)
    message("FASTA REF validation skipped: ", conditionMessage(e)))
} else {
  message("NOTE: Rsamtools not available; FASTA REF validation skipped.")
}

# ── 6. GTF import (gene records only) ────────────────────────────────────────
gtf_attr <- function(attr, key) {
  out <- rep(NA_character_, length(attr))
  hit <- grepl(paste0(key, ' "'), attr, fixed = TRUE)
  out[hit] <- sub(paste0('^.*', key, ' "([^"]*)".*$'), "\\1", attr[hit])
  out
}

import_gtf <- function(p) {
  if (requireNamespace("rtracklayer", quietly = TRUE)) {
    df <- as.data.frame(rtracklayer::import(p, format = "gtf", feature.type = "gene"))
    pick <- function(cols) {
      cols <- intersect(cols, names(df))
      if (!length(cols)) return(rep(NA_character_, nrow(df)))
      do.call(fcoalesce, lapply(cols, function(cl) as.character(df[[cl]])))
    }
    dt <- data.table(chrom   = sub("^chr", "", as.character(df$seqnames)),
                     start   = as.integer(df$start),
                     end     = as.integer(df$end),
                     strand  = as.character(df$strand),
                     gene    = pick(c("gene_name", "Name", "gene_id")),
                     biotype = pick(c("gene_type", "gene_biotype")))
  } else {
    # Skip the leading '#' comment block (fread's skip="#" would instead jump
    # to the first line *containing* '#')
    head_lines <- readLines(p, n = 1000L, warn = FALSE)
    n_skip <- match(FALSE, startsWith(head_lines, "#"), nomatch = 1L) - 1L
    gtf <- fread(p, sep = "\t", header = FALSE, quote = "", skip = n_skip,
                 select = c(1L, 3L, 4L, 5L, 7L, 9L),
                 col.names = c("seqname", "type", "start", "end", "strand", "attr"))
    gtf <- gtf[type == "gene"]
    dt <- gtf[, .(chrom   = sub("^chr", "", seqname),
                  start   = as.integer(start),
                  end     = as.integer(end),
                  strand,
                  # Some Ensembl genes lack gene_name: fall back to gene_id
                  gene    = fcoalesce(gtf_attr(attr, "gene_name"),
                                      gtf_attr(attr, "gene_id")),
                  biotype = fcoalesce(gtf_attr(attr, "gene_type"),
                                      gtf_attr(attr, "gene_biotype")))]
  }
  dt
}

genes <- import_gtf(gtf_path)
genes <- genes[chrom %in% unique(vcf$CHROM) &
                 !is.na(start) & !is.na(end) & start > 0 & end >= start]
setkey(genes, chrom, start, end)
cat(sprintf("GTF: %d genes on VCF chromosome(s)\n", nrow(genes)))

# ── 7. Gene + XCI annotation (vectorised, one row per SNV x gene) ────────────
snv <- vcf[, .(i = .I, chrom = CHROM, start = POS, end = POS)]
setkey(snv, chrom, start, end)

ov <- foverlaps(snv, genes, type = "within", nomatch = 0L)
setorder(ov, i, start)
in_gene <- seq_len(nrow(vcf)) %in% ov$i

# Gene and biotype stay paired per row (pasting unique() of each separately
# would misalign them when two overlapping genes share a biotype)
long <- unique(ov[!is.na(gene), .(i, gene, gtf_bio = biotype)], by = c("i", "gene"))
long[, `:=`(label = gene, xci_bio = NA_character_, xci_status = NA_character_)]

if (use_xci_ref) {
  # Pass 1: gene-name join
  long[xci_name, on = "gene", `:=`(xci_bio = i.xci_bio, xci_status = i.xci_status)]

  # Pass 2: coordinate overlap for SNVs without any status after pass 1
  unk_i <- setdiff(seq_len(nrow(vcf)), long[!is.na(xci_status), i])
  cat(sprintf("Pass 2: %d SNVs still unknown after name lookup\n", length(unk_i)))
  if (length(unk_i) && nrow(xci_coord)) {
    ov2 <- foverlaps(snv[i %in% unk_i], xci_coord, type = "within", nomatch = 0L)
    best2 <- ov2[order(i, -prio)][!duplicated(i),
                                  .(i, xci_gene = gene, xci_bio2 = biotype,
                                    xci_status2 = xci_status)]
    cat(sprintf("  resolved %d SNVs via coordinate overlap\n", nrow(best2)))
    # GTF-annotated SNVs: "gtfGene|xciGene" when names differ
    long[best2, on = "i", `:=`(
      label      = fifelse(gene == i.xci_gene, gene, paste0(gene, "|", i.xci_gene)),
      xci_bio    = i.xci_bio2,
      xci_status = i.xci_status2)]
    # Intergenic SNVs: the xci_ref gene is the only source
    long <- rbind(long,
                  best2[!i %in% long$i,
                        .(i, gene = xci_gene, gtf_bio = NA_character_, label = xci_gene,
                          xci_bio = xci_bio2, xci_status = xci_status2)])
    setorder(long, i)   # stable: keeps gene order within an SNV
  }
  long[is.na(xci_status), xci_status := "unknown"]
}
long[, bio := merge_bio(gtf_bio, xci_bio)]

snv_ann <- long[, .(gene_merged    = paste(label, collapse = ";"),
                    biotype_merged = paste0(label, ":", bio, collapse = ";"),
                    xci_status     = xci_status[which.max(prio_of(xci_status))]),
                by = i]

vcf[, `:=`(gene_merged = "unannotated", biotype_merged = "unknown")]
vcf[snv_ann$i, `:=`(gene_merged = snv_ann$gene_merged,
                    biotype_merged = snv_ann$biotype_merged)]
if (use_xci_ref) {
  vcf[, xci_status := "unknown"]
  vcf[snv_ann$i, xci_status := snv_ann$xci_status]
  cat("XCI status (per SNV): ",
      paste(sprintf("%s=%d", c(xci_levels, "unknown"),
                    sapply(c(xci_levels, "unknown"), function(s) sum(vcf$xci_status == s))),
            collapse = ", "), "\n", sep = "")
}

# ── 8. Nearest gene for intergenic SNVs (optional) ───────────────────────────
if (opt$nearest) {
  vcf[, `:=`(nearest_gene = NA_character_, nearest_strand = NA_character_,
             nearest_side = NA_character_, nearest_dist_bp = NA_integer_)]
  miss <- snv[!i %in% ov$i, .(i, chrom, pos = start)]   # snv is keyed (re-ordered): select by i
  if (nrow(miss)) {
    # x.end / x.start: after a join the join column holds i's values, so the
    # gene coordinate must be taken explicitly from x
    up <- genes[, .(chrom, end, gene, strand)][
      miss, on = .(chrom, end = pos), roll = Inf, mult = "last",
      .(i, gene_up = gene, strand_up = strand, dist_up = i.pos - x.end)]
    down <- genes[, .(chrom, start, gene, strand)][
      miss, on = .(chrom, start = pos), roll = -Inf, mult = "first",
      .(i, gene_down = gene, strand_down = strand, dist_down = x.start - i.pos)]
    nr <- merge(up, down, by = "i")
    nr[, use_up := !is.na(dist_up) & (is.na(dist_down) | dist_up <= dist_down)]
    vcf[nr$i, `:=`(
      nearest_gene    = fifelse(nr$use_up, nr$gene_up,   nr$gene_down),
      nearest_strand  = fifelse(nr$use_up, nr$strand_up, nr$strand_down),
      nearest_side    = fcase(nr$use_up, "upstream",
                              !is.na(nr$dist_down), "downstream",
                              default = NA_character_),
      nearest_dist_bp = as.integer(fifelse(nr$use_up, nr$dist_up, nr$dist_down)))]
  }
}

# ── 9. Haplotype table: VCF (X1A + AD) merged with excise1.5 haplotype TSV ───
#  Per SNV, the phase decides which AD count goes to which haplotype:
#    X1 allele = REF  ->  x1_umi = ref_umi, x2_umi = alt_umi
#    X1 allele = ALT  ->  x1_umi = alt_umi, x2_umi = ref_umi
#  The phase is taken from the haplotype TSV when the SNV is in it, otherwise
#  from the VCF X1A tag. phase_check records how the two sources agree:
#    concordant       both phased, same orientation
#    discordant       both phased, opposite orientation
#    hap_phased_only  in both files but VCF X1A=Unk
#    unassigned       VCF X1A=Unk and absent from TSV (escape candidate)
#    vcf_only         VCF phased but absent from TSV
#    hap_only         in TSV but absent from VCF (no UMI counts)
#    allele_mismatch  REF/ALT differ between the files (phase not trusted)
hap <- vcf[, .(CHROM, POS, ID, REF, ALT, X1A, ref_umi, alt_umi, is_mono,
               in_vcf = TRUE)]

if (use_hap) {
  h <- fread(hap_path, sep = "\t", header = TRUE, colClasses = "character")
  req <- c("POS", "X1_allele", "X2_allele")
  if (!all(req %in% names(h)))
    stop("hap file must contain columns ", paste(req, collapse = ", "),
         " (found: ", paste(names(h), collapse = ", "), ")")
  h <- h[, .(CHROM     = if ("CHROM" %in% names(h)) sub("^chr", "", CHROM)
                         else unique(vcf$CHROM)[1],
             POS       = as.integer(POS),
             hap_REF   = if ("REF" %in% names(h)) REF else NA_character_,
             hap_ALT   = if ("ALT" %in% names(h)) ALT else NA_character_,
             X1_allele, X2_allele,
             hap_mono  = if ("mono" %in% names(h)) as.integer(mono) else NA_integer_,
             in_hap    = TRUE)]

  hap <- merge(hap, h, by = c("CHROM", "POS"), all = TRUE)
  hap[is.na(in_vcf), in_vcf := FALSE]
  hap[is.na(in_hap), in_hap := FALSE]
  hap[in_vcf == FALSE, `:=`(REF = hap_REF, ALT = hap_ALT)]
  hap[, allele_ok := !in_vcf | !in_hap | is.na(hap_REF) |
                     (hap_REF == REF & hap_ALT == ALT)]
  hap[, hap_x1a := fcase(X1_allele == REF & X2_allele == ALT, "Ref",
                         X1_allele == ALT & X2_allele == REF, "Alt",
                         default = NA_character_)]
  hap[, phase_check := fcase(
    in_hap & !allele_ok,                          "allele_mismatch",
    in_vcf & in_hap & hap_x1a == X1A,             "concordant",
    in_vcf & in_hap & X1A %chin% c("Ref", "Alt"), "discordant",
    in_vcf & in_hap,                              "hap_phased_only",
    in_vcf & X1A %chin% "Unk",                    "unassigned",
    in_vcf,                                       "vcf_only",
    default =                                     "hap_only")]
  hap[, `:=`(phase       = fifelse(in_hap, hap_x1a, X1A),
             phase_source = fifelse(in_hap & allele_ok, "hap_tsv", "vcf_X1A"))]
  hap[phase_check == "allele_mismatch", phase := NA_character_]
  hap[, mono := fcoalesce(hap_mono, as.integer(is_mono | ref_umi == 0L | alt_umi == 0L))]

  pc <- hap[, .N, by = phase_check]
  cat("Haplotype merge (VCF x TSV): ",
      paste(sprintf("%s=%d", pc$phase_check, pc$N), collapse = ", "), "\n", sep = "")
  n_conc <- sum(hap$phase_check == "concordant")
  n_disc <- sum(hap$phase_check == "discordant")
  if (n_disc > n_conc)
    warn("most phased SNVs are discordant between VCF and haplotype TSV. ",
         "X1/X2 labels are arbitrary per run - the two files probably come ",
         "from different runs with swapped X1/X2 labels.")
  else if (n_disc > 0)
    warn(n_disc, " SNV(s) phased in opposite orientation by VCF and TSV ",
         "(TSV orientation used).")
} else {
  hap[, `:=`(phase = X1A, phase_source = "vcf_X1A",
             mono = as.integer(is_mono | ref_umi == 0L | alt_umi == 0L))]
}

hap[!phase %chin% c("Ref", "Alt"), `:=`(phase = NA_character_, phase_source = NA_character_)]
hap[, `:=`(x1_allele = fcase(phase == "Ref", REF, phase == "Alt", ALT),
           x2_allele = fcase(phase == "Ref", ALT, phase == "Alt", REF),
           x1_umi    = fcase(phase == "Ref", ref_umi, phase == "Alt", alt_umi),
           x2_umi    = fcase(phase == "Ref", alt_umi, phase == "Alt", ref_umi))]
hap[, total_umi := x1_umi + x2_umi]
hap[, x1_frac   := fifelse(total_umi > 0, round(x1_umi / total_umi, 4), NA_real_)]
setorder(hap, CHROM, POS)

# Attach gene / XCI annotation to the haplotype table
hap[vcf, on = c("CHROM", "POS"), gene_merged := i.gene_merged]
if (use_xci_ref) hap[vcf, on = c("CHROM", "POS"), xci_status := i.xci_status]

if (use_hap) {
  vcf[hap, on = c("CHROM", "POS"), `:=`(x1_allele = i.x1_allele, x2_allele = i.x2_allele)]
  merged_file <- paste0(out_base, ".haplotype_merged.tsv")
  merged_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "X1A", "X1_allele", "X2_allele",
                   "phase_check", "phase_source", "mono",
                   "ref_umi", "alt_umi", "x1_umi", "x2_umi", "total_umi", "x1_frac",
                   "gene_merged", if (use_xci_ref) "xci_status")
  fwrite(hap[, ..merged_cols], merged_file, sep = "\t", na = "NA", quote = FALSE)
  cat(sprintf("Wrote merged haplotype table (%d SNVs) to %s\n", nrow(hap), merged_file))
}

# ── 10. Write annotated TSV ──────────────────────────────────────────────────
out_cols <- c(vcf_cols,
              if (use_hap) c("x1_allele", "x2_allele"),
              "gene_merged", "biotype_merged",
              if (use_xci_ref) "xci_status",
              if (opt$nearest) c("nearest_gene", "nearest_strand",
                                 "nearest_side", "nearest_dist_bp"))
fwrite(vcf[, ..out_cols], opt$out, sep = "\t", na = "NA", quote = FALSE)
cat(sprintf("Wrote %d annotated variants to %s\n", nrow(vcf), opt$out))

# ── 11. Summary (optional) ───────────────────────────────────────────────────
# Fixed-width text table: numbers right-aligned with thousands separators,
# text left-aligned, NA shown as "-".
fmt_table <- function(df, right = character(0), indent = "  ") {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  is_right <- vapply(df, is.numeric, logical(1)) | names(df) %in% right
  cells <- lapply(df, function(v) {
    s <- if (is.integer(v)) formatC(v, format = "d", big.mark = ",")
         else if (is.numeric(v))   # 1 decimal if that is exact, else 3
           formatC(v, format = "f",
                   digits = if (isTRUE(all(abs(v * 10 - round(v * 10)) < 1e-9, na.rm = TRUE))) 1 else 3)
         else as.character(v)
    s[is.na(v)] <- "-"
    s
  })
  w   <- mapply(function(h, s) max(nchar(c(h, s))), names(df), cells)
  pad <- function(s, w, r) formatC(s, width = if (r) w else -w)
  c(paste0(indent, paste(mapply(pad, names(df), w, is_right), collapse = "  ")),
    paste0(indent, paste(strrep("-", w), collapse = "  ")),
    if (nrow(df)) paste0(indent, do.call(paste, c(Map(pad, cells, w, is_right), sep = "  "))))
}
display_names <- function(x) {
  map <- c(gene = "Gene", snvs = "SNVs", first_pos = "First_pos", last_pos = "Last_pos",
           span_kb = "Span_kb", biotype = "Biotype", xci_status = "XCI_status",
           x1_umi = "X1_UMIs", x2_umi = "X2_UMIs", x1_frac = "X1_frac")
  fcoalesce(unname(map[x]), x)
}
pct <- function(x, n) if (n > 0) sprintf("%.1f%%", 100 * x / n) else "-"
kv  <- function(k, v) sprintf("  %-26s %s", k, v)
section <- function(title) c("", strrep("-", 78), paste0("  ", title), strrep("-", 78))

if (opt$summary) {
  summary_file <- paste0(out_base, ".summary.txt")
  gene_file    <- paste0(out_base, ".gene_summary.tsv")
  n_total      <- nrow(vcf)
  has_umi      <- any(!is.na(hap$x1_umi))

  # Per-gene table (an SNV overlapping two genes counts for both)
  long[, `:=`(CHROM = vcf$CHROM[i], POS = vcf$POS[i])]
  if (has_umi) long[hap, on = c("CHROM", "POS"), `:=`(x1_umi = i.x1_umi, x2_umi = i.x2_umi)]
  gene_sum <- long[, c(list(
      snvs      = uniqueN(i),
      first_pos = min(POS),
      last_pos  = max(POS),
      span_kb   = round((max(POS) - min(POS)) / 1e3, 1),
      biotype   = paste(unique(bio), collapse = ";")),
    if (use_xci_ref) list(xci_status = xci_status[1]),
    if (has_umi) list(x1_umi = sum(x1_umi, na.rm = TRUE),
                      x2_umi = sum(x2_umi, na.rm = TRUE))),
    by = .(gene = label)]
  if (has_umi) gene_sum[, x1_frac := fifelse(x1_umi + x2_umi > 0,
                                             x1_umi / (x1_umi + x2_umi), NA_real_)]
  setorder(gene_sum, -snvs, first_pos)
  fwrite(gene_sum, gene_file, sep = "\t", na = "NA", quote = FALSE)

  n_in_gene <- sum(in_gene)
  rep_lines <- c(
    strrep("=", 78),
    "  EXCISE ANNOTATION SUMMARY",
    strrep("=", 78),
    kv("Generated",       format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    kv("Input VCF",       vcf_path),
    kv("Haplotype TSV",   if (use_hap) hap_path else "none (phase from VCF X1A tag)"),
    kv("Reference FASTA", fasta_path),
    kv("Gene GTF",        gtf_path),
    kv("XCI reference",   if (use_xci_ref) xci_path else "none (GTF-only annotation)"),
    kv("Annotated table", opt$out),
    kv("Per-gene table",  gene_file),

    section("1. SNV OVERVIEW"),
    fmt_table(data.frame(
      Category = c("Total SNVs", "Inside a GTF gene", "Intergenic (no GTF gene)",
                   if (use_xci_ref) "  - placed by XCI-ref coordinates",
                   "Unique genes / loci"),
      Count = as.integer(c(n_total, n_in_gene, n_total - n_in_gene,
                           if (use_xci_ref) sum(!in_gene & vcf$gene_merged != "unannotated"),
                           nrow(gene_sum))),
      Percent = c(pct(n_total, n_total), pct(n_in_gene, n_total),
                  pct(n_total - n_in_gene, n_total),
                  if (use_xci_ref) pct(sum(!in_gene & vcf$gene_merged != "unannotated"), n_total),
                  "")),
      right = "Percent")
  )

  if (use_xci_ref) {
    st <- c(xci_levels, "unknown")
    rep_lines <- c(rep_lines, section("2. XCI STATUS"),
      fmt_table(data.frame(
        Status  = st,
        SNVs    = vapply(st, function(s) sum(vcf$xci_status == s), integer(1)),
        Percent = vapply(st, function(s) pct(sum(vcf$xci_status == s), n_total), ""),
        Genes   = vapply(st, function(s) sum(gene_sum$xci_status == s), integer(1))),
        right = "Percent"),
      "  SNV status = highest-priority status of its genes (escape > variable > inactive > par).")
    for (s in c("escape", "variable")) {
      g <- gene_sum[xci_status == s]
      if (nrow(g))
        rep_lines <- c(rep_lines, "", sprintf("  %s genes (SNVs):", tools::toTitleCase(s)),
                       strwrap(paste(sprintf("%s (%d)", g$gene, g$snvs), collapse = ", "),
                               width = 76, indent = 4, exdent = 4))
    }
  }

  bt <- long[, .(SNVs = uniqueN(i), Genes = uniqueN(label)), by = .(Biotype = bio)]
  setorder(bt, -SNVs)
  rep_lines <- c(rep_lines, section("3. BIOTYPE (top 10)"),
                 if (nrow(bt)) fmt_table(head(bt, 10)) else "  No SNV inside an annotated gene.")

  # Haplotype phasing
  x1a <- fcoalesce(vcf$X1A, "missing")
  ph_rows <- data.frame(
    Category = c("X1A = Ref  (X1 carries REF)", "X1A = Alt  (X1 carries ALT)",
                 "X1A = Unk  (unassigned / escape candidate)", "X1A tag missing",
                 "Monoallelic (one allele observed)"),
    SNVs = as.integer(c(sum(x1a == "Ref"), sum(x1a == "Alt"), sum(x1a == "Unk"),
                        sum(x1a == "missing"), sum(hap$mono[hap$in_vcf] == 1L, na.rm = TRUE))))
  ph_rows$Percent <- vapply(ph_rows$SNVs, pct, "", n = n_total)
  rep_lines <- c(rep_lines, section("4. HAPLOTYPE PHASING"),
                 fmt_table(ph_rows, right = "Percent"))
  if (use_hap) {
    pc_all <- c("concordant", "discordant", "hap_phased_only", "unassigned",
                "vcf_only", "hap_only", "allele_mismatch")
    rep_lines <- c(rep_lines, "", "  VCF X1A vs haplotype TSV:",
                   fmt_table(data.frame(
                     Check = pc_all,
                     SNVs  = vapply(pc_all, function(s) sum(hap$phase_check == s), integer(1)))))
  }
  if (has_umi) {
    ph   <- hap[!is.na(x1_umi)]
    bial <- ph[mono == 0L & total_umi > 0]
    sx1 <- sum(ph$x1_umi); sx2 <- sum(ph$x2_umi)
    rep_lines <- c(rep_lines, "", "  Allele-specific UMIs on phased SNVs:",
      kv("  Phased SNVs", formatC(nrow(ph), format = "d", big.mark = ",")),
      kv("  X1 UMIs (sum)", sprintf("%s  (%s)", formatC(sx1, format = "d", big.mark = ","),
                                    pct(sx1, sx1 + sx2))),
      kv("  X2 UMIs (sum)", sprintf("%s  (%s)", formatC(sx2, format = "d", big.mark = ","),
                                    pct(sx2, sx1 + sx2))),
      kv("  Median X1 fraction", if (nrow(bial))
        sprintf("%.3f  (biallelic SNVs, n=%d)", median(bial$x1_frac), nrow(bial)) else "-"),
      "  Note: X1/X2 labels are set per run so that X2 is active in more cells;",
      "  X2 > X1 is therefore expected and is not an independent finding.")
  }

  top_n <- min(30L, nrow(gene_sum))
  rep_lines <- c(rep_lines,
                 section(sprintf("5. PER-GENE TABLE (top %d of %d by SNV count; full table: %s)",
                                 top_n, nrow(gene_sum), basename(gene_file))),
                 if (top_n) fmt_table(setnames(head(gene_sum, top_n), display_names(names(gene_sum))))
                 else "  No annotated genes.",
                 "", strrep("=", 78))

  writeLines(rep_lines, summary_file)
  cat(sprintf("Summary written to %s\n          and %s\n", summary_file, gene_file))
}

# ── 12. chrX haplotype UMI plots (optional) ──────────────────────────────────
if (!opt$plot) {
  cat("Haplotype plot skipped (use --plot=TRUE to enable).\n")
  quit(save = "no", status = 0)
}
if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("ggplot2 required for plotting: install.packages('ggplot2')")
suppressPackageStartupMessages(library(ggplot2))

xchr <- hap[CHROM == "X" & !is.na(x1_umi) & total_umi > 0]
if (nrow(xchr) == 0L) {
  message("No phased chrX SNVs with UMI counts - skipping plots.")
  quit(save = "no", status = 0)
}
if (!use_xci_ref) xchr[, xci_status := "unknown"]
xchr[, snv_type := factor(fifelse(mono == 1L, "monoallelic", "biallelic"),
                          levels = c("biallelic", "monoallelic"))]
unk_dt <- hap[CHROM == "X" & is.na(phase) & in_vcf, .(pos_Mb = POS / 1e6)]

# Landmarks: XIST/XACT from the GTF when present; PAR/centromere by build,
# inferred from the XIST position (GRCh38 ~73.8 Mb, GRCh37 ~73.0 Mb)
gtf_pos <- function(sym) genes[chrom == "X" & gene == sym, start[1] / 1e6]
xist_gtf <- gtf_pos("XIST")
build <- if (length(xist_gtf) && !is.na(xist_gtf) && xist_gtf < 73.5) "GRCh37" else "GRCh38"
landmarks <- list(   # Mb
  GRCh38 = list(par1 = c(0.010, 2.781), par2 = c(155.701, 156.031),
                cen = c(58.10, 61.00), xist = 73.821, xact = 113.616),
  GRCh37 = list(par1 = c(0.060, 2.700), par2 = c(154.931, 155.261),
                cen = c(58.63, 61.63), xist = 73.040, xact = 112.983))
lm <- landmarks[[build]]
if (length(xist_gtf) && !is.na(xist_gtf)) lm$xist <- xist_gtf
xact_gtf <- gtf_pos("XACT")
if (length(xact_gtf) && !is.na(xact_gtf)) lm$xact <- xact_gtf

label_x1 <- "X1 (minor haplotype)"
label_x2 <- "X2 (major haplotype)"
pal <- setNames(c("#C0392B", "#2471A3"), c(label_x1, label_x2))

plot_dt <- melt(xchr[, .(POS, snv_type, x1_umi, x2_umi)],
                id.vars = c("POS", "snv_type"), measure.vars = c("x1_umi", "x2_umi"),
                variable.name = "haplotype", value.name = "umi")
plot_dt[, haplotype := factor(fifelse(haplotype == "x1_umi", label_x1, label_x2),
                              levels = c(label_x1, label_x2))]
plot_dt[, pos_Mb := POS / 1e6]

base_theme <- theme_bw(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey40", size = 9),
        plot.caption  = element_text(colour = "grey50", size = 7.5, hjust = 0),
        legend.position = "top",
        panel.grid.minor = element_line(colour = "grey94"))
src_str <- paste(basename(c(vcf_path, if (use_hap) hap_path)), collapse = " + ")

# 12a. UMI along chrX
y_ceiling <- max(5, ceiling(quantile(plot_dt$umi, 0.99) * 1.2))
n_clipped <- sum(plot_dt$umi > y_ceiling)
if (n_clipped) message(sprintf("Note: %d point(s) above y = %d not shown.", n_clipped, y_ceiling))
escape_dt <- xchr[xci_status == "escape", .(pos_Mb = POS / 1e6)]
span_v    <- max(0.05, min(0.35, 300 / nrow(xchr)))
n_mono    <- sum(xchr$snv_type == "monoallelic")

lab_text <- function(x, label, colour = "grey35", hjust = 0.5)
  annotate("text", x = x, y = Inf, vjust = 1.6, hjust = hjust, size = 2.6,
           label = label, colour = colour, fontface = "italic")

p_pos <- ggplot(plot_dt, aes(pos_Mb, umi, colour = haplotype)) +
  annotate("rect", xmin = lm$par1[1], xmax = lm$par1[2], ymin = -Inf, ymax = Inf,
           fill = "#F5CBA7", alpha = 0.35) +
  annotate("rect", xmin = lm$par2[1], xmax = lm$par2[2], ymin = -Inf, ymax = Inf,
           fill = "#F5CBA7", alpha = 0.35) +
  annotate("rect", xmin = lm$cen[1], xmax = lm$cen[2], ymin = -Inf, ymax = Inf,
           fill = "#D5D8DC", alpha = 0.45) +
  lab_text(mean(lm$par1), "PAR1") + lab_text(mean(lm$par2), "PAR2") +
  lab_text(mean(lm$cen), "CEN") +
  geom_vline(xintercept = lm$xist, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = lm$xact, linetype = "dotted", colour = "#8E44AD", linewidth = 0.5) +
  lab_text(lm$xist, "XIST", "grey40", hjust = 1.1) +
  lab_text(lm$xact, "XACT", "#8E44AD", hjust = -0.1) +
  { if (nrow(escape_dt))
      geom_rug(data = escape_dt, aes(x = pos_Mb), inherit.aes = FALSE, sides = "b",
               colour = "#27AE60", linewidth = 0.5, length = unit(0.04, "npc")) } +
  { if (nrow(unk_dt))
      geom_rug(data = unk_dt, aes(x = pos_Mb), inherit.aes = FALSE, sides = "t",
               colour = "grey55", linewidth = 0.3, length = unit(0.03, "npc")) } +
  geom_point(aes(shape = snv_type), alpha = 0.25, size = 0.8) +
  geom_smooth(method = "loess", formula = y ~ x, span = span_v,
              se = TRUE, linewidth = 1.1, alpha = 0.12) +
  scale_shape_manual(values = c(biallelic = 16, monoallelic = 17), name = NULL, drop = FALSE) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_x_continuous(name = "chrX position (Mb)", breaks = seq(0, 160, 10),
                     minor_breaks = seq(0, 160, 5), expand = c(0.005, 0)) +
  scale_y_continuous(name = "UMIs per SNV", breaks = pretty(c(0, y_ceiling), n = 6)) +
  coord_cartesian(ylim = c(0, y_ceiling)) +
  labs(title = "Haplotype-specific expression along chromosome X",
       subtitle = sprintf("%d phased SNVs (%d biallelic, %d monoallelic)  |  loess span %.2f  |  y clipped at %d  |  %s",
                          nrow(xchr), nrow(xchr) - n_mono, n_mono, span_v, y_ceiling, src_str),
       caption = paste0(
         sprintf("Landmarks %s. Dashed: XIST (%.2f Mb). Dotted: XACT (%.2f Mb). ",
                 build, lm$xist, lm$xact),
         "Orange: PAR1/PAR2. Grey band: centromere. ",
         if (nrow(escape_dt)) "Green rug: escape genes. " else "",
         if (nrow(unk_dt)) sprintf("Grey top rug: %d unphased SNVs (X1A=Unk). ", nrow(unk_dt)) else "",
         "Circles: biallelic; triangles: monoallelic.")) +
  base_theme

# 12b. Distribution of per-SNV UMIs on X1 vs X2
med <- plot_dt[, .(med = median(umi), n = .N), by = .(snv_type, haplotype)]
med[, lab := sprintf("median %g\nn = %d", med, n)]
bial <- xchr[snv_type == "biallelic"]
wt_str <- if (nrow(bial) >= 5) {
  p <- suppressWarnings(wilcox.test(bial$x2_umi, bial$x1_umi, paired = TRUE)$p.value)
  p_str <- format.pval(p, digits = 2)
  sprintf("Paired Wilcoxon (biallelic, X2 vs X1): p %s",
          if (startsWith(p_str, "<")) p_str else paste("=", p_str))
} else "Too few biallelic SNVs for a paired test"

p_dist <- ggplot(plot_dt, aes(haplotype, umi, fill = haplotype)) +
  geom_violin(data = plot_dt[, if (.N >= 3) .SD, by = .(snv_type, haplotype)],
              scale = "width", alpha = 0.35, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85) +
  geom_text(data = med, aes(label = lab, y = Inf), vjust = 1.3, size = 3,
            colour = "grey25") +
  facet_wrap(~ snv_type, drop = TRUE) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(name = "UMIs per SNV (log1p scale)", trans = "log1p",
                     breaks = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000),
                     expand = expansion(mult = c(0.02, 0.25))) +
  labs(x = NULL, title = "Per-SNV UMI distribution on the two X haplotypes",
       subtitle = sprintf("%s  |  %s", wt_str, src_str),
       caption = paste("Each phased SNV contributes one X1 and one X2 value.",
                       "Monoallelic SNVs have one haplotype at 0 by definition.")) +
  base_theme

# 12c. Per-SNV X1 fraction (biallelic only)
xci_cols <- c(escape = "#27AE60", variable = "#E67E22", inactive = "#7F8C8D",
              par = "#A04000", unknown = "#BDC3C7")
p_frac <- ggplot(bial, aes(x1_frac, fill = factor(xci_status, levels = names(xci_cols)))) +
  geom_histogram(breaks = seq(0, 1, 0.05), colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = xci_cols, name = "XCI status", drop = TRUE) +
  scale_x_continuous(name = "X1 fraction = X1 UMIs / (X1 + X2)", breaks = seq(0, 1, 0.1)) +
  labs(y = "Biallelic SNVs",
       title = "Allelic balance of phased biallelic SNVs",
       subtitle = sprintf("n = %d  |  median X1 fraction %.3f  |  %s", nrow(bial),
                          if (nrow(bial)) median(bial$x1_frac) else NA_real_, src_str)) +
  base_theme
if (!use_xci_ref) p_frac <- p_frac + guides(fill = "none")

# ── 13. Save plots ───────────────────────────────────────────────────────────
save_plot <- function(p, suffix, w, h) {
  f <- paste0(out_base, suffix, c(".pdf", ".png"))
  ggsave(f[1], p, width = w, height = h, device = "pdf")
  ggsave(f[2], p, width = w, height = h, dpi = 200, device = "png")
  cat(sprintf("  %s\n  %s\n", f[1], f[2]))
}
cat("chrX haplotype plots saved:\n")
save_plot(p_pos, "_chrX_haplotype", 16, 5)
save_plot(p_dist, "_chrX_haplotype_umi_dist", if (n_mono) 9 else 6, 5.5)
if (nrow(bial)) save_plot(p_frac, "_chrX_haplotype_x1frac", 8, 5)
