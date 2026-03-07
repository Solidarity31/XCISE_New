#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
# Arguments: must specify vcf, fasta, and gtf manually
opt <- list(vcf=NULL, fasta=NULL, gtf=NULL, out=NULL, nearest="TRUE")
for (a in args) {
  kv <- strsplit(a, "=", fixed=TRUE)[[1]]
  if (length(kv)==2) opt[[sub("^--","",kv[1])]] <- kv[2]
}

# Check required arguments
stopifnot(!is.null(opt$vcf), !is.null(opt$fasta), !is.null(opt$gtf))
opt$nearest <- toupper(opt$nearest) %in% c("1","TRUE","T","YES","Y")

# Verify files exist
fasta_path <- path.expand(opt$fasta)
gtf_path <- path.expand(opt$gtf)
if (!file.exists(fasta_path)) stop("FASTA not found: ", fasta_path)
if (!file.exists(gtf_path)) stop("GTF not found: ", gtf_path)

vcf_path <- opt$vcf
hdr <- grep("^#CHROM", readLines(vcf_path, warn = FALSE))[1]
if (is.na(hdr)) stop("VCF missing #CHROM header")
vcf <- fread(vcf_path, sep="\t", quote="", skip=hdr, header=FALSE, fill=TRUE,
             col.names=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"),
             na.strings=c("","NA"))
vcf[, CHROM := sub("^chr","", CHROM)]
vcf[, POS := as.integer(POS)]

# import GTF (prefer rtracklayer, else fread)
# import GTF (prefer rtracklayer, else fread)
import_gtf <- function(p) {
  if (requireNamespace("rtracklayer", quietly=TRUE)) {
    g <- rtracklayer::import(p)
    df <- as.data.frame(g)
    nm <- names(df)
    name_col <- intersect(c("gene_name","Name","gene_id"), nm)[1]
    type_col <- intersect(c("gene_type","gene_biotype"), nm)[1]
    dt <- data.table(
      chrom   = sub("^chr","", as.character(df$seqnames)),
      start   = as.integer(df$start),
      end     = as.integer(df$end),
      strand  = as.character(df$strand),
      type    = if ("type" %in% nm) df[["type"]] else NA_character_,
      gene    = if (!is.na(name_col)) df[[name_col]] else NA_character_,
      biotype = if (!is.na(type_col)) df[[type_col]] else NA_character_
    )
    dt[type == "gene", .(chrom, start, end, strand, gene, biotype)]
  } else {
    gtf <- fread(p, sep="\t", header=FALSE, skip="#",
                 col.names=c("seqname","source","type","start","end","score","strand","frame","attr"))
    gtf <- gtf[type=="gene"]
    gtf[, gene_name := sub('.*gene_name "([^"]+)".*','\\1', attr)]
    gtf[, gene_type := fifelse(grepl('gene_type "', attr),
                               sub('.*gene_type "([^"]+)".*','\\1', attr),
                               fifelse(grepl('gene_biotype "', attr),
                                       sub('.*gene_biotype "([^"]+)".*','\\1', attr),
                                       NA_character_))]
    gtf[, .(chrom=sub("^chr","", seqname),
            start=as.integer(start), end=as.integer(end), strand,
            gene=gene_name, biotype=gene_type)]
  }
}

genes <- import_gtf(gtf_path)
genes <- genes[!is.na(start) & !is.na(end) & start>0 & end>=start]
setkey(genes, chrom, start, end)

snv <- data.table(i=seq_len(nrow(vcf)), chrom=vcf$CHROM, start=vcf$POS, end=vcf$POS)
setkey(snv, chrom, start, end)

# exact overlaps
ov <- foverlaps(snv, genes, type="within", nomatch=0L)
ann <- ov[, .(gene=paste(unique(na.omit(gene)), collapse=";"),
              biotype=paste(unique(na.omit(biotype)), collapse=";")), by=i]

# Initialize all variants as "unannotated"
vcf[, `:=`(gene="unannotated", biotype="unannotated")]
# Update with actual annotations where overlaps exist
if (nrow(ann)) vcf[ann$i, `:=`(gene=ann$gene, biotype=ann$biotype)]

# nearest gene for intergenic sites
if (opt$nearest) {
  nohit <- setdiff(snv$i, ann$i)
  if (length(nohit)) {
    snv_miss <- snv[i %in% nohit, .(i, chrom, pos = start)]
    
    # Prepare per-gene start/end tables
    g_end   <- genes[, .(chrom, end, gene, strand)]
    g_start <- genes[, .(chrom, start, gene, strand)]
    
    # Keys for rolling joins
    setkey(g_end,   chrom, end)
    setkey(g_start, chrom, start)
    setkey(snv_miss, chrom, pos)
    
    # Nearest upstream gene end: end <= pos  (roll to the last end before pos)
    up <- g_end[snv_miss, on = .(chrom, end = pos), roll = Inf,
                mult = "last", nomatch = NA]
    # After this join, 'up' has columns from g_end plus snv_miss's 'i' & 'pos'
    setnames(up, c("gene","strand"), c("gene_up","strand_up"))
    up[, dist_up := pos - end]
    
    # Nearest downstream gene start: start >= pos (roll to the first start after pos)
    down <- g_start[snv_miss, on = .(chrom, start = pos), roll = -Inf,
                    mult = "first", nomatch = NA]
    setnames(down, c("gene","strand"), c("gene_down","strand_down"))
    down[, dist_down := start - pos]
    
    # Merge candidates
    nearest <- merge(
      up[,   .(i, pos, gene_up,   strand_up,   dist_up)],
      down[, .(i, pos, gene_down, strand_down, dist_down)],
      by = c("i","pos"), all = TRUE
    )
    
    # Choose nearer side (prefer upstream on ties to be deterministic)
    choose <- nearest[, {
      if (!is.na(dist_up) && (is.na(dist_down) || dist_up <= dist_down)) {
        .(nearest_gene = gene_up,
          nearest_strand = strand_up,
          nearest_side = "upstream",
          nearest_dist_bp = as.integer(dist_up))
      } else if (!is.na(dist_down)) {
        .(nearest_gene = gene_down,
          nearest_strand = strand_down,
          nearest_side = "downstream",
          nearest_dist_bp = as.integer(dist_down))
      } else {
        .(nearest_gene = NA_character_,
          nearest_strand = NA_character_,
          nearest_side = NA_character_,
          nearest_dist_bp = NA_integer_)
      }
    }, by = i]
    
    # Attach to VCF table
    vcf[, c("nearest_gene","nearest_strand","nearest_side","nearest_dist_bp") :=
          .(NA_character_, NA_character_, NA_character_, NA_integer_)]
    if (nrow(choose)) {
      vcf[choose$i, `:=`(
        nearest_gene = choose$nearest_gene,
        nearest_strand = choose$nearest_strand,
        nearest_side = choose$nearest_side,
        nearest_dist_bp = choose$nearest_dist_bp
      )]
    }
  }
}

if (is.null(opt$out)) {
  base <- tools::file_path_sans_ext(basename(vcf_path))
  opt$out <- file.path(dirname(vcf_path), paste0(base, ".annot.tsv"))
}
# Generate summary statistics
summary_file <- sub("\\.tsv$", ".summary.txt", opt$out)

# Overall statistics
total_snvs <- nrow(vcf)
annotated_snvs <- sum(vcf$gene != "unannotated")
unannotated_snvs <- sum(vcf$gene == "unannotated")

# Per-gene statistics (for annotated SNVs only)
gene_stats <- vcf[gene != "unannotated", {
  # Handle cases where multiple genes are annotated (semicolon-separated)
  genes_list <- unlist(strsplit(gene, ";"))
  biotypes_list <- unlist(strsplit(biotype, ";"))
  
  # Create a row for each gene (in case of multi-gene annotations)
  data.table(
    gene = genes_list,
    biotype = biotypes_list,
    chrom = CHROM,
    pos = POS
  )
}]

# Aggregate by gene
gene_summary <- gene_stats[, .(
  n_snvs = .N,
  chromosome = paste(unique(chrom), collapse=","),
  min_pos = min(pos),
  max_pos = max(pos),
  interval_bp = max(pos) - min(pos),
  biotype = paste(unique(biotype), collapse=";")
), by = gene]

# Sort by number of SNVs (descending)
setorder(gene_summary, -n_snvs)

# Write summary file
sink(summary_file)
cat("=" , rep("=", 70), "\n", sep="")
cat("  VARIANT ANNOTATION SUMMARY STATISTICS\n")
cat("=" , rep("=", 70), "\n", sep="")
cat("\n")
cat("Input VCF:", vcf_path, "\n")
cat("Reference GTF:", gtf_path, "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")
cat("=" , rep("=", 70), "\n", sep="")
cat("  OVERALL STATISTICS\n")
cat("=" , rep("=", 70), "\n", sep="")
cat(sprintf("Total SNVs:        %d\n", total_snvs))
cat(sprintf("Annotated SNVs:    %d (%.2f%%)\n", annotated_snvs, 100*annotated_snvs/total_snvs))
cat(sprintf("Unannotated SNVs:  %d (%.2f%%)\n", unannotated_snvs, 100*unannotated_snvs/total_snvs))
cat(sprintf("Unique genes:      %d\n", nrow(gene_summary)))
cat("\n")
cat("=" , rep("=", 70), "\n", sep="")
cat("  PER-GENE STATISTICS\n")
cat("=" , rep("=", 70), "\n", sep="")
cat("\n")

if (nrow(gene_summary) > 0) {
  cat(sprintf("%-20s %8s %6s %12s %12s %12s %s\n", 
              "Gene", "N_SNVs", "Chr", "Min_Pos", "Max_Pos", "Interval_bp", "Biotype"))
  cat(rep("-", 100), "\n", sep="")
  
  for (i in 1:nrow(gene_summary)) {
    cat(sprintf("%-20s %8d %6s %12d %12d %12d %s\n",
                gene_summary$gene[i],
                gene_summary$n_snvs[i],
                gene_summary$chromosome[i],
                gene_summary$min_pos[i],
                gene_summary$max_pos[i],
                gene_summary$interval_bp[i],
                gene_summary$biotype[i]))
  }
} else {
  cat("No annotated SNVs found.\n")
}

cat("\n")
cat("=" , rep("=", 70), "\n", sep="")
cat("  TOP 20 GENES BY SNV COUNT\n")
cat("=" , rep("=", 70), "\n", sep="")
cat("\n")

if (nrow(gene_summary) > 0) {
  top20 <- head(gene_summary, 20)
  for (i in 1:nrow(top20)) {
    cat(sprintf("%2d. %-20s: %d SNVs (Chr %s, %d-%d bp, %d bp interval)\n",
                i,
                top20$gene[i],
                top20$n_snvs[i],
                top20$chromosome[i],
                top20$min_pos[i],
                top20$max_pos[i],
                top20$interval_bp[i]))
  }
}

sink()

cat(sprintf("Annotated %d SNVs; wrote %s\n", nrow(vcf), opt$out))
cat(sprintf("Summary statistics written to %s\n", summary_file))
fwrite(vcf, opt$out, sep="\t", na="NA")
cat(sprintf("Annotated %d SNVs; wrote %s\n", nrow(vcf), opt$out))