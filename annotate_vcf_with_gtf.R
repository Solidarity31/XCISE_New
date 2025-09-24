#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
# tiny arg parser
opt <- list(vcf=NULL, gtf="gencode_v46_GRCh38", cache="~/proj_xci/refcache",
            out=NULL, nearest="TRUE")
for (a in args) {
  kv <- strsplit(a, "=", fixed=TRUE)[[1]]
  if (length(kv)==2) opt[[sub("^--","",kv[1])]] <- kv[2]
}
stopifnot(!is.null(opt$vcf))
opt$nearest <- toupper(opt$nearest) %in% c("1","TRUE","T","YES","Y")

# registry
url_map <- list(
  gencode_v46_GRCh38 = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz",
  ensembl_113_GRCh38 = "https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz"
)

vcf_path <- opt$vcf
hdr <- grep("^#CHROM", readLines(vcf_path, warn = FALSE))[1]
if (is.na(hdr)) stop("VCF missing #CHROM header")
vcf <- fread(vcf_path, sep="\t", quote="", skip=hdr, header=FALSE, fill=TRUE,
             col.names=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"),
             na.strings=c("","NA"))
vcf[, CHROM := sub("^chr","", CHROM)]
vcf[, POS := as.integer(POS)]

# acquire GTF (cache if URL token)
dir.create(path.expand(opt$cache), recursive=TRUE, showWarnings=FALSE)
if (opt$gtf %in% names(url_map)) {
  gtf_path <- file.path(path.expand(opt$cache), basename(url_map[[opt$gtf]]))
  if (!file.exists(gtf_path)) download.file(url_map[[opt$gtf]], gtf_path, mode="wb")
} else if (grepl("://", opt$gtf)) {
  gtf_path <- file.path(path.expand(opt$cache), basename(opt$gtf))
  if (!file.exists(gtf_path)) download.file(opt$gtf, gtf_path, mode="wb")
} else {
  gtf_path <- path.expand(opt$gtf)
  if (!file.exists(gtf_path)) stop("GTF not found: ", gtf_path)
}

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
    gtf <- fread(p, sep="\t", header=FALSE, comment.char="#",
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
vcf[, `:=`(gene=NA_character_, biotype=NA_character_)]
if (nrow(ann)) vcf[ann$i, `:=`(gene=ann$gene, biotype=ann$biotype)]

# nearest gene for intergenic sites
if (opt$nearest) {
  nohit <- setdiff(snv$i, ann$i)
  if (length(nohit)) {
    snv_miss <- snv[i %in% nohit, .(i, chrom, pos=start)]
    g_end <- genes[, .(chrom, start, end, gene, strand)]; setkey(g_end, chrom, end)
    up   <- g_end[snv_miss, on=.(chrom, end<=pos), mult="last", nomatch=NA][, dist_up := pos - end]
    g_start <- genes[, .(chrom, start, end, gene, strand)]; setkey(g_start, chrom, start)
    down <- g_start[snv_miss, on=.(chrom, start>=pos), mult="first", nomatch=NA][, dist_down := start - pos]
    
    nearest <- merge(
      up[, .(i, gene_up=gene, strand_up=strand, dist_up)],
      down[, .(i, gene_down=gene, strand_down=strand, dist_down)],
      by="i", all=TRUE
    )
    
    choose <- nearest[, {
      if (!is.na(dist_up) && (is.na(dist_down) || dist_up <= dist_down))
        .(nearest_gene=gene_up, nearest_strand=strand_up, nearest_side="upstream", nearest_dist_bp=as.integer(dist_up))
      else if (!is.na(dist_down))
        .(nearest_gene=gene_down, nearest_strand=strand_down, nearest_side="downstream", nearest_dist_bp=as.integer(dist_down))
      else
        .(nearest_gene=NA_character_, nearest_strand=NA_character_, nearest_side=NA_character_, nearest_dist_bp=NA_integer_)
    }, by=i]
    
    vcf[, c("nearest_gene","nearest_strand","nearest_side","nearest_dist_bp") :=
          .(NA_character_, NA_character_, NA_character_, NA_integer_)]
    if (nrow(choose))
      vcf[choose$i, `:=`(nearest_gene=choose$nearest_gene,
                         nearest_strand=choose$nearest_strand,
                         nearest_side=choose$nearest_side,
                         nearest_dist_bp=choose$nearest_dist_bp)]
  }
}

if (is.null(opt$out)) {
  base <- tools::file_path_sans_ext(basename(vcf_path))
  opt$out <- file.path(dirname(vcf_path), paste0(base, ".annot.tsv"))
}
fwrite(vcf, opt$out, sep="\t", na="NA")
cat(sprintf("Annotated %d SNVs; wrote %s\n", nrow(vcf), opt$out))
