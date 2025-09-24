#!/usr/bin/env Rscript

# Minimal deps: base R + parallel + tools (for file_ext)
suppressWarnings(suppressMessages({
  # no external packages required
}))

args <- commandArgs(trailingOnly = TRUE)

die <- function(msg, usage = NULL) {
  cat(msg, "\n", file = stderr())
  if (!is.null(usage)) cat(usage, "\n", file = stderr())
  quit(status = 1)
}

usage <- "Usage: Rscript xcise_runner.R -o <output_prefix> -s <vcf_file> -b <bam1> [ -b <bam2> ... ] [-r <chromosome>] [-t <tries>] [-u <min_umis>] [-m <min_maf>] [-p <penalty>] [-j <jobs>] [-samthreads <N>] [-seed <int>] [-pretrain] [-pre_cap <N>] [-pre_wmin <W>] [-pre_tieskip] [-sm] [-ms <N>] [-x] [-i]"

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------- Parse CLI (manual, Perl-like) ----------------------
i <- 1
bams <- character()
sample <- NULL
snv_file <- NULL
chromosome <- "X"
tries <- 100L
min_maf <- 0
min_umis <- 10L
discordant_penalty <- 5L
scramble_alleles <- FALSE
improve_existing <- FALSE
no_softmasking <- FALSE
monostretch <- 15L

jobs <- 1L
samthreads <- 1L
seed <- NA_integer_

pretrain <- FALSE
pre_cap_pairs <- 5000L
pre_min_weight <- 3L
pre_tie_skip <- TRUE

expect_next <- function() {
  if (i >= length(args)) die("Unexpected/incomplete parameter", usage)
  args[i + 1]
}

while (i <= length(args)) {
  a <- args[i]
  if (a == "-b") {
    j <- i + 1
    while (j <= length(args) && !grepl("^-", args[j])) {
      if (!file.exists(args[j])) die(paste("BAM file", args[j], "does not exist"))
      bams <- c(bams, args[j]); j <- j + 1
    }
    i <- j
  } else if (a == "-o") { sample <- expect_next(); i <- i + 2
  } else if (a == "-s") { snv_file <- expect_next(); i <- i + 2
  } else if (a == "-r") { chromosome <- expect_next(); i <- i + 2
  } else if (a == "-t") { tries <- as.integer(expect_next()); i <- i + 2
  } else if (a == "-u") { min_umis <- as.integer(expect_next()); i <- i + 2
  } else if (a == "-m") { min_maf <- as.numeric(expect_next()); i <- i + 2
  } else if (a == "-p") { discordant_penalty <- as.integer(expect_next()); i <- i + 2
  } else if (a == "-x") { scramble_alleles <- TRUE; i <- i + 1
  } else if (a == "-i") { improve_existing <- TRUE; i <- i + 1
  } else if (a == "-sm") { no_softmasking <- TRUE; i <- i + 1
  } else if (a == "-ms") { monostretch <- as.integer(expect_next()); i <- i + 2
  } else if (a == "-j") { jobs <- max(1L, as.integer(expect_next())); i <- i + 2
  } else if (a == "-samthreads") { samthreads <- max(1L, as.integer(expect_next())); i <- i + 2
  } else if (a == "-seed") { seed <- as.integer(expect_next()); i <- i + 2
  } else if (a == "-pretrain") { pretrain <- TRUE; i <- i + 1
  } else if (a == "-pre_cap") { pre_cap_pairs <- as.integer(expect_next()); i <- i + 2
  } else if (a == "-pre_wmin") { pre_min_weight <- as.integer(expect_next()); i <- i + 2
  } else if (a == "-pre_tieskip") { pre_tie_skip <- TRUE; i <- i + 1
  } else {
    die(paste("Unexpected/incomplete parameter:", a), usage)
  }
}

if (is.null(sample)) die("No Sample name", usage)
message("Sample: ", sample)
if (is.null(snv_file)) die("No VCF(SNV) file", usage)
message("VCF file: ", snv_file)
if (length(bams) == 0) die("No BAM file", usage)
message("BAM files: ", length(bams), " (", paste(bams, collapse = ","), ")")
message("Chromosome: ", chromosome)

if (is.na(tries) || tries <= 0) tries <- 100L
message("Number of tries: ", tries)
if (is.na(min_umis) || min_umis <= 0) min_umis <- 10L
message("Min number of UMIs per SNV: ", min_umis)
if (is.na(min_maf) || min_maf < 0 || min_maf > 0.5) min_maf <- 0
message("Min MAF per SNV: ", min_maf)
if (is.na(discordant_penalty) || discordant_penalty < 1) discordant_penalty <- 5L
message("Discordant penalty: ", discordant_penalty)

if (pre_cap_pairs < 100) pre_cap_pairs <- 5000L
if (pre_min_weight < 1) pre_min_weight <- 3L
if (!is.na(seed)) set.seed(seed)

# ---------------------- Read VCF positions (by chromosome) ------------------
read_snv_positions <- function(vcf_file, chrom) {
  con <- if (grepl("\\.gz$", vcf_file)) {
    pipe(sprintf("gunzip -c %s", shQuote(vcf_file)), "r")
  } else {
    file(vcf_file, "r")
  }
  on.exit(close(con))
  snp_info <- new.env(parent = emptyenv())
  n <- 0L
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (length(lines) == 0) break
    for (ln in lines) {
      if (startsWith(ln, "#")) next
      arr <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      if (length(arr) < 5) next
      if (arr[1] != chrom) next
      pos <- as.integer(arr[2])
      snp_info[[as.character(pos)]] <- list(rs = arr[3], ref = arr[4], alt = arr[5])
      n <- n + 1L
    }
  }
  list(env = snp_info, n = n)
}

vcf <- read_snv_positions(snv_file, chromosome)
snp_info <- vcf$env
message(vcf$n, " SNVs were loaded from VCF file ", snv_file)
if (length(ls(snp_info)) == 0) die("Zero SNVs, run terminated")

# ---------------------- Read BAMs via samtools view ------------------------
# Data structures (list-of-lists of named integer vectors):
# umis[[pos]][["1"]] and umis[[pos]][["2"]] are named integer vectors keyed by "UMI\tCB"
# cb2allele[[cb]][[pos]] is named integer vector keyed by UMI with value 1 or 2
umis <- new.env(parent = emptyenv())
cb2allele <- new.env(parent = emptyenv())
inf_reads <- 0L
inf_alleles <- 0L
hard_umis <- 0L
soft_umis <- 0L

has_mono_stretch <- function(seqstr, k) {
  # e.g., G{15} or A{15} etc.
  grepl(paste0("(G{", k, "}|A{", k, "}|T{", k, "}|C{", k, "})"), seqstr, perl = TRUE)
}

extract_tag <- function(line, tag) { # return first match after <TAB>tag:TYPE:VALUE
  m <- regexpr(paste0("\\t", tag, ":[A-Za-z]:([^\\t]+)"), line, perl = TRUE)
  if (m[1] == -1) return(NA_character_)
  sub("^.*:[A-Za-z]:", "", regmatches(line, m)[[1]])
}

get_cb <- function(line) {
  cb <- extract_tag(line, "CB")
  if (is.na(cb)) cb <- extract_tag(line, "RG")
  if (is.na(cb) || cb == "-") return(NA_character_)
  cb
}

get_ub <- function(line) {
  ub <- extract_tag(line, "UB")
  if (is.na(ub) || ub == "-") return(NA_character_)
  ub
}

ensure_vec <- function(env, key) {
  k <- as.character(key)
  if (!exists(k, envir = env, inherits = FALSE)) assign(k, list(`1` = integer(0), `2` = integer(0)), envir = env)
  get(k, envir = env, inherits = FALSE)
}
set_vec <- function(env, key, val) assign(as.character(key), val, envir = env)

ensure_cb <- function(env, cb) {
  if (!exists(cb, envir = env, inherits = FALSE)) assign(cb, new.env(parent = emptyenv()), envir = env)
  get(cb, envir = env, inherits = FALSE)
}

inc_name <- function(v, nm) {
  if (length(v) == 0L) {
    v <- structure(1L, names = nm); return(v)
  }
  pos <- match(nm, names(v), 0L)
  if (pos == 0L) {
    v <- c(v, structure(1L, names = nm))
  } else {
    v[pos] <- v[pos] + 1L
  }
  v
}

message("Reading BAM files ...")
for (file in bams) {
  message("    Reading ", file, " ...")
  cmd <- if (samthreads > 1) {
    sprintf("samtools view -@ %d -F 256 %s %s", samthreads, shQuote(file), shQuote(chromosome))
  } else {
    sprintf("samtools view -F 256 %s %s", shQuote(file), shQuote(chromosome))
  }
  con <- pipe(cmd, "r")
  repeat {
    lines <- readLines(con, n = 10000L, warn = FALSE)
    if (length(lines) == 0) break
    for (line in lines) {
      if (!grepl("\\tvG:B:i,\\d", line, perl = TRUE)) next
      if (!grepl("\\tvW:i:1", line, perl = TRUE)) next
      parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
      cigar <- parts[6]
      seqstr <- parts[10]
      if (no_softmasking && grepl("\\dS", cigar, perl = TRUE)) next
      if (has_mono_stretch(seqstr, monostretch)) next
      
      vGs <- extract_tag(line, "vG")
      vAs <- extract_tag(line, "vA")
      if (is.na(vGs) || is.na(vAs)) next
      vG <- as.integer(strsplit(vGs, ",", fixed = TRUE)[[1]])
      vA <- suppressWarnings(as.integer(strsplit(vAs, ",", fixed = TRUE)[[1]]))
      if (length(vG) == 0 || length(vA) == 0) next
      
      cb <- get_cb(line)
      if (is.na(cb)) die(paste("Found no CB, no RG tags in BAM file in read:", line))
      ub <- get_ub(line)
      if (!is.na(ub)) {
        if (ub == "-") next
        hard_umis <- hard_umis + 1L
      } else {
        # synthesize UMI: read,flag,chr,pos,cigar,tlen
        # read line: qname flag rname pos mapq cigar rnext pnext tlen
        read <- parts[1]; flag <- parts[2]; chr <- parts[3]; pos <- parts[4]; tlen <- parts[9]
        ub <- paste(cb, chr, pos, flag, cigar, tlen, sep = "_")
        soft_umis <- soft_umis + 1L
      }
      inf_reads <- inf_reads + 1L
      
      for (k in seq_along(vG)) {
        pos0 <- vG[k] + 1L
        pos_key <- as.character(pos0)
        if (!exists(pos_key, envir = snp_info, inherits = FALSE)) next
        alle <- vA[k]
        if (!(alle == 1L || alle == 2L)) next
        if (scramble_alleles) alle <- 1L + sample.int(2L, 1) - 1L
        
        # Update umis
        rec <- ensure_vec(umis, pos_key)
        nm <- paste0(ub, "\t", cb)
        if (alle == 1L) rec[["1"]] <- inc_name(rec[["1"]], nm) else rec[["2"]] <- inc_name(rec[["2"]], nm)
        set_vec(umis, pos_key, rec)
        
        # Update cb2allele
        cb_env <- ensure_cb(cb2allele, cb)
        if (!exists(pos_key, envir = cb_env, inherits = FALSE)) assign(pos_key, integer(0), envir = cb_env)
        v <- get(pos_key, envir = cb_env, inherits = FALSE)
        # store by UMI (without CB) exactly like Perl (per-cell UMI)
        ukey <- ub
        pos_idx <- match(ukey, names(v), 0L)
        if (pos_idx == 0L) {
          v <- c(v, structure(alle, names = ukey))
        } else {
          # if duplicate same UMI at same pos, last write wins (OK for majority later)
          v[pos_idx] <- alle
        }
        assign(pos_key, v, envir = cb_env)
        
        inf_alleles <- inf_alleles + 1L
      }
    }
  }
  close(con)
  message("    Reads with allelic info: ", inf_reads,
          " Alleles: ", inf_alleles,
          " UMIs:", hard_umis, "/", soft_umis, " (hard/soft)")
}

covered_pos <- ls(umis)
n_cells <- length(ls(cb2allele))
message(length(covered_pos), " SNV positions were covered in ", n_cells, " cells were loaded from BAM files")
if (length(covered_pos) == 0) die("Zero SNVs, run terminated")
if (n_cells == 0) die("Zero cells, run terminated")

# ---------------------- Filter SNVs by UMI & MAF ---------------------------
low_umis <- 0L
low_maf <- 0L
ratios <- numeric()

keep_pos <- integer(0)
for (pos_key in covered_pos) {
  rec <- get(pos_key, envir = umis, inherits = FALSE)
  n1 <- length(rec[["1"]]); n2 <- length(rec[["2"]])
  if ((n1 + n2) < min_umis) { low_umis <- low_umis + 1L; next }
  af <- if ((n1 + n2) > 0) n1 / (n1 + n2) else 0
  mn <- min(n1, n2); mx <- max(n1, n2)
  if ((mn + mx) > 0) ratios <- c(ratios, 100 * mn / (mn + mx))
  if (af < min_maf || af > (1 - min_maf)) { low_maf <- low_maf + 1L; next }
  keep_pos <- c(keep_pos, as.integer(pos_key))
}
if (length(keep_pos) == 0) die("No SNVs to phase")
if (low_umis > 0) message("Excluded ", low_umis, " SNVs having less than ", min_umis, " reads/UMIs")
if (low_maf > 0) message("Excluded ", low_maf, " SNVs having less than ", min_maf, " minor allele frequency")
ratios <- sort(ratios)
median_maf <- if (length(ratios) == 0) 0 else {
  if (length(ratios) %% 2 == 1) ratios[(length(ratios) + 1)/2] else mean(ratios[(length(ratios)/2):(length(ratios)/2 + 1)])
}
message(sprintf("Median MAF: %0.2f", median_maf))

phased_pos <- sort(keep_pos)
# initial directions: {-1, 0, 1}; pretraining may overwrite some
phased_dir <- sample.int(3L, length(phased_pos), replace = TRUE) - 2L

# ---------------------- Pretraining (unsupervised seeding) -----------------
cell_majority_calls <- function(cb2allele_env, phased_positions, tie_skip = TRUE) {
  is_pos <- new.env(parent = emptyenv())
  for (p in phased_positions) assign(as.character(p), TRUE, envir = is_pos)
  out <- new.env(parent = emptyenv()) # cb -> env(pos -> allele)
  cbs <- ls(cb2allele_env)
  for (cb in cbs) {
    cb_env <- get(cb, envir = cb2allele_env, inherits = FALSE)
    pos_keys <- ls(cb_env)
    for (pk in pos_keys) {
      if (!exists(pk, envir = is_pos, inherits = FALSE)) next
      v <- get(pk, envir = cb_env, inherits = FALSE) # named int (umi -> 1/2)
      if (length(v) == 0) next
      c1 <- sum(v == 1L, na.rm = TRUE)
      c2 <- sum(v == 2L, na.rm = TRUE)
      if ((c1 + c2) == 0) next
      if (tie_skip && c1 == c2) next
      alle <- if (c1 >= c2) 1L else 2L
      if (!exists(cb, envir = out, inherits = FALSE)) assign(cb, new.env(parent = emptyenv()), envir = out)
      out_cb <- get(cb, envir = out, inherits = FALSE)
      assign(pk, structure(alle, names = NULL), envir = out_cb)
    }
  }
  out
}

add_vote <- function(W, deg, pi, pj, same) {
  s <- if (same) 1L else -1L
  if (is.null(W[[pi]])) W[[pi]] <- integer()
  if (is.null(W[[pj]])) W[[pj]] <- integer()
  W[[pi]][as.character(pj)] <- (W[[pi]][as.character(pj)] %||% 0L) + s
  W[[pj]][as.character(pi)] <- (W[[pj]][as.character(pi)] %||% 0L) + s
  deg[[as.character(pi)]] <- (deg[[as.character(pi)]] %||% 0L) + 1L
  deg[[as.character(pj)]] <- (deg[[as.character(pj)]] %||% 0L) + 1L
  list(W = W, deg = deg)
}

build_graph_from_cells <- function(cell_calls_env, pre_cap) {
  W <- list(); deg <- list()
  cbs <- ls(cell_calls_env)
  for (cb in cbs) {
    pos_env <- get(cb, envir = cell_calls_env, inherits = FALSE)
    pos_keys <- ls(pos_env)
    L <- length(pos_keys)
    if (L < 2) next
    all_pairs <- L * (L - 1) / 2
    if (all_pairs <= pre_cap) {
      for (i in 1:(L - 1)) for (j in (i + 1):L) {
        pi <- as.integer(pos_keys[i]); pj <- as.integer(pos_keys[j])
        ai <- get(pos_keys[i], envir = pos_env); aj <- get(pos_keys[j], envir = pos_env)
        tmp <- add_vote(W, deg, pi, pj, ai == aj); W <- tmp$W; deg <- tmp$deg
      }
    } else {
      seen <- new.env(parent = emptyenv())
      need <- pre_cap; tries <- 0L
      while (length(ls(seen)) < need && tries < 5L * need) {
        tries <- tries + 1L
        i <- sample.int(L, 1); j <- sample.int(L, 1); if (i == j) next
        ii <- min(i, j); jj <- max(i, j)
        key <- paste0(ii, "#", jj)
        if (exists(key, envir = seen, inherits = FALSE)) next
        assign(key, TRUE, envir = seen)
        pi <- as.integer(pos_keys[ii]); pj <- as.integer(pos_keys[jj])
        ai <- get(pos_keys[ii], envir = pos_env); aj <- get(pos_keys[jj], envir = pos_env)
        tmp <- add_vote(W, deg, pi, pj, ai == aj); W <- tmp$W; deg <- tmp$deg
      }
    }
  }
  list(W = W, deg = deg)
}

bfs_seed_from_graph <- function(W, deg, phased_positions, wmin) {
  dir <- new.env(parent = emptyenv())
  cands <- as.integer(phased_positions)
  # order by degree desc
  degv <- sapply(as.character(cands), function(k) deg[[as.character(k)]] %||% 0L)
  ord <- order(-degv)
  cands <- cands[ord]
  for (start in cands) {
    k <- as.character(start)
    if (exists(k, envir = dir, inherits = FALSE)) next
    has_neighbors <- !is.null(W[[start]]) && length(W[[start]]) > 0
    if (!has_neighbors) next
    assign(k, 1L, envir = dir)
    q <- list(start)
    while (length(q) > 0) {
      u <- q[[1]]; q <- q[-1]
      Wu <- W[[u]]
      if (is.null(Wu) || length(Wu) == 0) next
      for (vkey in names(Wu)) {
        w <- Wu[[vkey]]
        if (abs(w) < wmin) next
        want <- if (w >= 0) (get(as.character(u), envir = dir)) else -(get(as.character(u), envir = dir))
        if (!exists(vkey, envir = dir, inherits = FALSE)) {
          assign(vkey, as.integer(want), envir = dir)
          q[[length(q) + 1]] <- as.integer(vkey)
        }
      }
    }
  }
  sapply(as.character(phased_positions), function(p) get(as.character(p), envir = dir, inherits = FALSE) %||% 0L, USE.NAMES = FALSE)
}

run_pretraining_seed <- function(cb2allele_env, phased_positions, pre_cap_pairs, pre_min_weight, pre_tie_skip) {
  cc <- cell_majority_calls(cb2allele_env, phased_positions, tie_skip = pre_tie_skip)
  GG <- build_graph_from_cells(cc, pre_cap_pairs)
  seed_dir <- bfs_seed_from_graph(GG$W, GG$deg, phased_positions, pre_min_weight)
  # stats
  edges <- 0L; strong <- 0L
  for (i in names(GG$W)) {
    edges <- edges + length(GG$W[[i]])
    strong <- strong + sum(abs(unlist(GG$W[[i]])) >= pre_min_weight)
  }
  assigned <- sum(seed_dir != 0)
  message("Pretrain: nodes=", length(phased_positions),
          " edges~", edges, " strong-edges~", strong,
          " assigned=", assigned, " (", sprintf("%0.1f", 100 * assigned / length(phased_positions)), "%)")
  seed_dir
}

if (pretrain) {
  message("Running unsupervised pretraining to seed directions ...")
  pre_seed <- run_pretraining_seed(cb2allele, phased_pos, pre_cap_pairs, pre_min_weight, pre_tie_skip)
  phased_dir <- pre_seed
}

# ---------------------- Core scoring helpers -------------------------------
build_cb2phase <- function(phased_pos, phased_dir, umis_env) {
  cb2phase <- new.env(parent = emptyenv()) # cb -> c(c1,c2)
  for (idx in seq_along(phased_pos)) {
    dir <- phased_dir[idx]
    if (dir == 0L) next
    pos_key <- as.character(phased_pos[idx])
    rec <- get(pos_key, envir = umis_env, inherits = FALSE)
    for (alle in c("1", "2")) {
      if (length(rec[[alle]]) == 0) next
      nms <- names(rec[[alle]])
      # bc key: "UMI\tCB"
      cbs <- sub("^.*\t", "", nms, perl = TRUE)
      for (cb in cbs) {
        if (!exists(cb, envir = cb2phase, inherits = FALSE)) assign(cb, c(0L, 0L), envir = cb2phase)
        v <- get(cb, envir = cb2phase, inherits = FALSE)
        if (dir == 1L) {
          # allele stays
          if (alle == "1") v[1] <- v[1] + 1L else v[2] <- v[2] + 1L
        } else if (dir == -1L) {
          # swap
          if (alle == "1") v[2] <- v[2] + 1L else v[1] <- v[1] + 1L
        }
        assign(cb, v, envir = cb2phase)
      }
    }
  }
  cb2phase
}

score_from_cb2phase <- function(cb2phase_env, penalty) {
  cbs <- ls(cb2phase_env)
  total <- 0L; disc <- 0L; conc <- 0L
  for (cb in cbs) {
    v <- get(cb, envir = cb2phase_env, inherits = FALSE)
    mn <- min(v[1], v[2]); mx <- max(v[1], v[2])
    total <- total + mn + mx
    if (mx == 0L) next
    disc <- disc + mn
    conc <- conc + (mx - 1L)
  }
  list(score = conc - penalty * disc, total = total, disc = disc, conc = conc)
}

# Apply delta for flipping one SNV to new state, update cb2phase in place
apply_flip <- function(cb2phase_env, umis_env, pos, old_dir, new_dir) {
  if (old_dir == new_dir) return(invisible(NULL))
  pos_key <- as.character(pos)
  rec <- get(pos_key, envir = umis_env, inherits = FALSE)
  adj <- function(cb, alle, s) {
    if (!exists(cb, envir = cb2phase_env, inherits = FALSE)) assign(cb, c(0L, 0L), envir = cb2phase_env)
    v <- get(cb, envir = cb2phase_env, inherits = FALSE)
    if (alle == 1L) v[1] <- v[1] + s else v[2] <- v[2] + s
    assign(cb, v, envir = cb2phase_env)
  }
  # Remove old_dir contributions
  if (old_dir != 0L) {
    for (alle in c(1L, 2L)) {
      vec <- rec[[as.character(alle)]]
      if (length(vec) == 0) next
      cbs <- sub("^.*\t", "", names(vec), perl = TRUE)
      if (old_dir == 1L) {
        # direct
        for (cb in cbs) adj(cb, alle, -1L)
      } else if (old_dir == -1L) {
        # swapped
        for (cb in cbs) adj(cb, 3L - alle, -1L)
      }
    }
  }
  # Add new_dir contributions
  if (new_dir != 0L) {
    for (alle in c(1L, 2L)) {
      vec <- rec[[as.character(alle)]]
      if (length(vec) == 0) next
      cbs <- sub("^.*\t", "", names(vec), perl = TRUE)
      if (new_dir == 1L) {
        for (cb in cbs) adj(cb, alle, +1L)
      } else if (new_dir == -1L) {
        for (cb in cbs) adj(cb, 3L - alle, +1L)
      }
    }
  }
  invisible(NULL)
}

# ---------------------- Single try (greedy improvement) --------------------
run_single_try <- function(try_id, phased_pos, phased_dir, umis_env, penalty, pretrain, seed) {
  if (!is.na(seed)) set.seed(seed + try_id * 1337L)
  # Shuffle order; keep pretrain directions as-is; otherwise randomize dirs
  ord <- sample.int(length(phased_pos))
  pos <- phased_pos[ord]
  dir <- phased_dir[ord]
  if (!pretrain) dir <- sample(-1L:1L, length(dir), replace = TRUE)
  
  message("Try #", try_id, ", initializing order and XCI status ...")
  cb2phase <- build_cb2phase(pos, dir, umis_env)
  sc <- score_from_cb2phase(cb2phase, penalty)
  best_score <- sc$score
  message("    SNVs: ", length(pos), " Initial score: ", best_score,
          ", discordance rate : ",
          if ((sc$disc + sc$conc) > 0) sprintf("%0.2f", 100 * sc$disc / (sc$disc + sc$conc)) else "0.00")
  
  pass <- 0L; imps <- 1L
  while (imps > 0L) {
    imps <- 0L; pass <- pass + 1L
    last_imp <- -1L
    for (e in seq_along(dir)) {
      old <- dir[e]
      # compute baseline score (sanity)
      base <- score_from_cb2phase(cb2phase, penalty)$score
      if (base != best_score) stop("score_before ne best_score")
      for (new in c(-1L, 0L, 1L)) {
        if (new == old) next
        apply_flip(cb2phase, umis_env, pos[e], old, new)
        new_sc <- score_from_cb2phase(cb2phase, penalty)$score
        if (new_sc > best_score || (new_sc == best_score && new == 0L)) {
          imps <- imps + 1L
          last_imp <- e
          dir[e] <- new
          best_score <- new_sc
          # keep updated cb2phase
        } else {
          # restore
          apply_flip(cb2phase, umis_env, pos[e], new, old)
          chk <- score_from_cb2phase(cb2phase, penalty)$score
          if (chk != base) stop("score_restored ne score_before")
        }
      }
    }
    message("    Pass: ", pass, " Improvements: ", imps, " Score: ", best_score)
  }
  list(score = best_score, pos = pos, dir = dir)
}

# ---------------------- Parallel tries -------------------------------------
library(parallel)
res <- mclapply(seq_len(tries), function(tt)
  run_single_try(tt, phased_pos, phased_dir, umis, discordant_penalty, pretrain, seed),
  mc.cores = jobs
)

best_idx <- which.max(vapply(res, `[[`, numeric(1), "score"))
global_best <- res[[best_idx]]$score
# realign best dir to parent order (phased_pos)
pos2dir <- structure(res[[best_idx]]$dir, names = as.character(res[[best_idx]]$pos))
global_best_dir <- vapply(as.character(phased_pos), function(p) pos2dir[[p]] %||% 0L, integer(1))

# Optionally compare with previous run summary (-i)
if (improve_existing && file.exists(paste0(sample, "_chr", chromosome, "_XCISE_summary.txt"))) {
  first_line <- readLines(paste0(sample, "_chr", chromosome, "_XCISE_summary.txt"), n = 1L)
  m <- regexpr("^Best\\s+score\\s+:\\s+(\\d+)", first_line, perl = TRUE)
  if (m[1] != -1) {
    prev <- as.numeric(sub(".*:\\s+", "", regmatches(first_line, m)[[1]]))
    if (!is.na(prev) && prev > global_best) {
      global_best <- prev
      message("    Updated global best score from previous run: ", global_best)
      # (Like Perl, we don't recover previous best directions—score only.)
    }
  }
}

# ---------------------- Outputs (VCF, bc2xci, summary) ---------------------
# rebuild cb2phase with best directions
cb2phase_best <- build_cb2phase(phased_pos, global_best_dir, umis)

write_summary <- function(path, best_score, cb2phase_env, global_best_dir, median_maf) {
  cbs <- ls(cb2phase_env)
  total <- 0L; disc <- 0L; conc <- 0L
  for (cb in cbs) {
    v <- get(cb, envir = cb2phase_env, inherits = FALSE)
    mn <- min(v[1], v[2]); mx <- max(v[1], v[2])
    total <- total + mn + mx
    if (mx == 0L) next
    disc <- disc + mn
    conc <- conc + (mx - 1L)
  }
  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  writeLines(sprintf("Best score  : %d", best_score), con)
  writeLines(sprintf("Total UMIs  : %d", total), con)
  writeLines(sprintf("Concordant  : %d", conc), con)
  writeLines(sprintf("Discordant  : %d", disc), con)
  if ((conc + disc) > 0)
    writeLines(sprintf("Discordancy : %f", disc / (conc + disc)), con)
  else
    writeLines("Discordancy : N/A", con)
  writeLines(sprintf("Total SNVs  : %d", length(global_best_dir)), con)
  writeLines(sprintf("Median MAF  : %3.2f %%", median_maf), con)
  # counts by direction
  ref1 <- sum(global_best_dir == 1L)
  alt1 <- sum(global_best_dir == -1L)
  noninf <- sum(global_best_dir == 0L)
  writeLines(sprintf("SNV Ref/Alt : %d", ref1), con)
  writeLines(sprintf("SNV Alt/Ref : %d", alt1), con)
  writeLines(sprintf("XCI-inf SNVs: %d", ref1 + alt1), con)
  writeLines(sprintf("Non-inf SNVs: %d", noninf), con)
  invisible(NULL)
}

message("    Outputting VCF...")
vcf_out <- file(paste0(sample, "_chr", chromosome, "_XCISE.vcf"), "w"); on.exit(close(vcf_out), add = TRUE)
writeLines("##fileformat=VCFv4.2", vcf_out)
writeLines(paste("#CHROM","POS","ID","REF","ALT","QUAL","INFO", sep = "\t"), vcf_out)

ord_idx <- order(phased_pos)
for (e in ord_idx) {
  pos <- phased_pos[e]
  pos_key <- as.character(pos)
  rec <- get(pos_key, envir = umis, inherits = FALSE)
  n1 <- length(rec[["1"]]); n2 <- length(rec[["2"]])
  info <- get(pos_key, envir = snp_info, inherits = FALSE)
  x1alle <- if (global_best_dir[e] == 0L) "Unk" else if (global_best_dir[e] == 1L) "Ref" else if (global_best_dir[e] == -1L) "Alt" else "Err"
  line <- paste(chromosome, pos, info$rs, info$ref, info$alt,
                n1 + n2, "PASS",
                paste0("X1A=", x1alle, ";AD=", n1, ",", n2),
                sep = "\t")
  writeLines(line, vcf_out)
}

message("    Outputting barcode to phase data...")
bc_out <- file(paste0(sample, "_chr", chromosome, "_XCISE_bc2xci.txt"), "w"); on.exit(close(bc_out), add = TRUE)
total0 <- total1 <- total2 <- total3 <- total4 <- 0L

cbs <- sort(ls(cb2allele))
for (cb in cbs) {
  cb_env <- get(cb, envir = cb2allele, inherits = FALSE)
  hap1 <- 0L; hap2 <- 0L
  for (e in seq_along(global_best_dir)) {
    d <- global_best_dir[e]; if (d == 0L) next
    pos <- phased_pos[e]; pk <- as.character(pos)
    if (!exists(pk, envir = cb_env, inherits = FALSE)) next
    v <- get(pk, envir = cb_env, inherits = FALSE) # named int (umi -> 1/2)
    # count per-UMI contributions like Perl (each UMI once)
    if (d == 1L) {
      hap1 <- hap1 + sum(v == 1L) + sum(v == 2L & FALSE) # explicit
      hap2 <- hap2 + sum(v == 2L & TRUE) * 0L             # noop
      hap2 <- hap2 + sum(v == 2L & FALSE)
      # effectively: hap1 += (v==1); hap2 += (v==2) handled below:
      hap2 <- hap2 + sum(v == 2L) * 0L
    } else if (d == -1L) {
      # swap
      # hap1 += (v==2), hap2 += (v==1)
      # do them explicitly:
      hap1 <- hap1 + sum(v == 2L)
      hap2 <- hap2 + sum(v == 1L)
    }
    if (d == 1L) {
      hap1 <- hap1 + 0L + sum(v == 1L)
      hap2 <- hap2 + 0L + sum(v == 2L) * 0L
    }
  }
  # The above blocks ensure alignment with Perl’s two cases;
  # simplify to exact Perl logic:
  hap1 <- 0L; hap2 <- 0L
  for (e in seq_along(global_best_dir)) {
    d <- global_best_dir[e]; if (d == 0L) next
    pos <- phased_pos[e]; pk <- as.character(pos)
    if (!exists(pk, envir = cb_env, inherits = FALSE)) next
    v <- get(pk, envir = cb_env, inherits = FALSE)
    if (d == 1L) {
      hap1 <- hap1 + sum(v == 1L)
      hap2 <- hap2 + sum(v == 2L)
    } else if (d == -1L) {
      hap1 <- hap1 + sum(v == 2L)
      hap2 <- hap2 + sum(v == 1L)
    }
  }
  lab <- "?"
  if (hap1 == 0L && hap2 == 0L) { total0 <- total0 + 1L; lab <- "Unknown"
  } else if (hap1 >= 2L && hap1 / (hap1 + hap2) >= 0.9) { total1 <- total1 + 1L; lab <- "X1"
  } else if (hap2 >= 2L && hap2 / (hap1 + hap2) >= 0.9) { total2 <- total2 + 1L; lab <- "X2"
  } else if (hap1 > 0L && hap2 > 0L) { total3 <- total3 + 1L; lab <- "Both"
  } else { total4 <- total4 + 1L; lab <- "Low_coverage" }
  writeLines(paste(cb, hap1, hap2, lab, sep = "\t"), bc_out)
}

grand_total <- total0 + total1 + total2 + total3 + total4

# Summary (first block)
summary_path <- paste0(sample, "_chr", chromosome, "_XCISE_summary.txt")
write_summary(summary_path, global_best, cb2phase_best, global_best_dir, median_maf)

# Append category counts
append_cat <- function(path, total1, total2, total3, total4, total0, grand) {
  con <- file(path, "a"); on.exit(close(con), add = TRUE)
  w <- function(fmt, ...) writeLines(sprintf(fmt, ...), con)
  w("X1 cells    : %5d( %3.2f %% )", total1, if (grand > 0) 100 * total1 / grand else 0)
  w("X2 cells    : %5d( %3.2f %% )", total2, if (grand > 0) 100 * total2 / grand else 0)
  w("Both X      : %5d( %3.2f %% )", total3, if (grand > 0) 100 * total3 / grand else 0)
  w("LowCoverage : %5d( %3.2f %% )", total4, if (grand > 0) 100 * total4 / grand else 0)
  w("Unknown     : %5d( %3.2f %% )", total0, if (grand > 0) 100 * total0 / grand else 0)
}

append_cat(summary_path, total1, total2, total3, total4, total0, grand_total)

message("    X1/X2/Both/LowC/Unknown cells: ",
        paste(total1, total2, total3, total4, total0, sep = " / "))
if (total1 > total2) message("    Switching X1 and X2, so that there are more X2 cells than X1 cells ...")

# Flip all directions until X1 <= X2 (like Perl: do { ... } until total1 <= total2)
while (total1 > total2) {
  global_best_dir <- -global_best_dir
  # recompute bc categories quickly
  total0 <- total1 <- total2 <- total3 <- total4 <- 0L
  # (we do NOT re-write VCF or bc2xci again; Perl only appends summary counts)
  # but we need to append updated summary block:
  # rebuild cb2phase_best for reporting
  cb2phase_best <- build_cb2phase(phased_pos, global_best_dir, umis)
  write_summary(summary_path, global_best, cb2phase_best, global_best_dir, median_maf)
  
  # recount cats for logging/appending
  cbs <- sort(ls(cb2allele))
  for (cb in cbs) {
    cb_env <- get(cb, envir = cb2allele, inherits = FALSE)
    hap1 <- 0L; hap2 <- 0L
    for (e in seq_along(global_best_dir)) {
      d <- global_best_dir[e]; if (d == 0L) next
      pos <- phased_pos[e]; pk <- as.character(pos)
      if (!exists(pk, envir = cb_env, inherits = FALSE)) next
      v <- get(pk, envir = cb_env, inherits = FALSE)
      if (d == 1L) { hap1 <- hap1 + sum(v == 1L); hap2 <- hap2 + sum(v == 2L) }
      else if (d == -1L) { hap1 <- hap1 + sum(v == 2L); hap2 <- hap2 + sum(v == 1L) }
    }
    if (hap1 == 0L && hap2 == 0L) total0 <- total0 + 1L
    else if (hap1 >= 2L && hap1 / (hap1 + hap2) >= 0.9) total1 <- total1 + 1L
    else if (hap2 >= 2L && hap2 / (hap1 + hap2) >= 0.9) total2 <- total2 + 1L
    else if (hap1 > 0L && hap2 > 0L) total3 <- total3 + 1L
    else total4 <- total4 + 1L
  }
  grand_total <- total0 + total1 + total2 + total3 + total4
  append_cat(summary_path, total1, total2, total3, total4, total0, grand_total)
}

message("Done.")
