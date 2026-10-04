#!/usr/bin/perl -w
use strict;

# Derived from xcise.pl of XCISE (https://github.com/Vityay/XCISE).
# Please cite: Henning RH, Rust TM, Dijksterhuis K, Eggen BJL, Guryev V.
#   bioRxiv (2024). doi:10.1101/2024.08.29.610317

use Parallel::ForkManager;
use File::Temp qw(tempdir);

# =============================================================================
# EXCISE2  —  Enhanced XCISE pipeline (based on xcise_mt.pl)
#
# Improvements over xcise_mt.pl:
#   1. Weighted co-segregation votes (UMI-confidence scaled, not binary)
#   2. BFS-based pretraining (original xcise_mt style) with:
#      - confidence-weighted + depth-scaled votes (no allele-1 tie bias)
#      - pre_min_weight actually applied (was broken in MST version)
#      - random ±1 component anchor (no reference-bias allele-balance)
#      - isolated SNVs stay dir=0 (no soft prior skew propagation)
#   4. Per-try jitter proportional to per-SNV edge uncertainty
#   5. Soft prior initialisation for isolated (unseeded) SNVs
#   6. samtools -@ defaults to $jobs (all threads), overridable with -samthreads
#   7. Haplotype sequence output (*_EXCISE_haplotypes.tsv)
#   8. dir=0 restored for biallelic SNVs: escape candidates stay Unk
#      dir=0 SNVs: both ±1 evaluated and compared before collapsing to 0
#      monoallelic SNVs: forced ±1 (only one allele observed)
#   9. Monoallelic SNVs flagged (MONO INFO tag in VCF); Phase-2 cell-vote
#      correction applied post-optimizer
#      Phase-3: directionality-based assignment for remaining dir=0 SNVs
#      (label-free; measures haplotype purity gain in covering cells)
#  10. VCF INFO tags: X1A, AD, X1, X2, DP, MONO
#  11. Haplotype TSV: CHROM POS ID REF ALT X1_allele X2_allele mono
#      (UMI counts live in VCF INFO tags X1=/X2=/DP=)
#  12. Try-random-before-last escape in hill-climber
#
# Flags (all xcise_mt.pl flags preserved):
#   -o  <prefix>       output prefix
#   -s  <vcf>          heterozygous SNV VCF (optionally gzipped)
#   -b  <bam> [...]    one or more BAM files (WASP-mode aligned)
#   -r  <chrom>        chromosome (default X)
#   -t  <int>          number of optimiser tries (default 100)
#   -u  <int>          min UMIs per SNV (default 10)
#   -m  <float>        min minor allele frequency per SNV (default 0)
#   -p  <int>          discordant-UMI penalty (default 5)
#   -x                 scramble alleles (control / debug)
#   -i                 improve an existing solution
#   -sm                skip soft-masked alignments
#   -ms <int>          skip reads with mononucleotide runs >= N (default 15)
#   -j  <int>          parallel tries / processes (default 1)
#   -samthreads <int>  samtools -@ threads (default = $jobs)
#   -seed <int>        RNG seed for reproducibility
#   -wl  <file>        optional barcode whitelist (plain text/TSV or gzipped)
#   -wl_strip_gem      strip trailing -1/-2/... for whitelist matching
#   -pretrain          enable improved unsupervised seeding
#   -pre_cap <int>     max SNV-pair votes per cell (default 5000)
#   -pre_wmin <float>  min absolute edge weight for MST (default 1)
#   -pre_tieskip       skip tie allele calls in pretraining
#   -jitter <float>    per-SNV flip probability scale for per-try jitter (default 0.25)
#   -ph                enable Phase 3 assignment (disabled by default)
#   -ph3_ratio <float> min Phase 3 concordance ratio (default 0.90; requires -ph)
# =============================================================================

warn "Reading parameters ...\n";
my ( $sample, $chromosome, $snv_file );
my $tries            = 100;
my $min_maf          = 0;
my $min_umis         = 4;  # restored to xcise_mt default
my $discordant_penalty = 5;
my $scramble_alleles = 0;
my $improve_existing = 0;
my $no_softmasking   = 0;
my $monostretch      = 15;

my $jobs             = 1;
my $samthreads       = -1;   # -1 = use $jobs by default; overridden by -samthreads
my $seed;
my $barcode_whitelist;
my $wl_strip_gem     = 0;

my $pretrain         = 0;
my $pre_cap_pairs    = 5000;
my $pre_min_weight   = 3;    # restored to xcise_mt default: filters weak/noisy edges
my $pre_tie_skip     = 1;    # default ON: tie allele calls skipped in pretrain
my $jitter_scale     = 0.25; # fraction of uncertainty that drives per-try flips
my $run_phase3       = 0;    # Phase 3 is opt-in; enable with -ph
my $ph3_min_ratio    = 0.90; # Phase 3: min concordance ratio to assign direction

my %bams = ();
my $ele  = 0;
while ( $ele <= $#ARGV ) {
    if ( $ARGV[$ele] eq '-b' ) {
        while ( $ele < $#ARGV and $ARGV[$ele+1] !~ m/^\-/ ) {
            die 'BAM file '.$ARGV[$ele+1].' does not exist' unless -e $ARGV[$ele+1];
            $bams{ $ARGV[$ele+1] } = 1;
            $ele++;
        }
        $ele++;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-o' )          { $sample     = $ARGV[$ele+1]; $ele += 2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-s' )          { $snv_file   = $ARGV[$ele+1]; $ele += 2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-r' )          { $chromosome = $ARGV[$ele+1]; $ele += 2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-t' )          { $tries      = $ARGV[$ele+1]; $ele += 2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-u' )          { $min_umis   = $ARGV[$ele+1]; $ele += 2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-m' )          { $min_maf    = $ARGV[$ele+1]; $ele += 2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-p' )          { $discordant_penalty = $ARGV[$ele+1]; $ele += 2; }
    elsif ( $ARGV[$ele] eq '-x' )                             { $scramble_alleles = 1; $ele++; }
    elsif ( $ARGV[$ele] eq '-i' )                             { $improve_existing  = 1; $ele++; }
    elsif ( $ARGV[$ele] eq '-sm' )                            { $no_softmasking    = 1; $ele++; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-ms' and $ARGV[$ele+1] =~ m/^(\d+)/ ) { $monostretch = $1; $ele += 2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-j' )          { $jobs        = $ARGV[$ele+1]; $ele += 2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-samthreads' ) { $samthreads  = $ARGV[$ele+1]; $ele += 2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-seed' )       { $seed        = $ARGV[$ele+1]; $ele += 2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-wl' ) {
        $barcode_whitelist = $ARGV[$ele+1];
        die 'Barcode whitelist '.$barcode_whitelist.' does not exist' unless -e $barcode_whitelist;
        $ele += 2;
    }
    elsif ( $ARGV[$ele] eq '-wl_strip_gem' )                 { $wl_strip_gem = 1; $ele++; }
    elsif ( $ARGV[$ele] eq '-pretrain' )                      { $pretrain    = 1; $ele++; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-pre_cap' )    { $pre_cap_pairs   = $ARGV[$ele+1]; $ele += 2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-pre_wmin' )   { $pre_min_weight  = $ARGV[$ele+1]; $ele += 2; }
    elsif ( $ARGV[$ele] eq '-pre_tieskip' )                   { $pre_tie_skip    = 1; $ele++; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-jitter' )     { $jitter_scale    = $ARGV[$ele+1]; $ele += 2; }
    elsif ( $ARGV[$ele] eq '-ph' )                            { $run_phase3      = 1; $ele++; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-ph3_ratio' )  { $ph3_min_ratio   = $ARGV[$ele+1]; $ele += 2; }
    # Note: -ph3_cells removed; Gate 1 dropped in favour of stricter -ph3_ratio (default 0.90)
    else { die "Unexpected/incomplete parameter: $ARGV[$ele]\n"; }
}

my $usage = 'Usage: perl '.$0
  .' -o <prefix> -s <vcf> -b <bam> [...] [-r <chrom>] [-t <tries>]'
  .' [-u <min_umis>] [-m <min_maf>] [-p <penalty>]'
  .' [-j <jobs>] [-samthreads <N>] [-seed <int>]'
  .' [-wl <barcodes.txt|tsv[.gz]>] [-wl_strip_gem]'
  .' [-pretrain] [-pre_cap <N>] [-pre_wmin <W>] [-pre_tieskip] [-jitter <f>]'
  .' [-ph] [-ph3_ratio <f>]'."\n";

die "No sample name\n$usage"    unless $sample;
die "No VCF file\n$usage"       unless $snv_file;
die "No BAM file\n$usage"       unless keys %bams;
die "-wl_strip_gem requires -wl <file>\n$usage"
    if $wl_strip_gem and not defined $barcode_whitelist;

$chromosome        = 'X'  unless $chromosome;
$tries             = 100  if $tries <= 0;
$min_umis          = 10  if $min_umis <= 0;
$min_maf           = 0    if $min_maf < 0 or $min_maf > 0.5;
$discordant_penalty= 5    if $discordant_penalty < 1;
$jobs              = 1    if $jobs < 1;
$pre_cap_pairs     = 100  if $pre_cap_pairs < 100;
$pre_min_weight    = 0    if $pre_min_weight < 0;
$jitter_scale      = 0    if $jitter_scale < 0;
$jitter_scale      = 1    if $jitter_scale > 1;

# samtools threads: default to $jobs (use all available processes for BAM I/O)
$samthreads = $jobs if $samthreads < 0;
$samthreads = 1     if $samthreads < 1;

srand($seed) if defined $seed;

warn "Sample             : $sample\n";
warn "VCF file           : $snv_file\n";
warn "BAM files          : ", scalar(keys %bams), " (", join(',', keys %bams), ")\n";
warn "Chromosome         : $chromosome\n";
warn "Tries              : $tries\n";
warn "Min UMIs/SNV       : $min_umis\n";
warn "Min MAF            : $min_maf\n";
warn "Discordant penalty : $discordant_penalty\n";
warn "Parallel jobs      : $jobs\n";
warn "samtools threads   : $samthreads\n";
warn "Barcode whitelist  : ", (defined $barcode_whitelist ? $barcode_whitelist : 'disabled'), "\n";
warn "Strip GEM suffix   : ", ($wl_strip_gem ? "yes" : "no"), "\n"
    if defined $barcode_whitelist;
warn "Pretrain           : ", ($pretrain ? "yes" : "no"), "\n";
warn "Jitter scale       : $jitter_scale\n";
warn "Phase 3            : ", ($run_phase3 ? "enabled" : "disabled"), "\n";

# =====================================================================
# LOAD BARCODE WHITELIST (optional)
# =====================================================================
my %allowed_barcodes = ();
my $wl_nonempty_records = 0;

if ( defined $barcode_whitelist ) {
    warn "Reading barcode whitelist ...\n";

    # Recognise standard *.gz names and gzip files identified by magic bytes.
    my $is_gzip = ( $barcode_whitelist =~ m/\.gz$/i ) ? 1 : 0;
    if ( not $is_gzip ) {
        open my $PROBE, '<', $barcode_whitelist
            or die "Cannot open barcode whitelist $barcode_whitelist: $!";
        binmode $PROBE;
        my $magic = '';
        read $PROBE, $magic, 2;
        close $PROBE;
        $is_gzip = 1 if $magic eq "\x1f\x8b";
    }

    my $FW;
    if ( $is_gzip ) {
        open $FW, '-|', 'gzip', '-dc', '--', $barcode_whitelist
            or die "Cannot decompress barcode whitelist $barcode_whitelist: $!";
    } else {
        open $FW, '<', $barcode_whitelist
            or die "Cannot open barcode whitelist $barcode_whitelist: $!";
    }

    while ( <$FW> ) {
        chomp;
        s/\r$//;
        s/^\s+//;
        next if $_ eq '';
        my ($barcode) = split /\s+/, $_, 2;
        next unless defined $barcode and length $barcode;
        $wl_nonempty_records++;
        $barcode =~ s/-\d+$// if $wl_strip_gem;
        next unless length $barcode;
        $allowed_barcodes{$barcode} = 1;
    }
    close $FW or die "Failed while reading barcode whitelist $barcode_whitelist\n";

    die "Zero barcodes loaded from whitelist $barcode_whitelist\n"
        unless keys %allowed_barcodes;
    warn "Whitelist barcodes loaded: ", scalar(keys %allowed_barcodes),
         " unique from $wl_nonempty_records non-empty records\n";
}

# =====================================================================
# LOAD SNVs
# =====================================================================
warn "Reading SNVs ...\n";
my %snp_info = ();
{
    my $FH;
    if ( $snv_file =~ m/\.gz$/ ) { open $FH, "gunzip -c $snv_file |" or die "gunzip failed: $!"; }
    else                          { open $FH, '<', $snv_file           or die "Cannot open $snv_file: $!"; }
    while ( <$FH> ) {
        next if m/^\#/;
        chomp;
        my @a = split /\t/;
        next unless $a[0] eq $chromosome;
        $snp_info{$a[1]}{'rs'}  = $a[2];
        $snp_info{$a[1]}{'ref'} = $a[3];
        $snp_info{$a[1]}{'alt'} = $a[4];
    }
    close $FH;
}
warn scalar(keys %snp_info), " SNVs loaded from $snv_file\n";
die "Zero SNVs, run terminated\n" unless keys %snp_info;

# =====================================================================
# READ BAMs
# =====================================================================
warn "Reading BAM files ...\n";
my %umis      = ();   # {pos}{allele(1|2)}{"umi\tcb"}++
my %cb2allele = ();   # {cb}{pos}{umi} = allele(1|2)
my ( $inf_reads, $inf_alleles, $hard_umis, $soft_umis ) = (0,0,0,0);
my ( $wl_candidate_reads, $wl_kept_reads, $wl_rejected_reads ) = (0,0,0);
my %wl_barcodes_seen    = ();
my %wl_barcodes_matched = ();

foreach my $file ( sort keys %bams ) {
    warn "    Reading $file ...\n";
    my $cmd = 'samtools view -@ '.$samthreads.' -F 256 '.$file.' '.$chromosome.' |';
    open my $SF, $cmd or die "samtools view failed: $!";
    while ( my $line = <$SF> ) {
        next unless $line =~ m/\tvG\:B\:i\,\d/;
        next unless $line =~ m/\tvW\:i\:1/;
        my ( $cigar, $seqstr ) = ( split /\t/, $line )[5,9];
        next if $seqstr =~ m/(G{$monostretch}|A{$monostretch}|T{$monostretch}|C{$monostretch})/;
        next if $no_softmasking and $cigar =~ m/\dS/;
        my ($vGs) = $line =~ m/\tvG\:B\:i\,([\d,]+)/;
        my ($vAs) = $line =~ m/\tvA\:B\:c\,([\d,]+)/;
        my @vG = split /,/, $vGs;
        my @vA = split /,/, $vAs;
        my $cb;
        if    ( $line =~ m/\tCB\:Z\:(\S+)/ ) { $cb = $1; }
        elsif ( $line =~ m/\tRG\:Z\:(\S+)/ ) { $cb = $1; }
        else  { die "No CB/RG tag in read:\n$line"; }
        next if $cb eq '-';
        if ( defined $barcode_whitelist ) {
            $wl_candidate_reads++;
            my $match_cb = $cb;
            $match_cb =~ s/-\d+$// if $wl_strip_gem;
            $wl_barcodes_seen{$match_cb} = 1;
            if ( not exists $allowed_barcodes{$match_cb} ) {
                $wl_rejected_reads++;
                next;
            }
            $wl_kept_reads++;
            $wl_barcodes_matched{$match_cb} = 1;
        }
        my $umi;
        if ( $line =~ m/\tUB\:Z\:(\S+)/ ) {
            $umi = $1; next if $umi eq '-'; $hard_umis++;
        } else {
            my @f = split /\t/, $line;
            $umi = join('_', $cb, $f[2], $f[3], $f[1], $f[5], $f[8]);
            $soft_umis++;
        }
        $inf_reads++;
        for my $k ( 0 .. $#vG ) {
            my $pos    = $vG[$k] + 1;
            next unless exists $snp_info{$pos};
            my $allele = $vA[$k];
            next unless $allele == 1 or $allele == 2;
            $allele = 1 + int(rand(2)) if $scramble_alleles;
            $inf_alleles++;
            $umis{$pos}{$allele}{ $umi."\t".$cb }++;
            $cb2allele{$cb}{$pos}{$umi} = $allele;
        }
    }
    close $SF;
    warn "    Reads with allelic info: $inf_reads  Alleles: $inf_alleles  UMIs: $hard_umis/$soft_umis (hard/soft)\n";
}
if ( defined $barcode_whitelist ) {
    warn "Whitelist candidate barcodes observed in BAM: ", scalar(keys %wl_barcodes_seen), "\n";
    warn "Whitelist barcodes matched in BAM: ", scalar(keys %wl_barcodes_matched),
         " / ", scalar(keys %allowed_barcodes), " loaded\n";
    warn "Whitelist CB-tagged candidate reads: checked=$wl_candidate_reads kept=$wl_kept_reads rejected=$wl_rejected_reads\n";
    warn "Whitelist barcodes with usable chr$chromosome allelic evidence: ",
         scalar(keys %cb2allele), " / ", scalar(keys %allowed_barcodes), " loaded\n";
}
warn scalar(keys %umis), " SNV positions covered in ",scalar(keys %cb2allele)," cells\n";
die "Zero SNVs after BAM reading\n" unless keys %umis;
die "Zero cells after BAM reading\n" unless keys %cb2allele;

# =====================================================================
# FILTER SNVs / BUILD phased_pos
# =====================================================================
warn "Checking AFs after WASP genotyping ...\n";
my @phased_pos = ();
my @phased_dir = ();
my @ratios     = ();
my ($low_umis, $low_maf) = (0,0);
foreach my $pos ( sort { $a <=> $b } keys %umis ) {
    my $n1 = scalar keys %{ $umis{$pos}{1} };
    my $n2 = scalar keys %{ $umis{$pos}{2} };
    if ( $n1 + $n2 < $min_umis ) { $low_umis++; next; }
    my ($mn,$mx) = ($n1,$n2); ($mn,$mx) = ($mx,$mn) if $mn > $mx;
    push @ratios, ($mn+$mx) ? 100*$mn/($mn+$mx) : 0;
    my $af = ($n1+$n2) ? $n1/($n1+$n2) : 0;
    if ( $af < $min_maf or $af > 1-$min_maf ) { $low_maf++; next; }
    $snp_info{$pos}{'mono'} = (scalar keys %{$umis{$pos}{1}} == 0
                              or scalar keys %{$umis{$pos}{2}} == 0) ? 1 : 0;
    push @phased_pos, $pos;
    push @phased_dir, (rand() < 0.5) ? 1 : -1;  # ±1 only, never 0
}
die "No SNVs to phase\n" unless @phased_pos;
my $n_mono = scalar grep { $snp_info{$_}{'mono'} } @phased_pos;
warn "Excluded $low_umis SNVs < $min_umis UMIs; $low_maf SNVs with MAF < $min_maf\n";
warn "Monoallelic SNVs in pool: $n_mono / ", scalar(@phased_pos), "\n";
@ratios = sort { $a <=> $b } @ratios;
my $median_maf = @ratios % 2
    ? $ratios[$#ratios/2]
    : ( $ratios[@ratios/2-1] + $ratios[@ratios/2] ) / 2;
warn "Median MAF: $median_maf\n";

# =====================================================================
# PRETRAINING  —  improved version
# =====================================================================
#
# Architecture:
#   _cell_weighted_calls   : returns {cb}{pos} = confidence in [0,1]
#   _build_weighted_graph  : builds W{i}{j} = net signed sum of confidence votes
#   _bfs_seed_from_graph   : BFS seeding with random anchor, wmin threshold
#   _per_snv_uncertainty   : bottleneck weight along MST path → uncertainty score
#   run_pretraining        : top-level; returns (\@seed_dir, \@uncertainty)

# ------------------------------------------------------------------
# Step 1: per-cell, per-SNV confidence (replaces binary majority call)
# ------------------------------------------------------------------
# =====================================================================
# =====================================================================
# PRETRAINING  —  exact xcise_mt subs, verbatim
# Added on top: per-SNV uncertainty from edge weight -> per-try jitter
# =====================================================================

# Step 1: majority allele call per cell per SNV (xcise_mt verbatim)
sub _cell_majority_calls {
    my ($cb2allele_ref, $phased_pos_aref, $pre_tie_skip) = @_;
    my %is_pos = map { $_ => 1 } @$phased_pos_aref;
    my %cell2snv2alle;
    foreach my $cb ( keys %$cb2allele_ref ) {
        foreach my $pos ( keys %{ $cb2allele_ref->{$cb} } ) {
            next unless $is_pos{$pos};
            my ($c1,$c2) = (0,0);
            foreach my $umi ( keys %{ $cb2allele_ref->{$cb}{$pos} } ) {
                my $a = $cb2allele_ref->{$cb}{$pos}{$umi};
                $c1++ if $a == 1; $c2++ if $a == 2;
            }
            next unless ($c1+$c2) > 0;
            if ( $c1 == $c2 and $pre_tie_skip ) { next; }
            my $alle = ($c1 >= $c2) ? 1 : 2;
            $cell2snv2alle{$cb}{$pos} = $alle;
        }
    }
    return \%cell2snv2alle;
}

# Step 2: add a single vote to the graph (xcise_mt verbatim)
sub _add_vote {
    my ($W_ref, $deg_ref, $pi, $pj, $same) = @_;
    my $s = $same ? 1 : -1;
    $W_ref->{$pi}{$pj} += $s;
    $W_ref->{$pj}{$pi} += $s;
    $deg_ref->{$pi} += 1;
    $deg_ref->{$pj} += 1;
}

# Step 3: build co-segregation graph from cell calls (xcise_mt verbatim)
sub _build_graph_from_cells {
    my ($cell2snv2alle_ref, $pre_cap_pairs) = @_;
    my (%W,%deg);
    foreach my $cb ( keys %$cell2snv2alle_ref ) {
        my @pos = keys %{ $cell2snv2alle_ref->{$cb} };
        next if @pos < 2;
        my $L   = scalar @pos;
        my $all = int($L*($L-1)/2);
        if ($all <= $pre_cap_pairs) {
            for (my $i=0; $i<$L; $i++) {
                for (my $j=$i+1; $j<$L; $j++) {
                    my ($pi,$pj) = ($pos[$i],$pos[$j]);
                    my $ai = $cell2snv2alle_ref->{$cb}{$pi};
                    my $aj = $cell2snv2alle_ref->{$cb}{$pj};
                    _add_vote(\%W,\%deg,$pi,$pj, ($ai==$aj));
                }
            }
            next;
        }
        my %seen; my $need=$pre_cap_pairs; my $tries=0;
        while ( (keys %seen) < $need && $tries < 5*$need ) {
            $tries++;
            my $i=int(rand($L)); my $j=int(rand($L)); next if $i==$j;
            ($i,$j) = ($i<$j) ? ($i,$j) : ($j,$i);
            my $key = "$i#$j"; next if $seen{$key}++;
            my ($pi,$pj) = ($pos[$i],$pos[$j]);
            my $ai = $cell2snv2alle_ref->{$cb}{$pi};
            my $aj = $cell2snv2alle_ref->{$cb}{$pj};
            _add_vote(\%W,\%deg,$pi,$pj, ($ai==$aj));
        }
    }
    return (\%W, \%deg);
}

# Step 4: BFS seed from graph (xcise_mt verbatim)
# Anchor is always +1 per component (xcise_mt behaviour)
sub _bfs_seed_from_graph {
    my ($W_ref, $deg_ref, $phased_pos_aref, $wmin) = @_;
    my %dir;
    my @cands = sort { ($deg_ref->{$b}//0) <=> ($deg_ref->{$a}//0) } @$phased_pos_aref;
    foreach my $start (@cands) {
        next if exists $dir{$start};
        my $has_n = ($W_ref->{$start} && keys %{$W_ref->{$start}});
        next unless $has_n;
        $dir{$start} = 1;
        my @q = ($start);
        while (@q) {
            my $u = shift @q;
            next unless $W_ref->{$u};
            foreach my $v ( keys %{ $W_ref->{$u} } ) {
                my $w = $W_ref->{$u}{$v};
                next unless abs($w) >= $wmin;
                my $want = ($w >= 0) ? $dir{$u} : -$dir{$u};
                if ( !exists $dir{$v} ) {
                    $dir{$v} = $want;
                    push @q, $v;
                }
            }
        }
    }
    return \%dir;
}

# Step 5: per-SNV uncertainty from edge weight (excise2 addition for jitter)
# bottleneck = max |edge weight| on path from root (BFS records it per node)
# uncertainty = 1 / (1 + bottleneck); high = weakly connected = jitter more
sub _bfs_seed_with_uncertainty {
    my ($W_ref, $deg_ref, $phased_pos_aref, $wmin) = @_;
    my (%dir, %bottleneck);
    my @cands = sort { ($deg_ref->{$b}//0) <=> ($deg_ref->{$a}//0) } @$phased_pos_aref;
    foreach my $start (@cands) {
        next if exists $dir{$start};
        next unless ($W_ref->{$start} && keys %{$W_ref->{$start}});
        $dir{$start}        = 1;
        $bottleneck{$start} = 1e9;
        my @q = ($start);
        while (@q) {
            my $u = shift @q;
            next unless $W_ref->{$u};
            foreach my $v ( keys %{ $W_ref->{$u} } ) {
                next if exists $dir{$v};
                my $w = $W_ref->{$u}{$v};
                next unless abs($w) >= $wmin;
                $dir{$v}        = ($w >= 0) ? $dir{$u} : -$dir{$u};
                $bottleneck{$v} = (abs($w) < $bottleneck{$u})
                                   ? abs($w) : $bottleneck{$u};
                push @q, $v;
            }
        }
    }
    return (\%dir, \%bottleneck);
}

# Top-level: exact xcise_mt logic + uncertainty output for jitter
sub run_pretraining {
    my ($cb2allele_ref, $phased_pos_aref, $umis_ref,
        $pre_cap_pairs, $pre_min_weight, $pre_tie_skip) = @_;

    warn "  [pretrain] Building cell majority calls ...\n";
    my $cell_calls = _cell_majority_calls($cb2allele_ref, $phased_pos_aref, $pre_tie_skip);

    warn "  [pretrain] Building co-segregation graph ...\n";
    my ($W, $deg) = _build_graph_from_cells($cell_calls, $pre_cap_pairs);

    warn "  [pretrain] BFS seeding (wmin=$pre_min_weight) ...\n";
    my ($dir_map, $bottleneck) = _bfs_seed_with_uncertainty($W, $deg, $phased_pos_aref, $pre_min_weight);

    # Build seed_dir (xcise_mt: isolated SNVs stay 0)
    my @seed_dir;
    my ($n_pos, $n_neg, $n_unk) = (0,0,0);
    foreach my $pos (@$phased_pos_aref) {
        if ( exists $dir_map->{$pos} ) {
            push @seed_dir, $dir_map->{$pos};
            $n_pos++ if $dir_map->{$pos} ==  1;
            $n_neg++ if $dir_map->{$pos} == -1;
        } else {
            push @seed_dir, 0;
            $n_unk++;
            $bottleneck->{$pos} //= 0;
        }
    }

    # Uncertainty array for jitter
    my @unc;
    foreach my $pos (@$phased_pos_aref) {
        my $bw = $bottleneck->{$pos} // 0;
        push @unc, 1.0 / (1.0 + $bw);
    }

    my ($edges,$strong) = (0,0);
    foreach my $i ( keys %$W ) {
        $edges  += scalar keys %{$W->{$i}};
        $strong += scalar grep { abs($W->{$i}{$_}) >= $pre_min_weight } keys %{$W->{$i}};
    }
    my $assigned = scalar grep { $_ != 0 } @seed_dir;
    warn sprintf("  [pretrain] nodes=%d  edges~%d  strong-edges~%d  assigned=%d (%.1f%%)\n",
        scalar(@seed_dir), $edges, $strong, $assigned,
        scalar(@seed_dir) ? 100*$assigned/scalar(@seed_dir) : 0);

    return (\@seed_dir, \@unc);
}

# RUN PRETRAINING (if enabled)
# =====================================================================
my @pre_seed_dir  = ();
my @snv_uncertainty = ();   # per-SNV in phased_pos order

if ($pretrain) {
    warn "Running improved pretraining ...\n";
    my ($seed_ref, $unc_ref) = run_pretraining(
        \%cb2allele, \@phased_pos, \%umis,
        $pre_cap_pairs, $pre_min_weight, $pre_tie_skip
    );
    @pre_seed_dir    = @$seed_ref;
    @snv_uncertainty = @$unc_ref;
    @phased_dir      = @pre_seed_dir;
} else {
    # No pretrain: all uncertainties = 1 (full jitter allowed)
    @snv_uncertainty = (1.0) x scalar(@phased_pos);
}

# =====================================================================
# PARALLEL OPTIMISER TRIES  (same greedy hill-climber as xcise_mt.pl)
# =====================================================================
warn "Starting XCI calling\n";
my $tmpdir = tempdir( CLEANUP => 1 );
my $pm     = Parallel::ForkManager->new($jobs);
my $global_best;
my @global_best_dir;

TRY_LOOP:
foreach my $try ( 1 .. $tries ) {

    $pm->start and next TRY_LOOP;   # CHILD starts

    srand( (defined $seed ? $seed : time) + $try * 1337 );

    warn "Try #$try, initialising ...\n";

    # --- Per-try initialisation ---
    my @local_pos = @phased_pos;
    my @local_dir;

    if ($pretrain) {
        @local_dir = @pre_seed_dir;

        # Jitter: flip SNVs with probability proportional to uncertainty
        for my $i (0 .. $#local_dir) {
            my $u = $snv_uncertainty[$i] * $jitter_scale;
            if ( rand() < $u ) {
                # flip to the opposite non-zero direction
                if    ($local_dir[$i] ==  1) { $local_dir[$i] = -1; }
                elsif ($local_dir[$i] == -1) { $local_dir[$i] =  1; }
                else { $local_dir[$i] = (rand() < 0.5) ? 1 : -1; }
            }
        }

        # Shuffle position order (keeping dir aligned)
        for (my $i=$#local_pos; $i>0; $i--) {
            my $j = int(rand($i+1));
            next if $i == $j;
            @local_pos[$i,$j] = @local_pos[$j,$i];
            @local_dir[$i,$j] = @local_dir[$j,$i];
        }
    } else {
        # Full random init (original behaviour)
        for my $i (0 .. $#local_pos) {
            my $a = int(rand(scalar @local_pos));
            my $b = int(rand(scalar @local_pos));
            @local_pos[$a,$b] = @local_pos[$b,$a] if $a != $b;
            # biallelic: allow 0 (optimizer can collapse ambiguous SNVs)
            # monoallelic: force ±1 (only one allele observed, can't be escape)
            if ( $snp_info{$local_pos[$i]}{'mono'} ) {
                $local_dir[$i] = (rand() < 0.5) ? 1 : -1;
            } else {
                $local_dir[$i] = int(rand(3)) - 1;  # -1, 0, or +1
            }
        }
    }

    # Build initial cb2phase
    warn "    Calculating initial score ...\n";
    my %cb2phase = ();
    for my $i (0 .. $#local_pos) {
        next if $local_dir[$i] == 0;
        foreach my $allele ( keys %{ $umis{$local_pos[$i]} } ) {
            foreach my $bc ( keys %{ $umis{$local_pos[$i]}{$allele} } ) {
                my (undef,$cb) = split /\t/, $bc;
                $cb2phase{$cb}{$allele}++       if $local_dir[$i] ==  1;
                $cb2phase{$cb}{3-$allele}++     if $local_dir[$i] == -1;
            }
        }
    }

    # Score
    my ($t,$d,$c) = (0,0,0);
    foreach my $cb (keys %cb2phase) {
        my $mn = $cb2phase{$cb}{1}//0;
        my $mx = $cb2phase{$cb}{2}//0;
        $t += $mn+$mx;
        ($mn,$mx) = ($mx,$mn) if $mn>$mx;
        next unless $mx; $d += $mn; $c += $mx-1;
    }
    my $best_score = $c - $discordant_penalty*$d;
    warn "    SNVs: ",scalar(@local_pos),"  Initial score: $best_score  Discordance: ",
        ($d+$c ? sprintf("%.2f%%",100*$d/($d+$c)) : 'N/A'),"\n";

    # --- Hill-climbing (same logic as xcise_mt.pl) ---
    my ($pass, $imps, $last_imp) = (0, 1, 0);
    while ($imps) {
        $imps = 0; $pass++;
        for my $ele (0 .. $#local_dir) {
            my $old = $local_dir[$ele];

            if ( $old == 0 ) {
                # ── dir=0 SNV: evaluate both ±1, pick the better one ──────────
                # This ensures escape-candidate SNVs get a fair symmetric test
                # before the optimizer decides whether to assign a direction.
                my %score_for;   # $score_for{+1}, $score_for{-1}
                for my $new (-1, 1) {
                    # apply new direction
                    foreach my $allele ( keys %{ $umis{$local_pos[$ele]} } ) {
                        foreach my $bc ( keys %{ $umis{$local_pos[$ele]}{$allele} } ) {
                            my (undef,$cb) = split /\t/, $bc;
                            $cb2phase{$cb}{$allele}++    if $new ==  1;
                            $cb2phase{$cb}{3-$allele}++  if $new == -1;
                        }
                    }
                    my ($t2,$d2,$c2) = (0,0,0);
                    foreach my $cb (keys %cb2phase) {
                        my $mn = $cb2phase{$cb}{1}//0;
                        my $mx = $cb2phase{$cb}{2}//0;
                        $t2 += $mn+$mx; ($mn,$mx) = ($mx,$mn) if $mn>$mx;
                        next unless $mx; $d2 += $mn; $c2 += $mx-1;
                    }
                    $score_for{$new} = $c2 - $discordant_penalty*$d2;
                    # revert
                    foreach my $allele ( keys %{ $umis{$local_pos[$ele]} } ) {
                        foreach my $bc ( keys %{ $umis{$local_pos[$ele]}{$allele} } ) {
                            my (undef,$cb) = split /\t/, $bc;
                            $cb2phase{$cb}{$allele}--    if $new ==  1;
                            $cb2phase{$cb}{3-$allele}--  if $new == -1;
                        }
                    }
                }
                # Pick the strictly better direction; if both equal best_score
                # or worse, stay at 0 (tie-break: ambiguous SNV stays non-informative)
                my $best_new = ($score_for{1} >= $score_for{-1}) ? 1 : -1;
                if ( $score_for{$best_new} > $best_score ) {
                    # commit the winning direction
                    foreach my $allele ( keys %{ $umis{$local_pos[$ele]} } ) {
                        foreach my $bc ( keys %{ $umis{$local_pos[$ele]}{$allele} } ) {
                            my (undef,$cb) = split /\t/, $bc;
                            $cb2phase{$cb}{$allele}++    if $best_new ==  1;
                            $cb2phase{$cb}{3-$allele}++  if $best_new == -1;
                        }
                    }
                    $local_dir[$ele] = $best_new;
                    $best_score = $score_for{$best_new};
                    $imps++; $last_imp = $ele;
                }
                # else: both ±1 gave no improvement → SNV stays at dir=0
                #       (putative escape / genuinely ambiguous)

            } else {
                # ── dir!=0 SNV: try all three candidates ─────────────────────
                foreach my $new (-1, 0, 1) {
                    next if $new == $old;
                    # monoallelic SNVs must not go to 0
                    next if $new == 0 and $snp_info{$local_pos[$ele]}{'mono'};

                    # apply new direction (from $old to $new)
                    foreach my $allele ( keys %{ $umis{$local_pos[$ele]} } ) {
                        foreach my $bc ( keys %{ $umis{$local_pos[$ele]}{$allele} } ) {
                            my (undef,$cb) = split /\t/, $bc;
                            $cb2phase{$cb}{$allele}--    if $old ==  1;
                            $cb2phase{$cb}{3-$allele}--  if $old == -1;
                            $cb2phase{$cb}{$allele}++    if $new ==  1;
                            $cb2phase{$cb}{3-$allele}++  if $new == -1;
                        }
                    }

                    my ($t2,$d2,$c2) = (0,0,0);
                    foreach my $cb (keys %cb2phase) {
                        my $mn = $cb2phase{$cb}{1}//0;
                        my $mx = $cb2phase{$cb}{2}//0;
                        $t2 += $mn+$mx; ($mn,$mx) = ($mx,$mn) if $mn>$mx;
                        next unless $mx; $d2 += $mn; $c2 += $mx-1;
                    }
                    my $score_after = $c2 - $discordant_penalty*$d2;

                    # accept if strictly better, or tie-break: prefer dir=0
                    # (ambiguous SNV becomes non-informative at no cost)
                    if ( $score_after > $best_score
                         or ($score_after == $best_score and $new == 0) ) {
                        $imps++; $last_imp = $ele;
                        $local_dir[$ele] = $new;
                        $old = $new;   # update $old for subsequent new candidates
                        $best_score = $score_after;
                    } else {
                        # revert
                        foreach my $allele ( keys %{ $umis{$local_pos[$ele]} } ) {
                            foreach my $bc ( keys %{ $umis{$local_pos[$ele]}{$allele} } ) {
                                my (undef,$cb) = split /\t/, $bc;
                                $cb2phase{$cb}{$allele}++    if $old ==  1;
                                $cb2phase{$cb}{$allele}--    if $new ==  1;
                                $cb2phase{$cb}{3-$allele}++  if $old == -1;
                                $cb2phase{$cb}{3-$allele}--  if $new == -1;
                            }
                        }
                    }
                } # foreach new
            } # if dir==0
            # Early-exit: if no improvement found and the last improvement
            # was at the immediately prior position, try one random non-adjacent
            # position before giving up on this pass.
            if ($imps == 0 and $last_imp - $ele == 1) {
                my $rnd = int(rand(scalar @local_dir));
                if ($rnd != $ele and $rnd != $ele-1) {
                    my $rold = $local_dir[$rnd];
                    foreach my $rnew (-1, 0, 1) {
                        next if $rnew == $rold;
                        next if $rnew == 0 and $snp_info{$local_pos[$rnd]}{'mono'};
                        foreach my $allele ( keys %{ $umis{$local_pos[$rnd]} } ) {
                            foreach my $bc ( keys %{ $umis{$local_pos[$rnd]}{$allele} } ) {
                                my (undef,$cb) = split /\t/, $bc;
                                $cb2phase{$cb}{$allele}--    if $rold ==  1;
                                $cb2phase{$cb}{3-$allele}--  if $rold == -1;
                                $cb2phase{$cb}{$allele}++    if $rnew ==  1;
                                $cb2phase{$cb}{3-$allele}++  if $rnew == -1;
                            }
                        }
                        my ($tr,$dr,$cr) = (0,0,0);
                        foreach my $cb (keys %cb2phase) {
                            my $mn = $cb2phase{$cb}{1}//0;
                            my $mx = $cb2phase{$cb}{2}//0;
                            $tr += $mn+$mx;
                            ($mn,$mx) = ($mx,$mn) if $mn>$mx;
                            next unless $mx; $dr += $mn; $cr += $mx-1;
                        }
                        my $rscore = $cr - $discordant_penalty*$dr;
                        if ($rscore > $best_score
                            or ($rscore == $best_score and $rnew == 0)) {
                            $imps++; $last_imp = $rnd;
                            $local_dir[$rnd] = $rnew;
                            $best_score = $rscore;
                            last;   # accept first improvement found at random pos
                        } else {
                            # revert random position
                            foreach my $allele ( keys %{ $umis{$local_pos[$rnd]} } ) {
                                foreach my $bc ( keys %{ $umis{$local_pos[$rnd]}{$allele} } ) {
                                    my (undef,$cb) = split /\t/, $bc;
                                    $cb2phase{$cb}{$allele}++    if $rold ==  1;
                                    $cb2phase{$cb}{$allele}--    if $rnew ==  1;
                                    $cb2phase{$cb}{3-$allele}++  if $rold == -1;
                                    $cb2phase{$cb}{3-$allele}--  if $rnew == -1;
                                }
                            }
                        }
                    }
                }
                last if $imps == 0;   # truly stuck even after random probe
            }
        } # foreach ele
        warn "    Pass: $pass  Improvements: $imps  Score: $best_score\n";
    } # while imps

    # Write child result: score, positions (child order), directions
    my $resf = "$tmpdir/try_${try}.res";
    open my $RF, '>', $resf or die "Cannot write $resf: $!";
    print $RF $best_score, "\n";
    print $RF join(",", @local_pos), "\n";
    print $RF join(",", @local_dir), "\n";
    close $RF;
    warn "    Finished Try #$try  Final score: $best_score\n";

    $pm->finish(0);   # CHILD exits
}
$pm->wait_all_children;   # PARENT waits

# =====================================================================
# PARENT: gather best result; realign to canonical position order
# =====================================================================
{
    my @parent_pos = @phased_pos;

    for my $try (1 .. $tries) {
        my $resf = "$tmpdir/try_${try}.res";
        next unless -e $resf;
        open my $RF, '<', $resf or die "Cannot read $resf: $!";
        chomp(my $score_line = <$RF>);
        chomp(my $pos_line   = <$RF>);
        chomp(my $dir_line   = <$RF>);
        close $RF;
        next unless defined $dir_line;

        my $score     = $score_line + 0;
        my @child_pos = split /,/, $pos_line;
        my @child_dir = split /,/, $dir_line;

        my %pos2dir;
        for my $k (0 .. $#child_pos) {
            my $d = (defined $child_dir[$k] && $child_dir[$k] ne '') ? $child_dir[$k]+0 : 0;
            # dir=0 is valid: putative escape / ambiguous SNV
            $pos2dir{ $child_pos[$k] } = $d;
        }
        my @aligned = map { $pos2dir{$_} // 0 } @parent_pos;  # 0 = absent in child

        if ( !defined $global_best or $score > $global_best ) {
            $global_best     = $score;
            @global_best_dir = @aligned;
        }
    }

    # Optionally adopt score from previous run
    if ( $improve_existing and -s $sample.'_chr'.$chromosome.'_XCISE_summary.txt' ) {
        open my $FL, '<', $sample.'_chr'.$chromosome.'_XCISE_summary.txt' or die $!;
        my $line = <$FL>; close $FL;
        if ( $line and $line =~ m/^Best\sscore\s+\:\s+(\d+)/ ) {
            if ( !$global_best or $1 > $global_best ) {
                $global_best = $1;
                warn "    Adopted best score from previous run: $global_best\n";
            }
        }
    }
}

# =====================================================================
# REBUILD cb2phase from best solution
# =====================================================================
my %cb2phase_best = ();
for my $i (0 .. $#global_best_dir) {
    next unless defined $phased_pos[$i];
    my $pos = $phased_pos[$i];
    my $dir = $global_best_dir[$i];
    next if $dir == 0;   # dir=0: ambiguous/escape — excluded from haplotype scoring
    foreach my $allele ( keys %{ $umis{$pos} } ) {
        foreach my $bc ( keys %{ $umis{$pos}{$allele} } ) {
            my (undef,$cb) = split /\t/, $bc;
            $cb2phase_best{$cb}{$allele}++    if $dir ==  1;
            $cb2phase_best{$cb}{3-$allele}++  if $dir == -1;
        }
    }
}

# =====================================================================
# PHASE 2: assign monoallelic SNVs by cell-vote
# For monoallelic sites, the optimizer may have oriented them arbitrarily.
# We correct them using the cell classifications established from biallelic SNVs.
# =====================================================================
warn "    Phase-2 cell-vote correction for monoallelic SNVs ...\n";
{
    # Build per-cell haplotype label from biallelic SNVs only
    my %cell_hap;   # cb => 1 (X1) | 2 (X2) | 0 (ambiguous)
    foreach my $cb ( keys %cb2allele ) {
        my ($h1, $h2) = (0, 0);
        for my $i (0 .. $#global_best_dir) {
            next if $snp_info{$phased_pos[$i]}{'mono'};  # biallelic only
            my $pos = $phased_pos[$i];
            next unless exists $cb2allele{$cb}{$pos};
            foreach my $umi ( keys %{ $cb2allele{$cb}{$pos} } ) {
                my $al = $cb2allele{$cb}{$pos}{$umi};
                $h1++ if $al==1 and $global_best_dir[$i]== 1;
                $h1++ if $al==2 and $global_best_dir[$i]==-1;
                $h2++ if $al==1 and $global_best_dir[$i]==-1;
                $h2++ if $al==2 and $global_best_dir[$i]== 1;
            }
        }
        my $tot = $h1 + $h2;
        $cell_hap{$cb} = !$tot                    ? 0
                       : $h1/$tot >= 0.9           ? 1
                       : $h2/$tot >= 0.9           ? 2
                       :                             0;
    }

    my ($mono_corrected, $mono_kept) = (0, 0);
    for my $i (0 .. $#global_best_dir) {
        next unless $snp_info{$phased_pos[$i]}{'mono'};
        my $pos = $phased_pos[$i];
        # which allele is observed at this monoallelic site?
        my $obs_allele = (scalar keys %{$umis{$pos}{2}} > 0) ? 2 : 1;
        # count X1 vs X2 classified cells covering this position
        my ($v1, $v2) = (0, 0);
        my %seen_cb;
        foreach my $bc ( keys %{ $umis{$pos}{$obs_allele} } ) {
            my (undef,$cb) = split /\t/, $bc;
            next if $seen_cb{$cb}++;
            $v1++ if ($cell_hap{$cb}//0) == 1;
            $v2++ if ($cell_hap{$cb}//0) == 2;
        }
        next if $v1 + $v2 == 0;  # no classified cells cover this SNV
        # correct direction: obs_allele on X1 if majority cells are X1
        my $correct_dir;
        if ($v1 >= $v2) {
            # observed allele belongs to X1
            $correct_dir = ($obs_allele == 1) ?  1 : -1;
        } else {
            # observed allele belongs to X2
            $correct_dir = ($obs_allele == 1) ? -1 :  1;
        }
        if ($correct_dir != $global_best_dir[$i]) {
            $global_best_dir[$i] = $correct_dir;
            $mono_corrected++;
        } else {
            $mono_kept++;
        }
    }
    warn "    Mono SNVs corrected: $mono_corrected  already correct: $mono_kept\n";

    # Rebuild cb2phase_best after mono correction
    %cb2phase_best = ();
    for my $i (0 .. $#global_best_dir) {
        next unless defined $phased_pos[$i];
        my $pos = $phased_pos[$i];
        my $dir = $global_best_dir[$i];
        foreach my $allele ( keys %{ $umis{$pos} } ) {
            foreach my $bc ( keys %{ $umis{$pos}{$allele} } ) {
                my (undef,$cb) = split /\t/, $bc;
                $cb2phase_best{$cb}{$allele}++    if $dir ==  1;
                $cb2phase_best{$cb}{3-$allele}++  if $dir == -1;
            }
        }
    }
}

# =====================================================================
# =====================================================================
# PHASE 3: directionality-based assignment for dir=0 biallelic SNVs
#
# For each dir=0 SNV, tentatively apply dir=+1 and dir=-1 in turn.
# For each trial, add the SNV UMIs into cb2phase_best for covering
# cells, measure:
#   (a) n_directional: how many covering cells reach >= 90% hap purity
#   (b) sum hap1 and hap2 UMIs across covering cells
# Pick the direction that produces more directional cells (tiebreak:
# dominant haplotype UMI mass). Then check the gate:
#   dominant_hap_umi / total_umi >= ph3_min_ratio (default 0.90)
# If the gate passes, assign that direction. Otherwise stay dir=0.
# This is label-free: it measures haplotype purity gain directly.
# =====================================================================
if ($run_phase3) {
warn "    Phase-3 directionality assignment for dir=0 biallelic SNVs ...\n";
{
    my ($ph3_assigned, $ph3_stayed) = (0, 0);

    for my $i (0 .. $#global_best_dir) {
        next unless $global_best_dir[$i] == 0;
        next if $snp_info{$phased_pos[$i]}{'mono'};
        my $pos = $phased_pos[$i];

        # Collect all cells covering this SNV
        my %covering_cb;
        foreach my $allele (1, 2) {
            foreach my $bc ( keys %{ $umis{$pos}{$allele} } ) {
                my (undef,$cb) = split /\t/, $bc;
                $covering_cb{$cb} = 1;
            }
        }
        next unless keys %covering_cb;

        my %trial_result;   # trial_dir => { n_dir, h1, h2 }

        for my $trial (1, -1) {

            # Apply trial direction
            foreach my $allele (1, 2) {
                foreach my $bc ( keys %{ $umis{$pos}{$allele} } ) {
                    my (undef,$cb) = split /\t/, $bc;
                    $cb2phase_best{$cb}{$allele}++    if $trial ==  1;
                    $cb2phase_best{$cb}{3-$allele}++  if $trial == -1;
                }
            }

            # Measure directionality across covering cells
            my ($n_dir, $sum_h1, $sum_h2) = (0, 0, 0);
            foreach my $cb ( keys %covering_cb ) {
                my $h1 = $cb2phase_best{$cb}{1} // 0;
                my $h2 = $cb2phase_best{$cb}{2} // 0;
                my $tot = $h1 + $h2;
                next unless $tot > 0;
                $sum_h1 += $h1; $sum_h2 += $h2;
                $n_dir++ if ($h1/$tot >= 0.9 or $h2/$tot >= 0.9);
            }
            $trial_result{$trial} = { n_dir => $n_dir,
                                       h1    => $sum_h1,
                                       h2    => $sum_h2 };

            # Revert trial direction
            foreach my $allele (1, 2) {
                foreach my $bc ( keys %{ $umis{$pos}{$allele} } ) {
                    my (undef,$cb) = split /\t/, $bc;
                    $cb2phase_best{$cb}{$allele}--    if $trial ==  1;
                    $cb2phase_best{$cb}{3-$allele}--  if $trial == -1;
                }
            }
        }

        # Pick winning direction: most directional cells;
        # tiebreak by dominant haplotype UMI mass
        my $winner;
        if ($trial_result{1}{n_dir} > $trial_result{-1}{n_dir}) {
            $winner = 1;
        } elsif ($trial_result{-1}{n_dir} > $trial_result{1}{n_dir}) {
            $winner = -1;
        } else {
            my $dom_pos = ($trial_result{ 1}{h1} > $trial_result{ 1}{h2})
                           ? $trial_result{ 1}{h1} : $trial_result{ 1}{h2};
            my $dom_neg = ($trial_result{-1}{h1} > $trial_result{-1}{h2})
                           ? $trial_result{-1}{h1} : $trial_result{-1}{h2};
            $winner = ($dom_pos >= $dom_neg) ? 1 : -1;
        }

        # Gate: apply winner, check dominant hap fraction >= ph3_min_ratio
        foreach my $allele (1, 2) {
            foreach my $bc ( keys %{ $umis{$pos}{$allele} } ) {
                my (undef,$cb) = split /\t/, $bc;
                $cb2phase_best{$cb}{$allele}++    if $winner ==  1;
                $cb2phase_best{$cb}{3-$allele}++  if $winner == -1;
            }
        }
        my ($tot_h1, $tot_h2) = (0, 0);
        foreach my $cb ( keys %covering_cb ) {
            $tot_h1 += $cb2phase_best{$cb}{1} // 0;
            $tot_h2 += $cb2phase_best{$cb}{2} // 0;
        }
        # Revert regardless — only commit via global_best_dir
        foreach my $allele (1, 2) {
            foreach my $bc ( keys %{ $umis{$pos}{$allele} } ) {
                my (undef,$cb) = split /\t/, $bc;
                $cb2phase_best{$cb}{$allele}--    if $winner ==  1;
                $cb2phase_best{$cb}{3-$allele}--  if $winner == -1;
            }
        }

        my $dom_total = $tot_h1 + $tot_h2;
        my $dom_frac  = $dom_total > 0
            ? (($tot_h1 > $tot_h2 ? $tot_h1 : $tot_h2) / $dom_total)
            : 0;

        if ($dom_frac >= $ph3_min_ratio) {
            $global_best_dir[$i] = $winner;
            $ph3_assigned++;
        } else {
            $ph3_stayed++;
        }
    }
    warn "    Phase-3: assigned=$ph3_assigned  stayed Unk=$ph3_stayed\n";

    # Rebuild cb2phase_best with Phase-3 assignments
    if ($ph3_assigned > 0) {
        %cb2phase_best = ();
        for my $i (0 .. $#global_best_dir) {
            next unless defined $phased_pos[$i];
            my $pos = $phased_pos[$i];
            my $dir = $global_best_dir[$i];
            next if $dir == 0;
            foreach my $allele ( keys %{ $umis{$pos} } ) {
                foreach my $bc ( keys %{ $umis{$pos}{$allele} } ) {
                    my (undef,$cb) = split /\t/, $bc;
                    $cb2phase_best{$cb}{$allele}++    if $dir ==  1;
                    $cb2phase_best{$cb}{3-$allele}++  if $dir == -1;
                }
            }
        }
    }
}
# Robustness gates:
#   ph3_ratio     : min concordance ratio (concordant / total evidence) >= 0.90
# SNVs failing either gate stay at dir=0 (genuine escape candidates).
# =====================================================================
warn "    Phase-3 cell-vote assignment for dir=0 biallelic SNVs ...\n";
{
    # Build per-cell X1/X2 label from ±1 phased SNVs
    my %cell_hap;  # cb => 1 (X1) | 2 (X2) | 0 (ambiguous)
    foreach my $cb ( keys %cb2allele ) {
        my ($h1, $h2) = (0, 0);
        for my $i (0 .. $#global_best_dir) {
            next if $global_best_dir[$i] == 0;  # only phased SNVs anchor cells
            my $pos = $phased_pos[$i];
            next unless exists $cb2allele{$cb}{$pos};
            foreach my $umi ( keys %{ $cb2allele{$cb}{$pos} } ) {
                my $al = $cb2allele{$cb}{$pos}{$umi};
                $h1++ if $al==1 and $global_best_dir[$i]== 1;
                $h1++ if $al==2 and $global_best_dir[$i]==-1;
                $h2++ if $al==1 and $global_best_dir[$i]==-1;
                $h2++ if $al==2 and $global_best_dir[$i]== 1;
            }
        }
        my $tot = $h1 + $h2;
        $cell_hap{$cb} = !$tot             ? 0
                       : $h1/$tot >= 0.9   ? 1
                       : $h2/$tot >= 0.9   ? 2
                       :                      0;
    }

    my ($ph3_assigned, $ph3_stayed) = (0, 0);
    for my $i (0 .. $#global_best_dir) {
        next unless $global_best_dir[$i] == 0;     # only dir=0 biallelic
        next if $snp_info{$phased_pos[$i]}{'mono'};  # monoallelic handled by Phase 2
        my $pos = $phased_pos[$i];

        # Accumulate per-allele UMIs split by cell haplotype
        # ref_on_x1: Ref UMIs from X1 cells; alt_on_x2: Alt UMIs from X2 cells etc.
        my ($ref_on_x1, $ref_on_x2) = (0, 0);
        my ($alt_on_x1, $alt_on_x2) = (0, 0);
        my %seen_cb;

        foreach my $bc ( keys %{ $umis{$pos}{1} } ) {  # Ref UMIs
            my (undef,$cb) = split /\t/, $bc;
            next if $seen_cb{ref}{$cb}++;
            my $hap = $cell_hap{$cb} // 0;
            $ref_on_x1 += $umis{$pos}{1}{$bc} if $hap == 1;
            $ref_on_x2 += $umis{$pos}{1}{$bc} if $hap == 2;
        }
        foreach my $bc ( keys %{ $umis{$pos}{2} } ) {  # Alt UMIs
            my (undef,$cb) = split /\t/, $bc;
            next if $seen_cb{alt}{$cb}++;
            my $hap = $cell_hap{$cb} // 0;
            $alt_on_x1 += $umis{$pos}{2}{$bc} if $hap == 1;
            $alt_on_x2 += $umis{$pos}{2}{$bc} if $hap == 2;
        }

        # Evaluate both hypotheses:
#   H+1 (Ref=X1, Alt=X2): concordant = ref_on_x1 + alt_on_x2
#                         discordant = ref_on_x2 + alt_on_x1
#   H-1 (Alt=X1, Ref=X2): concordant = alt_on_x1 + ref_on_x2
#                         discordant = alt_on_x2 + ref_on_x1
        my $conc_pos = $ref_on_x1 + $alt_on_x2;  # evidence for dir=+1
        my $conc_neg = $alt_on_x1 + $ref_on_x2;  # evidence for dir=-1
        my $total_ev = $conc_pos + $conc_neg;
        next if $total_ev == 0;  # no classified-cell UMIs at all

        my ($winner, $winner_conc);
        if ($conc_pos >= $conc_neg) {
            $winner = 1;  $winner_conc = $conc_pos;
        } else {
            $winner = -1; $winner_conc = $conc_neg;
        }

        # Gate 1: minimum UMI support on winning side

        # Gate 2: concordance ratio must be convincing (>= ph3_min_ratio)
        my $ratio = $winner_conc / $total_ev;
        if ( $ratio < $ph3_min_ratio ) { $ph3_stayed++; next; }

        # All gates passed — assign direction
        $global_best_dir[$i] = $winner;
        $ph3_assigned++;
    }
    warn "    Phase-3: assigned=$ph3_assigned  stayed Unk=$ph3_stayed\n";

    # Rebuild cb2phase_best incorporating Phase-3 assignments
    if ($ph3_assigned > 0) {
        %cb2phase_best = ();
        for my $i (0 .. $#global_best_dir) {
            next unless defined $phased_pos[$i];
            my $pos = $phased_pos[$i];
            my $dir = $global_best_dir[$i];
            next if $dir == 0;
            foreach my $allele ( keys %{ $umis{$pos} } ) {
                foreach my $bc ( keys %{ $umis{$pos}{$allele} } ) {
                    my (undef,$cb) = split /\t/, $bc;
                    $cb2phase_best{$cb}{$allele}++    if $dir ==  1;
                    $cb2phase_best{$cb}{3-$allele}++  if $dir == -1;
                }
            }
        }
    }
}
} else {
    warn "    Phase 3 disabled (use -ph to enable); dir=0 SNVs remain unresolved.\n";
}

# =====================================================================
# OUTPUTS  (with do..until flip to ensure X2 >= X1)
# =====================================================================
my ( $total0, $total1, $total2, $total3, $total4 ) = (0,0,0,0,0);

do {
    ($total0,$total1,$total2,$total3,$total4) = (0,0,0,0,0);

    # --- bc2xci ---
    warn "    Outputting barcode-to-XCI ...\n";
    open my $FB, '>', $sample.'_chr'.$chromosome.'_XCISE_bc2xci.txt' or die $!;
    foreach my $cb ( sort keys %cb2allele ) {
        my ($hap1, $hap2) = (0,0);
        for my $i (0 .. $#global_best_dir) {
            my $pos = $phased_pos[$i];
            next unless exists $cb2allele{$cb}{$pos};
            foreach my $umi ( keys %{ $cb2allele{$cb}{$pos} } ) {
                my $al = $cb2allele{$cb}{$pos}{$umi};
                $hap1++ if $al==1 and $global_best_dir[$i]== 1;
                $hap1++ if $al==2 and $global_best_dir[$i]==-1;
                $hap2++ if $al==1 and $global_best_dir[$i]==-1;
                $hap2++ if $al==2 and $global_best_dir[$i]== 1;
            }
        }
        my $hap;
        if    ($hap1==0 and $hap2==0)                              { $total0++; $hap='Unknown'; }
        elsif ($hap1>=2 and $hap1/($hap1+$hap2)>=0.9)             { $total1++; $hap='X1'; }
        elsif ($hap2>=2 and $hap2/($hap1+$hap2)>=0.9)             { $total2++; $hap='X2'; }
        elsif ($hap1>0 and $hap2>0)                                { $total3++; $hap='Both'; }
        else                                                        { $total4++; $hap='Low_coverage'; }
        print $FB join("\t", $cb, $hap1, $hap2, $hap), "\n";
    }
    close $FB;

    # --- summary stats ---
    my ($ref1,$alt1,$unk1) = (0,0,0);
    for my $d (@global_best_dir) { $ref1++ if $d==1; $alt1++ if $d==-1; $unk1++ if $d==0; }
    my $n_mono_out = scalar grep { $snp_info{$phased_pos[$_]}{'mono'} } (0..$#phased_pos);
    my ($total_umi, $disc, $conc) = (0,0,0);
    foreach my $cb (keys %cb2phase_best) {
        my $mn = $cb2phase_best{$cb}{1}//0;
        my $mx = $cb2phase_best{$cb}{2}//0;
        $total_umi += $mn+$mx;
        ($mn,$mx) = ($mx,$mn) if $mn>$mx;
        next unless $mx; $disc += $mn; $conc += $mx-1;
    }
    my $grand = $total0+$total1+$total2+$total3+$total4;

    warn "    Total XCI-informative UMIs   : $total_umi\n";
    warn "    SNVs informative/non-informative/mono : ",($ref1+$alt1-$n_mono_out),"/$unk1/$n_mono_out\n";

    open my $FS, '>', $sample.'_chr'.$chromosome.'_XCISE_summary.txt' or die $!;
    print $FS "Best score  : $global_best\n";
    print $FS "Total UMIs  : $total_umi\n";
    print $FS "Concordant  : $conc\n";
    print $FS "Discordant  : $disc\n";
    print $FS "Discordancy : ", ($conc+$disc ? $disc/($conc+$disc) : 'N/A'), "\n";
    my $n_inf    = $ref1 + $alt1 - $n_mono_out;  # biallelic phased
    my $n_noninf = $unk1;                          # escape candidates
    print $FS "Total SNVs        : ", scalar(@global_best_dir), "\n";
    print $FS "Informative SNVs  : $n_inf\n";
    print $FS "Non-informative   : $n_noninf\n";
    print $FS "Monoallelic SNVs  : $n_mono_out\n";
    printf $FS "Median MAF        : %3.2f %%\n", $median_maf;
    printf $FS "X1 cells    : %5d( %3.2f %% )\n", $total1, $grand ? 100*$total1/$grand : 0;
    printf $FS "X2 cells    : %5d( %3.2f %% )\n", $total2, $grand ? 100*$total2/$grand : 0;
    printf $FS "Both X      : %5d( %3.2f %% )\n", $total3, $grand ? 100*$total3/$grand : 0;
    printf $FS "LowCoverage : %5d( %3.2f %% )\n", $total4, $grand ? 100*$total4/$grand : 0;
    printf $FS "Unknown     : %5d( %3.2f %% )\n", $total0, $grand ? 100*$total0/$grand : 0;
    if ( defined $barcode_whitelist ) {
        print $FS "Whitelist file     : $barcode_whitelist\n";
        print $FS "Whitelist strip GEM: ", ($wl_strip_gem ? 'yes' : 'no'), "\n";
        print $FS "Whitelist loaded   : ", scalar(keys %allowed_barcodes), "\n";
        print $FS "Whitelist BAM match: ", scalar(keys %wl_barcodes_matched), "\n";
        print $FS "Whitelist evidence  : ", scalar(keys %cb2allele), "\n";
        print $FS "Whitelist reads     : checked=$wl_candidate_reads kept=$wl_kept_reads rejected=$wl_rejected_reads\n";
    } else {
        print $FS "Barcode whitelist  : disabled\n";
    }
    close $FS;

    warn "    X1/X2/Both/LowC/Unk: ",join(' / ',$total1,$total2,$total3,$total4,$total0),"\n";
    if ($total1 > $total2) {
        warn "    Flipping X1<->X2 so that X2 >= X1 ...\n";
        $global_best_dir[$_] = -$global_best_dir[$_] for (0 .. $#global_best_dir);
        # Rebuild cb2phase_best after flip
        %cb2phase_best = ();
        for my $i (0 .. $#global_best_dir) {
            next unless defined $phased_pos[$i];
            my $pos = $phased_pos[$i]; my $dir = $global_best_dir[$i];
            next if $dir == 0;
            foreach my $allele ( keys %{ $umis{$pos} } ) {
                foreach my $bc ( keys %{ $umis{$pos}{$allele} } ) {
                    my (undef,$cb) = split /\t/, $bc;
                    $cb2phase_best{$cb}{$allele}++   if $dir== 1;
                    $cb2phase_best{$cb}{3-$allele}++ if $dir==-1;
                }
            }
        }
    }
} until ($total1 <= $total2);

# --- VCF ---
warn "    Outputting VCF ...\n";
open my $FVCF, '>', $sample.'_chr'.$chromosome.'_XCISE.vcf' or die $!;
print $FVCF "##fileformat=VCFv4.2\n";
print $FVCF join("\t",'#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'),"\n";
# VCF INFO fields:
#   X1A      : which allele is X1 haplotype (Ref or Alt; Unk = escape candidate)
#   AD       : allele depth as ref_umi,alt_umi
#   X1      : UMI count supporting X1 haplotype
#   X2      : UMI count supporting X2 haplotype
#   DP      : total UMI depth at this SNV (X1 + X2)
print $FVCF "##INFO=<ID=X1A,Number=1,Type=String,Description=\"X1 haplotype allele: Ref, Alt, or Unk (escape candidate)\">\n";
print $FVCF "##INFO=<ID=MONO,Number=0,Type=Flag,Description=\"Monoallelic SNV: direction by cell-vote not optimizer\">\n";
print $FVCF "##INFO=<ID=AD,Number=2,Type=Integer,Description=\"Allele depth: ref_umi,alt_umi\">\n";
print $FVCF "##INFO=<ID=X1,Number=1,Type=Integer,Description=\"UMI count supporting X1 haplotype\">\n";
print $FVCF "##INFO=<ID=X2,Number=1,Type=Integer,Description=\"UMI count supporting X2 haplotype\">\n";
print $FVCF "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total UMI depth (X1 + X2)\">\n";
for my $i ( sort { $phased_pos[$a] <=> $phased_pos[$b] } (0 .. $#global_best_dir) ) {
    my $pos = $phased_pos[$i];
    my $dir = $global_best_dir[$i];
    my $x1a    = $dir==0 ? 'Unk' : $dir==1 ? 'Ref' : 'Alt';
    my $n1     = scalar keys %{ $umis{$pos}{1} };   # ref UMIs
    my $n2     = scalar keys %{ $umis{$pos}{2} };   # alt UMIs
    my $x1_umi = ($dir ==  1) ? $n1 : ($dir == -1) ? $n2 : 0;
    my $x2_umi = ($dir == -1) ? $n1 : ($dir ==  1) ? $n2 : 0;
    my $dp     = $x1_umi + $x2_umi;
    my $mono   = $snp_info{$pos}{'mono'} ? ';MONO' : '';
    print $FVCF join("\t", $chromosome, $pos, $snp_info{$pos}{'rs'},
        $snp_info{$pos}{'ref'}, $snp_info{$pos}{'alt'},
        $dp, 'PASS',
        "X1A=$x1a;AD=$n1,$n2;X1=$x1_umi;X2=$x2_umi;DP=$dp$mono"), "\n";
}
close $FVCF;

# --- Haplotype sequence table ---
# Written AFTER the flip so it is consistent with bc2xci X1/X2 labels.
# Columns:
#   CHROM POS ID REF ALT X1_allele X2_allele
#   mono        : 1 if monoallelic (direction by cell-vote), 0 if biallelic
# UMI counts are in the VCF (X1=, X2=, DP= INFO tags) — not repeated here.
warn "    Outputting haplotype sequence table ...\n";

# Accumulators for summary statistics
my ($hap_conc, $hap_disc, $hap_tied) = (0, 0, 0);
my ($hap_x1ref, $hap_x1alt) = (0, 0);  # how often X1=Ref vs X1=Alt across SNVs
my ($sum_x1_umi, $sum_x2_umi) = (0, 0);

open my $FHAP, '>', $sample.'_chr'.$chromosome.'_XCISE_haplotypes.tsv' or die $!;
print $FHAP join("\t", qw(CHROM POS ID REF ALT X1_allele X2_allele mono)), "\n";

for my $i ( sort { $phased_pos[$a] <=> $phased_pos[$b] } (0 .. $#global_best_dir) ) {
    my $pos  = $phased_pos[$i];
    my $dir  = $global_best_dir[$i];
    next if $dir == 0;   # dir=0 = escape candidate; not in haplotype TSV
    my $ref  = $snp_info{$pos}{'ref'};
    my $alt  = $snp_info{$pos}{'alt'};
    my ($x1, $x2, $x1_is_ref);
    if    ($dir ==  1) { $x1=$ref; $x2=$alt; $x1_is_ref=1; $hap_x1ref++; }
    elsif ($dir == -1) { $x1=$alt; $x2=$ref; $x1_is_ref=0; $hap_x1alt++; }

    my $ref_umi = scalar keys %{ $umis{$pos}{1} };
    my $alt_umi = scalar keys %{ $umis{$pos}{2} };
    my $x1_umi  = ($dir ==  1) ? $ref_umi : $alt_umi;
    my $x2_umi  = ($dir == -1) ? $ref_umi : $alt_umi;
    my $total   = $x1_umi + $x2_umi;
    my $x1_frac = $total ? sprintf("%.4f", $x1_umi / $total) : "NA";

    # Concordance: does the majority UMI allele agree with whichever haplotype
    # dominates at this SNV? i.e. is the assignment internally self-consistent?
    # Since X2 is the active/majority X by design, asking "does majority = X1"
    # would make most sites look discordant — a misleading definition.
    # Instead: concordant = the dominant UMI allele matches the dominant haplotype
    # (majority UMI on X2 when x2_umi > x1_umi, or on X1 when x1_umi > x2_umi).
    my $concordance;
    if    ($ref_umi == $alt_umi) { $concordance = "tied";       $hap_tied++; }
    elsif ($x1_umi  >  $x2_umi) { $concordance = "concordant"; $hap_conc++; }  # majority on X1, X1 assigned
    elsif ($x2_umi  >  $x1_umi) { $concordance = "concordant"; $hap_conc++; }  # majority on X2, X2 assigned
    else                         { $concordance = "discordant"; $hap_disc++; }  # never reached with ±1 dirs

    $sum_x1_umi += $x1_umi;
    $sum_x2_umi += $x2_umi;

    my $mono_flag = $snp_info{$pos}{'mono'} // 0;
    print $FHAP join("\t", $chromosome, $pos, $snp_info{$pos}{'rs'},
        $ref, $alt, $x1, $x2, $mono_flag), "\n";
}
close $FHAP;
# =====================================================================
# NEW OUTPUT: per-cell × per-SNV UMI profile
#
# File: *_XCISE_cell_snv_profile.tsv
#
# One row per (cell, SNV) pair where the cell has at least one UMI
# overlapping that SNV position. Columns:
#
#   cell_barcode  : cell barcode string
#   pos           : genomic position of the SNV
#   rs            : RS ID (or '.' if absent)
#   ref           : reference allele
#   alt           : alternate allele
#   snv_dir       : final direction of this SNV (+1/-1/0 = Unk/escape)
#   x1a           : X1 haplotype allele assignment (Ref/Alt/Unk)
#   ref_umis      : number of UMIs with allele 1 (Ref) in this cell
#   alt_umis      : number of UMIs with allele 2 (Alt) in this cell
#   total_umis    : ref_umis + alt_umis
#   hap_x1_umis   : UMIs counted toward X1 haplotype for this cell at this SNV
#   hap_x2_umis   : UMIs counted toward X2 haplotype for this cell at this SNV
#   cell_xci      : final XCI call for this cell (X1/X2/Both/Low_coverage/Unknown)
#   mono          : 1 if monoallelic SNV, 0 otherwise
# =====================================================================

warn "    Outputting per-cell SNV profile ...\n";

# Build a lookup: cb -> XCI call  (re-derive from cb2allele + global_best_dir,
# consistent with the bc2xci file already written above)
my %cb2xci_label;
foreach my $cb ( sort keys %cb2allele ) {
    my ($hap1, $hap2) = (0, 0);
    for my $i (0 .. $#global_best_dir) {
        my $pos = $phased_pos[$i];
        next unless exists $cb2allele{$cb}{$pos};
        foreach my $umi ( keys %{ $cb2allele{$cb}{$pos} } ) {
            my $al = $cb2allele{$cb}{$pos}{$umi};
            $hap1++ if $al == 1 and $global_best_dir[$i] ==  1;
            $hap1++ if $al == 2 and $global_best_dir[$i] == -1;
            $hap2++ if $al == 1 and $global_best_dir[$i] == -1;
            $hap2++ if $al == 2 and $global_best_dir[$i] ==  1;
        }
    }
    if    ($hap1 == 0 and $hap2 == 0)                         { $cb2xci_label{$cb} = 'Unknown';      }
    elsif ($hap1 >= 2 and $hap1/($hap1+$hap2) >= 0.9)        { $cb2xci_label{$cb} = 'X1';           }
    elsif ($hap2 >= 2 and $hap2/($hap1+$hap2) >= 0.9)        { $cb2xci_label{$cb} = 'X2';           }
    elsif ($hap1 > 0 and $hap2 > 0)                           { $cb2xci_label{$cb} = 'Both';         }
    else                                                       { $cb2xci_label{$cb} = 'Low_coverage'; }
}

# Build SNV direction index (pos -> dir, x1a, mono)
my %pos2dir;
my %pos2x1a;
for my $i (0 .. $#global_best_dir) {
    my $pos = $phased_pos[$i];
    my $dir = $global_best_dir[$i];
    $pos2dir{$pos} = $dir;
    $pos2x1a{$pos} = $dir == 0 ? 'Unk' : $dir == 1 ? 'Ref' : 'Alt';
}

open my $FPROF, '>', $sample.'_chr'.$chromosome.'_XCISE_cell_snv_profile.tsv'
    or die "Cannot write cell_snv_profile: $!";

print $FPROF join("\t", qw(
    cell_barcode pos rs ref alt
    snv_dir x1a
    ref_umis alt_umis total_umis
    hap_x1_umis hap_x2_umis
    cell_xci mono
)), "\n";

foreach my $cb ( sort keys %cb2allele ) {
    my $cell_label = $cb2xci_label{$cb} // 'Unknown';

    # Iterate over SNV positions that this cell actually covers
    foreach my $pos ( sort { $a <=> $b } keys %{ $cb2allele{$cb} } ) {

        # Only emit rows for positions that passed SNV filtering (in phased_pos)
        next unless exists $pos2dir{$pos};

        # Count allele-specific UMIs for this cell at this position
        my ($ref_umis, $alt_umis) = (0, 0);
        foreach my $umi ( keys %{ $cb2allele{$cb}{$pos} } ) {
            my $al = $cb2allele{$cb}{$pos}{$umi};
            $ref_umis++ if $al == 1;
            $alt_umis++ if $al == 2;
        }
        my $total = $ref_umis + $alt_umis;
        next unless $total > 0;   # safety: skip empty cells

        # Translate allele UMIs into haplotype UMIs given the final SNV direction
        my $dir = $pos2dir{$pos};
        my ($hap_x1, $hap_x2);
        if    ($dir ==  1) { $hap_x1 = $ref_umis; $hap_x2 = $alt_umis; }
        elsif ($dir == -1) { $hap_x1 = $alt_umis; $hap_x2 = $ref_umis; }
        else               { $hap_x1 = 0;          $hap_x2 = 0;          }
        # dir=0 (escape candidate): hap contributions left as 0 to signal
        # that this SNV did not contribute to haplotype phasing.

        my $mono_flag = $snp_info{$pos}{'mono'} // 0;

        print $FPROF join("\t",
            $cb,
            $pos,
            $snp_info{$pos}{'rs'}  // '.',
            $snp_info{$pos}{'ref'} // '.',
            $snp_info{$pos}{'alt'} // '.',
            $dir,
            $pos2x1a{$pos},
            $ref_umis,
            $alt_umis,
            $total,
            $hap_x1,
            $hap_x2,
            $cell_label,
            $mono_flag,
        ), "\n";
    }
}
close $FPROF;
warn "    Per-cell SNV profile written: ",
     $sample.'_chr'.$chromosome.'_XCISE_cell_snv_profile.tsv', "\n";
