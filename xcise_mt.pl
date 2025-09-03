#!/usr/bin/perl -w
use strict;

# Minimal deps
use Parallel::ForkManager;       # parallel tries and intra-try reductions
use File::Temp qw(tempdir);      # temp files/dirs
use List::Util qw(min);          # small helper

# ----------------------------
# CLI args & defaults
# ----------------------------
warn "Reading parameters ...\n";
my ( $sample, $chromosome, $snv_file );
my $tries = 100;
my $min_maf = 0;
my $min_umis = 10;
my $discordant_penalty = 5;
my $scramble_alleles = 0;
my $improve_existing = 0;
my $no_softmasking = 0;
my $monostretch = 15;

# Parallelism / IO knobs
my $jobs = 1;            # -j (concurrent tries)
my $samthreads = 1;      # -samthreads (samtools -@ threads for BAM read)
my $seed;                # -seed (base RNG seed)
my $k_cell = 1;          # -k (cell-handling workers per try)

my %bams = ();
my $ele = 0;
while ( $ele <= $#ARGV ) {
    if ( $ARGV[$ele] eq '-b' ) { # BAM files
        while ( $ele < $#ARGV and $ARGV[$ele+1] !~ m/^\-/) {
            die 'BAM file '.$ARGV[$ele+1].' does not exist' unless -e $ARGV[$ele+1];
            $bams{$ARGV[$ele+1]}=1;
            $ele++;
        }
        $ele++;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-o' ) { $sample = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-s' ) { $snv_file = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-r' ) { $chromosome = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-t' ) { $tries = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-u' ) { $min_umis = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-m' ) { $min_maf = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-p' ) { $discordant_penalty = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele <= $#ARGV and $ARGV[$ele] eq '-x' ) { $scramble_alleles = 1; $ele++; }
    elsif ( $ele <= $#ARGV and $ARGV[$ele] eq '-i' ) { $improve_existing = 1; $ele++; }
    elsif ( $ele <= $#ARGV and $ARGV[$ele] eq '-sm' ) { $no_softmasking = 1; $ele++; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-ms' and $ARGV[$ele+1] =~ m/^(\d+)/ ) { $monostretch = $1; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-j' ) { $jobs = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-samthreads' ) { $samthreads = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-seed' ) { $seed = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-k' ) { $k_cell = $ARGV[$ele+1]; $ele+=2; }
    else { die "Unexpected/incomplete parameter:$ARGV[$ele]"; }
}
my $usage = 'Usage: perl '.$0.' -o <output_prefix> -s <vcf_file> -b <bam_file1> [bam_file2] [-r <chromosome>] [-t <tries>] [-u <min_umis>] [-m <min_maf>] [-p <penalty>] [-j <jobs>] [-samthreads <N>] [-seed <int>] [-k <cell_threads>]'."\n";
die "No Sample name\n".$usage unless $sample;
warn "Sample: ", $sample, "\n";
die "No VCF(SNV) file\n".$usage unless $snv_file;
warn "VCF file: ", $snv_file, "\n";
die "No BAM file\n".$usage unless keys %bams;
warn "BAM files: ", scalar(keys %bams),' (',join(',', keys %bams), ")\n";
$chromosome = 'X' unless $chromosome;
warn "Chromosome: ", $chromosome, "\n";
$tries = 100 if $tries <= 0;
warn "Number of tries: ", $tries, "\n";
$min_umis = 10 if $min_umis <= 0;
warn "Min number of UMIs per SNV: ", $min_umis, "\n";
$min_maf = 0 if $min_maf < 0 or $min_maf > 0.5;
warn "Min MAF per SNV: ", $min_maf, "\n";
$discordant_penalty = 5 if $discordant_penalty < 1;
warn "Discordant penalty: ", $discordant_penalty, "\n";
$jobs = 1 if $jobs < 1;
$samthreads = 1 if $samthreads < 1;
$k_cell = 1 if $k_cell < 1;
srand($seed) if defined $seed;

# ----------------------------
# SNVs from VCF
# ----------------------------
warn "Reading SNVs ...\n";
my %snp_info = ();
open F, $snv_file;
open F, 'gunzip -c '.$snv_file.' |' if $snv_file =~ m/\.gz$/;
while ( <F> ) {
    next if m/^\#/;
    chomp;
    my @arr = split /\t/;
    next unless $arr[0] eq $chromosome;
    $snp_info{$arr[1]}{'rs'}  = $arr[2];
    $snp_info{$arr[1]}{'ref'} = $arr[3];
    $snp_info{$arr[1]}{'alt'} = $arr[4];
}
close F;
warn scalar( keys %snp_info), " SNVs were loaded from VCF file $snv_file\n";
die "Zero SNVs, run terminated" unless keys %snp_info;

# ----------------------------
# Read BAMs, build UMIs/cell maps
# ----------------------------
warn "Reading BAM files ...\n";
my %umis = ();        # pos -> allele(1/2) -> "UMI\tCB" -> count
my %cb2allele = ();   # CB -> pos -> UMI -> allele(1/2)
my $inf_reads = 0;
my $inf_alleles = 0;
my ( $hard_umis, $soft_umis ) = ( 0, 0 );
foreach my $file ( sort keys %bams ) {
    warn "    Reading $file ...\n";
    my $cmd = ($samthreads > 1)
        ? 'samtools view -@ '.$samthreads.' -F 256 '.$file.' '.$chromosome.' |'
        : 'samtools view -F 256 '.$file.' '.$chromosome.' |';
    open F, $cmd or die "samtools view failed";
    while ( my $line = <F> ) {
        next unless $line =~ m/\tvG\:B\:i\,\d/; # allele-specific tag exists
        next unless $line =~ m/\tvW\:i\:1/;     # passed WASP unambiguous
        my ( $cigar, $seqstr ) = ( split /\t/, $line )[5,9];
        next if $seqstr =~ m/(G{$monostretch}|A{$monostretch}|T{$monostretch}|C{$monostretch})/; # long mono-nucleotide stretch
        next if $no_softmasking and $cigar =~ m/\dS/; # drop soft-clipped if requested

        my ( $vGs ) = $line =~ m/\tvG\:B\:i\,([\d+\,]+)/;
        my @vG = split /\,/, $vGs;
        my ( $vAs ) = $line =~ m/\tvA\:B\:c\,([\d+\,]+)/;
        my @vA = split /\,/, $vAs;

        my $cb;
        if ( $line =~ m/\tCB\:Z\:(\S+)/ or $line =~ m/\tRG\:Z\:(\S+)/ ) {
            $cb = $1;
            next if $cb eq '-';
        } else {
            die "Found no CB, no RG tags in BAM file in read $line";
        }

        my $umi;
        if ( $line =~ m/\tUB\:Z\:(\S+)/ ) {
            $umi = $1;
            next if $umi eq '-';
            $hard_umis++;
        } else {
            my ( $read, $flag, $chr, $pos, $mapq, $cigar2, $chr2, $pos2, $tlen ) = split /\t/, $line;
            $umi = join('_', $cb, $chr, $pos, $flag, $cigar2, $tlen );
            $soft_umis++;
        }

        $inf_reads++;
        foreach my $ele ( 0 .. $#vG ) {
            my $pos = $vG[$ele]+1;
            next unless exists($snp_info{$pos});
            my $allele = $vA[$ele];
            next unless $allele == 1 or $allele == 2;
            $allele = 1+int(rand(2)) if $scramble_alleles;
            $inf_alleles++;
            $umis{$pos}{$allele}{$umi."\t".$cb}++;
            $cb2allele{$cb}{$pos}{$umi}=$allele;
        }
    }
    close F;
    warn "    Reads with allelic info: $inf_reads Alleles: $inf_alleles UMIs:$hard_umis\/$soft_umis (hard/soft)\n";
}
warn scalar keys %umis, " SNV positions were covered in ", scalar keys %cb2allele, " cells were loaded from BAM files\n";
die "Zero SNVs, run terminated" unless keys %umis;
die "Zero cells, run terminated" unless keys %cb2allele;

# ----------------------------
# Filter SNVs by coverage/MAF; seed directions
# ----------------------------
warn "Checking AFs after WASP genotyping ...\n";
my @phased_pos = ();
my @phased_dir = ();
my @ratios = ();
my $low_umis = 0;
my $low_maf = 0;
foreach my $pos ( sort {$a<=>$b} keys %umis ) {
    if ( scalar(keys %{$umis{$pos}{1}}) + scalar(keys %{$umis{$pos}{2}} ) < $min_umis ) {
        $low_umis++;
        next;
    }
    my $af = scalar(keys %{$umis{$pos}{1}}) / ( scalar(keys %{$umis{$pos}{1}}) + scalar(keys %{$umis{$pos}{2}} ) );
    my ( $minv, $maxv ) = ( scalar(keys %{$umis{$pos}{1}}), scalar(keys %{$umis{$pos}{2}}));
    ( $minv, $maxv ) = ( $maxv, $minv ) if $minv > $maxv;
    push @ratios, 100*$minv/($minv+$maxv);
    if ( $af < $min_maf or $af > ( 1 - $min_maf ) ) {
        $low_maf++;
        next;
    }
    push @phased_pos, $pos;
    push @phased_dir, int(rand(3))-1;  # -1,0,1
}
die "No SNVs to phase\n" unless @phased_pos;
warn "Excluded $low_umis SNVs having less than $min_umis reads/UMIs\n" if $low_umis;
warn "Excluded $low_maf SNVs having less than $min_maf minor allele frequency\n" if $low_maf;
@ratios = sort {$a<=>$b} @ratios;
my $median_maf = @ratios % 2 ? $ratios[$#ratios/2] : ($ratios[@ratios/2-1]+$ratios[@ratios/2])/2;
warn "Median MAF:", $median_maf, "\n";

warn "Starting XCI calling\n";

# ----------------------------
# Helpers: parallel reduction over cells and scoring
# ----------------------------
sub _reduce_chunk_cb2phase {
    my ($cb2phase_ref, $slice_ref) = @_;
    my ($total, $disc, $conc) = (0, 0, 0);
    foreach my $cb (@$slice_ref) {
        my $minv = exists($cb2phase_ref->{$cb}{1}) ? $cb2phase_ref->{$cb}{1} : 0;
        my $maxv = exists($cb2phase_ref->{$cb}{2}) ? $cb2phase_ref->{$cb}{2} : 0;
        $total += $minv + $maxv;
        ($minv, $maxv) = ($maxv, $minv) if $minv > $maxv;
        next unless $maxv;
        $disc += $minv;
        $conc += $maxv - 1;
    }
    return ($total, $disc, $conc);
}

sub reduce_cb2phase {
    my ($cb2phase_ref, $k) = @_;
    my @cbs = keys %$cb2phase_ref;
    return (0,0,0) unless @cbs;
    if ($k <= 1 || @cbs < 2000) {          # avoid fork overhead for small jobs
        return _reduce_chunk_cb2phase($cb2phase_ref, \@cbs);
    }
    my $fm = Parallel::ForkManager->new($k);
    my $tmpd = tempdir(CLEANUP => 1);
    my $chunk = int(@cbs / $k) + 1;

    for (my $i=0; $i<$k; $i++) {
        my $start = $i * $chunk;
        last if $start >= @cbs;
        my $end   = min($start + $chunk - 1, $#cbs);
        my @slice = @cbs[$start .. $end];

        $fm->start and next;
        my ($t,$d,$c) = _reduce_chunk_cb2phase($cb2phase_ref, \@slice);
        open my $W, '>', "$tmpd/part_$i.txt" or die "write part_$i failed";
        print $W "$t\t$d\t$c\n";
        close $W;
        $fm->finish(0);
    }
    $fm->wait_all_children;

    my ($T,$D,$C) = (0,0,0);
    for (my $i=0; $i<$k; $i++) {
        my $f = "$tmpd/part_$i.txt";
        next unless -e $f;
        open my $R, '<', $f or die "read part_$i failed";
        my $line = <$R>; close $R;
        next unless defined $line;
        my ($t,$d,$c) = split /\t/, $line;
        $T += ($t//0); $D += ($d//0); $C += ($c//0);
    }
    return ($T,$D,$C);
}

sub score_from_cb2phase {
    my ($cb2phase_ref, $k, $penalty) = @_;
    my ($total, $disc, $conc) = reduce_cb2phase($cb2phase_ref, $k);
    my $score = $conc - $penalty * $disc;
    return ($total, $disc, $conc, $score);
}

# ----------------------------
# Parallel tries manager & tmp
# ----------------------------
my $tmpdir = tempdir( CLEANUP => 1 );
my $pm     = Parallel::ForkManager->new($jobs);

# Global best across children
my $global_best;            # score
my @global_best_dir;        # directions aligned to parent @phased_pos

# ----------------------------
# TRY LOOP (may be parallel with -j)
# ----------------------------
TRY_LOOP:
foreach my $try ( 1 .. $tries ) {

    $pm->start and next TRY_LOOP;  # CHILD begins

    # Child-local RNG seed
    srand( (defined $seed ? $seed : time) + $try * 1337 );

    warn "Try #$try, randomizing order and XCI status ...\n";
    foreach my $ele ( 0 .. $#phased_pos ) {
        my $ele1 = int(rand(@phased_pos));
        my $ele2 = int(rand(@phased_pos));
        @phased_pos[$ele1,$ele2] = @phased_pos[$ele2,$ele1] if $ele1 != $ele2;
        $phased_dir[$ele] = int(rand(3))-1;
    }

    warn "    Calculating initial score ...\n";
    my %cb2phase = ();
    foreach my $ele ( 0 .. $#phased_pos ) {
        foreach my $allele ( keys %{$umis{$phased_pos[$ele]}} ) {
            foreach my $bc ( keys %{$umis{$phased_pos[$ele]}{$allele}} ) {
                my ( $umi, $cb ) = split /\t/, $bc;
                $cb2phase{$cb}{$allele}++      if $phased_dir[$ele] == 1;
                $cb2phase{$cb}{3-$allele}++    if $phased_dir[$ele] == -1;
            }
        }
    }

    my ( $t, $d, $c, $score ) = score_from_cb2phase(\%cb2phase, $k_cell, $discordant_penalty);
    warn "    SNVs: ",scalar(@phased_pos)," Initial score: $score, discordance rate : ", ($d+$c ? 100*$d/($d+$c) : 0),"\n";

    my $best_score = $score;
    my ( $pass, $imps, $last_imp ) = ( 0, 1, 0 );

    while ( $imps ) {
        $imps = 0;
        $pass++;

        foreach my $ele ( 0 .. $#phased_dir ) {
            foreach my $new ( -1, 0, 1 ) {
                next if $new == $phased_dir[$ele];
                my $old = $phased_dir[$ele];

                # Score before (parallel reduce)
                my ($total_b, $disc_b, $conc_b, $score_before) =
                    score_from_cb2phase(\%cb2phase, $k_cell, $discordant_penalty);
                die 'score_before ne best_score' if $score_before != $best_score;

                # Apply candidate change for this SNV
                foreach my $allele ( keys %{$umis{$phased_pos[$ele]}} ) {
                    foreach my $bc ( keys %{$umis{$phased_pos[$ele]}{$allele}} ) {
                        my ( $umi, $cb ) = split /\t/, $bc;
                        $cb2phase{$cb}{$allele}--       if $old == 1;
                        $cb2phase{$cb}{$allele}++       if $new == 1;
                        $cb2phase{$cb}{3-$allele}--     if $old == -1;
                        $cb2phase{$cb}{3-$allele}++     if $new == -1;
                    }
                }

                # Score after (parallel reduce)
                my ($total_a, $disc_a, $conc_a, $score_after) =
                    score_from_cb2phase(\%cb2phase, $k_cell, $discordant_penalty);

                if ( $score_after > $best_score or ( $score_after == $best_score and $new == 0 ) ) {
                    $imps++;
                    $last_imp = $ele;
                    $phased_dir[$ele] = $new;
                    $best_score = $score_after;
                }
                else {
                    # Restore hash and verify
                    foreach my $allele ( keys %{$umis{$phased_pos[$ele]}} ) {
                        foreach my $bc ( keys %{$umis{$phased_pos[$ele]}{$allele}} ) {
                            my ( $umi, $cb ) = split /\t/, $bc;
                            $cb2phase{$cb}{$allele}++       if $old == 1;
                            $cb2phase{$cb}{$allele}--       if $new == 1;
                            $cb2phase{$cb}{3-$allele}++     if $old == -1;
                            $cb2phase{$cb}{3-$allele}--     if $new == -1;
                        }
                    }
                    my (undef, undef, undef, $score_restored) =
                        score_from_cb2phase(\%cb2phase, $k_cell, $discordant_penalty);
                    die 'score_restored ne score_before' if $score_before != $score_restored;
                }
            } # each new {-1,0,1}
            last if $imps == 0 and $last_imp - $ele == 1;
        } # each ele

        warn "    Pass: $pass Improvements: $imps Score: $best_score\n";
    } # while imps

    # CHILD: write results for parent reduction
    my $resf = "$tmpdir/try_${try}.res";
    open my $RF, '>', $resf or die "Cannot write $resf: $!";
    print $RF $best_score, "\n";
    print $RF join(",", @phased_pos), "\n";
    print $RF join(",", @phased_dir), "\n";
    close $RF;
    warn "    Finished Try #$try with Final score: $best_score\n";

    $pm->finish(0);  # CHILD exits
}
$pm->wait_all_children;  # PARENT waits for all tries

# ----------------------------
# PARENT: pick best try, rebuild, output
# ----------------------------
{
    my @parent_pos = @phased_pos;   # canonical parent order

    for my $try (1 .. $tries) {
        my $resf = "$tmpdir/try_${try}.res";
        next unless -e $resf;
        open my $RF, '<', $resf or die "Cannot read $resf: $!";
        my $line1 = <$RF>; my $line2 = <$RF>; my $line3 = <$RF>;
        close $RF;
        next unless defined $line1 and defined $line2 and defined $line3;
        chomp $line1; chomp $line2; chomp $line3;

        my $score = $line1 + 0;
        my @child_pos = split(/,/, $line2);
        my @child_dir = split(/,/, $line3);

        my %pos2dir;
        for (my $k=0; $k<=$#child_pos; $k++) {
            my $p = $child_pos[$k];
            my $d = (defined $child_dir[$k] && $child_dir[$k] ne '') ? $child_dir[$k]+0 : 0;
            $pos2dir{$p} = $d;
        }
        my @dir_parent_order = map { exists $pos2dir{$_} ? $pos2dir{$_} : 0 } @parent_pos;

        if ( !defined $global_best or $score > $global_best ) {
            $global_best = $score;
            @global_best_dir = @dir_parent_order;
        }
    }

    # Optionally compare with previous run summary (-i)
    if ( $improve_existing and -s $sample.'_chr'.$chromosome.'_XCISE_summary.txt' ) {
        open FLOG, $sample.'_chr'.$chromosome.'_XCISE_summary.txt';
        my $line = <FLOG>;
        die "Error reading summary file from previous run: $line" unless $line =~ m/^Best\sscore\s+\:\s+(\d+)/;
        if ( !$global_best or $1 > $global_best ) {
             $global_best = $1;
             warn "    Updated global best score from previous run: $global_best\n";
        }
        close FLOG;
    }

    # Rebuild cb2phase using best directions then output
    my %cb2phase_best = ();
    for my $ele ( 0 .. $#global_best_dir ) {
        next unless defined $phased_pos[$ele];
        foreach my $allele ( keys %{$umis{$phased_pos[$ele]}} ) {
            foreach my $bc ( keys %{$umis{$phased_pos[$ele]}{$allele}} ) {
                my ( $umi, $cb ) = split /\t/, $bc;
                $cb2phase_best{$cb}{$allele}++    if $global_best_dir[$ele] == 1;
                $cb2phase_best{$cb}{3-$allele}++  if $global_best_dir[$ele] == -1;
            }
        }
    }

    my ( $total0, $total1, $total2, $total3, $total4 );

    # Final totals via reducer (parallel)
    my ($total_um, $disc, $conc, undef) =
        score_from_cb2phase(\%cb2phase_best, $k_cell, $discordant_penalty);

    warn '    Total XCI-informative UMIs   : ', $total_um, "\n";

    # Summaries & files
    my $best_score = $global_best;
    my ($ref1, $alt1, $noninf) = (0,0,0);
    foreach my $ele ( 0 .. $#global_best_dir ) {
        if ( $global_best_dir[$ele] == 1 ) { $ref1++; }
        elsif ( $global_best_dir[$ele] == -1 ) { $alt1++; }
        else { $noninf++; }
    }
    warn '    SNVs XCI-informative/non-inf : ', $ref1+$alt1,'/',$noninf, "\n";

    open FLOG, '>', $sample.'_chr'.$chromosome.'_XCISE_summary.txt';
    print FLOG 'Best score  : ', $best_score, "\n";
    print FLOG 'Total UMIs  : ', $total_um, "\n";
    print FLOG 'Concordant  : ', $conc, "\n";
    print FLOG 'Discordant  : ', $disc, "\n";
    if ( $conc+$disc ) {
        print FLOG 'Discordancy : ', $disc/($conc+$disc), "\n";
        warn '    Final discordancy rate       : ', 100*$disc/($conc+$disc), "\n";
    } else {
        print FLOG 'Discordancy : N/A', "\n";
        warn '    Final discordancy rate       : N/A', "\n";
    }
    print FLOG 'Total SNVs  : ', scalar(@global_best_dir), "\n";
    printf FLOG "Median MAF  : %3.2f %%\n", $median_maf;
    print FLOG 'SNV Ref/Alt : ', $ref1, "\n";
    print FLOG 'SNV Alt/Ref : ', $alt1, "\n";
    print FLOG 'XCI-inf SNVs: ', $ref1+$alt1, "\n";
    print FLOG 'Non-inf SNVs: ', $noninf, "\n";

    warn "    Outputting VCF...\n";
    open FVCF, '>', $sample.'_chr'.$chromosome.'_XCISE.vcf';
    print FVCF "##fileformat=VCFv4.2\n";
    print FVCF join( "\t", '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'INFO' ), "\n";
    foreach my $ele ( sort {$phased_pos[$a]<=>$phased_pos[$b]} ( 0 .. $#global_best_dir ) ) {
        my $pos = $phased_pos[$ele];
        my $x1_allele = $global_best_dir[$ele] == 0 ? 'Unk'
                        : $global_best_dir[$ele] == 1 ? 'Ref'
                        : $global_best_dir[$ele] == -1 ? 'Alt' : 'Err';
        print FVCF join( "\t",
            $chromosome, $pos,
            $snp_info{$pos}{'rs'},
            $snp_info{$pos}{'ref'},
            $snp_info{$pos}{'alt'},
            scalar(keys %{$umis{$pos}{1}})+scalar(keys %{$umis{$pos}{2}}),
            'PASS',
            'X1A='.$x1_allele.';AD='.scalar(keys %{$umis{$pos}{1}}).','.scalar(keys %{$umis{$pos}{2}})
        ),"\n";
    }
    close FVCF;

    warn "    Outputting barcode to phase data...\n";
    ( $total0, $total1, $total2, $total3, $total4 ) = ( 0, 0, 0, 0, 0 );
    open F, '>', $sample.'_chr'.$chromosome.'_XCISE_bc2xci.txt';
    foreach my $cb ( sort keys %cb2allele ) {
        my ( $hap1, $hap2 ) = ( 0, 0 );
        foreach my $ele ( 0 .. $#global_best_dir ) {
            my $pos = $phased_pos[$ele];
            next unless exists( $cb2allele{$cb}{$pos} );
            foreach my $umi ( keys %{$cb2allele{$cb}{$pos}} ) {
                $hap1++ if $cb2allele{$cb}{$pos}{$umi} == 1 and $global_best_dir[$ele] == 1;
                $hap1++ if $cb2allele{$cb}{$pos}{$umi} == 2 and $global_best_dir[$ele] == -1;
                $hap2++ if $cb2allele{$cb}{$pos}{$umi} == 1 and $global_best_dir[$ele] == -1;
                $hap2++ if $cb2allele{$cb}{$pos}{$umi} == 2 and $global_best_dir[$ele] == 1;
            }
        }
        my $hap = '?';
        if ( $hap1 == 0 and $hap2 == 0 ) {
            $total0++; $hap = 'Unknown';
        }
        elsif ( $hap1 >= 2 and $hap1/($hap1+$hap2) >= 0.9 ) {
            $total1++; $hap = 'X1';
        }
        elsif ( $hap2 >= 2 and $hap2/($hap1+$hap2) >= 0.9 ) {
            $total2++; $hap = 'X2';
        }
        elsif ( $hap1 > 0 and $hap2 > 0 ) {
            $total3++; $hap = 'Both';
        }
        else {
            $total4++; $hap = 'Low_coverage';
        }
        print F join( "\t", $cb, $hap1, $hap2, $hap ), "\n";
    }
    close F;
    my $grand_total = $total0+$total1+$total2+$total3+$total4;

    warn "    Outputing summary...\n";
    open FLOG2, '>>', $sample.'_chr'.$chromosome.'_XCISE_summary.txt';
    printf FLOG2 "X1 cells    : %5d( %3.2f %% )\n", $total1, $grand_total ? 100*$total1/$grand_total : 0;
    printf FLOG2 "X2 cells    : %5d( %3.2f %% )\n", $total2, $grand_total ? 100*$total2/$grand_total : 0;
    printf FLOG2 "Both X      : %5d( %3.2f %% )\n", $total3, $grand_total ? 100*$total3/$grand_total : 0;
    printf FLOG2 "LowCoverage : %5d( %3.2f %% )\n", $total4, $grand_total ? 100*$total4/$grand_total : 0;
    printf FLOG2 "Unknown     : %5d( %3.2f %% )\n", $total0, $grand_total ? 100*$total0/$grand_total : 0;
    close FLOG2;

    warn "    X1/X2/Both/LowC/Unknown cells: ", join( ' / ', $total1, $total2, $total3, $total4, $total0 ), "\n";
    warn "    Switching X1 and X2, so that there are more X2 cells than X1 cells ...\n" if $total1 > $total2;

    # Flip all directions until X2 >= X1 (as in original)
    foreach my $ele ( 0 .. $#global_best_dir ) { $global_best_dir[$ele] = -$global_best_dir[$ele] }
    redo if $total1 > $total2;  # run the reporting block again after flip
}
