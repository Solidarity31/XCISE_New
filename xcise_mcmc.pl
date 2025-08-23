#!/usr/bin/perl -w
use strict;
use Parallel::ForkManager;       # for multi-process parallel tries (-j)
use File::Temp qw(tempdir);      # temp dir to collect child results
use List::Util qw(shuffle min max);

# XCISE (MCMC+MT edition) — drop-in replacement for xcise.pl
# Adds: multi-SNV proposals and parallel tries
# Outputs remain identical to original XCISE.
#
# New flags (in addition to original):
#   -cs <int>        MCMC chain steps per try         (default 10000)
#   -T0 <float>      Initial temperature              (default 1.0)
#   -Tmin <float>    Final temperature                (default 0.01)
#   -seed <int>      RNG seed                         (optional)
#   -j <int>         Parallel tries (processes)       (default 1)
#   -blocksize <int> Number of SNVs to flip per step  (default 1)
#   -blocktype <str> "random" or "window"             (default random)
#   -samthreads <int>Threads for samtools view (-@)   (default 1)
#
# Example:
# perl xcise_mcmc_multithread.pl -o P1 -s P1_hets.vcf -b P1_WASP.bam -t 80 -j 8 -cs 4000 -blocksize 3 -blocktype random -samthreads 8

warn "Reading parameters ...\n";
my ( $sample, $chromosome, $snv_file );
my $tries = 100;
my $min_maf = 0.0;
my $min_umis = 10;
my $discordant_penalty = 5;
my $scramble_alleles = 0;
my $improve_existing = 0;
my $no_softmasking = 0;
my $monostretch = 15;

# MCMC params
my $chain_steps = 10000;
my $T0 = 1.0;
my $Tmin = 0.01;
my $seed;

# New params
my $jobs = 1;                # parallel tries
my $blocksize = 1;           # SNVs per proposal
my $blocktype = 'random';    # or 'window'
my $samthreads = 1;          # samtools -@

my %bams = ();
my $ele = 0;
while ( $ele <= $#ARGV ) {
    if ( $ARGV[$ele] eq '-b' ) { # BAM files
        while ( $ele < $#ARGV and $ARGV[$ele+1] !~ m/^\-/) {
            die 'BAM file '.$ARGV[$ele+1].' does not exist' unless -e $ARGV[$ele+1];
            $bams{$ARGV[$ele+1]}=1; $ele++;
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
    elsif ( $ele <= $#ARGV and $ARGV[$ele] eq '-sm') { $no_softmasking = 1; $ele++; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-ms' and $ARGV[$ele+1] =~ m/^(\d+)/ ) { $monostretch = $1; $ele+=2; }

    # New
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-cs' )  { $chain_steps = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-T0' )  { $T0          = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-Tmin'){ $Tmin        = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-seed'){ $seed        = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-j' )   { $jobs        = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-blocksize' ) { $blocksize = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-blocktype' ) { $blocktype = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-samthreads' ) { $samthreads = $ARGV[$ele+1]; $ele+=2; }

    else { die "Unexpected/incomplete parameter:$ARGV[$ele]"; }
}
my $usage = 'Usage: perl '.$0.' -o <output_prefix> -s <vcf_file> -b <bam1> [bam2 ...] '
          . '[-r <chrom>] [-t <tries>] [-u <min_umis>] [-m <min_maf>] [-p <penalty>] '
          . '[-cs <steps>] [-T0 <T0>] [-Tmin <Tmin>] [-seed <int>] '
          . '[-j <jobs>] [-blocksize <k>] [-blocktype random|window] [-samthreads <N>]'."\n";
die "No Sample name\n".$usage unless $sample;
warn "Sample: ", $sample, "\n";
die "No VCF(SNV) file\n".$usage unless $snv_file;
warn "VCF file: ", $snv_file, "\n";
die "No BAM file\n".$usage unless keys %bams;
warn "BAM files: ", scalar(keys %bams),' (',join(',', keys %bams), ")\n";
$chromosome = 'X' unless $chromosome; warn "Chromosome: ", $chromosome, "\n";
$tries = 100 if $tries <= 0; warn "Number of tries: ", $tries, "\n";
$min_umis = 10 if $min_umis <= 0; warn "Min UMIs per SNV: ", $min_umis, "\n";
$min_maf = 0 if $min_maf < 0 or $min_maf > 0.5; warn "Min MAF per SNV: ", $min_maf, "\n";
$discordant_penalty = 5 if $discordant_penalty < 1; warn "Discordant penalty: ", $discordant_penalty, "\n";
$chain_steps = 10000 if $chain_steps <= 0; $T0 = 1.0 if $T0 <= 0; $Tmin = 0.01 if $Tmin <= 0;
$jobs = 1 if $jobs < 1; $blocksize = 1 if $blocksize < 1; $blocktype = lc($blocktype);
$samthreads = 1 if $samthreads < 1; srand($seed) if defined $seed;

# -------------------- Load SNVs --------------------
warn "Reading SNVs ...\n";
my %snp_info = ();
my $VF;  # VCF handle
if ( $snv_file =~ /\.gz$/ ) { open($VF, "-|", "gunzip", "-c", $snv_file) or die "Cannot gunzip $snv_file: $!"; }
else { open($VF, "<", $snv_file) or die "Cannot open $snv_file: $!"; }
while ( my $line = <$VF> ) {
    next if $line =~ m/^\#/; chomp $line;
    my @arr = split /\t/, $line; next unless $arr[0] eq $chromosome;
    $snp_info{$arr[1]}{'rs'}  = $arr[2];
    $snp_info{$arr[1]}{'ref'} = $arr[3];
    $snp_info{$arr[1]}{'alt'} = $arr[4];
}
close $VF;
warn scalar(keys %snp_info), " SNVs were loaded from VCF file $snv_file\n";
die "Zero SNVs, run terminated" unless keys %snp_info;

# -------------------- Read BAMs & collect UMIs --------------------
warn "Reading BAM files ...\n";
my %umis = ();          # %umis{pos}{allele}{"umi\tcb"}++
my %cb2allele = ();     # %cb2allele{cb}{pos}{umi}=allele(1/2)
my ( $inf_reads, $inf_alleles, $hard_umis, $soft_umis ) = (0,0,0,0);
foreach my $file ( sort keys %bams ) {
    warn "    Reading $file ...\n";
    my $cmd = ($samthreads>1)
        ? 'samtools view -@ '.$samthreads.' -F 256 ' . $file . ' ' . $chromosome . ' |'
        : 'samtools view -F 256 ' . $file . ' ' . $chromosome . ' |';
    open my $SF, $cmd or die "samtools view failed";
    while ( my $line = <$SF> ) {
        next unless $line =~ m/\tvG\:B\:i\,\d/;  # allele-specific positions
        next unless $line =~ m/\tvW\:i\:1/;        # WASP-passed
        my ( $cigar, $seqstr ) = ( split /\t/, $line )[5,9];
        next if $seqstr =~ m/(G{$monostretch}|A{$monostretch}|T{$monostretch}|C{$monostretch})/;
        next if $no_softmasking and $cigar =~ m/\dS/;
        my ( $vGs ) = $line =~ m/\tvG\:B\:i\,([\d+\,]+)/; my @vG = split /\,/, $vGs;
        my ( $vAs ) = $line =~ m/\tvA\:B\:c\,([\d+\,]+)/; my @vA = split /\,/, $vAs;
        my $cb;
        if ( $line =~ m/\tCB\:Z\:(\S+)/ or $line =~ m/\tRG\:Z\:(\S+)/ ) { $cb = $1; next if $cb eq '-'; }
        else { die "Found no CB/RG tag in BAM line:\n$line"; }
        my $umi;
        if ( $line =~ m/\tUB\:Z\:(\S+)/ ) { $umi = $1; next if $umi eq '-'; $hard_umis++; }
        else { my ( $read, $flag, $chr, $pos, $mapq, $cg, $chr2, $pos2, $tlen ) = split /\t/, $line; $umi = join('_', $cb, $chr, $pos, $flag, $cg, $tlen ); $soft_umis++; }
        $inf_reads++;
        for my $k ( 0 .. $#vG ) { my $pos = $vG[$k]+1; next unless exists $snp_info{$pos}; my $allele = $vA[$k]; next unless $allele == 1 or $allele == 2; $allele = 1+int(rand(2)) if $scramble_alleles; $inf_alleles++; $umis{$pos}{$allele}{$umi."\t".$cb}++; $cb2allele{$cb}{$pos}{$umi}=$allele; }
    }
    close $SF;
}
warn scalar keys %umis, " SNV positions were covered in ", scalar keys %cb2allele, " cells were loaded from BAM files\n";
die "Zero SNVs, run terminated" unless keys %umis;
die "Zero cells, run terminated" unless keys %cb2allele;

# -------------------- Select informative SNVs; median MAF --------------------
warn "Checking AFs after WASP genotyping ...\n";
my @phased_pos = (); my @phased_dir = (); my @ratios = (); my ( $low_umis, $low_maf ) = (0,0);
foreach my $pos ( sort {$a<=>$b} keys %umis ) {
    my $n1 = scalar(keys %{$umis{$pos}{1}}); my $n2 = scalar(keys %{$umis{$pos}{2}});
    if ( $n1 + $n2 < $min_umis ) { $low_umis++; next; }
    my ($minv,$maxv) = ($n1,$n2); ($minv,$maxv) = ($maxv,$minv) if $minv > $maxv;
    push @ratios, ($minv+$maxv) ? 100*$minv/($minv+$maxv) : 0;
    my $af = ($n1+$n2) ? $n1/($n1+$n2) : 0;
    if ( $af < $min_maf or $af > (1-$min_maf) ) { $low_maf++; next; }
    push @phased_pos, $pos; push @phased_dir, int(rand(3))-1;
}
die "No SNVs to phase\n" unless @phased_pos;
@ratios = sort {$a<=>$b} @ratios; my $median_maf = @ratios % 2 ? $ratios[$#ratios/2] : ($ratios[@ratios/2-1]+$ratios[@ratios/2])/2;
warn "Excluded $low_umis SNVs < $min_umis UMIs; Excluded $low_maf SNVs with minor AF < $min_maf\n";
warn "Median MAF: $median_maf\n";

# -------------------- Helpers --------------------
sub build_cb2phase {
    my ( $umis_ref, $pp_ref, $pd_ref ) = @_;
    my %cb2phase = ();
    for my $i ( 0 .. $#$pp_ref ) {
        my $pos = $pp_ref->[$i]; my $d   = $pd_ref->[$i]; next if $d == 0;
        foreach my $allele ( keys %{$umis_ref->{$pos}} ) {
            foreach my $bc ( keys %{$umis_ref->{$pos}{$allele}} ) {
                my ( $umi, $cb ) = split /\t/, $bc;
                if    ( $d == 1 ) { $cb2phase{$cb}{$allele}++ }
                elsif ( $d == -1){ $cb2phase{$cb}{3-$allele}++ }
            }
        }
    }
    return \%cb2phase;
}

sub score_cb2phase {
    my ( $cb2p_ref, $penalty ) = @_;
    my ( $disc, $conc ) = (0,0);
    foreach my $cb ( keys %$cb2p_ref ) {
        my $a1 = exists($cb2p_ref->{$cb}{1}) ? $cb2p_ref->{$cb}{1} : 0;
        my $a2 = exists($cb2p_ref->{$cb}{2}) ? $cb2p_ref->{$cb}{2} : 0;
        ($a1,$a2) = ($a2,$a1) if $a1 > $a2; next unless $a2;
        $disc += $a1; $conc += ($a2 - 1);
    }
    return ($conc - $penalty * $disc, $disc, $conc);
}

# Apply or revert a block of changes to cb2phase given indices and new states
sub apply_block_changes {
    my ( $cb2p_ref, $umis_ref, $pp_ref, $idxs_ref, $old_ref, $new_ref, $dir_ref ) = @_;
    for my $j ( 0 .. $#$idxs_ref ) {
        my $i = $idxs_ref->[$j]; my $old = $old_ref->[$j]; my $new = $new_ref->[$j];
        next if $old == $new; # nothing to do
        my $pos = $pp_ref->[$i];
        foreach my $allele ( keys %{$umis_ref->{$pos}} ) {
            foreach my $bc ( keys %{$umis_ref->{$pos}{$allele}} ) {
                my ( $umi, $cb ) = split /\t/, $bc;
                # remove old contribution
                $cb2p_ref->{$cb}{$allele}--   if $old == 1;
                $cb2p_ref->{$cb}{3-$allele}-- if $old == -1;
                # add new contribution
                $cb2p_ref->{$cb}{$allele}++   if $new == 1;
                $cb2p_ref->{$cb}{3-$allele}++ if $new == -1;
            }
        }
        $dir_ref->[$i] = $new; # update direction vector
    }
}

# -------------------- Parallel MCMC across tries --------------------
warn "Starting XCI calling (MCMC + multi-SNV proposals, parallel tries)\n";
my $tmpdir = tempdir( CLEANUP => 1 );
my $pm = Parallel::ForkManager->new($jobs);

TRY: for my $try ( 1 .. $tries ) {
    $pm->start and next;  # child process

    # Child-local RNG nudged by try id so runs differ even w/o -seed
    srand( (defined $seed ? $seed : time) + $try * 1337 );

    # initialize dirs randomly
    my @local_dir = map { int(rand(3))-1 } (0 .. $#phased_pos);
    my $cb2phase_ref = build_cb2phase(\%umis, \@phased_pos, \@local_dir);
    my ($curr_score) = score_cb2phase($cb2phase_ref, $discordant_penalty);
    my $best_score = $curr_score; my @best_dir = @local_dir;

    my $accepted = 0; # optional: acceptance counter

    for ( my $s = 1; $s <= $chain_steps; $s++ ) {
        my $T = $T0 * (($Tmin/$T0) ** ($s/$chain_steps)); $T = 1e-6 if $T <= 0;

        # ---- propose a block move ----
        my @idxs = ();
        if ( $blocktype eq 'window' ) {
            my $w = min($blocksize, scalar(@local_dir));
            my $start = int(rand( @local_dir - $w + 1 ));
            @idxs = ($start .. $start+$w-1);
        } else { # random unique indices
            my $k = min($blocksize, scalar(@local_dir));
            my @perm = shuffle( 0 .. $#local_dir );
            @idxs = @perm[0..$k-1];
        }

        my (@old_states, @new_states);
        for my $i ( @idxs ) {
            my $old = $local_dir[$i]; my @cand = grep { $_ != $old } (-1,0,1);
            my $new = $cand[ int(rand(@cand)) ]; push @old_states, $old; push @new_states, $new;
        }

        # compute score before
        my ($score_before) = score_cb2phase($cb2phase_ref, $discordant_penalty);

        # apply tentative block (update cb2phase + dir vector)
        apply_block_changes($cb2phase_ref, \%umis, \@phased_pos, \@idxs, \@old_states, \@new_states, \@local_dir);

        # score after
        my ($score_after) = score_cb2phase($cb2phase_ref, $discordant_penalty);
        my $delta = $score_after - $score_before;
        my $accept = ($delta >= 0) ? 1 : ( exp($delta / $T) > rand() );

        if ( $accept ) {
            $accepted++;
            if ( $score_after > $best_score ) { $best_score = $score_after; @best_dir = @local_dir; }
        } else {
            # revert: swap new/old and apply again
            apply_block_changes($cb2phase_ref, \%umis, \@phased_pos, \@idxs, \@new_states, \@old_states, \@local_dir);
        }

        # optional progress log (commented to keep it fast)
        # if ( $s % 2000 == 0 ) { warn "Try #$try step $s/$chain_steps best=$best_score acc=".(sprintf('%.3f',$accepted/$s))."\n"; }
    }

    # write child result
    my $resf = "$tmpdir/try_${try}.res";
    open my $RF, '>', $resf or die "Cannot write $resf: $!";
    print $RF $best_score, "\n";
    print $RF join(",", @best_dir), "\n";
    close $RF;

    $pm->finish(0);
}
$pm->wait_all_children;

# parent: gather best across tries
my ($global_best, @global_best_dir, @global_best_pos);
@global_best_pos = @phased_pos;
for my $try ( 1 .. $tries ) {
    my $resf = "$tmpdir/try_${try}.res"; next unless -e $resf;
    open my $RF, '<', $resf or die "Cannot read $resf: $!";
    my $line1 = <$RF>; my $line2 = <$RF>; close $RF; chomp $line1; chomp $line2;
    my $best_score = $line1 + 0; my @best_dir = split(/,/, $line2);
    if ( !defined $global_best or $best_score > $global_best ) { $global_best = $best_score; @global_best_dir = @best_dir; }
}

# optionally compare with previous run summary (-i)
if ( $improve_existing and -s $sample.'_chr'.$chromosome.'_XCISE_summary.txt' ) {
    open my $FLOG, $sample.'_chr'.$chromosome.'_XCISE_summary.txt';
    my $line = <$FLOG>;
    if ( $line and $line =~ m/^Best\sscore\s+\:\s+(\d+)/ ) {
        if ( !$global_best or $1 > $global_best ) { $global_best = $1; warn "    Adopted global best score from previous run: $global_best\n"; }
    }
    close $FLOG;
}

# lock in best solution
@phased_dir = @global_best_dir if @global_best_dir;

# -------------------- Outputs (same as original) --------------------
my $cb2phase_best = build_cb2phase(\%umis, \@phased_pos, \@phased_dir);
my ( $final_score, $disc, $conc ) = score_cb2phase($cb2phase_best, $discordant_penalty);
my $total_umis = $disc + $conc;

# 1) VCF
warn "    Outputting VCF...\n";
open my $FVCF, '>', $sample.'_chr'.$chromosome.'_XCISE.vcf';
print $FVCF "##fileformat=VCFv4.2\n";
print $FVCF join( "\t", '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'INFO' ), "\n";
for my $i ( sort { $phased_pos[$a] <=> $phased_pos[$b] } ( 0 .. $#phased_dir ) ) {
    my $pos = $phased_pos[$i];
    my $x1_allele = $phased_dir[$i] == 0 ? 'Unk' : $phased_dir[$i] == 1 ? 'Ref' : $phased_dir[$i] == -1 ? 'Alt' : 'Err';
    my $n1 = scalar(keys %{$umis{$pos}{1}}); my $n2 = scalar(keys %{$umis{$pos}{2}});
    print $FVCF join( "\t", $chromosome, $pos, $snp_info{$pos}{'rs'}, $snp_info{$pos}{'ref'}, $snp_info{$pos}{'alt'}, $n1+$n2, 'PASS', 'X1A='.$x1_allele.';AD='.$n1.','.$n2 ),"\n";
}
close $FVCF;

# 2) bc2xci & 3) summary (with flip loop to prefer more X2 than X1)
my ( $total0, $total1, $total2, $total3, $total4 ) = ( 0, 0, 0, 0, 0 );
do {
    ($total0, $total1, $total2, $total3, $total4) = (0,0,0,0,0);
    warn "    Outputting barcode to phase data...\n";
    open my $FB, '>', $sample.'_chr'.$chromosome.'_XCISE_bc2xci.txt';
    foreach my $cb ( sort keys %cb2allele ) {
        my ( $hap1, $hap2 ) = ( 0, 0 );
        for my $i ( 0 .. $#phased_dir ) {
            my $pos = $phased_pos[$i]; next unless exists( $cb2allele{$cb}{$pos} );
            foreach my $umi ( keys %{$cb2allele{$cb}{$pos}} ) {
                my $al = $cb2allele{$cb}{$pos}{$umi};
                $hap1++ if $al == 1 and $phased_dir[$i] == 1;
                $hap1++ if $al == 2 and $phased_dir[$i] == -1;
                $hap2++ if $al == 1 and $phased_dir[$i] == -1;
                $hap2++ if $al == 2 and $phased_dir[$i] == 1;
            }
        }
        my $hap = '?';
        if ( $hap1 == 0 and $hap2 == 0 ) { $total0++; $hap = 'Unknown'; }
        elsif ( $hap1 >= 2 and $hap1/($hap1+$hap2) >= 0.9 ) { $total1++; $hap = 'X1'; }
        elsif ( $hap2 >= 2 and $hap2/($hap1+$hap2) >= 0.9 ) { $total2++; $hap = 'X2'; }
        elsif ( $hap1 > 0 and $hap2 > 0 ) { $total3++; $hap = 'Both'; }
        else { $total4++; $hap = 'Low_coverage'; }
        print $FB join( "\t", $cb, $hap1, $hap2, $hap ), "\n";
    }
    close $FB;

    my $grand_total = $total0+$total1+$total2+$total3+$total4;
    warn "    Outputing summary...\n";
    open my $FS, '>', $sample.'_chr'.$chromosome.'_XCISE_summary.txt';
    print $FS 'Best score  : ', $final_score, "\n";
    print $FS 'Total UMIs  : ', $total_umis, "\n";
    print $FS 'Concordant  : ', $conc, "\n";
    print $FS 'Discordant  : ', $disc, "\n";
    if ( $conc+$disc ) { print $FS 'Discordancy : ', $disc/($conc+$disc), "\n"; }
    else { print $FS 'Discordancy : N/A', "\n"; }
    print $FS 'Total SNVs  : ', scalar(@phased_dir), "\n";
    printf $FS "Median MAF  : %3.2f %%\n", $median_maf;
    my $ref1 = 0; my $alt1 = 0; my $noninf = 0; for my $d ( @phased_dir ) { $ref1++ if $d==1; $alt1++ if $d==-1; $noninf++ if $d==0; }
    print $FS 'SNV Ref/Alt : ', $ref1, "\n";
    print $FS 'SNV Alt/Ref : ', $alt1, "\n";
    print $FS 'XCI-inf SNVs: ', $ref1+$alt1, "\n";
    print $FS 'Non-inf SNVs: ', $noninf, "\n";
    printf $FS "X1 cells    : %5d( %3.2f %% )\n", $total1, $grand_total ? 100*$total1/$grand_total : 0;
    printf $FS "X2 cells    : %5d( %3.2f %% )\n", $total2, $grand_total ? 100*$total2/$grand_total : 0;
    printf $FS "Both X      : %5d( %3.2f %% )\n", $total3, $grand_total ? 100*$total3/$grand_total : 0;
    printf $FS "LowCoverage : %5d( %3.2f %% )\n", $total4, $grand_total ? 100*$total4/$grand_total : 0;
    printf $FS "Unknown     : %5d( %3.2f %% )\n", $total0, $grand_total ? 100*$total0/$grand_total : 0;
    close $FS;

    warn "    X1/X2/Both/LowC/Unknown cells: ", join( ' / ', $total1, $total2, $total3, $total4, $total0 ), "\n";
    warn "    Switching X1 and X2 so there are more X2 cells than X1 ...\n" if $total1 > $total2;
    foreach my $i ( 0 .. $#phased_dir ) { $phased_dir[$i] = -$phased_dir[$i] }
} until ( $total1 <= $total2 );

warn "Done.\n";
