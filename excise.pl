#!/usr/bin/perl -w
use strict;

# Minimal deps for multi-process tries and temp dir
use Parallel::ForkManager;
use File::Temp qw(tempdir);

# ---------------------- CLI ARGS & DEFAULTS ----------------------
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

# Parallelism / IO
my $jobs = 1;            # -j
my $samthreads = 1;      # -samthreads
my $seed;                # -seed

# Pretraining knobs
my $pretrain       = 0;      # -pretrain
my $pre_cap_pairs  = 5000;   # -pre_cap
my $pre_min_weight = 3;      # -pre_wmin
my $pre_tie_skip   = 1;      # -pre_tieskip

# MCMC knobs
my $chain_steps = 10000;     # -cs
my $T0          = 1.0;       # -T0
my $Tmin        = 0.01;      # -Tmin
my $blocksize   = 1;         # -blocksize
my $blocktype   = 'random';  # -blocktype random|window

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
    elsif ( $ele <= $#ARGV and $ARGV[$ele] eq '-sm' ) { $no_softmasking = 1; $ele++; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-ms' and $ARGV[$ele+1] =~ m/^(\d+)/ ) { $monostretch = $1; $ele+=2; }

    # Parallel/IO/seed
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-j' ) { $jobs = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-samthreads' ) { $samthreads = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-seed' ) { $seed = $ARGV[$ele+1]; $ele+=2; }

    # Pretrain
    elsif ( $ele <= $#ARGV and $ARGV[$ele] eq '-pretrain' ) { $pretrain = 1; $ele++; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-pre_cap' ) { $pre_cap_pairs = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-pre_wmin' ) { $pre_min_weight = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele <= $#ARGV and $ARGV[$ele] eq '-pre_tieskip' ) { $pre_tie_skip = 1; $ele++; }

    # MCMC
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-cs' )        { $chain_steps = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-T0' )        { $T0          = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-Tmin' )      { $Tmin        = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-blocksize' ) { $blocksize   = $ARGV[$ele+1]; $ele+=2; }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-blocktype' ) { $blocktype   = lc($ARGV[$ele+1]); $ele+=2; }

    else { die "Unexpected/incomplete parameter:$ARGV[$ele]"; }
}

my $usage = 'Usage: perl '.$0.' -o <output_prefix> -s <vcf_file> -b <bam_file1> [bam_file2] '
          . '[-r <chromosome>] [-t <tries>] [-u <min_umis>] [-m <min_maf>] [-p <penalty>] '
          . '[-j <jobs>] [-samthreads <N>] [-seed <int>] '
          . '[-pretrain] [-pre_cap <N>] [-pre_wmin <W>] [-pre_tieskip] '
          . '[-cs <steps>] [-T0 <T0>] [-Tmin <Tmin>] [-blocksize <k>] [-blocktype random|window]'."\n";

die "No Sample name\n".$usage unless $sample;
warn "Sample: ", $sample, "\n";
die "No VCF(SNV) file\n".$usage unless $snv_file;
warn "VCF file: ", $snv_file, "\n";
die "No BAM file\n".$usage unless keys %bams;
warn "BAM files: ", scalar(keys %bams),' (',join(',', keys %bams), ")\n";
$chromosome = 'X' unless $chromosome; warn "Chromosome: ", $chromosome, "\n";
$tries = 100 if $tries <= 0; warn "Number of tries: ", $tries, "\n";
$min_umis = 10 if $min_umis <= 0; warn "Min number of UMIs per SNV: ", $min_umis, "\n";
$min_maf = 0 if $min_maf < 0 or $min_maf > 0.5; warn "Min MAF per SNV: ", $min_maf, "\n";
$discordant_penalty = 5 if $discordant_penalty < 1; warn "Discordant penalty: ", $discordant_penalty, "\n";
$jobs = 1 if $jobs < 1;
$samthreads = 1 if $samthreads < 1;
$pre_cap_pairs  = 5000 if $pre_cap_pairs < 100;
$pre_min_weight = 3 if $pre_min_weight < 1;
$chain_steps = 10000 if $chain_steps <= 0;
$T0    = 1.0  if $T0    <= 0;
$Tmin  = 0.01 if $Tmin  <= 0;
$blocksize = 1 if $blocksize < 1;
$blocktype = 'random' unless $blocktype =~ /^(random|window)$/;
srand($seed) if defined $seed;

# ---------------------- LOAD SNVs ----------------------
warn "Reading SNVs ...\n";
my %snp_info = ();
open F, $snv_file;
open F, 'gunzip -c '.$snv_file.' |' if $snv_file =~ m/\.gz$/;
while ( <F> ) {
    next if m/^\#/;
    chomp;
    my @arr = split /\t/;
    next unless $arr[0] eq $chromosome;
    $snp_info{$arr[1]}{'rs'} = $arr[2];
    $snp_info{$arr[1]}{'ref'} = $arr[3];
    $snp_info{$arr[1]}{'alt'} = $arr[4];
}
close F;
warn scalar( keys %snp_info), " SNVs were loaded from VCF file $snv_file\n";
die "Zero SNVs, run terminated" unless keys %snp_info;

# ---------------------- READ BAMs ----------------------
warn "Reading BAM files ...\n";
my %umis = ();          # %umis{pos}{allele}{"umi\tcb"}++
my %cb2allele = ();     # %cb2allele{cb}{pos}{umi}=allele(1/2)
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
        next unless $line =~ m/\tvG\:B\:i\,\d/;
        next unless $line =~ m/\tvW\:i\:1/;
        my ( $cigar, $seqstr ) = ( split /\t/, $line )[5,9];
        next if $seqstr =~ m/(G{$monostretch}|A{$monostretch}|T{$monostretch}|C{$monostretch})/;
        next if $no_softmasking and $cigar =~ m/\dS/;
        my ( $vGs ) = $line =~ m/\tvG\:B\:i\,([\d+\,]+)/;
        my @vG = split /\,/, $vGs;
        my ( $vAs ) = $line =~ m/\tvA\:B\:c\,([\d+\,]+)/;
        my @vA = split /\,/, $vAs;
        my $cb;
        if ( $line =~ m/\tCB\:Z\:(\S+)/ or $line =~ m/\tRG\:Z\:(\S+)/ ) {
            $cb = $1; next if $cb eq '-';
        } else {
            die "Found no CB, no RG tags in BAM file in read $line";
        }
        my $umi;
        if ( $line =~ m/\tUB\:Z\:(\S+)/ ) { $umi = $1; next if $umi eq '-'; $hard_umis++; }
        else {
            my ( $read, $flag, $chr, $pos, $mapq, $cigar2, $chr2, $pos2, $tlen ) = split /\t/, $line;
            $umi = join('_', $cb, $chr, $pos, $flag, $cigar2, $tlen ); $soft_umis++;
        }
        $inf_reads++;
        foreach my $i ( 0 .. $#vG ) {
            my $pos = $vG[$i]+1; next unless exists($snp_info{$pos});
            my $allele = $vA[$i]; next unless $allele == 1 or $allele == 2;
            $allele = 1+int(rand(2)) if $scramble_alleles;
            $inf_alleles++;
            $umis{$pos}{$allele}{$umi."\t".$cb}++;
            $cb2allele{$cb}{$pos}{$umi}=$allele;
        }
    }
    close F;
    warn "    Reads with allelic info: $inf_reads Alleles: $inf_alleles UMIs:$hard_umis/$soft_umis (hard/soft)\n";
}
warn scalar keys %umis, " SNV positions were covered in ", scalar keys %cb2allele, " cells were loaded from BAM files\n";
die "Zero SNVs, run terminated" unless keys %umis;
die "Zero cells, run terminated" unless keys %cb2allele;

# ---------------------- FILTER SNVs ----------------------
warn "Checking AFs after WASP genotyping ...\n";
my @phased_pos = ();
my @phased_dir = ();
my @ratios = ();
my $low_umis = 0;
my $low_maf = 0;
foreach my $pos ( sort {$a<=>$b} keys %umis ) {
    if ( scalar(keys %{$umis{$pos}{1}}) + scalar(keys %{$umis{$pos}{2}} ) < $min_umis ) { $low_umis++; next; }
    my $af = scalar(keys %{$umis{$pos}{1}}) / ( scalar(keys %{$umis{$pos}{1}}) + scalar(keys %{$umis{$pos}{2}} ) );
    my ( $min, $max ) = ( scalar(keys %{$umis{$pos}{1}}), scalar(keys %{$umis{$pos}{2}}));
    ( $min, $max ) = ( $max, $min ) if $min > $max;
    push @ratios, 100*$min/($min+$max);
    if ( $af < $min_maf or $af > ( 1 - $min_maf ) ) { $low_maf++; next; }
    push @phased_pos, $pos;
    push @phased_dir, int(rand(3))-1;   # placeholder; may be replaced by pretrain/MCMC init
}
die "No SNVs to phase\n" unless @phased_pos;
warn "Excluded $low_umis SNVs having less than $min_umis reads/UMIs\n" if $low_umis;
warn "Excluded $low_maf SNVs having less than $min_maf minor allele frequency\n" if $low_maf;
@ratios = sort {$a<=>$b} @ratios;
my $median_maf = @ratios % 2 ? $ratios[$#ratios/2] : ($ratios[@ratios/2-1]+$ratios[@ratios/2])/2;
warn "Median MAF:", $median_maf, "\n";

# ---------------------- PRETRAIN (optional) ----------------------
sub _cell_majority_calls {
    my ($cb2allele_ref, $phased_pos_aref, $pre_tie_skip) = @_;
    my %is_pos = map { $_ => 1 } @$phased_pos_aref;
    my %cell2snv2alle;
    foreach my $cb (keys %$cb2allele_ref) {
        foreach my $pos (keys %{$cb2allele_ref->{$cb}}) {
            next unless $is_pos{$pos};
            my ($c1,$c2) = (0,0);
            foreach my $umi (keys %{$cb2allele_ref->{$cb}{$pos}}) {
                my $a = $cb2allele_ref->{$cb}{$pos}{$umi};
                $c1++ if $a==1; $c2++ if $a==2;
            }
            next unless ($c1+$c2)>0;
            if ($c1==$c2 && $pre_tie_skip) { next; }
            my $alle = ($c1>=$c2) ? 1 : 2;
            $cell2snv2alle{$cb}{$pos} = $alle;
        }
    }
    return \%cell2snv2alle;
}
sub _add_vote {
    my ($W_ref, $deg_ref, $pi, $pj, $same) = @_;
    my $s = $same ? 1 : -1;
    $W_ref->{$pi}{$pj} += $s;
    $W_ref->{$pj}{$pi} += $s;
    $deg_ref->{$pi} += 1;
    $deg_ref->{$pj} += 1;
}
sub _build_graph_from_cells {
    my ($cell2snv2alle_ref, $pre_cap_pairs) = @_;
    my (%W,%deg);
    foreach my $cb (keys %$cell2snv2alle_ref) {
        my @pos = keys %{$cell2snv2alle_ref->{$cb}};
        next if @pos < 2;
        my $L = scalar @pos;
        my $all = int($L*($L-1)/2);
        if ($all <= $pre_cap_pairs) {
            for (my $i=0;$i<$L;$i++){
                for (my $j=$i+1;$j<$L;$j++){
                    my ($pi,$pj)=($pos[$i],$pos[$j]);
                    my $ai = $cell2snv2alle_ref->{$cb}{$pi};
                    my $aj = $cell2snv2alle_ref->{$cb}{$pj};
                    _add_vote(\%W,\%deg,$pi,$pj, ($ai==$aj));
                }
            }
            next;
        }
        my %seen; my $need=$pre_cap_pairs; my $tries=0;
        while ((keys %seen)<$need && $tries<5*$need){
            $tries++;
            my $i=int(rand($L)); my $j=int(rand($L)); next if $i==$j;
            ($i,$j)=($i<$j)?($i,$j):($j,$i);
            my $key="$i#$j"; next if $seen{$key}++;
            my ($pi,$pj)=($pos[$i],$pos[$j]);
            my $ai=$cell2snv2alle_ref->{$cb}{$pi};
            my $aj=$cell2snv2alle_ref->{$cb}{$pj};
            _add_vote(\%W,\%deg,$pi,$pj, ($ai==$aj));
        }
    }
    return (\%W,\%deg);
}
sub _bfs_seed_from_graph {
    my ($W_ref, $deg_ref, $phased_pos_aref, $wmin) = @_;
    my %dir;
    my @cands = sort { ($deg_ref->{$b}//0) <=> ($deg_ref->{$a}//0) } @$phased_pos_aref;
    foreach my $start (@cands) {
        next if exists $dir{$start};
        my $has_n = ($W_ref->{$start} && keys %{$W_ref->{$start}});
        next unless $has_n;
        $dir{$start} = 1;
        my @q=($start);
        while(@q){
            my $u=shift @q;
            next unless $W_ref->{$u};
            foreach my $v (keys %{$W_ref->{$u}}){
                my $w=$W_ref->{$u}{$v};
                next unless abs($w) >= $wmin;
                my $want = ($w>=0) ? $dir{$u} : -$dir{$u};
                if (!exists $dir{$v}) { $dir{$v} = $want; push @q,$v; }
            }
        }
    }
    return \%dir;
}
sub run_pretraining_seed {
    my ($cb2allele_ref, $phased_pos_aref, $pre_cap_pairs, $pre_min_weight, $pre_tie_skip) = @_;
    my $cell_calls = _cell_majority_calls($cb2allele_ref, $phased_pos_aref, $pre_tie_skip);
    my ($W,$deg)   = _build_graph_from_cells($cell_calls, $pre_cap_pairs);
    my $dir_map    = _bfs_seed_from_graph($W, $deg, $phased_pos_aref, $pre_min_weight);
    my @seed_dir   = map { exists $dir_map->{$_} ? $dir_map->{$_} : 0 } @$phased_pos_aref;

    my ($edges,$strong)=(0,0);
    foreach my $i (keys %$W) {
        $edges  += scalar(keys %{$W->{$i}});
        $strong += scalar grep { abs($W->{$i}{$_}) >= $pre_min_weight } keys %{$W->{$i}};
    }
    my $assigned = scalar(grep { $_ != 0 } @seed_dir);
    warn "Pretrain: nodes=",scalar(@$phased_pos_aref),
         " edges~",$edges," strong-edges~",$strong,
         " assigned=", $assigned, " (", sprintf("%0.1f", 100*$assigned/@$phased_pos_aref), "%)\n";
    return \@seed_dir;
}

my $pre_seed_ref;
if ($pretrain) {
    warn "Running unsupervised pretraining to seed directions ...\n";
    $pre_seed_ref = run_pretraining_seed(\%cb2allele, \@phased_pos, $pre_cap_pairs, $pre_min_weight, $pre_tie_skip);
    @phased_dir = @$pre_seed_ref; # seed replaces random init; unknowns remain 0
}

# ---------------------- MCMC HELPERS ----------------------
sub build_cb2phase {
    my ( $umis_ref, $pp_ref, $pd_ref ) = @_;
    my %cb2phase;
    for my $i ( 0 .. $#$pp_ref ) {
        my $pos = $pp_ref->[$i]; my $d = $pd_ref->[$i];
        next if $d == 0;
        foreach my $allele ( keys %{$umis_ref->{$pos}} ) {
            foreach my $bc ( keys %{$umis_ref->{$pos}{$allele}} ) {
                my ( $umi, $cb ) = split /\t/, $bc;
                if    ( $d == 1 )  { $cb2phase{$cb}{$allele}++ }
                elsif ( $d == -1 ) { $cb2phase{$cb}{3-$allele}++ }
            }
        }
    }
    return \%cb2phase;
}
sub score_cb2phase {
    my ( $cb2p_ref, $penalty ) = @_;
    my ( $disc, $conc ) = (0,0);
    foreach my $cb ( keys %$cb2p_ref ) {
        my $a1 = $cb2p_ref->{$cb}{1} // 0;
        my $a2 = $cb2p_ref->{$cb}{2} // 0;
        ($a1,$a2) = ($a2,$a1) if $a1 > $a2; # a2 max
        next unless $a2;
        $disc += $a1;
        $conc += ($a2 - 1);
    }
    return ($conc - $penalty * $disc, $disc, $conc);
}
sub apply_block_changes {
    my ( $cb2p_ref, $umis_ref, $pp_ref, $idxs_ref, $old_ref, $new_ref, $dir_ref ) = @_;
    for my $j ( 0 .. $#$idxs_ref ) {
        my $i   = $idxs_ref->[$j];
        my $old = $old_ref->[$j];
        my $new = $new_ref->[$j];
        next if $old == $new;
        my $pos = $pp_ref->[$i];
        foreach my $allele ( keys %{$umis_ref->{$pos}} ) {
            foreach my $bc ( keys %{$umis_ref->{$pos}{$allele}} ) {
                my ( $umi, $cb ) = split /\t/, $bc;
                # remove old
                $cb2p_ref->{$cb}{$allele}--     if $old == 1;
                $cb2p_ref->{$cb}{3-$allele}--   if $old == -1;
                # add new
                $cb2p_ref->{$cb}{$allele}++     if $new == 1;
                $cb2p_ref->{$cb}{3-$allele}++   if $new == -1;
            }
        }
        $dir_ref->[$i] = $new; # commit proposed state
    }
}
sub _rand_unique_indices {
    my ($n, $k) = @_;
    $k = $n if $k > $n;
    my %seen; my @idxs;
    while (@idxs < $k) {
        my $i = int(rand($n));
        next if $seen{$i}++;
        push @idxs, $i;
    }
    return @idxs;
}

# ---------------------- PARALLEL TRIES WITH MCMC ----------------------
warn "Starting XCI calling\n";
my $tmpdir = tempdir( CLEANUP => 1 );
my $pm = Parallel::ForkManager->new($jobs);

TRY_LOOP:
foreach my $try ( 1 .. $tries ) {
    $pm->start and next TRY_LOOP;  # CHILD

    # Child-local RNG
    srand( (defined $seed ? $seed : time) + $try * 1337 );

    warn "Try #$try, initializing seed and building cb2phase ...\n";

    # Local direction vector per child:
    my @local_dir;
    if ($pretrain) { @local_dir = @phased_dir; }
    else { @local_dir = map { int(rand(3)) - 1 } (0 .. $#phased_pos); }

    my $cb2phase_ref = build_cb2phase(\%umis, \@phased_pos, \@local_dir);
    my ($curr_score) = score_cb2phase($cb2phase_ref, $discordant_penalty);
    my $best_score = $curr_score;
    my @best_dir   = @local_dir;
    my $accepted   = 0;

    for ( my $s = 1; $s <= $chain_steps; $s++ ) {
        my $T = $T0 * (($Tmin/$T0) ** ($s/$chain_steps));
        $T = 1e-6 if $T <= 0;

        # propose a block
        my @idxs;
        if ($blocktype eq 'window') {
            my $w = $blocksize < @local_dir ? $blocksize : scalar(@local_dir);
            my $start = int(rand(@local_dir - $w + 1));
            @idxs = ($start .. $start + $w - 1);
        } else {
            @idxs = _rand_unique_indices(scalar(@local_dir), $blocksize);
        }

        my (@old_states, @new_states);
        for my $i (@idxs) {
            my $old = $local_dir[$i];
            my @cand = grep { $_ != $old } (-1,0,1);
            my $new = $cand[ int(rand(@cand)) ];
            push @old_states, $old; push @new_states, $new;
        }

        my ($score_before) = score_cb2phase($cb2phase_ref, $discordant_penalty);
        apply_block_changes($cb2phase_ref, \%umis, \@phased_pos, \@idxs, \@old_states, \@new_states, \@local_dir);
        my ($score_after)  = score_cb2phase($cb2phase_ref, $discordant_penalty);
        my $delta = $score_after - $score_before;

        my $accept = ($delta >= 0) ? 1 : ( exp($delta / $T) > rand() );
        if ($accept) {
            $accepted++;
            if ($score_after > $best_score) { $best_score = $score_after; @best_dir = @local_dir; }
        } else {
            apply_block_changes($cb2phase_ref, \%umis, \@phased_pos, \@idxs, \@new_states, \@old_states, \@local_dir); # revert
        }
    }

    # CHILD: write 3-line result
    my $resf = "$tmpdir/try_${try}.res";
    open my $RF, '>', $resf or die "Cannot write $resf: $!";
    print $RF $best_score, "\n";
    print $RF join(",", @phased_pos), "\n";
    print $RF join(",", @best_dir), "\n";
    close $RF;
    warn "    Finished Try #$try with Final score: $best_score\n";
    $pm->finish(0);
}
$pm->wait_all_children;  # PARENT

# ---------------------- REDUCE BEST & OUTPUTS (unchanged formats) ----------------------
{
    my $global_best;
    my @global_best_dir;
    my @parent_pos = @phased_pos;

    for my $try (1 .. $tries) {
        my $resf = "$tmpdir/try_${try}.res"; next unless -e $resf;
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
            $global_best = $score; @global_best_dir = @dir_parent_order;
        }
    }

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

    # Rebuild cb2phase with best dirs
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
    my $best_score = $global_best;

    do {
        my ( $ref1, $alt1, $noninf ) = ( 0, 0, 0 );
        foreach my $ele ( 0 .. $#global_best_dir ) {
            if    ( $global_best_dir[$ele] == 1 )  { $ref1++; }
            elsif ( $global_best_dir[$ele] == -1 ) { $alt1++; }
            else                                   { $noninf++; }
        }
        my ( $total, $disc, $conc ) = ( 0, 0, 0 );
        foreach my $cb ( keys %cb2phase_best ) {
            my $min = $cb2phase_best{$cb}{1} // 0;
            my $max = $cb2phase_best{$cb}{2} // 0;
            $total += $min + $max;
            ( $min, $max ) = ( $max, $min ) if $min > $max;
            next unless $max;
            $disc += $min;
            $conc += $max - 1;
        }

        warn '    Total XCI-informative UMIs   : ', $total, "\n";
        warn '    SNVs XCI-informative/non-inf : ', $ref1+$alt1,'/',$noninf, "\n";
        open FLOG, '>', $sample.'_chr'.$chromosome.'_XCISE_summary.txt';
        print FLOG 'Best score  : ', $best_score, "\n";
        print FLOG 'Total UMIs  : ', $total, "\n";
        print FLOG 'Concordant  : ', $conc, "\n";
        print FLOG 'Discordant  : ', $disc, "\n";
        if ( $conc+$disc ) { print FLOG 'Discordancy : ', $disc/($conc+$disc), "\n"; warn '    Final discordancy rate       : ', 100*$disc/($conc+$disc), "\n"; }
        else { print FLOG 'Discordancy : N/A', "\n"; warn '    Final discordancy rate       : N/A', "\n"; }
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
            my $x1_allele = $global_best_dir[$ele] == 0 ? 'Unk' : $global_best_dir[$ele] == 1 ? 'Ref' : $global_best_dir[$ele] == -1 ? 'Alt' : 'Err';
            print FVCF join( "\t", $chromosome, $pos, $snp_info{$pos}{'rs'}, $snp_info{$pos}{'ref'}, $snp_info{$pos}{'alt'}, scalar(keys %{$umis{$pos}{1}})+scalar(keys %{$umis{$pos}{2}}), 'PASS', 'X1A='.$x1_allele.';AD='.scalar(keys %{$umis{$pos}{1}}).','.scalar(keys %{$umis{$pos}{2}}) ),"\n";
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
            if ( $hap1 == 0 and $hap2 == 0 ) { $total0++; $hap = 'Unknown'; }
            elsif ( $hap1 >= 2 and $hap1/($hap1+$hap2) >= 0.9 ) { $total1++; $hap = 'X1'; }
            elsif ( $hap2 >= 2 and $hap2/($hap1+$hap2) >= 0.9 ) { $total2++; $hap = 'X2'; }
            elsif ( $hap1 > 0 and $hap2>0 ) { $total3++; $hap = 'Both'; }
            else { $total4++; $hap = 'Low_coverage'; }
            print F join( "\t", $cb, $hap1, $hap2, $hap ), "\n";
        }
        close F;
        my $grand_total = $total0+$total1+$total2+$total3+$total4;

        warn "    Outputing summary...\n";
        open FLOG, '>>', $sample.'_chr'.$chromosome.'_XCISE_summary.txt';
        printf FLOG "X1 cells    : %5d( %3.2f %% )\n", $total1, $grand_total ? 100*$total1/$grand_total : 0;
        printf FLOG "X2 cells    : %5d( %3.2f %% )\n", $total2, $grand_total ? 100*$total2/$grand_total : 0;
        printf FLOG "Both X      : %5d( %3.2f %% )\n", $total3, $grand_total ? 100*$total3/$grand_total : 0;
        printf FLOG "LowCoverage : %5d( %3.2f %% )\n", $total4, $grand_total ? 100*$total4/$grand_total : 0;
        printf FLOG "Unknown     : %5d( %3.2f %% )\n", $total0, $grand_total ? 100*$total0/$grand_total : 0;
        close FLOG;

        warn "    X1/X2/Both/LowC/Unknown cells: ", join( ' / ', $total1, $total2, $total3, $total4, $total0 ), "\n";
        warn "    Switching X1 and X2, so that there are more X2 cells than X1 cells ...\n" if $total1 > $total2;

        foreach my $ele ( 0 .. $#global_best_dir ) { $global_best_dir[$ele] = -$global_best_dir[$ele] }
    } until ( $total1 <= $total2 );
}
