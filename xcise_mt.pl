#!/usr/bin/perl -w
use strict;

# NEW: minimal deps for multi-process tries and temp dir
use Parallel::ForkManager;       # NEW
use File::Temp qw(tempdir);      # NEW

# Input parameters
# -b bam (can be multiple)
# -s snv
# -r chromosome
# -o output name
# -t number of tries
# -u min number of UMIs per SNV
# -m allowed MAF per SNV
# -p discordant penalty
# -x scramble alleles randomly
# -i improve existsing solution
# -sm Do not use alignments with softmasking
# -ms Skip alignments with long mononucleotide stretches
# NEW:
# -j parallel tries (processes)
# -samthreads threads for samtools -@
# -seed base RNG seed
# -pretrain enable unsupervised seeding (NEW)
# -pre_cap per-cell cap on SNV-SNV pairs (NEW)
# -pre_wmin abs edge weight threshold for propagation (NEW)
# -pre_tieskip skip cells’ SNV calls with 1:1 ties (NEW)

warn "Reading parameters ...\n"; 
my ( $sample, $chromosome, $snv_file);
my $tries = 100;
my $min_maf = 0;
my $min_umis = 10;
my $discordant_penalty = 5;
my $scramble_alleles = 0;
my $improve_existing = 0;
my $no_softmasking = 0;
my $monostretch = 15;

# NEW: parallelism / IO knobs
my $jobs = 1;            # NEW: -j
my $samthreads = 1;      # NEW: -samthreads
my $seed;                # NEW: -seed

# NEW: pretraining knobs
my $pretrain       = 0;      # -pretrain
my $pre_cap_pairs  = 5000;   # -pre_cap
my $pre_min_weight = 3;      # -pre_wmin
my $pre_tie_skip   = 1;      # -pre_tieskip

my %bams = ();
my $ele = 0;
while ( $ele <= $#ARGV ) {
    if ( $ARGV[$ele] eq '-b' ) { #BAM files
        while ( $ele < $#ARGV and $ARGV[$ele+1] !~ m/^\-/) {
            die 'BAM file '.$ARGV[$ele+1].' does not exist' unless -e $ARGV[$ele+1];
            $bams{$ARGV[$ele+1]}=1;
            $ele++;
        }
        $ele++;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-o' ) { #Output prefix
        $sample = $ARGV[$ele+1];
        $ele+=2;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-s' ) { #SNV file
        $snv_file = $ARGV[$ele+1];
        $ele+=2;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-r' ) { #Chromosome
        $chromosome = $ARGV[$ele+1];
        $ele+=2;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-t' ) { #Number of tries
        $tries = $ARGV[$ele+1];
        $ele+=2;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-u' ) { #Min number of UMIs per SNV
        $min_umis = $ARGV[$ele+1];
        $ele+=2;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-m' ) { #Allowed MAF per SNV
        $min_maf = $ARGV[$ele+1];
        $ele+=2;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-p' ) { #Discrodant penalty
        $discordant_penalty = $ARGV[$ele+1];
        $ele+=2;
    }
    elsif ( $ele <= $#ARGV and $ARGV[$ele] eq '-x' ) { #Scramble alleles
        $scramble_alleles = 1;
        $ele++;
    }
    elsif ( $ele <= $#ARGV and $ARGV[$ele] eq '-i' ) { #Improve existsing solution
        $improve_existing = 1;
        $ele++;
    }
    elsif ( $ele <= $#ARGV and $ARGV[$ele] eq '-sm' ) { #Do not use alignments with softmasking
        $no_softmasking = 1;
        $ele++;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-ms' and $ARGV[$ele+1] =~ m/^(\d+)/ ) { #Do not use reads with long mononucleotide stretches
        $monostretch = $1;
        $ele+=2;
    }
    # NEW args:
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-j' ) { $jobs = $ARGV[$ele+1]; $ele+=2; }             # NEW
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-samthreads' ) { $samthreads = $ARGV[$ele+1]; $ele+=2; } # NEW
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-seed' ) { $seed = $ARGV[$ele+1]; $ele+=2; }          # NEW
    elsif ( $ele <= $#ARGV and $ARGV[$ele] eq '-pretrain' ) { $pretrain = 1; $ele++; }              # NEW
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-pre_cap' ) { $pre_cap_pairs = $ARGV[$ele+1]; $ele+=2; }    # NEW
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-pre_wmin' ) { $pre_min_weight = $ARGV[$ele+1]; $ele+=2; }  # NEW
    elsif ( $ele <= $#ARGV and $ARGV[$ele] eq '-pre_tieskip' ) { $pre_tie_skip = 1; $ele++; }            # NEW
    else {
        die "Unexpected/incomplete parameter:$ARGV[$ele]";
    }
}
my $usage = 'Usage: perl '.$0.' -o <output_prefix> -s <vcf_file> -b <bam_file1> [bam_file2] [-r <chromosome>] [-t <tries>] [-u <min_umis>] [-m <min_maf>] [-p <penalty>] [-j <jobs>] [-samthreads <N>] [-seed <int>] [-pretrain] [-pre_cap <N>] [-pre_wmin <W>] [-pre_tieskip]'."\n"; # CHANGED
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
# NEW: clamps
$jobs = 1 if $jobs < 1;              # NEW
$samthreads = 1 if $samthreads < 1;  # NEW
$pre_cap_pairs  = 5000 if $pre_cap_pairs < 100;   # NEW
$pre_min_weight = 3 if $pre_min_weight < 1;       # NEW
srand($seed) if defined $seed;       # NEW

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

warn "Reading BAM files ...\n";
my %umis = ();
my %cb2allele = ();
my $inf_reads = 0;
my $inf_alleles = 0;
my ( $hard_umis, $soft_umis ) = ( 0, 0 );
foreach my $file ( sort keys %bams ) {
    warn "    Reading $file ...\n";
    # CHANGED: allow threaded decompression when -samthreads > 1
    my $cmd = ($samthreads > 1)
        ? 'samtools view -@ '.$samthreads.' -F 256 '.$file.' '.$chromosome.' |'
        : 'samtools view -F 256 '.$file.' '.$chromosome.' |';
    open F, $cmd or die "samtools view failed"; # CHANGED
    while ( my $line = <F> ) {
        next unless $line =~ m/\tvG\:B\:i\,\d/; # Should have alelle-specific tag
        next unless $line =~ m/\tvW\:i\:1/; # Only unambguously allele-specific reads that passed WASP
        my ( $cigar, $seqstr ) = ( split /\t/, $line )[5,9];
        next if $seqstr =~ m/(G{$monostretch}|A{$monostretch}|T{$monostretch}|C{$monostretch})/; # Skip reads with long mononucleotide stretches
        next if $no_softmasking and $cigar =~ m/\dS/; # if enabled, only use end-to-end aighned reads
        my ( $vGs ) = $line =~ m/\tvG\:B\:i\,([\d+\,]+)/;
        my @vG = split /\,/, $vGs;
        my ( $vAs ) = $line =~ m/\tvA\:B\:c\,([\d+\,]+)/;
        my @vA = split /\,/, $vAs;
        my $cb;
        if ( $line =~ m/\tCB\:Z\:(\S+)/ or $line =~ m/\tRG\:Z\:(\S+)/ ) {
            $cb = $1;
            next if $cb eq '-';
        }
        else {
            die "Found no CB, no RG tags in BAM file in read $line";
        }

        my $umi;
        if ( $line =~ m/\tUB\:Z\:(\S+)/ ) {
            $umi = $1;
            next if $umi eq '-';
            $hard_umis++;
        }
        else {
            my ( $read, $flag, $chr, $pos, $mapq, $cigar, $chr2, $pos2, $tlen ) = split /\t/, $line;
            $umi = join('_', $cb, $chr, $pos, $flag, $cigar, $tlen );
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
    my ( $min, $max ) = ( scalar(keys %{$umis{$pos}{1}}), scalar(keys %{$umis{$pos}{2}}));
    ( $min, $max ) = ( $max, $min ) if $min > $max;
    push @ratios, 100*$min/($min+$max); 
    if ( $af < $min_maf or $af > ( 1 - $min_maf ) ) {
        $low_maf++;
        next;
    }
    push @phased_pos, $pos;
    push @phased_dir, int(rand(3))-1;
}
die "No SNVs to phase\n" unless @phased_pos;
warn "Excluded $low_umis SNVs having less than $min_umis reads/UMIs\n" if $low_umis;
warn "Excluded $low_maf SNVs having less than $min_maf minor allele frequency\n" if $low_maf;
@ratios = sort {$a<=>$b} @ratios;
my $median_maf = @ratios % 2 ? $ratios[$#ratios/2] : ($ratios[@ratios/2-1]+$ratios[@ratios/2])/2;
warn "Median MAF:", $median_maf, "\n";

# ======== UNSUPERVISED PRETRAINING (optional, minimal-change) ========
# Helpers
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
                if (!exists $dir{$v}) {
                    $dir{$v} = $want;
                    push @q,$v;
                }
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
# ======== END PRETRAINING ========

# If enabled, compute seed directions now (kept minimal; only initialization changes)
my $pre_seed_ref;
if ($pretrain) {
    warn "Running unsupervised pretraining to seed directions ...\n";
    $pre_seed_ref = run_pretraining_seed(\%cb2allele, \@phased_pos, $pre_cap_pairs, $pre_min_weight, $pre_tie_skip);
    @phased_dir = @$pre_seed_ref; # seed replaces random init; unknowns remain 0
}

warn "Starting XCI calling\n";

# NEW: parallel manager + temp dir for child results
my $tmpdir = tempdir( CLEANUP => 1 );             # NEW
my $pm = Parallel::ForkManager->new($jobs);       # NEW

# We'll collect best score/dir from children here:
my $global_best;                                   # CHANGED: undef
my @global_best_dir;                               # NEW

TRY_LOOP:
foreach my $try ( 1 .. $tries ) {

    $pm->start and next TRY_LOOP;                 # NEW: CHILD starts here

    # NEW: child-local RNG seed (deterministic w/ -seed)
    srand( (defined $seed ? $seed : time) + $try * 1337 );  # NEW

    warn "Try #$try, initializing order and XCI status ...\n";
    if ($pretrain) {
        # Shuffle positions while keeping seeded directions aligned
        for (my $i=$#phased_pos; $i>0; $i--) {
            my $j = int(rand($i+1));
            next if $i==$j;
            @phased_pos[$i,$j] = @phased_pos[$j,$i];
            @phased_dir[$i,$j] = @phased_dir[$j,$i];
        }
        # DO NOT randomize directions; we keep the pretraining seed
    } else {
        foreach my $ele ( 0 .. $#phased_pos ) {
            my $ele1 = int(rand(@phased_pos));
            my $ele2 = int(rand(@phased_pos));
            @phased_pos[$ele1,$ele2] = @phased_pos[$ele2,$ele1] if $ele1 != $ele2;
            $phased_dir[$ele] = int(rand(3))-1;
        }
    }
    
    warn "    Calculating initial score ...\n";
    my ( $t, $d, $c, $hap2, $hap0, $hap1 ) = (0,0,0,0,0,0);
    my %cb2phase = ();
    foreach my $ele ( 0 .. $#phased_pos ) {
        foreach my $allele ( keys %{$umis{$phased_pos[$ele]}} ) {
            foreach my $bc ( keys %{$umis{$phased_pos[$ele]}{$allele}} ) {
                my ( $umi, $cb ) = split /\t/, $bc;
                $cb2phase{$cb}{$allele}++ if $phased_dir[$ele] == 1;
                $cb2phase{$cb}{3-$allele}++ if $phased_dir[$ele] == -1;
            }
        }
    }
    foreach my $cb ( keys %cb2phase ) {
        my $min = exists($cb2phase{$cb}{1}) ? $cb2phase{$cb}{1} : 0; 
        my $max = exists($cb2phase{$cb}{2}) ? $cb2phase{$cb}{2} : 0; 
        $t += $min + $max;
        ( $min, $max ) = ( $max, $min ) if $min > $max;
        next unless $max;
        $d += $min;
        $c += $max - 1;
    }
    my $score = $c - $discordant_penalty * $d;
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

                my ( $total, $disc, $conc ) = ( 0, 0, 0 );
                foreach my $cb ( keys %cb2phase ) {
                    my $min = exists($cb2phase{$cb}{1}) ? $cb2phase{$cb}{1} : 0; 
                    my $max = exists($cb2phase{$cb}{2}) ? $cb2phase{$cb}{2} : 0; 
                    $total += $min + $max;
                    ( $min, $max ) = ( $max, $min ) if $min > $max;
                    next unless $max;
                    $disc += $min;
                    $conc += $max - 1;
                }
                my $score_before = $conc - $discordant_penalty * $disc;
                die 'score_before ne best_score' if $score_before != $best_score;

                foreach my $allele ( keys %{$umis{$phased_pos[$ele]}} ) {
                    foreach my $bc ( keys %{$umis{$phased_pos[$ele]}{$allele}} ) {
                        my ( $umi, $cb ) = split /\t/, $bc;
                        $cb2phase{$cb}{$allele}-- if $old == 1;
                        $cb2phase{$cb}{$allele}++ if $new == 1;
                        $cb2phase{$cb}{3-$allele}-- if $old == -1;
                        $cb2phase{$cb}{3-$allele}++ if $new == -1;
                    }
                }

                ( $total, $disc, $conc ) = ( 0, 0, 0 );
                foreach my $cb ( keys %cb2phase ) {
                    my $min = exists($cb2phase{$cb}{1}) ? $cb2phase{$cb}{1} : 0; 
                    my $max = exists($cb2phase{$cb}{2}) ? $cb2phase{$cb}{2} : 0; 
                    $total += $min + $max;
                    ( $min, $max ) = ( $max, $min ) if $min > $max;
                    next unless $max;
                    $disc += $min;
                    $conc += $max - 1;
                }
                my $score_after = $conc - $discordant_penalty * $disc;

                if ( $score_after > $best_score or ( $score_after == $best_score and $new == 0 ) ) {
                    $imps++;
                    $last_imp = $ele;
                    $phased_dir[$ele] = $new;
                    $best_score = $score_after;
                }
                else {
                    #Restore hash
                    foreach my $allele ( keys %{$umis{$phased_pos[$ele]}} ) {
                        foreach my $bc ( keys %{$umis{$phased_pos[$ele]}{$allele}} ) {
                            my ( $umi, $cb ) = split /\t/, $bc;
                            $cb2phase{$cb}{$allele}++ if $old == 1;
                            $cb2phase{$cb}{$allele}-- if $new == 1;
                            $cb2phase{$cb}{3-$allele}++ if $old == -1;
                            $cb2phase{$cb}{3-$allele}-- if $new == -1;
                        }
                    }
                    my ( $total_, $disc_, $conc_ ) = ( 0, 0, 0 );
                    foreach my $cb ( keys %cb2phase ) {
                        my $min = exists($cb2phase{$cb}{1}) ? $cb2phase{$cb}{1} : 0; 
                        my $max = exists($cb2phase{$cb}{2}) ? $cb2phase{$cb}{2} : 0; 
                        $total_ += $min + $max;
                        ( $min, $max ) = ( $max, $min ) if $min > $max;
                        next unless $max;
                        $disc_ += $min;
                        $conc_ += $max - 1;
                    }
                    my $score_restored = $conc_ - $discordant_penalty * $disc_;
                    die 'score_restored ne score_before' if $score_before != $score_restored;
                } # if best score
            } # foreach -1,0,1 
            last if $imps == 0 and $last_imp - $ele == 1; #No imps so far and the ext element is where we made improvement in the last pass
        } # foreach ele
        warn "    Pass: $pass Improvements: $imps Score: $best_score\n";
    } # while imps

    # CHILD: write 3-line result: score, positions, directions  (Option B)
    my $resf = "$tmpdir/try_${try}.res";           # NEW
    open my $RF, '>', $resf or die "Cannot write $resf: $!";
    print $RF $best_score, "\n";
    print $RF join(",", @phased_pos), "\n";       # NEW: child's position order
    print $RF join(",", @phased_dir), "\n";       # NEW: child's directions aligned to positions
    close $RF;
    warn "    Finished Try #$try with Final score: $best_score\n";

    $pm->finish(0);                                # NEW: CHILD exits
}
$pm->wait_all_children;                            # NEW: PARENT waits

# PARENT: reduce best-scoring try, realign to parent order, then output once
{
    my @parent_pos = @phased_pos;   # parent’s canonical order
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

        # Build pos -> dir map and realign to parent order
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

    # Optional: compare with previous run summary (-i)
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

    # Rebuild cb2phase using the best directions, then produce outputs (same as original)
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
    my $best_score = $global_best;  # reuse original variable names below

    do {
        my ( $ref1, $alt1, $noninf ) = ( 0, 0, 0 );
        foreach my $ele ( 0 .. $#global_best_dir ) {
            if ( $global_best_dir[$ele] == 1 ) { $ref1++; }
            elsif ( $global_best_dir[$ele] == -1 ) { $alt1++; }
            else { $noninf++; }
        }
        my ( $total, $disc, $conc ) = ( 0, 0, 0 );
        foreach my $cb ( keys %cb2phase_best ) {
            my $min = exists($cb2phase_best{$cb}{1}) ? $cb2phase_best{$cb}{1} : 0; 
            my $max = exists($cb2phase_best{$cb}{2}) ? $cb2phase_best{$cb}{2} : 0; 
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
        if ( $conc+$disc ) {
            print FLOG 'Discordancy : ', $disc/($conc+$disc), "\n";
            warn '    Final discordancy rate       : ', 100*$disc/($conc+$disc), "\n";
        }
        else {
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
            if ( $hap1 == 0 and $hap2 == 0 ) {
                $total0++;
                $hap = 'Unknown';
            }
            elsif ( $hap1 >= 2 and $hap1/($hap1+$hap2) >= 0.9 ) {
                $total1++;
                $hap = 'X1';
            } 
            elsif ( $hap2 >= 2 and $hap2/($hap1+$hap2) >= 0.9 ) {
                $total2++;
                $hap = 'X2';
            } 
            elsif ( $hap1 > 0 and $hap2>0 ) {
                $total3++;
                $hap = 'Both';
            } 
            else {
                $total4++;
                $hap = 'Low_coverage';
            }
            print F join( "\t", $cb, $hap1, $hap2, $hap ), "\n";
        }
        close F;
        my $grand_total = $total0+$total1+$total2+$total3+$total4;

        warn "    Outputing summary...\n";
        open FLOG, '>>', $sample.'_chr'.$chromosome.'_XCISE_summary.txt'; # ensure totals appended if loop repeats
        printf FLOG "X1 cells    : %5d( %3.2f %% )\n", $total1, $grand_total ? 100*$total1/$grand_total : 0;
        printf FLOG "X2 cells    : %5d( %3.2f %% )\n", $total2, $grand_total ? 100*$total2/$grand_total : 0;
        printf FLOG "Both X      : %5d( %3.2f %% )\n", $total3, $grand_total ? 100*$total3/$grand_total : 0;
        printf FLOG "LowCoverage : %5d( %3.2f %% )\n", $total4, $grand_total ? 100*$total4/$grand_total : 0;
        printf FLOG "Unknown     : %5d( %3.2f %% )\n", $total0, $grand_total ? 100*$total0/$grand_total : 0;
        close FLOG;

        warn "    X1/X2/Both/LowC/Unknown cells: ", join( ' / ', $total1, $total2, $total3, $total4, $total0 ), "\n";
        warn "    Switching X1 and X2, so that there are more X2 cells than X1 cells ...\n" if $total1 > $total2;

        # Flip all directions if needed (same as original)
        foreach my $ele ( 0 .. $#global_best_dir ) { $global_best_dir[$ele] = -$global_best_dir[$ele] }
    } until ( $total1 <= $total2 );
}
