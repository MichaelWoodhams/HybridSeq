#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;
use Bio::Tree::Draw::Cladogram;
use Bio::TreeIO;

# Run a bunch of HybridSeq simulations all the way to
# tree diagram.

# Options:
# -n#: Number of runs
# -s#: initial seed. (Will be incremented by one per run)
# -f<base>: base name for files. Will read parameters from <base>.txt
# (except for seed), will write summary stats to <base>_<seed>.stats,
# Dollo data to <base>_<seed>_Dollo.phy, dist matrix to <base>_<seed>_Dist.phy,
# tree diagram to <base>_<seed>.eps.

our ($opt_n,$opt_s,$opt_f)=(0,0,"");

sub usage() {
    print "Usage:\n  perl $0 [options]\n";
    print "Command line options:\n";
    print "    -f<base>: use <base> as base name of input and output files\n";
    print "    -n#: number of runs\n";
    print "    -s#: seed for first run\n";
}

# Read from file named in first arg, copy to file in second arg,
# replacing seed with third arg
sub replaceSeed($$$) {
    my ($infile, $outfile, $seed) = @_;
    open (IN,"<$infile") or die "Could not open $infile";
    my @lines = <IN>;
    close (IN);
    map { s/seed = .*/seed = $seed/ } @lines;
    open (OUT,">$outfile") or die "Could not write $outfile";
    print OUT @lines;
    close OUT;
}

# main code starts here
getopts('n:s:f:');
# Defaults:
my $n = $opt_n || 1;
my $start_seed = $opt_s || 1;
my $base = $opt_f || "dollo";

for my $seed ($start_seed .. ($start_seed+$n-1)) {
    replaceSeed("$base.txt","tempParameters.txt",$seed);
    open (SIM, "java -jar HybridSeq.jar -ptempParameters.txt -d$base"."_$seed"."_Dollo.phy |") or die "Could not run HybridSeq";
    my @lines = <SIM>;
    close (SIM) or die "Error in running HybridSeq";
# Optionally trim @lines here
    open (STATS, ">$base"."_$seed.stats") or die "Could not open $base"."_$seed.stats";
    print STATS @lines;
    close STATS;
    system "add_Dollo_dist.pl -f phylip $base"."_$seed"."_Dollo.phy > $base"."_$seed"."_Dist.phy";
    system "fastme -i $base"."_$seed"."_Dist.phy -o $base"."_$seed.tree";
    my $treeio = Bio::TreeIO->new('-format' => 'newick', '-file' => "$base"."_$seed.tree");
    my $tree = $treeio->next_tree;
    $treeio->close();
    my $diagram = Bio::Tree::Draw::Cladogram->new(-tree => $tree);
    $diagram->print(-file => "$base"."_$seed.eps");
}
