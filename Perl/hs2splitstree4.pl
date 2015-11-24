#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

our ($opt_c, $opt_a) = ("","");

sub usage() {
    print 
"Usage: $0 [-c] [-a] <raw>.nex {<cooked.nex>}
where <raw.nex> is the name of a nexus file output by HybridSeq and
<cooked.nex> is the file name to write to. (If omitted, <raw>_st.nex is 
used. If input file is "-", standard input is used.) Produces an edited 
version of <raw>.nex which is acceptable  as input to SplitsTree4.
-c: Splitstree4 will use the coalescent trees. (By default, uses lineage
trees.)
-a: <cooked.nex> will contain all the blocks in <raw>.nex, but edited and
reoordered to work with Splitstree4. (By default, blocks which would 
cause Splitstree4 to issue warnings are omitted from <cooked.nex>.)
\n";
    exit 1;
}

getopts('ca');
usage() if (!@ARGV || $ARGV[0] !~ /\.nex$/ && $ARGV[0] ne "-");

my $inFile = shift;
open (IN,$inFile) or die;
my @inLines = <IN>;
close IN;

my @outLines = ();
my @delayedLines = (); # for now, these will just be deleted.
my @treesBlocks = (); # array of arrays
my @currentTreeBlock = ();

my $inDelayedBlock = 0;
my $inTreesBlock = 0;
my $foundNonComment = 0; # to correctly place comment saying this program has run.

for my $line (@inLines) {
    if (!$foundNonComment && $line !~ /^[\[#]/) {
	# Happens once only:
	push @outLines, "[Subsequently processed by $0 to make this file acceptable input to SplitsTree4]";
	$foundNonComment = 1;
    }
    $inDelayedBlock = 1 if ($line =~ /begin (filo|hybridseq|enewick)/i);
    if ($line =~ /begin (trees)/i) {
	$inTreesBlock=1;
	@currentTreeBlock=();
    }
    $line =~ s/(FORMAT .*) NOTOKENS(.*)/$1$2/;
    $line =~ s/(FORMAT .*) ITEMS=[^\s;]+(.*)/$1$2/;
    $line =~ s/(FORMAT .*) STATESFORMAT=[^\s;]+(.*)/$1$2/;
    
    if ($inTreesBlock) {
	push @currentTreeBlock, $line;
    } elsif ($inDelayedBlock) {
	push @delayedLines, $line;
    } else {
	push @outLines, $line;
    }
    if ($line =~ /end;/i) {
	if ($inTreesBlock) {
	    push @treesBlocks, [@currentTreeBlock];
	    $inTreesBlock = 0;
	}
	$inDelayedBlock = 0;
    }
}

# Default outfile name derived from input file name:
$inFile =~ s/\.nex/_st.nex/;
my $outFile = shift || $inFile;

open (OUT, ">$outFile") or die;
print OUT @outLines;
# Search @treesBlocks for the right one
my $goodTreeBlock;
my $match = $opt_c ? qr/Randomly selected coalescent/ : qr/Randomly selected lineage/;
my @badTreeBlocks;
for my $treeBlockPtr (@treesBlocks) {
    if ($treeBlockPtr->[1] =~ $match) {
	$goodTreeBlock = $treeBlockPtr;
    } else {
	push @badTreeBlocks, $treeBlockPtr;
    }
}
if ($opt_a) {
    # IF printing all: print the bad tree blocks first
    # (Splitstree4 uses the last trees block)
    for my $treeBlockPtr (@badTreeBlocks) {
	print OUT @$treeBlockPtr;
    }
}
print OUT @$goodTreeBlock;
if ($opt_a) {
    # If printing all: print the blocks unused by SplitsTree4
    # (some of which cause problems if they preceed the trees block.)
    print OUT @delayedLines;
}
close OUT;
