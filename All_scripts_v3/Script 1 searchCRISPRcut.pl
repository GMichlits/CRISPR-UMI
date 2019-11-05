#!/usr/bin/perl
use strict;
use warnings;

########################################################
# Search all CRISPR site (NGG)
# IMP/IMBA Bioinformatics department
########################################################


#function to search all cut sites
sub searchSeq
{
    my ($chr, $seq, $pam, $strand) = @_;
    my $offset = 0;

    my $result = index($seq, $pam, $offset);

#cut is between 16./17. or 17./18. nucleotide of the gDNA. Take 15-20

#0b           34567890
#XXXXXXXXXXXXXXXCCCXXNGG
#1b                 90

#0b  0123
#     CCNXXCCCXXXXXXXXXXXXXXX
#1b  0123456789

    while ($result != -1) 
    {
	if ($strand eq "+") 
	{
	    if (($result -7) >= 0)
	    {
		print "$chr\t". ($result - 7 ) . "\t". ($result - 1) . "\t1\t1\t$strand\n";
	    }
	} else {
	    print "$chr\t". ($result + 3 ) . "\t". ($result + 9) . "\t1\t1\t$strand\n";
	}


	$offset = $result + 1;
	$result = index($seq, $pam, $offset);
    }    
}

#CRISPR defining PAM
my $PAM = "GG";
my $PAM_rev = "CC";

#search in mm10 genome
my $genomeFile = "mm10.fa";

open(GENOME, "$genomeFile");

my $chr = <GENOME>;
$chr =~ s/\s.*//;
$chr =~ s/>//;

#iterate through chromosomes
my $seq = "";
while (<GENOME>)
{
    chomp;

    if (/^>/)
    {
	searchSeq($chr, $seq, $PAM, "+");
	searchSeq($chr, $seq, $PAM_rev, "-");
	$chr = $_;
	$chr =~ s/\s.*//;
	$chr =~ s/>//;

	$seq = "";
	next;
    }

    $seq .= $_;
}
    
searchSeq($chr, $seq, $PAM, "+");
searchSeq($chr, $seq, $PAM_rev, "-");
	
