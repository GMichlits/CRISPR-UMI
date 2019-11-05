#!/usr/bin/perl
use strict;

########################################################
# Correct the exon number of the GTF ((-)-strand) and sort
# IMP/IMBA Bioinformatics department
########################################################

#get exon count of transcripts
open(REF, "refGene.txt");
my $exonC = {};
while (<REF>)
{
    my @f = split "\t";
    $exonC->{$f[1]} = $f[8];
}
close REF;


my $out = {};
while (my $l=<>)
{
    next if ($l=~/start_codon/);
    next if ($l=~/stop_codon/);

    my @f = split "\t", $l;
    if ($l=~/transcript_id\s+"(.*?)";/)
    {
	my $tid = $1;
	my $ec = $exonC->{$tid} + 1;
	if ($l=~/exon_number\s+"(.*?)";/)
	{
	    my $en = $1;

	    #reverse exon number for "-" strand
	    if ($f[6] eq "-")
	    {
		$en = $ec - $en;
		$l =~ s/exon_number\s+".*?";/exon_number "$en";/;
	    }
	    if ($f[2] eq "exon")
	    {
		$out->{$f[0]}->{$tid}->{$en}->{"exon"} = $l;
	    } elsif ($f[2] eq "CDS")
	    {
		$out->{$f[0]}->{$tid}->{$en}->{"CDS"} = $l;
	    }
	}
    }
}

my $gNum = 0;
foreach my $chr (sort keys %$out)
{
    $gNum += scalar keys %{$out->{$chr}};
    foreach my $id (sort keys %{$out->{$chr}})
    {
	#sort output by exon number
	foreach my $en (sort { $a <=> $b } keys %{$out->{$chr}->{$id}})
	{
	    print $out->{$chr}->{$id}->{$en}->{exon};
	    print $out->{$chr}->{$id}->{$en}->{CDS};
	}
    }
}

print STDERR "Number of transcripts in GTF: $gNum\n";
