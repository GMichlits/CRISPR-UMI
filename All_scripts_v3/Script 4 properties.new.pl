#!/usr/bin/perl

########################################################
# Extract and combine properties
# IMP/IMBA Bioinformatics department
########################################################


use List::Util qw(min max);
use strict;

#read transcript to gene relationship
open(REF, "refGene.txt");
my $numTrans={};
while(<REF>)
{
    next if (! /\tNM_/);
    my @c = split "\t";
    $numTrans->{$c[12]}->{$c[1]}=1;
}
close(REF);

#read GTF properties
open(GTF, "tmp/mm10_refSeq.cdsMod3Len.gtf");
my $transToExonLength = {};
while (my $line = <GTF>)
{
    chomp $line;
    my @c = split "\t";
    my $tid = $1 if ($line =~ /transcript_id "(.*?)"/);
    my $eid = $1 if ($line =~ /exon_number "(.*?)"/);
    my $l = $1 if ($line =~ /length "(.*?)"/);
    $transToExonLength->{$tid}->{$eid} = $l;
}
close GTF;

#calculate transcript length and cumulative length of exon
foreach my $tid (keys %$transToExonLength)
{
    my $length=0;
    foreach my $eid (sort {$a <=> $b} %{$transToExonLength->{$tid}})
    {
	$transToExonLength->{$tid}->{$eid} += $length; #cummulative length
	$length = $transToExonLength->{$tid}->{$eid};
    }
    $transToExonLength->{$tid}->{whole} = $length;
}

#read 2nd ATG results
open(ATG, "tmp/ATG2nd.txt");
my $atg={};
while (<ATG>)
{
    chomp;
    my @c = split "\t";

    $atg->{$c[0]} = $c[1];
}
close ATG;

#read pfam results
open(PFAMSCAN, "pfamscan.txt");
my $pfam={};
while (<PFAMSCAN>)
{
    next if (! /^NM_/);

    chomp;
    s/\s+/\t/g;
    my @c = split "\t";
    $pfam->{$c[0]}->{$c[3]}->{$c[4]} = $c[5] . "|" . $c[7];
}
close PFAMSCAN;

#go through intersection CDS vs cut region
my $gene={};
while (my $line=<>)
{
    chomp $line;
    my @c = split "\t", $line;

    #extend to (4N-gDNA(20)-NGG-3N)
    my $start = 0;
    my $stop = 0;
    if ($c[5] eq "-")
    {
	$start = $c[1] - 6;
	$stop = $c[2] + 18;
    } else {
	$start = $c[1] - 18;
	$stop = $c[2] + 6;
    }

    next if ($start < 0);
    my $region = "$c[0]:$start-$stop($c[5])";

    my $gid = $line;
    $gid =~ s/.*gene_id "(.*?)".*/\1/;

    #CDS position
    my $tid = $1 if ($line =~ /transcript_id "(.*?)"/);
    my $eid = $1 if ($line =~ /exon_number "(.*?)"/);
    $eid -= 1;
    my $cutLen = 0;
    if ($c[12] eq "-") #transcript orientation
    {
	$cutLen = $c[10] - ($c[2] - 1) + $transToExonLength->{$tid}->{$eid};
    } else {
	$cutLen = ($c[1] + 1) - ($c[9] - 1) + $transToExonLength->{$tid}->{$eid};
    }

    #distance penalty
    my $d = 0.5 + 0.5*(1-($cutLen / $transToExonLength->{$tid}->{whole}));
    if ((! exists $gene->{$gid}->{$region}->{"relDist"}) || ($gene->{$gid}->{$region}->{"relDist"} > $d))
    {
	$gene->{$gid}->{$region}->{"relDist"} = $d;
	$gene->{$gid}->{$region}->{"dist"} = $cutLen;
	$gene->{$gid}->{$region}->{"transLen"} = $transToExonLength->{$tid}->{whole};
    } 
    #END - CDS position

    #pfam
    my $protPos = int($cutLen / 3);
    
    if (exists $pfam->{$tid})
    {
	foreach my $start (sort {$a <=> $b} keys %{$pfam->{$tid}})
	{
	    next if ($start > $protPos);
	    foreach my $end (sort {$a <=> $b} keys %{$pfam->{$tid}->{$start}})
	    {
		next if ($end < $protPos);
		push @{$gene->{$gid}->{$region}->{"pfam"}}, $pfam->{$tid}->{$start}->{$end};
	    }
	}
    }
    #END - pfam

    #ATG
    if (exists $atg->{$tid})
    {
	if ($atg->{$tid} > $cutLen)
	{
	    push @{$gene->{$gid}->{$region}->{"ATG"}}, $atg->{$tid};
	}
    }
    #END - ATG

    #further properties
    my $tid = $line;
    $tid =~ s/.*transcript_id "(.*?)".*/\1/;

    my $eid = $line;
    $eid =~ s/.*exon_id "(.*?)".*/\1/;

    my $mod3 = $line;
    $mod3 =~ s/.*mod3 "(.*?)".*/\1/;

    my $l = $line;
    $l =~ s/.*length "(.*?)".*/\1/;

    $gene->{$gid}->{$region}->{"tid"}->{$tid} = 1;
    $gene->{$gid}->{$region}->{"eid"}->{$eid} = 1;

    if (! exists $gene->{$gid}->{$region}->{"mod3"})
    {
	$gene->{$gid}->{$region}->{"mod3"} = $mod3;
    } elsif ($mod3 == 0)
    {
	$gene->{$gid}->{$region}->{"mod3"} = $mod3;
    }

    if (! exists $gene->{$gid}->{$region}->{"length"})
    {
	$gene->{$gid}->{$region}->{"length"} = $l;
    } elsif ($gene->{$gid}->{$region}->{"length"} < $l)
    {
	$gene->{$gid}->{$region}->{"length"} = $l;
    }
}

#sort and write
foreach my $gid (sort keys %$gene)
{
    foreach my $region (sort keys %{$gene->{$gid}})
    {
	print "$gid\t$region\t$gene->{$gid}->{$region}->{length}\t";
	
	if ($gene->{$gid}->{$region}->{mod3} == 0)
	{
	    print "TRUE\t";
	} else {
	    print "FALSE\t";
	}
	
	#common region
	my $maxTransNum = scalar keys %{$numTrans->{$gid}};
	my $transNum = scalar keys %{$gene->{$gid}->{$region}->{"tid"}};
	if ($transNum == $maxTransNum)
	{
	    print "TRUE\t";
	} else {
	    print "FALSE\t";
	}

	my $exon = join "|", sort keys $gene->{$gid}->{$region}->{"eid"};
	print "$exon\t";

	if (exists $gene->{$gid}->{$region}->{"pfam"})
	{
	    my $pfam = join "|", sort keys { map { $_ => 1 } @{$gene->{$gid}->{$region}->{"pfam"}} };
	    print "$pfam\t";
	} else {
	    print "\t";
	}

	print $gene->{$gid}->{$region}->{"relDist"} . "\t";
	print $gene->{$gid}->{$region}->{"dist"} . "\t";
	print $gene->{$gid}->{$region}->{"transLen"} . "\t";

	if (exists $gene->{$gid}->{$region}->{"ATG"})
	{
	    my $atg = join "|", @{$gene->{$gid}->{$region}->{"ATG"}};
	    print "$atg\n";
	} else {
	    print "\n";
	}

    }
}
