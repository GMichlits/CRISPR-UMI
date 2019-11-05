#!/usr/bin/perl
use strict;
use List::Util qw(min max);

########################################################
# Generate final gDNA table with all properties an filters
# IMP/IMBA Bioinformatics department
########################################################

#create alghtough file with subset
my $TOP=4;

#read and filter gDNA by restriction site, length, none ATGC letters
open(SCORE, "cut.guide30.score.txt");
my $guide = {};
while (<SCORE>)
{
    chomp;
    my @c = split " ";

    my $seq = $c[1];


    if ($seq =~ /[UWSMKRYBDHVN]/)
    {
	print STDERR "WARNING: excluded $c[0] ($seq)\n";
	next;
    }

    if (length($seq) != 30)
    {
	print STDERR "Warning: sequence not 30NT long! Excluded $c[0] ($seq)\n";
	next
    }
    my $gs = substr($seq, 4, 20);
    if ( $gs =~ /(GAAGAC)|(GTCTTC)|(CTCGAG)|(CGTCTC)|(GAGACG)|(^AAGAC)|(CTCGA$)/i ) #xho (CTCGAG); BsmBI (CGTCTC)|(GAGACG)
    {
	print STDERR "Warning: sequence contains BbsI/Xho/BsmBI site in guide $gs! Excluded $c[0] ($seq) due to $1$2$3$4$5$6$7\n";
	next
    }
    $guide->{$c[0]}->{seq} = $c[1];
    $guide->{$c[0]}->{score} = $c[2];
}   
close SCORE;

#print scalar(keys(%$guide)) . "\n";

#Mark and penalize score if offtarges with <=2 mismatches exist
open(NM1, "tmp/guide.2MM.NM2count.txt");
while (<NM1>)
{
    chomp;
    if (/(\d+)\s+>(\S+)/) # ">" from FASTA
    {
	my $pos = $2;
	my $cnt = $1;
	next if (! exists($guide->{$pos}));
	if ($cnt > 1)
	{
	    $guide->{$pos}->{score} *= 0.2; #offtarget 2MM
	    $guide->{$pos}->{MM2} = "TRUE";
	}  else {
	    $guide->{$pos}->{MM2} = "FALSE";
	}
    }
}
close NM1;


#read further properties
open(PROP, "tmp/cut.propertie.new.txt");
my $sort = {};
while (<PROP>)
{
    chomp;
    my @c = split "\t";

    next if (! exists($guide->{$c[1]}));


    push @{$sort->{$c[0]}->{$guide->{$c[1]}->{score}}}, $c[1]; #take score to sort & filter

    $guide->{$c[1]}->{gene}->{$c[0]}->{length} = $c[2];
    $guide->{$c[1]}->{gene}->{$c[0]}->{mod3} = $c[3];
    $guide->{$c[1]}->{gene}->{$c[0]}->{commonCDS} = $c[4]; #commonCDS
    $guide->{$c[1]}->{gene}->{$c[0]}->{exon} = $c[5];
    $guide->{$c[1]}->{gene}->{$c[0]}->{pfam} = $c[6];
    $guide->{$c[1]}->{gene}->{$c[0]}->{relDist} = $c[7];
    $guide->{$c[1]}->{gene}->{$c[0]}->{dist} = $c[8];
    $guide->{$c[1]}->{gene}->{$c[0]}->{transLen} = $c[9];
    $guide->{$c[1]}->{gene}->{$c[0]}->{ATG} = $c[10];
}
close PROP;

#read offtargets if exist (2nd call of this script)
open(OT, "tmp/offtarget.txt");
my $ot = {};
my $offtargetRead=0;
while (<OT>)
{
    chomp;
    my @c = split "\t";
    $ot->{$c[0]}->{OT95} = $c[1];
    $ot->{$c[0]}->{OT90} = $c[2];
    $ot->{$c[0]}->{OT80} = $c[3];
    $offtargetRead=1;
}
close OT;

#metascore calculation
my $metascore={};
foreach my $gene (keys %$sort)
{
    my $exon = {};
    foreach my $score (sort {$b <=> $a} keys %{$sort->{$gene}})
    {
	foreach my $pos (sort @{$sort->{$gene}->{$score}})
	{
	    #pfam
	    if ($guide->{$pos}->{gene}->{$gene}->{pfam} ne "")
	    {
		if ($guide->{$pos}->{gene}->{$gene}->{pfam} =~ /Domain/)
		{
		    $guide->{$pos}->{gene}->{$gene}->{metascore} = $guide->{$pos}->{score};
		} elsif ($guide->{$pos}->{gene}->{$gene}->{pfam} =~ /Family/)
		{
		    $guide->{$pos}->{gene}->{$gene}->{metascore} = $guide->{$pos}->{score} * 0.9;
		} else {
		    $guide->{$pos}->{gene}->{$gene}->{metascore} = $guide->{$pos}->{score} * 0.8;
		}
	    } else {
		$guide->{$pos}->{gene}->{$gene}->{metascore} = $guide->{$pos}->{score} * 0.5;

		#mod3 & length
		if ($guide->{$pos}->{gene}->{$gene}->{mod3} eq "TRUE")
		{
		    if ($guide->{$pos}->{gene}->{$gene}->{length} <= 100)
		    {
			$guide->{$pos}->{gene}->{$gene}->{metascore} *= 0.2;
		    } elsif ($guide->{$pos}->{gene}->{$gene}->{length} <= 300) {
			$guide->{$pos}->{gene}->{$gene}->{metascore} *= 0.4;
		    } elsif ($guide->{$pos}->{gene}->{$gene}->{length} <= 500) {
			$guide->{$pos}->{gene}->{$gene}->{metascore} *= 0.7;
		    } else {
			$guide->{$pos}->{gene}->{$gene}->{metascore} *= 0.9;
		    }
		} else {
		    if ($guide->{$pos}->{gene}->{$gene}->{length} <= 100)
		    {
			$guide->{$pos}->{gene}->{$gene}->{metascore} *= 0.4;
		    } elsif ($guide->{$pos}->{gene}->{$gene}->{length} <= 300) {
			$guide->{$pos}->{gene}->{$gene}->{metascore} *= 0.6;
		    } elsif ($guide->{$pos}->{gene}->{$gene}->{length} <= 500) {
			$guide->{$pos}->{gene}->{$gene}->{metascore} *= 0.8;
		    } else {
			$guide->{$pos}->{gene}->{$gene}->{metascore} *= 1;
		    }
		}
		
	    }
	    
	    #exon usage
	    $exon->{$guide->{$pos}->{gene}->{$gene}->{exon}} += 1;
	    $guide->{$pos}->{gene}->{$gene}->{exonUsage} = $exon->{$guide->{$pos}->{gene}->{$gene}->{exon}};
	    if ($exon->{$guide->{$pos}->{gene}->{$gene}->{exon}} >= 4)
	    {
		$guide->{$pos}->{gene}->{$gene}->{metascore} *= 0.3;
	    } elsif ($exon->{$guide->{$pos}->{gene}->{$gene}->{exon}} == 3) {
		$guide->{$pos}->{gene}->{$gene}->{metascore} *= 0.5;
	    } elsif ($exon->{$guide->{$pos}->{gene}->{$gene}->{exon}} == 2) {
		$guide->{$pos}->{gene}->{$gene}->{metascore} *= 0.7;
	    }

	    #commonCDS
	    if ($guide->{$pos}->{gene}->{$gene}->{commonCDS} eq "FALSE")
	    {
		$guide->{$pos}->{gene}->{$gene}->{metascore} *= 0.2;
	    }
	    
	    #transcript distance
	    $guide->{$pos}->{gene}->{$gene}->{metascore} *= $guide->{$pos}->{gene}->{$gene}->{relDist};

	    #ATG
	    if ($guide->{$pos}->{gene}->{$gene}->{ATG} ne "")
	    {
		$guide->{$pos}->{gene}->{$gene}->{metascore} *= 0.3;
	    }


	    if (exists $ot->{$pos})
	    {
	        #OT80
		if ($ot->{$pos}->{OT80} > 50)
		{
		    $guide->{$pos}->{gene}->{$gene}->{metascore} *= 0.8;
		}

	        #OT90
		if ($ot->{$pos}->{OT90} > 0)
		{
		    $guide->{$pos}->{gene}->{$gene}->{metascore} *= 0.1;
		}
	    } elsif ($offtargetRead == 1) {
		$guide->{$pos}->{gene}->{$gene}->{metascore} = 0; #if no offtarget was calculated for that gene set metascore 0!
	    }

	    push @{$metascore->{$gene}->{$guide->{$pos}->{gene}->{$gene}->{metascore}}}, $pos;
	    
	}
#	last if ($cnt > 50);
    }
}


#write table
print "genename\tpos\tgeneOverlap\tseq [4N-gDNA(20)-NGG-3N]\tscore\tmetascore\tdistanceScore\tmin CDS distance\tCDS length\texon length\tmod3\texon\texonUsage\tpfam\tATG\tMM2\tcommonCDS\tOT95\tOT90\tOT80\n";

open(SELECTED, ">mergeData.2MM.new.top" . $TOP . ".txt");
print SELECTED "genename\tpos\tgeneOverlap\tseq [4N-gDNA(20)-NGG-3N]\tscore\tmetascore\tdistanceScore\tmin CDS distance\tCDS length\texon length\tmod3\texon\texonUsage\tpfam\tATG\tMM2\tcommonCDS\tdistanceOK\tOT95\tOT90\tOT80\n";

system("echo track type=bedGraph name=score > mergeData.2MM.new.score.bedgraph");
system("echo track type=bedGraph name=metascore > mergeData.2MM.new.metascore.bedgraph");
open(SCORE, "| sort -k1,1 -k2,2n -k3,3n -u >> mergeData.2MM.new.score.bedgraph");
open(METASCORE, "| sort -k1,1 -k2,2n -k3,3n -u >> mergeData.2MM.new.metascore.bedgraph");

system("echo track type=bedGraph name=metascoreTop" . $TOP . " > mergeData.2MM.new.metascoreTop" . $TOP . ".bedgraph");
open(METASCORE4, "| sort -k1,1 -k2,2n -k3,3n -u >> mergeData.2MM.new.metascoreTop" . $TOP . ".bedgraph");

my $minDist=5;
my $selectedStart={};
foreach my $gene (sort keys %$metascore)
{
    my $cnt = 1;
    my $cntPos = 1;
    my @distanceExclPos = ();
    foreach my $s (sort {$b <=> $a} keys %{$metascore->{$gene}})
    {
	foreach my $pos (sort @{$metascore->{$gene}->{$s}})
	{
	    last if ($cntPos > 50);
	    $cntPos++;

	    print "$gene\t$pos\t" . 
		join("|", keys(%{$guide->{$pos}->{gene}})) . "\t" .
		$guide->{$pos}->{seq} . "\t" .
		$guide->{$pos}->{score} . "\t" .
		$guide->{$pos}->{gene}->{$gene}->{metascore} . "\t" .
		$guide->{$pos}->{gene}->{$gene}->{relDist} . "\t" .
		$guide->{$pos}->{gene}->{$gene}->{dist} . "\t" .
		$guide->{$pos}->{gene}->{$gene}->{transLen} . "\t" .
		$guide->{$pos}->{gene}->{$gene}->{length} . "\t" . 
		$guide->{$pos}->{gene}->{$gene}->{mod3} . "\t" . 
		$guide->{$pos}->{gene}->{$gene}->{exon} . "\t" . 
		$guide->{$pos}->{gene}->{$gene}->{exonUsage} . "\t" . 
		$guide->{$pos}->{gene}->{$gene}->{pfam} . "\t" .
		$guide->{$pos}->{gene}->{$gene}->{ATG} . "\t" .
		$guide->{$pos}->{MM2} . "\t" .
		$guide->{$pos}->{gene}->{$gene}->{commonCDS} . "\t";

	    if (exists $ot->{$pos})
	    {
		print $ot->{$pos}->{OT95} . "\t" .
		    $ot->{$pos}->{OT90} . "\t" .
		    $ot->{$pos}->{OT80} . "\n";
	    } else {
		print "NA\tNA\tNA\n";
	    }

	    my $distOk=1;
	    if ($pos =~ /(^.*?):(\d+)-(\d+)\((.)\)/)
	    {
		my $chr=$1;
		my $strand=$4;
		my $leftmost=$2;

		if (exists $selectedStart->{$gene}->{$chr}->{$strand})
		{
		    foreach my $l (@{$selectedStart->{$gene}->{$chr}->{$strand}})
		    {
			if (abs($leftmost-$l) < $minDist)
			{
			    push @distanceExclPos, $pos;
			    $distOk=0;
			    last;
			}
		    }
		}
	    
		if ($cnt <= $TOP && $distOk)
		{
		    push @{$selectedStart->{$gene}->{$chr}->{$strand}}, $leftmost;

		    print SELECTED "$gene\t$pos\t" . 
			join("|", keys(%{$guide->{$pos}->{gene}})) . "\t" .
			$guide->{$pos}->{seq} . "\t" .
			$guide->{$pos}->{score} . "\t" .
			$guide->{$pos}->{gene}->{$gene}->{metascore} . "\t" .
			$guide->{$pos}->{gene}->{$gene}->{relDist} . "\t" .
			$guide->{$pos}->{gene}->{$gene}->{dist} . "\t" .
			$guide->{$pos}->{gene}->{$gene}->{transLen} . "\t" .
			$guide->{$pos}->{gene}->{$gene}->{length} . "\t" . 
			$guide->{$pos}->{gene}->{$gene}->{mod3} . "\t" . 
			$guide->{$pos}->{gene}->{$gene}->{exon} . "\t" . 
			$guide->{$pos}->{gene}->{$gene}->{exonUsage} . "\t" . 
			$guide->{$pos}->{gene}->{$gene}->{pfam} . "\t" .
			$guide->{$pos}->{gene}->{$gene}->{ATG} . "\t" .
			$guide->{$pos}->{MM2} . "\t" .
			$guide->{$pos}->{gene}->{$gene}->{commonCDS} . "\t" .
			"TRUE\t"; #distance OK
		    
		    
		    if (exists $ot->{$pos})
		    {
			print SELECTED $ot->{$pos}->{OT95} . "\t" .
			    $ot->{$pos}->{OT90} . "\t" .
			    $ot->{$pos}->{OT80} . "\n";
		    } else {
			print SELECTED "NA\tNA\tNA\n";
		    }
		}
	    }

	    if ($pos =~ /(^.*?):(\d+)-(\d+)\((.)\)/)
	    {
		my $chr = $1;
		if ($4 eq "-")
		{
		    my $s = ($3-26);
		    my $e = ($2+4);
		    print SCORE "$chr\t$s\t$e\t" . $guide->{$pos}->{score} . "\n";
		    print METASCORE "$chr\t$s\t$e\t" . $guide->{$pos}->{gene}->{$gene}->{metascore} . "\n";
		    if ($cnt <= $TOP && $distOk)
		    {
			print METASCORE4 "$chr\t$s\t$e\t" . $guide->{$pos}->{gene}->{$gene}->{metascore} . "\n";
			$cnt++;
		    }


		} else {
		    my $s = ($2+25);
		    my $e = ($3-5);
		    print SCORE "$chr\t$s\t$e\t" . $guide->{$pos}->{score} . "\n";
		    print METASCORE "$chr\t$s\t$e\t" . $guide->{$pos}->{gene}->{$gene}->{metascore} . "\n";
		    if ($cnt <= $TOP && $distOk)
		    {
			print METASCORE4 "$chr\t$s\t$e\t" . $guide->{$pos}->{gene}->{$gene}->{metascore} . "\n";
			$cnt++;
		    }
		}
	    }
	}
	last if ($cntPos > 50);
    }

    if ($cnt <= $TOP)
    {
	foreach my $pos (@distanceExclPos)
	{
	    next if ($cnt > $TOP);
	    $cnt++;

	    print SELECTED "$gene\t$pos\t" . 
		join("|", keys(%{$guide->{$pos}->{gene}})) . "\t" .
		$guide->{$pos}->{seq} . "\t" .
		$guide->{$pos}->{score} . "\t" .
		$guide->{$pos}->{gene}->{$gene}->{metascore} . "\t" .
		$guide->{$pos}->{gene}->{$gene}->{relDist} . "\t" .
		$guide->{$pos}->{gene}->{$gene}->{dist} . "\t" .
		$guide->{$pos}->{gene}->{$gene}->{transLen} . "\t" .
		$guide->{$pos}->{gene}->{$gene}->{length} . "\t" . 
		$guide->{$pos}->{gene}->{$gene}->{mod3} . "\t" . 
		$guide->{$pos}->{gene}->{$gene}->{exon} . "\t" . 
		$guide->{$pos}->{gene}->{$gene}->{exonUsage} . "\t" . 
		$guide->{$pos}->{gene}->{$gene}->{pfam} . "\t" .
		$guide->{$pos}->{gene}->{$gene}->{ATG} . "\t" .
		$guide->{$pos}->{MM2} . "\t" .
		$guide->{$pos}->{gene}->{$gene}->{commonCDS} . "\t" .
		"FALSE\t"; #distance OK
	    if (exists $ot->{$pos})
	    {
		print SELECTED $ot->{$pos}->{OT95} . "\t" .
		    $ot->{$pos}->{OT90} . "\t" .
		    $ot->{$pos}->{OT80} . "\n";
	    } else {
		print SELECTED "NA\tNA\tNA\n";
	    }

	    if ($pos =~ /(^.*?):(\d+)-(\d+)\((.)\)/)
	    {
		my $chr = $1;
		if ($4 eq "-")
		{
		    my $s = ($3-26);
		    my $e = ($2+4);
		    print METASCORE4 "$chr\t$s\t$e\t" . $guide->{$pos}->{gene}->{$gene}->{metascore} . "\n";
		} else {
		    my $s = ($2+25);
		    my $e = ($3-5);
		    print METASCORE4 "$chr\t$s\t$e\t" . $guide->{$pos}->{gene}->{$gene}->{metascore} . "\n";
		}
	    }
	}
    }
}

close SCORE;
close METASCORE;
close METASCORE4;
close SELECTED;
