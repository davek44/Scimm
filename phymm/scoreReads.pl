#!/usr/bin/perl

use strict;

use Cwd;

$| = 1;

my $dataFile = shift;

# DK addition
my $ignoreFile = shift;

$dataFile =~ s/\(/\\\(/g;

$dataFile =~ s/\)/\\\)/g;

my $outputPrefix = $dataFile;

$outputPrefix =~ s/\//\_/g;

$outputPrefix =~ s/\./\_/g;

die("Usage: $0 <file containing query reads>\n") if not $dataFile;

# Set / check BLAST environment variables.

my $oldBlastDB = $ENV{'BLASTDB'};

my $oldBlastMat = $ENV{'BLASTMAT'};

die("BLAST error: environment variable 'BLASTMAT' not found; please check your local blast binary installation.\n") if not $oldBlastMat;

# Grab the list of ICMs.

opendir DOT, '.genomeData' or die("Can't open . for scanning.\n");

my @dirs = readdir DOT;

closedir DOT;

my @ICMs;

foreach my $dir ( @dirs ) {
   
   if ( -d ".genomeData/$dir" and $dir ne '.' and $dir ne '..' ) {
      
      &scanDir(".genomeData/$dir", \@ICMs);
   }
}

# Create a reverse-complement copy of the query file.

system(".scripts/revCompFASTA.pl $dataFile");

my $revDataFile = $dataFile;

$revDataFile =~ /(\.[^\.]+)$/;

my $extension = $1;

$revDataFile =~ s/$extension$/\.revComp$extension/;

my $fileCheck = $revDataFile;

$fileCheck =~ s/\\//g;

if ( not -e $fileCheck ) {
   
   die("File renaming problem [revCompFASTA.pl]: could not detect revComp file \"$fileCheck\".\n");
}


# Read organism names so we can generate results files comparable to what BLAST will output.

my $accFile = '.taxonomyData/.0_accessionMap/accessionMap.txt';

open IN, "<$accFile" or die("Can't open $accFile for reading.\n");

my $speciesDirName = {};

while ( my $line = <IN> ) {
   
   chomp $line;

   (my $orgName, my $prefix, my $seqType, my $desc) = split(/\t/, $line);
   
   $speciesDirName->{$prefix} = $orgName;
}

close IN;

# Scan for taxonomic metadata.

my $taxFile = '.taxonomyData/.3_parsedTaxData/distributionOfTaxa.txt';

open IN, "<$taxFile" or die("Can't open $taxFile for reading.\n");

my $tax = {};

while ( my $line = <IN> ) {
   
   if ( $line =~ /^\S/ ) {
      
      chomp $line;

      (my $taxType, my $taxVal, my $prefixAndSpecies, my $dirName) = split(/\t/, $line);

      if ( $taxType eq 'phylum' or $taxType eq 'class' or $taxType eq 'order' or $taxType eq 'family' or $taxType eq 'genus' ) {
	 
	 $tax->{$dirName}->{$taxType} = $taxVal;
      }
   }
}

close IN;

############################################################
# Make a list of ICMs to ignore from $ignoreFile
#
# Author: David Kelley
############################################################
open(IGNOREF, $ignoreFile);
my %ignore_icms = ();
#print "Ignoring ICMS:\n";
while(<IGNOREF>) {
    chomp;
    $ignore_icms{$_} = 1;
    #print $_, "\n";
}
############################################################


# Score the query data (using both forward and reverse directions) with the 1-periodic genome-sequence-trained ICMs.

print "Scoring reads with Phymm...";

my $score = {};

my $topScoringICM = {};

my $rawPhymmFile = 'rawPhymmOutput_' . $outputPrefix . '.txt';

$rawPhymmFile =~ s/\\//g;

open OUT, ">$rawPhymmFile" or die("Can't open $rawPhymmFile for writing.\n");

print OUT "QUERY_ID\tICM\tSCORE\n";

my $icmCount = @ICMs;

my $icmsFinished = 0;

my $logFile = $outputPrefix . '_progress.txt';

$logFile =~ s/\\//g;

open LOG, ">$logFile" or die("Can't open $logFile for writing.\n");

print LOG "Scoring reads with Phymm...\n\n";

foreach my $ICM ( sort { $a cmp $b } @ICMs ) {

   ########################################
   # Ignore ICM if its in our hash
   #
   # Author: David Kelley
   ########################################
   $ICM =~ /genomeData\/(\S+)\//;
   if(exists($ignore_icms{$1})) {
       next;
   }
   #print $1, "\n";
   ########################################
   
   my $icmPrefix = $ICM;

   $icmPrefix =~ s/.+\/([^\/]+)$/$1/;

   $icmPrefix =~ s/\.icm$//;

   my $fullScore = {};
   
   my $command = '(.scripts/.icmCode/bin/simple-score -N ' . $ICM . ' < ' . $dataFile . ' > tempFwd_' . $outputPrefix . '.txt) >& errFile_' . $outputPrefix . '.txt';
   #print $command,"\n";
   system($command);
   
   my $inFile = "tempFwd_${outputPrefix}.txt";

   $inFile =~ s/\\//g;

   open IN, "<$inFile" or die("Can't open $inFile for reading.\n");
   
   while ( my $line = <IN> ) {
      
      if ( $line =~ /(\S+)\s+(\S+)/ ) {
	 
	 my $queryID = $1;

	 my $queryScore = $2;
	 
	 $fullScore->{$queryID}->{$icmPrefix} = $queryScore;
	 
	 if ( not $score->{$queryID} or $queryScore > $score->{$queryID} ) {
	    
	    $score->{$queryID} = $queryScore;

	    $topScoringICM->{$queryID} = $icmPrefix;
	 }
      }
   }

   close IN;

   my $command = '(.scripts/.icmCode/bin/simple-score -N ' . $ICM . ' < ' . $revDataFile . ' > tempRev_' . $outputPrefix . '.txt) >& errFile_' . $outputPrefix . '.txt';
   
   system($command);
   
   $inFile = "tempRev_${outputPrefix}.txt";
   
   $inFile =~ s/\\//g;

   open IN, "<$inFile" or die("Can't open $inFile for reading.\n");
   
   while ( my $line = <IN> ) {
      
      if ( $line =~ /(\S+)\s+(\S+)/ ) {
	 
	 my $queryID = $1;

	 my $queryScore = $2;
	 
	 if ( $queryScore > $fullScore->{$queryID}->{$icmPrefix} ) {
	    
	    $fullScore->{$queryID}->{$icmPrefix} = $queryScore;
	 }
	 
	 if ( not $score->{$queryID} or $queryScore > $score->{$queryID} ) {
	    
	    $score->{$queryID} = $queryScore;

	    $topScoringICM->{$queryID} = $icmPrefix;
	 }
      }
   }

   close IN;
   
   foreach my $queryID ( sort { $a cmp $b } keys %$fullScore ) {
      
      foreach my $icmID ( keys %{$fullScore->{$queryID}} ) {
	 
	 print OUT "$queryID\t$speciesDirName->{$icmID}\t$fullScore->{$queryID}->{$icmID}\n";
      }
   }

   $icmsFinished++;

   if ( $icmsFinished % 50 == 0 ) {
      
      print LOG "   ...finished scoring read set with $icmsFinished / $icmCount IMMs...\n";
   }
}

print LOG "\n...done.\n";

close LOG;

close OUT;

# Write the Phymm-only results.

my $phymmOut = 'results.01.phymm_' . $outputPrefix . '.txt';

$phymmOut =~ s/\\//g;

open OUT, ">$phymmOut" or die("Can't open $phymmOut for writing.\n");

print OUT "QUERY_ID\tBEST_MATCH\tSCORE\tGENUS\tFAMILY\tORDER\tCLASS\tPHYLUM\n";

foreach my $queryID ( sort { $a cmp $b } keys %$score ) {
   
   my $spName = $speciesDirName->{$topScoringICM->{$queryID}};
   
   print OUT "$queryID\t$spName\t$score->{$queryID}\t$tax->{$spName}->{'genus'}\t$tax->{$spName}->{'family'}\t$tax->{$spName}->{'order'}\t$tax->{$spName}->{'class'}\t$tax->{$spName}->{'phylum'}\n";
}

close OUT;

system("rm tempFwd_${outputPrefix}.txt");

system("rm tempRev_${outputPrefix}.txt");

system("rm errFile_${outputPrefix}.txt");

print "done.\n\n";


############################################################
# Skip BLAST
#
# Author: David Kelley
############################################################
#exit 0;
############################################################

# Score the query data using BLAST.

print "Scoring reads with BLAST...";

my $newBlastDB = cwd;

$newBlastDB .= '/.blastData/';

$ENV{'BLASTDB'} = $newBlastDB;

my $blastCmd = "blastall -i $dataFile -o rawBlastOutput_${outputPrefix}.txt -d bacteriaAndArchaea -p blastn -a 2 -m 9";

system($blastCmd);

if ( -e "error.log" ) {
   
   system("rm error.log");
}

# Parse the raw BLAST results file to generate a best-hit list.

my $blastFile = 'rawBlastOutput_' . $outputPrefix . '.txt';

$blastFile =~ s/\\//g;

my $blastScore = {};

my $blastMatch = {};

open IN, "<$blastFile" or die("Can't open $blastFile for reading.\n");

while ( my $line = <IN> ) {
   
   if ( $line !~ /^#/ ) {
      
      chomp $line;

      my @fields = split(/\t/, $line);

      my $queryID = $fields[0];

      my $matchName = $fields[1];

      $matchName =~ s/\.\d+$//;

      my $currentScore = $fields[10];

      if ( not $blastScore->{$queryID} or $blastScore->{$queryID} > $currentScore ) {
	 
	 $blastScore->{$queryID} = $currentScore;

	 $blastMatch->{$queryID} = $matchName;
      }
   }
}

close IN;

my $blastOut = 'results.02.blast_' . $outputPrefix . '.txt';

$blastOut =~ s/\\//g;

open OUT, ">$blastOut" or die("Can't open $blastOut for writing.\n");

print OUT "QUERY_ID\tBEST_MATCH\tSCORE\tGENUS\tFAMILY\tORDER\tCLASS\tPHYLUM\n";

foreach my $queryID ( sort { $a cmp $b } keys %$blastScore ) {
   
   my $spName = $blastMatch->{$queryID};
   
   print OUT "$queryID\t$spName\t$blastScore->{$queryID}\t$tax->{$spName}->{'genus'}\t$tax->{$spName}->{'family'}\t$tax->{$spName}->{'order'}\t$tax->{$spName}->{'class'}\t$tax->{$spName}->{'phylum'}\n";
}

close OUT;

print "done.\n\n";

# Combine the scores to generate PhymmBL predictions.

print "Combining scores...";

# Read the raw BLAST E-values.

$blastFile = 'rawBlastOutput_' . $outputPrefix . '.txt';

$blastFile =~ s/\\//g;

$blastScore = {};

open IN, "<$blastFile" or die("Can't open $blastFile for reading.\n");

while ( my $line = <IN> ) {
   
   if ( $line !~ /^#/ ) {
      
      chomp $line;

      my @fields = split(/\t/, $line);

      my $queryID = $fields[0];

      my $matchName = $fields[1];

      $matchName =~ s/\.\d+$//;

      my $currentScore = $fields[10];

      if ( not $blastScore->{$queryID}->{$matchName} or $blastScore->{$queryID}->{$matchName} > $currentScore ) {
	 
	 $blastScore->{$queryID}->{$matchName} = $currentScore;
      }
   }
}

close IN;

# Read the raw Phymm scores.

my $phymmFile = 'rawPhymmOutput_' . $outputPrefix . '.txt';

$phymmFile =~ s/\\//g;

my $phymmScore = {};

open IN, "<$phymmFile" or die("Can't open $phymmFile for reading.\n");

while ( my $line = <IN> ) {
   
   if ( $line !~ /^QUERY/ ) {
      
      chomp $line;

      my @fields = split(/\t/, $line);

      my $queryID = $fields[0];

      my $matchName = $fields[1];

      my $currentScore = $fields[2];

      if ( not $phymmScore->{$queryID}->{$matchName} or $phymmScore->{$queryID}->{$matchName} < $currentScore ) {
	 
	 $phymmScore->{$queryID}->{$matchName} = $currentScore;
      }
   }
}

close IN;

# Iterate on each query ID.  Phymm will score everything;
# BLAST will (very occasionally) turn up no matches, so we
# use the $phymmScore keys as the reference index list.
# 
# Combine scores: the formula is [IMM score] + [1.2 * (4 - log(BLAST E-value))].

my $combinedScore = {};

my $combinedMatch = {};

foreach my $queryID ( keys %$phymmScore ) {
   
   # Check each match.

   foreach my $matchName ( keys %{$phymmScore->{$queryID}} ) {
      
      # Make sure there's a corresponding BLAST score.

      if ( not $blastScore->{$queryID} or not $blastScore->{$queryID}->{$matchName} ) {
	 
	 # If there isn't, just find the best Phymm score.

	 if ( not $combinedScore->{$queryID} ) {
	    
	    $combinedScore->{$queryID} = $phymmScore->{$queryID}->{$matchName};

	    $combinedMatch->{$queryID} = $matchName;

	 } elsif ( $combinedScore->{$queryID} < $phymmScore->{$queryID}->{$matchName} ) {
	    
	    $combinedScore->{$queryID} = $phymmScore->{$queryID}->{$matchName};

	    $combinedMatch->{$queryID} = $matchName;
	 }

      } else {
	 
	 # There is a corresponding BLAST score: compute the combined score and compare.
	 
	 my $tempScore = $phymmScore->{$queryID}->{$matchName};
	 
	 if ( $blastScore->{$queryID}->{$matchName} == 0 ) {
	    
	    $tempScore += 1000;
	    
	 } else {
	    
	    $tempScore += 1.2 * (4 - log($blastScore->{$queryID}->{$matchName}));
	 }
	 
	 # Update the best-hit record for this query if we've got
	 # something better than something we've seen so far.
	 
	 if ( not $combinedScore->{$queryID} or $tempScore > $combinedScore->{$queryID} ) {
	    
	    $combinedScore->{$queryID} = $tempScore;

	    $combinedMatch->{$queryID} = $matchName;
	 }

      } # end if ( there's a BLAST score for this query and this matching reference organism )

   } # end foreach ( reference organism to which this query was compared )

} # end foreach ( query key in the $phymmScore hash )

# Output the combined results.

my $combinedOut = 'results.03.combined_' . $outputPrefix . '.txt';

$combinedOut =~ s/\\//g;

open OUT, ">$combinedOut" or die("Can't open $combinedOut for writing.\n");

print OUT "QUERY_ID\tBEST_MATCH\tSCORE\tGENUS\tFAMILY\tORDER\tCLASS\tPHYLUM\n";

foreach my $queryID ( sort { $a cmp $b } keys %$combinedScore ) {
   
   my $matchName = $combinedMatch->{$queryID};
   
   print OUT "$queryID\t$matchName\t$combinedScore->{$queryID}\t$tax->{$matchName}->{'genus'}\t$tax->{$matchName}->{'family'}\t$tax->{$matchName}->{'order'}\t$tax->{$matchName}->{'class'}\t$tax->{$matchName}->{'phylum'}\n";
}

close OUT;

print "done.\n\n";


############################################################
# No plasmids!  Only train on largest fasta file in each
# strain directory
#
# Author: David Kelley
############################################################
sub scanDir {
   
   my $dir = shift;

   my $arrayRef = shift;
   
   opendir DOT, $dir or die("Can't open $dir for scanning.\n");

   my @files = readdir DOT;

   closedir DOT;

   # original
   #foreach my $file ( @files ) {      
   #   if ( $file =~ /\.icm$/ ) {
   #	 push @$arrayRef, "$dir/$file";
   #   }
   #}

   # DK
   my $filesize;
   my $maxsize = 0;
   my $maxfna;
   foreach my $file (@files) {
       if($file =~ /.fna$/) {
	   $filesize = -s "$dir/$file";
	   if($filesize > $maxsize) {
	       $maxsize = $filesize;
	       $maxfna = $file
	   }	   
       }
   }

   my $maxicm = substr($maxfna, 0, -3) . "icm";
   push @$arrayRef, "$dir/$maxicm";
}

############################################################
