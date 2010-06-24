#!/usr/bin/perl

use strict;

$| = 1;

##########################################################################
# 
# User interactivity prologue.  Figure out what the user wants to do, within the confines of what we're allowing.
# 
##########################################################################

my $fileLoc = '';

my @fileLocs = ();

# Public service announcement.

print "\nPLEASE NOTE -- IMPORTANT: This script can be used to build genome models\n";

print "for only one organism at a time.  You can specify multiple fasta files\n";

print "belonging to the same organism, and you'll be asked later whether or not you\n";

print "want to create one model per file or create one model combining all of the\n";

print "specified files, but either way, all models constructed in this way will\n";

print "be associated with the same organism.\n\n";

print "If you want to construct models for multiple organisms, run this script\n";

print "repeatedly, once per organism.\n\n";

print "------------------------------------------------------------------------\n\n";

##########################################################################
# 
# Get taxonomic data.
# 
##########################################################################

my $tax = {};

print "Please input the taxonomic data for the organism represented by your input files.\n";

print "Every level except genus and species may be left blank (if there's no existing\n";

print "classification for that level).\n\n";

print "Phylum: ";

my $response = <STDIN>;

chomp $response;

$tax->{'phylum'} = $response;

print "Class: ";

$response = <STDIN>;

chomp $response;

$tax->{'class'} = $response;

print "Order: ";

$response = <STDIN>;

chomp $response;

$tax->{'order'} = $response;

print "Family: ";

$response = <STDIN>;

chomp $response;

$tax->{'family'} = $response;

print "Genus [required]: ";

$response = <STDIN>;

chomp $response;

while ( length($response) == 0 ) {
   
   print "Genus [required]: ";

   $response = <STDIN>;

   chomp $response;
}

$tax->{'genus'} = $response;

print "Species (not the full binomial name with the genus, but just the species label) [required]: ";

$response = <STDIN>;

chomp $response;

while ( length($response) == 0 ) {
   
   print "Species (not the full binomial name with the genus, but just the species label) [required]: ";

   $response = <STDIN>;

   chomp $response;
}

$tax->{'species'} = $response;

print "Strain: ";

$response = <STDIN>;

chomp $response;

$tax->{'strain'} = $response;

print "\n";

my $displayBinomialAndStrain = "$tax->{'genus'} $tax->{'species'}";

if ( length($tax->{'strain'}) > 0 ) {
   
   $displayBinomialAndStrain .= " $tax->{'strain'}";
}

##########################################################################
# 
# Get input file information.
# 
##########################################################################

print "Adding custom genome for \"$displayBinomialAndStrain.\"  Genome data will be built from fasta files provided by you.\n\nHave you got one or more files? [Enter S for a single file (default), M for more]: ";

$response = <STDIN>;

chomp $response;

if ( length($response) > 1 or ( length($response) == 1 and $response !~ /[sSmM]/ ) ) {
   
   print "\n   Please select from among the options given.\n\n";

   while ( length($response) > 1 or ( length($response) == 1 and $response !~ /[sSmM]/ ) ) {
      
      print "[Enter S for a single file (default), M for more]: ";

      $response = <STDIN>;

      chomp $response;
   }
}

if ( length($response) == 0 or $response =~ /[sS]/ ) {
   
   # Just one input file.  Ask for its location and scan it to make sure it's a valid fasta file.

   my $foundCorrectFileType = 0;

   while ( $foundCorrectFileType != 1 ) {
      
      print "\nPlease enter the location of the fasta file you wish to add to the database: ";

      $fileLoc = <STDIN>;
      
      chomp $fileLoc;

      if ( not -e $fileLoc ) {
	 
	 print "\n   The file you specified doesn't appear to exist.\n";

      } elsif ( not -f $fileLoc or -B $fileLoc ) {
	 
	 print "\n   The location you entered isn't a regular file.  Please provide a fasta file.\n";

      } else {
	 
	 print "\n   Checking file format...";

	 $foundCorrectFileType = &checkFileFormat($fileLoc);

	 print "done.  File complies with standard nucleotide fasta spec.\n\n";
      }
   }

} elsif ( $response =~ /[mM]/ ) {
   
   # The user has multiple fasta input files.
   # 
   # Figure out whether they just want to provide a directory location which will then be scanned for fasta files,
   # or they want to enter file locations by hand one at a time.

   print "\nYou can specify a directory to scan for fasta files contained within it (all files\n";
   
   print "with the extensions \".fasta\", \".fa\" and \".fna\" [not case-sensitive] will be used;\n";
   
   print "all others in the specified directory will be ignored), or you can specify the files one at a time.\n\n";

   print "Please select D for a directory scan (default) or M for manual entry of fasta file locations, one at a time: ";

   $response = <STDIN>;

   chomp $response;
   
   if ( length($response) > 1 or ( length($response) == 1 and $response !~ /[dDmM]/ ) ) {
      
      print "\n   Please select from among the options given.\n\n";

      while ( length($response) > 1 or ( length($response) == 1 and $response !~ /[dDmM]/ ) ) {

	 print "[Enter D for directory scan (default), M for manual file entry]: ";

	 $response = <STDIN>;

	 chomp $response;
      }
   }

   if ( length($response) == 0 or $response =~ /[dD]/ ) {
      
      # User wants to specify a directory to scan.  Ask for its location, check for the presence of at least one fasta file,
      # and scan all contained fasta files for correct formatting.
      
      my $foundGoodDir = 0;

      while ( not $foundGoodDir ) {
	 
	 print "\nPlease enter the location of the directory to scan: ";

	 my $dirLoc = <STDIN>;

	 chomp $dirLoc;

	 if ( not -e $dirLoc ) {
	    
	    print "\n   The directory you specified doesn't appear to exist.\n";

	 } elsif ( not -d $dirLoc ) {
	    
	    print "\n   The location you specified isn't a directory.\n";

	 } else {
	    
	    print "\n   Directory located.  Scanning for fasta files...\n";

	    opendir DOT, $dirLoc or die("\nCan't open $dirLoc for scanning.  Aborting.\n\n");

	    my @files = readdir DOT;

	    closedir DOT;
	    
	    my $foundFile = 0;

	    foreach my $file ( @files ) {
	       
	       if ( $file =~ /\.fna$/i or $file =~ /\.fa$/i or $file =~ /\.fasta$/i ) {
		  
		  $foundFile = 1;
		  
		  print "\n      Found file $file\.  Checking file format...";

		  my $fileCheckResult = &checkFileFormat("$dirLoc/$file");

		  print "done.  File complies withs standard nucleotide fasta spec.\n";

		  push @fileLocs, "$dirLoc/$file";
	       }
	    }

	    if ( not $foundFile ) {
	       
	       die("\n   ...done.  No fasta files (extensions .fna, .fa, or .fasta) were found in the specified directory.  Aborting.\n\n");

	    } else {
	       
	       $foundGoodDir = 1;
	       
	       my $fileCount = $#fileLocs + 1;

	       my $plural = 's';

	       if ( $fileCount == 1 ) {
		  
		  $plural = '';
	       }

	       print "\n   ...done.  $fileCount fasta file$plural found and verified as valid.\n\n";
	    }
	 }
      }

   } elsif ( $response =~ /[mM]/ ) {
      
      # User would like to specify multiple fasta files, one at a time.  Ask for their locations, and scan each for
      # correct formatting.
      
      my $done = 0;

      while ( not $done ) {
	 
	 print "\nPlease enter the name of a file to use, or enter a blank line when you're done: ";

	 $response = <STDIN>;

	 chomp $response;

	 if ( length($response) == 0 ) {
	    
	    $done = 1;

	 } else {
	    
	    my $nextFile = $response;

	    if ( &alreadyEntered($nextFile, \@fileLocs) ) {
	       
	       print "\n   You already entered that file.\n";

	    } elsif ( not -e $nextFile ) {
	       
	       print "\n   The file you specified doesn't appear to exist.\n";

	    } elsif ( not -f $nextFile or -B $nextFile ) {
	       
	       print "\n   The location you entered isn't a regular file.  Please provide a fasta file.\n";

	    } else {
	       
	       print "\n   Checking file format...";

	       my $fileCheckResult = &checkFileFormat($nextFile);

	       print "done.  File complies with standard nucleotide fasta spec.\n";

	       push @fileLocs, $nextFile;
	    }
	 }
      }
      
      my $fileCount = $#fileLocs + 1;

      my $plural = 's';

      if ( $fileCount == 0 ) {
	 
	 die("\n   You didn't enter any files.  Aborting.\n\n");

      } elsif ( $fileCount == 1 ) {
	 
	 $plural = '';
      }

      print "\nDone.  $fileCount file$plural entered and verified as valid.\n\n"; 
   }
}

##########################################################################
# 
# We've now got either a single file to process or a list of multiple files
# to work on.  All files have been verified as valid fasta files.
# 
# If we've got multiple files, find out whether the user wants to create
# a single genome model per file, or to combine all the files into a single
# model.
# 
##########################################################################

my $combineMultipleFiles = 1;

if ( $fileLoc eq '' ) {
   
   # If $fileLoc never got set, we have a list of files instead.  Is there more
   # than one file in the list?  (User could've selected a directory with just
   # one fasta file in it; we need to check.)

   my $fileCount = $#fileLocs + 1;

   if ( $fileCount > 1 ) {
      
      # We've got multiple files in the list.  Poll the user on how to proceed.

      print "OK: you've entered multiple files.  We can proceed in one of two ways.\n\n";

      print "We can construct one genome model per file (recommended if you have, say,\n";

      print "one file per chromosome or plasmid), or we can combine the data in all\n";

      print "the files into a single genome model (recommended if you have, for example,\n";
      
      print "files containing multiple chromosomal contigs, sequenced and assembled\n";
      
      print "from the same organism.)\n\n";

      print "--- NOTE AGAIN that regardless of which option you choose, all models that are\n";

      print "built will be associated with the same species.  If you want to construct models\n";
      
      print "for multiple species, you need to run this script multiple times. ---\n\n";
      
      print "How would you like to proceed?\n\n";

      print "[Enter S for a single combined model (default), or M for multiple models]: ";

      $response = <STDIN>;

      chomp $response;

      while ( length($response) > 1 or ( length($response) == 1 and $response !~ /[sSmM]/ ) ) {
	 
	 print "[Enter S for a single combined model (default), or M for multiple models]: ";

	 $response = <STDIN>;

	 chomp $response;
      }

      if ( $response =~ /m/i ) {
	 
	 # One model per file.

	 $combineMultipleFiles = 0;

      } else {
	 
	 # Combine files into a single model.  This is the default setting for $combineMultipleFiles.
      }

      print "\n";
   }
}

#######################################################################################################
# 
# Prompt the user to regenerate the BLAST database at the end of the automated phase following this
# one.  If they're going to be running this script multiple times in one sitting, they may not want
# to do this bit yet, so we've provided a separate script to do that manually when they're ready.
# 
#######################################################################################################

my $regenBlastPlural = 's';

my $regenBlastVerb = 'have';

if ( $fileLoc != '' or $#fileLocs < 1 ) {
   
   # Just one input file.

   $regenBlastPlural = '';

   $regenBlastVerb = 'has';
}

print "Once your input file${regenBlastPlural} $regenBlastVerb been processed, the local BLAST database will be automatically\n";

print "rebuilt to reflect the new sequence data.  If for any reason you don't want this to happen\n";

print "automatically, you can use the \"rebuildBlastDB.pl\" script to do so manually;\n";

print "this would be useful if you intend to run this script multiple times to add multiple\n";

print "organisms to the database, and don't want to wait for the BLAST DB to be rebuilt each time.\n\n";

print "If you select manual regeneration, please be cautioned that if the DB is not regenerated\n";

print "when you've finished - either automatically on the last run of this script, or manually\n";

print "using \"rebuildBlastDB.pl\" - your PhymmBL results will be *wrong*, because Phymm and\n";

print "BLAST will be operating with different datasets.\n\n";

print "Given all that, would you like to automatically rebuild the BLAST DB when this script has\n";

print "finished (A, default) or rebuild the BLAST DB by hand later on (H)? : ";

$response = <STDIN>;

chomp $response;

while ( length($response) > 1 or ( length($response) == 1 and $response !~ /[aAhH]/ ) ) {
   
   print "\nPlease enter (A) for automatic regeneration when this script has completed,\n";

   print "or (H) for regeneration by hand using \"rebuildBlastDB.pl\" when you've finished\n";

   print "entering all your genome data: ";

   $response = <STDIN>;

   chomp $response;
}

my $rebuildBlastDB = 1;

if ( $response =~ /[hH]/ ) {
   
   # Waive the auto-regeneration.

   $rebuildBlastDB = 0;

} elsif ( $response =~ /[aA]/ or length($response) == 0 ) {
   
   # Default.  $rebuildBlastDB has already been set to 1.
}

#######################################################################################################
# 
# Make a directory for the ICM and training files in a special .userAdded subdirectory of .genomeData,
# and move all the FASTA input into that directory in the appropriate form.
# 
#######################################################################################################

print "\n\n   Thank you.  The user-input phase has now ended.  Feel free to go get some coffee or something. :)\n\n\n";

print "Configuring directory structure...";

if ( not -e '.genomeData/.userAdded' ) {
   
   system("mkdir .genomeData/.userAdded");
}

my $speciesName = $tax->{'genus'} . ' ' . $tax->{'species'};

if ( length($tax->{'strain'}) > 0 ) {
   
   $speciesName .= ' ' . $tax->{'strain'};
}

$speciesName = &dirify($speciesName);

my $newDirPath = ".genomeData/.userAdded/$speciesName";

my $nextFastaID = 0;

if ( not -e $newDirPath ) {
   
   system("mkdir $newDirPath");

} else {
   
   # We're adding to the existing species directory.  Scan any .fna files in there to find out
   # the current unique-prefix index so we don't muck up the existing data.

   opendir DOT, $newDirPath or die("\n\n   Can't open $newDirPath for scanning; it already exists, but there's something wrong...\n\n");

   my @subs = readdir DOT;

   closedir DOT;

   foreach my $sub ( @subs ) {
      
      if ( $sub =~ /\.(\d+)\.fna/ ) {
	 
	 my $index = $1;

	 if ( $index >= $nextFastaID ) {
	    
	    $nextFastaID = $index + 1;
	 }
      }
   }
}

# Save the first FASTA ID we actually wrote to, so that we can properly update the accession map later on.

my $firstFastaID_written = $nextFastaID;

print "done.\n";

my $plural = 's';

if ( $fileLoc ne '' or $#fileLocs == 0 ) {
   
   $plural = '';
}

my $verb = "Copying";

if ( $combineMultipleFiles ) {
   
   $verb = "Concatenating";
}

print "$verb fasta input file${plural}...";

if ( $fileLoc ne '' ) {
   
   # We just have a single genome file.  Rename it to give it a new prefix that we know is going to be unique,
   # and make sure any headers contain the assigned species name so BLAST results will correspond properly to
   # Phymm results.

   my $oldFile = $fileLoc;
   
   my $tempFile = 'temp__header__rewrite.fna';

   my $newFile = "$newDirPath/$speciesName\.$nextFastaID\.fna";

   system("cp $oldFile $tempFile");
   
   open IN, "<$tempFile" or die("Can't open $tempFile for reading.\n");

   open OUT, ">$newFile" or die("Can't open $newFile for writing.\n");

   while ( my $line = <IN> ) {
      
      if ( $line =~ /^>/ ) {
	 
	 $line =~ s/^>/>$speciesName \| /;
      }

      print OUT $line;
   }

   close OUT;

   close IN;

   system("rm $tempFile");

} elsif ( $#fileLocs == 0 ) {
   
   # Again, we just have a single genome file, but this time, it's in the fileLocs array.  Rename it to give
   # it a new prefix that we know is going to be unique.

   my $oldFile = $fileLocs[0];

   my $tempFile = 'temp__header__rewrite.fna';

   my $newFile = "$newDirPath/$speciesName\.$nextFastaID\.fna";

   system("cp $oldFile $tempFile");
   
   open IN, "<$tempFile" or die("Can't open $tempFile for reading.\n");

   open OUT, ">$newFile" or die("Can't open $newFile for writing.\n");

   while ( my $line = <IN> ) {
      
      if ( $line =~ /^>/ ) {
	 
	 $line =~ s/^>/>$speciesName \| /;
      }

      print OUT $line;
   }

   close OUT;

   close IN;

   system("rm $tempFile");

} else {
   
   # We have multiple files.

   if ( $combineMultipleFiles ) {
      
      # We've been asked to combine them into a single model.  Start by concatenating them and
      # dumping the result in the new genome directory.

      my $sysCmd = 'cat ';

      foreach my $file ( @fileLocs ) {
	 
	 $sysCmd .= "$file ";
      }

      my $newFile = "$newDirPath/$speciesName\.$nextFastaID\.fna";

      my $tempFile = 'temp__header__rewrite.fna';

      $sysCmd .= "> $tempFile";

      system($sysCmd);

      open IN, "<$tempFile" or die("Can't open $tempFile for reading.\n");

      open OUT, ">$newFile" or die("Can't open $newFile for writing.\n");

      while ( my $line = <IN> ) {
      
	 if ( $line =~ /^>/ ) {
	 
	    $line =~ s/^>/>$speciesName \| /;
	 }

	 print OUT $line;
      }

      close OUT;

      close IN;

      system("rm $tempFile");

   } else {
      
      # We've been asked to make a single model for each input file.  Dump them all
      # into the new genome directory, renaming as we go.

      foreach my $file ( @fileLocs ) {
	 
	 my $oldFile = $file;

	 my $tempFile = 'temp__header__rewrite.fna';

	 my $newFile = "$newDirPath/$speciesName\.$nextFastaID\.fna";

	 system("cp $oldFile $tempFile");
   
	 open IN, "<$tempFile" or die("Can't open $tempFile for reading.\n");

	 open OUT, ">$newFile" or die("Can't open $newFile for writing.\n");

	 while ( my $line = <IN> ) {
      
	    if ( $line =~ /^>/ ) {
	 
	       $line =~ s/^>/>$speciesName \| /;
	    }

	    print OUT $line;
	 }

	 close OUT;

	 close IN;

	 system("rm $tempFile");

	 $nextFastaID++;
      }
   }
}
   
print "done.\n\n";

# Save the last FASTA ID we actually wrote to, so that we can properly update the accession map later on.

my $lastFastaID_written = $nextFastaID - 1;

#######################################################################################################
# 
# Update the local accession map and taxonomy database.
# 
#######################################################################################################

print "Updating taxonomic data store...";

my $userAccFile = '.taxonomyData/.0_accessionMap/accessionMap_userAdded.txt';

open OUT, ">>$userAccFile" or die("\n\n   Can't open $userAccFile for appending.\n\n");

foreach my $fastaIndex ( $firstFastaID_written .. $lastFastaID_written ) {
   
   my $line = "$speciesName\t$speciesName\.$fastaIndex\n";

   print OUT $line;
}

close OUT;

my $userDistFile = '.taxonomyData/.3_parsedTaxData/distributionOfTaxa_userAdded.txt';

open OUT, ">>$userDistFile" or die("\n\n   Can't open $userDistFile for appending.\n\n");

foreach my $fastaIndex ( $firstFastaID_written .. $lastFastaID_written ) {
   
   my $accession = "$speciesName\.${fastaIndex}";

   foreach my $taxType ( 'phylum', 'class', 'order', 'family', 'genus', 'species' ) {
      
      print OUT "$taxType\t$tax->{$taxType}\t$accession \($tax->{'genus'} $tax->{'species'}\)\t$speciesName\n";
   }
}

close OUT;

print "done.\n\n";

#######################################################################################################
# 
# Construct training files for the FASTA input as needed.
# 
#######################################################################################################

print "Constructing genomic training files for IMMs as needed...";

system(".scripts/createTrainingFiles_userAdded.pl $newDirPath 500");

print "done.\n\n";

#######################################################################################################
# 
# Build ICMs as needed.
# 
#######################################################################################################

print "Building IMMs as needed...";

system(".scripts/buildICMs_userAdded.pl $newDirPath");

print "done.\n\n";

#######################################################################################################
# 
# Rebuild the BLAST DB unless the user asked us not to.
# 
#######################################################################################################

if ( $rebuildBlastDB ) {
   
   print "Rebuilding BLAST database...";

   &rebuildBlastDB();

   print "done.\n\n";
}










##########################################################################
# 
# Subroutines.
# 
##########################################################################

# Check the specified input file to make sure it's a valid [multi]fasta file.  Return 1 for success, 2 for disallowed
# characters in non-comment lines.

sub checkFileFormat {
   
   my $filePath = shift;
   
   my $formatOK = 1;

   open IN, "<$filePath" or die("\nCan't open $filePath for format checking.\n");
   
   my $lineCount = 0;

   while ( my $line = <IN> ) {
      
      $lineCount++;

      chomp $line;

      if ( $line =~ /^>/ or $line =~ /^\;/ ) {
	 
	 # Comment line.  These aren't the droids you're looking for.

      } elsif ( $line =~ /\S\s\S/ ) {
	 
	 # We'll accept leading and trailing spaces just to be helpful, but the sequence data itself can't contain spaces.

	 die("\n\n   There's a space character that isn't supposed to be where it is in file $filePath (line $lineCount).  Aborting.\n\n");

      } elsif ( $line =~ /(\S+)/ ) {
	 
	 # We've got the data from the current data line; check that it's limited to the proper nucleotide codes.

	 my $lineData = $1;

	 if ( $lineData !~ /^[aAcCgGtTuUrRyYkKmMsSwWbBdDhHvVnNxX]+$/ ) {
	    
	    $lineData =~ /([^aAcCgGtTuUrRyYkKmMsSwWbBdDhHvVnNxX])/;

	    my $badChar = $1;

	    die("\n\n   There's a disallowed nucleotide symbol (\"$badChar\") on data line $lineCount in file $filePath\.\n"
	    . "   Allowable symbols are [ACGTURYKMSWBDHVNX].  If you're using a protein fasta file, that won't work properly.\n"
	    . "   Please submit only nucleotide fasta files for processing.  Aborting.\n\n");
	 }
      }
   }
   
   close IN;

   return $formatOK;
}

# Check to see if the file argument is already contained in the array argument.  If it is, return 1, otherwise return 0.

sub alreadyEntered {
   
   my $file = shift;

   my $arrayRef = shift;

   my @arrayArg = @$arrayRef;

   my $maxIndex = $#arrayArg;
   
   if ( $maxIndex == -1 ) {
      
      return 0;

   } else {
      
      for ( my $i = 0; $i <= $maxIndex; $i++ ) {
	 
	 if ( $file eq $arrayArg[$i] ) {
	    
	    return 1;
	 }
      }
   }

   return 0;
}

# Transform an arbitrary ASCII string into a standardized directory name.

sub dirify {
   
   my $result = shift;

   $result =~ s/\_/UNDERSCORE/g;

   $result =~ s/\s+/\_/g;

   $result =~ s/\//SLASH/g;

   $result =~ s/\(/LPAREN/g;

   $result =~ s/\)/RPAREN/g;

   $result =~ s/\'/SINGLEQUOTE/g;

   $result =~ s/\"/DOUBLEQUOTE/g;

   $result =~ s/\:/COLONCOLON/g;

   $result =~ s/\;/SEMICOLON/g;

   return $result;
}

# Rebuild the BLAST DB.

sub rebuildBlastDB {
   
   my $locateBlastBinary = `which blastall 2>&1`;
   
   if ( $locateBlastBinary =~ /no\s+blastall\s+in/ or $locateBlastBinary =~ /command\s+not\s+found/i ) {
      
      die("\n\n   You must have a local copy of the standalone BLAST software installed.  Please see the README for details.\n\n");

   } else {
      
      # Standalone BLAST is installed.
      
      # Concatenate the genome files and build a local BLAST DB from them.
      # 
      # First, grab a list of relative paths to all the RefSeq genomic FASTA files.
      
      my @fastaFiles = ();
      
      opendir DOT, '.genomeData' or die("\n\n   Can't open .genomeData for scanning.\n\n");

      my @subs = readdir DOT;

      closedir DOT;

      foreach my $sub ( @subs ) {
	 
	 if ( -d ".genomeData/$sub" and $sub !~ /^\./ ) {
	    
	    opendir DOT, ".genomeData/$sub" or die("\n\n   Can't open .genomeData/$sub for scanning.\n\n");

	    my @files = readdir DOT;

	    closedir DOT;

	    foreach my $file ( @files ) {
	       
	       if ( $file =~ /\.fna$/ ) {
		  
		  push @fastaFiles, ".genomeData/$sub/$file";
	       }
	    }
	 }
      }
      
      # Next, go into .genomeData/.userAdded and grab those too.
      
      my $userFolder = '.genomeData/.userAdded';
      
      if ( -e $userFolder ) {
	 
	 opendir DOT, $userFolder or die("\n\n   Can't open $userFolder for scanning.\n\n");

	 @subs = readdir DOT;

	 closedir DOT;

	 foreach my $sub ( @subs ) {

	    if ( -d "$userFolder/$sub" and $sub ne '.' and $sub ne '..' ) {
	    
	       opendir DOT, "$userFolder/$sub" or die("\n\n   Can't open $userFolder/$sub for scanning.\n\n");

	       my @files = readdir DOT;

	       closedir DOT;

	       foreach my $file ( @files ) {
	       
		  if ( $file =~ /\.fna$/ ) {
		  
		     push @fastaFiles, "$userFolder/$sub/$file";
		  }
	       }
	    }
	 }

	 # Concatenate all the FASTA files together for database creation.
      
	 my $command = 'cat ';
	 
	 my $argCounter = 0;
	 
	 my $wrapped = 0;

	 foreach my $file ( @fastaFiles ) {
	 
	    $command .= "$file \\\n";

	    $argCounter++;

	    if ( $argCounter == 10 ) {
	       
	       if ( not $wrapped ) {
		  
		  $wrapped = 1;

		  $command .= "> blastDB_data.txt\n";

		  $command .= 'cat ';

	       } else {
		  
		  $command .= ">> blastDB_data.txt\n";

		  $command .= 'cat ';
	       }

	       $argCounter = 1;
	    }
	 }
	 
	 if ( not $wrapped ) {
	    
	    $command .= "> blastDB_data.txt\n";

	 } else {
	    
	    $command .= ">> blastDB_data.txt\n";
	 }
	 
	 # There's a problem with using system() for a command this long, apparently, so we use a hack.

	 my $shLoc = `which sh`;

	 chomp $shLoc;

	 my $firstLine = "#!$shLoc\n";

	 my $secondLine = "\n";

	 my $thirdLine = $command;

	 my $outFile = "command_temp.sh";

	 open OUT, ">$outFile" or die("Can't open $outFile for writing.\n");

	 print OUT $firstLine;

	 print OUT $secondLine;

	 print OUT $thirdLine;

	 close OUT;

	 system("chmod 775 command_temp.sh");
	 
	 system("./command_temp.sh");
	 
	 system("rm command_temp.sh");

	 # Create the local database copy.

	 system("formatdb -i blastDB_data.txt -p F -n phymm_BLAST_DB");

	 system("mv phymm_BLAST_DB.* .blastData");

	 system("rm blastDB_data.txt");

	 system("rm formatdb.log");
      
      } # end if ( there's a '.genomeData/.userAdded' directory )

   } # end if ( local copy of BLAST binary is installed )

}



