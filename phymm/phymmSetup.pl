#!/usr/bin/perl

use strict;

$| = 1;

use Net::FTP;

# Figure out whether to obtain a new copy of the RefSeq genomic database,
# use an existing one, or refresh existing data with new data from RefSeq.

my $newCopy = 0;

&checkDataSource();

if ( $newCopy == 0 or $newCopy == 1 ) {
   
   # USE EXISTING genome database ($newCopy == 0).  This is essentially
   # a cleanup operation.
   # 
   # --> OR <--
   # 
   # OVERWRITE genome database with new data ($newCopy == 1): We've been
   # told to completely overwrite existing genome data with current RefSeq
   # data, treating this as a first-time installation.
   
   # 0. At this point, the data's already been downloaded and moved to the proper place.

   # 1. Rename directories as necessary.

   &renameNewRefSeqDirs('.genomeData');
   
   # 2. Generate metafiles describing the RefSeq data, using NCBI's taxonomy data to help.

   &generateRefSeqMetafiles();

   # 3. Extract FASTA sequences from the GenBank files as necessary.

   &genbankToFasta();
   
   # 4. Create training files for the ICMs as necessary.

   &createTrainingFiles();
   
   # 5. Build/update the ICM suite.

   &buildICMs();
   
   # 6. Build the local BLAST database from the RefSeq files.

   &generateBlastDB();

} elsif ( $newCopy == 2 ) {
   
   # REFRESH existing genome database with changes from RefSeq since last update.
   
   # We're only doing two things here: adding newly-deposited genomes,
   # and updating existing ones.  We're /not/ deleting genomes which aren't
   # in the current RefSeq instance (for several reasons, including preserving
   # manually-added organisms).  If you want old organisms removed, just delete
   # the subdirectory of .genomeData corresponding to the offending organism(s).
   
   # 0. At this point, the new RefSeq data has been extracted
   #    to '.tempGenomeData', and the existing genome data remains
   #    untouched in '.genomeData'.
   
   # 1 Rename genome directories - in the existing genome database - as necessary.

   &renameNewRefSeqDirs('.genomeData');

   # 2 Rename genome directories - in the temporary holding bin containing the new RefSeq data - as necessary.

   &renameNewRefSeqDirs('.tempGenomeData');
   
   # 3. Perform genome data updates as needed.
   
   &updateGenomeData();
   
   # 4. Remove the temp dir.  We're done with it.
   
   print "Removing temporary genome data directory...";

   system("rm -rf .tempGenomeData");

   print "done.\n\n";

   # 5. Extract FASTA sequences from the GenBank files as necessary.

   &genbankToFasta();
   
   # 6. Generate metafiles describing the RefSeq data, using NCBI's taxonomy data to help.

   &generateRefSeqMetafiles();

   # 7. Create training files for the ICMs as necessary.

   &createTrainingFiles();
   
   # 8. Build/update the ICM suite.

   &buildICMs();
   
   # 9. Build the local BLAST database from the RefSeq files.

   &generateBlastDB();
}








###########################################################
# SUBROUTINES
###########################################################

sub checkDataSource {
   
   if ( &checkRefSeqData() == 1 ) {
      
      print "\nLocal genome database detected.  Please choose one of the following:\n\n";

      print "Use (E)xisting genome data without changing it\n\n";

      print "(R)efresh database with newest RefSeq data (will not alter genomes which haven't been\n";
      print "changed in the RefSeq database since Phymm was installed)\n\n";

      print "completely (O)verwrite existing genome data with fresh/current copy from RefSeq\n\n";

      print "[Default is (E)]: ";

      my $response = <STDIN>;

      chomp $response;

      while ( length($response) > 1 or ( length($response) == 1 and $response !~ /[eErRoO]/ ) ) {
	 
	 print "\nPlease enter E (default: use existing genome data), R (refresh with newest data) or O (overwrite with new copy): ";

	 $response = <STDIN>;

	 chomp $response;
      }
      
      print "\n";
      
      if ( $response =~ /[oO]/ ) {
	 
	 # Grab the newest data from the server.
	 
	 &downloadNewRefSeqData();
	 
	 # Wipe .genomeData.  If the .userAdded subdirectory's in there, this syntax should leave it alone.

	 system("rm -rf .genomeData/*");

	 # Extract downloaded tarball into .genomeData.
	 
	 &expandNewRefSeqData_genomeData();

	 $newCopy = 1;
	 
      } elsif ( $response =~ /[rR]/ ) {

	 # Refresh data with genomes added since install, and update genomes which've been altered since install.

	 &downloadNewRefSeqData();
	 
	 # Extract downloaded tarball into a temp directory for comparison/update.

	 &expandNewRefSeqData_tempDir();

	 $newCopy = 2;

      } else {
	 
	 # We're going to use existing data.

	 $newCopy = 0; # (reaffirming default value)
      }

   } else {
      
      # In this case, $checkRefSeqData() has returned 0 and there's no tarball and no data in .genomeData.  Do a fresh install.

      print "\nNo RefSeq genomic data detected; download now? [Y/n] ";

      my $response = <STDIN>;

      chomp $response;

      while ( length($response) > 1 or ( length($response) == 1 and $response !~ /[yYnN]/ ) ) {
	 
	 print "\nPlease enter Y (default: get a new copy of RefSeq genome data) or N (don't): ";

	 $response = <STDIN>;

	 chomp $response;
      }
      
      print "\n";
      
      if ( $response !~ /[nN]/ ) {
	 
	 # Grab the newest data from the server.
	 
	 &downloadNewRefSeqData();
	 
	 # Extract downloaded tarball into .genomeData.
	 
	 &expandNewRefSeqData_genomeData();

	 $newCopy = 1;

      } else {
	 
	 die("Can't continue Phymm setup without genomic data; exiting.\n\n");
      }
   }
}

# Check for the existence of a gzipped tarball from the RefSeq FTP server,
# as well as any content in the .genomeData directory.

sub checkRefSeqData {
   
   if ( not -e '.genomeData' ) {
      
      # For whatever reason, the .genomeData directory doesn't (or has ceased to) exist.

      die("\nDirectory \".genomeData\" is missing; can't continue with setup.\n\n");

   } else {
      
      # The .genomeData directory exists, as expected.

      # If there's *anything at all* in the '.genomeData' directory, assume it's a genomic library and return 1.
   
      opendir DATA, ".genomeData";

      my @subs = readdir DATA;

      closedir DATA;

      foreach my $sub ( @subs ) {

	 if ( -d ".genomeData/$sub" and $sub !~ /^\./ ) {
	    
	    return 1;
	 }
      }
      
      # The .genomeData directory is empty.  Return 0.

      return 0;
   }
}

# Download a current tarball from RefSeq.

sub downloadNewRefSeqData {
   
   my $ftpSite = 'ftp.ncbi.nih.gov';

   print "Connecting to $ftpSite...";

   my $ftp = Net::FTP->new($ftpSite, Passive => 1, Debug => 0, Timeout => 2400) or die "Can't connect to $ftpSite: $@";

   $ftp->login() or die "Can't login: " . $ftp->message;

   print "done.\n\n";

   my $targetDir = '/genomes/Bacteria';

   $ftp->cwd($targetDir) or die "Can't change remote directory: " . $ftp->message;

   my $targetFile = 'all.gbk.tar.gz';

   print "Downloading $targetFile (this will take a while; the file is over 2GB)...";

   $ftp->binary();
   
   $ftp->get($targetFile, "./$targetFile") or die "Failed to download file: " . $ftp->message;

   print "done.\n\n";

   $ftp->quit;
}

# Unzip/untar the downloaded genome-data file into the .genomeData directory.

sub expandNewRefSeqData_genomeData {
   
   print "Expanding downloaded genome data for database overwrite (this'll take a while)...";

   system("tar -zxvf all.gbk.tar.gz -C .genomeData >& /dev/null");
   
   print "done.\n\nRemoving compressed genome data...";
   
   system("rm -rf all.gbk.tar.gz");
   
   print "done.\n\n";
}

# Unzip/untar the downloaded genome-data file into a temp directory for comparison/update.

sub expandNewRefSeqData_tempDir {
   
   print "Expanding downloaded genome data for database update (this'll take a while)...";
   
   system("mkdir .tempGenomeData");

   system("tar -zxvf all.gbk.tar.gz -C .tempGenomeData >& /dev/null");
   
   print "done.\n\nRemoving tarball...";
   
   system("rm -rf all.gbk.tar.gz");
   
   print "done.\n\n";
}

# Standardize the names of the RefSeq organism directories to match the names given in the .gbk files.

sub renameNewRefSeqDirs {
   
   # Get the name of the genome data directory to scan.

   my $dataDir = shift;

   # Load known GenBank SOURCE-field inconsistencies for later processing.

   my $GBerrors = {};

   &loadGBerrors($GBerrors);

   print "Standardizing organism directory names and file data in $dataDir...";
   
   opendir DATA, "$dataDir" or die("Can't open $dataDir for scanning.\n");

   my @subs = readdir DATA;

   closedir DATA;
   
   # Iterate over subdirectories of .genomeData.

   foreach my $sub ( @subs ) {
      
      # Organism-wide variable.  Used for sanity checking.

      my $orgName = '';
      
      if ( -d "$dataDir/$sub" and $sub !~ /^\./ ) {
	 
	 opendir SUB, "$dataDir/$sub" or die("Can't open $dataDir/$sub for scanning.\n");

	 my @files = readdir SUB;

	 closedir SUB;
	 
	 # Delete unwanted files.

	 foreach my $file ( @files ) {
	    
	    if ( $file =~ /\%$/ ) {
	       
	       # This file's been obsolesced by NCBI, but for whatever reason remained in the distribution.  Kill it.
	       
	       system("rm $dataDir/$sub/$file");
	    }

	    # We only want accessions beginning with NC_ - all others are in varying stages of unreadiness.
	    
	    # Also, NC_008265 is a phage and not a bacterium; delete it.

	    if ( $file =~ /\.gbk$/ and ( $file !~ /NC_/ or $file =~ /NC_008265/ ) ) {
	       
	       system("rm $dataDir/$sub/$file");
	       
	    }

	 } # end ( loop deleting unwanted files )
	 
	 # Refresh the file list.  We may have deleted some.

	 opendir SUB, "$dataDir/$sub" or die("Can't open $dataDir/$sub for scanning.\n");

	 @files = readdir SUB;

	 closedir SUB;
	 
	 # Go through each GenBank file in the current subdirectory of .genomeData.
	 # 
	 # 1. Grab the (possibly multi-line) organism name from the SOURCE field.
	 # 
	 # 2. If $orgName hasn't been set yet, set it to the current (retrieved) organism name.
	 # 
	 # 3. If $orgName has been set, then compare it to the current (retrieved) organism name.
	 #    If they don't match, it means there are different listings in the SOURCE field
	 #    across different files in the same subdirectory.  Check against a list of known
	 #    GenBank inconsistencies and fix if possible.  Otherwise, complain with a warning.
	 
	 foreach my $file ( @files ) {
	    
	    if ( $file =~ /\.gbk$/ ) {
	       
	       my $fullFile = "$dataDir/$sub/$file";
	       
	       open IN, "<$fullFile" or die("Can't open $fullFile for reading.\n");
	       
	       my $found = 0;
	       
	       my $recording = 0;
	       
	       my $currentOrgName = '';
	       
	       while ( not $found and my $line = <IN> ) {
		  
		  chomp $line;

		  if ( $recording and $line !~ /^\s*ORGANISM/ ) {
		     
		     $line =~ /^\s+(.*)/;

		     $currentOrgName .= ' ' . $1;

		  } elsif ( $line =~ /^SOURCE\s+(.*)/ ) {
		     
		     $currentOrgName = $1;

		     $recording = 1;
		     
		  } elsif ( $line =~ /^\s*ORGANISM/ ) {
		     
		     $found = 1;
		  }
	       }

	       close IN;
	       
	       # Sanity check on organism name.
	       
	       $currentOrgName = &fixOrgName($currentOrgName, $GBerrors);

	       if ( $orgName eq '' ) {
		  
		  $orgName = $currentOrgName;
	       
	       } elsif ( $orgName ne $currentOrgName ) {
		  
		  print "\n\n   WARNING: Different source organisms given for different .gbk files in data directory \"$dataDir/$sub\": \"$orgName\" and \"$currentOrgName\".\n\n";
	       }
	       
	    } # end ( filename check for /\.gbk$/ )
	    
	 } # end ( file iterator on current subdirectory )
	 
	 # Check to ensure we haven't just emptied the current subdirectory.  If we have, delete it.

	 opendir SUB, "$dataDir/$sub" or die("Can't open $dataDir/$sub for scanning.\n");

	 @files = readdir SUB;

	 closedir SUB;

	 my $found = 0;

	 foreach my $file ( @files ) {
	    
	    if ( $file ne '.' and $file ne '..' ) {
	       
	       $found = 1;
	    }
	 }

	 if ( not $found ) {
	    
	    # No files left in this subdirectory.  Delete it.
	    
	    system("rm -rf $dataDir/$sub\n");
	    
	 }
	 
	 if ( $found ) {
	    
	    # The subdirectory isn't empty.  Continue.

	    # If we're not dealing with a completely alien directory - i.e., if we
	    # successfully found an $orgName - process the directory, otherwise ignore.

	    if ( $orgName ne '' ) {
	       
	       # Rename the current subdirectory to a standard name based on the species name quoted in its .gbk files.

	       my $newSubName = &toDirName($orgName);
	       
	       if ( $newSubName ne $sub ) {
		  
		  # We've changed the directory name.

		  if ( not -e "$dataDir/$newSubName" ) {
	       
		     # There's not yet a directory labeled with the new name.  Just move the old one over.
	       
		     my $tempSub = $sub;

		     $tempSub =~ s/\(/\\\(/g;

		     $tempSub =~ s/\)/\\\)/g;

		     $tempSub =~ s/\'/\\\'/g;

		     $tempSub =~ s/\"/\\\"/g;

		     my $command = "mv $dataDir/$tempSub $dataDir/$newSubName\n";

		     system($command);

		  } else {
	       
		     # We've made a change, but the target directory (labeled with the
		     # new name) already exists; move the /files/ from the old dir to
		     # the new one.
	       
		     my $tempSub = $sub;

		     $tempSub =~ s/\(/\\\(/g;

		     $tempSub =~ s/\)/\\\)/g;

		     $tempSub =~ s/\'/\\\'/g;

		     $tempSub =~ s/\"/\\\"/g;

		     my $command = "mv $dataDir/$tempSub/* $dataDir/$newSubName\n";

		     system($command);
	       
		     # Delete the now-empty directory labeled with the old name.
	       
		     $command = "rm -rf $dataDir/$tempSub\n";

		     system($command);
		  }

	       } # end ( check to see if we changed the directory name from what it was )

	    } # end ( check to discover whether or not we've been able to establish an $orgName )
	    
	 } # end ( check to ensure current subdirectory isn't empty )
	 
      } # end ( filter on /^\./ )
      
   } # end foreach ( subdirectory in .genomeData )
   
   print "done.\n\n";
}

# Generate all needed metadata for later Phymm processing.

sub generateRefSeqMetafiles {
   
   print "Generating taxonomic metafiles...";
   
   my $species = {};
   
   my $definition = {};
   
   my $dataDir = '.genomeData';

   opendir DATA, "$dataDir" or die("Can't open $dataDir for scanning.\n");

   my @subs = readdir DATA;

   closedir DATA;
   
   # Go through each subdirectory of .genomeData and grab the defline.  Reconstruct the organism name from the directory name.

   foreach my $sub ( @subs ) {
      
      my $defLine = '';
      
      my $orgName = '';
      
      if ( -d "$dataDir/$sub" and $sub !~ /^\./ ) {
	 
	 opendir SUB, "$dataDir/$sub" or die("Can't open $dataDir/$sub for scanning.\n");

	 my @files = readdir SUB;

	 closedir SUB;
	 
	 $orgName = &fromDirName($sub);

	 foreach my $file ( @files ) {
	    
	    if ( $file =~ /(.+)\.gbk$/ ) {
	       
	       my $prefix = $1;
	       
	       my $fullFile = "$dataDir/$sub/$file";

	       open IN, "<$fullFile" or die("Can't open $fullFile for reading.\n");
	       
	       my $found = 0;
	       
	       $defLine = '';
	       
	       my $recording = 0;
	       
	       while ( not $found and my $line = <IN> ) {
		  
		  chomp $line;

		  if ( $recording and $line !~ /^ACCESSION/ ) {
		     
		     $line =~ /^\s+(.*)/;

		     $defLine .= ' ' . $1;

		  } elsif ( $line =~ /^DEFINITION\s+(.*)/ ) {
		     
		     $defLine = $1;

		     $recording = 1;
		     
		  } elsif ( $line =~ /^ACCESSION/ ) {
		     
		     $found = 1;
		  }
	       }
	       
	       $species->{$prefix} = $orgName;

	       $definition->{$prefix} = $defLine;
	       
	    } # end ( filename check for /\.gbk$/ )
	    
	 } # end ( file iterator on current subdirectory )
	 
      } # end ( filter on /^\./ )
      
   } # end foreach ( subdirectory in .genomeData )
   
   # Output all the saved data to an accession map file.

   my $idMapFile = '.taxonomyData/.0_accessionMap/accessionMap.txt';

   open OUT, ">$idMapFile" or die("Can't open $idMapFile for writing.\n");
   
   foreach my $prefix ( sort { $species->{$a} cmp $species->{$b} || $definition->{$a} cmp $definition->{$b} } keys %$species ) {
      
      my $orgName = $species->{$prefix};

      my $defLine = $definition->{$prefix};

      my $seqType = '';

      if ( $defLine =~ /chromosom/ or $defLine =~ /complete\s+genome/ ) {
	 
	 $seqType = 'chromosome';

      } elsif ( $defLine =~ /plasmid/ ) {
	 
	 $seqType = 'plasmid';

      } elsif ( $defLine =~ /complete\s+sequence/ ) {
	 
	 $seqType = 'chromosome';

      } else {
	 
	 die("Bad 'DEFINITION' line: $prefix \/\/ $orgName \/\/ $defLine\n");
      }

      my $newSubName = &toDirName($orgName);

      print OUT "$newSubName\t$prefix\t$seqType\t$defLine\n";
   }
   
   close OUT;
   
   # Grab the NCBI taxonomy, parse, and save in a format we can cope with.

   system(".scripts/.taxonomyInfo/0_getNCBIpages.pl");
   
   system(".scripts/.taxonomyInfo/1_stripHTML_bacteria.pl");
   
   system(".scripts/.taxonomyInfo/2_getTree_bacteria.pl");
   
   system(".scripts/.taxonomyInfo/3_stripHTML_archaea.pl");

   system(".scripts/.taxonomyInfo/4_getTree_archaea.pl");

   system(".scripts/.taxonomyInfo/5_configureTaxonomicData.pl");
   
   print "done.\n\n";
}

# Convert the .gbk files that came with the RefSeq distribution to corresponding .fna files which contain only the nucleotide sequence from the .gbk file.

sub genbankToFasta {
   
   print "Extracting nucleotide sequences from .gbk files into .fna files as necessary...";
      
   system(".scripts/genbankToFasta.pl .genomeData");

   print "done.\n\n";
}

# Create training files for the ICMs if they don't already exist.

sub createTrainingFiles {
   
   print "Creating training files for IMMs as necessary...";
      
   system(".scripts/createTrainingFiles.pl .genomeData 500");

   print "done.\n\n";
}

# Build the ICMs which form the core of Phymm's classification system

sub buildICMs {
   
   # First, check to see if the ICM code from the Glimmer package has been untarred and/or compiled yet.  Do whatever needs doing.

   if ( not -e ".scripts/.icmCode/LICENSE" ) {
      
      # The package hasn't been uncompressed yet.  Do so.
      
      print "Uncompressing IMM code...";
      
      system("tar -zxvf .scripts/.icmCode/glimmer3_plusSS.tgz -C .scripts/.icmCode >& /dev/null");
      
      print "done.\n\n";
   }
   
   # Next, resolve a conflict based on mutually exclusive syntax requirements across the gcc 4.3/4.4 barrier.

   open (CMD, "(gcc -v) 2>&1 |");

   while ( my $line = <CMD> ) {
   
      if ( $line =~ /gcc\s+version\s+(\S+)/ ) {
      
	 my $versionString = $1;
      
	 my @versionNums = ();

	 while ( $versionString =~ /^(\d+)\./ ) {
	 
	    my $nextNum = $1;

	    push @versionNums, $nextNum;

	    $versionString =~ s/^$nextNum\.//;
	 }
      
	 if ( $versionString =~ /^(\d+)/ ) {
	 
	    push @versionNums, $1;
	 }

	 my $topNum = $versionNums[0];
      
	 my $secondNum;
      
	 if ( $#versionNums > 0 ) {
	 
	    $secondNum = $versionNums[1];
	 }

	 my $newGCC = 0;

	 if ( $topNum > 4 or ( $topNum == 4 and $secondNum >= 4 ) ) {
	 
	    $newGCC = 1;
	 }
	 
	 # Rename the appropriate files, if we haven't already done this.

	 if ( $newGCC == 0 ) {
	    
	    # We've got an older version of gcc.

	    if ( not -e '.scripts/.icmCode/SimpleMake/icm.cc' ) {
	       
	       system("mv .scripts/.icmCode/SimpleMake/icm_oldgcc.cc .scripts/.icmCode/SimpleMake/icm.cc");
	       system("mv .scripts/.icmCode/SimpleMake/gene_oldgcc.cc .scripts/.icmCode/SimpleMake/gene.cc");
	       system("mv .scripts/.icmCode/SimpleMake/delcher_oldgcc.hh .scripts/.icmCode/SimpleMake/delcher.hh");
	    }

	 } else {
	    
	    # We've got a newer version of gcc.

	    if ( not -e '.scripts/.icmCode/SimpleMake/icm.cc' ) {
	       
	       system("mv .scripts/.icmCode/SimpleMake/icm_newgcc.cc .scripts/.icmCode/SimpleMake/icm.cc");
	       system("mv .scripts/.icmCode/SimpleMake/gene_newgcc.cc .scripts/.icmCode/SimpleMake/gene.cc");
	       system("mv .scripts/.icmCode/SimpleMake/delcher_newgcc.hh .scripts/.icmCode/SimpleMake/delcher.hh");
	    }
	 }
      }
   }

   close CMD;

   # Now we can compile.

   if ( not -e ".scripts/.icmCode/bin/build-icm" ) {
      
      # The code hasn't been compiled yet.  Do so.
      
      print "Compiling IMM code...";
      
      chdir(".scripts/.icmCode/SimpleMake");

      system("make >& /dev/null");

      chdir("../../..");
      
      print "done.\n\n";
   }
   
   # If there was a problem with the compilation, we just
   # redirected any messages to that effect to a black hole.
   # Check for one of the more critical binaries.

   if ( not -e ".scripts/.icmCode/bin/build-icm" ) {
      
      die("The IMM code does not appear to have compiled successfully.\n\nPlease try to compile manually (cd .scripts/.icmCode/SimpleMake; make);\nif that doesn't work, check to make sure you have\nthe appropriate version of gcc (gcc -v) for the version of Phymm\nthat you've downloaded (see website for details).\n");
   }

   # If ICMs haven't yet been built, build them.
   
   system(".scripts/buildICMs.pl .genomeData");
}

# Count the number of ICM training files in the .genomeData hierarchy.

sub countTrainingFiles {
   
   my $retVal = 0;
   
   my $root = '.genomeData';

   opendir ROOT, $root or die("Can't open $root for scanning.\n");

   my @subs = readdir ROOT;

   closedir ROOT;

   foreach my $sub ( @subs ) {
      
      if ( -d "$root/$sub" and $sub !~ /^\./ ) {
	 
	 opendir SUB, "$root/$sub" or die("Can't open $root/$sub for scanning.\n");

	 my @files = readdir SUB;

	 closedir SUB;

	 foreach my $file ( @files ) {
	    
	    if ( $file =~ /GenomicTrainingSeqs/ ) {
	       
	       $retVal++;
	    }
	 }
      }
   }
   
   return $retVal;
}

# Build a local BLAST database from the RefSeq FASTA files.

sub generateBlastDB {
   
   my $locateBlastBinary = `which blastall 2>&1`;
   
   if ( $locateBlastBinary =~ /no\s+blastall\s+in/ or $locateBlastBinary =~ /command\s+not\s+found/i ) {
      
      die("You must have a local copy of the standalone BLAST software installed.  Please see the README for details.\n\n");

   } else {
      
      # Standalone BLAST is installed.  Do we need to regenerate the BLAST DB?  If there's no DB, generate one.
      
      my $regenDB = 0;

      if ( not -e ".blastData/phymm_BLAST_DB.nsq" ) {
	 
	 $regenDB = 1;

      } else {
	 
	 # If there is one, search the lastmod timestamps of the .fna files in .genomeData.
	 # If we find one that's newer than the blast DB, regenerate the DB.
	 # 
	 # We won't check .fna files in .genomeData/.userAdded here, because if the user has
	 # used addCustomGenome.pl, then they've been repeatedly instructed to remember to
	 # rebuild the BLAST DB when they were done (if they selected the option to do it by
	 # hand in the first place, that is; by default, it's rebuilt automatically).  If they
	 # decided not to do it despite repeated cautions, well then who knows, they might
	 # have their reasons.

	 my $blastTimeStamp = (stat('.blastData/phymm_BLAST_DB.nsq'))[9];
	 
	 my $genomeDir = '.genomeData';

	 opendir DOT, $genomeDir or die("Can't open $genomeDir for scanning.\n");

	 my @subs = readdir DOT;

	 closedir DOT;

	 OUTER: foreach my $sub ( @subs ) {
	    
	    my $fullSub = "$genomeDir/$sub";

	    if ( -d $fullSub and $sub !~ /^\./ ) {
	       
	       opendir DOT, $fullSub or die("Can't open $fullSub for scanning.\n");

	       my @files = readdir DOT;

	       closedir DOT;

	       foreach my $file ( @files ) {
		  
		  my $fullFile = "$fullSub/$file";

		  if ( $file =~ /\.fna$/ ) {
		     
		     my $currentTimeStamp = (stat($fullFile))[9];

		     if ( $currentTimeStamp > $blastTimeStamp ) {
			
			# Found one.  Stop looking.

			$regenDB = 1;

			last OUTER;
		     }
		  }
	       }
	    }
	 }
      }

      if ( $regenDB ) {
	 
	 # Concatenate the genome files and build a local BLAST DB from them.
	 
	 print "Generating local BLAST database from RefSeq data...";
      
	 # First, grab a list of relative paths to all the RefSeq genomic FASTA files.
      
	 my @fastaFiles = ();
      
	 opendir DOT, '.genomeData' or die("Can't open .genomeData for scanning.\n");

	 my @subs = readdir DOT;

	 closedir DOT;

	 foreach my $sub ( @subs ) {
	 
	    if ( -d ".genomeData/$sub" and $sub !~ /^\./ ) {
	    
	       opendir DOT, ".genomeData/$sub" or die("Can't open .genomeData/$sub for scanning.\n");

	       my @files = readdir DOT;

	       closedir DOT;

	       foreach my $file ( @files ) {
	       
		  if ( $file =~ /\.fna$/ ) {
		  
		     push @fastaFiles, ".genomeData/$sub/$file";
		  }
	       }
	    }
	 }
	 
	 # Next, grab a list of relative paths to any user-added FASTA files.
	 
	 my $userDir = '.genomeData/.userAdded';
	 
	 if ( -e $userDir ) {
	    
	    opendir DOT, $userDir or die("Can't open $userDir for scanning.\n");

	    @subs = readdir DOT;

	    closedir DOT;

	    foreach my $sub ( @subs ) {
	       
	       if ( -d "$userDir/$sub" and $sub !~ /^\./ ) {
		  
		  opendir DOT, "$userDir/$sub" or die("Can't open $userDir/$sub for scanning.\n");

		  my @files = readdir DOT;

		  closedir DOT;

		  foreach my $file ( @files ) {
		     
		     if ( $file =~ /\.fna$/ ) {
			
			push @fastaFiles, "$userDir/$sub/$file";
		     }
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
      
	 print "done.\n\n";

      } # end if ( we're generating a DB at all )

   } # end if ( local copy of BLAST binary is installed )
}

sub toDirName {
   
   my $name = shift;

   my $result = $name;

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

sub fromDirName {
   
   my $name = shift;

   my $result = $name;

   $result =~ s/\_/ /g;

   $result =~ s/UNDERSCORE/\_/g;

   $result =~ s/SLASH/\//g;

   $result =~ s/LPAREN/\(/g;

   $result =~ s/RPAREN/\)/g;

   $result =~ s/SINGLEQUOTE/\'/g;

   $result =~ s/DOUBLEQUOTE/\"/g;
	    
   $result =~ s/COLONCOLON/\:/g;

   $result =~ s/SEMICOLON/\;/g;

   return $result;
}

sub loadGBerrors {
   
   my $errors = shift;

   my $errorFile = '.taxonomyData/knownGBerrors.txt';

   open IN, "<$errorFile" or die("Can't open $errorFile for reading.\n");

   while ( my $line = <IN> ) {
      
      chomp $line;

      (my $canonName, my $pattern) = split(/\t/, $line);

      $errors->{$pattern} = $canonName;
   }

   close IN;
}

sub fixOrgName {
   
   my $orgName = shift;

   my $errors = shift;

   foreach my $pattern ( keys %$errors ) {
      
      my $correctName = $errors->{$pattern};

      if ( $orgName =~ /$pattern/ ) {
	 
	 $orgName = $correctName;

	 return $orgName;
      }
   }

   return $orgName;
}

sub updateGenomeData {

   # Classify each genome directory in the new RefSeq set as one of the following:
   #    [new organism added since existing installation was last updated],
   #    [organism contained in existing installation, but /hasn't/ been changed], or
   #    [organism contained in existing installation and /has/ been changed].

   
   # Load the existing genome directory list, along with .fna file sizes.

   my $extantGenomes = {};

   opendir GENOMEDATA, '.genomeData' or die("Can't open .genomeData for scanning.\n");

   my @subs = readdir GENOMEDATA;

   closedir GENOMEDATA;

   foreach my $sub ( @subs ) {
      
      if ( -d ".genomeData/$sub" and $sub !~ /^\./ ) {
	 
	 opendir DOT,  ".genomeData/$sub" or die("Can't open .genomeData/$sub for scanning.\n");

	 my @files = readdir DOT;

	 closedir DOT;

	 foreach my $file ( @files ) {
	    
	    if ( $file =~ /\.fna$/ ) {
	       
	       my $fileSize = -s ".genomeData/$sub/$file";

	       $extantGenomes->{$sub}->{$file} = $fileSize;
	    }
	 }
      }
   }
   
   # Create .fna files for the new download.

   print "Creating temporary FASTA files from new .gbk files for update check...";
      
   system(".scripts/genbankToFasta.pl .tempGenomeData");

   print "done.\n\n";

   # Load the new genome directory list, along with .fna file sizes.

   my $newGenomes = {};

   opendir GENOMEDATA, '.tempGenomeData' or die("Can't open .tempGenomeData for scanning.\n");

   my @subs = readdir GENOMEDATA;

   closedir GENOMEDATA;

   foreach my $sub ( @subs ) {
      
      if ( -d ".tempGenomeData/$sub" and $sub !~ /^\./ ) {
	 
	 opendir DOT,  ".tempGenomeData/$sub" or die("Can't open .tempGenomeData/$sub for scanning.\n");

	 my @files = readdir DOT;

	 closedir DOT;

	 foreach my $file ( @files ) {
	    
	    if ( $file =~ /\.fna$/ ) {
	       
	       my $fileSize = -s ".tempGenomeData/$sub/$file";

	       $newGenomes->{$sub}->{$file} = $fileSize;
	    }
	 }
      }
   }

   # First, find the completely new critters, if any.
   
   my $addedGenomeHash = {};

   foreach my $newGenome ( keys %$newGenomes ) {
      
      if ( not exists($extantGenomes->{$newGenome}) ) {
	 
	 $addedGenomeHash->{$newGenome} = 1;
      }
   }
   
   my @addedGenomes = sort { $a cmp $b } keys %$addedGenomeHash;

   if ( $#addedGenomes == -1 ) {
      
      print "No new genomes to add.\n\n";

   } else {
      
      print "Adding new genomes:\n\n";

      foreach my $newGenome ( @addedGenomes ) {
	 
	 print "   $newGenome\n";
      }

      print "\n";
   }

   # Next, identify existing genomes which need to be updated, if any.
   
   my $updatedGenomeHash = {};

   foreach my $newGenome ( keys %$newGenomes ) {
      
      foreach my $file ( keys %{$newGenomes->{$newGenome}} ) {
	 
	 if ( not exists($extantGenomes->{$newGenome})
	   or not exists($extantGenomes->{$newGenome}->{$file})
	   or $extantGenomes->{$newGenome}->{$file} != $newGenomes->{$newGenome}->{$file} ) {
	    
	    if ( not exists($addedGenomeHash->{$newGenome}) ) {
	       
	       $updatedGenomeHash->{$newGenome}->{$file} = 1;
	    }
	 }
      }
   }

   my @updatedGenomes = sort { $a cmp $b } keys %$updatedGenomeHash;

   if ( $#updatedGenomes == -1 ) {
      
      print "No updates found for existing genomes.\n\n";

   } else {
      
      print "Updating existing genomes with new data:\n\n";

      foreach my $updateGenome ( @updatedGenomes ) {
	 
	 print "   $updateGenome\n";
      }

      print "\n";
   }

   # Process the additions: just move the directories over.

   foreach my $addGenome ( @addedGenomes ) {
      
      system("cp -R .tempGenomeData/$addGenome .genomeData");
   }
   
   # Process the updates: remove old .gbk files, as well
   # as any files generated from them, and copy the new
   # .gbk files over.
   
   foreach my $updateGenome ( @updatedGenomes ) {
      
      my @updateFiles = keys %{$updatedGenomeHash->{$updateGenome}};

      foreach my $updateFile ( @updateFiles ) {
	 
	 $updateFile =~ /^(.+)\.fna/;

	 my $prefix = $1;

	 my $command = "rm -rf .genomeData/$updateGenome/$prefix\*";

	 system($command);

	 $command = "cp .tempGenomeData/$updateGenome/$prefix\.gbk .genomeData/$updateGenome/$prefix\.gbk";

	 system($command);
      }
   }
}
