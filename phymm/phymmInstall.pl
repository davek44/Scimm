#!/usr/bin/perl

use strict;

$| = 1;

use Net::FTP;

# Figure out whether to obtain a new copy of the RefSeq genomic database or to use an existing one.

my $newCopy = 0;

&checkDataSource();

# If we have a new (or old but unexpanded) tarball,

if ( $newCopy ) {
   
   # Unzip & untar,

   &expandNewRefSeqData();

   # rename directories if necessary,

   &renameNewRefSeqDirs();
   
   # create FASTA sequences from the GenBank files if necessary,

   &genbankToFasta();
   
   # create training files for the ICMs if necessary,

   &createTrainingFiles();
   
   # generate metafiles describing the RefSeq data, using NCBI's taxonomy data to help,

   &generateRefSeqMetafiles();

   # build the ICM suite,

   &buildICMs();
   
   # build local BLAST database from RefSeq files,

   &generateBlastDB();

} else {
   
   # rename directories if necessary,

   &renameNewRefSeqDirs();
   
   # create FASTA sequences from the GenBank files if necessary,

   &genbankToFasta();
   
   # create training files for the ICMs if necessary,

   &createTrainingFiles();
   
   # generate metafiles describing the RefSeq data, using NCBI's taxonomy data to help,

   &generateRefSeqMetafiles();

   # build the ICM suite,

   &buildICMs();
   
   # build local BLAST database from RefSeq files,

   &generateBlastDB();
}






###########################################################
# SUBROUTINES
###########################################################

sub checkDataSource {
   
   if ( &checkRefSeqData() ) {
      
      print "\nLocal copy of RefSeq genome data detected.\n";

      print "(U)se existing data or (R)eplace with new copy from RefSeq? [U] ";

      my $response = <STDIN>;

      chomp $response;

      while ( length($response) > 1 or ( length($response) == 1 and $response !~ /[uUrR]/ ) ) {
	 
	 print "\nPlease enter U (default: use existing genome data) or R (replace with new copy): ";

	 $response = <STDIN>;

	 chomp $response;
      }
      
      print "\n";
      
      if ( $response =~ /[rR]/ ) {
	 
	 # Grab the newest data from the server.
	 
	 &downloadNewRefSeqData();

	 $newCopy = 1;
	 
      } else {
	 
	 # We're going to use existing data.  If it's just the tarball, though, we want to flag it for expansion.

	 if ( &checkRefSeqData() == 1 ) {
	    
	    $newCopy = 1;
	 }
      }

   } else {
      
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
	 
	 $newCopy = 1;

      } else {
	 
	 die("Can't continue Phymm setup without genomic data; exiting.\n\n");
      }
   }
}

# Check for the existence of a gzipped tarball from the RefSeq FTP server.

sub checkRefSeqData {
   
   # If there's anything in the .genomeData directory, return 2.
   
   # If we have a tarball but nothing in the .genomeData directory, return 1.
   
   # If the tarball is missing and the .genomeData directory is empty, return 0.
   
   if ( not -e './.genomeData' ) {
      
      if ( -e './all.gbk.tar.gz' ) {
	 
	 return 1;

      } else {
	 
	 return 0;
      }

   } else {
      
      # If there's *anything at all* in the '.genomeData' directory, assume it's a genomic library.

      opendir DATA, ".genomeData";

      my @subs = readdir DATA;

      closedir DATA;

      foreach my $sub ( @subs ) {

	 if ( $sub ne '.' and $sub ne '..' ) {
	    
	    return 2;
	 }
      }
      
      if ( -e './all.gbk.tar.gz' ) {
	 
	 return 1;

      } else {
	 
	 return 0;
      }
   }
}

# Replace any existing RefSeq genome data with a fresh copy.

sub downloadNewRefSeqData {
   
   my $ftpSite = 'ftp.ncbi.nih.gov';

   print "Connecting to $ftpSite...";

   my $ftp = Net::FTP->new($ftpSite, Passive => 1, Debug => 0) or die "Can't connect to $ftpSite: $@";

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

# Unzip/untar the downloaded genome-data file.

sub expandNewRefSeqData {
   
   print "Expanding downloaded genome data (this'll take a while)...";

   system("tar -zxvf all.gbk.tar.gz -C .genomeData >& /dev/null");
   
   print "done.\n\nRemoving compressed genome data...";
   
   system("rm -rf all.gbk.tar.gz");
   
   print "done.\n\n";
}

# Standardize the names of the RefSeq organism directories to match the names given in the .gbk files.

sub renameNewRefSeqDirs {
   
   print "Standardizing organism directory names...";
   
   my $dataDir = '.genomeData';

   opendir DATA, "$dataDir" or die("Can't open $dataDir for scanning.\n");

   my @subs = readdir DATA;

   closedir DATA;

   foreach my $sub ( @subs ) {
      
      my $orgName = '';
      
      if ( $sub ne '.' and $sub ne '..' ) {
	 
	 opendir SUB, "$dataDir/$sub" or die("Can't open $dataDir/$sub for scanning.\n");

	 my @files = readdir SUB;

	 closedir SUB;

	 foreach my $file ( @files ) {
	    
	    if ( $file =~ /\%$/ ) {
	       
	       # This file's been obsolesced by NCBI, but for whatever reason remained in the distribution.  Kill it.
	       
	       system("rm $dataDir/$sub/$file");
	    }

	    # We only want accessions beginning with NC_ - all others are in varying stages of unreadiness.
	    
	    # Also, NC_008265 is a phage and not a bacterium; delete it.

	    if ( $file =~ /\.gbk$/ and ( $file !~ /NC_/ or $file =~ /NC_008265/ ) ) {
	       
	       system("rm $dataDir/$sub/$file");
	       
	    } elsif ( $file =~ /\.gbk$/ ) {
	       
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
	       
	       # Fix GenBank's known inconsistencies.
	       
	       if ( $currentOrgName =~ /^Nostoc\s+sp\.\s+PCC\s+7120/ ) {
		  
		  $orgName = 'Nostoc sp. PCC 7120';
		  
	       } elsif ( $currentOrgName =~ /^Synechocystis\s+sp\.\s+PCC\s+6803/ ) {
		  
		  $orgName = 'Synechocystis sp. PCC 6803';

	       } elsif ( $currentOrgName =~ /^Buchnera\s+aphidicola/ and $currentOrgName =~ /cedri/ ) {
		  
		  $orgName = 'Buchnera aphidicola str. Cc (Cinara cedri)';

	       } elsif ( $sub eq 'Bifidobacterium_longum_DJO10A' ) {
		  
		  $orgName = 'Bifidobacterium longum DJO10A';
		  
	       } elsif ( $currentOrgName =~ /^Bacillus\s+licheniformis\s+ATCC\s+14580/ ) {
		  
		  $orgName = 'Bacillus licheniformis ATCC 14580';

	       } elsif ( $currentOrgName =~ /^Alkalilimnicola\s+ehrlichei/ ) {
		  
		  $orgName = 'Alkalilimnicola ehrlichii';

	       } elsif ( $currentOrgName =~ /^Silicibacter\s+pomeroyi/ ) {
		  
		  $orgName = 'Ruegeria pomeroyi';
		  
	       } elsif ( $currentOrgName =~ /^Acidovorax/ and $currentOrgName =~ /citrulli/ ) {
		  
		  $orgName = 'Acidovorax citrulli AAC00-1';
		  
	       } elsif ( $currentOrgName =~ /^Candidatus\s+Methanosphaerula\s+palustris/ ) {
		  
		  $orgName = 'Methanosphaerula palustris E1-9c';
		  
	       } elsif ( $currentOrgName =~ /^Corynebacterium\s+aurimucosum/ ) {
		  
		  $orgName = 'Corynebacterium aurimucosum ATCC 700975';
		  
	       } elsif ( $currentOrgName =~ /^Yersinia\s+pestis\s+KIM/ ) {
		  
		  $orgName = 'Yersinia pestis KIM D27';
		  
	       } elsif ( $orgName eq '' ) {
		  
		  $orgName = $currentOrgName;
	       
	       } elsif ( $orgName ne $currentOrgName ) {
		  
		  print "\n\n   WARNING: Different source organisms given for different .gbk files in data directory \"$dataDir/$sub\".\n\n";
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
	    
	 } else {
	    
	    # If we're not dealing with a completely alien directory - i.e., if we
	    # successfully found an $orgName - process the directory, otherwise ignore.

	    if ( $orgName ne '' ) {
	       
	       # Rename the current subdirectory to a standard name based on the species name quoted in its .gbk files.

	       my $newSubName = $orgName;
	    
	       $newSubName =~ s/\_/UNDERSCORE/g;

	       $newSubName =~ s/\s+/\_/g;

	       $newSubName =~ s/\//SLASH/g;

	       $newSubName =~ s/\(/LPAREN/g;

	       $newSubName =~ s/\)/RPAREN/g;

	       $newSubName =~ s/\'/SINGLEQUOTE/g;

	       $newSubName =~ s/\"/DOUBLEQUOTE/g;
	    
	       $newSubName =~ s/\:/COLONCOLON/g;

	       if ( $newSubName ne $sub and not -e "$dataDir/$newSubName" ) {
	       
		  # We've made a change to the name, and there's not yet a
		  # directory labeled with the new name.  Just move the old one over.
	       
		  my $tempSub = $sub;

		  $tempSub =~ s/\(/\\\(/g;

		  $tempSub =~ s/\)/\\\)/g;

		  my $command = "mv $dataDir/$tempSub $dataDir/$newSubName\n";

		  system($command);

	       } elsif ( $newSubName ne $sub ) {
	       
		  # We've made a change, but the target directory (labeled with the
		  # new name) already exists; move the files from the old dir to
		  # the new one.
	       
		  my $tempSub = $sub;

		  $tempSub =~ s/\(/\\\(/g;

		  $tempSub =~ s/\)/\\\)/g;

		  my $command = "mv $dataDir/$tempSub/* $dataDir/$newSubName\n";

		  system($command);
	       
		  # Delete the now-empty directory labeled with the old name.
	       
		  $command = "rm -rf $dataDir/$tempSub\n";

		  system($command);
	       }

	    } # end ( check to discover whether or not we've been able to establish an $orgName )
	    
	 } # end ( check to ensure current subdirectory isn't empty )
	 
      } # end ( filter on '.' and '..' )
      
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

   foreach my $sub ( @subs ) {
      
      my $defLine = '';
      
      my $orgName = '';
      
      if ( $sub ne '.' and $sub ne '..' ) {
	 
	 opendir SUB, "$dataDir/$sub" or die("Can't open $dataDir/$sub for scanning.\n");

	 my @files = readdir SUB;

	 closedir SUB;

	 foreach my $file ( @files ) {
	    
	    if ( $file =~ /(.+)\.gbk$/ ) {
	       
	       my $prefix = $1;
	       
	       my $fullFile = "$dataDir/$sub/$file";

	       open IN, "<$fullFile" or die("Can't open $fullFile for reading.\n");
	       
	       my $found = 0;
	       
	       $orgName = '';
	       
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
	       
	       $recording = 0;
	       
	       $found = 0;

	       while ( not $found and my $line = <IN> ) {
		  
		  chomp $line;

		  if ( $recording and $line !~ /^\s*ORGANISM/ ) {
		     
		     $line =~ /^\s+(.*)/;

		     $orgName .= ' ' . $1;

		  } elsif ( $line =~ /^SOURCE\s+(.*)/ ) {
		     
		     $orgName = $1;

		     $recording = 1;
		     
		  } elsif ( $line =~ /^\s*ORGANISM/ ) {
		     
		     $found = 1;
		  }
	       }
	       
	       close IN;
	       
	       # Exceptions for poorly-curated GenBank files.
	       
	       if ( $orgName =~ /^Nostoc\s+sp\.\s+PCC\s+7120/ ) {
		  
		  $orgName = 'Nostoc sp. PCC 7120';
		  
	       } elsif ( $orgName =~ /^Synechocystis\s+sp\.\s+PCC\s+6803/ ) {
		  
		  $orgName = 'Synechocystis sp. PCC 6803';

	       } elsif ( $orgName =~ /^Buchnera\s+aphidicola/ and $orgName =~ /cedri/ ) {
		  
		  $orgName = 'Buchnera aphidicola str. Cc (Cinara cedri)';

	       } elsif ( $orgName =~ /^Silicibacter\s+pomeroyi/ ) {
		  
		  $orgName = 'Ruegeria pomeroyi';

	       } elsif ( $sub eq 'Bifidobacterium_longum_DJO10A' ) {
		  
		  $orgName = 'Bifidobacterium longum DJO10A';
	       
	       } elsif ( $orgName =~ /^Alkalilimnicola\s+ehrlichei/ ) {
		  
		  $orgName = 'Alkalilimnicola ehrlichii';

	       } elsif ( $orgName =~ /^Silicibacter\s+pomeroyi/ ) {
		  
		  $orgName = 'Ruegeria pomeroyi';
		  
	       } elsif ( $orgName =~ /^Acidovorax/ and $orgName =~ /citrulli/ ) {
		  
		  $orgName = 'Acidovorax citrulli AAC00-1';
		  
	       } elsif ( $orgName =~ /^Candidatus\s+Methanosphaerula\s+palustris/ ) {
		  
		  $orgName = 'Methanosphaerula palustris E1-9c';
		  
	       }

	       $species->{$prefix} = $orgName;

	       $definition->{$prefix} = $defLine;
	       
	    } # end ( filename check for /\.gbk$/ )
	    
	 } # end ( file iterator on current subdirectory )
	 
      } # end ( filter on '.' and '..' )
      
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

      my $newSubName = $orgName;

      $newSubName =~ s/\_/UNDERSCORE/g;

      $newSubName =~ s/\s+/\_/g;

      $newSubName =~ s/\//SLASH/g;

      $newSubName =~ s/\(/LPAREN/g;

      $newSubName =~ s/\)/RPAREN/g;

      $newSubName =~ s/\'/SINGLEQUOTE/g;

      $newSubName =~ s/\"/DOUBLEQUOTE/g;

      $newSubName =~ s/\:/COLONCOLON/g;
      
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
      
   system("./.scripts/genbankToFasta.pl .genomeData");

   print "done.\n\n";
}

# Create training files for the ICMs if they don't already exist.

sub createTrainingFiles {
   
   print "Creating training files for ICMs as necessary...";
      
   system("./.scripts/createTrainingFiles.pl .genomeData 500");

   print "done.\n\n";
}

# Build the ICMs which form the core of Phymm's classification system

sub buildICMs {
   
   # First, check to see if the C code from Glimmer has been untarred and/or compiled yet.  Do whatever needs doing.

   if ( not -e ".scripts/.icmCode/LICENSE" ) {
      
      # The package hasn't been uncompressed yet.  Do so.
      
      print "Uncompressing Glimmer code...";
      
      system("tar -zxvf .scripts/.icmCode/glimmer3_plusSS.tgz -C .scripts/.icmCode >& /dev/null");
      
      print "done.\n\n";
   }

   if ( not -e ".scripts/.icmCode/bin/build-icm" ) {
      
      # The code hasn't been compiled yet.  Do so.
      
      print "Compiling Glimmer code...";
      
      chdir(".scripts/.icmCode/SimpleMake");

      system("make >& /dev/null");

      chdir("../../..");

      print "done.\n\n";
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
      
      if ( -d "$root/$sub" and $sub ne '.' and $sub ne '..' ) {
	 
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
      
      # Standalone BLAST is installed.  Do we need to regenerate the BLAST DB?  Just check to see whether or not something's there.
      
      if ( not -e ".blastData/bacteriaAndArchaea.nsq" ) {
	 
	 # Concatenate the genome files and build a local BLAST DB from them.
	 
	 print "Generating local BLAST database from RefSeq data...";
      
	 # First, grab a list of relative paths to all the genomic FASTA files.
      
	 my @fastaFiles = ();
      
	 opendir DOT, '.genomeData' or die("Can't open .genomeData for scanning.\n");

	 my @subs = readdir DOT;

	 closedir DOT;

	 foreach my $sub ( @subs ) {
	 
	    if ( -d ".genomeData/$sub" and $sub ne '.' and $sub ne '..' ) {
	    
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
      
	 # Concatenate all the FASTA files together for database creation.
      
	 my $command = 'cat ';

	 foreach my $file ( @fastaFiles ) {
	 
	    $command .= "$file ";
	 }

	 $command .= "> blastDB_data.txt";

	 system($command);
      
	 # Create the local database copy.

	 system("formatdb -i blastDB_data.txt -p F -n bacteriaAndArchaea");

	 system("mv bacteriaAndArchaea.* .blastData");

	 system("rm blastDB_data.txt");

	 system("rm formatdb.log");
      
	 print "done.\n\n";

      } # end if ( there's already a BLAST database which has been built )

   } # end if ( local copy of BLAST binary is installed )

}
