#!/usr/bin/perl

use strict;

$| = 1;

#######################################################################################################
# 
# Rebuild the BLAST DB.
# 
#######################################################################################################

print "\nAre you sure you want to rebuild the BLAST database (Y/N): ";

my $response = <STDIN>;

chomp $response;

while ( length($response) != 1 or $response !~ /[yYnN]/ ) {
   
   print "Please enter Y or N: ";

   $response = <STDIN>;

   chomp $response;
}

if ( $response =~ /[yY]/ ) {
   
   print "Rebuilding BLAST database...";

   &rebuildBlastDB();

   print "done.\n\n";

} else {
   
   print "\nOkay; exiting without rebuildling BLAST database.\n\n";
}



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

