#!/usr/bin/perl
use strict;
use warnings;
#unhash line below and edit if you don't have new enough bioperl and for some reason
#don't want to upgrade your bioperl version
#use lib '/path/to/BioDirectory/on/your/system';
use Bio::DB::Taxonomy;
use Bio::Taxon;
use Getopt::Long;
use Time::Local;
use Date::Calc;
use Net::FTP;

=head1 Name

findSpeciesInGenus.pl

=head1 Synopsis

This script is used to find all the species or subspecies below a certain level of the NCBI taxonomy. Input is a taxon name or a taxonomic id. The output is a tab delimited text file containing three columns. The first column contains the taxonomic ids of the species or subspecies, the second column contains their name and the third column contains whether they are a species or subspecies. 

To do its job, indices of the NCBI taxonomy must be available locally. If this script cannot find them, or the indices you have are more than 180 days old, it will download the taxonomy file from the NCBI, unpack the necessary parts and index the files. (The indexing can take a while to run!)

See the Assumptions section below for more about where your files will download to and where indices are built by default. 

This script is intended as a precursor to tax2Blast.pl, which can take the tab delimited output file as input. It then fetches all the nucleotide and peptides sequences for those organisms. Optionally, it can also create blast databases for these sequences. Note that tax2Blast.pl is intended for use with organisms without vast quantities of sequence data available. Further details about tax2Blast.pl and the script itself are available at http://darwin.nerc-oxford.ac.uk/pgp-wiki/index.php/Database_download_and_formating#Script:_tax2Blast.pl.

=head1 Description and Usage

The following commands all allow you to find all the species or subspecies below a certain level of the NCBI taxonomy. Write these to a tab delimited file with column 1 containing taxonomic ids, column 2 containing taxon names and column three containing "species" or "subspecies". 

findSpeciesInGenus.pl -taxid NCBI_taxid  -outfile outfileName

findSpeciesInGenus.pl -taxon genusName  -outfile outFileName

To run this script, the nodes.dmp and names.dmp files of the NCBI taxonomy must be held locally. If you just wish to update these, or download and index them for the first time, then the script can be run as:

findSpeciesInGenus.pl -updateonly

See the section below for tips regarding this.

The output file from this script can be used as input to the tax2Blast.pl script using the -taxfile option. This will cause all sequences in the ncbi for those organisms to be retrieved. A copy of that script, and documentation about using it, can be found at http://darwin.nerc-oxford.ac.uk/pgp-wiki/index.php/Database_download_and_formating#Script:_tax2Blast.pl


=head2 Before you start...

If you already have the NCBI taxonomy indices locally, then please edit the $idx_dir location within this script.

If this script has not been run successfully on your system before, then check that the directories given in the script for $idx_dir and $downloadDir exist and are writable by you. If they are not, then either change the directory locations these variables are set to, or find a system administrator with permissions to write to these directories and get them to run the script using the -updateonly flag.

If your version of Bioperl is not new enough, you won't have Bio::Taxon. And your version of Bio::DB::Taxonomy will probably be out of date. A good solution is to upgrade your version of Bioperl. If for some reason that is not what you want to do, grab the folder biotax.tar from our site and unpack it somewher suitable. This will give you a folder called Bio. It contains only two files, in the appropriate directory hierarchy, so that you can run this script. Before you run the script, edit the line near the top of the file to give the full path to where the Bio directory you just unpacked is:

use lib '/path/to/BioDirectory/on/your/system';

=head1 Assumptions

This script was written specifically for users of NEBC Bio-Linux (http://nebc.nox.ac.uk/biolinux.html), and thus there are assumptions made about default file locations and what programs area already available on the system. 

For example, if indices need to be built or rebuilt on your system, the NCBI taxonomy file taxdump.tar.gz will be downloaded into /tmp and indices will be built into  /home/db/taxonomy. To change these settings, edit $idx_dir and $downloadDir in the script. 

This script requires a recent version of Bioperl as you must have the Bio::DB::Taxonomy and Bio::Taxon modules installed. 


=cut

#note that the text below is necessary because the official Bioperl release
#is missing this function! It has been put into the CVS in 2007, but is still not in release

sub Bio::Taxon::each_Descendent {
    my ($self) = shift;
    my $db ||= $self->db_handle || return;
    return $db->each_Descendent($self);
}

my $failureMsg =  "\nUsage:\n$0 -taxon <taxonName> -outfile <outfile> \nOR\n$0 -taxid <taxid> -outfile <outfile>\nOR\n$0 -updateonly\n\nIf the -updateonly flag is set, the ncbi taxonomy indices are updated for you, but no search is run.\n\n\n";

if(scalar @ARGV < 1)
{
	die $failureMsg;
}


my $taxon = "NULL";
my $taxid = 0;
my $rank = "species";
my $updateonly;
my $outfile;
my $gettaxid=0; 
my $gettaxon=0; 
my $continue = 0;
my $idx_dir = '/tmp/taxonomy';
my $downloadDir = '/tmp';
my $indexAge = 180;  #number of days after which you should rebuild your taxonomy index

sub checkDirExists;
sub checkDirGood;
sub timeSinceIndexBuild;
sub checkIndices;
sub databaseChecks;
sub unpackNCBITax;
sub downloadNCBITax;

GetOptions (	"taxon:s" =>\$taxon,
				"taxid:i" =>\$taxid,
				"outfile:s" =>\$outfile,
				"updateonly!" =>\$updateonly,
				"rank:s" =>\$rank,
				"idx:s" =>\$idx_dir,
				"download:s" =>\$downloadDir
           );


#add in functions from tax2Blast, but use bioperl functions to run them. (e.g. getting taxid from taxon)
#add in error messages when no children found.

#note storage of taxonomy files in $idx_dir, which is also where indices will be built
my ($nodefile,$namesfile) = ("$idx_dir/nodes.dmp","$idx_dir/names.dmp");


if ($updateonly) 
{
	#do database checks here.
	my $dirExists = checkDirExists($idx_dir);
	my $dirGood = 0;
	if ($dirExists)
	{
		$dirGood = checkDirGood($idx_dir);
	}
	unless ($dirGood) { die "\nDirectory given for NCBI taxonomy indices is $idx_dir, and insufficient privileges prevented its creation.\n\n"; }
	
	#download NCBI taxonomy files to /tmp
	my $downloadFail = downloadNCBITax($downloadDir);
 	if ( $downloadFail > 0 ) {die "Problem downloading taxonomic files: $!\nExiting program.\n";  }

	my $success = unpackNCBITax($downloadDir, $idx_dir, $updateonly);

	if ($success) 
	{
		print "About to build NCBI Taxonomy indices. Please wait....\n\n";
		my $db = new Bio::DB::Taxonomy(-source    => 'flatfile',
                               -nodesfile => $nodefile,
                               -namesfile => $namesfile,
                               -directory => $idx_dir  , 
				-force => '1'
						);
		print "Finished building NCBI Taxonomy indices. Script can now be run in search mode.\n\n";
		exit;
	}
	else  
	{ 
		print "\nDownloading or moving new taxonomy files did not succeed.\n\n"; 
		exit; 
	}   #note that this else should never get called and so is probably redundant.

}

if ($taxid && $taxon) 
   {
     print "\n\nBoth a taxon id and a taxon name were given. Taxonomic id $taxid will be used.\n\n"; 
   } 
elsif (($taxid == 0) && ($taxon eq "NULL" )) 
   {
    print "\nThis script requires a valid taxon id or a taxon name to be given after the appropriate command line flag.\n\n";
    print $failureMsg;
    exit;
   } 


unless ($outfile)
{
        print "\nAn outfile was not specified with the flag -outfile.\nUnless by immediately using Ctrl-c now, an output file, tmp.txt, will be written to the current directory.\n";
        sleep("2");

		#print "\nContinuing now...\n";
        $outfile = "tmp.txt";
}



my ($force, $needsUpdating) = databaseChecks($idx_dir, $indexAge);

if ( ($force == 1) && ( $needsUpdating == 1) ) 
{ 
	print "NCBI taxonomy indices need to be updated. Be patient as this will take some time before the search commences...\n"; 
        #download NCBI taxonomy files to /tmp
        my $downloadFail = downloadNCBITax($downloadDir);
        if ( $downloadFail > 0 ) {die "Encountered problem downloading taxonomic files: $!\nnow exiting.\n";  }

        my $force = unpackNCBITax($downloadDir, $idx_dir, $updateonly); #set to 1 if we got files alright. Otherwise 0. 
}

elsif ( $force == 1 ) 
{  #I don't think we can ever get into this else clause...leaving it for good measure, although bad form.
	print "NCBI taxonomy indices need to be updated. Be patient as this will take some time before the search commences...\n";
	#download NCBI taxonomy files to /tmp
        my $downloadFail = downloadNCBITax($downloadDir);
        if ( $downloadFail > 0 ) {die "Encountered problem downloading taxonomic files: $!\nnow exiting.\n";  }
	my $force = unpackNCBITax($downloadDir, $idx_dir, $updateonly); #set to 1 if we got files alright. Otherwise 0.

	if ($force == 0)  { print "Encountered a problem updating the indices. Previous indices will be used instead.\n"; }
}

##############

#hopefully we have exited by now if there is a problem. 

my $db = new Bio::DB::Taxonomy(	-source    => 'flatfile',
								-nodesfile => $nodefile,
								-namesfile => $namesfile,
								-directory => $idx_dir  , 
								-force => $force
						);
if($taxon)
   { 
	my @taxonids = $db->get_taxonids($taxon);
	my $sizeList = @taxonids;
	if ( $sizeList == 0 ) 
	{
		print "\nNo taxonomic id associated with $taxon.\n\n";
		exit;
	}
	elsif ($sizeList > 1)
	{
		print "\nMore than one taxonomic id encountered for $taxon; only the first will be used: $taxonids[0].\n\n";
		#should do a for loop here to print out what the other ids are.
		print "The extra taxonomic ids are:\n";
		for(my $excess=1;$excess<=$sizeList;$excess++)
		{
			print "$taxonids[$excess]\n" ;
		}
		$taxid = $taxonids[0];
	}
	else
	{
		$taxid = $taxonids[0];
	}
   }

my $node = $db->get_taxon(-taxonid => $taxid);
#print "$node\n";

#print $node->id, " ", $node->scientific_name, " ", $node->rank, "\n";
# to only get children that are of a particular rank in the taxonomy test if their rank is 'species' for example
#my @species= grep { $_->rank eq 'species' } $node->get_all_Descendents;
#my @subspecies = grep { $_->rank eq 'subspecies' } $node->get_all_Descendents;

my @species= grep { $_->rank eq $rank } $node->get_all_Descendents;

my $childSize = @species;   #note that this contains species AND subspecies. See regexp above for why.

print "\nTaxa in $rank is $childSize and are written to $outfile.\n\n";  

open (OUTFILE, ">$outfile")  or die "\nUnable to write to $outfile.\n\n";

for my $child ( @species ) {
    print OUTFILE $child->id, "\t"; # NCBI taxa id
    print OUTFILE $child->scientific_name, "\t"; # scientific name
    print OUTFILE $child->rank, "\n"; # e.g. species
}

close OUTFILE;



################################
sub checkDirExists
{
	my $dir = $_[0];
	my $dirExists = 0;
	my $dirGood = 0;
 	if ( -e $dir )  {  $dirExists = 1;  }

	else
	{
		my $check = mkdir("$dir");
		if ( $check == 0 )
		{
			print "\nFailed to create $dir and so will not be able to build NCBI indices with current settings. Ensure that the function call is made with sufficient privileges.\n\n ";
			exit;
			
		} 
		else 
		{ 
			print "$dir has been created and now stores NCBI taxonomy indices.\n";
			$dirExists = 1; 

		}
	}

	return $dirExists;  #need to catch a return value in updateonly checks. 
}


###############
sub checkDirGood
{
	my $dir = $_[0];
	my $dirGood = 0;
	#presumably if user has permission to write to the index directory then it exists
	if (( -w $dir ) && (-x $dir)) { $dirGood = 1; } ### COMMENTED JME
	#$dirGood = 1;

	return $dirGood;

}

############
sub timeSinceIndexBuild
{
	my $file = $_[0];
	my $fileDate = (stat($file))[9];  #get mtime of file
	my $today = time();

	my ($year1,$month1,$day1) = Date::Calc::Time_to_Date($fileDate);
	my ($year2,$month2,$day2) = Date::Calc::Time_to_Date($today);


	my $indexAge = Date::Calc::Delta_Days($year1,$month1,$day1, $year2,$month2,$day2);

	return($indexAge);
}


#############

sub checkIndices
{
	my $idx_dir = $_[0];
	my $indexAge = $_[1]; #given in days
	my $force = 0;
	my $needUpdating = 0;
	#check for all files. 

	if (( -e "$idx_dir/parents" ) && ( -e "$idx_dir/nodes" ) && ("$idx_dir/names2id") && ("$idx_dir/id2names"))
	{
		my $myIndexAge = timeSinceIndexBuild("$idx_dir/id2names");
		if ($myIndexAge > $indexAge ) 
		{
			print "NCBI taxonomy indices are more than $indexAge days old. Preparing to rebuild...\n";
			my $havePermissions = checkDirGood($idx_dir);
			$needUpdating = 1;
			if ($havePermissions) { $force = 1; }
			else { print "Unable to rebuild the indices to $idx_dir. Ensure that the function call with -updateonly is made with sufficient privileges.\n\n"; }

		}
		#else, $force should remain 0 and search will just run.

	}
	else    #index files don't already exist
	{ 
		my $havePermissions = checkDirGood($idx_dir);
		if ($havePermissions) 
		{ 
			$force = 1; 
			print "NCBI taxonomy indices need to be updated. Be patient as this will take some time...\n"; 
		}
		else 
		{ 
			print "NCBI taxonomy indices need to be updated. Be patient as this will take some time before the search commences...\n\n";
			#
		}  
		
	}  
	

	return ($force, $needUpdating);  #force is 1 when no indices, force is 0 if indices there needsUpdating should be set to 1 if indices are more than $indexAge days old. 
#These are separate because I'd like to add a function in that prompts the user for 
#whether they want to update the indices or not. But right now, I'll set the default 
#to set force to 1 if indices are over six months old. 

}

#############
sub databaseChecks
{
	my $idx_dir = $_[0];
	my $indexAge = $_[1];  #given in days
 	checkDirExists($idx_dir);  #will exit if they can't create this. Don't need to catch a return value.
	my $dirGood = 0;
	

#want to replace error below with nice function that asks user about making the directory


	my ($force,$needUpdating) = checkIndices($idx_dir, $indexAge); 
	if ($needUpdating) { $force = 1; }

	#if force is set to 1, then check dir is good
	if ($force == 1) { $dirGood = checkDirGood($idx_dir); }

	if (( $dirGood == 0) && ($force == 1))  #probably an unnecessary check. Not sure you can have this condition at this point
	{
		print "Unable to rebuild the indices to $idx_dir. Ensure that the function call with -updateonly is made with sufficient privileges.\n\n";
		exit;
	}

	return ($force,$needUpdating);


}

###############

sub downloadNCBITax
{
    my $downloadDir = $_[0];
    my $downloadFail = 1;
    #problem with ftp. It seems to be missing the tail end of the gzip file so uncompressing fails. 
    #works with command line ftp. Try wget:
	if ( -e "$downloadDir/taxdump.tar.gz")
	{
		print "\nRemoving $downloadDir/taxdump.tar.gz...\n";
		unlink("$downloadDir/taxdump.tar.gz");  #should there be a check on success here?
		print "Downloading updated taxdump.tar.gz to $downloadDir...\n";
	}
	$downloadFail = system("wget -P $downloadDir ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz");
	
   # my $ftp;

   # $ftp = Net::FTP->new("ftp.ncbi.nlm.nih.gov", Debug => 0)
   #   or die "Cannot connect to ftp.ncbi.nlm.nih.gov: $@";

   # $ftp->login("anonymous",'-anonymous@')
   #   or die "Cannot login ", $ftp->message;

   # $ftp->cwd("/pub/taxonomy")
   #   or die "Cannot change working directory ", $ftp->message;

   # $ftp->get("taxdump.tar.gz", "$downloadDir/taxdump.tar.gz")
   #   or die "Failed to get taxdump.tar.gz ", $ftp->message;
#	wait;
#    $ftp->quit;

	return $downloadFail;

}

############

sub unpackNCBITax
{
	my $downloadDir = $_[0];
	my $idx_dir = $_[1];
	my $updateonly = $_[2];
	my $success = 0;
#unpack the nodes.dmp and names.dmp files and put them where they're expected
        #puts files into working directory
	#tried the following, but didn't work. (Works on command line.)
	my $nodes = system("gunzip -c $downloadDir/taxdump.tar.gz | tar -C $idx_dir -xf - nodes.dmp");
	wait;
	my $names = system("gunzip -c $downloadDir/taxdump.tar.gz | tar -C $idx_dir -xf - names.dmp");
	wait;
	#new try:
#	system("gunzip $downloadDir/taxdump.tar.gz");
#	system("tar xf $downloadDir/taxdump.tar  nodes.dmp");
	

	if (($nodes == 0) && ($names == 0)) 
	{ 
		#$nodes = system ("mv nodes.dmp $idx_dir\n"); 
		#$names = system ("mv nodes.dmp $idx_dir\n");
		$success = 1;
		print "taxdump.tar.gz file downloaded fine. nodes.dmp and names.dmp unpacked successfully\n";
	} 
	else 
	{ 
		print "Unable to extract nodes.dmp and/or names.dmp file from the taxonomy tar file in $downloadDir\n";
		if ($updateonly) 
		{
			print "\nEncountered problem fetching nodes.dmp and/or names.dmp file.\n\n";
			exit;
		}
		else #if your index was too old for example
		{
			print "\nUnable to rebuild the indices to $idx_dir. Ensure that the function call with -updateonly is made with sufficient privileges to update indices.\n\n"; 
			$success = 0;  #make this into $force = 0 back in the main body of script.
		}
		#ideally should have an elsif clause here that picks up on whether you just need to updat the index or not. If so, it continues running on old indices but reports problem. 
	}	
	return $success;   
	
}