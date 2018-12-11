#!/usr/bin/perl

## The script was used to run QuantiSNP on Minerva high performance cluster.
## You need to modifiy it according to the system you are using if you would like to use it.
## Please refer to original QuantiSNP documents (https://sites.google.com/site/quantisnp/) for more information 

use strict;
use Getopt::Long;

my $in_dir="";   ## input directory
my $out_dir="";  ## output directory

GetOptions("in_dir=s" => \$in_dir,
		   "out_dir=s" => \$out_dir);

my $out_file=$out_dir."/quantisnp.cnv";

opendir(DIR, $in_dir) or die "cannot open $in_dir: $!"; ## 2018-12-11
open(OUT1, ">", $out_file) or die $!;

my $flag = 1;
while (defined(my $folder = readdir(DIR))) {
	
	next if ($folder=~/^\./);
	##next if ($folder=~/^INTERNAL/);
	##next if ($folder=~/^CONTROL/);
	
	my $filename=$folder.".cnv"; ## $folder is Sample ID
	my $file=$in_dir."/".$folder."/".$filename;
	
	#print "$flag", "$file\n";
	open(IN1, "<$file") or die $!;
	while (my $line=<IN1>) {
		next if ($line=~/^Sample/);
		print OUT1 $line;
	}

	close IN1;
	$flag = $flag + 1;
}

close OUT1;

print "Analysis completed!\n";
