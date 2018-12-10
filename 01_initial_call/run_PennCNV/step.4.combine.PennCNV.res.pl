#!/usr/bin/env perl

## The script was used to run PennCNV on Minerva high performance cluster.
## You need to modifiy it according to the system you are using if you would like to use it.
## Please refer to original PennCNV documents (http://penncnv.openbioinformatics.org/en/latest/) for more information

use strict;
use Getopt::Long;

my $in_dir="";   ## input directory
my $out_dir="";  ## output directory

GetOptions("in_dir=s" => \$in_dir,
	   "out_dir=s" => \$out_dir);

my $out_file=$out_dir."/"."CNV.PennCNV.rawcnv";
my $out_log=$out_dir."/"."CNV.PennCNV.log";

opendir(DIR, $in_dir) or die "cannot open $in_dir: $!"; ## 2018-12-10
open(OUT1, ">", $out_file) or die $!;
open(OUT2, ">", $out_log) or die $!;

while (defined(my $folder = readdir(DIR))) {
	
	next if ($folder=~/^\./);

	my $filename=$folder.".rawcnv";
	my $logname=$folder.".log";
	my $file=$in_dir."/".$folder."/".$filename;
	my $logfile=$in_dir."/".$folder."/".$logname;
	## print "$file\n";
	open(IN1, "<$file") or die $!;
	open(IN2, "<$logfile") or die $!;
	while (my $line=<IN1>) {
		print OUT1 $line;
	}

	while (my $line=<IN2>) {
		print OUT2 $line;
	}

	close IN1;
	close IN2;
}

close OUT1;
close OUT2;

print "Analysis completed!\n";
