#!/usr/bin/env perl

## The script was used to run PennCNV on Minerva high performance cluster.
## You need to modifiy it according to the system you are using if you would like to use it.
## Please refer to original PennCNV documents (http://penncnv.openbioinformatics.org/en/latest/) for more information

use Getopt::Long;

$in_dir="";
$out_dir="";

GetOptions("in_dir=s" => \$in_dir,
		   "out_dir=s" => \$out_dir);

$out_file=$out_dir."CNV.PennCNV.rawcnv";
$out_log=$out_dir."CNV.PennCNV.log";

opendir(DIR, $in_dir) or "cannot open $in_dir: $!";
open(OUT1, ">", $out_file) or die $!;
open(OUT2, ">", $out_log) or die $!;

while (defined($folder = readdir(DIR))) {
	
	next if ($folder=~/^\./);

	$filename=$folder.".rawcnv";
	$logname=$folder.".log";
	$file=$in_dir.$folder."/".$filename;
	$logfile=$in_dir.$folder."/".$logname;
	print "$file\n";
	open(IN1, "<$file") or die $!;
	open(IN2, "<$logfile") or die $!;
	while ($line=<IN1>) {
		print OUT1 $line;
	}

	while ($line=<IN2>) {
		print OUT2 $line;
	}

	close IN1;
	close IN2;
}

close OUT1;
close OUT2;
