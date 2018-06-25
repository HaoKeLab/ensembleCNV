#!/usr/bin/env perl

use Carp;

$file_snps_selected = "/sc/orga/projects/haok01a/chengh04/MEGA/MegaEX_Inga/analysis_Part1/batch_effect/dat/snps.randomly.select.txt";
$reportfile = "/sc/orga/projects/haok01a/chengh04/MEGA/MegaEX_Inga/FinalReport/Plates_Part1/MegaEX_Inga_Part1_FinalReport.txt";
$file_marix_LRR="/sc/orga/projects/haok01a/chengh04/MEGA/MegaEX_Inga/analysis_Part1/batch_effect/res/matrix.LRR.snps.randomly.select.txt";

## read in selected snps
open(IN, "< $file_snps_selected") or die "Can't open snps file $in_snps: $!";
%snps=();
while ($line=<IN>) {
	chomp $line;
	print "$line\n";
	$snps{$line}++;
}

close IN;
@snps=(keys %snps);
print "total number of snps:".scalar(@snps)."\n";

## build matrix file for all samples using 100000 selected snps
%samples=(); ## for all samples
%hash=();
open(REPORT, "< $reportfile") or die "can't open finalreport $in_file: $!";

my (@field);
my ($count_line, $sample_index, $name_index, $LRR_index) = 0;

while (<REPORT>) {
	$count_line++;
	m/^\[Data\]/ and last;
	$count_line > 1000 and confess "Error: after reading 1000 lines in $reportfile, still cannot find [Data] section. The $reportfile file may not be in Illumina report format.\n";
}

$_ = <REPORT>;
s/[\r\n]+$//;
$count_line++;
@field = split (/\t/, $_);
@field >= 3 or confess confess "Error: invalid header line (at least 3 tab-delimited fields, including 'SNP Name', 'Sample ID', 'Log R Ratio' expected) in report file $reportfile: <$_>\n";

for my $i (0 .. @field-1) {
	$field[$i] eq 'SNP Name' and $name_index = $i;
	$field[$i] eq 'Sample ID' and $sample_index = $i;
	$field[$i] eq 'Log R Ratio' and $LRR_index = $i;
}

defined $name_index or confess "Error: the 'SNP Name' field is not found in header line in report file $reportfile: <$_>\n";
defined $sample_index or confess "Error: the 'Sample ID' field is not found in header line in report file $reportfile: <$_>\n";
defined $LRR_index or confess "Error: the 'Log R Ratio' field is not found in header line in report file $reportfile: <$_>\n";

my $flagsample=0; ##  lrr save flag
my $lrrsample=(); ##
my $SampleIDraw=();
my $total=0;
my $flageof=0;

while ($line = <REPORT>) {
	$flageof=1 if eof; ## add file eof flag
	chomp $line;

	@line=split(/\t/, $line);	

	## tansform Log R Ratio
	if (exists($samples{$line[$sample_index]})&&exists($snps{$line[$name_index]})) {

		$lrrvalue=$line[$LRR_index];
		$lrrvalue=~ tr/\015//d;
		$lrrsample=$lrrsample."\t".$lrrvalue;
		$flagsample=1;
		$SampleIDraw=$line[$sample_index];	
		$total++;

		if ($flageof == 1) {
			$hash{$SampleIDraw}=$lrrsample;
			print "SampleID:$SampleIDraw\t".scalar(keys %samples)."\t$total\n";
			last;
		}

	} elsif (exists($samples{$line[$sample_index]})) {
			
		if ($flageof == 1) {
			$hash{$SampleIDraw}=$lrrsample;
			print "SampleID:$SampleIDraw\t".scalar(keys %samples)."\t$total\n";
			last;
		} else {
			next;
		}
			
	} else {

		if ($flagsample==0) {
			## first sample
			if (exists($snps{$line[$name_index]})) {
				$samples{$line[$sample_index]}++;
				$lrrvalue=$line[$LRR_index];
				$lrrvalue=~ tr/\015//d;
				$lrrsample=$lrrvalue;
				$total++;
			}
		} elsif ($flagsample==1) {
			if (exists($snps{$line[$name_index]})) {

				$hash{$SampleIDraw}=$lrrsample;
				print "SampleID:$SampleIDraw\t".scalar(keys %samples)."\t$total\n";

				$samples{$line[$sample_index]}++;
				$lrrsample=();
				$lrrvalue=$line[$LRR_index];
				$lrrvalue=~ tr/\015//d;
				$lrrsample=$lrrvalue;
				$total=1;
			}
		}

	}

}

close IN;

## save LRR matrix
open(OUT, ">", $file_marix_LRR) or die $!;
foreach my $item (keys %hash) {
	print OUT "$item\t$hash{$item}\n";
}
close OUT;
