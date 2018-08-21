#!/usr/bin/env perl

## build matrix file for all samples using 100000 selected SNPs

use Carp;

## input
my $file_snps_selected  = $ARGV[0];   ## seleted SNPs file from the first step
my $reportfile          = $ARGV[1];   ## finalreport from Genome Studio
my $file_matrix_LRR     = $ARGV[2];   ## output LRR matrix file

## read in selected snps
open(IN, "< $file_snps_selected") or die "Error: can't open snps file $in_snps: $!";
%snps=();
while ($line=<IN>) {
	chomp $line;
	#print "$line\n";
	$snps{$line}++;
}
close IN;

@snps=(keys %snps);
print "total number of SNPs: ".scalar(@snps)."\n";

## parse the header of final report
open(REPORT, "< $reportfile") or die "Error: can't open finalreport $in_file: $!";

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

## parse data part of final report
%samples = (); ## hash for sample ID 
%hash = ();    ## hash of LRR values at selected SNPs for one sample

my $flagsample = 0;    ## indicator of the first sample =0; following samples =1
my $lrrsample = ();    ## tab-delimited LRR values for one sample
my $SampleIDraw = ();  ## temporary sample ID of one sample
my $total = 0;         ## counter of current number of LRR values recorded in $lrrsample
my $flageof = 0;       ## indicaotr of eof =0 not EOF; =1 EOF

while ($line = <REPORT>) {
	$flageof = 1 if eof; ## add file eof flag
	chomp $line;

	@line=split(/\t/, $line);	

	## tansform Log R Ratio
	if (exists($samples{$line[$sample_index]}) && exists($snps{$line[$name_index]})) { ##%snps has been converted to @snps in line 22??

		$lrrvalue = $line[$LRR_index];
		$lrrvalue =~ tr/\015//d;
		$lrrsample = $lrrsample."\t".$lrrvalue;
		$flagsample = 1;
		$SampleIDraw = $line[$sample_index];	
		$total++;

		if ($flageof == 1) {
			$hash{$SampleIDraw} = $lrrsample;
			print "SampleID: $SampleIDraw\t".scalar(keys %samples)."\t$total\n";
			last;
		}

	} elsif (exists($samples{$line[$sample_index]})) {
			
		if ($flageof == 1) {
			$hash{$SampleIDraw} = $lrrsample;
			print "SampleID: $SampleIDraw\t".scalar(keys %samples)."\t$total\n";
			last;
		} else {
			next;
		}
			
	} else {

		if ($flagsample == 0) {
			
			## initialize the first sample
			if (exists($snps{$line[$name_index]})) {
				$samples{$line[$sample_index]}++;
				$lrrvalue = $line[$LRR_index];
				$lrrvalue =~ tr/\015//d;
				$lrrsample = $lrrvalue;
				$total++;
			}
		} elsif ($flagsample == 1) {

			if (exists($snps{$line[$name_index]})) {

				## complete the previous sample
				$hash{$SampleIDraw} = $lrrsample;
				print "SampleID: $SampleIDraw\t".scalar(keys %samples)."\t$total\n";

				## initialize another new sample
				$samples{$line[$sample_index]}++;
				$lrrsample = ();
				$lrrvalue = $line[$LRR_index];
				$lrrvalue =~ tr/\015//d;
				$lrrsample = $lrrvalue;
				$total = 1;
			}
		}
	}
}

close IN;

## save LRR matrix
open(OUT, ">", $file_matrix_LRR) or die "Error: can't open file $file_matrix_LRR: $!";
foreach my $item (keys %hash) {
	print OUT "$item\t$hash{$item}\n";
}
close OUT;
