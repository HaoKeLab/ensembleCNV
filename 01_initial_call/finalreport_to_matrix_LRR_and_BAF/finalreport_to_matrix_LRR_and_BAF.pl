#!/usr/bin/env perl

use Data::Dumper;
use Carp;

## This script is used to generate chromosome-wise LRR and BAF matrices from GenomeStudio finalreport.
## The LRR and BAF matrices are used in CNV genotyping.
## Finalreport from GenomeStudio is supposed to include the following columns: 
## "Sample ID", "Chr", "Position", "SNP Name", "Log R Ratio", "B Allele Freq"

## module load perl5

## input
my $reportfile = $ARGV[0];     ## finalreport from Genome Studio
my $path_output = $ARGV[1];    ## path to save results (LRR and BAF chr-matrix)

print "Finalreport:", $reportfile, "\n";
print "Path_output:", $path_output, "\n";

## mkdir LRR folder
mkdir join('/', $path_output, "LRR");
mkdir join('/', $path_output, "BAF");

my ($count_line, $name_index, $sample_index, $LRR_index, $BAF_index, $chr_index, $position_index) = 0;	
my (@field);
open(REPORT, $reportfile) or confess "Error: cannot read from input report file $reportfile: $!\n";
while (<REPORT>) {
	$count_line++;
	m/^\[Data\]/ and last;
	$count_line > 1000 and confess "Error: after reading 1000 lines in $reportfile, still cannot find [Data] section. The $reportfile file may not be in Illumina report format.\n";
}

$_ = <REPORT>;
s/[\r\n]+$//;
$count_line++;
@field = split (/\t/, $_);
@field >= 6 or confess "Error: invalid header line (at least 6 tab-delimited fields, including 'SNP Name', 'Sample ID', 'B Allele Freq', 'Log R Ratio', 'Chr', 'Position' expected) in report file $reportfile: <$_>\n";

for my $i (0 .. @field-1) {
	$field[$i] eq 'SNP Name' and $name_index = $i;
	$field[$i] eq 'Sample ID' and $sample_index = $i;
	$field[$i] eq 'B Allele Freq' and $BAF_index = $i;
	$field[$i] eq 'Log R Ratio' and $LRR_index = $i;
	$field[$i] eq 'Chr' and $chr_index = $i;
	$field[$i] eq 'Position' and $position_index = $i;
}


defined $name_index or confess "Error: the 'SNP Name' field is not found in header line in report file $reportfile: <$_>\n";
defined $sample_index or confess "Error: the 'Sample ID' field is not found in header line in report file $reportfile: <$_>\n";
defined $BAF_index or confess "Error: the 'B Allele Freq' field is not found in header line in report file $reportfile: <$_>\n";
defined $LRR_index or confess "Error: the 'Log R Ratio' field is not found in header line in report file $reportfile: <$_>\n";
defined $position_index or confess "Error: the 'Position' field is not found in header line in report file $reportfile: <$_>\n";
defined $chr_index or confess "Error: the 'Chr' field is not found in header line in report file $reportfile: <$_>\n";

## output file name
my $path_LRR = join("/", $path_output, "LRR/");
my $path_BAF = join("/", $path_output, "BAF/");
my @LRR_ends = ($path_LRR, ".tab");
my @BAF_ends = ($path_BAF, ".tab");
my $file_samples_order = $path_output."/"."samples_order.txt";
my $file_snps_number = $path_output."/"."snps_number.txt";
my $file_snps_name = $path_output."/"."snps_name.txt";
my $file_snps_position = $path_output."/"."snps_position.txt";
my %samples = ();
my $sample_order = 1;
my $sample_before = ();

my %handle = ();
for (1..22) {
	$handle{$_} = ();
}

my %handle_LRR = ();
my %handle_BAF = ();
for (1..22) {
	$handle_LRR{$_} = ();
	$handle_BAF{$_} = ();
}

my %handle_save_LRR = ();
my %handle_save_BAF = ();

my %snp_chr_number = ();
for (1..22) {
	$snp_chr_number{$_} = ();
}

my %snp_chr_name = ();
for (1..22) {
	$snp_chr_name{$_} = ();
}

my %snp_chr_position = ();  
my $flag_snp = "init";

my $flag_sampleID = "old";
my $flag_snp_save = "yes";

while($line = <REPORT>) {

	if ( eof ) {

		chomp $line;

		@line = split(/\t/, $line);
		$sample1 = $line[$sample_index];
		$chr1 = $line[$chr_index];
		$position1 = $line[$position_index];
		$snp1 = $line[$name_index];
		$LRR1 = $line[$LRR_index];
		$BAF1 = $line[$BAF_index];

		if ( exists $handle_LRR{$chr1} ) {
			$handle_LRR{$chr1}{$snp1} = $LRR1;
			$handle_BAF{$chr1}{$snp1} = $BAF1;
		}

		foreach my $item (sort keys %handle) {
			if (! exists $handle_save_LRR{$item}) {

				$file_LRR = join($item, @LRR_ends);
				$file_BAF = join($item, @BAF_ends);

				open my $fh_LRR, q{>>}, join($item, @LRR_ends);
				open my $fh_BAF, q{>>}, join($item, @BAF_ends);

				$handle_save_LRR{$item} = $fh_LRR;
				$handle_save_BAF{$item} = $fh_BAF;
			}

			my @array_chr_snp = ();
			my @array_chr_LRR = ();
			my @array_chr_BAF = ();

			foreach my $snpname (sort keys $handle_LRR{$item}) {
				push @array_chr_LRR, $handle_LRR{$item}{$snpname};
				push @array_chr_BAF, $handle_BAF{$item}{$snpname};
				push @array_chr_snp, $snpname;
			}

			$combine_chr_snp = join("___", @array_chr_snp);
						
			die "chr $item snp number are not same $snp_chr_number{$item}" unless (scalar(@array_chr_snp) == $snp_chr_number{$item});
			die "chr $item snp name are not same" unless ($combine_chr_snp eq $snp_chr_name{$item}); 

			my $chr_LRR = join("\t", $sample_before, join("\t", @array_chr_LRR));
			my $chr_BAF = join("\t", $sample_before, join("\t", @array_chr_BAF));

			print { $handle_save_LRR{$item} } qq($chr_LRR\n);
			print { $handle_save_BAF{$item} } qq($chr_BAF\n);
		}

		last;
	}

	chomp $line;

	@line = split(/\t/, $line);
	$sample1 = $line[$sample_index];
	$chr1 = $line[$chr_index];
	$position1 = $line[$position_index];
	$snp1 = $line[$name_index];
	$LRR1 = $line[$LRR_index];
	$BAF1 = $line[$BAF_index];

	if ((! exists $samples{$sample1}) && ($flag_sampleID eq "old")) {
		$samples{$sample1} = $sample_order;
		$sample_order = $sample_order + 1;
		print "$sample1\n";
		# print "condition_1\n";

		$flag_sampleID = "new";
		if ( exists $handle_LRR{$chr1} ) {
			# $handle{$chr1}{$snp1} = $LRR1;
			$handle_LRR{$chr1}{$snp1} = $LRR1;
			$handle_BAF{$chr1}{$snp1} = $BAF1;

			if ($flag_snp_save eq "yes") {
				$snp_chr_position{$snp1} = $position1;
			}
		}
	} 

	if ( (exists $samples{$sample1}) && ($flag_sampleID eq "new") ) {

		if ( exists $handle_LRR{$chr1} ) {
			# $handle{$chr1}{$snp1} = $LRR1;
			$handle_LRR{$chr1}{$snp1} = $LRR1;
			$handle_BAF{$chr1}{$snp1} = $BAF1;
			$sample_before = $sample1;

			if ($flag_snp_save eq "yes") {
				$snp_chr_position{$snp1} = $position1;
			}
		}

	}

	if ( (! exists $samples{$sample1}) && ($flag_sampleID eq "new") ) {

		foreach my $item (sort keys %handle_LRR) {
			if (! exists $handle_save_LRR{$item}) {

				print join($item, @LRR_ends)."\n";
				print join($item, @BAF_ends)."\n";

				open my $fh_LRR, q{>>}, join($item, @LRR_ends);
				open my $fh_BAF, q{>>}, join($item, @BAF_ends);

				$handle_save_LRR{$item} = $fh_LRR;
				$handle_save_BAF{$item} = $fh_BAF;
			}

			@array_chr_LRR = ();
			@array_chr_BAF = ();
			@array_chr_snp = ();
			foreach my $snpname (sort keys $handle_LRR{$item}) {
				push @array_chr_LRR, $handle_LRR{$item}{$snpname};
				push @array_chr_BAF, $handle_BAF{$item}{$snpname};
				push @array_chr_snp, $snpname;			
			}

			$combine_chr_snp = join("___", @array_chr_snp);

			if ($flag_snp eq  "init") {
				$snp_chr_number{$item} = scalar(@array_chr_snp);
				$snp_chr_name{$item} = $combine_chr_snp;
			} 
			
			print "$item\t".scalar(@array_chr_snp)."\t".scalar(keys %handle)."\n";
			die "chr $item snp number are not same $snp_chr_number{$item}" unless (scalar(@array_chr_snp) == $snp_chr_number{$item});
			die "chr $item snp name are not same" unless ($combine_chr_snp eq $snp_chr_name{$item}); 

			my $chr_LRR = join("\t", $sample_before, join("\t", @array_chr_LRR));
			my $chr_BAF = join("\t", $sample_before, join("\t", @array_chr_BAF));

			print { $handle_save_LRR{$item} } qq($chr_LRR\n);
			print { $handle_save_BAF{$item} } qq($chr_BAF\n);
		}

		if ($flag_snp_save eq "yes") {

			open($fh_snp_number, ">$file_snps_number");
			foreach my $item_chr (keys %snp_chr_number) {
				print $fh_snp_number "$item_chr\t$snp_chr_number{$item_chr}\n";
			}
			close $fh_snp_number;

			open($fh_snp_name, ">$file_snps_name");
			foreach my $item_chr (keys %snp_chr_name) {
				print $fh_snp_name "$item_chr\t$snp_chr_name{$item_chr}\n";
			}
			close $fh_snp_name;

			open($fh_snp_position, ">$file_snps_position");
			foreach my $item_snp (keys %snp_chr_position) {
				print $fh_snp_position "$item_snp\t$snp_chr_position{$item_snp}\n";
			}
			close $fh_snp_position;

			$flag_snp_save = "no";
		}
		$flag_snp = "pass";  ## use the first sample information
		my %handle_LRR = ();
		my %handle_BAF = ();
		for (1..22) {
			$handle_LRR{$_} = ();
			$handle_BAF{$_} = ();
		}

		$samples{$sample1} = $sample_order;
		$sample_order  = $sample_order + 1;
		print "$sample1\n";
		if ( exists $handle_LRR{$chr1} ) {
			# $handle{$chr1}{$snp1} = $LRR1;

			$handle_LRR{$chr1}{$snp1} = $LRR1;
			$handle_BAF{$chr1}{$snp1} = $BAF1;
		}
	}
}

map { close $handle_save_LRR{$_} } keys %handle_save_LRR;
map { close $handle_save_BAF{$_} } keys %handle_save_BAF;

open(OUT, ">$file_samples_order");

foreach my $item (keys %samples) {
	print OUT "$item\t$samples{$item}\n";
}

close OUT;

