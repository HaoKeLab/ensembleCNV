#!/usr/bin/env perl

## Acknowlegement: This script was adapted from split_illumina_report.pl in PennCNV package
##                 originally developed by Kai Wang, PhD <kw2701@cumc.columbia.edu>.

use warnings;
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2010-11-28 12:47:21 -0800 (Sun, 28 Nov 2010) $';

our ($verbose, $help, $man);
our ($reportfile, $prefix, $suffix, $numeric_name, $comma, $tolerate, $revised_file);

GetOptions ('verbose'=>\$verbose, 'help'=>\$help, 'man'=>\$man, 'prefix=s'=>\$prefix, 'suffix=s'=>\$suffix, 'numeric_name'=>\$numeric_name, 'comma'=>\$comma, 'tolerate'=>\$tolerate, 
	'revised_file=s'=>\$revised_file) or pod2usage ();

# If used, %revised is populated with content from $revised_file.
our %revised = (); # Key = sample ID in Illumina report; value = sample ID to use.

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error");

($reportfile) = @ARGV;
$prefix ||= '';
$suffix ||= '';


#
# Populate %revised if necessary:

if ( $revised_file ) {
	open (REV, $revised_file) or confess "Error: cannot read file $revised_file :$!\n";
	while ( my $line = <REV> ) {
		chomp($line);
		my ($orig_ID, $revised_ID) = split /\s+/, $line;
		if ( exists $revised{$orig_ID} ) {
			warn "\nOriginal ID $orig_ID already encountered in revised_file; keeping old value.\n";
		}
		else {
			$revised{$orig_ID} = $revised_ID;
		}
	}
}



splitIlluminaReport ($reportfile, $prefix, $suffix, $numeric_name, $comma);



sub splitIlluminaReport {
	my ($reportfile, $prefix, $suffix, $numeric_name, $comma) = @_; ## add chr_index
	my ($count_line, $count_file, $pre, $name_index, $sample_index, $lrr_index, $baf_index, $chr_index) = (0, 0, 'NA');
	my (@field);
	
	open (REPORT, $reportfile) or confess "Error: cannot read from input report file $reportfile: $!\n";
	while (<REPORT>) {
		$count_line++;
		m/^\[Data\]/ and last;
		$count_line > 1000 and confess "Error: after reading 1000 lines in $reportfile, still cannot find [Data] section. The $reportfile file may not be in Illumina report format.\n";
	}
	
	$_ = <REPORT>;
	s/[\r\n]+$//;
	$count_line++;
	@field = $comma ? (split (/,/, $_)) : (split (/\t/, $_));
	@field >= 4 or confess "Error: invalid header line (at least 4 tab-delimited fields, including 'SNP Name', 'Sample ID', 'B Allele Freq', 'Log R Ratio' expected) in report file $reportfile: <$_>\n";
	
	for my $i (0 .. @field-1) {
		$field[$i] eq 'SNP Name' and $name_index = $i;
		$field[$i] eq 'Sample ID' and $sample_index = $i;
		$field[$i] eq 'B Allele Freq' and $baf_index = $i;
		$field[$i] eq 'Log R Ratio' and $lrr_index = $i;
		$field[$i] eq 'Chr' and $chr_index = $i;
	}

	defined $name_index or confess "Error: the 'SNP Name' field is not found in header line in report file $reportfile: <$_>\n";
	defined $sample_index or confess "Error: the 'Sample ID' field is not found in header line in report file $reportfile: <$_>\n";
	defined $baf_index or confess "Error: the 'B Allele Freq' field is not found in header line in report file $reportfile: <$_>\n";
	defined $lrr_index or confess "Error: the 'Log R Ratio' field is not found in header line in report file $reportfile: <$_>\n";
	
	my $got_sample_ID = 0; # Use this so we don't report on multiple adjacent lines

	while (<REPORT>) {
		s/[\r\n]+$//;
		$count_line++;
		@field = $comma ? (split (/,/, $_)) : (split (/\t/, $_));
		@field >= 4 or confess "Error: invalid data line (at least 4 tab- or comma-delimited fields expected) in report file $reportfile: <$_>\n";
		defined $field[$name_index] or confess "Error: the 'SNP Name' field is not found in data line in report file $reportfile: <$_>\n";

		next if ($field[$chr_index]=~/^0|XY|Y|MT/); ##/X|Y|MT/
		#defined $field[$sample_index] or confess "Error: the 'Sample ID' field is not found in data line in report file $reportfile line $.; line is:\n<$_>\n";
	
		# Sometimes a "blank" sampleID can be produced in BeadStuidio files
		# (i,.e. just whitespace). Be resilient to this.

		if ( $field[$sample_index] !~ /\S+/ ) {
		  if ( $got_sample_ID == 0 ) {
			# already reported this; keep quiet.
		  }
		  else { # Just found this bad'un.
			$got_sample_ID = 0;
			print "\n***************\nWARNING: No sample ID found at line $.; line is:\n<$_>\nSkipping...\n***************\n\n";
		  }
		  next;
		}
		else {
		  $got_sample_ID = 1; # OK!
		}

		if (not defined $field[$baf_index]) {
			if ($tolerate) {
				print STDERR "WARNING: Skipping marker $field[$name_index] for $field[$sample_index] due to lack of BAF information\n";
				next;
			} else {
				confess "Error: the 'B Allele Freq' field is not found in data line in report file $reportfile: <$_>\n";
			}
		}
		if (not defined $field[$lrr_index]) {
			if ($tolerate) {
				print STDERR "WARNING: Skipping marker $field[$name_index] for $field[$sample_index] due to lack of LRR information\n";
				next;
			} else {
				confess "Error: the 'Log R Ratio' field is not found in data line in report file $reportfile: <$_>\n";
				next;
			}
		}
		
		if ($field[$sample_index] eq $pre) {
			print OUT join ("\t", @field[$name_index, $lrr_index, $baf_index]), "\n";
   		} else {
			$count_file++;
			my $outname;
			my $revised_sample;
			if ($numeric_name) {
				$outname = "$prefix"."split$count_file$suffix";
			} else {
				$outname = "$prefix$field[$sample_index]$suffix";
				if ($revised_file) {
					my $key = $field[$sample_index];
					$revised_sample = $revised{$key};
					$outname = "$prefix$revised_sample$suffix";
					# print "\noutname is $outname.";
				}
			}
			print STDERR "NOTICE: Writing to output signal file $outname\n";
			open (OUT, ">$outname") or confess "Error: cannot write to output file $outname:$!\n";
			if ( $revised_file ) {
				print OUT "Name\t$revised_sample.Log R Ratio\t$revised_sample.B Allele Freq\n";
			}
			else {
				print OUT "Name\t$field[$sample_index].Log R Ratio\t$field[$sample_index].B Allele Freq\n";
			}
			print OUT join("\t",@field[$name_index, $lrr_index, $baf_index]), "\n";
		}
		$pre = $field[$sample_index];
	}
	print STDERR "NOTICE: Finished processing $count_line lines in report file $reportfile, and generated $count_file output signal intensity files\n";
}


=head1 SYNOPSIS
 finalreport_to_iPattern.pl [arguments] <reportfile>
 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	-p, --prefix <string>		prefix of output file name
 	-s, --suffix <string>		suffix of output file name
 	-n, --numeric_name		use numeric file name (default: Sample ID is file name)
 	-c, --comma			fields are comma-delimited (default: fields are tab-delimited)
 	-t, --tolerate			tolerate records without LRR/BAF information
 	-r, --revised_file <file>	path to "revised" file of alternate sample IDs
 	    --tolerate			tolerate when LRR/BAF do not exist in input file
 Function: split the Illumina report file to individual signal intensity files
 Example: finalreport_to_iPattern.pl -prefix signal/ -suffix .txt HapMap.report
          finalreport_to_iPattern.pl -num HapMap.report
 Version: $LastChangedDate: 2018-02-06 11:44:00 -0800 (Mon, 06 Feb 2018) $
=head1 OPTIONS
=over 8
=item B<--help>
print a brief help message and exit
=item B<--man>
print the complete manual of how to use the program
=item B<--verbose>
use verbose output
=item B<--prefix>
specify the prefix for output file name. By default the file name is the "Sample 
ID" field in the report file.
=item B<--suffix>
specify the suffix for output file name. By default the file name is the "Sample 
ID" field in the report file.
=item B<--numeric_name>
specify that the output file name should be arranged in a numeric manner (such 
as split1, split2, split3, etc.).
=item B<--comma>
specify that the data fields in report file is comma-delimited. By default, the 
tab-delimited format is assumed.
=item B<--tolerate>
tolerate the problem when Log R Ratio (LRR) or B Allele Frequency (BAF) values 
do not exist in input file. In this case, the program will still run, but the 
output will not be immediately useful for CNV inference. For example, sometimes 
one may get a Report file with X and Y values, and can generate splitted 
individual files with X/Y values, and then use a simple script to convert the 
X/Y to LRR/BAF for each individual.
=item B<--revised_file>
specify path to file listing sampleIDs to use instead of those appearing in
the Illumina report file. The revised file must have 2 whitespace-separated
columns, without a header; the first column is the original IDs (i.e. those
in the Illumina report file), and the second the corresponding (revised) ID
to use instead. You must provide one such mapping for each sampleID in the
Illumina file.
=back
=head1 DESCRIPTION
This program is used to convert the "Final Report" file from Illumina GenomeStudio 
software to a format that can be used by the iPattern program for CNV 
detection.
The report file should contain at least ten tab-delimited fields, including 
'SNP Name', 'Sample ID', 'B Allele Freq', 'Log R Ratio', 'Chr', 'Position', 
'Allele1 - Forward', 'Allele1 - Forward',  'X', 'Y', (To run PennCNV, only the first 
4 fields are needed, but the last 6 are required to run iPattern). If not, try to generate 
the report file again in GenomeStudio, and make sure that these fields are 
selected to be exported.
It is okay for the final report files to contain additional columns, but they will be 
ignored and not used in the output files. For example, the first a few lines of 
a typical Illumina final report file is presented below:

	[Header]
	GSGT Version    2.0.3
	Processing Date 9/26/2017 10:45 AM
	Content         MEGA_Consortium_v2_15070954_A2.bpm
	Num SNPs        2036060
	Total SNPs      2036060
	Num Samples     1859
	Total Samples   1859
	[Data]
	Sample ID       SNP Name        Chr     Position        Allele1 - Forward       Allele2 - Forward       X       Y      B Allele Freq   Log R Ratio
	EFZ80734A       1:10001102-G-T  1       10001102        T       T       0.595   0.035   0.0164  -0.1132
	EFZ80734A       1:100011159-T-G 1       100011159       G       G       0.034   0.912   1.0000  0.0512
	EFZ80734A       1:10002775-GA   1       10002775        G       G       0.013   1.298   1.0000  0.0318
	EFZ80734A       1:100122796-C-T 1       100122796       T       T       1.947   0.000   0.0000  0.1345
	EFZ80734A       1:100152282-CT  1       100152282       C       C       0.077   1.469   0.9968  0.0205
	EFZ80734A       1:100154376-GA  1       100154376       G       G       0.014   0.838   1.0000  0.0918
	EFZ80734A       1:100154844-CA  1       100154844       C       C       0.097   1.599   0.9970  0.0027
	EFZ80734A       1:100155035-AC  1       100155035       A       A       1.812   0.038   0.0000  0.0425
	EFZ80734A       1:100155084-CT  1       100155084       C       C       0.078   1.680   0.9912  0.0309
	EFZ80734A       1:100182985-CA  1       100182985       C       C       0.022   1.652   1.0000  0.0928

Each line in the file include measurements for one marker for one sample. The 
'SNP Name', 'Sample ID', 'B Allele Freq', 'Log R Ratio' fields will be examined 
and be written to output files.
=over 8
=item * B<Explanation of the -revised_file argument>
Sometimes users may want to use a different sample ID in the splitted output 
file, and the --revised_file argument can be useful. The revised file must have 
2 whitespace-separated columns, without a header; the first column is the 
original IDs (i.e. those in the Illumina final report file), and the second the 
corresponding (revised) ID to use instead. For example, given an Illumina file 
with lines such as:
	...
	EFZ80734A       1:100155084-CT  1       100155084       C       C       0.078   1.680   0.9912  0.0309
	EFZ80734A       1:100182985-CA  1       100182985       C       C       0.022   1.652   1.0000  0.0928
	...
by default the original program will create a file called "EFZ80734A", with header line:
	Name    EFZ80734A.Log R Ratio   EFZ80734A.B Allele Freq
However, alternative ID can be used in the output. Suppose this is 
internal_ID_123. Then, given a file with entries such as:
	...
	EFZ80734A   internal_ID_123
	...
One can run the script with the -r option and get a file called "internal_ID_123", with header:
	Name    internal_ID_123.Log R Ratio   internal_ID_123.B Allele Freq
=item * B<File format issues>
It is possible for Illumina final call report files to be output such that they 
are grouped by SNP and not sample. This means that the results from different 
samples would be interleaved with the results from other samples. Currently the 
script cannot cope with this; all it will do is to output the data from the last 
occurrence of each sample.
=back
This script was adapted from split_illumina_report.pl in PennCNV package
originally developed by Kai Wang, PhD <kw2701@cumc.columbia.edu>.
For questions, comments or bug reports, please contact us at 
Haoxiang Cheng <haoxiang.cheng@mssm.edu>, or
Zhongyang Zhang <zhongyang.zhang@mssm.edu>, or
Ke Hao <ke.hao@mssm.edu>.
=cut