#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Carp;

## This script is used to generate chromosome-wise LRR and BAF matrices from GenomeStudio finalreport.
## The LRR and BAF matrices are used in CNV genotyping.
## Finalreport from GenomeStudio is supposed to include the following columns: 
## "Sample ID", "Chr", "Position", "SNP Name", "Log R Ratio", "B Allele Freq", required by PennCNV, QuantiSNP, and iPattern
## ("Allele1 - Forward and Allele2 - Forward") or ("Allele1 - Top" and "Allele2 - Top"), "X" and "Y", required by iPattern, 
## report warnings if these 4 fields are not found

## input
my $reportfile = $ARGV[0];     ## finalreport from Genome Studio
my $path_output = $ARGV[1];    ## path to save results (LRR and BAF chr-matrix)

print "Finalreport:", $reportfile, "\n";
print "Path_output:", $path_output, "\n";

my @warnings = ();

## haoxiang 2018-11-26
my ($count_line, $name_index, $sample_index, $LRR_index, $BAF_index, $chr_index, $position_index, $X_index, $Y_index, $Allele1Forward_index, $Allele2Forward_index, $Allele1Top_index, $Allele2Top_index) = (0);	
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

my $n_field = scalar(@field); ## haoxiang
for my $i (0 .. @field-1) {
	$field[$i] eq 'SNP Name' and $name_index = $i;
	$field[$i] eq 'Sample ID' and $sample_index = $i;
	$field[$i] eq 'B Allele Freq' and $BAF_index = $i;
	$field[$i] eq 'Log R Ratio' and $LRR_index = $i;
	$field[$i] eq 'Chr' and $chr_index = $i;
	$field[$i] eq 'Position' and $position_index = $i;
	
	$field[$i] eq 'X' and $X_index = $i;
	$field[$i] eq 'Y' and $Y_index = $i;
	
	$field[$i] eq "Allele1 - Top" and $Allele1Top_index = $i;
	$field[$i] eq "Allele2 - Top" and $Allele2Top_index = $i;
	$field[$i] eq "Allele1 - Forward" and $Allele1Forward_index = $i;
	$field[$i] eq "Allele2 - Forward" and $Allele2Forward_index = $i;
}

## check whether all columns are included ## haoxiang
defined $name_index or confess "Error: the 'SNP Name' field is not found in header line of report file $reportfile: <$_>\n";
defined $sample_index or confess "Error: the 'Sample ID' field is not found in header line of report file $reportfile: <$_>\n";
defined $BAF_index or confess "Error: the 'B Allele Freq' field is not found in header line of report file $reportfile: <$_>\n";
defined $LRR_index or confess "Error: the 'Log R Ratio' field is not found in header line of report file $reportfile: <$_>\n";
defined $position_index or confess "Error: the 'Position' field is not found in header line of report file $reportfile: <$_>\n";
defined $chr_index or confess "Error: the 'Chr' field is not found in header line of report file $reportfile: <$_>\n";

# add warning 2018-11-26
if ( !defined $X_index ) {
    my $warn="Warning: the 'X' field is not found which is required by iPattern\n"; 
    print $warn; 
    push @warnings, $warn;
}
if ( !defined $Y_index ) {
    my $warn="Warning: the 'Y' field is not found which is required by iPattern\n"; 
    print $warn; 
    push @warnings, $warn;
}

my $flag_AllelesTop = defined $Allele1Top_index && defined $Allele2Top_index;
my $flag_AllelesForward = defined $Allele1Forward_index && defined $Allele2Forward_index;

my $flags_Alleles = $flag_AllelesTop || $flag_AllelesForward;

if ( !$flags_Alleles ) {
	my $warn="Warning: the paired ('Allele1 - Top' and 'Allele2 - Top') or ('Allele1 - Forward' and 'Allele2 - Forward') fields are not found which are required by iPattern\n";
    print $warn;
    push @warnings, $warn;
}

print "======== start checking integrity of finalreport: $reportfile ========\n";

print "number of fields in header line: $n_field\n";
## check if all rows have the same number of fields ####################################
my %sample_each_chr_nSNP = ();
my $content_line = 0;
while( my $line = <REPORT>) {
    
    $line =~ s/[\r\n]+$//;
    $content_line++;
    
    my @line = split(/\t/, $line);
    
    my $n1 = scalar(@line); ## each row must have the same number of fields as the header line
    confess "Error: different number of fields: $n1 than header line at line: $count_line\n" if $n_field != $n1 ;
    
    my $sample = $line[$sample_index];
    my $chr = $line[$chr_index];
    
    if ( !exists($sample_each_chr_nSNP{$sample}) ) {
        $sample_each_chr_nSNP{$sample}{$chr} = 1;
    
    } elsif ( exists($sample_each_chr_nSNP{$sample})) {
        if (exists($sample_each_chr_nSNP{$sample}{$chr})) {
            $sample_each_chr_nSNP{$sample}{$chr} = $sample_each_chr_nSNP{$sample}{$chr} + 1;
        } else {
            $sample_each_chr_nSNP{$sample}{$chr} = 1;
        }
    }
}

my $flag = 0;
my %chr_nSNP = ();
my $nSNP = 0;
my $nSample = 0;
foreach my $s (sort keys %sample_each_chr_nSNP) { ## sample
    my $nSNP_each_sample = 0;
    foreach my $c (sort keys %{$sample_each_chr_nSNP{$s}}) { ## chr
        
        if ( $flag == 0) {
            $chr_nSNP{$c} = $sample_each_chr_nSNP{$s}{$c};
            $nSNP = $nSNP + $sample_each_chr_nSNP{$s}{$c};
        }
        
        $nSNP_each_sample = $nSNP_each_sample + $sample_each_chr_nSNP{$s}{$c};
        
        ## each chr in each sample must have the same number of probes
        confess "Error: Sample: $s has a different number of probes: $sample_each_chr_nSNP{$s}{$c} on chr: $c than other samples.\nPlease check the integrity of finalreport $reportfile. You may need to regenerate the finalreport by Genome Studio.\n" if $sample_each_chr_nSNP{$s}{$c} != $chr_nSNP{$c};
    }

    ## each sample must have the same total number of probes
    confess "Error: Sample: $s has a diffrent total number of probes: $nSNP_each_sample than other samples.\nPlease check the integrity of finalreport $reportfile. You may need to regenerate the finalreport by Genome Studio.\n" if $nSNP_each_sample != $nSNP;

    $flag = 1;
    $nSample = $nSample + 1;
}

close(REPORT);

print "total number of probes: $nSNP\n";
print "total number of records: $content_line\n";
print "number of samples: $nSample\n";
if ( $nSample < 100 ) {
    my $warn = "Warning: Too few samples (<100) will lead to invalid results!!!\n";
    print $warn;
    push @warnings, $warn;
}

# summary information on each chr
my @array_chrs = (); # init following hash tables
print "summary information for each chr:\n";
for my $cc (sort keys %chr_nSNP) {
    print "chr: $cc\t\tnumber of probes: $chr_nSNP{$cc}\n";
    push @array_chrs, $cc;
}

print "======== end checking integrity of finalreport: $reportfile ========\n\n";
## end check ##################################################


print "======== start generating LRR and BAF matrices from finalreport: $reportfile ========\n";
## mkdir LRR BAF folder
my $path_output_LRR = join('/', $path_output, "LRR");
my $path_output_BAF = join('/', $path_output, "BAF");
`mkdir -p $path_output_LRR`;
`mkdir -p $path_output_BAF`;

## output file name
my $path_LRR = join("/", $path_output, "LRR/");
my $path_BAF = join("/", $path_output, "BAF/");
my @LRR_ends = ($path_LRR, ".tab");
my @BAF_ends = ($path_BAF, ".tab");
my $file_samples_order = $path_output."/"."samples_order.txt";
my $file_snps_number = $path_output."/"."snps_number.txt";
my $file_snps_name = $path_output."/"."snps_name.txt";
my $file_snps_position = $path_output."/"."SNP_pos.txt";
my %samples = ();
my $sample_order = 1;
my $sample_before = ();

open(REPORT, $reportfile) or confess "Error: cannot read from input report file $reportfile: $!\n";

while (<REPORT>) {
	m/^\[Data\]/ and last; 
}
my $line = <REPORT>;

# init
my %handle = ();
my %handle_LRR = ();
my %handle_BAF = ();
my %snp_chr_number = ();
my %snp_chr_name = ();

foreach my $chr (@array_chrs) {
    $handle{$chr} = ();
    $handle_LRR{$chr} = ();
    $handle_BAF{$chr} = ();
    $snp_chr_number{$chr} = ();
    $snp_chr_name{$chr} = ();
}


my %handle_save_LRR = ();
my %handle_save_BAF = ();
my %snp_chr_position = ();

my $flag_snp = "init";
my $flag_sampleID = "old";
my $flag_snp_save = "yes";

while ($line = <REPORT>) {
    
    $line =~ s/[\r\n]+$//;
    
    if ( eof ) {
        
        my @line = split(/\t/, $line);
	my $sample1 = $line[$sample_index];
	my $chr1 = $line[$chr_index];
	my $position1 = $line[$position_index];
	my $snp1 = $line[$name_index];
	my $LRR1 = $line[$LRR_index];
	my $BAF1 = $line[$BAF_index];
        
        $handle_LRR{$chr1}{$snp1} = $LRR1;
	$handle_BAF{$chr1}{$snp1} = $BAF1;
        
	foreach my $item (sort keys %handle) {
	    if (! exists $handle_save_LRR{$item}) {

		open my $fh_LRR, q{>}, join($item, @LRR_ends);
		open my $fh_BAF, q{>}, join($item, @BAF_ends);

		$handle_save_LRR{$item} = $fh_LRR;
		$handle_save_BAF{$item} = $fh_BAF;
	    }

	    my @array_chr_snp = ();
            my @array_chr_LRR = ();
            my @array_chr_BAF = ();
            
            foreach my $snpname (sort keys %{$handle_LRR{$item}}) {
    		push @array_chr_LRR, $handle_LRR{$item}{$snpname};
		push @array_chr_BAF, $handle_BAF{$item}{$snpname};
		push @array_chr_snp, $snpname;
	    }

	    my $combine_chr_snp = join("___", @array_chr_snp);

            ## 2018-11-05 haoxiang ##
	    my $n_sample = keys %samples;
	    if ( $n_sample > 1) {
		die "chr $item snp number are not same $snp_chr_number{$item}" unless (scalar(@array_chr_snp) == $snp_chr_number{$item});
		die "chr $item snp name are not same" unless ($combine_chr_snp eq $snp_chr_name{$item}); 
	    }			
			
	    if ( $n_sample == 1) {
                
                ## last snp
                if ( $flag_snp_save eq "yes" ) {
                    $snp_chr_position{$snp1}{$chr1} = $position1;
                } 

		$combine_chr_snp = join("___", @array_chr_snp);

		$snp_chr_number{$item} = scalar(@array_chr_snp);
		$snp_chr_name{$item} = $combine_chr_snp;

		open(my $fh_snp_number, ">$file_snps_number");
		foreach my $item_chr (keys %snp_chr_number) {
		    print $fh_snp_number "$item_chr\t$snp_chr_number{$item_chr}\n";
		}
		close $fh_snp_number;

		open(my $fh_snp_name, ">$file_snps_name");
		foreach my $item_chr (keys %snp_chr_name) {
		    print $fh_snp_name "$item_chr\t$snp_chr_name{$item_chr}\n";
		}
		close $fh_snp_name;

		open(my $fh_snp_position, ">$file_snps_position") or die "$!";
                print $fh_snp_position "Name\tChr\tPosition\n";
                foreach my $item_snp (keys %snp_chr_position) {
                    foreach my $item_chr (keys %{$snp_chr_position{$item_snp}}) {
                        print $fh_snp_position "$item_snp\t$item_chr\t$snp_chr_position{$item_snp}{$item_chr}\n";
                    }
                }
                close($fh_snp_position);
                
	    }
	    ## end ##

            my $chr_LRR = join("\t", $sample_before, join("\t", @array_chr_LRR));
	    my $chr_BAF = join("\t", $sample_before, join("\t", @array_chr_BAF));

            print { $handle_save_LRR{$item} } qq($chr_LRR\n);
	    print { $handle_save_BAF{$item} } qq($chr_BAF\n);
	}

	last;
    }


    my @line = split(/\t/, $line);
    
    my $sample1 = $line[$sample_index];
    my $chr1 = $line[$chr_index];
    my $position1 = $line[$position_index];
    my $snp1 = $line[$name_index];
    my $LRR1 = $line[$LRR_index];
    my $BAF1 = $line[$BAF_index];
    
    if ( (!exists($samples{$sample1})) && ($flag_sampleID eq "old") ) {
        $samples{$sample1} = $sample_order;
        
        print "order: $sample_order\tSample_ID: $sample1\n";
        $sample_order = $sample_order + 1;
        
        $flag_sampleID = "new";
        
        $handle_LRR{$chr1}{$snp1} = $LRR1;
        $handle_BAF{$chr1}{$snp1} = $BAF1;
            
        if ( $flag_snp_save eq "yes" ) {
            $snp_chr_position{$snp1}{$chr1} = $position1;
        }
    }
    
    if ( (exists($samples{$sample1})) && ($flag_sampleID eq "new") ) {
        
        $handle_LRR{$chr1}{$snp1} = $LRR1;
	$handle_BAF{$chr1}{$snp1} = $BAF1;
	$sample_before = $sample1;
            
        if ( $flag_snp_save eq "yes" ) {
            $snp_chr_position{$snp1}{$chr1} = $position1;
        } 
    }
    
    if ( (!exists($samples{$sample1})) && ($flag_sampleID eq "new") ) {
        
    ## starting with a new sample, so save data from the previous sample
	foreach my $item (sort keys %handle_LRR) {
	    if (! exists $handle_save_LRR{$item}) {

        #print join($item, @LRR_ends)."\n";
		#print join($item, @BAF_ends)."\n";

		open my $fh_LRR, q{>}, join($item, @LRR_ends);
		open my $fh_BAF, q{>}, join($item, @BAF_ends);

		$handle_save_LRR{$item} = $fh_LRR;
		$handle_save_BAF{$item} = $fh_BAF;
	    }

	    my @array_chr_LRR = ();
	    my @array_chr_BAF = ();
	    my @array_chr_snp = ();
	    foreach my $snpname (sort keys %{$handle_LRR{$item}}) {
		push @array_chr_LRR, $handle_LRR{$item}{$snpname};
		push @array_chr_BAF, $handle_BAF{$item}{$snpname};
		push @array_chr_snp, $snpname;			
	    }

	    my $combine_chr_snp = join("___", @array_chr_snp);

            ## only apply for the first sample
	    if ($flag_snp eq  "init") {
    		$snp_chr_number{$item} = scalar(@array_chr_snp);
                $snp_chr_name{$item} = $combine_chr_snp;
            } 
			
	    #print "$item\t".scalar(@array_chr_snp)."\t".scalar(keys %handle)."\n";
	    die "chr $item snp number are not same $snp_chr_number{$item}" unless (scalar(@array_chr_snp) == $snp_chr_number{$item});
	    die "chr $item snp name are not same" unless ($combine_chr_snp eq $snp_chr_name{$item}); 

	    my $chr_LRR = join("\t", $sample_before, join("\t", @array_chr_LRR));
	    my $chr_BAF = join("\t", $sample_before, join("\t", @array_chr_BAF));

	    print { $handle_save_LRR{$item} } qq($chr_LRR\n);
	    print { $handle_save_BAF{$item} } qq($chr_BAF\n);
	}

    ## only apply for the first sample
	if ($flag_snp_save eq "yes") {

	    open(my $fh_snp_number, ">$file_snps_number") or die "$!";
	    foreach my $item_chr (keys %snp_chr_number) {
                print $fh_snp_number "$item_chr\t$snp_chr_number{$item_chr}\n";
	    }
	    close $fh_snp_number;

	    open(my $fh_snp_name, ">$file_snps_name") or die "$!";
	    foreach my $item_chr (keys %snp_chr_name) {
		print $fh_snp_name "$item_chr\t$snp_chr_name{$item_chr}\n";
	    }
	    close $fh_snp_name;

            open(my $fh_snp_position, ">$file_snps_position") or die "$!";
            print $fh_snp_position "Name\tChr\tPosition\n";
            foreach my $item_snp (keys %snp_chr_position) {
                foreach my $item_chr (keys %{$snp_chr_position{$item_snp}}) {
                    print $fh_snp_position "$item_snp\t$item_chr\t$snp_chr_position{$item_snp}{$item_chr}\n";
                }
            }
            close($fh_snp_position);

	    $flag_snp_save = "no";
	}
		
	$flag_snp = "pass";  ## use the first sample information
	%handle_LRR = ();
	%handle_BAF = ();
	
	%handle_LRR = ();
	%handle_BAF = ();

    $samples{$sample1} = $sample_order;
    print "order: $sample_order\tSample_ID: $sample1\n";
	$sample_order  = $sample_order + 1;

	if ( !exists $handle_LRR{$chr1} ) {
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

close(REPORT);

print "======== end generating LRR and BAF matrices from finalreport: $reportfile ========\n";

print "\n";
print "=============== start checking consistency in output files ====================\n";

# number of samples
print "check $file_samples_order\n";
my $size = keys %samples;
confess "Error: inconsistent sample size in $file_samples_order\n" if $size != $nSample;

# number of SNP in each chr @array_chrs %chr_nSNP
print "check $file_snps_number file\n";
foreach my $chr (@array_chrs) {
    confess "Error: inconsistent number of probes in chr:$chr\n" if $chr_nSNP{$chr} != $snp_chr_number{$chr};  
}

print "check $file_snps_name\n";
foreach my $chr (@array_chrs) {
    my $n1_chr = $chr_nSNP{$chr};
    my $name = $snp_chr_name{$chr};
    my @name = split(/___/, $name);
    my $n2_chr = scalar(@name);
    confess "Error: inconsistent number of probes in chr:$chr\n" if $n1_chr != $n2_chr;
}

print "check $file_snps_position file \n";
my %nSNP_pos = ();
foreach my $chr (@array_chrs) {
    $nSNP_pos{$chr} = 0;
}

foreach my $item_snp (keys %snp_chr_position) {
    foreach my $item_chr (keys %{$snp_chr_position{$item_snp}}) {
        $nSNP_pos{$item_chr}++;
    }
}

foreach my $chr (@array_chrs) {
    confess "Error: inconsistent number of probes in chr:$chr\n" if $chr_nSNP{$chr} != $nSNP_pos{$chr};
}

print "check probe names in $file_snps_name and $file_snps_position\n";

my %snps_name = ();
foreach my $chr (keys %snp_chr_name) {
    my $name = $snp_chr_name{$chr};
    my @name = split(/___/, $name);
    foreach my $snp (@name) {
        $snps_name{$chr}{$snp}++;
    }
}

my %snps_pos = ();
# both side checked
foreach my $item_snp (keys %snp_chr_position) {
    foreach my $item_chr (keys %{$snp_chr_position{$item_snp}}) {
        $snps_pos{$item_chr}{$item_snp}++;
        confess "Error: probe: $item_snp in chr:$item_chr does not exist in $file_snps_position\n" if !exists($snps_name{$item_chr}{$item_snp});
    }
}

foreach my $chr (keys my %snp_name) {
    foreach my $snp (keys %{$snp_name{$chr}}) {
        confess "Error: probe: $snp in chr:$chr does not exist in $file_snps_position\n" if !exists($snps_pos{$chr}{$snp});
    }
}

print "================ end checking consistency in output files ===================\n\n\n";

my $nwarn = scalar(@warnings);
print "The processing of finalreport $reportfile has been successfully finished with $nwarn warnings.\n\n";
if ( @warnings ) {
    for my $i (0 .. @warnings-1) {
        print $warnings[$i];
    }
}

print "The output files are at $path_output\n";
