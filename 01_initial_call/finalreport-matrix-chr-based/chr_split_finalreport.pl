

use Data::Dumper;
my $finalreport = $ARGV[0];
my $sample_index = $ARGV[1];
my $chr_index = $ARGV[2];
my $position_index = $ARGV[3];
my $snpName_index = $ARGV[4];
my $LRR_index = $ARGV[5];
my $path_save_chr = $ARGV[6];
my $path_save_summary = $ARGV[7];


my @file_ends = ($path_save_chr, ".tab");
my $file_sample = $path_save_summary."sample_order.txt";
my $file_snp_number = $path_save_summary."snp_number.txt";
my $file_snp_name = $path_save_summary."snp_name.txt";
my $file_snp_position = $path_save_summary."snp_position.txt";
my %samples = ();
my $sample_order = 1;
my $sample_before = ();

# for test
# my $finalreport = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/raw_data_batch/batch1/FA_batch1_recluster_FinalReport.txt";
# my $finalreport = "/sc/orga/projects/haok01a/chengh04/Code_test/perl/tmp.finalreport";
# my $sample_index = 1;
# my $chr_index = 18;
# my $position_index = 19;
# my $snpName_index = 0;
# my $LRR_index = 34;

## save
# my @file_ends = qw(/sc/orga/projects/haok01a/chengh04/Food_Allergy/matrix_chr_batch/batch1/LRR/res/chr .tab);
# my $file_sample = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/matrix_chr_batch/batch1/LRR/sample_order.txt";
# my $file_snp_number = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/matrix_chr_batch/batch1/LRR/file_snp_number.txt";
# my $file_snp_name = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/matrix_chr_batch/batch1/LRR/file_snp_name.txt";
# my $file_snp_position = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/matrix_chr_batch/batch1/LRR/file_snp_position.txt";
# my %samples = ();
# my $sample_order = 1;
# my $sample_before = ();
# my $BAF_index = 33;

my %handle = ();
for (1..22) {
	$handle{$_} = ();
}

my %handle_save = ();

my %snp_chr_number = ();
for (1..22) {
	$snp_chr_number{$_} = ();
}

my %snp_chr_name = ();
for (1..22) {
	$snp_chr_name{$_} = ();
}

my %snp_chr_position = ();  
$flag_snp = "init";

open(IN, "<$finalreport") || die "can not open file $finalreport: $!";

my $flag_header = 0;
my $flag_sampleID = "old";
my $flag_snp_save = "yes";

while($line = <IN>) {

	if ( eof ) {

		chomp $line;

		@line = split(/\t/, $line);
		$sample1 = $line[$sample_index];
		$chr1 = $line[$chr_index];
		$position1 = $line[$position_index];
		$snp1 = $line[$snpName_index];
		$LRR1 = $line[$LRR_index];

		if ( exists $handle{$chr1} ) {
			$handle{$chr1}{$snp1} = $LRR1;
		}

		# print "condition_eof\n";
		# print "$line\n";
		foreach my $item (sort keys %handle) {
			if (! exists $handle_save{$item}) {

				$out_file = join($item, @file_ends);
				# print "$item\t$out_file\n";

				open my $fh, q{>>}, join($item, @file_ends);
				$handle_save{$item} = $fh;
			}

			# print Dumper(\$handle{$item});
			# print ref($handle{$item});
			my @array_chr_LRR = ();
			my @array_chr_snp = ();
			foreach my $snpname (sort keys $handle{$item}) {
				push @array_chr_LRR, $handle{$item}{$snpname};
				push @array_chr_snp, $snpname;
			}

			$combine_chr_snp = join("___", @array_chr_snp);
						
			die "chr $item snp number are not same $snp_chr_number{$item}" unless (scalar(@array_chr_snp) == $snp_chr_number{$item});
			die "chr $item snp name are not same" unless ($combine_chr_snp eq $snp_chr_name{$item}); 

			my $chr_LRR = join("\t", $sample_before, join("\t", @array_chr_LRR));
			print { $handle_save{$item} } qq($chr_LRR\n);
		}

		last;
	}

	if ( $line =~ /^SNP Name/) {
		$flag_header = 1;
		next;
	}
	next unless ($flag_header);
	chomp $line;

	@line = split(/\t/, $line);
	$sample1 = $line[$sample_index];
	$chr1 = $line[$chr_index];
	$position1 = $line[$position_index];
	$snp1 = $line[$snpName_index];
	$LRR1 = $line[$LRR_index];

	# print "$sample1\t$chr1\t$position1\t$snp1\t$LRR1\n";

 	## for first SampleID special
	if ((! exists $samples{$sample1}) && ($flag_sampleID eq "old")) {
		$samples{$sample1} = $sample_order;
		$sample_order = $sample_order + 1;
		print "$sample1\n";
		# print "condition_1\n";

		$flag_sampleID = "new";
		if ( exists $handle{$chr1} ) {
			$handle{$chr1}{$snp1} = $LRR1;

			if ($flag_snp_save eq "yes") {
				$snp_chr_position{$snp1} = $position1;
			}
		}
	} 

	if ( (exists $samples{$sample1}) && ($flag_sampleID eq "new") ) {

		if ( exists $handle{$chr1} ) {
			$handle{$chr1}{$snp1} = $LRR1;
			$sample_before = $sample1;
			# print "$sample1\n";
			if ($flag_snp_save eq "yes") {
				$snp_chr_position{$snp1} = $position1;
			}
		}

	}

	if ( (! exists $samples{$sample1}) && ($flag_sampleID eq "new") ) {

		foreach my $item (sort keys %handle) {
			if (! exists $handle_save{$item}) {

				print join($item, @file_ends)."\n";

				open my $fh, q{>>}, join($item, @file_ends);
				$handle_save{$item} = $fh;
			}

			@array_chr_LRR = ();
			@array_chr_snp = ();
			foreach my $snpname (sort keys $handle{$item}) {
				push @array_chr_LRR, $handle{$item}{$snpname};
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
			print { $handle_save{$item} } qq($chr_LRR\n);
		}

		if ($flag_snp_save eq "yes") {

			open($fh_snp_number, ">$file_snp_number");
			foreach my $item_chr (keys %snp_chr_number) {
				print $fh_snp_number "$item_chr\t$snp_chr_number{$item_chr}\n";
			}
			close $fh_snp_number;

			open($fh_snp_name, ">$file_snp_name");
			foreach my $item_chr (keys %snp_chr_name) {
				print $fh_snp_name "$item_chr\t$snp_chr_name{$item_chr}\n";
			}
			close $fh_snp_name;

			open($fh_snp_position, ">$file_snp_position");
			foreach my $item_snp (keys %snp_chr_position) {
				print $fh_snp_position "$item_snp\t$snp_chr_position{$item_snp}\n";
			}
			close $fh_snp_position;

			$flag_snp_save = "no";
		}
		$flag_snp = "pass";  ## use the first sample information
		## init %handle to save all chr in one sample 
		# print "$sample1\n";
		my %handle = ();
		for (1..22) {
			$handle{$_} = ();
		}
		$samples{$sample1} = $sample_order;
		$sample_order  = $sample_order + 1;
		print "$sample1\n";
		if ( exists $handle{$chr1} ) {
			$handle{$chr1}{$snp1} = $LRR1;
		}
	}
}

close IN;
map { close $handle_save{$_} } keys %handle_save;

open(OUT, ">$file_sample");
foreach my $item (keys %samples) {
	print OUT "$item\t$samples{$item}\n";
}
close OUT;

