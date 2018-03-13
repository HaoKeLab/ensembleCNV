

## read in selected snps
$in_dir_snps="/sc/orga/projects/haok01a/chengh04/MEGA/MegaEX_Inga/analysis_Part1/batch_effect/dat/";
$file_snps="snps.randomly.select.txt";

$in_snps=$in_dir_snps.$file_snps;

open(IN, "<$in_snps") or die "can't open snps file $in_snps: $!";
%snps=();
while ($line=<IN>) {
	chomp $line;
	print "$line\n";
	$snps{$line}++;
}

close IN;
@snps=(keys %snps);
print "total number of snps:".scalar(@snps)."\n";

## build matrix file for all samples
$in_dir="/sc/orga/projects/haok01a/chengh04/MEGA/MegaEX_Inga/FinalReport/Plates_Part1/";
$file_name="MegaEX_Inga_Part1_FinalReport.txt";

$in_file=$in_dir.$file_name;

%samples=(); ## for all samples
%hash=();
open(IN, "<$in_file") or die "can't open finalreport $in_file: $!";

$flag=0;
$LRR=0; ## [Data]
$idx=0;
$SampleID=0;
$SNPName=0;
$flagsample=0; ##  lrr save flag
$lrrsample=(); ##
$SampleIDraw=();
$total=0;
$flageof=0;

while ($line=<IN>) {

	$flageof=1 if eof; ## add file eof flag
	chomp $line; 
	$flag=1 if $line=~/^SNP Name/;
	next unless ($flag);

	if ($flag==1) {

		chomp $line;
		@line=split(/\t/, $line);
		print scalar(@line)."\n";
		# foreach my $item (@line) {
		# 	chomp $item;
		# 	print $item."\n";
		# 	$LRR=$idx if ($item=~/^Log R Ratio$/);
		# 	$SampleID=$idx if ($item=~/^Sample ID$/);
		# 	$SNPName=$idx if ($item=~/^SNP Name$/);
		# 	$idx++;
		# }
		$LRR=34;
		$SampleID=1;
		$SNPName=0;
		print "$line[$LRR]\n$LRR\n$line[$SampleID]\n$SampleID\n$line[$SNPName]\n$SNPName\n";
		$flag=2;
	} elsif ($flag==2) {

		@line=split(/\t/, $line);

		## tansform Log R Ratio

		if (exists($samples{$line[$SampleID]})&&exists($snps{$line[$SNPName]})) {

			$lrrvalue=$line[$LRR];
			$lrrvalue=~ tr/\015//d;
			# $lrrsample=$lrrsample."\t".$line[$LRR];
			$lrrsample=$lrrsample."\t".$lrrvalue;
			$flagsample=1;
			$SampleIDraw=$line[$SampleID];	
			$total++;

			if ($flageof == 1) {
				$hash{$SampleIDraw}=$lrrsample;
				print "SampleID:$SampleIDraw\t".scalar(keys %samples)."\t$total\n";
				last;
			}

		} elsif (exists($samples{$line[$SampleID]})) {
			
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
				if (exists($snps{$line[$SNPName]})) {
					$samples{$line[$SampleID]}++;
					$lrrvalue=$line[$LRR];
					$lrrvalue=~ tr/\015//d;
					$lrrsample=$lrrvalue;
					$total++;
				}
			} elsif ($flagsample==1) {
				if (exists($snps{$line[$SNPName]})) {

					$hash{$SampleIDraw}=$lrrsample;
					print "SampleID:$SampleIDraw\t".scalar(keys %samples)."\t$total\n";

					$samples{$line[$SampleID]}++;
					$lrrsample=();
					$lrrvalue=$line[$LRR];
					$lrrvalue=~ tr/\015//d;
					$lrrsample=$lrrvalue;
					$total=1;
				}
			}

		}

	}

}

close IN;

## save LRR matrix
$out_file="/sc/orga/projects/haok01a/chengh04/MEGA/MegaEX_Inga/analysis_Part1/batch_effect/res/matrix.LRR.snps.randomly.select.txt";

open(OUT, ">", $out_file) or die $!;

foreach my $item (keys %hash) {
	print OUT "$item\t$hash{$item}\n";
}

close OUT;
