
$in_dir="";
$out_dir="";
$out_file=$out_dir."quantisnp.cnv";

opendir(DIR, $in_dir) or "cannot open $in_dir:$!";
open(OUT1, ">", $out_file) or die $!;

$flag = 1;
while (defined($folder = readdir(DIR))) {
	
	next if ($folder=~/^\./);
	##next if ($folder=~/^INTERNAL/);
	##next if ($folder=~/^CONTROL/);
	$filename=$folder.".cnv";
	$file=$in_dir.$folder."/".$filename;
	print "$flag", "$file\n";
	open(IN1, "<$file") or die $!;
	while ($line=<IN1>) {
		next if ($line=~/^Sample/);
		print OUT1 $line;
	}

	close IN1;
	$flag = $flag + 1;
}

close OUT1;

