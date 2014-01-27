use Getopt::Std;
my %opts;
getopt('c', \%opts );
my $chromosome = $opts{'c'};
my $GENOMELOCATION='/home/NKI/b.evers/mm10/genome-fasta';
{
	open (IN, "$GENOMELOCATION/chr$chromosome.fa") or die "Could not open file\n";
	open (OUT, ">", "$GENOMELOCATION/chr$chromosome.stripped.fa") or die "Could not open output file\n";;
	read IN,$temp,length('chr' . $chromosome)+1;
	while (1) {
		if(eof(IN)) {
			last;
		}
		$i++;
	#	print "$i\n";
#		print "$i\n"  if(!($i % 1000000));
		read IN, $temp, 1;
		if(ord($temp) != 10) {
			syswrite OUT, $temp;
		}
	}
	close(IN);
	close(OUT);
}
