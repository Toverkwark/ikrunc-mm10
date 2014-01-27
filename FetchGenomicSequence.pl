sub FetchGenomicSequence($$$) {
	my $GENOMELOCATION='/home/NKI/b.evers/mm10/genome-fasta';
	my($Chromosome, $StartSite, $EndSite) = @_;
	my $GenomicSequence;
	open (IN, "$GENOMELOCATION/chr$Chromosome.stripped.fa") or die "Could not open genomic sequence file $GENOMELOCATION/chr$Chromosome.stripped.fa\n";
	if($EndSite>=$StartSite)
	{
		seek(IN,$StartSite-1,1);
		read IN,$GenomicSequence,$EndSite-$StartSite+1;
	}
	close (IN);
	$GenomicSequence =~ tr/a-z/A-Z/;
	return $GenomicSequence;
}
1;
