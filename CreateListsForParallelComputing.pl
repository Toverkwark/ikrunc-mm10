$NumberOfCores=64;
my $RefSeqFile = '/home/NKI/b.evers/mm10/RefSeq/refGene.txt';
my $PARALLELDIR = '/home/NKI/b.evers/mm10/scripts/parallel';
my $OUTPUTDIR = '/home/NKI/b.evers/mm10/output';
open (IN, $RefSeqFile) or die;
$CurrentCore=0;
for (my $i=0;$i<$NumberOfCores;$i++) {
	$OutputFileHandle = 'OUT.' . $i;
	$OutputFile = "$PARALLELDIR/CRISPRSearchChunk." . $i;
	#print "Opening $OutputFileHandle with file $OutputFile\n";
	open ($OutputFileHandle, ">", $OutputFile) or die "ERROR: Cannot open outputfile $OutputFile\n";
}

while (defined(my $Line = <IN>)) {
	chomp($Line);
	my @RefSeqValues = split( /\t/, $Line );
	$RefSeqID = $RefSeqValues[1];
	#unlink "/home/NKI/b.evers/output/$RefSeqID.qualities.4";
	unless (-e '$OUTPUTDIR/$RefSeqID.qualities.4') {
		$CurrentCore = 0 if ($CurrentCore==$NumberOfCores);
		$OutputFileHandle = 'OUT.' . $CurrentCore;
		select $OutputFileHandle;
		print "/home/NKI/b.evers/mm10/scripts/RunRefSeq.sh $RefSeqID\n";
		$CurrentCore++;
	}
}
close (IN) or die "ERROR: Cannot close inputfile";
for (my $i=0;$i<$NumberOfCores;$i++) {
	$OutputFileHandle = 'OUT.' . $i;
	$OutputFile = '$PARALLELDIR/CRISPRSearchChunk.' . $CurrentCore;
	close ($OutputFileHandle) or die "ERROR: Cannot close outputfile $OutputFile\n";
}
