use Getopt::Std;
use File::Basename;
use warnings;
use strict;
my $RefSeqFile = '/home/NKI/b.evers/mm10/RefSeq/refGene.txt';
my $ScriptName="FilterAndSelect.pl";
my $SelectNumberOfProtospacers;
my $CDSSearchFraction;
my $StartNeighbourhood;
my %opts;
my %Protospacers;
my %ProtospacerSequences;
getopt( 'iocsnr', \%opts );
print "FilterAndSelect.pl -iocsnr\n-i=InputFile\n-o=OutputFile\n-c=Fraction of the coding region to consider, default 1\n-s=#nucleotides around the start codon to consider, default 0\n-n Select the top # of protospacers, default 20, infinite when -1\n-r RefSeq ID (optional, if not given it's determined from input filename)\n";
die "ERROR in $ScriptName: No inputfile given.\n" unless my $InputFile = $opts{'i'};
die "ERROR in $ScriptName: No Outputfile given.\n" unless my $OutputFile = $opts{'o'};
$CDSSearchFraction = 1 unless $CDSSearchFraction = $opts{'c'};
$StartNeighbourhood = 0 unless $StartNeighbourhood = $opts{'s'};
die "ERROR in $ScriptName: Start neighbourhood can maximally be 150.\n" if ($StartNeighbourhood > 150); 
$SelectNumberOfProtospacers = 20 unless $SelectNumberOfProtospacers = $opts{'n'};

open (OUT, ">", $OutputFile) or die "ERROR in $ScriptName: Cannot open outputfile $OutputFile\n";
open (OUTSEQ, ">", $OutputFile . ".seq") or die "ERROR in $ScriptName: Cannot open outputfile with sequences " . $OutputFile . ".seq\n";

#Find and extract gene information
my $RefSeq;
$RefSeq = substr(basename($InputFile),0,index(basename($InputFile),'.')) unless $RefSeq = $opts{'r'};
my $RefSeqInfo = `grep -P '$RefSeq\t' $RefSeqFile`;
die "ERROR in $ScriptName: RefSeq $RefSeq cannot be found in the database.\n" if !$RefSeqInfo;
my @RefSeqValues = split( /\t/, $RefSeqInfo );
my $Chromosome = $RefSeqValues[2];
my $GeneOrientation = 0;
$GeneOrientation = 1 if($RefSeqValues[3] eq '+');
my $GeneStart    = $RefSeqValues[4];
my $GeneEnd      = $RefSeqValues[5];
my $ProteinStart = $RefSeqValues[6];
my $ProteinEnd   = $RefSeqValues[7];
my $NumberOfExons  = $RefSeqValues[8];
my @ExonStartSites = split( /,/, $RefSeqValues[9] );
my @ExonEndSites   = split( /,/, $RefSeqValues[10] );

#Write browser and track information to the output file
print OUT "browser position " .$Chromosome . ":" . $GeneStart . "-" . $GeneEnd . "\n";
print OUT "browser full refGene\n";
print OUT "track name='CRISPRs' description='Valid CRISPR protospacers' visibility='pack' itemRgb='On' db='hg19'\n";

#Determine the size of the protein coding region
my $CDSSize=0;
for (my $Exon = 0;$Exon < $NumberOfExons; $Exon++) {
	my $ExonStart= ($ExonStartSites[$Exon]);
	my $ExonEnd = ($ExonEndSites[$Exon]); 
	
	#Check if the current exon has coding sequence. If it has, set StartSite and EndSite according to coordinates that fall within protein coding sequence
	if($ProteinStart<$ExonEnd && $ProteinEnd>$ExonStart) {
		my $StartSite = ($ProteinStart > $ExonStart) ? $ProteinStart : $ExonStart;
		my $EndSite = ($ProteinEnd < $ExonEnd) ? $ProteinEnd : $ExonEnd;
		$CDSSize = $CDSSize + ($EndSite - $StartSite); 
	}
}

#Set the limit of how far to look into the CDS
my $CDSLimit = int(($CDSSize * $CDSSearchFraction) + 0.5);
$CDSLimit = $CDSSize - $CDSLimit if ($GeneOrientation == 0);

#Loop over the input file containing all target sites with off-target information
open (IN, $InputFile) or die "ERROR in $ScriptName: Cannot open inputfile $InputFile\n";
while (defined(my $Line = <IN>)) {
	chomp($Line);
	#For every line, extract the necessary information
	my @SAMValues = split( /\t/, $Line );
	my $ProtospacerSequence = $SAMValues[9];
	my $Orientation = '+';
	if($SAMValues[1] & 16) {
		$Orientation = '-';
		$ProtospacerSequence =~ tr/ACTG/TGAC/;
		$ProtospacerSequence = reverse($ProtospacerSequence);
	}
	$Chromosome = $SAMValues[2];
	my $Start = $SAMValues[3]-1;
	my $End = $Start+23;
	my $CutLocation;
	if ($Orientation eq '+') {
		$CutLocation = $End-5;
	}
	else {
		$CutLocation = $Start+6;
	} 

	my $NumberOfIdentical3PrimeTargets=$SAMValues[19];
	my $NumberOfIdentical3PrimeTargetsNearExons=$SAMValues[20];
	my $Degree=1;
	my $ClosestRelatives=$SAMValues[21];
	my $ClosestRelativesNearExons=$SAMValues[22];
	for (my $i=23;$i<=27;$i=$i+2) {
		unless ($ClosestRelatives > 0) {
			$ClosestRelatives=$SAMValues[$i];
			$ClosestRelativesNearExons=$SAMValues[$i+1];
			$Degree=$Degree+1;
		}
	}
	my $Name = $NumberOfIdentical3PrimeTargets . "(" . $NumberOfIdentical3PrimeTargetsNearExons . "):". $Degree . "MM:" . $ClosestRelatives . "(" . $ClosestRelativesNearExons . ")";
	
	#Determine the color and the score of any protospacer
	my $Score =0;
	my $ColorRed=0;
	my $ColorGreen=0;
	my $ColorBlue=1;
	
	#The score of any protospacer will be determined as follows
	#There will be one point added for every Degree, meaning the degree below that did not have relatives
	#3 Points will be added for protospacers that do not have OTHER identical 3' targets
	#2 Points will be added for protospacers that have 1 OTHER identical 3' target, which is not near Exons
	#1 Points will be added for protospacers that have more than 1 OTHER identical 3' target unless any of them near Exons
	#0.0001 points will be subtracted for the number of relatives in the closest degree
	#0.01 points will be substracted for the number of relatives in the closest degree that are near Exons
	$Score=$Degree;
	$Score=$Score + 3 if ($NumberOfIdentical3PrimeTargets == 1);
	$Score=$Score + 2 if ($NumberOfIdentical3PrimeTargets == 2 && $NumberOfIdentical3PrimeTargetsNearExons == 0);
	$Score=$Score + 1 if ($NumberOfIdentical3PrimeTargets > 2 && $NumberOfIdentical3PrimeTargetsNearExons == 0);
	$Score=$Score - (0.0001 * $ClosestRelatives);
	$Score=$Score - (0.01 * $ClosestRelativesNearExons);
	
	#Now, see if the protospacer is within the set criteria for CDS fraction or start codon neighbourhood
	my $WithinStartCodonNeighbourHood = 0;
	my $WithinChosenCDSFraction = 0;
	my $StartCodon=$ProteinStart;
	if(!$GeneOrientation) {
		$StartCodon=$ProteinEnd
	}
	if(abs($CutLocation - $StartCodon) <= $StartNeighbourhood) {
		$WithinStartCodonNeighbourHood = 1;
	}
	
	#Find what location in the CDS the cut location is
	my $CutLocationCDS=0;
	for (my $Exon = 0;$Exon < $NumberOfExons; $Exon++) {
		my $ExonStart= ($ExonStartSites[$Exon]);
		my $ExonEnd = ($ExonEndSites[$Exon]); 
	
		#Check if the current exon has coding sequence. If it has, set StartSite and EndSite according to coordinates that fall within protein coding sequence
		if($ProteinStart<$ExonEnd && $CutLocation>$ExonStart) {
			my $StartSite = ($ProteinStart > $ExonStart) ? $ProteinStart : $ExonStart;
			my $EndSite = ($CutLocation < $ExonEnd) ? $CutLocation : $ExonEnd;
			$CutLocationCDS = $CutLocationCDS + ($EndSite - $StartSite); 
		}
	}
	if($GeneOrientation == 1) {
		$WithinChosenCDSFraction = 1 if ($CutLocationCDS<$CDSLimit && $CutLocationCDS>0); 
	}
	else {
		$WithinChosenCDSFraction = 1 if ($CutLocationCDS>$CDSLimit && $CutLocationCDS<$CDSSize); 
	}
	#print "Protospacer info found for $RefSeq:$Chromosome:$Start\t$End\t$Orientation\t$CutLocation\t$WithinStartCodonNeighbourHood\t$WithinChosenCDSFraction\n";
	#print "$StartNeighbourhood\t$CutLocation\t$ProteinStart\t$WithinStartCodonNeighbourHood\t$WithinChosenCDSFraction\n";
	if ($WithinChosenCDSFraction || $WithinStartCodonNeighbourHood) {
		$Protospacers{"$Chromosome\t$Start\t$End\t$Name\t$Score\t$Orientation\t$Start\t$End\t$ColorRed,$ColorGreen,$ColorBlue\n"} = $Score;
		$ProtospacerSequences{"$Chromosome\t$Start\t$End\t$Name\t$Score\t$Orientation\t$Start\t$End\t$ColorRed,$ColorGreen,$ColorBlue\t$ProtospacerSequence\t$RefSeq"} = $Score;
	}
}

my $NumberOfProtospacers = 0;
foreach my $Protospacer (sort {$Protospacers{$b} <=> $Protospacers{$a}} keys %Protospacers) {
	print OUT $Protospacer;
	$NumberOfProtospacers++;
	if ($SelectNumberOfProtospacers>0) {
		last if ($NumberOfProtospacers >= $SelectNumberOfProtospacers);	
	}
}
my $MaxNumberOfProtospacers = $NumberOfProtospacers;
$NumberOfProtospacers = 0;
foreach my $ProtospacerSequence (sort {$ProtospacerSequences{$b} <=> $ProtospacerSequences{$a}} keys %ProtospacerSequences) {
	print OUTSEQ $ProtospacerSequence . "\t" . (($MaxNumberOfProtospacers - $NumberOfProtospacers) / $MaxNumberOfProtospacers) . "\n";
	$NumberOfProtospacers++;
	last if ($NumberOfProtospacers >= $MaxNumberOfProtospacers);
}

close (IN) or die "ERROR in $ScriptName: Cannot close inputfile $InputFile\n";
close (OUT) or die "ERROR in $ScriptName: Cannot close outputfile $OutputFile\n";
close (OUTSEQ) or die "ERROR in $ScriptName: Cannot close outputfile with sequences " . $OutputFile . ".seq\n";
