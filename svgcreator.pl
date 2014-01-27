use Getopt::Std;
use warnings;
use strict;
my $RefSeqFile = '/home/NKI/b.evers/mm10/RefSeq/refGene.txt';
my $HeaderFile = '/home/NKI/b.evers/mm10/scripts/header';
my $FooterFile = '/home/NKI/b.evers/mm10/scripts/footer';
my $ScriptName="svgcreator.pl";

#This script creates an svg graphical representation of the coding sequence of a gene with CRISPR location mapped onto it
#As input arguments, it takes the refseq id (r) and a SAM file containing the CRISPR locations to be mapped (i).
#As output argument, it takes a name of the svg file to create
my $YOffset=-20;
my $AnimationDistance = 30;
my $TriangleHeight = 15;
my $TriangleWidth = 4;
my $CollisionWidth=4;
my $LevelOffset = 15;

my %opts;
print "Usage:perl svgcreator.pl -i InputFile -o OutputFile -r RefSeqID\n";
getopt( 'roi', \%opts );
die "ERROR in $ScriptName: No RefSeq ID given.\n" unless my $RefSeq = $opts{'r'};
die "ERROR in $ScriptName: No Outputfile given.\n" unless my $OutputFile = $opts{'o'};
die "ERROR in $ScriptName: No Inputfile given.\n" unless my $InputFile = $opts{'i'};
open (IN, $InputFile) or die "Cannot open inputfile $InputFile\n";
my $RefSeqInfo = `grep -P "$RefSeq\t" $RefSeqFile`;
die "ERROR in $ScriptName: RefSeq $RefSeq cannot be found in the database.\n" if !$RefSeqInfo;
my @RefSeqValues = split( /\t/, $RefSeqInfo );

#First read in all intron/exon information of the gene we're plotting
my $Chromosome = $RefSeqValues[2];
#The refGene.txt file start sites are always one nt 5' of the actual nt, so compensate for that
my $GeneStart    = $RefSeqValues[4];
my $GeneEnd      = $RefSeqValues[5];
my $ProteinStart = $RefSeqValues[6];
my $ProteinEnd   = $RefSeqValues[7];
my $GeneOrientation = $RefSeqValues[3];
my $NumberOfExons  = $RefSeqValues[8];
my @ExonStartSites = split( /,/, $RefSeqValues[9] );
my @ExonEndSites   = split( /,/, $RefSeqValues[10] );

#Write the svg file
#First, write the header
`cat $HeaderFile > $OutputFile`;
#Determine the size of the mRNA
my $mRNASize = 0;
for (my $i=0;$i<$NumberOfExons;$i++) {
	$mRNASize =  $mRNASize + $ExonEndSites[$i] - $ExonStartSites[$i];
}

#Read in all CRISPRS to be mapped
my %DisplayObjects;
while (defined(my $Line = <IN>)) {
	chomp($Line);
	my @TargetSites = split( /\t/, $Line );
	my $TargetChromosome = $TargetSites[0];
	my $TargetOrientation = $TargetSites[5];
	my $TargetCutSite = $TargetSites[1]+23;
	$TargetCutSite = $TargetSites[1] if ($TargetOrientation eq '-');
	my $TargetLabel = $TargetSites[3];
	my $TargetScore = $TargetSites[11];
	my $TargetSequence=$TargetSites[9];

	#Verify that the target is in the gene, or within 250nt of the TSS
	if ($TargetChromosome eq $Chromosome && ($TargetCutSite >= $GeneStart-250 && $TargetCutSite <= $GeneEnd + 250)) {
		$DisplayObjects{$TargetOrientation}->{$TargetCutSite}->[0]=$TargetLabel . " - "  . $TargetSequence;
		$DisplayObjects{$TargetOrientation}->{$TargetCutSite}->[1]=$TargetScore;
	}
}

#Determine coding sequence position of cutsites
foreach my $Orientation (keys %DisplayObjects) {
	foreach my $TargetCutSite (keys $DisplayObjects{$Orientation}) {
		my $CodingSequencePosition = 0;
		if($GeneOrientation eq '+') {
			for (my $i=0;$i<$NumberOfExons;$i++) {
				if ($TargetCutSite <= $ExonEndSites[$i]) {
					$CodingSequencePosition = $CodingSequencePosition + ($TargetCutSite-$ExonStartSites[$i]);
					last;
				} 
				else
				{
					$CodingSequencePosition = $CodingSequencePosition + ($ExonEndSites[$i] - $ExonStartSites[$i]);
				}
			}
		}
		else {
			for (my $i=$NumberOfExons-1;$i>=0;$i--) {
				if ($TargetCutSite >= $ExonStartSites[$i]) {
					$CodingSequencePosition = $CodingSequencePosition + ($ExonEndSites[$i]-$TargetCutSite);
					last;
				} 
				else
				{
					$CodingSequencePosition = $CodingSequencePosition + ($ExonEndSites[$i] - $ExonStartSites[$i]);
				}
			}
		}
		my $RelativeMarkerPosition=1400*$CodingSequencePosition/$mRNASize;
		$RelativeMarkerPosition=$RelativeMarkerPosition+100;
		$DisplayObjects{$Orientation}->{$TargetCutSite}->[2]=$RelativeMarkerPosition;
	}
}

#Determine and assign colors
foreach my $Orientation (keys %DisplayObjects) {
	foreach my $TargetCutSite (keys $DisplayObjects{$Orientation}) {
		my $TriangleValue=$DisplayObjects{$Orientation}->{$TargetCutSite}->[1];
		$TriangleValue = int (255*($TriangleValue));
		my $TriangleColor = (255-$TriangleValue) . ",$TriangleValue,0"; 
		$DisplayObjects{$Orientation}->{$TargetCutSite}->[3]=$TriangleColor;
	}
}
		
#Do collision detection and resolution
my $MaxCollisionLevel=1;
foreach my $Orientation (keys %DisplayObjects) {
	my $CollisionsDetected=1;
	my $CollisionLevel=1;
	my %AllPositions;
	while ($CollisionsDetected) {
		$CollisionsDetected = 0;
		foreach my $TargetCutSite (keys $DisplayObjects{$Orientation}) {
			my $CollisionDetectedForThisTarget = 0;
			my $RelativeMarkerPosition=$DisplayObjects{$Orientation}->{$TargetCutSite}->[2];
			if(!($DisplayObjects{$Orientation}->{$TargetCutSite}->[4])) {
				if($AllPositions{$CollisionLevel}) {
					foreach my $Position (keys $AllPositions{$CollisionLevel}) {
						if ($RelativeMarkerPosition <= ($Position+$CollisionWidth) && $RelativeMarkerPosition >= ($Position-$CollisionWidth)) {
						#	print "Detected a collision with $RelativeMarkerPosition at position $Position at level $CollisionLevel\n";
							$CollisionsDetected++;
							$CollisionDetectedForThisTarget++;
							last;
						}
					}
				}
				if (!$CollisionDetectedForThisTarget) {
					$AllPositions{$CollisionLevel}->{$RelativeMarkerPosition}++;
					$DisplayObjects{$Orientation}->{$TargetCutSite}->[4]=$CollisionLevel;
				}
			}
		}
		$CollisionLevel++;
	}
	if($CollisionLevel > $MaxCollisionLevel && $Orientation eq '+') {
		$MaxCollisionLevel=$CollisionLevel;
	}
}
$YOffset=$YOffset+$AnimationDistance+($MaxCollisionLevel-2)*$LevelOffset;
print $RefSeq . ":" . $MaxCollisionLevel . "\t has offset $YOffset\n";

#Now print rectangles for all exons
open (OUT, ">>", $OutputFile) or die "Cannot open outpufile $OutputFile\n";
my $CurrentExon=0;
#Start 100 away from the beginning to accomodate for CRISPRs that are 5' of the TSS
my $CurrentPosition=0;
my $BoxSize;
$CurrentPosition=1500 if ($GeneOrientation eq '-');
#Start by drawing a piece of DNA upstream from the TSS
$BoxSize =  100;
print OUT "  <rect x=\"" . 0 . "\" y=\"" . (41+$YOffset) . "\" width=\"" . $BoxSize . "\" height=\"8\" fill=\"black\" stroke=\"black\" stroke-width=\"1\" />\'\n";
$CurrentPosition = $CurrentPosition + $BoxSize if ($GeneOrientation eq '+');

#Until you meet the proteinstart, plot full exons as UTR
while(!($ProteinStart >= $ExonStartSites[$CurrentExon] && $ProteinStart <= $ExonEndSites[$CurrentExon])) {
	$BoxSize =  1400*($ExonEndSites[$CurrentExon] - $ExonStartSites[$CurrentExon])/$mRNASize;
	$CurrentPosition=$CurrentPosition-$BoxSize if ($GeneOrientation eq '-');
	print OUT "  <rect x=\"" . $CurrentPosition . "\" y=\"" . (37+$YOffset) . "\" width=\"" . $BoxSize . "\" height=\"16\" fill=\"black\" stroke=\"black\" stroke-width=\"1\" onmousemove=\"ShowTooltip(evt, \'Exon " . ($GeneOrientation eq '+' ? ($CurrentExon+1) : $NumberOfExons-$CurrentExon) . "\')\" onmouseout=\"HideTooltip(evt)\"/>\'\n";	
	$CurrentPosition = $CurrentPosition + $BoxSize if ($GeneOrientation eq '+');
	$CurrentExon++;
}
#Print the piece of exon until the proteinstart as UTR
$BoxSize =  1400*($ProteinStart - $ExonStartSites[$CurrentExon])/$mRNASize;
$CurrentPosition=$CurrentPosition-$BoxSize if ($GeneOrientation eq '-');
print OUT "  <rect x=\"" . $CurrentPosition . "\" y=\"" . (37+$YOffset) . "\" width=\"" . $BoxSize . "\" height=\"16\" fill=\"black\" stroke=\"black\" stroke-width=\"1\" onmousemove=\"ShowTooltip(evt, \'Exon " . ($GeneOrientation eq '+' ? ($CurrentExon+1) : $NumberOfExons-$CurrentExon) . "\')\" onmouseout=\"HideTooltip(evt)\"/>\'\n";	
$CurrentPosition = $CurrentPosition + $BoxSize if ($GeneOrientation eq '+');
my $UseAlternativeStart=1;
#Until you meet the proteinend, plot exons as exon
while(!($ProteinEnd >= $ExonStartSites[$CurrentExon] && $ProteinEnd <= $ExonEndSites[$CurrentExon])) {
	$BoxSize =  1400*($ExonEndSites[$CurrentExon] - $ExonStartSites[$CurrentExon])/$mRNASize;
	$BoxSize =  1400*($ExonEndSites[$CurrentExon] - $ProteinStart)/$mRNASize if $UseAlternativeStart;
	$CurrentPosition=$CurrentPosition-$BoxSize if ($GeneOrientation eq '-');
	print OUT "  <rect x=\"" . $CurrentPosition . "\" y=\"" . (30+$YOffset) . "\" width=\"" . $BoxSize . "\" height=\"30\" fill=\"url(#grad1)\" stroke=\"black\" stroke-width=\"1\" onmousemove=\"ShowTooltip(evt, \'Exon " . ($GeneOrientation eq '+' ? ($CurrentExon+1) : $NumberOfExons-$CurrentExon) . "\')\" onmouseout=\"HideTooltip(evt)\"/>\'\n";	
	$CurrentPosition = $CurrentPosition + $BoxSize if ($GeneOrientation eq '+');
	$CurrentExon++;
	$UseAlternativeStart=0;
}
#Print the piece of exon until the proteinend as exon
$BoxSize =  1400*($ProteinEnd - $ExonStartSites[$CurrentExon])/$mRNASize;
$BoxSize =  1400*($ProteinEnd - $ProteinStart)/$mRNASize if $UseAlternativeStart;
$CurrentPosition=$CurrentPosition-$BoxSize if ($GeneOrientation eq '-');
print OUT "  <rect x=\"" . $CurrentPosition . "\" y=\"" . (30+$YOffset) . "\" width=\"" . $BoxSize . "\" height=\"30\" fill=\"url(#grad1)\" stroke=\"black\" stroke-width=\"1\" onmousemove=\"ShowTooltip(evt, \'Exon " . ($GeneOrientation eq '+' ? ($CurrentExon+1) : $NumberOfExons-$CurrentExon) . "\')\" onmouseout=\"HideTooltip(evt)\"/>\'\n";	
$CurrentPosition = $CurrentPosition + $BoxSize if ($GeneOrientation eq '+');
$UseAlternativeStart=1;	
#Print the rest as UTR
while($CurrentExon<$NumberOfExons) {
	$BoxSize =  1400*($ExonEndSites[$CurrentExon] - $ExonStartSites[$CurrentExon])/$mRNASize;
	$BoxSize =  1400*($ExonEndSites[$CurrentExon] - $ProteinEnd)/$mRNASize if $UseAlternativeStart;
	$CurrentPosition=$CurrentPosition-$BoxSize if ($GeneOrientation eq '-');
	print OUT "  <rect x=\"" . $CurrentPosition . "\" y=\"" . (37+$YOffset) . "\" width=\"" . $BoxSize . "\" height=\"16\" fill=\"black\" stroke=\"black\" stroke-width=\"1\" onmousemove=\"ShowTooltip(evt, \'Exon " . ($GeneOrientation eq '+' ? ($CurrentExon+1) : $NumberOfExons-$CurrentExon) . "\')\" onmouseout=\"HideTooltip(evt)\"/>\'\n";	
	$CurrentPosition = $CurrentPosition + $BoxSize if ($GeneOrientation eq '+');
	$CurrentExon++;
	$UseAlternativeStart=0;
}

#Perform printing of cut sites to output file
foreach my $Orientation (keys %DisplayObjects) {
	foreach my $TargetCutSite (keys $DisplayObjects{$Orientation}) {		
		my $RelativeMarkerPosition=$DisplayObjects{$Orientation}->{$TargetCutSite}->[2];
		my $TriangleColor=$DisplayObjects{$Orientation}->{$TargetCutSite}->[3];
		my $LayerOffset = $LevelOffset*($DisplayObjects{$Orientation}->{$TargetCutSite}->[4] - 1);
		if($Orientation eq '+') {
			print OUT "<polygon points=\"" . ($RelativeMarkerPosition - $TriangleWidth) . "\," . (30-$TriangleHeight+$YOffset-$LayerOffset) . " " . $RelativeMarkerPosition . "\," . (30+$YOffset-$LayerOffset) . " " . ($RelativeMarkerPosition+$TriangleWidth) . "\," . (30-$TriangleHeight+$YOffset-$LayerOffset) . "\" style=\"fill:rgb(" . $TriangleColor . ");stroke:black;stroke-width:1\" onmousemove=\"ShowTooltip(evt, \'" . $DisplayObjects{$Orientation}->{$TargetCutSite}->[0] . "\')\" onmouseout=\"HideTooltip(evt)\"\>\n";
			my $RandomTime=0.5*rand();
			print OUT "<animateTransform attributeName=\"transform\" attributeType=\"XML\" type=\"translate\" from=\"0 -" . $AnimationDistance . "\" to=\"0 0\" dur=\"" . $RandomTime . "s\"/>\n";
			print OUT "</polygon>\n";
		}
		else {
			print OUT "<polygon points=\"" . ($RelativeMarkerPosition - $TriangleWidth) . "\," . (60+$TriangleHeight+$YOffset+$LayerOffset) . " " . $RelativeMarkerPosition . "\," . (60+$YOffset+$LayerOffset) . " " . ($RelativeMarkerPosition+$TriangleWidth) . "\," . (60+$TriangleHeight+$YOffset+$LayerOffset) . "\" style=\"fill:rgb(" . $TriangleColor . ");stroke:black;stroke-width:1\" onmousemove=\"ShowTooltip(evt, \'" . $DisplayObjects{$Orientation}->{$TargetCutSite}->[0] . "\')\" onmouseout=\"HideTooltip(evt)\">\n";
			my $RandomTime=0.5*rand();
			print OUT "<animateTransform attributeName=\"transform\" attributeType=\"XML\" type=\"translate\" from=\"0 " . $AnimationDistance . "\" to=\"0 0\" dur=\"" . $RandomTime . "s\"/>\n";
			print OUT "</polygon>\n";
		}	
	}
} 
close (OUT) or die "ERROR in $ScriptName: Cannot close outputfile $OutputFile\n";

#Lastly, print the footer
`cat $FooterFile >> $OutputFile`;

