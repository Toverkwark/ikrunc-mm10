#!/bin/perl
$OutputFile="RunAll.sh";
$NumberOfParallelProcesses=64;
open (OUT, ">", $OutputFile) or die "ERROR in $ScriptName: Cannot open outpufile $OutputFile\n";
print OUT "#!/bin/bash\n";
for ($i=0;$i<$NumberOfParallelProcesses;$i++) {
	print OUT "screen -d -m nice -n 20 parallel/CRISPRSearchChunk.$i\n";
}

