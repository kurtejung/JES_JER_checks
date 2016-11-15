#!/usr/bin/perl

#sub produce_small_file_segment;

$filelist = $ARGV[0];
$njobs = $ARGV[1];
$nfiles = $ARGV[2];

$queue="cmscaf1nd";
$cmsswBase=$ENV{'CMSSW_BASE'};
$script="Validate_Jets.C";
$curWorkDir=$ENV{'PWD'};

open(FILE, "$filelist");
@lines = <FILE>;
close(FILE);

$itotal = @lines;

print "number of files: $itotal \n";

$ntot = 0;

for ($ifile = 0; $ifile < $njobs ; $ifile ++) {

    $istart = $ifile*$nfiles;
    $iend = $nfiles + $istart;

    if($iend >$itotal) {
        $iend = $itotal;
    }

    if($istart < $start_from) {
        next;
    }

    #-- produce condor configuration file ...
    $config_file = "candtree\_$istart\_$iend.sh";
    &produce_config_file;

    system("chmod u+x $config_file");
    print "subbing $config_file to $queue...\n";
    system("bsub -q $queue -N -u null@null -R \"swp>1000 && pool>30000\" -J validation_$ifile \'$config_file\' ");

    if($iend == $itotal) {
        last;
    }
}

#-----------------------------------------------------
# produce condor configuation file for each file list  

sub produce_config_file {
	open(OUTPUTFILE,">$config_file");
	print OUTPUTFILE ("cd $cmsswBase/src\n");
	print OUTPUTFILE ("eval `scram r -sh`\n");
	print OUTPUTFILE ("cd -\n");
	print OUTPUTFILE ("cp $curWorkDir/$script .\n");
	print OUTPUTFILE ("cp $curWorkDir/ForestContainer.h . \n");
	print OUTPUTFILE ("cp $curWorkDir/boundaries.h .\n");
	print OUTPUTFILE ("cp $curWorkDir/$filelist .\n");
	print OUTPUTFILE ("root -l -b -q $script\\+\\($istart,$iend,\\\"$filelist\\\"\\) \n");
	print OUTPUTFILE ("cp *.root $curWorkDir \n");
	print OUTPUTFILE ("exit 0\n");
}	
