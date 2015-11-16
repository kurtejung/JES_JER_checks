#!/bin/bash
cd /afs/cern.ch/work/r/rkunnawa/Run2_ECAL/JES_checks/CMSSW_7_5_5/src/JES_JER_checks/
#cmsenv
eval `scramv1 runtime -sh`

cd /afs/cern.ch/work/r/rkunnawa/Run2_ECAL/JES_checks/CMSSW_7_5_5/src/JES_JER_checks/

nJobs=97
i=0
while [ $i -le $nJobs ];
do 
   let "start=i*1"
   let "end=(i+1)*1"
  export FIRST=$start  
  export LAST=$end 
  echo "First = $FIRST and last file = $LAST"   
  bsub -R "pool>5000" -M 3000000 -q 1nh -J submit_job_${i} < /afs/cern.ch/work/r/rkunnawa/Run2_ECAL/JES_checks/CMSSW_7_5_5/src/JES_JER_checks/submit.sh
  let "i++"
done

echo "submit all jobs!"

#echo "Copying output files to " $destination
