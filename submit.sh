#!/bin/bash

cd /afs/cern.ch/work/r/rkunnawa/Run2_ECAL/JES_checks/CMSSW_7_5_5/src/JES_JER_checks/
#cmsenv
eval `scramv1 runtime -sh`


#Added by Ian
export X509_USER_PROXY=~/x509_user_proxy/proxy
voms-proxy-init --noregen
#</Ian>

cd /afs/cern.ch/work/r/rkunnawa/Run2_ECAL/JES_checks/CMSSW_7_5_5/src/JES_JER_checks/

echo "root -l -b -q runForest.C++"
echo "First = $FIRST and last file = $LAST"   

root -b > macro${FIRST}.log <<EOF
.x runForest_histosJESJER.C($FIRST,$LAST,"Pu","Calo","pthat80_globalfit_marta_${LAST}.root");
.q
EOF

echo "Done all jobs!"

#destination="/path/to/output/file/"

#echo "Copying output files to " $destination

#mv test_outputFile_$LAST.root $destination 
