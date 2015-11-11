#!bin/bash

root -l <<EOF
.L runForest_histosJESJER.C++
.q
EOF

rm run_MC.tar
tar -zcvf run_MC.tar mergedfile.txt runForest*.*
