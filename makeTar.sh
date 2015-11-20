#!bin/bash

root -l <<EOF
.L Validate_Jets.C++
.q
EOF

rm run_MC.tar
tar -zcvf run_MC.tar *.txt Validate_Je*.* *.h
