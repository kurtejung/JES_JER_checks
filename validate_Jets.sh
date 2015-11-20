#!/bin/bash

echo "going to run for PF jets R=0.4"
root -b -l <<EOF
.x Validate_Jets.C+(0, 1, 4, "PP", "Data", "PF", "Pu", "PromptValidation")
.q
EOF

echo "going to run for Calo jets R=0.4"
root -b -l <<EOF
.x Validate_Jets.C+(0, 1, 4, "PP", "Data", "Calo", "Pu","PromptValidation")
.q
EOF
