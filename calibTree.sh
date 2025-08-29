#!/bin/bash

homeDir=/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf

# for submitting jobs with slurm
runNb=$1
segNb=$2

fileNAME=$homeDir/temp/calibTree_${runNb}_${segNb}

root -b <<EOF
    .L $fileNAME.C+
        calibTree($runNb,$segNb);
    .q
EOF

# for running locally
#run_list=$1 # Path to the file list

# Loop through each file in the file list
#while infile= read -r file; do
#    root -b <<EOF
#    .L calibTree.C+
#        calibTree($file, -1);
#    .q
#EOF
#done < "$run_list"