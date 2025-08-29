#!/bin/bash

# This program is used for pi0 calibration
# execute with ./pi0Calib.sh <kinematics name> <target flag> <iteration to start>
# see /group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/MakeProdList/<kinematics name>_<target flag>_ProdList.txt
# Target flag: 0 for LH2, 1 for LD2, -1 for both

Kine=$1
Tar=$2
minIter=$3
maxIter=8

#run the calibration
for i_Iter in $(seq $minIter $maxIter)
do  
    echo "Start iteration " $i_Iter " for KinC_"$Kine
    root -b <<EOF
        .L pi0Calib.C+
            pi0Calib("$Kine",$Tar,$i_Iter)
        .q
EOF
done

#make plots and the coefficients file
root -b <<EOF
    .L Draw_pi0Calib.C
        Draw_pi0Calib("$Kine",$Tar)
    .q
EOF

#output coefficients of each run
root -b <<EOF
    .L outputCoef_PerRun.C
        outputCoef_PerRun("$Kine",$Tar)
    .q
EOF