#!/bin/bash

Kine=$1
Tar=$2
minIter=$3
maxIter=8

#create folder
# if [[ -e 'Result/'$runNb ]]; then
#     echo "Start the pi0 calibration for run" $runNb
# else
#     mkdir Result/$runNb
#     echo "Folder of run" $runNb " have been created"
#     echo "Start the pi0 calibration"
# fi

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