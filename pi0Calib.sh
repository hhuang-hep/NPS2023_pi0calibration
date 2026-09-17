#!/bin/bash

# This program is used for pi0 calibration
# execute with ./pi0Calib.sh <kinematics name (w/o KinC_x)> <target flag> <number of cycles to run> <cycle to start> <iteration to start>
# for example: ./pi0Calib.sh 36_5_3 0 2 0 1 execute for 
#   - KinC_x36_5_3, 
#   - LH2 target
#   - combining 2 cycles in a group
#   - starting from cycle 0
#   - starting from iteration 1
# see /group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/MakeProdList/<kinematics name>_<target flag>_ProdList.txt
# Target flag: 0 for LH2, 1 for LD2, -1 for both

Kine=$1
Tar=$2
group=$3
minCycle=$4
minIter=$5
maxIter=8
runlistDir="/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/calibList/ListPerCycle"

if [ "$group" -le 0 ]; then
    echo "ERROR: group must be > 0, got group=$group"
    exit 1
fi

# Find the number of cycles
nCycle=$(ls ${runlistDir}/x${Kine}_${Tar}_cycle*.txt | wc -l)

if [ -z "$nCycle" ]; then
    echo "ERROR: No cycle files found for $Kine"
    exit 1
fi

echo "Number of cycle for KinC_x$Kine is $nCycle"

# uncomment 2 lines below for statistic test__________________________________
# nCycle=$(($minCycle+$group)) # this line replace the number of cycles to run
# echo "For test, run only 1 group with $nCycle cycles"

# run the calibration
for (( firstcycle=minCycle; firstcycle<nCycle; firstcycle+=group ));do
    # echo "DEBUG: firstcycle=$firstcycle / nCycle=$nCycle / group=$group"
    for i_Iter in $(seq $minIter $((maxIter+1)));do
        echo "Start iteration $i_Iter for KinC_$Kine with $group cycles in a group from cycle $firstcycle"
        root -b <<EOF
            .L pi0Calib.C
                pi0Calib("${Kine}", ${Tar}, ${firstcycle}, ${group}, ${i_Iter})
            .q
EOF
done

# make plots and the coefficients file
root -b <<EOF
    .L Draw_pi0Calib.C
        Draw_pi0Calib("${Kine}", ${Tar}, ${firstcycle}, ${group})
    .q
EOF
done