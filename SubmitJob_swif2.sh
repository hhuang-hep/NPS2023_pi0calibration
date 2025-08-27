#!/bin/bash

# Run list for the jobs ===========================================================
listFile=x60_4b_-1_ProdList.txt
listDir="/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/MakeProdList"
listFileName=$listDir/$listFile

if [[ ! -f "$listFileName" ]]; then
  echo "Can't find the file: $listFileName"
  exit 1
fi

# Workflow settings ===========================================================

# Which workflow to run the jobs in
Workflow=pi0CalibTree_x60_4b_pass2

# Input, output rootfiles and slurm job names for each run and segment 
declare -A local_input_1
declare -A remote_input_1
declare -A local_output_1
declare -A remote_output_1
declare -A Slurm_job_name

# Remote input file directories
remote_input_1_Dir=mss:/mss/hallc/c-nps/users/kerver/WF/Kin60_4b_Kin60_2
# Note: the default input directory is $PWD for the macros in the same place of this script

# Remote output file directories
remote_output_1_Dir=$PWD/calibTree/x60_4b

while read -r run nseg; do
    for ((iseg=0; iseg<nseg; iseg++)); do
        key="${run},${iseg}"
        # Input file for each job
        local_input_1[$key]=nps_production_${run}_${iseg}_4_wf.root
        remote_input_1[$key]=${remote_input_1_Dir}/${local_input_1[$key]}
        # Output file for each job
        local_output_1[$key]=prodTree_pass2_${run}_${iseg}.root
        remote_output_1[$key]=${remote_output_1_Dir}/${local_output_1[$key]}
        # Slurm job name for each job
        Slurm_job_name[$key]=${Workflow}_${run}_${iseg}
    done
done < "$listFileName"

# Input scripts
local_input_2=rootlogon.C
remote_input_2=$PWD/$local_input_2

local_input_3=calibTree.C
remote_input_3=$PWD/$local_input_3


# The job script to run
Job_script=$PWD/JobScript_swif2.sh


# Globel settings for slurm jobs ===========================================================
Slurm_account=hallc
Slurm_partition=production
Slurm_constraint=el9
Slurm_cores=1
Slurm_disk_scratch=41000000000
Slurm_ram=2500000000
Slurm_time=2days
Slurm_shell=/bin/bash

# The followint commend lines are used to submit and run the jobs
# "./SubmitJob_swif2.sh": check the commends
# "./SubmitJob_swif2.sh add": add jobs to the workflow
# "./SubmitJob_swif2.sh run": start the jobs
# "./SubmitJob_swif2.sh status": check the job and workflow status
# Note: if the workflow does not exist, it will be created automatically

# Check the settings and add the jobs
ECHO="echo"   ## Used to echo swif commands for testing
if [ ${1:-unset} == "add" ]; then
    ECHO=""  ## disables the 'echo' and actually runs the swif commands
fi
# Run the jobs
if [ ${1:-unset} == "run" ]; then
    echo "Start the jobs in workflow $Workflow"
    swif2 run -workflow $Workflow
    swif2 notify $Workflow -when done -email h19901027@gmail.com
    swif2 notify $Workflow -when stalled -email h19901027@gmail.com
    exit 0
fi
# Check the status
if [ ${1:-unset} == "status" ]; then
    swif2 status -workflow $Workflow -jobs
    exit 0
fi

# Abandon ALL the jobs in the workflow
if [ ${1:-unset} == "abandon" ]; then
    swif2 abandon-jobs $Workflow -names -regexp '.*'
    exit 0
fi

# Delete the workflow (you should abandon all jobs first to clean up)
if [ ${1:-unset} == "delete" ]; then
    if [ -e "${Workflow}_status.txt" ]; then
        rm ${Workflow}_status.txt # Remove the status file if it exists
    fi
    swif2 status -workflow $Workflow -jobs >> ${remote_output_1_Dir}/${Workflow}_status.txt # Dump the status for debugging. Can be found in the output directory
    swif2 abandon-jobs $Workflow -names -regexp '.*'
    swif2 cancel $Workflow -delete
    exit 0
fi

# Commends to add the jobs to the workflow
while read -r run nseg; do
    for ((iseg=0; iseg<nseg; iseg++)); do
        key="${run},${iseg}"
        $ECHO \
        swif2 add-job \
            -create \
            -workflow "$Workflow" \
            -account "$Slurm_account" \
            -partition "$Slurm_partition" \
            -constraint "$Slurm_constraint"  \
            -cores "$Slurm_cores"  \
            -disk-scratch "$Slurm_disk_scratch"  \
            -ram "$Slurm_ram"  \
            -time "$Slurm_time"  \
            -shell "$Slurm_shell"  \
            -name "${Slurm_job_name[$key]}" \
            -input "${local_input_1[$key]}" "${remote_input_1[$key]}" \
            -input "$local_input_2" "$remote_input_2" \
            -input "$local_input_3" "$remote_input_3" \
            -output "${local_output_1[$key]}" "${remote_output_1[$key]}" \
            -stdout "/farm_out/hhuang/$Workflow/${Slurm_job_name[$key]}.out" \
            -stderr "/farm_out/hhuang/$Workflow/${Slurm_job_name[$key]}.err" \
            "$Job_script" "$run" "$iseg"
        
        if [ "$ECHO" == "echo" ]; then
            echo ""
        fi
    done
done < "$listFileName"

if [ ${1:-unset} == "submit" ]; then
    swif2 status -workflow $Workflow -jobs
fi