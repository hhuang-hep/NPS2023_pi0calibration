#!/bin/bash

#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --job-name=pi0Calib
#SBATCH --output=/farm_out/%u/out/%x-%A-%a.out
#SBATCH --error=/farm_out/%u/err/%x-%A-%a.err
#SBATCH --partition=production
#SBATCH --account=hallc
#SBATCH --mail-user=hhuang@jlab.org
#SBATCH --gres=disk:5120
#SBATCH --mem-per-cpu=10000
#SBATCH --time=24:00:00

#Environment setting
source /group/nps/hhuang/software/SetEnvi.sh
homeDir=/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf
jobDir=/farm_out/hhuang

runNb=$1
segNb=$2

#Copy files and run the script
cp $homeDir/calibTree.C $homeDir/temp/calibTree_${runNb}_${segNb}.C
srun $homeDir/calibTree.sh $runNb $segNb