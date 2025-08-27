#!/bin/bash

listDir="/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/MakeProdList"
listFileName=$listDir/$1

if [[ ! -f "$listFileName" ]]; then
  echo "Can't find the file: $listFileName"
  exit 1
fi

while read -r run nseg; do
  for ((iseg=0; iseg<nseg; iseg++)); do
    sbatch JobScript.sh $run $iseg
  done
done < "$listFileName"