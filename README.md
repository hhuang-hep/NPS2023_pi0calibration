# pi0 calibration scripts (NPS2023 rg1a)
## Introduction
This is the scripts for the energy calibration for the NPS calorimeter. \
The calibration method is based on Nucl. Instrum. Methods A, 566, 366–374 (2006) using the measured pi0.

## Required enviroment
 - root-6.30.04 in /group/halla/modulefiles
 - NPS software for clustering and photon reconstruction

## Download and compiling of NPS software
- Setup environment\
`module purge`\
`module use /group/halla/modulefiles`\
`module load root/6.30.04`

- Compile NPS software\
`git clone https://github.com/hhuang-hep/NPS_SOFT.git`\
`cd NPS_SOFT`\
`make biglib`\
`cd ..`

## Environment setup for the calibration
`module purge`\
`module use /group/halla/modulefiles`\
`module load root/6.30.04`

`setenv NPS_SOFT <path to your NPS software folder>` \
`setenv LD_LIBRARY_PATH ${NPS_SOFT}:${LD_LIBRARY_PATH}` \
`setenv PATH ${NPS_SOFT}:${PATH}`


## Calibration procedure
1. Saving small trees (per run per segment)
    #### Script and output file
    - Scripts: calibTree.C
    - Clustering using pulses within (-3, 3) ns and a 2x2 clustering threshold of 0.2 GeV
    - Output branches
        | Branch name | Type | Discription |
        |:------|:------|:------|
        | g.runnum | Double_t | Run number from T tree|
        | g.evnum | Double_t | Global event number from T tree|
        | H_react_x | Double_t | Vertex X positions|
        | H_react_y | Double_t | Vertex Y positions|
        | H_react_z | Double_t | Vertex Z positions|
        | H_gtr_px | Double_t | Scattered electron momentum X|
        | H_gtr_py | Double_t | Scattered electron momentum Y|
        | H_gtr_pz | Double_t | Scattered electron momentum Z|
        | m | Double_t | Reconstructed two-photon invariant mass|
        | mm2 | Double_t | Reconstructed two-photon missing mass|
        | caloev | TCaloEvent | Cluster information (Cluster energy & block energy in GeV)|

    #### Required modification to calibTree.C
    - Modify the following two lines on the top
        - #include "<absolute path to your pi0 calibration folder\>/calibHeader/Analysis.h"
        - #include "<absolute path to your NPS software folder\>/TDVCSDB.h"

    #### Test before running on the farm
    1. Change the directory pointed to the input file by modifying the line TString dataDir = "." in calibTree.C
    2. Currently (2026/09), the files after waveform-fitting are under /cache/hallc/c-nps/analysis/pass2/WF
    3. Check the available run and segment under that directory. 
    4. Excute `root -b -q 'calibTree.C(<run number>,<segment number>)` and check the output ROOT file

    #### Note
    - The 0.2 GeV clustering threshold can be change by modifying the line "Double_t clusThr = 0.2" in calibTree.C.
    - A boolean flag in calibTree.C, "bool use_wf = true" is set to use the waveform-fitted amplitude.
    - Only output the events with 2 clusters with a two-photon invariant mass within (0.05, 0.2).
    - 10 trees with different photon energy cut (larger than) from 0.5 to 1.5 GeV in step of 0.1 GeV are output. 
    - Names of the output trees are t_prod_<iEcut\>, where <iEcut\> = (photon energy cut - 0.5)/0.1 is a integer.