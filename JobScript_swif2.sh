#!/bin/bash

source /etc/skel/.bashrc
source /group/nps/hhuang/software/SetEnvi.sh

runNb=$1
segNb=$2

root -b <<EOF
    .L calibTree.C+
        calibTree($runNb,$segNb);
    .q
EOF