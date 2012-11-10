#!/bin/bash

MAXFILES=50
DSFILES=`\ls ../data/dataset-*.txt`

if [ ! -d ntuple/merged ]; then
  mkdir ntuple/merged
fi

for DSFILE in $DSFILES; do
  DATASET=`basename $DSFILE | sed -e 's/dataset-//g' -e 's/.txt//g'`

  hadd ntuple/merged/ntuple_$DATASET.root ntuple/unmerged/ntuple_${DATASET}_*.root
done
