#!/bin/bash

## Workarount to retrive PD name and menu name
## This workarount should be changed if script_arguments option 
## in the crab is working properly
#DATASETNAME=$1
#MENUNAME=$2

DATASETNAME=`cat CMSSW.sh | grep -m1 -o 'PrimaryDataset=\(.*\)' | sed -e 's/PrimaryDataset=//g'`
MENUNAME=NewMenuTest

cmsRun -j $RUNTIME_AREA/crab_fjr_$NJob.xml -p pset.py
echo =====================================================
ls -l output*.root
echo =====================================================

mkdir online

for FNAME in outputA outputHLTMON outputExpress; do
  ln -sf $FNAME.root inputfile.root
  FTYPE=`echo $FNAME | sed 's/output//g'`
  export WORKFLOW="/${DATASETNAME}/${MENUNAME}/online${FTYPE}"

  cmsRun onlineDQM_cfg.py

  mv DQM_*.root online/
#  rm -f $FNAME.root
done

for FNAME in outputHLTDQM; do
  ln -sf $FNAME.root inputfile.root
  FTYPE=`echo $FNAME | sed 's/output//g'`
  export WORKFLOW="/${DATASETNAME}/${MENUNAME}/online${FTYPE}"

  cmsRun convert_cfg.py
  cmsRun onlineDQM_dat_cfg.py

  mv DQM_*.root online/
#  rm -f $FNAME.root
  rm -f inputfile.dat
done

rm -f inputfile.root
tar czf online.tgz online
rm -f *.root

