#!/bin/bash

## Workarount to retrive PD name and menu name
## This workarount should be changed if script_arguments option 
## in the crab is working properly
#DATASETNAME=$1
#MENUNAME=$2
DATASETNAME=`cat ${RUNTIME_AREA}/CMSSW.sh | grep -m1 -o 'PrimaryDataset=\(.*\)' | sed -e 's/PrimaryDataset=//g'`
MENUNAME=NewMenuTest

mkdir log

## Start HLT-Rerunning
cmsRun -j $RUNTIME_AREA/crab_fjr_$NJob.xml -p pset.py >& log/hlt.log
echo =====================================================
ls -l output*.root
echo =====================================================

## Test online DQM
mkdir online

for FNAME in outputA outputHLTMON outputExpress; do
  ln -sf $FNAME.root inputfile.root
  FTYPE=`echo $FNAME | sed 's/output//g'`
  export WORKFLOW="/${DATASETNAME}/${MENUNAME}/online${FTYPE}"

  cmsRun onlineDQM_cfg.py >& log/onlineDQM_${FTYPE}.log

  ## Keep DQM output to the safe place
  mv DQM_*.root online/
done

for FNAME in outputHLTDQM; do
  ln -sf $FNAME.root inputfile.root
  FTYPE=`echo $FNAME | sed 's/output//g'`
  export WORKFLOW="/${DATASETNAME}/${MENUNAME}/online${FTYPE}"

  cmsRun convert_cfg.py >& log/onlineDQM_convert_${FTYPE}.log
  cmsRun onlineDQM_dat_cfg.py >& log/onlineDQM_${FTYPE}.log
  rm -f inputfile.dat

  ## Keep DQM output to the safe place
  mv DQM_*.root online/
done

## Do the reco+offlineDQM step
mkdir offline

for FNAME in outputA outputHLTMON outputExpress; do
  FTYPE=`echo $FNAME | sed 's/output//g'`
  cmsDriver.py step2 -s RAW2DIGI,RECO,DQM -n -1 --eventcontent DQM \
    --data --conditions auto:com10 --geometry Ideal \
    --filein file:${FNAME}.root --fileout step2_${FTYPE}.root --no_exec --python_filename=step2_cfg.py
  ## Temporary remove jetMETHLTOfflineSource due to DQM crash
  cat >> step2_cfg.py <<EOF
process.DQMOffline.remove(process.jetMETHLTOfflineSource)
EOF
  cmsRun step2_cfg.py >& log/offlineDQM_${FTYPE}.log
#  cmsDriver.py step3 -s HARVESTING:dqmHarvesting \
#    --data --conditions auto:com10 \
#    --filein file:step2_${NAME}.root
  mv step2_${FTYPE}.root offline/
done

## Everything done. Clear up output
rm -f inputfile.root
tar czf log.tgz log && rm -rf log
tar czf online.tgz online && rm -rf online
tar czf offline.tgz offline && rm -rf offline

## Make dummy files to save disk space
for FNAME in *.root; do
  rm -f $FNAME
  echo "DUMMY" > $FNAME
done

