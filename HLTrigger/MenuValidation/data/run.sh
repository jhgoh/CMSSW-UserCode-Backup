#!/bin/bash

cmsRun -j $RUNTIME_AREA/crab_fjr_$NJob.xml -p pset.py
echo =====================================================
ls -l output*.root
echo =====================================================

mkdir online

for FNAME in outputA outputHLTMON outputExpress; do
  ln -sf $FNAME.root inputfile.root
  cmsRun onlineDQM_cfg.py
  FTYPE=`echo $FNAME | sed 's/output//g'`

  for DQMOUTPUT in DQM_*.root; do
    RUNNR=`echo $DQMOUTPUT | sed 's/DQM.\+_R\([0-9]\+\).root/\1/g'`
    mv $DQMOUTPUT online/DQM_V0001_${FTYPE}_R${RUNNR}.root
  done

#  rm -f $FNAME.root
done

for FNAME in outputHLTDQM; do
  ln -sf $FNAME.root inputfile.root
  cmsRun convert_cfg.py
  cmsRun onlineDQM_dat_cfg.py
  FTYPE=`echo $FNAME | sed 's/output//g'`

  for DQMOUTPUT in DQM_*.root; do
    RUNNR=`echo $DQMOUTPUT | sed 's/DQM.\+_R\([0-9]\+\).root/\1/g'`
    mv $DQMOUTPUT online/DQM_V0001_${FTYPE}_R${RUNNR}.root
  done

#  rm -f $FNAME.root
  rm -f inputfile.dat
done

rm -f inputfile.root
tar czf online.tgz online
rm -f *.root

