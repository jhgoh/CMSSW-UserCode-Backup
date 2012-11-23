#!/bin/bash

mkdir log ntuple

cat > run.sh <<EOF
#!/bin/bash

export SAMPLE=\$1
export SECTION=\$2
export NSECTION=\$3

cd `pwd`
eval \`scram runtime -sh\`

cmsRun leptonAnalysis_cfg.py
EOF
chmod +x run.sh

MAXFILES=100
#for SAMPLE in Run2012AElEl Run2012BElEl Run2012CElEl; do
#for SAMPLE in Run2012AMuMu Run2012BMuMu Run2012CMuMu; do
for DATASET in Run2012-Run2012AMuEl Run2012-Run2012BMuEl Run2012-Run2012CMuEl Summer12-ZJets Summer12-TTbarTuneZ2 Summer12-Tbartw Summer12-Ttw Summer12-WJets Summer12-WW Summer12-WZ Summer12-ZJets10To50 Summer12-ZZ; do
  PROD=`echo $DATASET | cut -d- -f1`
  SAMPLE=`echo $DATASET | cut -d- -f2`
  SOURCE="$CMSSW_BASE/src/KoPFA/CommonTools/python/Sources/CMG/V5_10_0/${PROD}/cmgTuple_${SAMPLE}_cff.py"
  NFILES=`cat $SOURCE | grep .root | grep -v '^#' | wc -l`
  NSECTION=$(($NFILES/$MAXFILES))
  if [ $(($NSECTION*$MAXFILES)) == $NFILES ]; then
    NSECTION=$(($NSECTION-1))
  fi

  for SECTION in `seq 0 $NSECTION`; do
    echo bsub -oo log/${SAMPLE}_${SECTION}.log -q 8nh run.sh $SAMPLE $SECTION $(($NSECTION+1))
  done
done
