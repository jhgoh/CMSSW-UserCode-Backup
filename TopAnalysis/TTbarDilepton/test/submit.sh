#!/bin/bash

MAXFILES=50
DSFILES=`\ls ../data/dataset-*.txt`

cat > run.sh <<EOF
#!/bin/bash

cd `pwd`
eval \`scram runtime -sh\`

export DATASET=\$1
export SECTION=\$2
export MAXFILES=\$3

cmsRun ntuple_cfg.py

EOF

chmod +x run.sh

if [ ! -d log ]; then
  mkdir log
fi

if [ ! -d ntuple ]; then
  mkdir ntuple
fi

if [ ! -d ntuple/unmerged ]; then
  mkdir ntuple/unmerged
fi

for DSFILE in $DSFILES; do
  DATASET=`basename $DSFILE | sed -e 's/dataset-//g' -e 's/.txt//g'`
  NFILES=`cat $DSFILE | grep .root | grep -v '^#' | grep -v '^ *$' | wc -l`
  NJOBS=$(($NFILES/$MAXFILES))
  if [ $(($NJOBS*$MAXFILES)) -lt $NFILES ]; then
    NJOBS=$(($NJOBS+1))
  fi

  echo Submitting $DATASET,
  echo     NFILES=$NFILES
  echo     NJOBS=$NJOBS
  echo     MAXFILES=$MAXFILES

  for SECTION in `seq 0 $(($NJOBS-1))`; do
    #DATASET=$DATASET SECTION=$SECTION MAXFILES=$MAXFILES python ntuple_cfg.py
    bsub -q 8nh -oo log/${DATASET}_${SECTION}.log run.sh $DATASET $SECTION $MAXFILES
  done
  echo
done
