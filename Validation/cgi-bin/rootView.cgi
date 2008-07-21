#!/bin/bash

ROOTVER=5.12.00
if [ uname -a | grep -q x86_64 ]; then
	ARCH=slc4_amd64_gcc345
else
	ARCH=slc4_ia32_gcc345
fi

export ROOTSYS=/afs/cern.ch/sw/lcg/external/root/$ROOTVER/$ARCH/root
export PATH=/bin:/usr/bin:$ROOTSYS/bin
export LD_LIBRARY_PATH=$ROOTSYS/lib
export PYTHONPATH=$LD_LIBRARY_PATH

python rootView.py
RETVAL=$?

return $RETVAL

