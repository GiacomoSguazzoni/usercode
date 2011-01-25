#!/bin/bash

export SUFFIX=$2
export CRABDIR=$1
export INPUTFILENAME=Tracking_PFG_*_*.root
export OUTPUTFILENAME=Tracking_PFG_${SUFFIX}.root

hadd rootfiles/${OUTPUTFILENAME} ${CRABDIR}/res/${INPUTFILENAME}
#hadd /tmp/venturia/${OUTPUTFILENAME} ${CRABDIR}/res/${INPUTFILENAME}

