#!/bin/bash
cmsenv
source /cvmfs/cms.cern.ch/crab3/crab.sh
#voms-proxy-init -voms cms
for i in `cat datasets.txt` ; do
  export DATASET=$i 
  crab submit -c heppy_crab_config_env.py 
done 
