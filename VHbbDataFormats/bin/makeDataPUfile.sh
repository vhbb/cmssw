#!/bin/sh
JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-202305_8TeV_PromptReco_Collisions12_JSON.txt
PUJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-202478_corr.txt
#INBIAS=69400
MINBIAS=71000

pileupCalc.py -i ${JSON} --inputLumiJSON $PUJSON --calcMode true --minBiasXsec $MINBIAS --maxPileupBin 60 --numPileupBins 60  outputfile.root
