#!/bin/sh
#JSON=newJSON_max196509.txt
JSON=Cert_190456-196531_8TeV_13Jul2012ReReco_ert_190782-190949_8TeV_06Aug2012ReReco_Cert_198022-198523_8TeV_24Aug2012ReReco_Cert_198941-203002_8TeV_PromptReco_Recover201191-201191_8TeV_11Dec2012ReReco__203002-208686_8TeV_PromptReco_Collisions12_JSON.txt
#/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-202305_8TeV_PromptReco_Collisions12_JSON.txt
#PUJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-202478_corr.txt
PUJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-208686_corr.txt

MINBIAS=69400
MINBIASP=66475
MINBIASM=72450
#MINBIAS=71000

pileupCalc.py -i ${JSON} --inputLumiJSON $PUJSON --calcMode true --minBiasXsec $MINBIAS --maxPileupBin 60 --numPileupBins 60  outputfile.root
pileupCalc.py -i ${JSON} --inputLumiJSON $PUJSON --calcMode true --minBiasXsec $MINBIASP --maxPileupBin 60 --numPileupBins 60  outputfileP.root
pileupCalc.py -i ${JSON} --inputLumiJSON $PUJSON --calcMode true --minBiasXsec $MINBIASM --maxPileupBin 60 --numPileupBins 60  outputfileM.root
