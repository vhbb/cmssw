#!/bin/sh
#JSON=newJSON_max196509.txt
JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-202305_8TeV_PromptReco_Collisions12_JSON.txt
PUJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-202478_corr.txt
MINBIAS=69400
MINBIASP=66475
MINBIASM=72450
#MINBIAS=71000

pileupCalc.py -i ${JSON} --inputLumiJSON $PUJSON --calcMode true --minBiasXsec $MINBIAS --maxPileupBin 60 --numPileupBins 60  outputfile.root
pileupCalc.py -i ${JSON} --inputLumiJSON $PUJSON --calcMode true --minBiasXsec $MINBIASP --maxPileupBin 60 --numPileupBins 60  outputfileP.root
pileupCalc.py -i ${JSON} --inputLumiJSON $PUJSON --calcMode true --minBiasXsec $MINBIASM --maxPileupBin 60 --numPileupBins 60  outputfileM.root
