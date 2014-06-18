import FWCore.ParameterSet.Config as cms
import os
import re

process = cms.Process("FWLitePlots")

baseAddFiles = "/gpfs/gpfsddn/cms/user/arizzi/V11_bis/CMSSW_4_2_8_patch3/src/VHbbAnalysis/VHbbDataFormats/bin/" 

f = open(os.environ.get("FILETOPROCESS"))
lines = f.readlines()
f.close()

#fileNames   = cms.vstring('file:2l2bMetEdmNtuples.root'),         ## mandatory
process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(),
    PUmcfileName2011B= cms.string(baseAddFiles+"Fall11_Generated.root"),
    PUdatafileName2011B = cms.string(baseAddFiles+"Cert_175832-180252_PromptReco_JSON.pileupTruth_v2_finebin.root"),
    PUmcfileName = cms.string(baseAddFiles+"ttbarPU_35bins_fall11.root"),
    PUdatafileName = cms.string(baseAddFiles+"Pileup_2011_to_173692_LPLumiScale_68mb_35bins.root"),
    Weight3DfileName = cms.string(baseAddFiles+"Weight3D_Fall.root"),
    runMin  = cms.int32(-1),
    runMax  = cms.int32(-1),
    maxEvents   = cms.int32(-1),                             ## optional
    outputEvery = cms.uint32(0),                            ## optional
    skipEvents   = cms.int32(0),                             ## optional
    )


channel =  re.sub(".txt","",os.environ.get("FILETOPROCESS"))


for l in lines :
   process.fwliteInput.fileNames.append(re.sub("\n","",l))   
#   process.fwliteInput.fileNames.append("dcap://cmsdcache" + dirname + "/" + basename)
print "Number of files to process is %s" %(len(process.fwliteInput.fileNames)) 
    
    



#

fname = 'DiJetPt_' + channel + '.root'

process.fwliteOutput = cms.PSet(
    
    fileName  = cms.string(fname),## mandatory
    )

process.Analyzer = cms.PSet(
    triggers = cms.vstring(
        "HLT_IsoMu17_v.*" , #0
        "HLT_DoubleMu7_v.*", #1
        "HLT_Mu13_Mu8_v.*", #2
        "HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v.*", #3
        "HLT_Ele27_WP80_PFMHT50_v.*", #4
        "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v.*", #5
        "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v.*", #6
        "HLT_DiCentralJet20_BTagIP_MET65_v.*", #7
        "HLT_MET120_v.*", #8
        "HLT_CentralJet80_MET80_v.*", #9
        "HLT_PFMHT150_v.*", #10
        "HLT_DiCentralJet20_MET80_v.*", #11
        "HLT_DiCentralJet20_MET100_HBHENoiseFiltered_v.*", #12
        "HLT_IsoMu20_v.*", #13
        "HLT_IsoMu24_v.*", #14
        "HLT_IsoMu30_eta2p1_v.*", #15
        "HLT_Mu17_Mu8_v.*", #16
        "HLT_Ele17_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25_PFMHT15_v.*", #17
        "HLT_Ele22_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25_PFMHT20_v.*", #18
        "HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25_PFMHT20_v.*", #19
        "HLT_Mu30_v.*", #20 
        "HLT_Mu40_v.*", #21
        "HLT_Mu40_eta2p1_v.*", #22
        "HLT_IsoMu24_eta2p1_v.*", #23
        "HLT_IsoMu17_eta2p1_DiCentralJet30_v.*", #24
        "HLT_IsoMu17_eta2p1_DiCentralPFJet25_PFMHT15_v.*", #25
        "HLT_Ele30_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_DiCentralJet30_PFMHT25_v.*", #26
        "HLT_Ele27_WP80_DiCentralPFJet25_PFMHT15_v.*" #27 

  ),
    isMC =     cms.bool(True),
    verbose = cms.bool(False),
    readFromCandidates = cms.bool(False),
    jetPtThresholdZ = cms.double(20),
    jetPtThresholdW = cms.double(20),
    bJetCountThreshold = cms.double(0.898),
    useHighestPtHiggsW = cms.bool(True),
    useHighestPtHiggsZ = cms.bool(True),
    idMuFileName = cms.string(baseAddFiles+"ScaleEffs42.root"),
    hltMuFileName = cms.string(baseAddFiles+"ScaleFactor_muonEffsOnlyIsoToHLT2.2fb_efficiency.root"),
    hltEle1FileName = cms.string(baseAddFiles+"Ele17.root"),
    hltEle2FileName = cms.string(baseAddFiles+"Ele8NotEle17.root"),
    hltEle1AugFileName = cms.string(baseAddFiles+"Ele17Aug5PromptRecoV6.root"),
    hltEle2AugFileName = cms.string(baseAddFiles+"Ele8NotEle17Aug5PromptRecoV6.root"),
    idEle80FileName = cms.string(baseAddFiles+"PFElectronToWP80.root"),
    idEle95FileName = cms.string(baseAddFiles+"PFElectronToWP95.root"),
    hltJetEle1FileName = cms.string(baseAddFiles+"TriggerEfficiency_Jet30_PromptV4Aug05PromptV6.root"),
    hltJetEle2FileName = cms.string(baseAddFiles+"TriggerEfficiency_JetNo30_Jet25_PromptV4Aug05PromptV6.root"),
    recoEleFileName = cms.string(baseAddFiles+"EleReco.root"),
    hltSingleEleMayFileName = cms.string(baseAddFiles+"TriggerEfficiency_Electrons_May10.root"),
    hltSingleEleV4FileName = cms.string(baseAddFiles+"TriggerEfficiency_Electrons_PromptV4Aug05PromptV6.root"),
    idEleFileName = cms.string(baseAddFiles+"ScaleFactor_PFElectrons_DataMontecarlo.root"),
    hltMuOr30FileName =  cms.string(baseAddFiles+"ScaleFactor_muonEffsIsoToHLTMu30NotIso_efficiency.root"),
    btagEffFileName = cms.string(baseAddFiles+"btag_generic.txt")
    )





    
  
    

