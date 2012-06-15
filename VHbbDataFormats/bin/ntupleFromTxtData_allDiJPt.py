import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes

import FWCore.ParameterSet.Config as cms
import os
import re

process = cms.Process("FWLitePlots")

baseAddFiles = "/gpfs/gpfsddn/cms/user/arizzi/Hbb/CMSSW_5_2_5/src/VHbbAnalysis/VHbbDataFormats/bin/" 

f = open(os.environ.get("FILETOPROCESS"))
lines = f.readlines()
f.close()

#fileNames   = cms.vstring('file:2l2bMetEdmNtuples.root'),         ## mandatory
process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    PUmcfileName2011B= cms.string(baseAddFiles+"Summer12MCObserved.root"),
    PUdatafileName2011B = cms.string(baseAddFiles+"MyDataPileupHistogramObserved.root"),
    PUmcfileName = cms.string(baseAddFiles+"Summer12MCTrue.root"),
    PUdatafileName = cms.string(baseAddFiles+"Summer12DataTrue.root"),
#    Weight3DfileName = cms.string(baseAddFiles+"Weight3D_Summer12.root"),
    Weight3DfileName = cms.string(""),
    runMin  = cms.int32(-1),
    runMax  = cms.int32(-1),
    maxEvents   = cms.int32(-1),                             ## optional
    outputEvery = cms.uint32(0),                            ## optional
    skipEvents   = cms.int32(0),                             ## optional
    )


JSONfile = '/gpfs/gpfsddn/cms/user/arizzi/Hbb/CMSSW_5_2_5/src/VHbbAnalysis/VHbbDataFormats/bin/Cert_190456-195775_8TeV_PromptReco_Collisions12_JSON.txt'
lumiList = LumiList.LumiList (filename = JSONfile).getCMSSWString().split(',')
process.fwliteInput.lumisToProcess.extend(lumiList)

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
        "HLT_Ele27_WP80_DiCentralPFJet25_PFMHT15_v.*", #27
        "HLT_IsoPFTau35_Trk20_v.*", #28
        "HLT_IsoPFTau35_Trk20_MET45_v.*", #29
        "HLT_IsoPFTau35_Trk20_MET60_v.*", #30
        "HLT_IsoPFTau45_Trk20_MET60_v.*", #31
        "HLT_IsoPFTau35_Trk20_MET70_v.*", #32
        "HLT_MediumIsoPFTau35_Trk20_v.*", #33
        "HLT_MediumIsoPFTau35_Trk20_MET60_v.*", #34
        "HLT_MediumIsoPFTau35_Trk20_MET70_v.*", #35
        "HLT_LooseIsoPFTau35_Trk20_v.*", #36
        "HLT_LooseIsoPFTau35_Trk20_MET70_v.*", #37
        "HLT_LooseIsoPFTau35_Trk20_MET75_v.*", #38
        "HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v*", #39
        "HLT_DiCentralJet20_CaloMET65_BTagCSV07_PFMHT80*", #40
        "HLT_DiCentralPFJet30_PFMET80_BTagCSV07*", #41
        "HLT_PFMET150_v*", #42
        "HLT_L1ETM40_v*", #43
        "HLT_Ele27_WP80_v.*", #44
        "HLT_Ele27_WP80_WCandPt80_v.*", #45
        "HLT_IsoMu20_eta2p1_WCandPt80_v.*", #46
        "HLT_IsoMu20_WCandPt80_v.*", #47
        "HLT_Mu17_TkMu8_v.*", #48
	"HLT_DiCentralPFJet30_PFMET80_v*", #49

  ),
    isMC =     cms.bool(False),
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
    hltSingleEle2012Awp95 = cms.string(baseAddFiles+"triggerRootFiles/SingleEle.TrigEff.wp95.2012A.root"),
    hltSingleEle2012Awp80 = cms.string(baseAddFiles+"triggerRootFiles/SingleEle.TrigEff.wp80.2012A.root"),
    hltSingleMuon2012A = cms.string(baseAddFiles+"triggerRootFiles/SingleMu24OR40.TrigEff.2012A.root"),
    hltDoubleEle2012A_leg8 = cms.string(baseAddFiles+"triggerRootFiles/DoubleEle8.TrigEff.wp95.2012A.root"),
    hltDoubleEle2012A_leg17 = cms.string(baseAddFiles+"triggerRootFiles/DoubleEle17.TrigEff.wp95.2012A.root"),
    hltDoubleMuon2012A_leg8 = cms.string(baseAddFiles+"triggerRootFiles/DoubleMu8.TrigEff.2012A.root"),
    hltDoubleMuon2012A_leg17 = cms.string(baseAddFiles+"triggerRootFiles/DoubleMu17.TrigEff.2012A.root"),
    hltMuPlusWCandPt2012A_legMu = cms.string(baseAddFiles+"triggerRootFiles/SingleMu20Not24Or40.TrigEff.2012A.root"),
    hltMuPlusWCandPt2012A_legW = cms.string(baseAddFiles+"triggerRootFiles/WCandPt.TrigEff.2012A.root"),
    hltDoubleMuon2012A_dZ = cms.string(baseAddFiles+"triggerRootFiles/DoubleMuDz.TrigEff.2012AB.root"),
    hltDoubleEle2012A_dZ = cms.string(baseAddFiles+"triggerRootFiles/DoubleEleDz.TrigEff.2012AB.root"),
    csvDiscr = cms.string(baseAddFiles+"csvdiscr.root"),
    jecFolder = cms.string(baseAddFiles+"jec"),
    btagEffFileName = cms.string(baseAddFiles+"btag_generic.txt")
    )





    
  
    

