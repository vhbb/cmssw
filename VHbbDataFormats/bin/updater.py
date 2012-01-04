import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("FWLitePlots")

#fileNames   = cms.vstring('file:2l2bMetEdmNtuples.root'),         ## mandatory
fname =  os.environ.get("FILETOPROCESS")

#fname = "DiJetPt_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root"
#name = "DiJetPt_DYJetsToLL_PtZ-100_TuneZ2_7TeV-madgraph-tauola.root"

process.fwliteInput = cms.PSet(
    fileName   = cms.string(fname),
    PUmcfileName2011B= cms.string("Fall11_Generated.root"),
    PUdatafileName2011B = cms.string("Cert_175832-180252_PromptReco_JSON.pileupTruth_v2_finebin.root"),
    PUmcfileName = cms.string("ttbarPU_35bins.root"),
    PUdatafileName = cms.string("Pileup_2011_to_173692_LPLumiScale_68mb_35bins.root"),
    Weight3DfileName = cms.string(""),

#    Weight3DfileName = cms.string(""),
    maxEvents   = cms.int32(-1),                             ## optional
    skipEvents   = cms.int32(0),                             ## optional
    outputEvery = cms.uint32(0),                            ## optional
    )

fnameOut = "Updated_"+fname
process.fwliteOutput = cms.PSet(
    fileName  = cms.string(fnameOut),## mandatory
)

process.Analyzer = cms.PSet(
    replaceWeights = cms.bool(True),
    redoPU = cms.bool(True),
    redoTrigger = cms.bool(False),
    redoHiggs = cms.bool(False),
    idMuFileName = cms.string("ScaleEffs42.root"),
    hltMuFileName = cms.string("ScaleFactor_muonEffsOnlyIsoToHLT2.2fb_efficiency.root"),
    hltEle1FileName = cms.string("Ele17.root"),
    hltEle2FileName = cms.string("Ele8NotEle17.root"),
    hltEle1AugFileName = cms.string("Ele17Aug5PromptRecoV6.root"),
    hltEle2AugFileName = cms.string("Ele8NotEle17Aug5PromptRecoV6.root"),
    idEle80FileName = cms.string("PFElectronToWP80.root"),
    idEle95FileName = cms.string("PFElectronToWP95.root"),
    hltJetEle1FileName = cms.string("TriggerEfficiency_Jet30_PromptV4Aug05PromptV6.root"),
    hltJetEle2FileName = cms.string("TriggerEfficiency_JetNo30_Jet25_PromptV4Aug05PromptV6.root"),
    recoEleFileName = cms.string("EleReco.root"),
    hltSingleEleMayFileName = cms.string("TriggerEfficiency_Electrons_May10.root"),
    hltSingleEleV4FileName = cms.string("TriggerEfficiency_Electrons_PromptV4Aug05PromptV6.root"),
    idEleFileName = cms.string("ScaleFactor_PFElectrons_DataMontecarlo.root"),
    hltMuOr30FileName =  cms.string("ScaleFactor_muonEffsIsoToHLT2.2fb_efficiency.root"),
    btagEffFileName = cms.string("btag_generic.txt")
    )

    
  
    

