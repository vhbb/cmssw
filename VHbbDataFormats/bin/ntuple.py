import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

process = cms.Process("FWLitePlots")

#fileNames   = cms.vstring('file:2l2bMetEdmNtuples.root'),         ## mandatory
process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(
"root://cmsdcache//pnfs/pi.infn.it/data/cms/store/user/tboccali/WH_WToLNu_HToBB_M-120_8TeV-powheg-herwigpp/HBB_EDMNtupleV30_ProcV1_WH_WToLNu_HToBB_M-120/14fe2b624ddea84f5c39709f51bf546f/PAT.edm_51_1_3LJ.root"
#/pnfs/pi.infn.it/data/cms/store/user/tboccali/ZH_ZToNuNu_HToBB_M-120_8TeV-powheg-herwigpp/HBB_EDMNtupleV30_ProcV1_ZH_ZToNuNu_HToBB_M-120/14fe2b624ddea84f5c39709f51bf546f/PAT.edm_51_1_ukR.root

),


    PUmcfileName = cms.string("ttbarPU_35bins_fall11.root"),
    PUmcfileName2011B= cms.string("Fall11_Generated.root"),
    PUdatafileName2011B = cms.string("Cert_175832-180252_PromptReco_JSON.pileupTruth_v2_finebin.root"),
    PUdatafileName = cms.string("Pileup_2011_to_173692_LPLumiScale_68mb_35bins.root"),
    Weight3DfileName = cms.string(""),
    maxEvents   = cms.int32(-1),                             ## optional
    runMin  = cms.int32(-1),
    runMax  = cms.int32(-1),
    skipEvents   = cms.int32(0),                             ## optional
    outputEvery = cms.uint32(0),                            ## optional
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    )

# get JSON file 
#JSONfile = 'Input.json'
#lumiList = LumiList.LumiList (filename = JSONfile).getCMSSWString().split(',')

#Uncomment to run with JSON
#process.fwliteInput.lumisToProcess.extend(lumiList)


channel =  "WH120"
import os
dirnameOld = "//pnfs/pi.infn.it/data/cms/store/user/bortigno/DoubleMu/HBB_EDMNtupleV3_ProcV1_may/07fb60889166b64f474d8d0aa162db69/"




#for i in range(len(channels)):
 

#dirname =  dirnameOld 
#dirlist = os.listdir(dirname)
#basenamelist = os.listdir(dirname + "/")
#for basename in basenamelist:
#   process.fwliteInput.fileNames.append("root://cmsdcache" + dirname + "/" + basename)
print "Number of files to process is %s" %(len(process.fwliteInput.fileNames)) 
    
    



#


fname = 'Test' + channel + '.root'

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
        "HLT_LooseIsoPFTau35_Trk20_MET75_v.*" #38
        
   ),
    isMC =     cms.bool(True),
    verbose = cms.bool(False),
    readFromCandidates = cms.bool(False),
    jetPtThresholdZ = cms.double(20),
    jetPtThresholdW = cms.double(20),
    bJetCountThreshold = cms.double(0.898),
    useHighestPtHiggsW = cms.bool(True),
    useHighestPtHiggsZ = cms.bool(True),
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

    
  
    

