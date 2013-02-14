import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

process = cms.Process("FWLitePlots")

#fileNames   = cms.vstring('file:2l2bMetEdmNtuples.root'),         ## mandatory
process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(

#"file:/gpfs/ddn/srm/cms/store/user/lpchbb/dlopes/WH125_Summer12_53X_V42b/dlopes/WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp/HBB_EDMNtupleV42/9803889241b1fc304f795d3b3875632d//PAT.edm_10_1_lT1.root"
#"file:/gpfs/ddn/srm/cms/store/user/lpchbb/degrutto/TT_CT10_TuneZ2star_8TeV-powheg-tauolaSummer12_DR53X-PU_S10_START53_V7A-v2/degrutto/TT_CT10_TuneZ2star_8TeV-powheg-tauola/HBB_EDMNtupleV42/9803889241b1fc304f795d3b3875632d/PAT.edm_693_1_Sys.root"      
###"dcache:///pnfs/cms/WAX/11/store/user/lpchbb/msegala/WH_WToLNu_HToBB_M-120_8TeV-powheg-herwigpp/HBB_EDMNtupleV30_May2012Prod_2/e63abc9239312b3f728ddbd5ef32b303/"dcap:///pnfs/cms/WAX/11/store/user/degrutto//testMET2012/PAT.edm_11_1_oQQ.root" ,
##    "../../HbbAnalyzer/test/"dcap:///pnfs/cms/WAX/11/store/user/degrutto//testMET2012/PAT.edm.root",
#"dcache:///cmsdcache//pnfs/pi.infn.it/data/cms/store/user/tboccali/WH_WToLNu_HToBB_M-120_8TeV-powheg-herwigpp/HBB_EDMNtupleV30_ProcV1_WH_WToLNu_HToBB_M-120/14fe2b624ddea84f5c39709f51bf546f/"dcap:///pnfs/cms/WAX/11/store/user/degrutto//testMET2012/PAT.edm_51_1_3LJ.root"
#/pnfs/pi.infn.it/data/cms/store/user/tboccali/ZH_ZToNuNu_HToBB_M-120_8TeV-powheg-herwigpp/HBB_EDMNtupleV30_ProcV1_ZH_ZToNuNu_HToBB_M-120/14fe2b624ddea84f5c39709f51bf546f/"dcap:///pnfs/cms/WAX/11/store/user/degrutto//testMET2012/PAT.edm_51_1_ukR.root
"dcache:///pnfs/cms/WAX/11/store/user/lpchbb/msegala/WH_WToLNu_HToBB_M-125_TAU_ID_TEST_PromptV1EdmV42/PAT.edm_1_1_xCI.root"

),

    PUmcfileName2011B= cms.string("Summer12MCObserved.root"),
    PUdatafileName2011B = cms.string("MyDataPileupHistogramObserved.root"),
    PUmcfileName = cms.string("MC_S10_fromTwiki_60bins.root"),
    PUdatafileNameAB = cms.string("data_PU_60bins_190456-196509_69.4mb.root"),
    PUdatafileNameMinus = cms.string("data_PU_60bins_190456-208686_66.5mb.root"),
    PUdatafileNamePlus = cms.string("data_PU_60bins_190456-208686_72.4mb.root"),
    PUdatafileName = cms.string("data_PU_60bins_190456-208686_69.4mb.root"),
    badEventsFileName = cms.string("HCALLaser2012AllDatasets.txt"),
    Weight3DfileName = cms.string(""),
    maxEvents   = cms.int32(-1),                             ## optional
    runMin  = cms.int32(-1),
    runMax  = cms.int32(-1),
    skipEvents   = cms.int32(0),                             ## optional
    outputEvery = cms.uint32(0),                            ## optional
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    )

# get JSON file 
JSONfile = 'Cert_190456-203002_8TeV_PromptReco_Collisions12_JSON.txt'
lumiList = LumiList.LumiList (filename = JSONfile).getCMSSWString().split(',')

#Uncomment to run with JSON
process.fwliteInput.lumisToProcess.extend(lumiList)



channel =  "TEST"
import os
#dirnameOld = "//pnfs/pi.infn.it/data/cms/store/user/bortigno/DoubleMu/HBB_EDMNtupleV3_ProcV1_may/07fb60889166b64f474d8d0aa162db69/"
#dirnameOld = "//pnfs/cms/WAX/11/store/user/lpchbb/degrutto/METRun2011APromptV1EdmV31"
##dirnameOld = "//pnfs/cms/WAX/11/store/user/lpchbb/degrutto/METRun2011APromptV1EdmV31/degrutto/MET/HBB_EDMNtupleV31_May2012Prod/f3f16f0a8d4ba1fd82b25c3d99de1c78/"
#dirnameOld = "//pnfs/cms/WAX/11/store/user/lpchbb/degrutto/METRun2011BPromptV1EdmV31/degrutto/MET/HBB_EDMNtupleV31_May2012Prod/f3f16f0a8d4ba1fd82b25c3d99de1c78/"
#dirnameOld = "//pnfs/cms/WAX/11/store/user/lpchbb/degrutto/SingleMuRun2012APromptV1EdmV31/degrutto/SingleMu/HBB_EDMNtupleV31_May2012Prod/f3f16f0a8d4ba1fd82b25c3d99de1c78/"
#dirnameOld = "//pnfs/cms/WAX/11/store/user/lpchbb/degrutto/SingleMuRun2012BPromptV1EdmV31/degrutto/SingleMu/HBB_EDMNtupleV31_May2012Prod/f3f16f0a8d4ba1fd82b25c3d99de1c78/"
#dirnameOld = "//pnfs/cms/WAX/11/store/user/lpchbb/dlopes/SingleMuRun2012APromptV1EdmV33/dlopes/SingleMu/HBB_EDMNtupleV33/fc9898cbaf0c3f28aabf20b09408f841/"
#dirnameOld = "//pnfs/cms/WAX/11/store/user/lpchbb/dlopes/SingleMuRun2012BPromptV1EdmV33/dlopes/SingleMu/HBB_EDMNtupleV33/fc9898cbaf0c3f28aabf20b09408f841/"
##dirnameOld = "//pnfs/cms/WAX/11/store/user/lpchbb/degrutto/METRun2012APromptJul12/"
#dirnameOld = "//pnfs/cms/WAX/11/store/user/lpchbb/degrutto/METRun2012BPromptJul12/"
dirnameOld = "/"

#for i in range(len(channels)):
 

dirname =  dirnameOld 
dirlist = os.listdir(dirname)
basenamelist = os.listdir(dirname + "/")
for basename in basenamelist:
#   process.fwliteInput.fileNames.append("dcache:/" + dirname + "/" + basename)
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
        "HLT_LooseIsoPFTau35_Trk20_MET75_v.*", #38
        "HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v.*", #39
        "HLT_DiCentralJet20_CaloMET65_BTagCSV07_PFMHT80_v.*", #40
        "HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v.*", #41
        "HLT_PFMET150_v.*", #42
        "HLT_L1ETM40_v.*", #43
        "HLT_Ele27_WP80_v.*", #44
        "HLT_Ele27_WP80_WCandPt80_v.*", #45
        "HLT_IsoMu20_eta2p1_WCandPt80_v.*", #46
        "HLT_IsoMu20_WCandPt80_v.*", #47
        "HLT_Mu17_TkMu8_v.*", #48
        "HLT_DiCentralPFJet30_PFMHT80_v.*", #49 ## run2012A
        "HLT_DiCentralPFJet30_PFMET80_v.*", #50 ## run2012B prescaled       
        "HLT_LooseIsoPFTau35_Trk20_Prong1_v.*", #51
        "HLT_LooseIsoPFTau35_Trk20_Prong1_MET70_v.*", #52
        "HLT_LooseIsoPFTau35_Trk20_Prong1_MET75_v.*", #53
       
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
    hltSingleEle2012Awp95 = cms.string("triggerRootFiles/SingleEle.TrigEff.wp95.2012ABC.root"),
    hltSingleEle2012Awp80 = cms.string("triggerRootFiles/SingleEle.TrigEff.wp80.2012ABC.root"),
    hltSingleMuon2012A = cms.string("triggerRootFiles/SingleMu24OR40.TrigEff.2012ABC.root"),
    hltDoubleEle2012A_leg8 = cms.string("triggerRootFiles/DoubleEle8.TrigEff.wp95.2012ABC.root"),
    hltDoubleEle2012A_leg17 = cms.string("triggerRootFiles/DoubleEle17.TrigEff.wp95.2012ABC.root"),
    hltDoubleMuon2012A_leg8 = cms.string("triggerRootFiles/DoubleMu8.TrigEff.2012ABC.root"),
    hltDoubleMuon2012A_leg17 = cms.string("triggerRootFiles/DoubleMu17.TrigEff.2012ABC.root"),
    hltMuPlusWCandPt2012A_legMu = cms.string("triggerRootFiles/SingleMu20Not24Or40.TrigEff.2012ABC.root"),
    hltMuPlusWCandPt2012A_legW = cms.string("triggerRootFiles/WCandPt.TrigEff.2012ABC.root"),
    hltDoubleMuon2012A_dZ = cms.string("triggerRootFiles/DoubleMuDz.TrigEff.2012ABC.root"),
    hltDoubleEle2012A_dZ = cms.string("triggerRootFiles/DoubleEleDz.TrigEff.2012ABC.root"),
    idMu2012A = cms.string("triggerRootFiles/MuRecoId.ScaleFactor.2012ABC.root"),
    idEle2012A = cms.string("triggerRootFiles/EleRecoId.ScaleFactor.wp95.2012ABC.root"),
    idEle2012Awp80 = cms.string("triggerRootFiles/EleRecoId.ScaleFactor.wp80.2012ABC.root"),

    jecFolder = cms.string("jec"),
    csvDiscr = cms.string("csvdiscr.root"),
    btagEffFileName = cms.string("btag_generic.txt")
    )

    
  
    
process.Analyzer2012ABOnly = cms.PSet(
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
    hltSingleEle2012Awp95 = cms.string("triggerRootFiles/SingleEle.TrigEff.wp95.2012AB.root"),
    hltSingleEle2012Awp80 = cms.string("triggerRootFiles/SingleEle.TrigEff.wp80.2012AB.root"),
    hltSingleMuon2012A = cms.string("triggerRootFiles/SingleMu24OR40.TrigEff.2012AB.root"),
    hltDoubleEle2012A_leg8 = cms.string("triggerRootFiles/DoubleEle8.TrigEff.wp95.2012AB.root"),
    hltDoubleEle2012A_leg17 = cms.string("triggerRootFiles/DoubleEle17.TrigEff.wp95.2012AB.root"),
    hltDoubleMuon2012A_leg8 = cms.string("triggerRootFiles/DoubleMu8.TrigEff.2012AB.root"),
    hltDoubleMuon2012A_leg17 = cms.string("triggerRootFiles/DoubleMu17.TrigEff.2012AB.root"),
    hltMuPlusWCandPt2012A_legMu = cms.string("triggerRootFiles/SingleMu20Not24Or40.TrigEff.2012AB.root"),
    hltMuPlusWCandPt2012A_legW = cms.string("triggerRootFiles/WCandPt.TrigEff.2012AB.root"),
    hltDoubleMuon2012A_dZ = cms.string("triggerRootFiles/DoubleMuDz.TrigEff.2012AB.root"),
    hltDoubleEle2012A_dZ = cms.string("triggerRootFiles/DoubleEleDz.TrigEff.2012AB.root"),
    idMu2012A = cms.string("triggerRootFiles/MuRecoId.ScaleFactor.2012AB.root"),
    idEle2012A = cms.string("triggerRootFiles/EleRecoId.ScaleFactor.wp95.2012AB.root"),
    idEle2012Awp80 = cms.string("triggerRootFiles/EleRecoId.ScaleFactor.wp80.2012AB.root"),
)


