import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

process = cms.Process("FWLitePlots")

#fileNames   = cms.vstring('file:2l2bMetEdmNtuples.root'),         ## mandatory
process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_10_0_0dz.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_11_0_syf.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_12_0_98x.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_1_2_fZC.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_13_0_d07.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_14_0_gTV.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_15_0_k8w.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_16_0_Gyt.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_17_0_Sbb.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_18_0_08B.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_19_0_FPm.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_20_0_K9m.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_2_0_o7r.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_21_0_KC0.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_22_0_fym.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_23_0_PdH.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_24_0_gFX.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_25_0_Jta.root",
"file:PAT.edm_26_0_pL3.root",
#"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_26_0_pL3.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_27_0_2RN.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_28_0_ILe.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_30_0_lyu.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_3_0_JbI.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_31_0_l4y.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_32_0_jyU.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_33_0_CwC.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_34_0_aWH.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_35_0_6hS.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_36_0_3cS.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_37_0_akJ.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_38_0_YeV.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_39_1_Yve.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_40_0_ifl.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_4_0_q1y.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_41_0_hHm.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_5_0_hty.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_6_0_5vG.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_7_0_9A5.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_8_0_f1B.root",
"dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store/user/arizzi/TTJets_TuneZ2_7TeV-madgraph-tauola/HBB_EDMNtupleV3_ProcV1/f0063de7ad3f08dda539b5259b13f60a//PAT.edm_9_1_juk.root"

),
    PUmcfileName = cms.string(""),
    PUdatafileName = cms.string(""),
    maxEvents   = cms.int32(-1),                             ## optional
    outputEvery = cms.uint32(0),                            ## optional
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    )

# get JSON file 
JSONfile = 'Input.json'
lumiList = LumiList.LumiList (filename = JSONfile).getCMSSWString().split(',')

#Uncomment to run with JSON
#process.fwliteInput.lumisToProcess.extend(lumiList)


channel =  "TTbar"
import os
dirnameOld = "/pnfs/pi.infn.it/data/cms/store/user/bortigno/DoubleMu/HBB_EDMNtupleV3_ProcV1_may/07fb60889166b64f474d8d0aa162db69/"




#for i in range(len(channels)):
 

#dirname =  dirnameOld 
#dirlist = os.listdir(dirname)
#basenamelist = os.listdir(dirname + "/")
#for basename in basenamelist:
#   process.fwliteInput.fileNames.append("dcap://cmsdcache" + dirname + "/" + basename)
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
	"HLT_Ele27_CaloIdVT_CaloIsoT_TrkId_TrkIsoT_v.*", #3
	"HLT_Ele27_WP80_PFMHT50_v.*", #4
        "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v.*", #5
        "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v.*", #6
        "HLT_DiCentralJet20_BTagIP_MET65_v.*", #7
	"HLT_MET120_v.*" #8
	"HLT_CentralJet80_MET80_v.*", #9
	"HLT_PFMHT150_v.*", #10
	"HLT_DiCentralJet20_MET80_v.*", #11
   ),
    isMC =     cms.bool(False),
    verbose = cms.bool(False),
    readFromCandidates = cms.bool(True),
    jetPtThresholdZ = cms.double(20),
    jetPtThresholdW = cms.double(30),
    bJetCountThreshold = cms.double(0.898),
    useHighestPtHiggsW = cms.bool(True),
    useHighestPtHiggsZ = cms.bool(False),
    idMuFileName = cms.string("ScaleEffs42.root"),
    hltMuFileName = cms.string("IsoToHLT42.root"),

    )

    
  
    

