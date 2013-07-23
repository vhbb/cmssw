import FWCore.ParameterSet.Config as cms

process = cms.Process("BCA")

from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.trigTools import *

process.setName_("BCA")


#process = cms.Process("BCA")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = "ERROR"
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START38_V13::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'/store/user/leo/JetMETTau/BJetsPatDumpV8-CMSSW_3_7_0_patch2-Run2010A-PromptReco-v4-Runs_139460-139965/0409b26ad4f90c0ce2038da6d13e58f2/reducedPATV8_JetMETTau_Run2010A-PromptReco-v4-Runs_139460-139965_7_1_LHW.root'
    #'/store/user/leo/QCD_Pt80-herwig/BJetsPatDumpV8b-CMSSW_3_7_0_patch2-Spring10-START3X_V26-v1/fe983afb43c677d8daa662736af06b64/reducedPATV8b-QCD_Pt80-herwig_Spring10-START3X_V26_S09-v1_17_1_Emw.root'
    #'/store/user/leo/JetMETTau/BJetsPatDumpV8-CMSSW_3_7_0_patch2-Run2010A-PromptReco-v4-Runs_139460-139965/0409b26ad4f90c0ce2038da6d13e58f2/reducedPATV8_JetMETTau_Run2010A-PromptReco-v4-Runs_139460-139965_7_1_LHW.root'
    #'/store/user/leo/JetMET/BJetsPatDumpV8-CMSSW_3_7_0_patch2-Run2010A-PromptReco-v4_Runs141950-143731-2/95f977e9e64c4e6df50e2f0afa9a2a31/reducedPAT_30_2_KE8.root'
    'file:/shome/leo/Installations/CMSSW_3_8_6/src/Pattuplizer/reducedPAT_12_1_HKv.root'
)
                            
)


#process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1011")
#process.load ("RecoBTag.PerformanceDB.BTagPerformanceDB1011")
process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDBMC36X")
process.load ("RecoBTag.PerformanceDB.BTagPerformanceDBMC36X")

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#to be removed
process.ak5PFL2Relative.useCondDB = False
process.ak5PFL3Absolute.useCondDB = False
process.ak5PFResidual.useCondDB = False
process.ak5CaloL2Relative.useCondDB = False
process.ak5CaloL3Absolute.useCondDB = False
process.ak5CaloResidual.useCondDB = False


process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
    )


process.TFileService = cms.Service("TFileService",
                                  fileName = cms.string("miniTreePAT.root")
)

process.bhadrons = cms.EDProducer('MCBHadronProducer',
                                  quarkId = cms.uint32(5)
                                  )

process.bcandidates = cms.EDProducer('BCandidateProducer',
                                             src = cms.InputTag('selectedVertices','',''),
                                             primaryVertices = cms.InputTag('offlinePrimaryVerticesWithBS','',''),
                                             minDRUnique = cms.untracked.double(0.4),
                                             minvecSumIMifsmallDRUnique = cms.untracked.double(5.5), #5.5
                                             minCosPAtomerge = cms.untracked.double(0.99),
                                             maxPtreltomerge = cms.untracked.double(7777.0)


                                     )


process.bcanalyzer = cms.EDAnalyzer('BCorrAnalyzer',
                                    src = cms.InputTag('bcandidates','',''),
                                    simbsrc = cms.InputTag('bhadrons','',''),
                                    primaryVertices = cms.InputTag('offlinePrimaryVerticesWithBS','',''),
                                    trigressrc = cms.InputTag('TriggerResults','','HLT'),
                                    minInvM = cms.untracked.double(1.4),
                                    mindistSig3D = cms.untracked.double(5.0),
                                    maxEtaVertex = cms.untracked.double(2.0),
                                    
                                    #Jet analysis
                                    JetCollection      = cms.untracked.InputTag("selectedPatJets"),
                                    GenJetCollection      = cms.untracked.InputTag("selectedPatJets:genJets"),
                                    PFJetCollection      = cms.untracked.InputTag("selectedPatJets"),

                                    BCorrMethod = cms.untracked.string("MC"),
                                    
                                    isData       = cms.untracked.int32(0),

                                    JEC_PATH = cms.untracked.string(""),
                                    JEC_RES_FILE = cms.untracked.string("START38_V13_AK5PF_L2L3Residual.txt"),
                                    JEC_UNC_FILE = cms.untracked.string("START38_V13_AK5PF_Uncertainty.txt"),
                                    
                                    PFJetSelection_minPt = cms.untracked.double(8),
                                    PFJetSelection_maxEta= cms.untracked.double(5),
                                    
                                    BTag_tchel = cms.untracked.double(1.9),
                                    BTag_tchem = cms.untracked.double(3.99),
                                    BTag_tchpm = cms.untracked.double(2.17),
                                    BTag_tchpt = cms.untracked.double(4.31),
                                    BTag_ssvhem = cms.untracked.double(2.02),
                                    BTag_ssvhet = cms.untracked.double(3.4),
                                    BTag_ssvhpt = cms.untracked.double(2.),
                                    
                                    jetID = cms.PSet(useRecHits = cms.bool(True),
                                                     hbheRecHitsColl = cms.InputTag("hbhereco"),
                                                     hoRecHitsColl   = cms.InputTag("horeco"),
                                                     hfRecHitsColl   = cms.InputTag("hfreco"),
                                                     ebRecHitsColl   = cms.InputTag("ecalRecHit", "EcalRecHitsEB"),
                                                     eeRecHitsColl   = cms.InputTag("ecalRecHit", "EcalRecHitsEE")
                                                     ),

                                    impactParameterTagInfos = cms.untracked.InputTag("impactParameterTagInfosAOD")
                                    
                                    )



process.dump=cms.EDAnalyzer('EventContentAnalyzer')

##Vertex producer
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexing_cff')
process.inclusiveMergedVertices = process.vertexMerger.clone()
process.inclusiveMergedVertices.secondaryVertices = cms.InputTag("inclusiveVertices")
process.inclusiveMergedVertices.maxFraction = 0.2
process.inclusiveMergedVertices.minSignificance = 10.

process.load("RecoBTag/SecondaryVertex/bVertexFilter_cfi")
process.selectedVertices = process.bVertexFilter.clone()
process.selectedVertices.secondaryVertices = cms.InputTag("inclusiveMergedVertices")
process.selectedVertices.minVertices = 0
process.selectedVertices.vertexFilter.multiplicityMin = 3

process.inclusiveVertexFinder.clusterScale = 1.
process.inclusiveVertexFinder.clusterMinAngleCosine = 0.5

process.bcandidates = cms.EDProducer('BCandidateProducer',
                                                                          src = cms.InputTag('selectedVertices','',''),
                                                                          primaryVertices =
                                                                          cms.InputTag('offlinePrimaryVerticesWithBS','',''),
                                                                          minDRUnique = cms.untracked.double(0.4),
                                                                          minvecSumIMifsmallDRUnique = cms.untracked.double(5.5),
                                                                          minCosPAtomerge = cms.untracked.double(0.99),
                                                                          maxPtreltomerge = cms.untracked.double(7777.0)
                                                                          )




process.p = cms.Path(  #process.matchPATPF*
                       #process.matchPATCALO*
                       #process.selectedJetTriggerMatchHLTJet15U*
                       #process.selectedJetTriggerMatchHLTJet30U*
                       # patTriggerEvent*
    process.inclusiveVertexing
    *process.inclusiveMergedVertices
    *process.selectedVertices*process.bhadrons *process.bcandidates *
    #process.dump
    process.bcanalyzer
    )


