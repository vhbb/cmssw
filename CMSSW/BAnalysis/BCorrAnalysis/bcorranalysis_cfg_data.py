import FWCore.ParameterSet.Config as cms

process = cms.Process("BCA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "ERROR"
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR10_P_V6::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'/store/user/leo/JetMETTau/BJetsPatDumpV8-CMSSW_3_7_0_patch2-Run2010A-PromptReco-v4-Runs_139460-139965/0409b26ad4f90c0ce2038da6d13e58f2/reducedPATV8_JetMETTau_Run2010A-PromptReco-v4-Runs_139460-139965_7_1_LHW.root'
    #'/store/user/leo/QCD_Pt80-herwig/BJetsPatDumpV8b-CMSSW_3_7_0_patch2-Spring10-START3X_V26-v1/fe983afb43c677d8daa662736af06b64/reducedPATV8b-QCD_Pt80-herwig_Spring10-START3X_V26_S09-v1_17_1_Emw.root'
    #'/store/user/leo/JetMETTau/BJetsPatDumpV8-CMSSW_3_7_0_patch2-Run2010A-PromptReco-v4-Runs_139460-139965/0409b26ad4f90c0ce2038da6d13e58f2/reducedPATV8_JetMETTau_Run2010A-PromptReco-v4-Runs_139460-139965_7_1_LHW.root'
    #'/store/user/leo/JetMET/BJetsPatDumpV8-CMSSW_3_7_0_patch2-Run2010A-PromptReco-v4_Runs141950-143731-2/95f977e9e64c4e6df50e2f0afa9a2a31/reducedPAT_30_2_KE8.root'
    'file:/shome/leo/Installations/CMSSW_3_7_0_patch2/src/BAnalysis/BCorrAnalysis/reducedPAT.root'
)
                            
)


process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDBMC36X")
process.load ("RecoBTag.PerformanceDB.BTagPerformanceDBMC36X")

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
    )


process.TFileService = cms.Service("TFileService",
                                  fileName = cms.string("miniTreePAT.root")
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
                                    GenJetCollection      = cms.untracked.InputTag("ak5GenJets"),
                                    PFJetCollection      = cms.untracked.InputTag("selectedPatJetsAK5PF"),

                                    BCorrMethod = cms.untracked.string("MC"),
                                    
                                    isData       = cms.untracked.int32(1),
                                    CaloJetSelection_minPt = cms.untracked.double(10),
                                    CaloJetSelection_maxEta = cms.untracked.double(5),
                                    CaloJetSelection_EMF = cms.untracked.double(0.01),
                                    CaloJetSelection_fHPD = cms.untracked.double(0.98),
                                    CaloJetSelection_n90Hits = cms.untracked.double(1),
                                    
                                    PFJetSelection_minPt = cms.untracked.double(8),
                                    PFJetSelection_maxEta= cms.untracked.double(5),
                                    PFJetSelection_NHF= cms.untracked.double(1),
                                    PFJetSelection_NEF= cms.untracked.double(1),
                                    PFJetSelection_CEF= cms.untracked.double(1),
                                    PFJetSelection_CHF= cms.untracked.double(0.0),
                                    PFJetSelection_CM= cms.untracked.double(0.0),
                                    PFJetSelection_NCONST= cms.untracked.double(1),
                                    
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
                                                     )
                                    
                                    )


###for data
#process.GlobalTag.globaltag = cms.string('GR10_P_V6::All')
###for MC
#process.GlobalTag.globaltag = cms.string('START3X_V26B::All')

# firing trigger objects used in succeeding HLT path 'HLT_Jet15U'
process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
#process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cff" )
#from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi import *

process.load("BAnalysis.BCorrAnalysis.triggermatch_PF_cfi")
process.load("BAnalysis.BCorrAnalysis.triggermatch_CALO_cfi")

process.patTriggerMatcher += process.selectedJetTriggerMatchHLTL1Jet10U
process.patTriggerMatcher += process.selectedJetTriggerMatchHLTJet15U
process.patTriggerMatcher += process.selectedJetTriggerMatchHLTJet30U
process.patTriggerMatcher += process.selectedJetTriggerMatchHLTJet50U
process.patTriggerMatcher += process.selectedJetTriggerMatchHLTJet70U
process.patTriggerMatcher += process.selectedJetTriggerMatchHLTJet100U
process.patTriggerMatcher += process.selectedJetTriggerMatchHLTL1Jet10UCALO
process.patTriggerMatcher += process.selectedJetTriggerMatchHLTJet15UCALO
process.patTriggerMatcher += process.selectedJetTriggerMatchHLTJet30UCALO
process.patTriggerMatcher += process.selectedJetTriggerMatchHLTJet50UCALO
process.patTriggerMatcher += process.selectedJetTriggerMatchHLTJet70UCALO
process.patTriggerMatcher += process.selectedJetTriggerMatchHLTJet100UCALO


process.patTriggerEvent.patTriggerMatches = [ "selectedJetTriggerMatchHLTL1Jet10U","selectedJetTriggerMatchHLTJet15U","selectedJetTriggerMatchHLTJet30U",
                                              "selectedJetTriggerMatchHLTJet50U","selectedJetTriggerMatchHLTJet70U","selectedJetTriggerMatchHLTJet100U",
                                              "selectedJetTriggerMatchHLTL1Jet10UCALO","selectedJetTriggerMatchHLTJet15UCALO","selectedJetTriggerMatchHLTJet30UCALO",
                                              "selectedJetTriggerMatchHLTJet50UCALO","selectedJetTriggerMatchHLTJet70UCALO","selectedJetTriggerMatchHLTJet100UCALO"
                                              ]


process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.p = cms.Path(  process.matchPATPF*
                       process.matchPATCALO*
                       #process.selectedJetTriggerMatchHLTJet15U*
                       #process.selectedJetTriggerMatchHLTJet30U*
                        patTriggerEvent*
                       process.bcandidates *
                       #process.dump
                       process.bcanalyzer
                       )


