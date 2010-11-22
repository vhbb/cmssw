import FWCore.ParameterSet.Config as cms

#### V8
### using bcandproducer v1.3
### filtering at minVertices = 0


process = cms.Process("PAT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

#process.MessageLogger.cerr.threshold = "WARNING"
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration/EventContent/EventContent_cff')


process.GlobalTag.globaltag = 'START38_V12::All'

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
    #'/store/data/Run2010A/EG/RECO/v4/000/144/112/FA1804F7-D9B3-DF11-9D71-001D09F26C5C.root'
    '/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/146/713/E2ACD01C-05CA-DF11-8DF1-003048F024DE.root',
    '/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/146/713/38FCC647-07CA-DF11-90F4-001D09F254CE.root',
    '/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/146/713/243890B9-08CA-DF11-9886-001D09F241B9.root',
    '/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/146/713/8A2505BD-0ACA-DF11-A1B5-003048F024FA.root',
    '/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/146/713/22304547-07CA-DF11-86E4-001D09F28F11.root',
    '/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/146/713/5C8A6AB9-08CA-DF11-A6CB-001D09F2437B.root',
    '/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/146/713/5445D44B-0ECA-DF11-B888-001D09F2441B.root',
    '/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/146/713/94303675-FDC9-DF11-A6B8-001D09F2A465.root',
    '/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/146/713/726A87FC-02CA-DF11-AB1E-001617C3B76A.root',
    '/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/146/713/D2585EB4-16CA-DF11-9E6C-0019B9F581C9.root',
    '/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/146/713/5A66CC84-0BCA-DF11-B63E-0030487CD7EA.root',
    '/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/146/713/0A39850E-0ACA-DF11-852F-001D09F23A34.root'
    )) 


## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.trigTools import *

## Output Module Configuration (expects a path 'p')
#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('reducedPAT.root'),
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('filter') ),
                               outputCommands = cms.untracked.vstring('drop *' )
                               )

from PhysicsTools.PatAlgos.tools.coreTools import *
## remove MC matching from the default sequence
#removeMCMatching(process, ['All'])

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )

from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
switchJetCollection(process,cms.InputTag('ak5PFJets'),
                    doJTA        = True,
                    doBTagging   = True,
                    jetCorrLabel = ('AK5PF', cms.vstring(['L2Relative', 'L3Absolute'])),
                    doType1MET   = True,
                    genJetCollection=cms.InputTag("ak5GenJets"),
                    doJetID      = True
                    )
process.selectedPatJets.cut = cms.string("pt>8.")
process.patJets.embedPFCandidates = cms.bool(True)
process.selectedPatJets.embedPFCandidates = cms.bool(True)

###Trigger Matching
#process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
process.load("BAnalysis.BCorrAnalysis.triggermatch_PF_cfi")
from PhysicsTools.PatAlgos.tools.coreTools import removeCleaning
removeCleaning( process )
## Switch on PAT trigger
#from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process ) # This is optional and can be omitted.
switchOnTriggerMatching( process, ['selectedJetTriggerMatchHLTL1Jet6U','selectedJetTriggerMatchHLTL1Jet10U',
                                   'selectedJetTriggerMatchHLTJet15U','selectedJetTriggerMatchHLTJet30U',
                                   'selectedJetTriggerMatchHLTJet50U','selectedJetTriggerMatchHLTJet70U',
                                   'selectedJetTriggerMatchHLTJet100U'] )
# Switch to selected PAT objects in the trigger matching
removeCleaningFromTriggerMatching( process )


### Event Filter
process.goodPFJets = cms.EDFilter("PtMinCandViewSelector",
                                  src = cms.InputTag("selectedPatJets"),
                                  ptMin = cms.double(30)
                                  )

process.EventPFJetFilter = cms.EDFilter("CandViewCountFilter",
                                        src = cms.InputTag("goodPFJets"),
                                        minNumber = cms.uint32(1),
                                        filter = cms.bool(True)   # otherwise it won't filter the events, just produce an empty vertex collection.
                            )

###Sec vertex
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

####### Data cleaning
# require physics declared
process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter = cms.untracked.bool(True),
                                    debugOn = cms.untracked.bool(False),
                                    numtrack = cms.untracked.uint32(10),
                                    thresh = cms.untracked.double(0.25)
                                    )

process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed
#Good Bunch Crossings
process.bptxAnd = hltLevel1GTSeed.clone(L1TechTriggerSeeding = cms.bool(True), L1SeedsLogicalExpression = cms.string('0'))
#BSCNOBEAMHALO
process.bit40 = hltLevel1GTSeed.clone(L1TechTriggerSeeding = cms.bool(True), L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39) AND NOT ((42 AND NOT 43) OR (43 AND NOT 42))'))


#Require a good vertex
process.oneGoodVertexFilter = cms.EDFilter("VertexSelector",
                                           src = cms.InputTag("offlinePrimaryVertices"),
                                           cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
                                           filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
                                           )

process.dump=cms.EDAnalyzer('EventContentAnalyzer')
### rerunning btag (for 35X samples)
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertex3TrkES_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighPurBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighEffBJetTags_cfi")


process.filter = cms.Path(
    ##process.dump*
    ##    process.hltLevel1GTSeed *
    #process.bit40 *
    #process.bptxAnd *
    #process.scrapingVeto *
    ##    process.hltPhysicsDeclared *
    #process.oneGoodVertexFilter*
    #process.HBHENoiseFilter*
    process.simpleSecondaryVertexHighPurBJetTags*
    process.simpleSecondaryVertexHighEffBJetTags*
    process.patDefaultSequence*
    #process.goodPFJets*
    #process.EventPFJetFilter*
    #process.trackFilter*
    process.inclusiveVertexing*process.inclusiveMergedVertices*process.selectedVertices*process.bcandidates
    
    )

#process.schedule = cms.Schedule(process.filter, process.sv, process.pat)
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('reducedPAT.root'),
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('filter') ),
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_selectedPatElectrons*_*_*',
                                                                      'keep *_selectedPatMuons*_*_*',
                                                                      'keep *_selectedPatJets*_*_*',
                                                                      'keep *_patMETs*_*_*',
                                                                      #'keep *_selectedPatPFParticles*_*_*',
                                                                      "keep *_impactParameterTagInfos*_*_*",
                                                                      "keep *_secondaryVertexTagInfos_*_*",
                                                                      # GEN
                                                                      'keep recoGenParticles_genParticles*_*_*',
                                                                      'keep GenEventInfoProduct_*_*_*',
                                                                      'keep GenRunInfoProduct_*_*_*',
                                                                      # RECO
                                                                      'keep recoTracks_generalTracks*_*_*',
                                                                      #'keep *_towerMaker_*_*',
                                                                      'keep *_offlineBeamSpot_*_*',
                                                                      'keep *_offlinePrimaryVertices*_*_*',
                                                                      # TRIGGER
                                                                      'keep edmTriggerResults_TriggerResults*_*_*',
                                                                      'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',
                                                                      'keep *_hltTriggerSummaryAOD_*_*',
                                                                      'keep patTriggerAlgorithms_patTrigger_*_*',
                                                                      'keep patTriggerObjects_patTrigger_*_*',
                                                                      'keep patTriggerFilters_patTrigger_*_*',
                                                                      'keep patTriggerPaths_patTrigger_*_*',
                                                                      'keep patTriggerEvent_patTriggerEvent_*_*',
                                                                      'keep patTriggerObjectStandAlones_patTrigger_*_*',
                                                                      'keep *_patTrigger*_*_*',
                                                                      'keep patTrigger*_patTrigger_*_*',
                                                                      'keep patTriggerObjectStandAlonesedmAssociation_*_*_*',
                                                                      'keep *_*TriggerMatch_*_*',
                                                                      #SV
                                                                      'keep *_bcandidates_*_*',
                                                                      'keep *_selectedVertices*_*_*'
                                                                      )
                               )

process.outpath = cms.EndPath(process.out)

