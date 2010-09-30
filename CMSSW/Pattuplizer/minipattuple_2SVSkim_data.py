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


#process.GlobalTag.globaltag = 'GR_R_35X_V8B::All'
#process.GlobalTag.globaltag = 'GR10_P_V5::All'
#process.GlobalTag.globaltag = 'GR_R_37X_V6A::All'
#from 137436 
process.GlobalTag.globaltag = 'GR10_P_V7::All'
#process.GlobalTag.globaltag = 'GR_R_36X_V11A::All'
#process.GlobalTag.globaltag = 'GR_R_36X_V12B::All'
#Jun14th
#process.GlobalTag.globaltag = 'GR_R_36X_V12A::All'

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
    #'/store/data/Run2010A/EG/RECO/v4/000/144/112/FA1804F7-D9B3-DF11-9D71-001D09F26C5C.root'
'/store/data/Commissioning10/MinimumBias/RECO/v9/000/135/735/04AFE369-DD63-DF11-8CFE-0025B3E063F0.root',
'/store/data/Commissioning10/MinimumBias/RECO/v9/000/135/794/B205E0DD-4865-DF11-896F-001A64789E40.root',
'/store/data/Commissioning10/MinimumBias/RECO/v9/000/135/793/3658106B-4B65-DF11-9BB7-00E081B08CC7.root',
'/store/data/Commissioning10/MinimumBias/RECO/v9/000/135/787/786E4A71-4965-DF11-B987-003048D45FD6.root',
'/store/data/Commissioning10/MinimumBias/RECO/v9/000/135/777/D6CBB081-DC63-DF11-A67F-0025B3E05DB6.root',
'/store/data/Commissioning10/MinimumBias/RECO/v9/000/135/771/EC6C2805-DF63-DF11-9D21-00E081791867.root',
'/store/data/Commissioning10/MinimumBias/RECO/v9/000/135/735/98DC5094-D963-DF11-83A2-0025B3E05D5C.root',
'/store/data/Commissioning10/MinimumBias/RECO/v9/000/135/735/8ADB2D7C-DD63-DF11-8761-002481E14F1E.root'
#    'file:/shome/leo/Installations/CMSSW_3_7_0_patch2/src/Analysis/SimpleBAnalysis/test/3CE6FE16-F29E-DF11-8BD3-003048F118C4.root'
    )) 


##HLT Filter
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.singleJetHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.singleJetHLTFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.singleJetHLTFilter.HLTPaths = ["HLT_L1Jet6U", "HLT_L1Jet10U", "HLT_Jet15U"]
process.singleJetHLTFilter.andOr = cms.bool(True) # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true


## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")
### RBX noise filter
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')


## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('reducedPAT.root'),
                               #SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('filter ', 'vertexFilter') ),
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('filter ') ),
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
                                                                      'keep patTriggerObjectStandAlonesedmAssociation_*_*_*',
                                                                      #SV
                                                                      'keep *_bcandidates_*_*',
                                                                      'keep *_selectedVertices_*_*'
                                                                      ) 
                               )

from PhysicsTools.PatAlgos.tools.coreTools import *
## remove MC matching from the default sequence
removeMCMatching(process, ['All'])

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )

from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')

from PhysicsTools.PatAlgos.tools.jetTools import *
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
addJetCollection(process,cms.InputTag('ak5PFJets'),
                    'AK5', 'PF',
                    doJTA        = True,
                    doBTagging   = True,
                    jetCorrLabel = ('AK5', 'PF'),
                    doType1MET   = False,
                    doL1Cleaning = False,
                    doL1Counters = False,
                    genJetCollection=cms.InputTag("ak5PFJets"),
                    doJetID      = True,
                    jetIdLabel   = "ak5"
                    )

process.patJets.embedCaloTowers = cms.bool(True)
process.patJetsAK5PF.embedPFCandidates = cms.bool(True)
process.selectedPatJets.embedCaloTowers = cms.bool(True)
process.selectedPatJetsAK5PF.embedPFCandidates = cms.bool(True)
process.patJetCorrFactorsAK5PF.corrSample  = "Spring10"
process.patJetCorrFactors.corrSample  = "Spring10"
process.selectedPatJets.cut = cms.string("pt>10.")
process.selectedPatJetsAK5PF.cut = cms.string("pt>8.")



### Event Filter
process.goodPFJets = cms.EDFilter("PtMinCandViewSelector",
                                  src = cms.InputTag("selectedPatJetsAK5PF"),
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
                                           cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"),
                                           filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
                                           )

### rerunning btag (for 35X samples)
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertex3TrkES_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighPurBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighEffBJetTags_cfi")

process.dump=cms.EDAnalyzer('EventContentAnalyzer')

#process.vertexFilter = cms.Path(
#    #process.hltLevel1GTSeed *
#    process.bit40 *
#    process.bptxAnd *
#    process.scrapingVeto *
#    #process.hltPhysicsDeclared *
#    process.oneGoodVertexFilter* 
#    #process.singleJetHLTFilter+
#    process.HBHENoiseFilter*
#    process.goodPFJets*
#    process.EventPFJetFilter*
#    process.inclusiveVertexing*process.inclusiveMergedVertices*process.selectedVertices*process.bcandidates
#)

process.filter = cms.Path(
    #    process.hltLevel1GTSeed *
    process.bit40 *
    process.bptxAnd *
    process.scrapingVeto *
    #    process.hltPhysicsDeclared *
    process.oneGoodVertexFilter*
    #process.singleJetHLTFilter+
    process.HBHENoiseFilter*
    process.simpleSecondaryVertexHighPurBJetTags*
    process.simpleSecondaryVertexHighEffBJetTags*
    #process.dump
    process.patDefaultSequence*
    process.goodPFJets*
    process.EventPFJetFilter*
    process.inclusiveVertexing*process.inclusiveMergedVertices*process.selectedVertices*process.bcandidates
    
    )

#process.schedule = cms.Schedule(process.filter, process.sv, process.pat)

process.outpath = cms.EndPath(process.out)

