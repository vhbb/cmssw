import FWCore.ParameterSet.Config as cms

#### V8
### using bcandproducer v1.3
### filtering at minVertices = 0


process = cms.Process("PAT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration/EventContent/EventContent_cff')

#process.GlobalTag.globaltag = 'START3X_V26::All'
#process.GlobalTag.globaltag = 'START36_V9::All'
process.GlobalTag.globaltag = 'START36_V10::All'

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'file:/scratch/leo/QCD_Pt-15_7TeV-pythia8_START36_V10_SP10-v1.root'
    )) 

##ak5GenJets missing...
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

##HLT Filter
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.singleJetHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.singleJetHLTFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","REDIGI36")
process.singleJetHLTFilter.HLTPaths = ["HLT_L1Jet6U", "HLT_L1Jet10U", "HLT_Jet15U"]
process.singleJetHLTFilter.andOr = cms.bool(True) # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true


## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('reducedPATV8b.root'),
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('filter ', 'vertexFilter') ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_selectedPatElectrons*_*_*',
                                                                      'keep *_selectedPatMuons*_*_*',
                                                                      'keep *_selectedPatJets*_*_*',
                                                                      'keep *_patMETs*_*_*',
                                                                      'keep recoGenJets_ak5*_*_*',
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
                                                                      'keep *_hltTriggerSummaryAOD_*_*',
                                                                      'keep patTriggerObjects_patTrigger_*_*',
                                                                      'keep patTriggerFilters_patTrigger_*_*',
                                                                      'keep patTriggerPaths_patTrigger_*_*',
                                                                      'keep patTriggerEvent_patTriggerEvent_*_*',
                                                                      #SV
                                                                      'keep *_bcandidates_*_*',
                                                                      'keep *_selectedVertices_*_*'
                                                                      ) 
                               )

from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *

## uncomment the following line to add tcMET to the event content
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')

## uncomment the following line to add different jet collections
## to the event content
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
                    genJetCollection=cms.InputTag("ak5GenJets"),
                    doJetID      = True,
                    jetIdLabel   = "ak5"
                    )

process.patJets.embedCaloTowers = cms.bool(True)
process.patJetsAK5PF.embedPFCandidates = cms.bool(True)
process.patJetCorrFactorsAK5PF.corrSample  = "Spring10"
process.patJetCorrFactors.corrSample  = "Spring10"
process.selectedPatJets.cut = cms.string("pt>10.")
process.selectedPatJetsAK5PF.cut = cms.string("pt>8.")



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

process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertex3TrkES_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighPurBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighEffBJetTags_cfi")

process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.vertexFilter = cms.Path(
#    process.singleJetHLTFilter+
    process.inclusiveVertexing*process.inclusiveMergedVertices*process.selectedVertices*process.bcandidates
    )

process.filter = cms.Path(
#    process.singleJetHLTFilter+
    process.genJetParticles *
    process.ak5GenJets *
    process.simpleSecondaryVertexHighPurBJetTags*
    process.simpleSecondaryVertexHighEffBJetTags*
    #process.dump
    process.patDefaultSequence
    )

#process.schedule = cms.Schedule(process.filter, process.sv, process.pat)

process.outpath = cms.EndPath(process.out)

