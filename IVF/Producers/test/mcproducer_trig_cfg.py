import FWCore.ParameterSet.Config as cms

process = cms.Process("QCD")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.Tracer = cms.Service("Tracer")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    #skipEvents = cms.untracked.uint32(2867),                
    fileNames = cms.untracked.vstring(
    #'file:/shome/wehrlilu/SecVtx/SVProd/Zbb/CMSSW_3_1_2/src/step2_RAW2DIGI_RECO_100ev.root'
    'file:/shome/wehrlilu/SecVtx/SVProd/Zbb/F605AE2E-D990-DE11-BFCB-001D09F2545B.root'
    ),
    #duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

)

##############
process.bfilter = cms.EDFilter("PythiaFilter",
    ParticleID = cms.untracked.int32(5)
)
##############

##############
process.BHadrons = cms.EDProducer('MCBHadronProducer',
    jetsrc  = cms.InputTag('iterativeCone5CaloJets')
)
##############

##############
process.SimSecVertices = cms.EDProducer('MCSecVertexProducer',
    jetsrc  = cms.InputTag('iterativeCone5CaloJets')
    #jetsrc  = cms.InputTag('L2L3CorJetIC5Calo')
)
##############

##############
process.prunedGenParticles = cms.EDProducer(
    'GenParticlePruner',
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    'keep abs(pdgId)/100 == 5','keep abs(pdgId)/1000 == 5','keep abs(pdgId) == 5','keep status == 3','keep status == 2 & pt > 0.5'
    )
)
##############

process.out = cms.OutputModule("PoolOutputModule",
    #fileName = cms.untracked.string('trigtest.root'),
    fileName = cms.untracked.string('QCDpt80_trig.root'),
    #fileName = cms.untracked.string('qcdPt300_test.root'),
    SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('p')
    ),	    
    outputCommands = cms.untracked.vstring(
'drop *',
'keep *_hltTriggerSummaryAOD_*_HLT*',
'keep *_TriggerResults_*_HLT*',
'keep *_*_*_QCD',
'keep recoTracks_generalTracks_*_RECO',
'keep recoMuons_*_*_*',
#'keep recoTrackExtras_*_*_*',
'keep *_secondaryVertexTagInfos_*_*',
'keep recoVertexs_*_*_RECO',
'keep recoCaloJets_*_*_RECO',                          
'keep recoPFJets_*_*_RECO',                          
'keep recoGenParticles_genParticles_*_*',
'keep recoPFCandidates_particleFlow_*_RECO',
'keep CaloTowersSorted_towerMaker_*_RECO' 
))

##############
process.bfilter = cms.EDFilter("PythiaFilter",
    ParticleID = cms.untracked.int32(5)
)
##############


##############
#impact parameter, jettrackassociation
process.load("RecoBTag.ImpactParameter.impactParameter_cfi")
process.load("RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi")
process.ic5JetTracksAssociatorAtVertex.coneSize = 1.0
#GONE IN 225: process.impactParameterTagInfos.maximumDistanceToJetAxis = 9999.0
#GONE IN 225: process.impactParameterTagInfos.maximumDecayLength = 9999.0
process.impactParameterTagInfos.jetTracks = cms.InputTag("ic5JetTracksAssociatorAtVertex","","QCD")
##############

##############
#secondary vertices
#C10
from RecoBTag.SecondaryVertex.secondaryVertexTagInfos_cfi import *
#secondaryVertexTagInfos.trackSelection.jetDeltaRMax = 1.0
#secondaryVertexTagInfos.trackSelection.maxDistToAxis = 99999.0
#secondaryVertexTagInfos.vertexCuts.maxDeltaRToJetAxis = 99999.0
#C08
process.secondaryVertexTagInfosC08 = secondaryVertexTagInfos.clone()
process.secondaryVertexTagInfosC08.trackSelection.jetDeltaRMax = 0.8
#process.secondaryVertexTagInfosC08.trackSelection.maxDistToAxis = 99999.0
process.secondaryVertexTagInfosC08.vertexCuts.maxDeltaRToJetAxis = 0.8
#C03
#process.secondaryVertexTagInfosC03 = secondaryVertexTagInfos.clone()
#process.secondaryVertexTagInfosC03.trackSelection.jetDeltaRMax = 0.3
#process.secondaryVertexTagInfosC03.trackSelection.maxDistToAxis = 99999.0
#process.secondaryVertexTagInfosC03.vertexCuts.maxDeltaRToJetAxis = 99999.0

#NEW
#process.secondaryVertexTagInfosC03.vertexCuts.distVal2dMin = -99999
#process.secondaryVertexTagInfosC03.vertexCuts.distSig2dMin = -99999
#process.secondaryVertexTagInfosC08.vertexCuts.distVal2dMin = -99999
#process.secondaryVertexTagInfosC08.vertexCuts.distSig2dMin = -99999
#secondaryVertexTagInfos.vertexCuts.distVal2dMin = -99999
#secondaryVertexTagInfos.vertexCuts.distSig2dMin = -99999###

#process.secondaryVertexTagInfosC03.vertexCuts.fracPV = 2
#process.secondaryVertexTagInfosC08.vertexCuts.fracPV = 2
#secondaryVertexTagInfos.vertexCuts.fracPV = 2
#process.secondaryVertexTagInfosC03.vertexCuts.massMax = 99999
#process.secondaryVertexTagInfosC08.vertexCuts.massMax = 99999
#secondaryVertexTagInfos.vertexCuts.massMax = 99999

#process.secondaryVertexTagInfosC08.useBeamConstraint = False
#process.secondaryVertexTagInfosC03.useBeamConstraint = False
#secondaryVertexTagInfos.useBeamConstraint = False

##############

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "IDEAL_V9::All"
process.GlobalTag.globaltag = "MC_31X_V3::All"

#process.p = cms.Path(process.ic5JetTracksAssociatorAtVertex*process.impactParameterTagInfos*process.secondaryVertexTagInfosC08*process.secondaryVertexTagInfosC03*secondaryVertexTagInfos)
#process.p = cms.Path(process.ic5JetTracksAssociatorAtVertex*process.impactParameterTagInfos*secondaryVertexTagInfos)

##with prunedGenParticles
#process.p = cms.Path(process.prunedGenParticles*process.SimSecVertices*process.BHadrons*process.ic5JetTracksAssociatorAtVertex*process.impactParameterTagInfos*process.secondaryVertexTagInfosC08)

##without prunedGenParticles
process.p = cms.Path(process.SimSecVertices*process.BHadrons*process.ic5JetTracksAssociatorAtVertex*process.impactParameterTagInfos*process.secondaryVertexTagInfosC08)

process.e = cms.EndPath(process.out)
