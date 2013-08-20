import FWCore.ParameterSet.Config as cms

process = cms.Process("MCJP")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:/shome/wehrlilu/SecVtx/SVProd/twoBAllStatFiltered.root',
    )
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('BandTrackjets.root'),
    SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('p')
    ),	    
    outputCommands = cms.untracked.vstring(
'drop *',
'keep *_*_*_MCP',
'keep *_*_*_MCJP',
'keep recoTracks_generalTracks_*_RECO',
'keep recoMuons_*_*_*',
#'keep recoTrackExtras_*_*_*',
'keep *_secondaryVertexTagInfos_*_*',
'keep recoVertexs_*_*_RECO',
'keep recoCaloJets_*_*_RECO',                          
'keep CaloTowersSorted_towerMaker_*_RECO') 
)

##############
process.bfilter = cms.EDFilter("PythiaFilter",
    ParticleID = cms.untracked.int32(5)
)
##############



##############
#TRACK JETS

#from Configuration.StandardSequences.Reconstruction_cff import *
process.selectTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("generalTracks"),
    cut = cms.string('pt > 0.3 & pt<500 & numberOfValidHits > 0')
)
process.trackCandidates = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("selectTracks"),
    particleType = cms.string('pi+')
)
icTrackJetParameters = cms.PSet(
    src = cms.InputTag("trackCandidates"),
    verbose = cms.untracked.bool(False),
    jetPtMin = cms.double(1.0),
    inputEtMin = cms.double(0.3),
    seedThreshold = cms.double(1.0),
    debugLevel = cms.untracked.int32(0),
    jetType = cms.untracked.string('BasicJet'),
    inputEMin = cms.double(0.0)
)

process.ic3TrackJets = cms.EDProducer("IterativeConeJetProducer",
    icTrackJetParameters,
    alias = cms.untracked.string('ic3TrackJets'),
    coneRadius = cms.double(0.3)
                           
)
#JetTrackAssociation
process.TrackJetTracksAssociatorAtVertex = cms.EDFilter("JetTracksAssociatorAtVertex",
    #j2tParametersVX,
    tracks = cms.InputTag("generalTracks"),
    #DEF:
    coneSize = cms.double(0.3),
    #coneSize = cms.double(1.0), 
    jets = cms.InputTag("ic3TrackJets")
)
#impact parameter
process.load("RecoBTag.ImpactParameter.impactParameter_cfi")
process.TimpactParameterTagInfos = process.impactParameterTagInfos.clone();
process.TimpactParameterTagInfos.jetTracks = cms.InputTag("TrackJetTracksAssociatorAtVertex")
#secondary vertices
from RecoBTag.SecondaryVertex.secondaryVertexTagInfos_cfi import *
process.secondaryVertexTagInfosTrackJetsDef = secondaryVertexTagInfos.clone()
process.secondaryVertexTagInfosTrackJetsDef.trackIPTagInfos = cms.InputTag("TimpactParameterTagInfos","","MCJP")
###process.secondaryVertexTagInfosTrackJetsDef.trackSelection.jetDeltaRMax = 1.0
#secondaryVertexTagInfosTrackJetsDef.trackSelection.maxDistToAxis = 99999.0
###process.secondaryVertexTagInfosTrackJetsDef.vertexCuts.maxDeltaRToJetAxis = 1.0
#secondaryVertexTagInfosTrackJetsDef.vertexCuts.distVal2dMin = -99999
#secondaryVertexTagInfosTrackJetsDef.vertexCuts.distSig2dMin = -99999
#secondaryVertexTagInfosTrackJetsDef.vertexCuts.fracPV = 2
#secondaryVertexTagInfosTrackJetsDef.vertexCuts.massMax = 99999
#secondaryVertexTagInfosTrackJetsDef.useBeamConstraint = False
##############


##############
#B JETS


process.BJets = cms.EDProducer('MCBJetProducer',
)
from RecoJets.JetAssociationProducers.j2tParametersVX_cfi import *
process.BJetTracksAssociatorAtVertex = cms.EDFilter("JetTracksAssociatorAtVertex",
    j2tParametersVX,
    jets = cms.InputTag("BJets")
)
#impact parameter
#process.load("RecoBTag.ImpactParameter.impactParameter_cfi")
process.BimpactParameterTagInfos = process.impactParameterTagInfos.clone();
process.BimpactParameterTagInfos.jetTracks = cms.InputTag("BJetTracksAssociatorAtVertex")
#what about maximum transverse impact parameter? why?

#secondary vertices
from RecoBTag.SecondaryVertex.secondaryVertexTagInfos_cfi import *
process.secondaryVertexTagInfosBJetsDef = secondaryVertexTagInfos.clone()
#process.secondaryVertexTagInfosBJetsDef.trackSelection.jetDeltaRMax = 1.0
#process.secondaryVertexTagInfosBJetsDef.trackSelection.maxDistToAxis = 99999.0
#process.secondaryVertexTagInfosBJetsDef.vertexCuts.maxDeltaRToJetAxis = 1.0
process.secondaryVertexTagInfosBJetsDef.trackIPTagInfos = cms.InputTag("BimpactParameterTagInfos","","MCJP")
##############

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "IDEAL_V9::All"


process.p = cms.Path(process.selectTracks*process.trackCandidates*process.ic3TrackJets*process.TrackJetTracksAssociatorAtVertex*process.TimpactParameterTagInfos*process.secondaryVertexTagInfosTrackJetsDef*process.BJets*process.offlineBeamSpot*process.BJetTracksAssociatorAtVertex*process.BimpactParameterTagInfos*process.secondaryVertexTagInfosBJetsDef)


process.e = cms.EndPath(process.out)
