import FWCore.ParameterSet.Config as cms

process = cms.Process("MCPA")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.Tracer = cms.Service("Tracer")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    skipEvents = cms.untracked.uint32(2867),                
    fileNames = cms.untracked.vstring(
    'file:/shome/arizzi/twoBCone030810_95.root'
    )
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('newCutsC030810_4.root'),
    SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('p')
    ),	    
    outputCommands = cms.untracked.vstring(
'drop *',
'keep *_*_*_MCP',
'keep *_*_*_MCJP',
'keep *_*_*_MCPA',
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
#impact parameter, jettrackassociation
process.load("RecoBTag.ImpactParameter.impactParameter_cfi")
process.load("RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi")
process.ic5JetTracksAssociatorAtVertex.coneSize = 1.0
#GONE IN 225: process.impactParameterTagInfos.maximumDistanceToJetAxis = 9999.0
#GONE IN 225: process.impactParameterTagInfos.maximumDecayLength = 9999.0
#process.impactParameterTagInfos.jetTracks = cms.InputTag("BJetTracksAssociatorAtVertex")
##############

##############
#secondary vertices
#C10
from RecoBTag.SecondaryVertex.secondaryVertexTagInfos_cfi import *
secondaryVertexTagInfos.trackSelection.jetDeltaRMax = 1.0
secondaryVertexTagInfos.trackSelection.maxDistToAxis = 99999.0
secondaryVertexTagInfos.vertexCuts.maxDeltaRToJetAxis = 99999.0
#C08
process.secondaryVertexTagInfosC08 = secondaryVertexTagInfos.clone()
process.secondaryVertexTagInfosC08.trackSelection.jetDeltaRMax = 0.8
process.secondaryVertexTagInfosC08.trackSelection.maxDistToAxis = 99999.0
process.secondaryVertexTagInfosC08.vertexCuts.maxDeltaRToJetAxis = 99999.0
#C03
process.secondaryVertexTagInfosC03 = secondaryVertexTagInfos.clone()
#process.secondaryVertexTagInfosC03.trackSelection.jetDeltaRMax = 0.3
process.secondaryVertexTagInfosC03.trackSelection.maxDistToAxis = 99999.0
process.secondaryVertexTagInfosC03.vertexCuts.maxDeltaRToJetAxis = 99999.0

#NEW
process.secondaryVertexTagInfosC03.vertexCuts.distVal2dMin = -99999
process.secondaryVertexTagInfosC03.vertexCuts.distSig2dMin = -99999
process.secondaryVertexTagInfosC08.vertexCuts.distVal2dMin = -99999
process.secondaryVertexTagInfosC08.vertexCuts.distSig2dMin = -99999
secondaryVertexTagInfos.vertexCuts.distVal2dMin = -99999
secondaryVertexTagInfos.vertexCuts.distSig2dMin = -99999

process.secondaryVertexTagInfosC03.vertexCuts.fracPV = 2
process.secondaryVertexTagInfosC08.vertexCuts.fracPV = 2
secondaryVertexTagInfos.vertexCuts.fracPV = 2
process.secondaryVertexTagInfosC03.vertexCuts.massMax = 99999
process.secondaryVertexTagInfosC08.vertexCuts.massMax = 99999
secondaryVertexTagInfos.vertexCuts.massMax = 99999

process.secondaryVertexTagInfosC08.useBeamConstraint = False
process.secondaryVertexTagInfosC03.useBeamConstraint = False
secondaryVertexTagInfos.useBeamConstraint = False

##############

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "IDEAL_V9::All"


#process.p = cms.Path(process.ic5JetTracksAssociatorAtVertex*process.impactParameterTagInfos*process.secondaryVertexTagInfosC08*process.secondaryVertexTagInfosC03*secondaryVertexTagInfos)
#process.p = cms.Path(process.ic5JetTracksAssociatorAtVertex*process.impactParameterTagInfos*secondaryVertexTagInfos)
process.p = cms.Path(process.ic5JetTracksAssociatorAtVertex*process.impactParameterTagInfos*process.secondaryVertexTagInfosC08*secondaryVertexTagInfos)

process.e = cms.EndPath(process.out)
