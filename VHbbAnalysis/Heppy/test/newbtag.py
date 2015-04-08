import FWCore.ParameterSet.Config as cms


process = cms.Process("EX")
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:../ZLL-8A345C56-6665-E411-9C25-1CC1DE04DF20.root")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.OUT = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    outputCommands = cms.untracked.vstring(['drop *'])
)
process.endpath= cms.EndPath(process.OUT)

#process.OUT.outputCommands.append("keep *_ak5PFJetsCHS_*_EX")
# As tracks are not stored in miniAOD, and b-tag fwk for CMSSW < 72X does not accept candidates
# we need to recreate tracks and pv for btagging in standard reco format:
process.load('RecoBTag.Configuration.RecoBTag_cff')
process.load('RecoJets.Configuration.RecoJetAssociations_cff')
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.GlobalTag.globaltag = 'PLS170_V7AN1::All'

process.ak4JetTracksAssociatorAtVertexPF.jets = cms.InputTag("slimmedJets")
process.ak4JetTracksAssociatorAtVertexPF.tracks = cms.InputTag("unpackedTracksAndVertices")
process.impactParameterTagInfos.primaryVertex = cms.InputTag("unpackedTracksAndVertices")
process.inclusiveSecondaryVertexFinderTagInfos.extSVCollection = cms.InputTag("unpackedTracksAndVertices","secondary","")
process.combinedSecondaryVertex.trackMultiplicityMin = 1

process.combinedSecondaryVertexV2.calibrationRecords = cms.vstring(
'CombinedSVV2MVA_RecoVertex',
'CombinedSVV2MVA_PseudoVertex',
'CombinedSVV2MVA_NoVertex'
)


process.OUT.outputCommands.append("keep *_combinedInclusiveSecondaryVertexV2BJetTags_*_EX")
process.OUT.outputCommands.append("keep *_combinedSecondaryVertexBJetTags_*_EX")

process.options = cms.untracked.PSet( 
        wantSummary = cms.untracked.bool(True), # while the timing of this is not reliable in unscheduled mode, it still helps understanding what was actually run 
        allowUnscheduled = cms.untracked.bool(True)
)
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.BTauMVAJetTagComputerRecord = cms.ESSource("PoolDBESSource",
process.CondDBSetup,
timetype = cms.string('runnumber'),
toGet = cms.VPSet(cms.PSet(
record = cms.string('BTauGenericMVAJetTagComputerRcd'),
                tag = cms.string('MVAJetTags_620SLHCX')
)),
connect = cms.string('sqlite_file:MVAJetTags_620SLHCX_Phase1And2Upgrade.db'),
BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService')
)
process.es_prefer_BTauMVAJetTagComputerRecord = cms.ESPrefer("PoolDBESSource","BTauMVAJetTagComputerRecord")

