import FWCore.ParameterSet.Config as cms

process = cms.Process("MCP")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:/data2/wehrlilu/rootfiles/00ABE731-9382-DD11-9605-001D09F28D4A.root',
    'file:/data2/wehrlilu/rootfiles/027E628C-9382-DD11-A912-001617DBD472.root',
    'file:/data2/wehrlilu/rootfiles/0407E1C3-9382-DD11-92B9-001617C3B710.root',
    'file:/data2/wehrlilu/rootfiles/068C57CA-9E82-DD11-9186-001617C3B6DE.root',
    'file:/data2/wehrlilu/rootfiles/0E4F3055-9382-DD11-BC63-0019DB29C614.root',
    'file:/data2/wehrlilu/rootfiles/145FADDE-9482-DD11-ADE6-000423D99F1E.root',
    'file:/data2/wehrlilu/rootfiles/189471C1-9B82-DD11-B87E-001617E30D4A.root',
    'file:/data2/wehrlilu/rootfiles/1A8A1572-9582-DD11-BD78-000423D6C8E6.root',
    'file:/data2/wehrlilu/rootfiles/1AF94229-9782-DD11-88C4-000423D98E54.root',
    'file:/data2/wehrlilu/rootfiles/1C304281-9382-DD11-A456-001617E30D52.root',
    'file:/data2/wehrlilu/rootfiles/3A4B8C8A-9382-DD11-BEF5-000423D33970.root',
    'file:/data2/wehrlilu/rootfiles/3E648968-9682-DD11-9BDD-001617E30D52.root',
    'file:/data2/wehrlilu/rootfiles/44044E1A-9382-DD11-BC1B-001D09F24FEC.root',
    'file:/data2/wehrlilu/rootfiles/46337288-9382-DD11-A710-001617C3B6CE.root',
    'file:/data2/wehrlilu/rootfiles/4653B665-9382-DD11-9032-001617C3B6C6.root',
    'file:/data2/wehrlilu/rootfiles/469E0459-9382-DD11-BE67-001617E30F56.root',
    'file:/data2/wehrlilu/rootfiles/46DFCBDB-9482-DD11-9837-000423D99658.root',
    'file:/data2/wehrlilu/rootfiles/50C2C758-9382-DD11-9526-000423D991F0.root',
    'file:/data2/wehrlilu/rootfiles/52E50267-9382-DD11-B765-000423D9A212.root'
    )
)

###################
#Include the jet corrections
#process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08_cff")

# set the record's IOV. Must be defined once. Choose ANY correction service. #
#process.prefer("L2L3JetCorrectorIC5Calo")
##############


##############
process.bfilter = cms.EDFilter("PythiaFilter",
    ParticleID = cms.untracked.int32(5)
)
##############

##############
process.BHadrons = cms.EDProducer('MCBHadronProducer',
    jetsrc  = cms.InputTag('iterativeCone5CaloJets')
    #jetsrc  = cms.InputTag('L2L3CorJetIC5Calo')
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
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printList = cms.EDFilter("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint = cms.untracked.int32(-1)
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('twoBCone08.root'),
    SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('p')
    ),	    
    outputCommands = cms.untracked.vstring(
'drop *',
'keep *_*_*_MCP',
'keep recoTracks_generalTracks_*_RECO',
'keep recoMuons_*_*_*',
#'keep recoTrackExtras_*_*_*',
'keep *_secondaryVertexTagInfos_*_*',
'keep recoVertexs_*_*_RECO',
'keep recoCaloJets_*_*_RECO',                          
'keep CaloTowersSorted_towerMaker_*_RECO') 
)

##############
#Cone 08
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "IDEAL_V9::All"

process.load("RecoBTag.SecondaryVertex.secondaryVertexTagInfos_cfi")

#process.load("RecoBTag.SecondaryVertex.vertexCuts_cfi")
process.load("RecoBTag.ImpactParameter.impactParameter_cfi")
process.load("RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
process.secondaryVertexTagInfos.trackSelection.jetDeltaRMax = 0.8
process.ic5JetTracksAssociatorAtVertex.coneSize = 0.8
##############


process.p = cms.Path(process.prunedGenParticles*process.SimSecVertices*process.BHadrons*process.offlineBeamSpot*process.ic5JetTracksAssociatorAtVertex*process.impactParameterTagInfos*process.secondaryVertexTagInfos)


process.e = cms.EndPath(process.out)
