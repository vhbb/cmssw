import FWCore.ParameterSet.Config as cms

process = cms.Process("SSV")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

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

##############
#Include the jet corrections
#process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08_cff")

# set the record's IOV. Must be defined once. Choose ANY correction service. #
#process.prefer("L2L3JetCorrectorIC5Calo")
##############


##############
process.myProducerLabel = cms.EDProducer('MCSecVertexProducer',
    jetsrc  = cms.InputTag('iterativeCone5CaloJets')
    #jetsrc  = cms.InputTag('L2L3CorJetIC5Calo')
)
##############

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('output.root')
)

process.p = cms.Path(process.myProducerLabel)
#process.p = cms.Path(process.L2L3CorJetIC5Calo*process.myProducerLabel)

process.e = cms.EndPath(process.out)
