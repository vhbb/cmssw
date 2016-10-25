import FWCore.ParameterSet.Config as cms
process = cms.Process("runRivetAnalysis")

process.options   = cms.untracked.PSet(                           
    allowUnscheduled = cms.untracked.bool(False)
) 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",  fileNames = cms.untracked.vstring(

#'/store/mc/RunIISpring16reHLT80/VBFHToGG_M-125_13TeV_powheg_pythia8/AODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/00000/001631C3-4D3D-E611-AAB5-002590821180.root',
#'/store/mc/RunIISpring16reHLT80/VBFHToGG_M-125_13TeV_powheg_pythia8/AODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/00000/161298D1-3E3D-E611-9A56-0025904B893A.root',
#'/store/mc/RunIISpring16reHLT80/VBFHToGG_M-125_13TeV_powheg_pythia8/AODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/00000/18F2540A-493D-E611-8BCB-002590821180.root'
'/store/mc/RunIISpring16MiniAODv2/VBFHToGG_M-125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/00000/36A162AC-073E-E611-BD6A-0025901F96F4.root'
),

)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.rivetProducerHTXS = cms.EDProducer('HTXSRivetProducer',
  HepMCCollection = cms.InputTag('myGenerator','unsmeared'),
  ProductionMode = cms.string('VBF'),
)

#AOD
#process.myGenerator = cms.EDProducer("GenParticles2HepMCConverterHTXS",
#    genParticles = cms.InputTag("genParticles"),
#    genEventInfo = cms.InputTag("generator"),
#)
#process.p = cms.Path(process.myGenerator*process.rivetProducerHTXS)

#MINIAOD
process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
    inputPruned = cms.InputTag("prunedGenParticles"),
    inputPacked = cms.InputTag("packedGenParticles"),
)
process.myGenerator = cms.EDProducer("GenParticles2HepMCConverterHTXS",
    genParticles = cms.InputTag("mergedGenParticles"),
    genEventInfo = cms.InputTag("generator"),
)
process.p = cms.Path(process.mergedGenParticles*process.myGenerator*process.rivetProducerHTXS)

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *','keep *_*_*_runRivetAnalysis','keep *_generator_*_*','keep *_externalLHEProducer_*_*'),
    #fileName = cms.untracked.string('testHTXSRivet_VBFHgg_AOD_testMerger.root')
    fileName = cms.untracked.string('testHTXSRivet_VBFHgg_MINIAOD_testMerger.root')
)
process.o = cms.EndPath( process.out )
