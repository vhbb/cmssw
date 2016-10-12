import FWCore.ParameterSet.Config as cms
process = cms.Process("runRivetAnalysis")

process.options   = cms.untracked.PSet(                           
    allowUnscheduled = cms.untracked.bool(False)
) 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source = cms.Source("PoolSource",  fileNames = cms.untracked.vstring(
'/store/mc/RunIISpring16MiniAODv2/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/40000/2455EB86-8F38-E611-8E4B-0090FAA58D84.root'
#'/store/mc/RunIISpring16MiniAODv2/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/00000/10F007BD-8446-E611-A2DD-0025909091AA.root'
) 
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.myGenerator = cms.EDProducer("GenParticles2HepMCConverterHTXS",
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    genEventInfo = cms.InputTag("generator"),
)

process.rivetProducerHTXS = cms.EDProducer('HTXSRivetProducer',
  HepMCCollection = cms.InputTag('myGenerator','unsmeared'),
  ProductionMode = cms.string('GGF'),
  #ProductionMode = cms.string('QQ2ZH'),
)

process.p = cms.Path(process.myGenerator*process.rivetProducerHTXS)

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('testHTXSRivet.root')
)
process.o = cms.EndPath( process.out )
