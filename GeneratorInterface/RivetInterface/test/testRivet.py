import FWCore.ParameterSet.Config as cms
process = cms.Process("runRivetAnalysis")

process.options   = cms.untracked.PSet(                           
    allowUnscheduled = cms.untracked.bool(False)
) 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.source = cms.Source("PoolSource",  fileNames = cms.untracked.vstring(

# '/store/mc/RunIISpring16reHLT80/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/AODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/00000/042AE626-B644-E611-ADC8-20CF3056170E.root'
'/store/mc/RunIISpring16MiniAODv2/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/00000/10F007BD-8446-E611-A2DD-0025909091AA.root'

) 

)
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.myGenerator = cms.EDProducer("GenParticles2HepMCConverterHTXS",
    # genParticles = cms.InputTag("genParticles"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    genEventInfo = cms.InputTag("generator"),
)

#process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")
process.rivetAnalyzer = cms.EDAnalyzer('RivetAnalyzer',
  AnalysisNames = cms.vstring('HiggsTemplatCrossSections'),
  HepMCCollection = cms.InputTag('myGenerator','unsmeared'),
  UseExternalWeight = cms.bool(False),
  GenEventInfoCollection = cms.InputTag('generator'),
  useLHEweights = cms.bool(False),
  LHEweightNumber = cms.int32(0),
  LHECollection = cms.InputTag('source'),
  CrossSection = cms.double(1000),
  DoFinalize = cms.bool(True),
  ProduceDQMOutput = cms.bool(False),
  OutputFile = cms.string('out.aida')
)

process.rivetAnalyzerHTXS = cms.EDAnalyzer('HTXSRivetAnalyzer',
  AnalysisNames = cms.vstring('HiggsTemplateCrossSections'),
  HepMCCollection = cms.InputTag('myGenerator','unsmeared'),
  UseExternalWeight = cms.bool(False),
  GenEventInfoCollection = cms.InputTag('generator'),
  useLHEweights = cms.bool(False),
  LHEweightNumber = cms.int32(0),
  LHECollection = cms.InputTag('source'),
  CrossSection = cms.double(1000),
  DoFinalize = cms.bool(True),
  ProduceDQMOutput = cms.bool(False),
  OutputFile = cms.string('out.aida')
)

# process.p = cms.Path(process.myGenerator*process.rivetAnalyzer)
process.p = cms.Path(process.myGenerator
*process.rivetAnalyzerHTXS
)
