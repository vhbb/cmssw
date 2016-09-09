import FWCore.ParameterSet.Config as cms
process = cms.Process("runRivetAnalysis")

process.options   = cms.untracked.PSet(                           
    allowUnscheduled = cms.untracked.bool(False)
) 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.source = cms.Source("PoolSource",  fileNames = cms.untracked.vstring(

# '/store/mc/RunIISpring16reHLT80/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/AODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/00000/042AE626-B644-E611-ADC8-20CF3056170E.root'
# '/store/mc/RunIISpring16MiniAODv2/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/00000/10F007BD-8446-E611-A2DD-0025909091AA.root'
# '/store/mc/RunIISpring16MiniAODv2/ttHTobb_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/40000/0089CC67-6338-E611-947D-0025904C4E2A.root'
# '/store/mc/RunIISpring16MiniAODv2/GluGluHToBB_M125_13TeV_powheg_pythia8_CUETP8M1Down/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/20000/1E80D500-7A38-E611-A659-1418774124DE.root'
# '/store/mc/RunIISpring16MiniAODv2/GluGluHToBB_M125_13TeV_powheg_herwigpp/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/20000/0A8D2E42-6938-E611-81E8-00266CFFA23C.root'
# '/store/mc/RunIISpring16MiniAODv2/ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_herwigpp/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/40000/2AACBB94-8138-E611-9160-008CFA197FAC.root'
# '/store/mc/RunIISpring16MiniAODv2/WplusH_HToBB_WToLNu_M125_13TeV_powheg_herwigpp/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/20000/3C858D7E-2539-E611-A599-20CF3027A57D.root'
'/store/mc/RunIISpring16MiniAODv2/ZH_HToBB_ZToLL_M125_13TeV_powheg_herwigpp/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/70000/2AA52F8E-A238-E611-B098-00259020084C.root'
# '/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_herwigpp/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/040F30E6-9C38-E611-B0F7-0CC47A1DF80C.root'
) 

)
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.myGenerator = cms.EDProducer("GenParticles2HepMCConverterHTXS",
    # genParticles = cms.InputTag("genParticles"),
    genParticles = cms.InputTag("prunedGenParticles"),
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

#process.p = cms.Path(process.myGenerator*process.rivetAnalyzer)
# process.p = cms.Path(process.myGenerator*process.rivetAnalyzerHTXS)
process.p = cms.Path(process.myGenerator)
