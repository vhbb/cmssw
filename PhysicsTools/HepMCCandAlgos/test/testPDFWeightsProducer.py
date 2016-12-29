# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: testpdf -s NONE --no_exec --conditions auto:run2_mc -n -1 --filein file:~/work/CMSSWpdfweight/events_orig.lhe
import FWCore.ParameterSet.Config as cms

process = cms.Process('PDFweights')

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# Input source
process.source = cms.Source("PoolSource",
    # fileNames = cms.untracked.vstring('/store/mc/RunIISpring15FSPremix/SMS-T1bbbb_mGluino-1150_mLSP-400to975-1100to1125_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/MCRUN2_74_V9-v1/50000/320B2CE3-BE5D-E511-8C9E-B083FED76C6C.root')
    # fileNames = cms.untracked.vstring('/store/mc/RunIISpring16MiniAODv2/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/30000/021605F2-751A-E611-B237-0CC47AC33AB2.root')
    fileNames = cms.untracked.vstring('/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/00000/0064B539-803A-E611-BDEA-002590D0B060.root')
    # fileNames = cms.untracked.vstring('/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00000/E8090432-8628-E611-8713-001EC9ADFDC9.root')
)

process.options = cms.untracked.PSet(

)

process.out = cms.OutputModule( "PoolOutputModule",
  fileName = cms.untracked.string( "PDFWeightsProducer.root" ),
  outputCommands= cms.untracked.vstring(
    # "drop *",
    # "keep *_genParticles_*_*"
  )
)


# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.testpdf = cms.EDProducer(
                                 "PDFWeightsProducer",
                                 # mc2hessianCSV = cms.FileInPath('PhysicsTools/HepMCCandAlgos/data/NNPDF30_lo_as_0130_hessian_60.csv'), #MC2Hessian transformation matrix
                                 mc2hessianCSV = cms.FileInPath('PhysicsTools/HepMCCandAlgos/data/NNPDF30_nlo_as_0118_hessian_60.csv'), #MC2Hessian transformation matrix
                                 src = cms.InputTag("source"),
                                 # putting to 0 at least one of the parameters below will enforce hard-coded settings based on the mc2hessianCSV input file
                                 pdfWeightLHAnumber = cms.uint32(0), #index of first mc replica weight (careful, this should not be the nominal weight, which is repeated in some mc samples).  The majority of run2 LO madgraph_aMC@NLO samples with 5fs matrix element and pdf would use index 10, corresponding to pdf set 263001, the first alternate mc replica for the nominal pdf set 263000 used for these samples
                                 nPdfWeights = cms.uint32(0), #number of input weights
                                 nPdfEigWeights = cms.uint32(0), #number of output weights
                                 )

process.p = cms.Path(
    process.testpdf
)

process.o = cms.EndPath(
    process.out 
)

