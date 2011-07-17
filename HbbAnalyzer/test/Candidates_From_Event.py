#import configurations
import FWCore.ParameterSet.Config as cms

import os 


# define the process
process = cms.Process("Candidates")
process.source = cms.Source("PoolSource",
 fileNames = cms.untracked.vstring('file:/gpfs/gpfsddn/cms/user/boccali/hbb/ttjets/VHbbPAT.edm_5_0_eqA.root') ) 


process.out = cms.OutputModule(
    'PoolOutputModule',
    fileName       = cms.untracked.string('Cand.root'),
#    outputCommands = cms.untracked.vstring('drop *','keep *_HbbAnalyzerNew_*_*', 'keep *_hbbCandidates_*_*'),
    outputCommands = cms.untracked.vstring('drop *','keep *_*_*_Candidates'),
    dropMetaData = cms.untracked.string('ALL'),
    splitLevel = cms.untracked.int32(0)
    )


process.hbbCandidates = cms.EDProducer("HbbCandidateFinder",
				       VHbbEventLabel = cms.InputTag(""),
				       verbose = cms.bool(False) ,
				       jetPtThreshold = cms.double(30.)
				      )


# drop the meta data for dropped data
#process.out.dropMetaData = cms.string("DROPPED")

#process.out.fileName = '/tigress-hsm/dlopes/PatEDM.root'


# define path 'p'
process.p = cms.Path(	     process.hbbCandidates                     )


process.e = cms.EndPath(process.out)

#
