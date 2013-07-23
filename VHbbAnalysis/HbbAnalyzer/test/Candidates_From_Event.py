#import configurations
import FWCore.ParameterSet.Config as cms

import os 


# define the process
process = cms.Process("Candidates")
process.source = cms.Source("PoolSource",
 fileNames = cms.untracked.vstring('file:PAT.edm_194_0_mAx.root') ) 


process.out = cms.OutputModule(
    'PoolOutputModule',
    fileName       = cms.untracked.string('Cand.root'),
#    outputCommands = cms.untracked.vstring('drop *','keep *_HbbAnalyzerNew_*_*', 'keep *_hbbCandidates_*_*'),
    outputCommands = cms.untracked.vstring('keep *','drop VHbbEvent_*_*_*'),
    dropMetaData = cms.untracked.string('ALL'),
    splitLevel = cms.untracked.int32(0),
        SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('candidates')
    )
)


process.hbbHighPtCandidates = cms.EDProducer("HbbCandidateFinder",
				       VHbbEventLabel = cms.InputTag(""),
				       verbose = cms.bool(False) ,
				       useHighestPtHiggs = cms.bool(True)
				       jetPtThreshold = cms.double(30.),
             			       actAsAFilter = cms.bool(False)
				      )


# drop the meta data for dropped data
#process.out.dropMetaData = cms.string("DROPPED")

#process.out.fileName = '/tigress-hsm/dlopes/PatEDM.root'


# define path 'p'
process.p = cms.Path(	     process.hbbHighPtCandidates                     )


process.e = cms.EndPath(process.out)



#
