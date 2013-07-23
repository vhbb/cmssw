import FWCore.ParameterSet.Config as cms

process = cms.Process("BCP")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/CMSSW_3_5_6/src/filterSkim/filterQCD30.root'
        'file:/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/CMSSW_3_5_6/src/filterSkim/filterMC_DiJet_15_20.root'
        #'file:/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/CMSSW_3_5_6/src/filterSkim/filterMC_DiJet_20_30.root'


    )
)

process.bcandidates = cms.EDProducer('BCandidateProducer',
        src = cms.InputTag('selectedVertices','',''),
        primaryVertices = cms.InputTag('offlinePrimaryVerticesWithBS','',''),
	minDRUnique = cms.untracked.double(0.4),
        minvecSumIMifsmallDRUnique = cms.untracked.double(5.5),
        minCosPAtomerge = cms.untracked.double(0.99),
        maxPtreltomerge = cms.untracked.double(7777.0)

)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('filterMCQCD30_bcandidates.root'),
)


  
process.p = cms.Path(process.bcandidates)

process.e = cms.EndPath(process.out)
