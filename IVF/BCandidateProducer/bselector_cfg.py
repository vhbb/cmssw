import FWCore.ParameterSet.Config as cms

process = cms.Process("BCP")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/CMSSW_3_5_6/src/filterSkim/filterQCD30.root'
        #'file:/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/CMSSW_3_5_6/src/filterSkim/filterMC_DiJet_15_20.root'
        #'file:/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/CMSSW_3_5_6/src/filterSkim/filterMC_DiJet_20_30.root'


    )
)

process.bcandidates = cms.EDProducer('BSelector',
        src = cms.InputTag('selectedVertices','',''),
        primaryVertices = cms.InputTag('offlinePrimaryVerticesWithBS','',''),
        mininvMass = cms.untracked.double(0.9),
        mindist3D = cms.untracked.double(0.05),
        maxdist3D = cms.untracked.double(5.0),
        maxdist2D = cms.untracked.double(777.0),
        mindistSig3D = cms.untracked.double(10.0),
        cosFdTdsum = cms.untracked.double(0.98),
        minndof = cms.untracked.double(2.0),
        maxConeSize = cms.untracked.double(0.5),
        minpt = cms.untracked.double(10.0),
        maxSumScalarMass = cms.untracked.double(4.5),
        minSumCosFdTdsum = cms.untracked.double(1.995),
        #minvecSumIMifsmallDRUnique = cms.untracked.double(4.5) #was 6.0 and should go 3.5                                 

)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('filterMCQCD30_test.root'),
)


  
process.p = cms.Path(process.bcandidates)

process.e = cms.EndPath(process.out)
