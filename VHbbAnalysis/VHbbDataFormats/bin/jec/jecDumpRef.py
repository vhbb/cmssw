import FWCore.ParameterSet.Config as cms
process = cms.Process("jectxt")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# define your favorite global tag
process.GlobalTag.globaltag = 'GR_R_52_V7::All'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.source = cms.Source("EmptySource")
process.readAK5PFchs    = cms.EDAnalyzer('JetCorrectorDBReader',  
        # below is the communication to the database 
        payloadName    = cms.untracked.string('AK5PFchs'),
        # this is used ONLY for the name of the printed txt files. You can use any name that you like, 
        # but it is recommended to use the GT name that you retrieved the files from.
        globalTag      = cms.untracked.string('Reference'),
        printScreen    = cms.untracked.bool(False),
        createTextFile = cms.untracked.bool(True)
)
process.p = cms.Path(process.readAK5PFchs)


#files
# wget http://faculty.physics.tamu.edu/eusebi/jec/52X/SQLfiles/Summer12_V3_DATA.db
# wget http://faculty.physics.tamu.edu/eusebi/jec/52X/SQLfiles/Summer12_V3_MC.db
