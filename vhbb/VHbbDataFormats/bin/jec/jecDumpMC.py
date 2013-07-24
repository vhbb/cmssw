import FWCore.ParameterSet.Config as cms
process = cms.Process("jectxt")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# define your favorite global tag
process.GlobalTag.globaltag = 'START52_V9::All'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.source = cms.Source("EmptySource")
process.readAK5PFchs    = cms.EDAnalyzer('JetCorrectorDBReader',  
        # below is the communication to the database 
        payloadName    = cms.untracked.string('AK5PFchs'),
        # this is used ONLY for the name of the printed txt files. You can use any name that you like, 
        # but it is recommended to use the GT name that you retrieved the files from.
        globalTag      = cms.untracked.string('Summer12V3MC'),
        printScreen    = cms.untracked.bool(False),
        createTextFile = cms.untracked.bool(True)
)
process.p = cms.Path(process.readAK5PFchs)

process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Summer12_V3_MC_AK5PFchs'),
            # tag    = cms.string('JetCorrectorParametersCollection_Summer12_V3_MC_AK5PF'),
            label  = cms.untracked.string('AK5PFchs')
            ),

      ## here you add as many jet types as you need
      ## note that the tag name is specific for the particular sqlite file 
      ), 
      connect = cms.string('sqlite:Summer12_V3_MC.db')
     # uncomment above tag lines and this comment to use MC JEC
     # connect = cms.string('sqlite:Summer12_V3_MC.db')
)
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

#files
# wget http://faculty.physics.tamu.edu/eusebi/jec/52X/SQLfiles/Summer12_V3_DATA.db
# wget http://faculty.physics.tamu.edu/eusebi/jec/52X/SQLfiles/Summer12_V3_MC.db
