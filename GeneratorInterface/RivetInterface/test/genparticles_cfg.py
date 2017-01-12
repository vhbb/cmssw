import FWCore.ParameterSet.Config as cms
process = cms.Process("MAIN")
import sys

################################################################
# Read Options
################################################################
import FWCore.ParameterSet.VarParsing as parser
opts = parser.VarParsing('analysis')
opts.register('file', '',
    parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.string, 'input file')
opts.register('globalTag', '80X_dataRun2_Prompt_ICHEP16JEC_v0', parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.string, 'global tag')
opts.register('isMini', 1, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, 'Process as data')
opts.parseArguments()
infile = opts.file.split(',')
tag  = opts.globalTag
isMini =opts.isMini

print 'globalTag    : '+str(tag)

################################################################
# Standard setup
################################################################
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("EventTree.root"),
    closeFileFast = cms.untracked.bool(True)
)

################################################################
# Message Logging, summary, and number of events
################################################################
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    # allowUnscheduled = cms.untracked.bool(True)
)


################################################################
# Input files and global tags
################################################################
process.load("CondCore.CondDB.CondDB_cfi")
from CondCore.CondDB.CondDB_cfi import *

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*infile)
    )

process.GlobalTag.globaltag = cms.string(tag)


################################################################
# Gen info
################################################################
process.mergeGenParticles = cms.EDProducer('MergedGenParticleProducer',
    inputPruned = cms.InputTag('prunedGenParticles'),
    inputPacked = cms.InputTag('packedGenParticles')
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printTree = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(1),
  printVertex = cms.untracked.bool(True),
  printOnlyHardInteraction = cms.untracked.bool(False), # Print only status=3 particles. This will not work for Pythia8, which does not have any such particles.
  src = cms.InputTag("prunedGenParticles")
)

process.printTreeMerged = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(1),
  printVertex = cms.untracked.bool(True),
  printOnlyHardInteraction = cms.untracked.bool(False), # Print only status=3 particles. This will not work for Pythia8, which does not have any such particles.
  src = cms.InputTag("mergeGenParticles")
)

process.mygenerator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("mergeGenParticles"),
    genEventInfo = cms.InputTag("generator"),
)

if not isMini:
    process.mygenerator.genParticles = cms.InputTag("genParticles")


process.rivetProducerHTXS = cms.EDProducer('HTXSRivetProducer',
  HepMCCollection = cms.InputTag('mygenerator','unsmeared'),
  ProductionMode = cms.string('VBF'),
  #ProductionMode = cms.string('QQ2ZH'),
)

import UserCode.ICHiggsTauTau.default_producers_cfi as producers

process.icEventInfoProducer = producers.icEventInfoProducer.clone(
    includeHT           = cms.bool(False),
    lheProducer         = cms.InputTag("externalLHEProducer"),
    includeJetRho       = cms.bool(True),
    inputJetRho         = cms.InputTag("fixedGridRhoFastjetAll"),
    includeLeptonRho    = cms.bool(False),
    inputLeptonRho      = cms.InputTag("fixedGridRhoFastjetAll"),
    includeVertexCount  = cms.bool(False),
    inputVertices       = cms.InputTag("offlineSlimmedPrimaryVertices"),
    includeCSCFilter    = cms.bool(False),
    inputCSCFilter      = cms.InputTag("BeamHaloSummary"),
)

# if not isData:
#     process.icEventInfoProducer.includeLHEWeights = cms.bool(True)

process.icEventProducer = producers.icEventProducer.clone()

if isMini:
    process.p = cms.Path(
        process.mergeGenParticles+
        process.printTree+
        process.printTreeMerged+
        process.mygenerator+
        process.rivetProducerHTXS+
        process.icEventInfoProducer+
        process.icEventProducer
    )
else:
    process.p = cms.Path(
        process.mygenerator+
        process.rivetProducerHTXS+
        process.icEventInfoProducer+
        process.icEventProducer
    )

process.schedule = cms.Schedule(process.p)
# print process.dumpPython()
