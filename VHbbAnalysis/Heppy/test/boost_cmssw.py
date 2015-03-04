import FWCore.ParameterSet.Config as cms

process = cms.Process("EX")
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:///scratch/gregor/TTJets_MSDecaysCKM_central_Tune4C_13TeV_MiniAOD.root")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.OUT = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    outputCommands = cms.untracked.vstring(['drop *'])
)
process.endpath= cms.EndPath(process.OUT)


# Import some jet clustering defaults
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.PFJetParameters_cfi import *

# Select candidates that would pass CHS requirements
process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))

# HEPTopTagger (MultiR)
process.looseMultiRHTT = cms.EDProducer(
    "HTTTopJetProducer",
    PFJetParameters,
    AnomalousCellParameters,
    multiR         = cms.bool(True),
    algorithm      = cms.int32(1),
    jetAlgorithm   = cms.string("CambridgeAachen"),
    rParam         = cms.double(1.5),
    mode           = cms.int32(4),
    minFatjetPt    = cms.double(200.),
    minCandPt      = cms.double(200.),
    minSubjetPt    = cms.double(30.),
    writeCompound  = cms.bool(True),
    minCandMass    = cms.double(0.),
    maxCandMass    = cms.double(1000),
    massRatioWidth = cms.double(100.),
    minM23Cut      = cms.double(0.),
    minM13Cut      = cms.double(0.),
    maxM13Cut      = cms.double(2.))
process.looseMultiRHTT.src = cms.InputTag("chs")
process.looseMultiRHTT.jetPtMin = cms.double(200.)

# CA, R=1.5, pT > 200 GeV
process.ca15PFJetsCHS = cms.EDProducer(
        "FastjetJetProducer",
        PFJetParameters,
        AnomalousCellParameters,
        jetAlgorithm = cms.string("CambridgeAachen"),
        rParam       = cms.double(1.5))
process.ca15PFJetsCHS.src = cms.InputTag("chs")
process.ca15PFJetsCHS.jetPtMin = cms.double(200.)


# Calculate tau1, tau2 and tau3 for ungroomed CA R=1.5 jets
process.ca15PFJetsCHSNSubjettiness  = cms.EDProducer("NjettinessAdder",
                                                     src=cms.InputTag("ca15PFJetsCHS"),
                                                     cone=cms.double(1.5),
                                                     Njets = cms.vuint32(1,2,3)
)

# Apply trimming to CA R=1.5 jets
process.ca15PFTrimmedJetsCHS = process.ca15PFJetsCHS.clone(
    useTrimming = cms.bool(True),
    rFilt = cms.double(0.2),
    trimPtFracMin = cms.double(0.06),
    useExplicitGhosts = cms.bool(True))

# Let CMSSW take care of scheduling 
process.options = cms.untracked.PSet(     
        allowUnscheduled = cms.untracked.bool(True)
)


# Write all outputs
process.OUT.outputCommands.append("keep *_ca15PFJetsCHS_*_EX")
process.OUT.outputCommands.append("keep *_ca15PFJetsCHSNSubjettiness_*_EX")
process.OUT.outputCommands.append("keep *_ca15PFTrimmedJetsCHS_*_EX")
process.OUT.outputCommands.append("keep *_looseMultiRHTT_*_EX")

