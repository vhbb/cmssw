import FWCore.ParameterSet.Config as cms

skim = cms.EDFilter("SkimGen",
    mcTruthCollection = cms.InputTag("genParticles"),
    genpart=cms.vint32(22), # photons
    genpt=cms.double(50.),
    Debug=cms.bool(False)
)
