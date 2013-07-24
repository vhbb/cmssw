import FWCore.ParameterSet.Config as cms

doubleElectronHLTMatching = cms.EDProducer("DoubleElectronHLTMatching",
                                           InputCollection = cms.InputTag("gsfElectrons"),
                                           TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                                           HLTTriggerSummaryAOD = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                                           TriggerPaths = cms.vstring("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30",
                                                                      "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v7",
                                                                      "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17",
                                                                      "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17",
                                                                      ),
                                           RecoCuts   = cms.string(""),
                                           HLTCuts    = cms.string(""),
                                           TagLeg     = cms.bool(True),
                                           DeltaR     = cms.double(0.2),
                                           DoMatching = cms.bool(True)
                                           )

