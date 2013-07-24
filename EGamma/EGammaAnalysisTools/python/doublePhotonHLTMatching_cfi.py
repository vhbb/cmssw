import FWCore.ParameterSet.Config as cms

doublePhotonHLTMatching = cms.EDProducer("DoublePhotonHLTMatching",
                                         InputCollection = cms.InputTag("photons"),
                                         TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                                         HLTTriggerSummaryAOD = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                                         TriggerPaths = cms.vstring("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL",
                                                                    ),
                                         RecoCuts   = cms.string(""),
                                         HLTCuts    = cms.string(""),
                                         TagLeg     = cms.bool(True),
                                         DeltaR     = cms.double(0.2),
                                         DoMatching = cms.bool(True)
                                         )

