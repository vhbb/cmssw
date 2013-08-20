import FWCore.ParameterSet.Config as cms

photonHLTMatching = cms.EDProducer("PhotonHLTMatching",
                                   InputCollection = cms.InputTag("photons"),
                                   TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                                   HLTTriggerSummaryAOD = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                                   TriggerPaths = cms.vstring("HLT_Photon20_CaloIdVL",),
                                   ModuleLabels = cms.vstring("hltEG20CaloIdVLHEFilter",),
                                   RecoCuts = cms.string(""),
                                   HLTCuts = cms.string(""),
                                   DeltaR = cms.double(0.2),
                                   DoMatching = cms.bool(True)
                                   )

