import FWCore.ParameterSet.Config as cms

mvaTrigV0 = cms.EDFilter("ElectronIdMVAProducer",
                            verbose = cms.untracked.bool(False),
                            vertexTag = cms.InputTag('offlinePrimaryVertices'),
                            electronTag = cms.InputTag('gsfElectrons'),
                            reducedEBRecHitCollection = cms.InputTag('reducedEcalRecHitsEB'),
                            reducedEERecHitCollection = cms.InputTag('reducedEcalRecHitsEE'),
                            method = cms.string("BDT"),
                            mvaWeightFile = cms.vstring(
                                  "EGamma/EGammaAnalysisTools/data/Electrons_BDTG_TrigV0_Cat1.weights.xml",
                                  "EGamma/EGammaAnalysisTools/data/Electrons_BDTG_TrigV0_Cat2.weights.xml",
                                  "EGamma/EGammaAnalysisTools/data/Electrons_BDTG_TrigV0_Cat3.weights.xml",
                                  "EGamma/EGammaAnalysisTools/data/Electrons_BDTG_TrigV0_Cat4.weights.xml",
                                  "EGamma/EGammaAnalysisTools/data/Electrons_BDTG_TrigV0_Cat5.weights.xml",
                                  "EGamma/EGammaAnalysisTools/data/Electrons_BDTG_TrigV0_Cat6.weights.xml",
                            ),
                            Trig = cms.bool(True),
)

mvaNonTrigV0 = cms.EDFilter("ElectronIdMVAProducer",
                            verbose = cms.untracked.bool(False),
                            vertexTag = cms.InputTag('offlinePrimaryVertices'),
                            electronTag = cms.InputTag('gsfElectrons'),
                            reducedEBRecHitCollection = cms.InputTag('reducedEcalRecHitsEB'),
                            reducedEERecHitCollection = cms.InputTag('reducedEcalRecHitsEE'),
                            method = cms.string("BDT"),
                            mvaWeightFile = cms.vstring(
                                  "EGamma/EGammaAnalysisTools/data/Electrons_BDTG_NonTrigV0_Cat1.weights.xml",
                                  "EGamma/EGammaAnalysisTools/data/Electrons_BDTG_NonTrigV0_Cat2.weights.xml",
                                  "EGamma/EGammaAnalysisTools/data/Electrons_BDTG_NonTrigV0_Cat3.weights.xml",
                                  "EGamma/EGammaAnalysisTools/data/Electrons_BDTG_NonTrigV0_Cat4.weights.xml",
                                  "EGamma/EGammaAnalysisTools/data/Electrons_BDTG_NonTrigV0_Cat5.weights.xml",
                                  "EGamma/EGammaAnalysisTools/data/Electrons_BDTG_NonTrigV0_Cat6.weights.xml",
                            ),
                            Trig = cms.bool(False),
)
