## Trigger menu in MC: /dev/CMSSW_7_4_0/GRun/V41 (default in CMSSW_7_4_1)
## http://cmslxr.fnal.gov/source/HLTrigger/Configuration/python/HLT_GRun_cff.py?v=CMSSW_7_4_1

triggerTable = {
    "ZnnHbbAll" : [
        "HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDLoose_BTagCSV0p7_v*",
        "HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDLoose_v*",
        "HLT_PFMET90_PFMHT90_IDLoose_v*",
        "HLT_PFMET100_PFMHT100_IDLoose_v*",
        "HLT_PFMET110_PFMHT110_IDLoose_v*",
        "HLT_PFMET120_PFMHT120_IDLoose_v*",
        "HLT_PFMET170_NoiseCleaned_v*",
        "HLT_DiCentralPFJet70_PFMET120_NoiseCleaned_v*",
        "HLT_PFHT350_PFMET120_NoiseCleaned_v*",
    ],
    "ZnnHbbHighLumi" : [
        "HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDLoose_BTagCSV0p7_v*",
        "HLT_PFMET120_PFMHT120_IDLoose_v*",
        "HLT_PFMET170_NoiseCleaned_v*",
    ],
    "ZnnHbbLowLumi" : [
        "HLT_PFMET90_PFMHT90_IDLoose_v*",
    ],

    "ZeeHbbAll" : [
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
        "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
        "HLT_DoubleEle24_22_eta2p1_WP75_Gsf_v*",
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
        "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
        "HLT_Ele32_eta2p1_WP75_Gsf_v*",
        "HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v*",
        "HLT_Ele27_eta2p1_WP85_Gsf_HT200_v*",
        "HLT_Ele27_WP85_Gsf_v*",
        "HLT_Ele27_eta2p1_WP75_Gsf_v*",
        "HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v*",
        "HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",
        "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*",
    ],
    "ZeeHbbHighLumi" : [
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
        "HLT_DoubleEle24_22_eta2p1_WP75_Gsf_v*"
    ],
    "ZeeHbbLowLumi" : [
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
        "HLT_DoubleEle24_22_eta2p1_WP75_Gsf_v*"
    ],


    "WenHbbAll" : [
        "HLT_Ele32_eta2p1_WP75_Gsf_v*",
        "HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v*",
        "HLT_Ele22_eta2p1_WP75_Gsf_v*",
        "HLT_Ele27_eta2p1_WP85_Gsf_HT200_v*",
        "HLT_Ele27_WP85_Gsf_v*",
        "HLT_Ele27_eta2p1_WP75_Gsf_v*",
        "HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v*",
        "HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",
        "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*",
    ],
    "WenHbbHighLumi" : [
        "HLT_Ele32_eta2p1_WP75_Gsf_v*",
        "HLT_Ele27_eta2p1_WP85_Gsf_HT200_v*",
    ],
    "WenHbbLowLumi" : [
        "HLT_Ele27_eta2p1_WP75_Gsf_v*",
        "HLT_Ele27_eta2p1_WP85_Gsf_HT200_v*",
    ],

    "ZmmHbbAll" : [
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
        "HLT_Mu17_TkMu8_DZ_v*",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
        "HLT_DoubleIsoMu17_eta2p1_v*",
        "HLT_IsoMu24_eta2p1_v*",
        "HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v*",
        "HLT_Mu24_eta2p1_v*",
        "HLT_TkMu24_eta2p1_v*",
        "HLT_Mu24_v*",
        "HLT_IsoMu27_v*",
        "HLT_IsoTkMu27_v*",
        "HLT_TkMu27_v*",
        "HLT_Mu27_v*",
        "HLT_IsoMu20_eta2p1_v*",
        "HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v*",
        "HLT_Mu20_v*",
        "HLT_TkMu20_v*",
        "HLT_IsoMu20_v*",
        "HLT_IsoTkMu20_v*",
        "HLT_Mu20_v*",
        "HLT_TkMu20_v*",
        "HLT_Mu40_eta2p1_PFJet200_PFJet50_v*",
    ],
    "ZmmHbbHighLumi" : [
        "HLT_IsoMu24_eta2p1_v*",
        "HLT_IsoMu27_v*",
        "HLT_IsoTkMu27_v*",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
    ],
    "ZmmHbbLowLumi" : [
        "HLT_IsoMu20_v*",
        "HLT_IsoTkMu20_v*",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
    ],

    "WmnHbbAll" : [
        "HLT_IsoMu24_eta2p1_v*",
        "HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v*",
        "HLT_Mu24_eta2p1_v*",
        "HLT_TkMu24_eta2p1_v*",
        "HLT_Mu24_v*",
        "HLT_IsoMu27_v*",
        "HLT_IsoTkMu27_v*",
        "HLT_TkMu27_v*",
        "HLT_Mu27_v*",
        "HLT_IsoMu20_eta2p1_v*",
        "HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v*",
        "HLT_Mu20_v*",
        "HLT_TkMu20_v*",
        "HLT_IsoMu20_v*",
        "HLT_IsoTkMu20_v*",
        "HLT_Mu20_v*",
        "HLT_TkMu20_v*",
        "HLT_Mu40_eta2p1_PFJet200_PFJet50_v*",
        "HLT_IsoMu16_eta2p1_CaloMET30_v*",
        "HLT_Mu16_eta2p1_CaloMET30_v*",
        "HLT_PFMET120_NoiseCleaned_Mu5_v*",
        "HLT_MonoCentralPFJet80_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v*",
        "HLT_MonoCentralPFJet80_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v*",
        "HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v*",
        "HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v*",
    ],
    "WmnHbbHighLumi" : [
        "HLT_IsoMu24_eta2p1_v*",
        "HLT_IsoMu27_v*",
        "HLT_IsoTkMu27_v*",
        "HLT_IsoMu16_eta2p1_CaloMET30_v*"
    ],
    "WmnHbbLowLumi" : [
        "HLT_IsoMu20_v*",
        "HLT_IsoTkMu20_v*",
        "HLT_IsoMu16_eta2p1_CaloMET30_v*",
    ],

    "WtaunHbbAll" : [
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v*",
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v*",
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_v*",
    ],
    "WtaunHbbHighLumi" : [
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v*",
    ],
    "WtaunHbbLowLumi" : [
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v*",
    ],

    "VBFHbbAll" : [
        "HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v*",
        "HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v*",
        "HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v*",
        "HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v*",
        "HLT_QuadPFJet_VBF_v*",
        "HLT_L1_TripleJet_VBF_v*",
        "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v*",
        "HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v*",
    ],
    "VBFHbbHighLumi" : [
        "HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v*",
        "HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v*",
    ],
    "VBFHbbLowLumi" : [
        "HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v*",
        "HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v*",
    ],

    "HH4bAll" : [
        "HLT_QuadJet45_TripleCSV0p5_v*",
        "HLT_QuadJet45_DoubleCSV0p5_v*",
        "HLT_DoubleJet90_Double30_TripleCSV0p5_v*",
        "HLT_DoubleJet90_Double30_DoubleCSV0p5_v*",
    ],
    "HH4bHighLumi" : [
        "HLT_QuadJet45_TripleCSV0p5_v*",
        "HLT_DoubleJet90_Double30_TripleCSV0p5_v*",
    ],
    "HH4bLowLumi" : [
        "HLT_QuadJet45_TripleCSV0p5_v*",
        "HLT_DoubleJet90_Double30_TripleCSV0p5_v*",
    ],
    "ttHleptonic" : [
        "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
        "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*",
    ],
    "ttHhardonicAll" : [
        "HLT_PFHT450_SixJet40_PFBTagCSV_v*",
        "HLT_PFHT450_SixJet40_v*",
        "HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v*",
        "HLT_PFHT400_SixJet30_v*",
        "HLT_PFHT350_v*", 
    ],
    "ttHhardonicHighLumi" : [
        "HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v*",
        "HLT_PFHT450_SixJet40_PFBTagCSV_v*",
    ],
    "ttHhardonicLowLumi" : [
        "HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v*",
        "HLT_PFHT450_SixJet40_PFBTagCSV_v*",
    ],

    "hadronic" : [
        "HLT_PFHT750_4Jet_v*",
        "HLT_PFHT900_v*",
        "HLT_PFJet40_v*",
        "HLT_PFJet60_v*",
        "HLT_PFJet80_v*",
        "HLT_PFJet140_v*",
        "HLT_PFJet200_v*",
        "HLT_PFJet260_v*",
        "HLT_PFJet320_v*",
        "HLT_PFJet400_v*",
        "HLT_PFJet450_v*",
        "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p41_v*",
        "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v*",
        "HLT_AK8PFJet360_TrimMass30_v*",
    ],
}
