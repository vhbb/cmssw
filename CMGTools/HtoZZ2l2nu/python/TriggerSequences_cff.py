import FWCore.ParameterSet.Config as cms

###
### trigger paths of interest
### Trigger evolution in 2011 cf. http://fwyzard.web.cern.ch/fwyzard/hlt/summary
###
def getTriggerPaths(version=2012) :

    #################
    # 2012 triggers #
    #################
    if(version==2012) :
        mcTrigs = []
        
        DoubleElectron = ['HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
                          'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v']
        
        DoubleMu = ['HLT_Mu17_Mu8_v',
                    'HLT_Mu17_TkMu8_v',
                    'HLT_Mu22_TkMu8_v',
                    'HLT_Mu22_TkMu22_v']
        
        MuEG = ['HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
                'HLT_Mu17_Mu8_v',
                'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
                'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v'
                ]
        
        SingleMu = []
        
        Photon = ['HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_v',
                  'HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_v',                  
                  'HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_v',
                  'HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_v',
                  'HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_v',
                  'HLT_Photon135_v4_v',
                  'HLT_Photon150_v1_v',
                  'HLT_Photon160_v1_v',
                  'HLT_Photon250_NoHE_v1_v',
                  'HLT_Photon300_NoHE_v1_v']

    #################
    # 2011 triggers #
    #################
    else:

        mcTrigs = ['HLT_IsoMu17_v5','HLT_DoubleMu7_v1', 'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2', 'HLT_Mu8_Ele17_CaloIdL_v2', 'HLT_Mu10_Ele10_CaloIdL_v3', 'HLT_Mu17_Ele8_CaloIdL_v2']

        DoubleElectron = ['HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1',                                   #start run 160404
                     'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2',
                     'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3',
                     'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4',
                     'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5',
                     'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6',
                     'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v7',
                     'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v8',
                     'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1',
                     'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2',
                     'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3',
                     'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4',  
                     'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5', #start run 167039
                     'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6',
                     'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7',
                     'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8',
                     'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9',
                     'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10'
                     ]
        DoubleMu = ['HLT_DoubleMu7_v1',              #start run 160404
                    'HLT_DoubleMu7_v2',       
                    'HLT_Mu13_Mu8_v1',               #start run 165088 
                    'HLT_Mu13_Mu8_v2',
                    'HLT_Mu13_Mu8_v3',  
                    'HLT_Mu13_Mu8_v4',
                    'HLT_Mu13_Mu8_v5',
                    'HLT_Mu13_Mu8_v6',
                    'HLT_Mu13_Mu8_v7',
                    'HLT_Mu13_Mu8_v8',
                    'HLT_Mu13_Mu8_v9',
                    'HLT_Mu13_Mu8_v10',
                    'HLT_Mu13_Mu8_v11',
                    'HLT_Mu17_Mu8_v1',
                    'HLT_Mu17_Mu8_v2',
                    'HLT_Mu17_Mu8_v3',
                    'HLT_Mu17_Mu8_v4',
                    'HLT_Mu17_Mu8_v5',
                    'HLT_Mu17_Mu8_v6',
                    'HLT_Mu17_Mu8_v7',
                    'HLT_Mu17_TkMu8_v1',
                    'HLT_Mu17_TkMu8_v2',
                    'HLT_Mu17_TkMu8_v3',
                    'HLT_Mu17_TkMu8_v4'
                    ]
        MuEG = ['HLT_Mu17_Ele8_CaloIdL_v1',  #start run 160404
                'HLT_Mu17_Ele8_CaloIdL_v2',
                'HLT_Mu17_Ele8_CaloIdL_v3',   
                'HLT_Mu17_Ele8_CaloIdL_v4',
                'HLT_Mu17_Ele8_CaloIdL_v5',
                'HLT_Mu17_Ele8_CaloIdL_v6',
                'HLT_Mu17_Ele8_CaloIdL_v7',
                'HLT_Mu17_Ele8_CaloIdL_v8',
                'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v1', #start run 167039
                'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v2',
                'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v3', 
                'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4',
                'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v5',
                'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v6',
                'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7',
                'HLT_Mu8_Ele17_CaloIdL_v1',   #start run 160404
                'HLT_Mu8_Ele17_CaloIdL_v2',
                'HLT_Mu8_Ele17_CaloIdL_v3',     
                'HLT_Mu8_Ele17_CaloIdL_v4',
                'HLT_Mu8_Ele17_CaloIdL_v5',
                'HLT_Mu8_Ele17_CaloIdL_v6',
                'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v1', #start run 167039
                'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v2',
                'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3', 
                'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4',
                'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v5',
                'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v6',
                'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7'
                ]
        SingleMu = ['HLT_IsoMu17_v1',  #start run 160404
                    'HLT_IsoMu17_v2',
                    'HLT_IsoMu17_v3',  
                    'HLT_IsoMu17_v4',
                    'HLT_IsoMu17_v5',
                    'HLT_IsoMu17_v6',
                    'HLT_IsoMu17_v7',
                    'HLT_IsoMu17_v8',
                    'HLT_IsoMu17_v9',
                    'HLT_IsoMu17_v10',
                    'HLT_IsoMu17_v11',
                    'HLT_IsoMu17_v12',
                    'HLT_IsoMu17_v13',
                    'HLT_IsoMu24_v1',   #start run 160404
                    'HLT_IsoMu24_v2',
                    'HLT_IsoMu24_v3',  
                    'HLT_IsoMu24_v4',
                    'HLT_IsoMu24_v5',
                    'HLT_IsoMu24_v6',
                    'HLT_IsoMu24_v7', 
                    'HLT_IsoMu24_v8',
                    'HLT_IsoMu30_eta2p1_v1',   #start run 173236
                    'HLT_IsoMu30_eta2p1_v2',
                    'HLT_IsoMu30_eta2p1_v3',
                    'HLT_IsoMu30_eta2p1_v4',
                    'HLT_IsoMu30_eta2p1_v5',
                    'HLT_IsoMu30_eta2p1_v6',
                    'HLT_IsoMu30_eta2p1_v7',
                    'HLT_IsoMu34_eta2p1_v1',   #start run 173236
                    'HLT_IsoMu34_eta2p1_v2',
                    'HLT_IsoMu34_eta2p1_v3',
                    'HLT_IsoMu34_eta2p1_v4',
                    'HLT_IsoMu34_eta2p1_v5'
                    ]

        Photon=['HLT_Photon20_CaloIdVL_IsoL_v1','HLT_Photon20_CaloIdVL_IsoL_v2','HLT_Photon20_CaloIdVL_IsoL_v3','HLT_Photon20_CaloIdVL_IsoL_v4','HLT_Photon20_CaloIdVL_IsoL_v5','HLT_Photon20_CaloIdVL_IsoL_v6','HLT_Photon20_CaloIdVL_IsoL_v7','HLT_Photon20_CaloIdVL_IsoL_v8','HLT_Photon20_CaloIdVL_IsoL_v9',
                'HLT_Photon30_CaloIdVL_IsoL_v1','HLT_Photon30_CaloIdVL_IsoL_v2','HLT_Photon30_CaloIdVL_IsoL_v3','HLT_Photon30_CaloIdVL_IsoL_v4','HLT_Photon30_CaloIdVL_IsoL_v5','HLT_Photon30_CaloIdVL_IsoL_v6','HLT_Photon30_CaloIdVL_IsoL_v7','HLT_Photon30_CaloIdVL_IsoL_v8','HLT_Photon30_CaloIdVL_IsoL_v9','HLT_Photon30_CaloIdVL_IsoL_v10',
                'HLT_Photon50_CaloIdVL_IsoL_v1','HLT_Photon50_CaloIdVL_IsoL_v2','HLT_Photon50_CaloIdVL_IsoL_v3','HLT_Photon50_CaloIdVL_IsoL_v4','HLT_Photon50_CaloIdVL_IsoL_v5','HLT_Photon50_CaloIdVL_IsoL_v6','HLT_Photon50_CaloIdVL_IsoL_v7','HLT_Photon50_CaloIdVL_IsoL_v8','HLT_Photon50_CaloIdVL_IsoL_v9',
                'HLT_Photon75_CaloIdVL_IsoL_v1','HLT_Photon75_CaloIdVL_IsoL_v2','HLT_Photon75_CaloIdVL_IsoL_v3','HLT_Photon75_CaloIdVL_IsoL_v4','HLT_Photon75_CaloIdVL_IsoL_v5','HLT_Photon75_CaloIdVL_IsoL_v6','HLT_Photon75_CaloIdVL_IsoL_v7','HLT_Photon75_CaloIdVL_IsoL_v8','HLT_Photon75_CaloIdVL_IsoL_v9',
                'HLT_Photon90_CaloIdVL_IsoL_v1','HLT_Photon90_CaloIdVL_IsoL_v2','HLT_Photon90_CaloIdVL_IsoL_v3','HLT_Photon90_CaloIdVL_IsoL_v4','HLT_Photon90_CaloIdVL_IsoL_v5','HLT_Photon90_CaloIdVL_IsoL_v6','HLT_Photon90_CaloIdVL_IsoL_v7','HLT_Photon90_CaloIdVL_IsoL_v8','HLT_Photon90_CaloIdVL_IsoL_v9',
                'HLT_Photon125_v1','HLT_Photon125_v2',
                'HLT_Photon125_NoSpikeFilter_v1','HLT_Photon125_NoSpikeFilter_v2','HLT_Photon125_NoSpikeFilter_v3',
                'HLT_Photon135_v1','HLT_Photon135_v2','HLT_Photon135_v3',
                'HLT_Photon200_NoHE_v1','HLT_Photon200_NoHE_v2','HLT_Photon200_NoHE_v3','HLT_Photon200_NoHE_v4'  
                ]
    
    return DoubleElectron, DoubleMu, MuEG, Photon, SingleMu, mcTrigs

