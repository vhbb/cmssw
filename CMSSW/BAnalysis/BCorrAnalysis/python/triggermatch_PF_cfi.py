import FWCore.ParameterSet.Config as cms

selectedJetTriggerMatchHLTL1Jet10U = cms.EDProducer(
            "PATTriggerMatcherDRDPtLessByR"                    # match by DeltaR only, best match by DeltaR
                        , src     = cms.InputTag( 'selectedPatJetsAK5PF' )
                        , matched = cms.InputTag( 'patTrigger' )          # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
                        , andOr                      = cms.bool( False )  # AND
                        , filterIdsEnum              = cms.vstring( '*' ) # wildcard, overlaps with 'filterIds'
                        , filterIds                  = cms.vint32( 0 )    # wildcard, overlaps with 'filterIdsEnum'
                        , filterLabels               = cms.vstring( '*' ) # wildcard
                        , pathNames                  = cms.vstring(
            'HLT_L1Jet10U'
                        )
                         , pathLastFilterAcceptedOnly = cms.bool( True )   # select only trigger objects used in last filters of succeeding paths
                         , collectionTags             = cms.vstring( '*' ) # wildcard
                         , maxDPtRel = cms.double( 3.0 )
                         , maxDeltaR = cms.double( 0.4 )
                         , resolveAmbiguities    = cms.bool( True )        # only one match per trigger object
                         , resolveByMatchQuality = cms.bool( True )        # take best match found per reco object: by DeltaR here (s. above)

                        )

selectedJetTriggerMatchHLTL1Jet6U = cms.EDProducer(
                "PATTriggerMatcherDRDPtLessByR"                    # match by DeltaR only, best match by DeltaR
                                        , src     = cms.InputTag( 'selectedPatJetsAK5PF' )
                                        , matched = cms.InputTag( 'patTrigger' )          # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
                                        , andOr                      = cms.bool( False )  # AND
                                        , filterIdsEnum              = cms.vstring( '*' ) # wildcard, overlaps with 'filterIds'
                                        , filterIds                  = cms.vint32( 0 )    # wildcard, overlaps with 'filterIdsEnum'
                                        , filterLabels               = cms.vstring( '*' ) # wildcard
                                        , pathNames                  = cms.vstring(
                'HLT_L1Jet6U'
                                        )
                                         , pathLastFilterAcceptedOnly = cms.bool( True )   # select only trigger objects used in last filters of succeeding paths
                                         , collectionTags             = cms.vstring( '*' ) # wildcard
                                         , maxDPtRel = cms.double( 3.0 )
                                         , maxDeltaR = cms.double( 0.4 )
                                         , resolveAmbiguities    = cms.bool( True )        # only one match per trigger object
                                         , resolveByMatchQuality = cms.bool( True )        # take best match found per reco object: by DeltaR here (s. above)

                                        )


selectedJetTriggerMatchHLTJet15U = cms.EDProducer(
        "PATTriggerMatcherDRDPtLessByR"                    # match by DeltaR only, best match by DeltaR
            , src     = cms.InputTag( 'selectedPatJetsAK5PF' )
            , matched = cms.InputTag( 'patTrigger' )          # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
            , andOr                      = cms.bool( False )  # AND
            , filterIdsEnum              = cms.vstring( '*' ) # wildcard, overlaps with 'filterIds'
            , filterIds                  = cms.vint32( 0 )    # wildcard, overlaps with 'filterIdsEnum'
            , filterLabels               = cms.vstring( '*' ) # wildcard
            , pathNames                  = cms.vstring(
        'HLT_Jet15U'
            )
             , pathLastFilterAcceptedOnly = cms.bool( True )   # select only trigger objects used in last filters of succeeding paths
             , collectionTags             = cms.vstring( '*' ) # wildcard
             , maxDPtRel = cms.double( 3.0 )
             , maxDeltaR = cms.double( 0.4 )
             , resolveAmbiguities    = cms.bool( True )        # only one match per trigger object
             , resolveByMatchQuality = cms.bool( True )        # take best match found per reco object: by DeltaR here (s. above)

            )

selectedJetTriggerMatchHLTJet30U = cms.EDProducer(
        "PATTriggerMatcherDRDPtLessByR"                    # match by DeltaR only, best match by DeltaR
            , src     = cms.InputTag( 'selectedPatJetsAK5PF' )
            , matched = cms.InputTag( 'patTrigger' )          # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
            , andOr                      = cms.bool( False )  # AND
            , filterIdsEnum              = cms.vstring( '*' ) # wildcard, overlaps with 'filterIds'
            , filterIds                  = cms.vint32( 0 )    # wildcard, overlaps with 'filterIdsEnum'
            , filterLabels               = cms.vstring( '*' ) # wildcard
            , pathNames                  = cms.vstring(
        'HLT_Jet30U'
            )
             , pathLastFilterAcceptedOnly = cms.bool( True )   # select only trigger objects used in last filters of succeeding paths
             , collectionTags             = cms.vstring( '*' ) # wildcard
             , maxDPtRel = cms.double( 3.0 )
             , maxDeltaR = cms.double( 0.4 )
             , resolveAmbiguities    = cms.bool( True )        # only one match per trigger object
             , resolveByMatchQuality = cms.bool( True )        # take best match found per reco object: by DeltaR here (s. above)

            )

selectedJetTriggerMatchHLTJet50U = cms.EDProducer(
            "PATTriggerMatcherDRDPtLessByR"                    # match by DeltaR only, best match by DeltaR
                        , src     = cms.InputTag( 'selectedPatJetsAK5PF' )
                        , matched = cms.InputTag( 'patTrigger' )          # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
                        , andOr                      = cms.bool( False )  # AND
                        , filterIdsEnum              = cms.vstring( '*' ) # wildcard, overlaps with 'filterIds'
                        , filterIds                  = cms.vint32( 0 )    # wildcard, overlaps with 'filterIdsEnum'
                        , filterLabels               = cms.vstring( '*' ) # wildcard
                        , pathNames                  = cms.vstring(
            'HLT_Jet50U'
                        )
                         , pathLastFilterAcceptedOnly = cms.bool( True )   # select only trigger objects used in last filters of succeeding paths
                         , collectionTags             = cms.vstring( '*' ) # wildcard
                         , maxDPtRel = cms.double( 3.0 )
                         , maxDeltaR = cms.double( 0.4 )
                         , resolveAmbiguities    = cms.bool( True )        # only one match per trigger object
                         , resolveByMatchQuality = cms.bool( True )        # take best match found per reco object: by DeltaR here (s. above)

                        )

selectedJetTriggerMatchHLTJet70U = cms.EDProducer(
            "PATTriggerMatcherDRDPtLessByR"                    # match by DeltaR only, best match by DeltaR
                        , src     = cms.InputTag( 'selectedPatJetsAK5PF' )
                        , matched = cms.InputTag( 'patTrigger' )          # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
                        , andOr                      = cms.bool( False )  # AND
                        , filterIdsEnum              = cms.vstring( '*' ) # wildcard, overlaps with 'filterIds'
                        , filterIds                  = cms.vint32( 0 )    # wildcard, overlaps with 'filterIdsEnum'
                        , filterLabels               = cms.vstring( '*' ) # wildcard
                        , pathNames                  = cms.vstring(
            'HLT_Jet70U'
                        )
                         , pathLastFilterAcceptedOnly = cms.bool( True )   # select only trigger objects used in last filters of succeeding paths
                         , collectionTags             = cms.vstring( '*' ) # wildcard
                         , maxDPtRel = cms.double( 3.0 )
                         , maxDeltaR = cms.double( 0.4 )
                         , resolveAmbiguities    = cms.bool( True )        # only one match per trigger object
                         , resolveByMatchQuality = cms.bool( True )        # take best match found per reco object: by DeltaR here (s. above)

                        )

selectedJetTriggerMatchHLTJet100U = cms.EDProducer(
            "PATTriggerMatcherDRDPtLessByR"                    # match by DeltaR only, best match by DeltaR
                        , src     = cms.InputTag( 'selectedPatJetsAK5PF' )
                        , matched = cms.InputTag( 'patTrigger' )          # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
                        , andOr                      = cms.bool( False )  # AND
                        , filterIdsEnum              = cms.vstring( '*' ) # wildcard, overlaps with 'filterIds'
                        , filterIds                  = cms.vint32( 0 )    # wildcard, overlaps with 'filterIdsEnum'
                        , filterLabels               = cms.vstring( '*' ) # wildcard
                        , pathNames                  = cms.vstring(
            'HLT_Jet100U'
                        )
                         , pathLastFilterAcceptedOnly = cms.bool( True )   # select only trigger objects used in last filters of succeeding paths
                         , collectionTags             = cms.vstring( '*' ) # wildcard
                         , maxDPtRel = cms.double( 3.0 )
                         , maxDeltaR = cms.double( 0.4 )
                         , resolveAmbiguities    = cms.bool( True )        # only one match per trigger object
                         , resolveByMatchQuality = cms.bool( True )        # take best match found per reco object: by DeltaR here (s. above)

                        )




matchPATPF = cms.Sequence(
    selectedJetTriggerMatchHLTL1Jet10U*
    selectedJetTriggerMatchHLTL1Jet6U*
    selectedJetTriggerMatchHLTJet15U*
    selectedJetTriggerMatchHLTJet30U*
    selectedJetTriggerMatchHLTJet50U*
    selectedJetTriggerMatchHLTJet70U*
    selectedJetTriggerMatchHLTJet100U
    )
