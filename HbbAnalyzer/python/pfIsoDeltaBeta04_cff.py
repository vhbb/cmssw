### assembled and modified from CommonTools/ ParticleFlow/ python/ ParticleSelectors/ pfCandsForIsolation_cff.py
import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.pfNoPileUp_cff  import *
from CommonTools.ParticleFlow.ParticleSelectors.pfSortByType_cff import *

#Create the PU candidates
pfPileUpCandidates = cms.EDProducer("TPPFCandidatesOnPFCandidates",
                                                                             enable =  cms.bool( True ),
                                                                             verbose = cms.untracked.bool( False ),
                                                                             name = cms.untracked.string("pileUpCandidates"),
                                                                             topCollection = cms.InputTag("pfNoPileUp"),
                                                                            ### bottomCollection = cms.InputTag("particleFlowTmp") ## ..Tmp in 44X
                                                                             bottomCollection = cms.InputTag("particleFlow"),
                                    
                                                                             )


 #Take the PU charged particles
pfPUChargedCandidates = cms.EDFilter("PdgIdPFCandidateSelector",
                                                                     src = cms.InputTag("pfPileUpCandidates"),
                                                                     pdgId = cms.vint32(211,-211,321,-321,999211,2212,-2212,11,-11,13,-13)
                                                                     )


#Create All Charged Particles
pfAllChargedCandidates = cms.EDFilter("PdgIdPFCandidateSelector",
       src = cms.InputTag("pfNoPileUp"),
       pdgId = cms.vint32(211,-211,321,-321,999211,2212,-2212,11,-11,13,-13)
       )


pfPileUpCandidatesSequence = cms.Sequence(pfPileUpCandidates+
                pfPUChargedCandidates+
                pfAllChargedCandidates)



pfCandsForIsolationSequence = cms.Sequence(
     pfNoPileUpSequence +
     pfSortByTypeSequence +
     pfPileUpCandidatesSequence
       )

##assembled and modified from /CMSSW/RecoMuon/MuonIsolation/python/muonPFIsolation_cff.py from V02-03-01 (44X) and adapting also for electrons

import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.Isolation.tools_cfi import *
from PhysicsTools.PatAlgos.patSequences_cff import patMuons
from PhysicsTools.PatAlgos.patSequences_cff import patElectrons


##### MUONS ##################
#Now prepare the iso deposits
muPFIsoDepositCharged=isoDepositReplace('pfSelectedMuons','pfAllChargedHadrons')
muPFIsoDepositChargedAll=isoDepositReplace('pfSelectedMuons','pfAllChargedCandidates')
muPFIsoDepositNeutral=isoDepositReplace('pfSelectedMuons','pfAllNeutralHadrons')
muPFIsoDepositGamma=isoDepositReplace('pfSelectedMuons','pfAllPhotons')
muPFIsoDepositPU=isoDepositReplace('pfSelectedMuons','pfPUChargedCandidates')



#Now create isolation values for those isolation deposits 

muPFIsoValueCharged03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("muPFIsoDepositCharged"),
            deltaR = cms.double(0.3),
            weight = cms.string('1'),
            vetos = cms.vstring('0.0001','Threshold(0.0)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
            )
     )
)

muPFIsoValueChargedAll03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("muPFIsoDepositChargedAll"),
            deltaR = cms.double(0.3),
            weight = cms.string('1'),
            vetos = cms.vstring('0.0001','Threshold(0.0)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
     )
   )
)

muPFIsoValueGamma03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("muPFIsoDepositGamma"),
            deltaR = cms.double(0.3),
            weight = cms.string('1'),
            vetos = cms.vstring('0.01','Threshold(0.5)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
      )
   )
)

muPFIsoValueNeutral03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("muPFIsoDepositNeutral"),
            deltaR = cms.double(0.3),
            weight = cms.string('1'),
            vetos = cms.vstring('0.01','Threshold(0.5)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
    )
 )
)

muPFIsoValuePU03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("muPFIsoDepositPU"),
            deltaR = cms.double(0.3),
            weight = cms.string('1'),
            vetos = cms.vstring('0.01','Threshold(0.5)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
      )
   )
)



muPFIsoValueCharged04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("muPFIsoDepositCharged"),
            deltaR = cms.double(0.4),
            weight = cms.string('1'),
            vetos = cms.vstring('0.0001','Threshold(0.0)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
            )
     )
)




muPFIsoValueChargedAll04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("muPFIsoDepositChargedAll"),
            deltaR = cms.double(0.4),
            weight = cms.string('1'),
            vetos = cms.vstring('0.0001','Threshold(0.0)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
     )
   )
)

muPFIsoValueGamma04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("muPFIsoDepositGamma"),
            deltaR = cms.double(0.4),
            weight = cms.string('1'),
            vetos = cms.vstring('0.01','Threshold(0.5)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
      )
   )
)


muPFIsoValueNeutral04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("muPFIsoDepositNeutral"),
            deltaR = cms.double(0.4),
            weight = cms.string('1'),
            vetos = cms.vstring('0.01','Threshold(0.5)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
    )
 )

)
muPFIsoValuePU04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("muPFIsoDepositPU"),
            deltaR = cms.double(0.4),
            weight = cms.string('1'),
            vetos = cms.vstring('0.01','Threshold(0.5)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
      )
   )
)

patMuons.isoDeposits = cms.PSet(
       pfAllParticles    = cms.InputTag("muPFIsoDepositChargedAll"),
       pfChargedHadrons = cms.InputTag("muPFIsoDepositCharged"),
       pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutral"),
       pfPhotons        = cms.InputTag("muPFIsoDepositGamma")
)

patMuons.isolationValues = cms.PSet(
           pfAllParticles   = cms.InputTag("muPFIsoValueChargedAll04"),
           pfChargedHadrons = cms.InputTag("muPFIsoValueCharged04"),
           pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral04"),
           pfPhotons        = cms.InputTag("muPFIsoValueGamma04"),
           user = cms.VInputTag(
                        cms.InputTag("muPFIsoValueChargedAll03"),            
                        cms.InputTag("muPFIsoValueCharged03"),               
                        cms.InputTag("muPFIsoValueNeutral03"),               
                        cms.InputTag("muPFIsoValueGamma03"),                 
                        cms.InputTag("muPFIsoValuePU03"),
                        cms.InputTag("muPFIsoValuePU04"),
                        
          )
 )





muonPFIsolationSequence =  cms.Sequence(muPFIsoDepositCharged+
                                        muPFIsoDepositChargedAll+
                                        muPFIsoDepositGamma+
                                        muPFIsoDepositNeutral+
                                        muPFIsoDepositPU+
                                        ##############################
                                        muPFIsoValueCharged03+
                                        muPFIsoValueChargedAll03+
                                        muPFIsoValueGamma03+
                                        muPFIsoValueNeutral03+
                                        muPFIsoValuePU03+
                                        ############################## 
                                        muPFIsoValueCharged04+
                                        muPFIsoValueChargedAll04+
                                        muPFIsoValueGamma04+
                                        muPFIsoValueNeutral04+
                                        muPFIsoValuePU04 
                                       ## *     patMuons    
                                      

)                                         





########  ELECTRONS ###########
#Now prepare the iso deposits
elePFIsoDepositCharged=isoDepositReplace('pfSelectedElectrons','pfAllChargedHadrons')
elePFIsoDepositChargedAll=isoDepositReplace('pfSelectedElectrons','pfAllChargedCandidates')
elePFIsoDepositNeutral=isoDepositReplace('pfSelectedElectrons','pfAllNeutralHadrons')
elePFIsoDepositGamma=isoDepositReplace('pfSelectedElectrons','pfAllPhotons')
elePFIsoDepositPU=isoDepositReplace('pfSelectedElectrons','pfPUChargedCandidates')



#Now create isolation values for those isolation deposits 


elePFIsoValueCharged03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("elePFIsoDepositCharged"),
            deltaR = cms.double(0.3),
            weight = cms.string('1'),
            vetos = cms.vstring('0.0001','Threshold(0.0)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
            )
     )
)

elePFIsoValueChargedAll03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("elePFIsoDepositChargedAll"),
            deltaR = cms.double(0.3),
            weight = cms.string('1'),
            vetos = cms.vstring('0.0001','Threshold(0.0)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
     )
   )
)

elePFIsoValueGamma03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("elePFIsoDepositGamma"),
            deltaR = cms.double(0.3),
            weight = cms.string('1'),
            vetos = cms.vstring('0.01','Threshold(0.5)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
      )
   )
)

elePFIsoValueNeutral03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("elePFIsoDepositNeutral"),
            deltaR = cms.double(0.3),
            weight = cms.string('1'),
            vetos = cms.vstring('0.01','Threshold(0.5)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
    )
 )
)

elePFIsoValuePU03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("elePFIsoDepositPU"),
            deltaR = cms.double(0.3),
            weight = cms.string('1'),
            vetos = cms.vstring('0.01','Threshold(0.5)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
      )
   )
)


elePFIsoValueCharged04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("elePFIsoDepositCharged"),
            deltaR = cms.double(0.4),
            weight = cms.string('1'),
            vetos = cms.vstring('0.0001','Threshold(0.0)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
            )
     )
)




elePFIsoValueChargedAll04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("elePFIsoDepositChargedAll"),
            deltaR = cms.double(0.4),
            weight = cms.string('1'),
            vetos = cms.vstring('0.0001','Threshold(0.0)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
     )
   )
)

elePFIsoValueGamma04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("elePFIsoDepositGamma"),
            deltaR = cms.double(0.4),
            weight = cms.string('1'),
            vetos = cms.vstring('0.01','Threshold(0.5)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
      )
   )
)


elePFIsoValueNeutral04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("elePFIsoDepositNeutral"),
            deltaR = cms.double(0.4),
            weight = cms.string('1'),
            vetos = cms.vstring('0.01','Threshold(0.5)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
    )
 )

)
elePFIsoValuePU04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("elePFIsoDepositPU"),
            deltaR = cms.double(0.4),
            weight = cms.string('1'),
            vetos = cms.vstring('0.01','Threshold(0.5)'),
            skipDefaultVeto = cms.bool(True),
            mode = cms.string('sum')
      )
   )
)

patElectrons.isoDeposits = cms.PSet(
       pfAllParticles    = cms.InputTag("elePFIsoDepositChargedAll"),
       pfChargedHadrons = cms.InputTag("elePFIsoDepositCharged"),
       pfNeutralHadrons = cms.InputTag("elePFIsoDepositNeutral"),
       pfPhotons        = cms.InputTag("elePFIsoDepositGamma")
)

patElectrons.isolationValues = cms.PSet(
           pfAllParticles   = cms.InputTag("elePFIsoValueChargedAll04"),
           pfChargedHadrons = cms.InputTag("elePFIsoValueCharged04"),
           pfNeutralHadrons = cms.InputTag("elePFIsoValueNeutral04"),
           pfPhotons        = cms.InputTag("elePFIsoValueGamma04"),
           user = cms.VInputTag(
                        cms.InputTag("elePFIsoValueChargedAll03"),            
                        cms.InputTag("elePFIsoValueCharged03"),               
                        cms.InputTag("elePFIsoValueNeutral03"),               
                        cms.InputTag("elePFIsoValueGamma03"),                 
                        cms.InputTag("elePFIsoValuePU03"),
                        cms.InputTag("elePFIsoValuePU04"),
                        
          )
 )



electronPFIsolationSequence =  cms.Sequence(elePFIsoDepositCharged+
                                        elePFIsoDepositChargedAll+
                                        elePFIsoDepositGamma+
                                        elePFIsoDepositNeutral+
                                        elePFIsoDepositPU+
                                        ##############################
                                        elePFIsoValueCharged03+
                                        elePFIsoValueChargedAll03+
                                        elePFIsoValueGamma03+
                                        elePFIsoValueNeutral03+
                                        elePFIsoValuePU03+
                                        ############################## 
                                        elePFIsoValueCharged04+
                                        elePFIsoValueChargedAll04+
                                        elePFIsoValueGamma04+
                                        elePFIsoValueNeutral04+
                                        elePFIsoValuePU04 
                                       ## patElectrons
                                       

)                                         
