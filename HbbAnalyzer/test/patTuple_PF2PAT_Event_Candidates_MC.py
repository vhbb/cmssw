#import configurations
import FWCore.ParameterSet.Config as cms

import os 

isMC = True

# define the process
process = cms.Process("VH")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
SkipEvent = cms.untracked.vstring('ProductNotFound') )

## import skeleton process
#from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# source
# on lxbuild151
process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("file:/tmp/tboccali/trigger.root"))
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("file:/build1/tboccali/7C74874C-CA8E-E011-9782-001D09F25401.root"))
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/relval/CMSSW_4_2_3/RelValTTbar/GEN-SIM-RECO/MC_42_V12-v2/0066/3026A5BD-D97B-E011-A9D7-001A92811736.root"))


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if isMC == False :
	process.GlobalTag.globaltag = cms.string('GR_R_42_V13::All')
else :
	process.GlobalTag.globaltag = cms.string('START42_V12::All')

process.out1 = cms.OutputModule(
    'PoolOutputModule',
    fileName       = cms.untracked.string('PAT.edm.root'),
    outputCommands = cms.untracked.vstring(
	'drop *',
					   'keep *_HbbAnalyzerNew_*_*',
					   'keep *_hbbCandidates_*_*',
					   'keep PileupSummaryInfo_*_*_*',
					   'keep edmTriggerResults_*_*_*',
					   ),
    dropMetaData = cms.untracked.string('ALL'),
    splitLevel = cms.untracked.int32(0)
    )


process.out = cms.OutputModule(
    'PoolOutputModule',
    fileName       = cms.untracked.string('fanculo.edm.root'),
    outputCommands = cms.untracked.vstring('drop *','keep *_HbbAnalyzerNew_*_*')
    )


# load the standard PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# note that you can use a bunch of core tools of PAT 
# to taylor your PAT configuration; for a few examples
# uncomment the following lines

from PhysicsTools.PatAlgos.tools.coreTools import *
#removeMCMatching(process, 'Muons')
#removeAllPATObjectsBut(process, ['Muons'])
#removeSpecificPATObjects(process, ['Electrons', 'Muons', 'Taus'])

if isMC == False :
	removeMCMatching(process, ['All'])
        RunOnData()
        process.makePatMuons.remove(process.muonMatch)
        process.makePatTaus.remove(process.tauMatch)
        process.makePatTaus.remove(process.tauGenJets)
        process.makePatTaus.remove(process.tauGenJetsSelectorAllHadrons)
        process.makePatTaus.remove(process.tauGenJetMatch)
        process.makePatPhotons.remove(process.photonMatch)
        process.makePatJets.remove(process.patJetPartonMatch)
        process.makePatJets.remove(process.patJetGenJetMatch)
        process.makePatJets.remove(process.patJetFlavourId)

# add the trigger information to the configuration
from PhysicsTools.PatAlgos.tools.trigTools import *
#####switchOnTrigger( process )
from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent

## want pfNoElectron

process.load('CommonTools.ParticleFlow.PF2PAT_cff')
from PhysicsTools.PatAlgos.tools.pfTools import *

if isMC == False :
	usePF2PAT(process,runPF2PAT=False,jetAlgo='AK5',runOnMC=False,postfix='')
else :
	usePF2PAT(process,runPF2PAT=False,jetAlgo='AK5',runOnMC=True,postfix='')

process.pfPileUp.Enable = True
process.pfIsolatedElectrons.combinedIsolationCut = 0.35
process.pfIsolatedMuons.combinedIsolationCut     = 0.35

##bbb  study
#####process.pfNoMuon.enable  = False
#####process.pfIsolatedElectrons.combinedIsolationCut = 9999.
#####process.pfIsolatedMuons.combinedIsolationCut     = 9999.
######

### suggested recipe
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(7.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )
process.pfJets.Vertices = cms.InputTag('goodOfflinePrimaryVertices')
process.pfPileUp.checkClosestZVertex = cms.bool(False)
process.pfJets.doAreaFastjet = True
process.pfJets.doRhoFastjet = False
process.patJetCorrFactors.rho = cms.InputTag("kt6PFJets", "rho")
############



#! PU DENSITY RHO + L1 FastJET to AK5PF
#!
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJets = kt4PFJets.clone(
    rParam = cms.double(0.6),
    src = cms.InputTag('pfNoElectron'),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True),
    voronoiRfact = cms.double(0.9)
    )

## 4.1.4.
#process.load('RecoJets.JetProducers.ak5PFJets_cfi')
#process.ak5PFJets.jetPtMin      = 1.0
#process.ak5PFJets.doAreaFastjet = True
#process.ak5PFJets.src           = 'pfNoElectron'

process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')
#process.ak5PFL1Fastjet.useCondDB = False
#process.ak5PFL1Fastjet.srcRho = cms.InputTag('kt6PFJets','rho')
#process.ak5PFJetsL1 = cms.EDProducer(
#    'PFJetCorrectionProducer',
#    src        = cms.InputTag('ak5PFJets'),
#    correctors = cms.vstring('ak5PFL1Fastjet')
#    )
####

from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
process.ak5PFJets = ak5PFJets
process.ak5PFJets.src           = 'pfNoElectron'
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.jetPtMin      = 1.0


#######

process.ak7PFJets = process.ak5PFJets.clone( rParam = cms.double(0.7) )

from RecoJets.JetProducers.ak5CaloJets_cfi import ak5CaloJets
process.ak5CaloJets = ak5CaloJets
process.ak7CaloJets = ak5CaloJets.clone( rParam = cms.double(0.7) )

from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
process.ak5GenJets = ak5GenJets


# CATopJets
process.load("RecoJets.Configuration.GenJetParticles_cff")
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.ca4GenJets = ca4GenJets
process.load("VHbbAnalysis.HbbAnalyzer.caTopJets_cff")
##process.load("TopQuarkAnalysis.TopPairBSM.CATopJetTagger_cfi")

process.ca8GenJets = ca4GenJets.clone( rParam = cms.double(0.8) )

## tune DeltaR for C/A

process.caTopCaloJets.rBins = cms.vdouble(1.2,1.2,1.2)
process.caTopPFJets.rBins = cms.vdouble(1.2,1.2,1.2)

process.caTopCaloJets.rParam = cms.double(1.2)
process.caTopPFJets.rParam = cms.double(1.2)

#process.caTopCaloJets.rBins = cms.vdouble(1.4,1.4,1.4)
#process.caTopPFJets.rBins = cms.vdouble(1.4,1.4,1.4)
#process.caTopCaloJets.rParam = cms.double(1.4)
#process.caTopPFJets.rParam = cms.double(1.4)

##test
process.caTopCaloJets.useAdjacency = cms.int32(0)
process.caTopPFJets.useAdjacency = cms.int32(0)
######

process.ak5PFJets.src = cms.InputTag('pfNoElectron')
process.ak7PFJets.src = cms.InputTag('pfNoElectron')
process.caTopPFJets.src = cms.InputTag('pfNoElectron')


# subJetProducer
process.CAsubJetsProducer = cms.EDProducer('SubProducer',
jetTag = cms.untracked.InputTag("caTopPFJets"))


# switch jet collection to our juets
from PhysicsTools.PatAlgos.tools.jetTools import *


## ==== Example with CaloJets
addJetCollection(process, 
        cms.InputTag('caTopCaloJets'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
        "CA","Calo",
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
# temp off        jetCorrLabel=('KT6Calo', ['L2Relative', 'L3Absolute', 'L2L3Residual']),
##        jetCorrLabel=('KT6Calo', ['L2Relative', 'L3Absolute']),
        jetCorrLabel=('KT4Calo', ['L2Relative', 'L3Absolute']),
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
#        genJetCollection = cms.InputTag("ca8GenJets"),
        doJetID=False
                 )

## ==== Example withPFJets
addJetCollection(process, 
        cms.InputTag('caTopPFJets'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
        'CA',"PF",
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
#        jetCorrLabel=('KT6PF',  ['L2Relative', 'L3Absolute']),   # example jet correction name; set to None for no JEC
# temp off        jetCorrLabel=('KT6PF',  ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']),
####        jetCorrLabel=('KT6PF', ['L1FastJet', 'L2Relative', 'L3Absolute']),
              jetCorrLabel=('KT4PF', ['L1FastJet','L2Relative', 'L3Absolute']),
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
#        genJetCollection = cms.InputTag("ca8GenJets"),
        doJetID = False
                 )

addJetCollection(process,cms.InputTag("CAsubJetsProducer"),"subCA","PF",
##                 doJTA=True,doBTagging=True,jetCorrLabel=('KT4PF',  ['L2Relative', 'L3Absolute']),
                 doJTA=True,doBTagging=True,
#temp off                 jetCorrLabel=('KT4PF',  ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']),
              jetCorrLabel=('KT4PF', ['L1FastJet','L2Relative', 'L3Absolute']),
##                 doJTA=True,doBTagging=True,jetCorrLabel=None,
                 doType1MET=False,doL1Cleaning = False,doL1Counters=False,
#genJetCollection = cms.InputTag("ca4GenJets"),
                 doJetID = False)

## add standard jets for di-jet analysis type



addJetCollection(process,cms.InputTag("ak7CaloJets"),"AK7","Calo",
                 doJTA=True,doBTagging=True,jetCorrLabel=('AK7Calo',  ['L2Relative', 'L3Absolute']),
                 doType1MET=False,doL1Cleaning = False,doL1Counters=False,
                 doJetID = False)

addJetCollection(process,cms.InputTag("ak5PFJets"),"AK5","PF",
#addJetCollection(process,cms.InputTag("ak5PFJetsL1"),"AK5","PF",
                 doJTA=True,doBTagging=True,
#temp off                 jetCorrLabel=('AK5PFchs',  ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']),
                 jetCorrLabel=('AK5PFchs', ['L1FastJet','L2Relative', 'L3Absolute']),
                 doType1MET=False,doL1Cleaning = False,doL1Counters=False,
                 doJetID = False)

addJetCollection(process,cms.InputTag("ak7PFJets"),"AK7","PF",
                 doJTA=True,doBTagging=True,jetCorrLabel=('AK7PF',  ['L2Relative', 'L3Absolute']),
                 doType1MET=False,doL1Cleaning = False,doL1Counters=False,
                 doJetID = False)


switchJetCollection(process, 
        cms.InputTag('ak5CaloJets'),     # Jet collection; must be already in the event when patLayer0 sequence is executed
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
#temp off        jetCorrLabel=('AK5Calo',  ['L2Relative', 'L3Absolute', 'L2L3Residual']), 
        jetCorrLabel=('AK5Calo', ['L2Relative', 'L3Absolute']),
        doType1MET=False,
#        genJetCollection = cms.InputTag("ak5GenJets"),
        doJetID = False,
        jetIdLabel = 'ak5'
                    )

from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

# Jets
#setTagInfos(process,
#            coll = "allLayer1Jets",
#            tagInfos = cms.vstring("secondaryVertex")
#            )

for jetcoll in (process.selectedPatJetsCACalo, process.selectedPatJetsCAPF):
	jetcoll.embedGenJetMatch = cms.bool(False)
	#getJetMCFlavour uses jetFlavourAssociation*, see below
	jetcoll.getJetMCFlavour = cms.bool(True)
	#those two use jetPartonMatch*, see below
	jetcoll.addGenPartonMatch = cms.bool(True)
	jetcoll.embedGenPartonMatch = cms.bool(True)
	# Add CATopTag info... piggy-backing on b-tag functionality
	jetcoll.discriminatorSources = cms.VInputTag()
	jetcoll.addBTagInfo = cms.bool(True)
	jetcoll.addTagInfos = cms.bool(True)
	jetcoll.addDiscriminators = cms.bool(False)

##process.selectedPatJetsCACalo.tagInfoSources = cms.VInputTag( cms.InputTag('CATopCaloJetTagInfos') )
##process.selectedPatJetsCAPF.tagInfoSources = cms.VInputTag( cms.InputTag('CATopPFJetTagInfos') )

# jet flavor stuff
#process.jetPartons.withTop = cms.bool(True)
#process.jetPartonAssociation.coneSizeToAssociate = cms.double(1.2)
#process.jetPartonAssociation.doPriority = cms.bool(True)
#process.jetPartonAssociation.priorityList = cms.vint32(6)
#process.jetFlavourAssociation.definition = cms.int32(4)
#process.jetFlavourAssociation.physicsDefinition = cms.bool(False)
#process.jetPartonMatch.mcPdgId = cms.vint32(1,2,3,4,5,6,21)
#process.jetPartonMatch.maxDeltaR = cms.double(1.2)

# Place appropriate jet cuts (NB: no cut on number of constituents)
process.selectedPatJets.cut = cms.string('pt > 15. & abs(eta) < 5.0')
process.selectedPatJetsAK5PF.cut = cms.string('pt > 15. & abs(eta) < 5.0')
process.selectedPatJetsCACalo.cut = cms.string('pt > 0. & abs(eta) < 5.0')
process.selectedPatJetsCAPF.cut = cms.string('pt > 0. & abs(eta) < 5.0')
##process.selectedLayer1Muons.cut = cms.string('pt > 20. & abs(eta) < 2.5 & muonID("TMLastStationLoose")')
##process.selectedLayer1Electrons.cut = cms.string('pt > 20. & abs(eta) < 2.5 & electronID("eidLoose")')
process.selectedPatMuons.cut = cms.string('isGlobalMuon || (isTrackerMuon && muonID("TrackerMuonArbitrated"))')
process.selectedPatElectrons.cut = cms.string('')

#obsolete
algorithm = 'ca'

if algorithm == 'kt' :
    process.caTopCaloJets.algorithm = cms.int32(0)
    process.caTopGenJets.algorithm = cms.int32(0)
    process.caTopPFJets.algorithm = cms.int32(0)
elif algorithm == 'ca' :
    process.caTopCaloJets.algorithm = cms.int32(1)
    process.caTopGenJets.algorithm = cms.int32(1)
    process.caTopPFJets.algorithm = cms.int32(1)    
elif algorithm == 'antikt' :
    process.caTopCaloJets.algorithm = cms.int32(2)
    process.caTopGenJets.algorithm = cms.int32(2)
    process.caTopPFJets.algorithm = cms.int32(2)
else:
    print "Error: algorithm '%s' unknown. Use one of kt, ca, antikt." % algorithm
    raise AttributeError()
####



# pythia output
#process.printList = cms.EDAnalyzer( "ParticleListDrawer",
#                                src = cms.InputTag( "genParticles" ),
#                                maxEventsToPrint = cms.untracked.int32( 1 )
#)


# prune gen particles

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")


from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector 
process.goodPatJetsAK5PF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                           filterParams = pfJetIDSelector.clone(), src = cms.InputTag("selectedPatJetsAK5PF") )

process.HbbAnalyzerNew = cms.EDProducer("HbbAnalyzerNew",
    runOnMC = cms.bool(isMC),
    electronTag = cms.InputTag("selectedPatElectrons"),
    hltResultsTag = cms.InputTag("TriggerResults::HLT1"),
#    tauTag = cms.InputTag("cleanPatTaus"),
    tauTag = cms.InputTag("patTaus"),
    muonTag = cms.InputTag("selectedPatMuons"),
    jetTag = cms.InputTag("selectedPatJetsCAPF"),
    subjetTag = cms.InputTag("selectedPatJetssubCAPF"),
##    jetTag = cms.InputTag("selectedPatJetsCACalo"),
    simplejet1Tag = cms.InputTag("selectedPatJets"),
###					    simplejet1Tag = cms.InputTag("selectedPatJetsAK7PF"),
    simplejet2Tag = cms.InputTag("goodPatJetsAK5PF"),
    simplejet3Tag = cms.InputTag("selectedPatJetsAK7Calo"),
    simplejet4Tag = cms.InputTag("selectedPatJetsAK7PF"),
    photonTag = cms.InputTag("selectedPatPhotons"),
    metTag = cms.InputTag("patMETs"),
    dimuTag = cms.InputTag("dimuons"),
    dielecTag = cms.InputTag("dielectrons"),
    verbose = cms.untracked.bool(True)
)

process.dimuons = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(False),
    cut = cms.string('mass > 40 && mass < 200'),
    decay = cms.string('selectedPatMuons@+ selectedPatMuons@-')
)

process.dielectrons = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(False),
    cut = cms.string('mass > 40 && mass < 200'),
    decay = cms.string('selectedPatElectrons@+ selectedPatElectrons@-')
)

#process.Tracer = cms.Service("Tracer")


process.load("VHbbAnalysis.HbbAnalyzer.simpleEleIdSequence_cff")
process.patElectronIDs = cms.Sequence(process.simpleEleIdSequence)

if isMC == False :
	process.makePatElectrons = cms.Sequence(process.patElectronIDs*process.patElectronIsolation*process.patElectrons)
else :
	process.makePatElectrons = cms.Sequence(process.patElectronIDs*process.patElectronIsolation*process.electronMatch*process.patElectrons)

process.patElectrons.electronIDSources = cms.PSet(
    simpleEleId95relIso= cms.InputTag("simpleEleId95relIso"),
    simpleEleId90relIso= cms.InputTag("simpleEleId90relIso"),
    simpleEleId85relIso= cms.InputTag("simpleEleId85relIso"),
    simpleEleId80relIso= cms.InputTag("simpleEleId80relIso"),
    simpleEleId70relIso= cms.InputTag("simpleEleId70relIso"),
    simpleEleId60relIso= cms.InputTag("simpleEleId60relIso"),
    simpleEleId95cIso= cms.InputTag("simpleEleId95cIso"),
    simpleEleId90cIso= cms.InputTag("simpleEleId90cIso"),
    simpleEleId85cIso= cms.InputTag("simpleEleId85cIso"),
    simpleEleId80cIso= cms.InputTag("simpleEleId80cIso"),
    simpleEleId70cIso= cms.InputTag("simpleEleId70cIso"),
    simpleEleId60cIso= cms.InputTag("simpleEleId60cIso"),
    eidTight = cms.InputTag("eidTight"),
    eidLoose = cms.InputTag("eidLoose"),
    eidRobustTight = cms.InputTag("eidRobustTight"),
    eidRobustHighEnergy = cms.InputTag("eidRobustHighEnergy"),
    eidRobustLoose = cms.InputTag("eidRobustLoose")
    )



process.hbbCandidates = cms.EDProducer("HbbCandidateFinder",
				       VHbbEventLabel = cms.InputTag(""),
				       verbose = cms.bool(True) ,
				       jetPtThreshold = cms.double(30.)
				      )


# drop the meta data for dropped data
#process.out.dropMetaData = cms.string("DROPPED")

#process.out.fileName = '/tigress-hsm/dlopes/PatEDM.root'


# define path 'p'
process.p = cms.Path(
#process.genJetParticles*
                     process.goodOfflinePrimaryVertices*
                     process.PF2PAT*
                     process.makePatElectrons*
                     process.ak5CaloJets*
                     process.ak7CaloJets*
#                     process.ak5GenJets*
##                     process.ak5PFJets*
######                     process.kt6PFJets*process.ak5PFJets*process.ak5PFJetsL1*
                     process.kt6PFJets*process.ak5PFJets*
                     process.ak7PFJets*
#                     process.ca4GenJets*
#                     process.ca8GenJets*
#                     process.ca12GenJets*
#                     process.caTopGenJets*
                     process.caTopCaloJets*
                     process.caTopPFJets*
                     process.CAsubJetsProducer*
                     process.patDefaultSequence*
                     process.patPF2PATSequence* # added with usePF2PAT
                     process.patCandidates*
                     process.dimuons*
                     process.dielectrons*
                     process.goodPatJetsAK5PF*
                     process.HbbAnalyzerNew*
		     process.hbbCandidates
                     )


process.e = cms.EndPath(process.out1)

#
