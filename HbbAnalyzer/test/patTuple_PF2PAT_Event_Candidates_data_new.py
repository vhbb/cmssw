
#import configurations
import FWCore.ParameterSet.Config as cms

import os 

isMC = False

# define the process
process = cms.Process("VH")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")


## import skeleton process
#from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# source
process.source = cms.Source("PoolSource",
			    fileNames=cms.untracked.vstring(
		'root://cmsdcache7.pi.infn.it:7070//store/data/Run2011A/SingleMu/AOD/PromptReco-v2/000/163/374/8ABA0558-C470-E011-B566-003048D2C108.root',
#		"/store/mc/CMSSW_4_2_3/RelValProdTTbar/GEN-SIM-RECO/MC_42_V12_JobRobot-v1/0000/B89A0B07-818C-E011-953E-0030487CD7E0.root"
#	"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/mc/Summer11/BdToMuMu_2MuPtFilter_7TeV-pythia6-evtgen//GEN-SIM-RECO//PU_S4_START42_V11-v1///0000//92DE42C9-CD8C-E011-A421-001F296B758E.root"
	)
			    )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if isMC == False :
	process.GlobalTag.globaltag = cms.string('GR_R_42_V19::All')
else :
	process.GlobalTag.globaltag = cms.string('START42_V13::All')

process.out1 = cms.OutputModule(
    'PoolOutputModule',
    fileName       = cms.untracked.string('/tmp/tboccali/PAT.edm.root'),
    outputCommands = cms.untracked.vstring(
	'drop *',
	'keep *_HbbAnalyzerNew_*_*',
	'keep VHbbCandidates_*_*_*',
#	'keep PileupSummaryInfos_*_*_*',
	'keep edmTriggerResults_*_*_*',
	),
    dropMetaData = cms.untracked.string('ALL'),
    splitLevel = cms.untracked.int32(99),
        SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )

    )


process.out = cms.OutputModule(
    'PoolOutputModule',
    fileName       = cms.untracked.string('fanculo.edm.root'),
    outputCommands = cms.untracked.vstring('drop *','keep *_HbbAnalyzerNew_*_*')
    )

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
process.load("CommonTools.ParticleFlow.PF2PAT_EventContent_cff")
process.out.outputCommands.extend(process.PATEventContent.outputCommands)

# HB + HE noise filtering
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')

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

process.patJets.addTagInfos  = True
#process.selectedPatJets.addTagInfos  = True
###process.selectedPatJetsAK5PF.addTagInfos  = True
#process.selectedPatJetsAK7Calo.addTagInfos  = True
#process.selectedPatJetsAK7PF.addTagInfos  = True

process.load('CommonTools.ParticleFlow.PF2PAT_cff')
from PhysicsTools.PatAlgos.tools.pfTools import *
usePF2PAT(process,runPF2PAT=True,jetAlgo='AK5',runOnMC=isMC,postfix='')

# Get a list of good primary vertices, in 42x, these are DAF vertices
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )

process.pfPileUp.Enable = True
process.pfIsolatedElectrons.combinedIsolationCut = 0.35
process.pfIsolatedMuons.combinedIsolationCut     = 0.35
process.pfPileUp.checkClosestZVertex = cms.bool(False)
process.pfPileUp.Vertices = cms.InputTag('goodOfflinePrimaryVertices')
process.pfJets.doAreaFastjet = True
process.pfJets.doRhoFastjet = False


# Compute the mean pt per unit area (rho) from the
# PFchs inputs
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJets = kt4PFJets.clone(
    rParam = cms.double(0.6),
    src = cms.InputTag('pfNoElectron'),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True)
    )
process.patJetCorrFactors.rho = cms.InputTag("kt6PFJets", "rho")


process.load('RecoJets.JetProducers.ak5PFJets_cfi')
process.ak5PFJets.jetPtMin      = 1.0
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.src           = 'pfNoElectron'


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


if isMC == False :
        inputJetCorrLabel = ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']
else :
        inputJetCorrLabel = ['L1FastJet','L2Relative', 'L3Absolute']

addJetCollection(process, 
        cms.InputTag('caTopCaloJets'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
        "CA","Calo",
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
#        jetCorrLabel=('KT6Calo', inputJetCorrLabel),
        jetCorrLabel=('AK7Calo', inputJetCorrLabel),
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        doJetID=False
                 )

addJetCollection(process, 
        cms.InputTag('caTopPFJets'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
        'CA',"PF",
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
#        jetCorrLabel=('KT6PF',  inputJetCorrLabel),
        jetCorrLabel=('AK7PF',  inputJetCorrLabel),
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        doJetID = False
                 )

addJetCollection(process,cms.InputTag("CAsubJetsProducer"),"subCA","PF",
                 doJTA=True,doBTagging=True,jetCorrLabel=('KT4PF',  inputJetCorrLabel),
                 doType1MET=False,doL1Cleaning = False,doL1Counters=False,
                 doJetID = False)

addJetCollection(process,cms.InputTag("ak7CaloJets"),"AK7","Calo",
                 doJTA=True,doBTagging=True,jetCorrLabel=('AK7Calo',  inputJetCorrLabel),
                 doType1MET=False,doL1Cleaning = False,doL1Counters=False,
                 doJetID = False)

if isMC == False :
        addJetCollection(process,cms.InputTag("ak5PFJets"),"AK5","PF",
                 doJTA=True,doBTagging=True,jetCorrLabel=('AK5PFchs',  inputJetCorrLabel),
                 doType1MET=False,doL1Cleaning = False,doL1Counters=False,
                 doJetID = False)
else :
        addJetCollection(process,cms.InputTag("ak5PFJets"),"AK5","PF",
                 doJTA=True,doBTagging=True,jetCorrLabel=('AK5PFchs',  inputJetCorrLabel),
                 doType1MET=False, genJetCollection=cms.InputTag("ak5GenJets"),doL1Cleaning = False,doL1Counters=False,
                 doJetID = False)
	
addJetCollection(process,cms.InputTag("ak7PFJets"),"AK7","PF",
                 doJTA=True,doBTagging=True,jetCorrLabel=('AK7PF',  inputJetCorrLabel),
                 doType1MET=False,doL1Cleaning = False,doL1Counters=False,
                 doJetID = False)


switchJetCollection(process, 
        cms.InputTag('ak5CaloJets'),     # Jet collection; must be already in the event when patLayer0 sequence is executed
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
        jetCorrLabel=('AK5Calo', inputJetCorrLabel), 
        doType1MET=False,
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
	jetcoll.embedGenJetMatch = cms.bool(isMC)
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


#BTAGGING SF FROM DATABASE
#How to from :https://twiki.cern.ch/twiki/bin/view/CMS/BtagPerformanceDBV2#Usage_Examples
process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")
process.load ("RecoBTag.PerformanceDB.BTagPerformanceDB1107")


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
#process.selectedPatElectrons.cut = cms.string('')

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



# prune gen particles
## AR: is it used???
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")


from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector 
process.goodPatJetsAK5PF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                           filterParams = pfJetIDSelector.clone(), src = cms.InputTag("selectedPatJetsAK5PF") )


process.HbbAnalyzerNew = cms.EDProducer("HbbAnalyzerNew",
    runOnMC = cms.bool(isMC),
    hltResultsTag = cms.InputTag("TriggerResults::HLT"),
    electronTag = cms.InputTag("selectedElectronsMatched"),
#    electronTag = cms.InputTag("selectedPatElectrons"),
    tauTag = cms.InputTag("patTaus"),
#   muonTag = cms.InputTag("selectedPatMuons"),
    muonTag = cms.InputTag("selectedMuonsMatched"),
    jetTag = cms.InputTag("selectedPatJetsCAPF"),
    subjetTag = cms.InputTag("selectedPatJetssubCAPF"),
    simplejet1Tag = cms.InputTag("selectedPatJets"),
    simplejet2Tag = cms.InputTag("selectedPatJetsAK5PF"),
    simplejet3Tag = cms.InputTag("selectedPatJetsAK7Calo"),
    simplejet4Tag = cms.InputTag("selectedPatJetsAK7PF"),
    photonTag = cms.InputTag("selectedPatPhotons"),
    metTag = cms.InputTag("patMETs"),
    dimuTag = cms.InputTag("dimuons"),
    dielecTag = cms.InputTag("dielectrons"),

    verbose = cms.untracked.bool(False)
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

import VHbbAnalysis.HbbAnalyzer.simpleCutBasedElectronIDSpring10_cfi as vbtfid

process.eidVBTFRel95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95relIso' )
process.eidVBTFRel85 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '85relIso' )
process.eidVBTFRel80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80relIso' )
process.eidVBTFRel70 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '70relIso' )
process.eidVBTFCom95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95cIso'   )
process.eidVBTFCom85 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '85cIso'   )
process.eidVBTFCom80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80cIso'   )
process.eidVBTFCom70 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '70cIso'   )

process.eidSequence = cms.Sequence(
    process.eidVBTFRel95 +
    process.eidVBTFRel85 +
    process.eidVBTFRel80 +
    process.eidVBTFRel70 +
    process.eidVBTFCom95 +
    process.eidVBTFCom85 +
    process.eidVBTFCom80 +
    process.eidVBTFCom70
)

# Electron Selection
process.patElectrons.electronIDSources = cms.PSet(
    eidVBTFRel95 = cms.InputTag("eidVBTFRel95"),
    eidVBTFRel85 = cms.InputTag("eidVBTFRel85"),
    eidVBTFRel80 = cms.InputTag("eidVBTFRel80"),
    eidVBTFRel70 = cms.InputTag("eidVBTFRel70"),
    eidVBTFCom95 = cms.InputTag("eidVBTFCom95"),
    eidVBTFCom85 = cms.InputTag("eidVBTFCom85"),
    eidVBTFCom80 = cms.InputTag("eidVBTFCom80"),
    eidVBTFCom70 = cms.InputTag("eidVBTFCom70")
)
###################### trigger matching if isMC=false
switchOnTrigger(process, outputModule="")

defaultTriggerMatch = cms.EDProducer(
  "PATTriggerMatcherDRDPtLessByR"                                               # match by DeltaR only, best match by DeltaR
, src     = cms.InputTag( "selectedPatMuons" )
, matched = cms.InputTag( "patTrigger" )                                        # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
, matchedCuts = cms.string( 'path( "*Mu*" )' )
, maxDPtRel = cms.double( 0.5 )
, maxDeltaR = cms.double( 0.3 )
, resolveAmbiguities    = cms.bool( True )                                      # only one match per trigger object
, resolveByMatchQuality = cms.bool( True )                                      # take best match found per reco object: by DeltaR here (s. above)
)

process.selectedMuonsTriggerMatch = defaultTriggerMatch.clone(
        src         = cms.InputTag( "selectedPatMuons" )
        , matchedCuts = cms.string('path("*Mu*",0,0)')
#	        , matchedCuts = cms.string('path("HLT_IsoMu17_v*")|| path("HLT_Mu15_v*")')
)
# trigger object embedders for the same collections
process.selectedMuonsMatched = cms.EDProducer( "PATTriggerMatchMuonEmbedder",
        src     = cms.InputTag(  "selectedPatMuons" ),
        matches = cms.VInputTag( cms.InputTag('selectedMuonsTriggerMatch') )
    )

process.selectedElectronsTriggerMatch = defaultTriggerMatch.clone(
        src         = cms.InputTag( "selectedPatElectrons" )
        , matchedCuts = cms.string('path("*Ele*",0,0)')
#	        , matchedCuts = cms.string('path("HLT_IsoMu17_v*")|| path("HLT_Mu15_v*")')
)
# trigger object embedders for the same collections
process.selectedElectronsMatched = cms.EDProducer( "PATTriggerMatchElectronEmbedder",
        src     = cms.InputTag(  "selectedPatElectrons" ),
        matches = cms.VInputTag( cms.InputTag('selectedElectronsTriggerMatch') )
    )


process.leptonTrigMatch = cms.Sequence (
	process.selectedElectronsTriggerMatch
	*process.selectedElectronsMatched
	*process.selectedMuonsTriggerMatch
	*process.selectedMuonsMatched
	)


process.hbbBestCSVPt20Candidates = cms.EDFilter("HbbCandidateFinder",
                                       VHbbEventLabel = cms.InputTag(""),
                                       verbose = cms.bool(False) ,
                                       useHighestPtHiggs = cms.bool(False),
                                       jetPtThreshold = cms.double(20.),
                                       actAsAFilter = cms.bool(False)
                                      )

process.hbbHighestPtHiggsPt30Candidates = cms.EDFilter("HbbCandidateFinder",
                                       VHbbEventLabel = cms.InputTag(""),
                                       verbose = cms.bool(False) ,
                                       useHighestPtHiggs = cms.bool(True),
                                       jetPtThreshold = cms.double(30.),
                                       actAsAFilter = cms.bool(False)
                                      )

process.hbbCandidates = cms.EDFilter("HbbCandidateFinder",
				       VHbbEventLabel = cms.InputTag(""),
				       verbose = cms.bool(False) ,
				       jetPtThreshold = cms.double(30.),
				       useHighestPtHiggs=cms.bool(False),
              			       actAsAFilter = cms.bool(False)
				      )


process.patMETsHT = cms.EDProducer("MHTProducer",
  JetCollection = cms.InputTag("patJets"),
  MinJetPt      = cms.double(30),
  MaxJetEta     = cms.double(5)
)


process.dump = cms.EDAnalyzer("EventContentAnalyzer")

# drop the meta data for dropped data
#process.out.dropMetaData = cms.string("DROPPED")

#process.out.fileName = '/tigress-hsm/dlopes/PatEDM.root'


if isMC == False :
        process.p = cms.Path(
                    process.goodOfflinePrimaryVertices*
                     process.PF2PAT*
                     process.ak5CaloJets*
                     process.ak7CaloJets*
                     process.kt6PFJets*
                     process.ak5PFJets*
                     process.ak7PFJets*
                     process.caTopCaloJets*
                     process.caTopPFJets*
                     process.CAsubJetsProducer*
                     process.eidSequence*
                     process.patDefaultSequence*
#                     process.patPF2PATSequence* # added with usePF2PAT
		     process.patMETsHT*
                     process.dimuons*
                     process.dielectrons*
                     process.leptonTrigMatch*
                     process.HbbAnalyzerNew
#process.hbbCandidates*process.hbbHighestPtHiggsPt30Candidates*process.hbbBestCSVPt20Candidates
                     )
else :
        process.p = cms.Path(
                    process.goodOfflinePrimaryVertices*
                     process.genParticlesForJets*
                     process.ak5GenJets*
                     process.PF2PAT*
                     process.ak5CaloJets*
                     process.ak7CaloJets*
                     process.kt6PFJets*
                     process.ak5PFJets*
                     process.ak7PFJets*
                     process.caTopCaloJets*
                     process.caTopPFJets*
                     process.CAsubJetsProducer*
                     process.eidSequence*
                     process.patDefaultSequence*
#                     process.patPF2PATSequence* # added with usePF2PAT
		     process.patMETsHT*
#		     process.dump*
                     process.dimuons*
                     process.dielectrons*
                     process.leptonTrigMatch*
                     process.HbbAnalyzerNew
#process.hbbCandidates*process.hbbHighestPtHiggsPt30Candidates*process.hbbBestCSVPt20Candidates
                     )

process.hbhepath = cms.Path(process.HBHENoiseFilter)



#process.candidates = cms.Path(process.hbbCandidates*process.hbbHighestPtHiggsPt30Candidates*process.hbbBestCSVPt20Candidates)


#process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )

process.options = cms.untracked.PSet(
	wantSummary = cms.untracked.bool(True),
	Rethrow = cms.untracked.vstring('ProductNotFound')
	)
process.e = cms.EndPath(process.out1)


process.schedule = cms.Schedule(process.p, process.hbhepath, process.e)

#
