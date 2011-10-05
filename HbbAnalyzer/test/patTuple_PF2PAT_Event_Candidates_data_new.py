
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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

# source
process.source = cms.Source("PoolSource",
			    fileNames=cms.untracked.vstring(
#"rfio:/castor/cern.ch/user/d/degrutto/test/DYJetsToLL_PtZ-100.root",

#'file:/gpfs/gpfsddn/cms/user/arizzi/Hbb/submit/CMSSW_4_2_8_patch1/src/VHbbAnalysis/VHbbDataFormats/bin/submissions/testbortigno/24233412-65AD-E011-B930-E0CB4E553667.root'

#	'root://cmsdcache7.pi.infn.it:7070//store/mc/Summer11/ZH_ZToLL_HToBB_M-115_7TeV-powheg-herwigpp/AODSIM/PU_S4_START42_V11-v1/0000/02E676EE-BDAD-E011-B9ED-E0CB4EA0A929.root',
	'root://cmsdcache7.pi.infn.it:7070//store/mc/Summer11/ZH_ZToLL_HToBB_M-115_7TeV-powheg-herwigpp/AODSIM/PU_S4_START42_V11-v1/0000/32CECED6-BFAD-E011-B08D-00261834B5A4.root',
#	'root://cmsdcache7.pi.infn.it:7070//store/mc/Summer11/ZH_ZToLL_HToBB_M-115_7TeV-powheg-herwigpp/AODSIM/PU_S4_START42_V11-v1/0000/3EE916B3-C4AD-E011-9159-90E6BA0D09B0.root',
#	'root://cmsdcache7.pi.infn.it:7070//store/mc/Summer11/ZH_ZToLL_HToBB_M-115_7TeV-powheg-herwigpp/AODSIM/PU_S4_START42_V11-v1/0000/42622800-C1AD-E011-AD65-485B39800BF2.root',
#		"/store/mc/CMSSW_4_2_3/RelValProdTTbar/GEN-SIM-RECO/MC_42_V12_JobRobot-v1/0000/B89A0B07-818C-E011-953E-0030487CD7E0.root"
#	"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/mc/Summer11/BdToMuMu_2MuPtFilter_7TeV-pythia6-evtgen//GEN-SIM-RECO//PU_S4_START42_V11-v1///0000//92DE42C9-CD8C-E011-A421-001F296B758E.root"
	),
)
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if isMC == False :
	process.GlobalTag.globaltag = cms.string('GR_R_42_V19::All')
else :
	process.GlobalTag.globaltag = cms.string('START42_V13::All')

process.out1 = cms.OutputModule(
    'PoolOutputModule',
    fileName       = cms.untracked.string('PAT.edm.root'),
    outputCommands = cms.untracked.vstring(
	'drop *',
	'keep *_HbbAnalyzerNew_*_*',
	'keep VHbbCandidates_*_*_*',
#	'keep PileupSummaryInfos_*_*_*',
	'keep edmTriggerResults_*_*_*',
#        'keep *_hltTriggerSummaryAOD_*_*',
#        'keep *_selectedVertices_*_*',
#        'keep *_hltTriggerSummaryAOD_*_*',
#        'keep *_TriggerResults_*_*',
	'keep *_bcandidates_*_*',
	'keep *_bhadrons_*_*',
       	"keep *_HLTDiCentralJet20MET80_*_*",
       	"keep *_HLTDiCentralJet20MET100HBHENoiseFiltered_*_*",
	"keep *_HLTPFMHT150_*_*",
	"keep *_HLTQuadJet40_*_*",
	"keep *_HLTDoubleMu7_*_*",

	),
    dropMetaData = cms.untracked.string('ALL'),
    splitLevel = cms.untracked.int32(99),
        SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )

    )


process.out = cms.OutputModule(
    'PoolOutputModule',
    fileName       = cms.untracked.string('fake.edm.root'),
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
process.pfNoMuon.enable = False
process.pfNoElectron.enable = False

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

### to compute rho to subtract to lepton isolation cone, so with particles up to eta 2.5  
process.kt6PFJets25 = process.kt4PFJets.clone(
	src = 'pfNoElectron',
	rParam = 0.6,
	doRhoFastjet = True,
	Ghost_EtaMax = 2.5,
	Rho_EtaMax = 2.5
	)



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
process.selectedPatJetsAK7PF.cut = cms.string('pt > 15. & abs(eta) < 5.0')
process.selectedPatJetsCACalo.cut = cms.string('pt > 0. & abs(eta) < 5.0')
process.selectedPatJetsCAPF.cut = cms.string('pt > 0. & abs(eta) < 5.0')
##process.selectedLayer1Muons.cut = cms.string('pt > 20. & abs(eta) < 2.5 & muonID("TMLastStationLoose")')
##process.selectedLayer1Electrons.cut = cms.string('pt > 20. & abs(eta) < 2.5 & electronID("eidLoose")')
process.selectedPatMuons.cut = cms.string('isGlobalMuon || (isTrackerMuon && muonID("TrackerMuonArbitrated"))')
#process.selectedPatElectrons.cut = cms.string('')

process.selectedPatMuonsWithIso = process.selectedPatMuons.clone(
    src = 'selectedPatMuons',
     cut =  (
        "pt > 20 && isGlobalMuon && isTrackerMuon && globalTrack().normalizedChi2 < 10 &&" +
            "innerTrack().hitPattern().numberOfValidTrackerHits > 10 && "                      +
            "innerTrack().hitPattern().numberOfValidPixelHits > 0 && "                         +
            "globalTrack().hitPattern().numberOfValidMuonHits > 0 && "                         +
            "dB < 0.2 && "                                                                     +
            " chargedHadronIso + neutralHadronIso + photonIso  < 0.15 * pt && "                                               +
            "numberOfMatches > 1 && abs(eta) < 2.4"
        )
)

process.selectedPatElectronsWithIso = process.selectedPatElectrons.clone(
    src = 'selectedPatElectrons',
    cut = (
   "pt > 15.0 && abs(eta) < 2.5 &&"
   "(isEE || isEB) && !isEBEEGap &&"
   "ecalDrivenSeed==1 && (abs(superCluster.eta)<2.5)"
   " && !(1.4442<abs(superCluster.eta)<1.566)"
   " && (ecalEnergy*sin(superClusterPosition.theta)>20.0)"
   " && (gsfTrack.trackerExpectedHitsInner.numberOfHits == 0)"
    " && ((chargedHadronIso + neutralHadronIso + photonIso < 0.15 * pt)) "
   " && ((isEB"
   " && (sigmaIetaIeta<0.01)"
   " && ( -0.8<deltaPhiSuperClusterTrackAtVtx<0.8 )"
   " && ( -0.007<deltaEtaSuperClusterTrackAtVtx<0.007 )"
   ")"
   " || (isEE"
   " && (sigmaIetaIeta<0.03)"
   " && ( -0.7<deltaPhiSuperClusterTrackAtVtx<0.7 )"
   " && ( -0.01<deltaEtaSuperClusterTrackAtVtx<0.01 )"
   "))"
   ))


process.cleanPatJetsAK5PF = cms.EDProducer("PATJetCleaner",
					   src = cms.InputTag("selectedPatJetsAK5PF"),
					   # preselection = cms.string('pt > 20.0  && ( (neutralEmEnergy/energy < 0.99) &&  (neutralHadronEnergy/energy < 0.99) && numberOfDaughters>1) '),
					   preselection = cms.string('abs(eta)< 5.0 && pt > 15.0'),
					   
					   
					   checkOverlaps = cms.PSet(
	muons = cms.PSet(
	src       = cms.InputTag("selectedPatMuonsWithIso"),
	algorithm = cms.string("byDeltaR"),
	preselection        = cms.string(""),
	deltaR              = cms.double(0.5),
	checkRecoComponents = cms.bool(False),
	pairCut             = cms.string(""),
	requireNoOverlaps   = cms.bool(True),
	),
	electrons = cms.PSet(
	src       = cms.InputTag("selectedPatElectronsWithIso"),
	algorithm = cms.string("byDeltaR"),
	preselection        = cms.string(""),
	deltaR              = cms.double(0.5),
	checkRecoComponents = cms.bool(False),
	pairCut             = cms.string(""),
	requireNoOverlaps   = cms.bool(True),
	),
	),
					   finalCut = cms.string('')
					   
					   )
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

### rerun isolation for muons and electrons
process.load("VHbbAnalysis.HbbAnalyzer.pfIsoDeltaBeta04_cff")

process.selectedPatElectrons.cut = (
    "(ecalDrivenSeed==1) &&"                                       +
    "pt > 15.0 && abs(eta) < 2.5 &&"                               +
    "(isEE || isEB) && !isEBEEGap"                              
)

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
    simplejet2Tag = cms.InputTag("cleanPatJetsAK5PF"),
    simplejet4Tag = cms.InputTag("selectedPatJetsAK7Calo"),
    simplejet3Tag = cms.InputTag("selectedPatJetsAK7PF"),
    photonTag = cms.InputTag("selectedPatPhotons"),
    metTag = cms.InputTag("met"), #this input tag is used to fill calo MET 

    verbose = cms.untracked.bool(False)
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
  JetCollection = cms.InputTag("patJetsAK5PF"),
  MinJetPt      = cms.double(30),
  MaxJetEta     = cms.double(5)
)

process.pfMETNoPU = process.pfMET.clone()
process.pfMETNoPU.src=cms.InputTag("pfNoPileUp")

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

# drop the meta data for dropped data
#process.out.dropMetaData = cms.string("DROPPED")
process.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexing_cff')
process.inclusiveMergedVertices = process.vertexMerger.clone()
process.inclusiveMergedVertices.secondaryVertices = cms.InputTag("inclusiveVertices")
process.inclusiveMergedVertices.maxFraction = 0.2
process.inclusiveMergedVertices.minSignificance = cms.double(10.)

process.load("RecoBTag/SecondaryVertex/bVertexFilter_cfi")
process.selectedVertices = process.bVertexFilter.clone()
process.selectedVertices.secondaryVertices = cms.InputTag("inclusiveMergedVertices")
process.selectedVertices.minVertices = 0
process.selectedVertices.vertexFilter.multiplicityMin = 3

process.inclusiveVertexFinder.clusterScale = 1.
process.inclusiveVertexFinder.clusterMinAngleCosine = 0.5


#process.out.fileName = '/tigress-hsm/dlopes/PatEDM.root'
process.bcandidates = cms.EDProducer('BCandidateProducer',
                                     src = cms.InputTag('selectedVertices','',''),
                                     primaryVertices =
                                     cms.InputTag('offlinePrimaryVerticesWithBS','',''),
                                     minDRUnique = cms.untracked.double(0.4),
                                     minvecSumIMifsmallDRUnique = cms.untracked.double(5.5),
                                     minCosPAtomerge = cms.untracked.double(0.99),
                                     maxPtreltomerge = cms.untracked.double(7777.0)
                                     )


#
# ntuplizer di souvik
#
#

process.HLTDiCentralJet20MET80 = cms.EDProducer("HLTInfoDumperGeneral",
  HLTPath = cms.untracked.string('HLT_DiCentralJet20_MET80_v'),
  TriggerResults = cms.untracked.InputTag('TriggerResults', '', 'HLT'),
  TriggerEvent = cms.untracked.InputTag('hltTriggerSummaryAOD', '', 'HLT'),
  FilterNames = cms.untracked.VInputTag(
    cms.InputTag('hlt2CenJet20CentralRegional', '', 'HLT'),
    cms.InputTag('hltMET80', '', 'HLT')
  )
)

process.HLTDiCentralJet20MET100HBHENoiseFiltered = cms.EDProducer("HLTInfoDumperGeneral",
  HLTPath = cms.untracked.string('HLT_DiCentralJet20_MET100_HBHENoiseFiltered_v'),
  TriggerResults = cms.untracked.InputTag('TriggerResults', '', 'HLT'),
  TriggerEvent = cms.untracked.InputTag('hltTriggerSummaryAOD', '', 'HLT'),
  FilterNames = cms.untracked.VInputTag(
    cms.InputTag('hlt2CenJet20CentralRegional', '', 'HLT'),
    cms.InputTag('hltMET100', '', 'HLT')
  )
)

process.HLTPFMHT150 = cms.EDProducer("HLTInfoDumperGeneral",
  HLTPath = cms.untracked.string('HLT_PFMHT150_v'),
  TriggerResults = cms.untracked.InputTag('TriggerResults', '', 'HLT'),
  TriggerEvent = cms.untracked.InputTag('hltTriggerSummaryAOD', '', 'HLT'),
  FilterNames = cms.untracked.VInputTag(
    cms.InputTag('hltMET80', '', 'HLT'),
    cms.InputTag('hltPFMHT150Filter', '', 'HLT')
  )
)

process.HLTQuadJet40 = cms.EDProducer("HLTInfoDumperGeneral",
  HLTPath = cms.untracked.string('HLT_QuadJet40_v'),
  TriggerResults = cms.untracked.InputTag('TriggerResults', '', 'HLT'),
  TriggerEvent = cms.untracked.InputTag('hltTriggerSummaryAOD', '', 'HLT'),
  FilterNames = cms.untracked.VInputTag(
    cms.InputTag('hltQuadJet40Central', '', 'HLT')
  )
)

process.HLTDoubleMu7 = cms.EDProducer("HLTInfoDumperGeneral",
  HLTPath = cms.untracked.string('HLT_DoubleMu7_v5'),
  TriggerResults = cms.untracked.InputTag('TriggerResults', '', 'HLT'),
  TriggerEvent = cms.untracked.InputTag('hltTriggerSummaryAOD', '', 'HLT'),
  FilterNames = cms.untracked.VInputTag(
    cms.InputTag('hltL1DoubleMuon3L1Filtered0', '', 'HLT'),
    cms.InputTag('hltDiMuon3L2PreFiltered0', '', 'HLT'),
    cms.InputTag('hltDiMuonL3PreFiltered7', '', 'HLT')
  )
)




#simBHadron producer                                                                                                                                                                 
process.bhadrons = cms.EDProducer('MCBHadronProducer',
                                 quarkId = cms.uint32(5)
                                 )

### Paths ###

process.nTuplizePath=cms.Path(process.HLTDiCentralJet20MET80 +
                              process.HLTDiCentralJet20MET100HBHENoiseFiltered +
                              process.HLTPFMHT150 +
                              process.HLTQuadJet40 +
                              process.HLTDoubleMu7)



if isMC == False :
        process.p = cms.Path(
                    process.goodOfflinePrimaryVertices*
                     process.PF2PAT*
		     process.pfCandsForIsolationSequence *
		     process.muonPFIsolationSequence *
		     process.electronPFIsolationSequence *
                     process.ak5CaloJets*
                     process.ak7CaloJets*
                     process.kt6PFJets*
		     process.kt6PFJets25* 
                     process.ak5PFJets*
                     process.ak7PFJets*
                     process.caTopCaloJets*
                     process.caTopPFJets*
                     process.CAsubJetsProducer*
                     process.eidSequence*
                     process.patDefaultSequence*
		    		    #
		    # add no iso stuff
		    #
		    process.selectedPatMuonsWithIso*
		    process.selectedPatElectronsWithIso*
		    process.cleanPatJetsAK5PF*

#                     process.patPF2PATSequence* # added with usePF2PAT
		     process.patMETsHT*
		     process.pfMETNoPU*
                     process.leptonTrigMatch*
                     process.inclusiveVertexing*
                     process.inclusiveMergedVertices*process.selectedVertices*
                     process.bcandidates*
                     process.HbbAnalyzerNew

#process.hbbCandidates*process.hbbHighestPtHiggsPt30Candidates*process.hbbBestCSVPt20Candidates
                     )
else :
        process.p = cms.Path(
                    process.goodOfflinePrimaryVertices*
                     process.genParticlesForJets*
                     process.ak5GenJets*
                     process.PF2PAT*
		    process.pfCandsForIsolationSequence *
		     process.muonPFIsolationSequence *
		     process.electronPFIsolationSequence *
                     process.ak5CaloJets*
                     process.ak7CaloJets*
                     process.kt6PFJets*
	      	     process.kt6PFJets25* 
                     process.ak5PFJets*
                     process.ak7PFJets*
                     process.caTopCaloJets*
                     process.caTopPFJets*
                     process.CAsubJetsProducer*
                     process.eidSequence*
                     process.patDefaultSequence*
		    #
		    # add no iso stuff
		    #
		    process.selectedPatMuonsWithIso*
		    process.selectedPatElectronsWithIso*
		    process.cleanPatJetsAK5PF*
#                     process.patPF2PATSequence* # added with usePF2PAT
		     process.patMETsHT*
		     process.pfMETNoPU*
#		     process.dump*
                     process.leptonTrigMatch*
                     process.inclusiveVertexing*
                     process.inclusiveMergedVertices*process.selectedVertices*
                     process.bcandidates*
		     process.bhadrons*
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


process.schedule = cms.Schedule(process.p, process.hbhepath, process.nTuplizePath ,process.e)

#
