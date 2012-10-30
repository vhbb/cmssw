#import configurations
import FWCore.ParameterSet.Config as cms

import os 

isMC = True

# define the process
process = cms.Process("VH")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")


## import skeleton process
#from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


# source
process.source = cms.Source("PoolSource",
			    fileNames=cms.untracked.vstring(
"file:/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/mc/Fall11/ZH_ZToLL_HToBB_M-125_7TeV-powheg-herwigpp/AODSIM/PU_S6_START42_V14B-v1/0000/4AC41B53-B60C-E111-B38B-20CF3027A577.root",
"file:/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/mc/Fall11/ZH_ZToLL_HToBB_M-125_7TeV-powheg-herwigpp/AODSIM/PU_S6_START42_V14B-v1/0000/4E4F81B7-C10C-E111-984A-485B39800BA0.root",
"file:/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/mc/Fall11/ZH_ZToLL_HToBB_M-125_7TeV-powheg-herwigpp/AODSIM/PU_S6_START42_V14B-v1/0000/7EA0B77B-FF0C-E111-A498-20CF3027A5E9.root",
"file:/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/mc/Fall11/ZH_ZToLL_HToBB_M-125_7TeV-powheg-herwigpp/AODSIM/PU_S6_START42_V14B-v1/0000/84F0DC03-B10C-E111-8757-E0CB4E553644.root",
"file:/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/mc/Fall11/ZH_ZToLL_HToBB_M-125_7TeV-powheg-herwigpp/AODSIM/PU_S6_START42_V14B-v1/0000/9C1AF259-B30C-E111-9449-90E6BA0D09D7.root",
"file:/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/mc/Fall11/ZH_ZToLL_HToBB_M-125_7TeV-powheg-herwigpp/AODSIM/PU_S6_START42_V14B-v1/0000/A83884EC-B90C-E111-BA3D-E0CB4E29C503.root",
"file:/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/mc/Fall11/ZH_ZToLL_HToBB_M-125_7TeV-powheg-herwigpp/AODSIM/PU_S6_START42_V14B-v1/0000/A8D1BB82-AD0C-E111-8F15-20CF3019DEF7.root",
"file:/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/mc/Fall11/ZH_ZToLL_HToBB_M-125_7TeV-powheg-herwigpp/AODSIM/PU_S6_START42_V14B-v1/0000/B88CF49F-AA0C-E111-9C8B-E0CB4E19F96D.root",
"file:/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/mc/Fall11/ZH_ZToLL_HToBB_M-125_7TeV-powheg-herwigpp/AODSIM/PU_S6_START42_V14B-v1/0000/E091B3C4-B20C-E111-BDAC-00261834B53C.root",
"file:/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/mc/Fall11/ZH_ZToLL_HToBB_M-125_7TeV-powheg-herwigpp/AODSIM/PU_S6_START42_V14B-v1/0000/F063894C-AF0C-E111-989D-485B39800B62.root",
"file:/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/mc/Fall11/ZH_ZToLL_HToBB_M-125_7TeV-powheg-herwigpp/AODSIM/PU_S6_START42_V14B-v1/0000/FA629A0A-B90C-E111-8325-90E6BA442F2C.root"

),
)
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if isMC == False :
	process.GlobalTag.globaltag = cms.string('GR_R_42_V23::All')
else :
	process.GlobalTag.globaltag = cms.string('START42_V17::All')


process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printSavedList = cms.EDAnalyzer("ParticleListDrawer",
					     src = cms.InputTag("savedGenParticles"),
					     maxEventsToPrint = cms.untracked.int32(-1)
					 )
process.printList = cms.EDAnalyzer("ParticleListDrawer",
					     src = cms.InputTag("genParticles"),
					     maxEventsToPrint = cms.untracked.int32(-1)
					 )

process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
					     src = cms.InputTag("genParticles"),
					     printP4 = cms.untracked.bool(False),
					     printPtEtaPhi = cms.untracked.bool(False),
					     printVertex = cms.untracked.bool(True),
					     printStatus = cms.untracked.bool(False),
					     printIndex = cms.untracked.bool(False),
					     status = cms.untracked.vint32(1, 2, 3)
					 )

process.out1 = cms.OutputModule(
    'PoolOutputModule',
    fileName       = cms.untracked.string('PAT.edm.root'),
    outputCommands = cms.untracked.vstring(
	'drop *',
        'keep *_savedGenParticles_*_*',
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
	'keep *_pfMEtSignCovMatrix_*_*',
       	"keep *_HLTDiCentralJet20MET80_*_*",
       	"keep *_HLTDiCentralJet20MET100HBHENoiseFiltered_*_*",
	"keep *_HLTPFMHT150_*_*",
	"keep *_HLTQuadJet40_*_*",
	"keep *_HLTDoubleMu7_*_*",
	"keep *_EcalDeadCellEventFilter_*_*",
        "keep *_patType1CorrectedPFMet*_*_*",
	"keep *_patType1p2CorrectedPFMet*_*_*",
	"keep LHEEventProduct_source_*_LHE"
	
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

### BEGIN of filters required by JETMET for  MET analysis 
# HB + HE noise filtering
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
# ecal filter
process.load('JetMETAnalysis.ecalDeadCellTools.EcalDeadCellEventFilter_cfi') ## still not in the release
# CSC beam halo.... 
process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')
# HCAL laser filter
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
# tracking failurefilter
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.trackingFailureFilter.JetSource = cms.InputTag('ak5PFJetsL2L3Residual')
process.trackingFailureFilter.VertexSource = cms.InputTag('goodOfflinePrimaryVertices')
### END of filters required by JETMET for  MET analysis 



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

#### Taus ####
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

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
### already included in JetMETCorrections/Type1MET/python/pfMETCorrections_cff.py, just change the inputs
process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJets.src = cms.InputTag('pfNoElectron')

#process.kt6PFJets = process.kt4PFJets.clone(
#    rParam = cms.double(0.6),
#    src = cms.InputTag('pfNoElectron'),
#    doAreaFastjet = cms.bool(True),
#    doRhoFastjet = cms.bool(True)
#    )
#process.patJetCorrFactors.rho = cms.InputTag("kt6PFJets", "rho")

### to compute rho to subtract to lepton isolation cone, so with particles up to eta 2.5  

process.kt6PFJets25 = process.kt4PFJets.clone(
	src = 'pfNoElectron',
	rParam = 0.6,
	doRhoFastjet = True,
	Ghost_EtaMax = 2.5,
	Rho_EtaMax = 2.5
	)




from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
process.ak5PFJets = ak5PFJets
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
###process.load("VHbbAnalysis.HbbAnalyzer.HbbFilterJets_cff")

process.ca8GenJets = ca4GenJets.clone( rParam = cms.double(0.8) )


process.load('RecoJets.JetProducers.caSubjetFilterPFJets_cfi')
from RecoJets.JetProducers.caSubjetFilterPFJets_cfi import caSubjetFilterPFJets
process.caVHPFJets = caSubjetFilterPFJets.clone()

process.load('RecoJets.JetProducers.caSubjetFilterCaloJets_cfi')
from RecoJets.JetProducers.caSubjetFilterCaloJets_cfi import caSubjetFilterCaloJets
process.caVHCaloJets = caSubjetFilterCaloJets.clone()

process.load('RecoJets.JetProducers.caSubjetFilterGenJets_cfi')
from RecoJets.JetProducers.caSubjetFilterGenJets_cfi import caSubjetFilterGenJets
process.caVHGenJets = caSubjetFilterGenJets.clone()


## tune DeltaR for C/A

##process.caTopCaloJets.rBins = cms.vdouble(1.2,1.2,1.2)
##process.caTopPFJets.rBins = cms.vdouble(1.2,1.2,1.2)

##process.caTopCaloJets.rParam = cms.double(1.2)
##process.caTopPFJets.rParam = cms.double(1.2)

#process.caTopCaloJets.rBins = cms.vdouble(1.4,1.4,1.4)
#process.caTopPFJets.rBins = cms.vdouble(1.4,1.4,1.4)
#process.caTopCaloJets.rParam = cms.double(1.4)
#process.caTopPFJets.rParam = cms.double(1.4)

##test
process.caVHCaloJets.useAdjacency = cms.int32(0)
process.caVHPFJets.useAdjacency = cms.int32(0)
######

process.ak5PFJets.src = cms.InputTag('pfNoElectron')
process.ak7PFJets.src = cms.InputTag('pfNoElectron')
process.caVHPFJets.src = cms.InputTag('pfNoElectron')


# subJetProducer
##process.CAsubJetsProducer = cms.EDProducer('SubProducer',
##jetTag = cms.untracked.InputTag("caTopPFJets"))


# switch jet collection to our juets
from PhysicsTools.PatAlgos.tools.jetTools import *


if isMC == False :
        inputJetCorrLabel = ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']
else :
        inputJetCorrLabel = ['L1FastJet','L2Relative', 'L3Absolute']

addJetCollection(process, 
        cms.InputTag('caVHCaloJets:fat'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
        "CAVHFat","Calo",
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
        cms.InputTag('caVHCaloJets:sub'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
        "CAVHSub","Calo",
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
        cms.InputTag('caVHCaloJets:filter'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
        "CAVHFilter","Calo",
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
#        jetCorrLabel=('KT6Calo', inputJetCorrLabel),
        jetCorrLabel=('AK7Calo', inputJetCorrLabel),
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        doJetID=False
                 )

if isMC == False :
	addJetCollection(process, 
        cms.InputTag('caVHPFJets',"fat","VH"),         # Jet collection; must be already in the event when patLayer0 sequence is executed
        'CAVHFat',"PF",
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
#        jetCorrLabel=('KT6PF',  inputJetCorrLabel),
        jetCorrLabel=('AK7PF',  inputJetCorrLabel),
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        doJetID = False
                 )

else :
        addJetCollection(process,
        cms.InputTag('caVHPFJets',"fat","VH"),         # Jet collection; must be already in the event when patLayer0 sequence is executed
        'CAVHFat',"PF",
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
#        jetCorrLabel=('KT6PF',  inputJetCorrLabel),
        jetCorrLabel=('AK7PF',  inputJetCorrLabel),
        doType1MET=False,
        genJetCollection=cms.InputTag("ak5GenJets"),
        doL1Cleaning=False,
        doL1Counters=False,
        doJetID = False
                 )

if isMC == False :
	addJetCollection(process,
        cms.InputTag('caVHPFJets',"sub","VH"),         # Jet collection; must be already in the event when patLayer0 sequence is executed
        'CAVHSub',"PF",
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
#        jetCorrLabel=('KT6PF',  inputJetCorrLabel),
        jetCorrLabel=('KT4PF',  inputJetCorrLabel),
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        doJetID = False
                 )
else :
        addJetCollection(process,
        cms.InputTag('caVHPFJets',"sub","VH"),         # Jet collection; must be already in the event when patLayer0 sequence is executed
        'CAVHSub',"PF",
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
#        jetCorrLabel=('KT6PF',  inputJetCorrLabel),
        jetCorrLabel=('KT4PF',  inputJetCorrLabel),
        doType1MET=False,
        genJetCollection=cms.InputTag('caVHGenJets',"sub","VH"),
        doL1Cleaning=False,
        doL1Counters=False,
        doJetID = False
                 )

if isMC == False :
	addJetCollection(process,
        cms.InputTag('caVHPFJets',"filter","VH"),         # Jet collection; must be already in the event when patLayer0 sequence is executed
        'CAVHFilter',"PF",
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
        jetCorrLabel=('KT4PF',  inputJetCorrLabel),
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        doJetID = False
                 )
else :
        addJetCollection(process,
        cms.InputTag('caVHPFJets',"filter","VH"),         # Jet collection; must be already in the event when patLayer0 sequence is executed
       'CAVHFilter',"PF",
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
        jetCorrLabel=('KT4PF',  inputJetCorrLabel),
        doType1MET=False,
        genJetCollection=cms.InputTag('caVHGenJets',"filter","VH"),
        doL1Cleaning=False,
        doL1Counters=False,
        doJetID = False
                 )



#addJetCollection(process,cms.InputTag("CAsubJetsProducer"),"subCA","PF",
#                 doJTA=True,doBTagging=True,jetCorrLabel=('KT4PF',  inputJetCorrLabel),
#                 doType1MET=False,doL1Cleaning = False,doL1Counters=False,
#                 doJetID = False)

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
process.patJetCorrFactors.rho = cms.InputTag("kt6PFJets", "rho")

from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')


# Jets
#setTagInfos(process,
#            coll = "allLayer1Jets",
#            tagInfos = cms.vstring("secondaryVertex")
#            )

for jetcoll in (process.selectedPatJetsCAVHFatCalo, process.selectedPatJetsCAVHFatPF, process.selectedPatJetsCAVHSubPF, process.selectedPatJetsCAVHFilterPF):
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
process.selectedPatJetsCAVHFatCalo.cut = cms.string('pt > 15. & abs(eta) < 5.0')
process.selectedPatJetsCAVHFatPF.cut = cms.string('pt > 15. & abs(eta) < 5.0')
process.selectedPatJetsCAVHSubPF.cut = cms.string('pt > 10. & abs(eta) < 5.0')
process.selectedPatJetsCAVHFilterPF.cut = cms.string('pt > 2. & abs(eta) < 5.0')
##process.selectedLayer1Muons.cut = cms.string('pt > 20. & abs(eta) < 2.5 & muonID("TMLastStationLoose")')
##process.selectedLayer1Electrons.cut = cms.string('pt > 20. & abs(eta) < 2.5 & electronID("eidLoose")')
process.selectedPatMuons.cut = cms.string('isGlobalMuon || (isTrackerMuon && muonID("TrackerMuonArbitrated"))')
process.selectedPatMuons.cut = cms.string('')
#process.selectedPatElectrons.cut = cms.string('')

process.selectedPatMuonsWithIso = process.selectedPatMuons.clone(
    src = 'selectedPatMuons',
     cut =  (
        "pt > 10 && isGlobalMuon && isTrackerMuon && globalTrack().normalizedChi2 < 10 &&" +
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
   "pt > 10.0 && abs(eta) < 2.5 &&"
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


process.cleanPatJetsAK5PF = cms.EDProducer("PATJetCleaner",					   src = cms.InputTag("selectedPatJetsAK5PF"),
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
#algorithm = 'ca'

##if algorithm == 'kt' :
#    process.caTopCaloJets.algorithm = cms.int32(0)
#    process.caTopGenJets.algorithm = cms.int32(0)
#    process.caTopPFJets.algorithm = cms.int32(0)
#elif algorithm == 'ca' :
#    process.caTopCaloJets.algorithm = cms.int32(1)
#    process.caTopGenJets.algorithm = cms.int32(1)
#    process.caTopPFJets.algorithm = cms.int32(1)    
#elif algorithm == 'antikt' :
#    process.caTopCaloJets.algorithm = cms.int32(2)
#    process.caTopGenJets.algorithm = cms.int32(2)
#    process.caTopPFJets.algorithm = cms.int32(2)
#else:
#    print "Error: algorithm '%s' unknown. Use one of kt, ca, antikt." % algorithm
#    raise AttributeError()
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
    "pt > 5.0 && abs(eta) < 2.5 &&"                               +
    "(isEE || isEB) && !isEBEEGap"                              
)

from JetMETCorrections.Type1MET.pfMETCorrections_cff import pfCandsNotInJet
process.pfCandsNotInJetForPFMEtSignCovMatrix = pfCandsNotInJet.clone()

from RecoMET.METProducers.METSigParams_cfi import *
process.pfMEtSignCovMatrix = cms.EDProducer("PFMEtSignCovMatrixProducer",
  METSignificance_params,
  src = cms.VInputTag(
      'selectedPatMuonsWithIso',
	  'selectedPatElectronsWithIso',
      'cleanPatJetsAK5PF',
      'pfCandsNotInJetForPFMEtSignCovMatrix'
  ),
  addJERcorr = cms.PSet(
      inputFileName = cms.FileInPath('PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'),
      lutName = cms.string('pfJetResolutionMCtoDataCorrLUT')
  )
)


process.HbbAnalyzerNew = cms.EDProducer("HbbAnalyzerNew",
    runOnMC = cms.bool(isMC),
    hltResultsTag = cms.InputTag("TriggerResults::HLT"),
    lep_ptCutForBjets = cms.double(5),
    electronNoCutsTag = cms.InputTag("pfAllElectrons"),
    electronTag = cms.InputTag("selectedElectronsMatched"),
#    electronTag = cms.InputTag("selectedPatElectrons"),
    tauTag = cms.InputTag("patTaus"),
#   muonTag = cms.InputTag("selectedPatMuons"),
    muonNoCutsTag = cms.InputTag("pfAllMuons"),
    muonTag = cms.InputTag("selectedMuonsMatched"),
    jetTag = cms.InputTag("selectedPatJetsCAVHFatPF"),
    subjetTag = cms.InputTag("selectedPatJetsCAVHSubPF"),
    filterjetTag = cms.InputTag("selectedPatJetsCAVHFilterPF"),
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

# trigger match for taus
process.selectedTausTriggerMatch = defaultTriggerMatch.clone(
        src         = cms.InputTag( "selectedPatTaus" ), 
        matchedCuts = cms.string('path("*Tau*",0,0)')

)
# trigger object embedders for the same collections
process.selectedTausMatched = cms.EDProducer( "PATTriggerMatchTauEmbedder",
        src     = cms.InputTag(  "selectedPatTaus" ),
        matches = cms.VInputTag( cms.InputTag('selectedTausTriggerMatch') )
    )



process.leptonTrigMatch = cms.Sequence (
	process.selectedElectronsTriggerMatch
	*process.selectedElectronsMatched
	*process.selectedMuonsTriggerMatch
	*process.selectedMuonsMatched
	*process.selectedTausTriggerMatch
	*process.selectedTausMatched
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
#
# for charged MET
#
process.pfNoPileUpCharge  = cms.EDFilter(
   "GenericPFCandidateSelector",
   src = cms.InputTag("pfNoPileUp"),
     cut = cms.string("charge !=0")
 )
 
process.pfMETNoPUCharge = process.pfMET.clone()
process.pfMETNoPUCharge.src=cms.InputTag("pfNoPileUpCharge")

# type 1 +2 MET corrected
#process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
from PhysicsTools.PatUtils.patPFMETCorrections_cff import *
process.patPFMet = patPFMet
process.pfCandsNotInJet = pfCandsNotInJet


process.selectedPatJetsForMETtype1p2Corr = selectedPatJetsForMETtype1p2Corr
process.selectedPatJetsForMETtype1p2Corr.src = cms.InputTag('selectedPatJetsAK5PF')
process.selectedPatJetsForMETtype2Corr = selectedPatJetsForMETtype2Corr
process.selectedPatJetsForMETtype2Corr.src = cms.InputTag('selectedPatJetsAK5PF')
process.patPFJetMETtype1p2Corr = patPFJetMETtype1p2Corr
process.patPFJetMETtype1p2Corr.type1JetPtThreshold = cms.double(10.0)
process.patPFJetMETtype1p2Corr.skipEM = cms.bool(False)
process.patPFJetMETtype1p2Corr.skipMuons = cms.bool(False)
if isMC == False :
	process.patPFJetMETtype1p2Corr.jetCorrLabel = cms.string("L2L3Residual")
 	process.patPFMet.addGenMET = cms.bool(False)
process.patPFJetMETtype2Corr= patPFJetMETtype2Corr
process.pfCandMETcorr = pfCandMETcorr
process.patType1CorrectedPFMet = patType1CorrectedPFMet
process.patType1p2CorrectedPFMet = patType1p2CorrectedPFMet




#--------------------------------------------------------------------------------
# define sequence to run all modules
process.producePatPFMETCorrections = cms.Sequence(
    process.patPFMet
   * process.kt6PFJets
   * process.ak5PFJets
   * process.pfCandsNotInJet
   * process.selectedPatJetsForMETtype1p2Corr
   * process.selectedPatJetsForMETtype2Corr 
   * process.patPFJetMETtype1p2Corr
   * process.patPFJetMETtype2Corr
   * process.pfCandMETcorr 
   * process.patType1CorrectedPFMet
   * process.patType1p2CorrectedPFMet
   
 
    
)
#--------------------------------------------------------------------------------

from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties

runMEtUncertainties(process, 'selectedPatElectrons', '', 'selectedPatMuons', 'selectedPatTaus' , 'selectedPatJetsAK5PF') 


### build type1/2 correction also on top of pfMETnoPU
process.patType1CorrectedPFMetNoPU = process.patType1CorrectedPFMet.clone()
from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import patMETs
process.patPFMetNoPU = patMETs.clone(
    metSource = cms.InputTag('pfMETNoPU'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue')
)

#process.patPFMetNoPU = process.patMets.clone()
#process.patPFMetNoPU.metSource = cms.InputTag('pfMETNoPU'),
process.patType1CorrectedPFMetNoPU.src = cms.InputTag('patPFMetNoPU')
process.patType1p2CorrectedPFMetNoPU = process.patType1p2CorrectedPFMet.clone()
process.patType1p2CorrectedPFMetNoPU.src = cms.InputTag('patPFMetNoPU')
from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties




process.dump = cms.EDAnalyzer("EventContentAnalyzer")

# drop the meta data for dropped data
#process.out.dropMetaData = cms.string("DROPPED")
process.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexing_cff')
#process.inclusiveMergedVertices = process.vertexMerger.clone()
#process.inclusiveMergedVertices.secondaryVertices = cms.InputTag("inclusiveVertices")
#process.inclusiveMergedVertices.maxFraction = 0.2
#process.inclusiveMergedVertices.minSignificance = cms.double(10.)

process.load("RecoBTag/SecondaryVertex/bVertexFilter_cfi")
process.selectedVertices = process.bVertexFilter.clone()
process.selectedVertices.secondaryVertices = cms.InputTag("inclusiveMergedVertices")
process.selectedVertices.minVertices = 0
process.selectedVertices.vertexFilter.multiplicityMin = 3



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

#saves leptons with pt > 5 GeV, b and c quarks, B and D hadrons (and Lambda b,c), Z,W,Higgs, all status 3 particles. Daugthers of Z,W,H.
process.savedGenParticles = cms.EDProducer(
    "GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *", # this is the default
    "keep++ pdgId >= 23 && pdgId <= 25", #keep W,Z,H and theirs products
    "keep++ pdgId == 22 && pt > 15", #keep gamma above 15 GeV
    "drop++   status == 2 ", #drop all non stable decay products (and daughters) [quarks are going to be added below]
    "keep++ abs(pdgId) == 15", #keep tau and its decay prods
    "keep  numberOfMothers() > 0 && abs(mother(0).pdgId) == 15", #keep first generation of tau daugthers (this is redundant I think)
    "drop  numberOfMothers() > 0 && abs(mother(0).pdgId) == {pi0}", #drop pi0 daugthers photons
    "keep  (abs(pdgId) ==13 || abs(pdgId) ==11 || abs(pdgId) ==15 ) &&  pt > 5.0", #keep leptons of decent pT
    "keep  (abs(pdgId) > 400 &&  abs(pdgId) < 600)    ||     (  (abs(pdgId) > 4000 &&  abs(pdgId) < 6000)  )",  # track-back the origin of B/D
    "keep  (  (abs(pdgId) >= 4 &&  abs(pdgId) <= 6)) ", #keep heavy quarks
    "keep ( status == 3)"  #keep event summary status3 (for pythia)
    )
)




if isMC == False :
        process.p = cms.Path(
		                process.softElectronCands*
                    process.goodOfflinePrimaryVertices*
                     process.inclusiveVertexing*
                     process.PF2PAT*
		    process.PFTau* # re-run PFTau sequence as per tau POG recommendation
		    process.pfCandsForIsolationSequence *
		     process.muonPFIsolationSequence *
		     process.electronPFIsolationSequence *
                     process.ak5CaloJets*
                     process.ak7CaloJets*
                     process.kt6PFJets*
		     process.kt6PFJets25* 
                     process.ak5PFJets*
                     process.ak7PFJets*
                     process.caVHCaloJets*
                     process.caVHPFJets*
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
		     process.pfNoPileUpCharge*
       		     process.pfMETNoPUCharge*
		    process.patPFMetNoPU *
		    process.patType1CorrectedPFMetNoPU *
		    process.patType1p2CorrectedPFMetNoPU *
#		     process.producePatPFMETCorrections *
		    process.leptonTrigMatch*
#                     process.inclusiveMergedVertices* #already included now in inclusiveVertexing sequence
		     process.selectedVertices*
                     process.bcandidates*
                     process.HbbAnalyzerNew*
					 			process.pfCandsNotInJetForPFMEtSignCovMatrix*
			process.pfMEtSignCovMatrix


#process.hbbCandidates*process.hbbHighestPtHiggsPt30Candidates*process.hbbBestCSVPt20Candidates
                     )
else :
        process.p = cms.Path(
		process.softElectronCands*
#		process.printTree *
#		    process.printList *
	             process.savedGenParticles *
#		    process.printSavedList *
		process.goodOfflinePrimaryVertices*
                     process.inclusiveVertexing*
                     process.genParticlesForJets*
                     process.ak5GenJets*
                     process.caVHGenJets*
                     process.PF2PAT*
		process.PFTau* # re-run PFTau sequence as per tau POG recommendation
		    process.pfCandsForIsolationSequence *
		     process.muonPFIsolationSequence *
		     process.electronPFIsolationSequence *
                     process.ak5CaloJets*
                     process.ak7CaloJets*
                     process.kt6PFJets*
	      	     process.kt6PFJets25* 
                     process.ak5PFJets*
                     process.ak7PFJets*
                     process.caVHCaloJets*
                     process.caVHPFJets*
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
		     process.pfNoPileUpCharge*
       		     process.pfMETNoPUCharge*
		  process.patPFMetNoPU *
		process.patType1CorrectedPFMetNoPU *
		process.patType1p2CorrectedPFMetNoPU *
		 #    process.producePatPFMETCorrections *

#		     process.dump*
                     process.leptonTrigMatch*
#                     process.inclusiveMergedVertices* #already included now in inclusiveVertexing sequence
                     process.selectedVertices*
                     process.bcandidates*
		     process.bhadrons*
                     process.HbbAnalyzerNew*
					 					 			process.pfCandsNotInJetForPFMEtSignCovMatrix*
			process.pfMEtSignCovMatrix
#process.hbbCandidates*process.hbbHighestPtHiggsPt30Candidates*process.hbbBestCSVPt20Candidates
                     )

### begin of  filters required by JETMET for  MET analysis 
process.hbhepath = cms.Path(process.HBHENoiseFilter)
process.ecalFilter = cms.Path(process.EcalDeadCellEventFilter)
process.cschaloFilter = cms.Path(process.CSCTightHaloFilter)
process.hcallaserFilter=cms.Path(process.hcalLaserEventFilter)
process.trackingfailureFilter = cms.Path( process.goodOfflinePrimaryVertices*process.ak5PFJetsL2L3Residual*process.trackingFailureFilter)
### end of  filters required by JETMET for  MET analysis 
	    
process.load('GeneratorInterface.GenFilters.TotalKinematicsFilter_cfi');
process.totalKinematics=cms.Path(process.totalKinematicsFilter)

#process.candidates = cms.Path(process.hbbCandidates*process.hbbHighestPtHiggsPt30Candidates*process.hbbBestCSVPt20Candidates)


#process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )

process.options = cms.untracked.PSet(
	wantSummary = cms.untracked.bool(True),
	Rethrow = cms.untracked.vstring('ProductNotFound')
	)
process.e = cms.EndPath(process.out1)


if isMC == False :
 process.schedule = cms.Schedule(process.p, process.hbhepath, process.nTuplizePath ,process.ecalFilter, process.cschaloFilter, process.hcallaserFilter, process.trackingfailureFilter,  process.e)
else :
 process.schedule = cms.Schedule(process.p, process.hbhepath, process.nTuplizePath ,process.ecalFilter, process.cschaloFilter, process.hcallaserFilter, process.trackingfailureFilter, process.totalKinematics,process.e)


#
temp = process.dumpPython()
outputfile = file("patexp.py",'w')
outputfile.write(temp)
outputfile.close()
