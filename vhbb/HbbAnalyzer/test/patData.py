import FWCore.ParameterSet.Config as cms
isMC = False

if isMC == False :
        inputJetCorrLabel = ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']
else :
        inputJetCorrLabel = ['L1FastJet','L2Relative', 'L3Absolute']

process = cms.Process("VH")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

## Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#"/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/193/575/D8D5A7EC-EF99-E111-AB1D-003048D2BE12.root",
#"/store/data/Run2012A/MuHad/AOD/PromptReco-v1/000/193/556/5E52D114-BA99-E111-AF69-5404A63886C1.root"
##	"/store/data/Run2012A/MultiJet/AOD/PromptReco-v1/000/193/556/D84A9140-AC99-E111-A808-001D09F28D4A.root"
	
	#"/store/data/Run2012A/MET/RECO/PromptReco-v1/000/190/782/A2D1E62F-8684-E111-B651-00237DDBE41A.root"
#"/store/data/Run2012C/SingleMu/AOD/PromptReco-v2/000/200/369/60F3A873-74E2-E111-8A04-003048D37538.root"
#"/store/data/Run2012C/MET/AOD/PromptReco-v1/000/197/770/6E668370-77C3-E111-991E-BCAEC5364C4C.root"
	'/store/data/Run2012C/MET/AOD/PromptReco-v1/000/198/913/4AE2FEC7-78CE-E111-9D10-485B3962633D.root',


 )
)


## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond

if  isMC  :
   process.GlobalTag.globaltag = cms.string("START53_V7E::All" )
else :
   process.GlobalTag.globaltag = cms.string("GR_P_V40_AN1::All")

process.load("Configuration.StandardSequences.MagneticField_cff")

## Test JEC from test instances of the global DB
#process.load("PhysicsTools.PatAlgos.patTestJEC_cfi")

## Test JEC from local sqlite file
#process.load("PhysicsTools.PatAlgos.patTestJEC_local_cfi")

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root'),
                               # save only events passing the full path
			       dropMetaData = cms.untracked.string('ALL'),
	                       splitLevel = cms.untracked.int32(99),
                               SelectEvents = cms.untracked.PSet(      SelectEvents = cms.vstring('p')     ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
                               outputCommands = cms.untracked.vstring('drop *' )
                               )

process.outpath = cms.EndPath(process.out)


# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")


# Configure PAT to use PF2PAT instead of AOD sources
# this function will modify the PAT sequences.
from PhysicsTools.PatAlgos.tools.pfTools import *



#PF2PAT
from PhysicsTools.PatAlgos.tools.pfTools import *

# An empty postfix means that only PF2PAT is run,
# otherwise both standard PAT and PF2PAT are run. In the latter case PF2PAT
# collections have standard names + postfix (e.g. patElectronPFlow)
postfix = "PFlow"
jetAlgo = "AK5"
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=isMC, postfix=postfix)
usePFIso( process )
switchJetCollection(process, jetCollection=cms.InputTag('pfJetsPFlow'),doJTA=True,doBTagging=True, jetCorrLabel=('AK5PFchs', inputJetCorrLabel),doType1MET=False, doJetID=False)

process.pfPileUpPFlow.Vertices = cms.InputTag('goodOfflinePrimaryVertices')
process.pfPileUpPFlow.checkClosestZVertex = cms.bool(False)


from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )

if not isMC:
    # removing MC matching for standard PAT sequence
    # for the PF2PAT+PAT sequence, it is done in the usePF2PAT function
    removeMCMatchingPF2PAT( process, '' )

# top projections in PF2PAT:
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = False #do old style cleaning
getattr(process,"pfNoElectron"+postfix).enable = False #do old style cleaning
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True


#Electron ID
process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
process.mvaID = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )


########## Hbb specific stuff starts here ########################
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
# we have already PF Met as default in PF2PAT (AR)

#A few generic options
process.patJets.addTagInfos  = True

# rho2.5 calculation
process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJetsForIsolation = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

process.kt6PFJets25 = process.kt4PFJets.clone( src = 'pfNoElectron'+postfix,rParam = 0.6,doRhoFastjet = True,Ghost_EtaMax = 2.5, Rho_EtaMax = 2.5 )
process.kt6PFJetsCentralNeutral = process.kt6PFJets.clone( src = cms.InputTag("pfAllNeutralHadronsAndPhotons"+postfix), Ghost_EtaMax = cms.double(3.1), Rho_EtaMax = cms.double(2.5), inputEtMin = cms.double(0.5) )

from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
process.ak7PFJets = ak5PFJets.clone( rParam = cms.double(0.7) )
addJetCollection(process,cms.InputTag("ak7PFJets"),"AK7","PF", doJTA=True,doBTagging=True,jetCorrLabel=('AK7PF',  inputJetCorrLabel), doType1MET=False,doL1Cleaning = False,doL1Counters=False,doJetID = False)



##### SUBJET stuff
process.load("RecoJets.Configuration.GenJetParticles_cff")
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.ca4GenJets = ca4GenJets
process.ca8GenJets = ca4GenJets.clone( rParam = cms.double(0.8) )
process.load('RecoJets.JetProducers.caSubjetFilterPFJets_cfi')
from RecoJets.JetProducers.caSubjetFilterPFJets_cfi import caSubjetFilterPFJets
process.caVHPFJets = caSubjetFilterPFJets.clone(src=cms.InputTag('pfNoElectron'+postfix),useAdjacency = cms.int32(0))
process.load('RecoJets.JetProducers.caSubjetFilterCaloJets_cfi')
from RecoJets.JetProducers.caSubjetFilterCaloJets_cfi import caSubjetFilterCaloJets
process.caVHCaloJets = caSubjetFilterCaloJets.clone(useAdjacency = cms.int32(0))
process.load('RecoJets.JetProducers.caSubjetFilterGenJets_cfi')
from RecoJets.JetProducers.caSubjetFilterGenJets_cfi import caSubjetFilterGenJets
process.caVHGenJets = caSubjetFilterGenJets.clone()

from PhysicsTools.PatAlgos.tools.jetTools import *
#calo subjets in pat
#addJetCollection(process, cms.InputTag('caVHCaloJets:fat'),"CAVHFat","Calo",doJTA=True,doBTagging=True, jetCorrLabel=('AK5Calo', inputJetCorrLabel),doType1MET=False, doL1Cleaning=False, doL1Counters=False,doJetID=False)
#addJetCollection(process, cms.InputTag('caVHCaloJets:sub'),"CAVHSub","Calo",doJTA=True,doBTagging=True, jetCorrLabel=('AK5Calo', inputJetCorrLabel),doType1MET=False, doL1Cleaning=False, doL1Counters=False,doJetID=False)
#addJetCollection(process, cms.InputTag('caVHCaloJets:filter'),"CAVHFilter","Calo",doJTA=True,doBTagging=True, jetCorrLabel=('AK5Calo', inputJetCorrLabel),doType1MET=False, doL1Cleaning=False, doL1Counters=False,doJetID=False)
addJetCollection(process, cms.InputTag('caVHPFJets:fat'),"CAVHFat","PF",doJTA=True,doBTagging=True, jetCorrLabel=('AK5PF', inputJetCorrLabel),doType1MET=False, doL1Cleaning=False, doL1Counters=False,doJetID=False)
if isMC :
   addJetCollection(process, cms.InputTag('caVHPFJets:sub'),"CAVHSub","PF",doJTA=True,doBTagging=True, jetCorrLabel=('AK5PF', inputJetCorrLabel),doType1MET=False, doL1Cleaning=False, doL1Counters=False,doJetID=False, genJetCollection=cms.InputTag('caVHGenJets',"sub"))
   addJetCollection(process, cms.InputTag('caVHPFJets:filter'),"CAVHFilter","PF",doJTA=True,doBTagging=True, jetCorrLabel=('AK5PF', inputJetCorrLabel),doType1MET=False, doL1Cleaning=False, doL1Counters=False,doJetID=False, genJetCollection=cms.InputTag('caVHGenJets',"filter"))
else :
   addJetCollection(process, cms.InputTag('caVHPFJets:sub'),"CAVHSub","PF",doJTA=True,doBTagging=True, jetCorrLabel=('AK5PF', inputJetCorrLabel),doType1MET=False, doL1Cleaning=False, doL1Counters=False,doJetID=False)
   addJetCollection(process, cms.InputTag('caVHPFJets:filter'),"CAVHFilter","PF",doJTA=True,doBTagging=True, jetCorrLabel=('AK5PF', inputJetCorrLabel),doType1MET=False, doL1Cleaning=False, doL1Counters=False,doJetID=False)

# Place appropriate jet cuts (NB: no cut on number of constituents)
defaultJetCut = cms.string('pt > 15. & abs(eta) < 5.0')
defaultFatJetCut = cms.string('pt > 100. & abs(eta) < 5.0')
process.selectedPatJets.cut = defaultJetCut
process.selectedPatJetsAK7PF.cut = defaultJetCut
#process.selectedPatJetsCAVHFatCalo.cut = defaultFatJetCut
process.selectedPatJetsCAVHFatPF.cut = defaultFatJetCut
process.selectedPatJetsCAVHSubPF.cut = cms.string('pt > 15. & abs(eta) < 5.0')
process.selectedPatJetsCAVHFilterPF.cut = cms.string('pt > 5. & abs(eta) < 5.0')

#load btag SF 
process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")
process.load ("RecoBTag.PerformanceDB.BTagPerformanceDB1107")





#do cleaning only against isolatted and ID'ed leptons (TODO: review the electron selection)

process.cleanPatJets.checkOverlaps.muons.requireNoOverlaps = cms.bool(True)
process.cleanPatJets.checkOverlaps.electrons.requireNoOverlaps = cms.bool(True)

process.cleanPatJets.checkOverlaps.muons.preselection = ("pt > 20 && isGlobalMuon && globalTrack().normalizedChi2 < 10 && isPFMuon && "
            "innerTrack().hitPattern().trackerLayersWithMeasurement > 5 && "
            "innerTrack().hitPattern().numberOfValidPixelHits > 0 && "                        
            "globalTrack().hitPattern().numberOfValidMuonHits > 0 && "                         
            "numberOfMatchedStations > 1 && " 
            "dB < 0.2 && abs(eta) < 2.4 "                                                                     

           "&& ( chargedHadronIso + neutralHadronIso + photonIso ) < 0.10 * pt  "                                               
)

process.cleanPatJets.checkOverlaps.electrons.preselection =  (  "pt > 15.0 && abs(eta) < 2.5 &&"
   "(isEE || isEB) && !isEBEEGap &&"
  " (chargedHadronIso + neutralHadronIso + photonIso)/pt <0.10 &&"
   "dB < 0.02 && "  #dB is computed wrt PV but is transverse only, no info about dZ(vertex) 
   "( "
   "(isEE && ("
   "abs(deltaEtaSuperClusterTrackAtVtx) < 0.005 &&  abs(deltaPhiSuperClusterTrackAtVtx) < 0.02 && sigmaIetaIeta < 0.03 && hadronicOverEm < 0.10 &&  abs(1./ecalEnergy*(1.-eSuperClusterOverP)) < 0.05 "
   ")) || " 
    "(isEB && (  "
    "abs(deltaEtaSuperClusterTrackAtVtx) < 0.004 &&  abs(deltaPhiSuperClusterTrackAtVtx) < 0.03 && sigmaIetaIeta < 0.01 && hadronicOverEm < 0.12 && abs(1./ecalEnergy*(1.-eSuperClusterOverP)) < 0.05"
    "))"
#or use mvaNonTrigV0 and mvaTrigV0
    ")" 
)
#   "ecalDrivenSeed==1 && (abs(superCluster.eta)<2.5)"
 #  " && !(1.4442<abs(superCluster.eta)<1.566)"
 #  " && (ecalEnergy*sin(superClusterPosition.theta)>20.0)"
 #  " && (gsfTrack.trackerExpectedHitsInner.numberOfHits == 0)"
#    " && ((chargedHadronIso + neutralHadronIso + photonIso < 0.15 * pt)) "
#   " && ((isEB"
#   " && (sigmaIetaIeta<0.01)"
#   " && ( -0.8<deltaPhiSuperClusterTrackAtVtx<0.8 )"
#   " && ( -0.007<deltaEtaSuperClusterTrackAtVtx<0.007 )"
#   ")"
#   " || (isEE"
#   " && (sigmaIetaIeta<0.03)"
#   " && ( -0.7<deltaPhiSuperClusterTrackAtVtx<0.7 )"
#   " && ( -0.01<deltaEtaSuperClusterTrackAtVtx<0.01 )"
#   "))")


# use cone 0.3 for electron isolation rather than 0.4 (0.3 seem to be recommended by EGamma in the cut based ID twiki page)

process.patElectrons.isolationValues = cms.PSet(
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral04PFIdPFIso"),
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll04PFIdPFIso"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU04PFIdPFIso"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma04PFIdPFIso"),
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged04PFIdPFIso")
    )
process.patElectrons.isolationValuesNoPFId = cms.PSet(
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral04NoPFIdPFIso"),
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll04NoPFIdPFIso"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU04NoPFIdPFIso"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma04NoPFIdPFIso"),
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged04NoPFIdPFIso")
    )


process.patMuonsPFlow.isolationValues.user = cms.VInputTag(
        cms.InputTag("muPFIsoValueChargedAll03"+postfix),
        cms.InputTag("muPFIsoValueCharged03"+postfix),
        cms.InputTag("muPFIsoValueNeutral03"+postfix),
        cms.InputTag("muPFIsoValueGamma03"+postfix),
        cms.InputTag("muPFIsoValuePU03"+postfix),
        cms.InputTag("muPFIsoValuePU04"+postfix),
        )


#sourceElectrons = 'gsfElectrons'
#process.elPFIsoDepositCharged.src = sourceElectrons
#process.elPFIsoDepositChargedAll.src = sourceElectrons
#process.elPFIsoDepositNeutral.src = sourceElectrons
#process.elPFIsoDepositGamma.src = sourceElectrons
#process.elPFIsoDepositPU.src = sourceElectrons
#process.elePFIso = cms.Sequence (process.elPFIsoDepositCharged * process.elPFIsoDepositChargedAll * process.elPFIsoDepositNeutral * process.elPFIsoDepositGamma * process.elPFIsoDepositPU * process.pfElectronIsolationSequence )


#process.patElectrons.isolationValues = cms.PSet()
#process.patElectrons.isolationValues.user = cms.VInputTag(
#        cms.InputTag("elPFIsoValueChargedAll03PFId"),
#        cms.InputTag("elPFIsoValueCharged03PFId"),
#        cms.InputTag("elPFIsoValueNeutral03PFId"),
#        cms.InputTag("elPFIsoValueGamma03PFId"),
#       cms.InputTag("elPFIsoValuePU03PFId"),
#        cms.InputTag("elPFIsoValuePU04PFId"),
#
#        )

process.selectedPatElectrons.cut = (
    "(ecalDrivenSeed==1) &&"                                       +
    "pt > 5.0 && abs(eta) < 2.5 &&"                               +
    "(isEE || isEB) && !isEBEEGap"
)

######## MET Correction TODO: fix it!

# type 1 +2 MET corrected
process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
#--------------------------------------------------------------------------------

from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties

runMEtUncertainties(process, 'selectedPatElectrons', '', 'selectedPatMuons', 'selectedPatTaus' , 'selectedPatJets')







process.HbbAnalyzerNew = cms.EDProducer("HbbAnalyzerNew",
    runOnMC = cms.bool(isMC),
    hltResultsTag = cms.InputTag("TriggerResults::HLT"),
    lep_ptCutForBjets = cms.double(5),
    electronNoCutsTag = cms.InputTag("gsfElectrons"),
#    electronTag = cms.InputTag("selectedElectronsMatched"),
    electronTag = cms.InputTag("selectedPatElectrons"),
    tauTag = cms.InputTag("patTaus"),

    muonNoCutsTag = cms.InputTag("muons"),
   muonTag = cms.InputTag("selectedPatMuons"),
#    muonTag = cms.InputTag("selectedMuonsMatched"),
    jetTag = cms.InputTag("selectedPatJetsCAVHFatPF"),
    subjetTag = cms.InputTag("selectedPatJetsCAVHSubPF"),
    filterjetTag = cms.InputTag("selectedPatJetsCAVHFilterPF"),
    simplejet2Tag = cms.InputTag("cleanPatJets"),
    simplejet3Tag = cms.InputTag("selectedPatJetsAK7PF"),
    photonTag = cms.InputTag("selectedPatPhotons"),
    metTag = cms.InputTag("met"), #this input tag is used to fill calo MET 
    verbose = cms.untracked.bool(False),

    #TODO: clean up the analyzer
    simplejet1Tag = cms.InputTag("UNUSED_WAS_selectedPatJets"),
    simplejet4Tag = cms.InputTag("UNUSED_WAS_selectedPatJetsAK7Calo"),

)

## TODO review this and modify the HbbAnalyzer
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
    #MVA
    mvaTrigV0 = cms.InputTag("mvaTrigV0"),
    mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0"),
    #old ones
    eidVBTFRel95 = cms.InputTag("eidVBTFRel95"),
    eidVBTFRel85 = cms.InputTag("eidVBTFRel85"),
    eidVBTFRel80 = cms.InputTag("eidVBTFRel80"),
    eidVBTFRel70 = cms.InputTag("eidVBTFRel70"),
    eidVBTFCom95 = cms.InputTag("eidVBTFCom95"),
    eidVBTFCom85 = cms.InputTag("eidVBTFCom85"),
    eidVBTFCom80 = cms.InputTag("eidVBTFCom80"),
    eidVBTFCom70 = cms.InputTag("eidVBTFCom70")
)   


## Special METs
process.patMETsHT = cms.EDProducer("MHTProducer",
  JetCollection = cms.InputTag("patJets"),
  MinJetPt      = cms.double(30),
  MaxJetEta     = cms.double(5)
)
process.patPFMetNoPU = process.patMETs.clone(
    metSource = cms.InputTag('pfMETNoPU'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue')
)

process.pfMETNoPU = process.pfMETPFlow.clone()
process.pfMETNoPU.src=cms.InputTag("pfNoPileUp"+postfix)
#rocess.pfMETNoPU.jets = cms.InputTag("pfJets"+postfix)


 
process.pfNoPileUpCharge  = cms.EDFilter(
   "GenericPFCandidateSelector",
   src = cms.InputTag("pfNoPileUp"+postfix),
     cut = cms.string("charge!=0" )
)
process.pfMETNoPUCharge = process.pfMET.clone()
process.pfMETNoPUCharge.src=cms.InputTag("pfNoPileUpCharge")
process.pfMETNoPUCharge.calculateSignificance = cms.bool(False)





# load the PU JetID sequence
process.load("CMGTools.External.pujetidsequence_cff")
process.puJetId.jets =  cms.InputTag("cleanPatJets")
process.puJetMva.jets =  cms.InputTag("cleanPatJets")



### B Hadron truth
process.bhadrons = cms.EDProducer('MCBHadronProducer',
                                 quarkId = cms.uint32(5)
                                 )

### Save some gen particles
# saves leptons with pt > 5 GeV, b and c quarks, B and D hadrons (and Lambda b,c), Z,W,Higgs, all status 3 particles. Daugthers of Z,W,H.
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

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

## IVF and BCandidate producer for Vbb cross check analysis
process.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexing_cff')
process.load("RecoBTag/SecondaryVertex/bVertexFilter_cfi")
process.selectedVertices = process.bVertexFilter.clone()
process.selectedVertices.secondaryVertices = cms.InputTag("inclusiveMergedVertices")
process.selectedVertices.minVertices = 0
process.selectedVertices.vertexFilter.multiplicityMin = 3
process.bcandidates = cms.EDProducer('BCandidateProducer',
                                     src = cms.InputTag('selectedVertices','',''),
                                     primaryVertices =
                                     cms.InputTag('offlinePrimaryVerticesWithBS','',''),
                                     minDRUnique = cms.untracked.double(0.4),
                                     minvecSumIMifsmallDRUnique = cms.untracked.double(5.5),
                                     minCosPAtomerge = cms.untracked.double(0.99),
                                     maxPtreltomerge = cms.untracked.double(7777.0)
                                     )

process.bhadrons = cms.EDProducer('MCBHadronProducer',
                                 quarkId = cms.uint32(5)
                                 )

process.gen = cms.Sequence(  process.genParticlesForJets* process.caVHGenJets * process.bhadrons * process.savedGenParticles)
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )
process.common = cms.Sequence( 
    process.goodOfflinePrimaryVertices *
    process.softElectronCands*
    process.inclusiveVertexing*
    process.eidSequence*
    process.patPF2PATSequencePFlow *
    process.pfMETNoPU*
    process.pfNoPileUpCharge*
    process.pfMETNoPUCharge*
    process.kt6PFJetsCentralNeutral *
    process.kt6PFJetsForIsolation *
    process.kt6PFJets25*
    process.ak7PFJets*
    process.caVHCaloJets*
    process.caVHPFJets*
    process.mvaID *
    process.patDefaultSequence*
   process.patMETsHT*
    process.patPFMetNoPU *
#                    process.patType1CorrectedPFMetNoPU *
#                    process.patType1p2CorrectedPFMetNoPU *
#                    process.leptonTrigMatch*
    process.selectedVertices*
    process.bcandidates*
    process.producePatPFMETCorrections *
    process.puJetIdSqeuence *   ### it is not a typo ;-)
    process.HbbAnalyzerNew
)

if isMC : 
   process.p = cms.Path(process.gen * process.common)
else :
   process.p = cms.Path(process.common)

########## Hbb specific finishes here #######################################
 
############# MET Filter flags
### from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/RecoMET/METFilters/test/exampleICHEPrecommendation_cfg.py?revision=1.1&view=markup&pathrev=V00-00-08
## The good primary vertex filter ____________________________________________||
process.primaryVertexFilter = cms.EDFilter(
	"VertexSelector",
	src = cms.InputTag("offlinePrimaryVertices"),
	cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
	filter = cms.bool(True)
	)

## The beam scraping filter __________________________________________________||
process.noscraping = cms.EDFilter(
	"FilterOutScraping",
	applyfilter = cms.untracked.bool(True),
	debugOn = cms.untracked.bool(False),
	numtrack = cms.untracked.uint32(10),
	thresh = cms.untracked.double(0.25)
	)

## The iso-based HBHE noise filter ___________________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')

## The CSC beam halo tight filter ____________________________________________||
process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')

## The HCAL laser filter _____________________________________________________||
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
process.hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
process.hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)

## The ECAL dead cell trigger primitive filter _______________________________||
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
## For AOD and RECO recommendation to use recovered rechits
process.EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")

## The EE bad SuperCrystal filter ____________________________________________||
process.load('RecoMET.METFilters.eeBadScFilter_cfi')

## The Good vertices collection needed by the tracking failure filter ________||
process.goodVertices = cms.EDFilter(
	"VertexSelector",
	filter = cms.bool(False),
	src = cms.InputTag("offlinePrimaryVertices"),
	cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
	)

   ## The tracking failure filter _______________________________________________||
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')

process.hbhepath = cms.Path(process.HBHENoiseFilter)
process.ecalFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter )
process.cschaloFilter = cms.Path(process.CSCTightHaloFilter)
process.hcallaserFilter=cms.Path(process.hcalLaserEventFilter)
process.trackingfailureFilter = cms.Path(    process.goodVertices * process.trackingFailureFilter )
process.eebadscFilter= cms.Path(  process.eeBadScFilter)  					      
##### END of filters required by JETMET for  MET analysis 



#   process.GlobalTag.globaltag =  ...    ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
process.out.fileName = 'PAT.edm.root'
process.out.outputCommands = cms.untracked.vstring(
        'drop *',
	"keep *_puJetId_*_*", # input variables
	"keep *_puJetMva_*_*", # final MVAs and working point flags
        'keep *_savedGenParticles_*_*',
        'keep *_HbbAnalyzerNew_*_*',
        'keep VHbbCandidates_*_*_*',
#       'keep PileupSummaryInfos_*_*_*',
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
        "keep *_EcalDeadCellEventFilter_*_*",
        "keep *_patType1CorrectedPFMet*_*_*",
        "keep *_patType1p2CorrectedPFMet*_*_*",
        "keep *_patMETsHT*_*_*",
        "keep LHEEventProduct_source_*_LHE",
       'keep patTriggerAlgorithms_patTrigger_*_*',
       'keep patTriggerConditions_patTrigger_*_*',
       'keep patTriggerObjects_patTrigger_*_*',
       'keep patTriggerFilters_patTrigger_*_*',
       'keep patTriggerPaths_patTrigger_*_*',
       'keep *_patTriggerEvent_*_*'
	
        )


temp = process.dumpPython()
outputfile = file("patexpstd.py",'w')
outputfile.write(temp)
outputfile.close()

