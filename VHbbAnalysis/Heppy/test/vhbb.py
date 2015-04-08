#! /usr/bin/env python
fastObjects=True

#Switch to True to produce x1, x2, id1, id2, pdf scale
doPDFVars = False

import ROOT
from DataFormats.FWLite import *
import PhysicsTools.HeppyCore.framework.config as cfg
from VHbbAnalysis.Heppy.vhbbobj import *
from PhysicsTools.HeppyCore.utils.deltar import deltaPhi
from PhysicsTools.Heppy.analyzers.core.AutoFillTreeProducer  import * 

import logging
logging.basicConfig(level=logging.ERROR)


cfg.Analyzer.nosubdir = True

treeProducer= cfg.Analyzer(
	class_object=AutoFillTreeProducer, 
	verbose=False, 
	vectorTree = True,
        globalVariables	= [
		 NTupleVariable("Vtype", lambda ev : ev.Vtype, help="Event classification"),
		 NTupleVariable("VtypeSim", lambda ev : ev.VtypeSim, help="Event classification",mcOnly=True),
		 NTupleVariable("VMt", lambda ev : ev.V.goodMt, help="Transverse mass of the vector boson"),
		 NTupleVariable("HVdPhi", lambda ev : deltaPhi(ev.V.phi(),ev.H.phi()), help="Delta phi between Higgs and Z/W"),
		 NTupleVariable("fakeMET_sumet", lambda ev : ev.fakeMET.sumet, help="Fake SumET from Zmumu events removing muons"),
		 NTupleVariable("rho",  lambda ev: ev.rho, float, help="kt6PFJets rho"),
		 NTupleVariable("deltaR_jj",  lambda ev: deltaR(ev.hJets[0].eta(),ev.hJets[0].phi(),ev.hJets[1].eta(),ev.hJets[1].phi()) if len(ev.hJets) > 1 else -1, float, help="deltaR higgsJets"),
                 NTupleVariable("minDr3",    lambda ev: ev.minDr3, help="dR of closest jets for 3 jest case"),
                 NTupleVariable("lheNj",  lambda ev: ev.lheNj, float,mcOnly=True, help="number of jets at LHE level"),
                 NTupleVariable("lheNb",  lambda ev: ev.lheNb, float,mcOnly=True, help="number of b-jets at LHE level"),
                 NTupleVariable("lheNc",  lambda ev: ev.lheNc, float,mcOnly=True, help="number of c-jets at LHE level"),
                 NTupleVariable("lheNg",  lambda ev: ev.lheNg, float,mcOnly=True, help="number of gluon jets at LHE level"),
                 NTupleVariable("lheNl",  lambda ev: ev.lheNl, float,mcOnly=True, help="number of light(uds) jets at LHE level"),
		 NTupleVariable("lheV_pt",  lambda ev: ev.lheV_pt, float,mcOnly=True, help="Vector pT at LHE level"),
                 NTupleVariable("lheHT",  lambda ev: ev.lheHT, float,mcOnly=True, help="HT at LHE level"),
                 NTupleVariable("genTTHtoTauTauDecayMode", lambda ev: ev.genTTHtoTauTauDecayMode, int,mcOnly=True, help="gen level ttH, H -> tautau decay mode"),        
                 NTupleVariable("totSoftActivityJets", lambda ev: len([ x for x in ev.softActivityJets if x.pt()> 2 ] ), int, help="number of jets from soft activity with pt>2Gev"),
	],
	globalObjects = {
          "met"    : NTupleObject("met",     metType, help="PF E_{T}^{miss}, after default type 1 corrections"),
          "fakeMET"    : NTupleObject("fakeMET", fourVectorType, help="fake MET in Zmumu event obtained removing the muons"),
          "H"    : NTupleObject("H", fourVectorType, help="higgs"),
          "HCSV"    : NTupleObject("HCSV", fourVectorType, help="higgs CSV selection"),
          "H3cj"    : NTupleObject("H3cj", fourVectorType, help="higgs 3 cen jets selection"),
          "V"    : NTupleObject("V", fourVectorType, help="z or w"),
        },
	collections = {
		#standard dumping of objects
   	        "selectedLeptons" : NTupleCollection("selLeptons", leptonTypeVHbb, 8, help="Leptons after the preselection"),
#   	        "inclusiveLeptons" : NTupleCollection("incLeptons", leptonTypeVHbb, 8, help="Leptons after the preselection"),
		#old style stuff, can be removed at some point
   	        "vLeptons" : NTupleCollection("vLeptons", leptonTypeVHbb, 8, help="Leptons after the preselection"),
   	        "aLeptons" : NTupleCollection("aLeptons", leptonTypeVHbb, 8, help="Additional leptons, not passing the preselection"),
# now store only indices, this lines are left commented for possible debugging
#	        "hJets"       : NTupleCollection("hJets",     jetTypeVHbb, 8, sortDescendingBy = lambda jet : jet.btag('combinedSecondaryVertexBJetTags'),help="Higgs jets"),
#	        "aJets"       : NTupleCollection("aJets",     jetTypeVHbb, 8, sortDescendingBy = lambda jet : jet.btag('combinedSecondaryVertexBJetTags'),help="Additional jets"),

                "hjidx"       : NTupleCollection("hJidx",    objectInt, 2,help="Higgs jet indices"),
                "hjidxDiJetPtByCSV"       : NTupleCollection("hJidx_sortcsv",    objectInt, 2,help="Higgs jet indices within hJets with CSV sorting "),
                "ajidx"       : NTupleCollection("aJidx",    objectInt, 8,help="additional jet indices"),
                "hjidxCSV"       : NTupleCollection("hJCidx",    objectInt, 2,help="Higgs jet indices CSV"),
                "ajidxCSV"       : NTupleCollection("aJCidx",    objectInt, 8,help="additional jet indices CSV"),
                "hjidx3cj"       : NTupleCollection("hJ3Cidx",    objectInt, 3,help="Higgs jet indices 3 cen jets"),
                "ajidx3cj"       : NTupleCollection("aJ3Cidx",    objectInt, 8,help="additional jet indices 3 cen  jets"),

                "cleanJetsAll"       : NTupleCollection("Jet",     jetTypeVHbb, 15, help="Cental+fwd jets after full selection and cleaning, sorted by b-tag"),
                "inclusiveTaus"  : NTupleCollection("TauGood", tauTypeVHbb, 25, help="Taus after the preselection"),
                "softActivityJets"    : NTupleCollection("softActivityJets", fourVectorType, 5, help="jets made for soft activity"),
                "goodVertices"    : NTupleCollection("primaryVertices", primaryVertexType, 4, help="first four PVs"),

		#dump of gen objects
                "genJets"    : NTupleCollection("GenJet",   genParticleType, 15, help="Generated top quarks from hard scattering",filter=lambda x: x.pt() > 20,mcOnly=True),
                "gentopquarks"    : NTupleCollection("GenTop",     genParticleType, 4, help="Generated top quarks from hard scattering"),
                "gennusFromTop"    : NTupleCollection("GenNuFromTop",     genParticleType, 4, help="Generated neutrino from t->W decay"),
                "genbquarksFromH"      : NTupleCollection("GenBQuarkFromH",  genParticleType, 4, help="Generated bottom quarks from Higgs decays"),
                "genbquarksFromTop"      : NTupleCollection("GenBQuarkFromTop",  genParticleType, 4, help="Generated bottom quarks from top decays"),
                "genbquarksFromHafterISR"      : NTupleCollection("GenBQuarkFromHafterISR",  genParticleType, 4, help="Generated bottom quarks from Higgs decays"),
                "genwzquarks"     : NTupleCollection("GenWZQuark",   genParticleType, 6, help="Generated quarks from W/Z decays"),
                "genleps"         : NTupleCollection("GenLep",     genParticleType, 4, help="Generated leptons from W/Z decays"),
                "genlepsFromTop"         : NTupleCollection("GenLepFromTop",     genParticleType, 4, help="Generated leptons from t->W decays"),
                "gentauleps"      : NTupleCollection("GenLepFromTau", genParticleType, 6, help="Generated leptons from decays of taus from W/Z decays"),
		"genHiggsBosons"   : NTupleCollection("GenHiggsBoson", genParticleType, 4, help="Generated Higgs boson "),
		#"genZbosonsToLL"  : NTupleCollection("GenZbosonsToLL", genParticleType, 6, help="Generated W or Z bosons decaying to LL"),
		#"genWbosonsToLL"  : NTupleCollection("GenWbosonsToLL", genParticleType, 6, help="Generated W or Z bosons decaying to LL"),
		"genvbosons"       : NTupleCollection("GenVbosons", genParticleType, 6, help="Generated W or Z bosons, mass > 30"),
              
	}
	)

#Create shifted MET Ntuples
metNames={y:x for x,y in ROOT.pat.MET.__dict__.items() if y >= 0 and y < 13 and (x[-2:]=="Up" or x[-4:]=="Down")}
shifted_met_keys = ["met_shifted_{0}".format(n) for n in range(12)] #we do not need noShift I gueess
shifted_met_names = ["met_shifted_%s"%metNames[n] for n in range(12)] #we do not need noShift I gueess
shifted_mets = {mk: NTupleObject(nm, shiftedMetType, help="PF E_{T}^{miss}, after default type 1 corrections, shifted with %s" %mk) for mk,nm in zip(shifted_met_keys,shifted_met_names)}
treeProducer.globalObjects.update(shifted_mets)

# Lepton Analyzer, take its default config
from PhysicsTools.Heppy.analyzers.objects.LeptonAnalyzer import LeptonAnalyzer
LepAna = LeptonAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.objects.VertexAnalyzer import VertexAnalyzer
VertexAna = VertexAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.objects.PhotonAnalyzer import PhotonAnalyzer
PhoAna = PhotonAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.objects.TauAnalyzer import TauAnalyzer
TauAna = TauAnalyzer.defaultConfig
TauAna.inclusive_ptMin = 18.
TauAna.inclusive_etaMax = 2.5
TauAna.inclusive_dxyMax = 1000.
TauAna.inclusive_dzMax = 0.4
TauAna.inclusive_vetoLeptons = False
TauAna.inclusive_leptonVetoDR = 0.4
TauAna.inclusive_decayModeID = "decayModeFindingNewDMs"
TauAna.inclusive_tauID = "decayModeFindingNewDMs"
TauAna.inclusive_vetoLeptonsPOG = False
TauAna.inclusive_tauAntiMuonID = ""
TauAna.inclusive_tauAntiElectronID = ""
TauAna.inclusive_tauLooseID = "decayModeFindingNewDMs"

from PhysicsTools.Heppy.analyzers.objects.JetAnalyzer import JetAnalyzer
JetAna = JetAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.gen.LHEAnalyzer import LHEAnalyzer 
LHEAna = LHEAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.gen.GeneratorAnalyzer import GeneratorAnalyzer 
GenAna = GeneratorAnalyzer.defaultConfig
from VHbbAnalysis.Heppy.VHGeneratorAnalyzer import GeneratorAnalyzer as  VHGeneratorAnalyzer
VHGenAna = VHGeneratorAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.objects.METAnalyzer import METAnalyzer
METAna = METAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.core.PileUpAnalyzer import PileUpAnalyzer
PUAna = PileUpAnalyzer.defaultConfig

from VHbbAnalysis.Heppy.VHbbAnalyzer import VHbbAnalyzer
JetAna.jetPt = 15
JetAna.doQG=True
JetAna.QGpath="pdfQG_AK4chs_antib_13TeV_v1.root"
JetAna.recalibrateJets=True
JetAna.jecPath="jec"
JetAna.mcGT="PHYS14_V4_MC"

VHbb = cfg.Analyzer(
    verbose=False,
    class_object=VHbbAnalyzer,
    wEleSelection = lambda x : x.pt() > 25 and x.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-tight"),
    wMuSelection = lambda x : x.pt() > 25 and x.muonID("POG_ID_Tight"),
    zEleSelection = lambda x : x.pt() > 15 and x.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-loose"),
    zMuSelection = lambda x : x.pt() > 10 and x.muonID("POG_ID_Loose"),
    zLeadingElePt = 20,
    zLeadingMuPt = 20,
    higgsJetsPreSelection = lambda x:  x.puJetId() > 0 and x.jetID('POG_PFID_Loose') and x.pt() >  15 ,
    passall=False,
    doSoftActivity=True
)

from VHbbAnalysis.Heppy.TTHtoTauTauAnalyzer import TTHtoTauTauAnalyzer
TTHtoTauTau = cfg.Analyzer(
    verbose = False,
    class_object = TTHtoTauTauAnalyzer,
)
from VHbbAnalysis.Heppy.TTHtoTauTauGeneratorAnalyzer import TTHtoTauTauGeneratorAnalyzer
TTHtoTauTauGen = cfg.Analyzer(
    verbose = False,
    class_object = TTHtoTauTauGeneratorAnalyzer,
)

#from VHbbAnalysis.Heppy.HeppyShell import HeppyShell
#sh = cfg.Analyzer( class_object=HeppyShell)

from PhysicsTools.Heppy.analyzers.core.TriggerBitAnalyzer import TriggerBitAnalyzer
TrigAna = cfg.Analyzer(
    verbose = False,
    class_object = TriggerBitAnalyzer,
    triggerBits = {
        "METBTAG" : [ "HLT_PFMET120_NoiseCleaned_BTagCSV07_v*" ],
        "MET"     : [ "HLT_PFMET170_NoiseCleaned_v*" ],
        "DIELE"   : [ "HLT_Ele23_Ele12_CaloId_TrackId_Iso_v*" ],
        "ELE"     : [ "HLT_Ele32_eta2p1_WP85_Gsf_v*", "HLT_Ele32_eta2p1_WP85_Gsf_v*" ],
        "DIMU"    : [ "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*" ],
        "MU"      : [ "HLT_IsoTkMu24_eta2p1_IterTrk02_v*", "HLT_IsoTkMu24_IterTrk02_v*" ],
        "TAU"     : ["HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v*"],
        "HH4bQuad":["HLT_QuadJet45_TripleCSV0p5_v*"],
        "HH4bDouble":["HLT_DoubleJet90_Double30_TripleCSV0p5_v*"],  
        "MUE"     : [ "HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v*" ],
        "EMU"     : [ "HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v*" ],
        "METTAU"  : [ "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v*" ],
        "ELETAU"  : [ "HLT_Ele22_eta2p1_WP85_Gsf_LooseIsoPFTau20_v*" ],
        "MUTAU"   : [ "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v*" ],
        "DiTAU"   : [ "HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v*" ],
   },
#   processName = 'HLT',
#   outprefix = 'HLT'
   )

from PhysicsTools.HeppyCore.framework.services.tfile import TFileService 
output_service = cfg.Service(
      TFileService,
      'outputfile',
      name="outputfile",
      fname='tree.root',
      option='recreate'
    )

from PhysicsTools.Heppy.analyzers.core.TriggerBitAnalyzer import TriggerBitAnalyzer
FlagsAna = TriggerBitAnalyzer.defaultEventFlagsConfig

from PhysicsTools.Heppy.analyzers.gen.PDFWeightsAnalyzer import PDFWeightsAnalyzer
PdfAna = cfg.Analyzer(PDFWeightsAnalyzer,
    PDFWeights = [],
    doPDFVars = doPDFVars
)

if doPDFVars:
    treeProducer.globalVariables += [
        NTupleVariable("pdf_x1",  lambda ev: ev.pdf_x1, float,mcOnly=True, help="PDF energy fraction of first parton"),
        NTupleVariable("pdf_x2",  lambda ev: ev.pdf_x2, float,mcOnly=True, help="PDF energy fraction of second parton"),
        NTupleVariable("pdf_id1",  lambda ev: ev.pdf_id1, int,mcOnly=True, help="PDF id of first parton"),
        NTupleVariable("pdf_id2",  lambda ev: ev.pdf_id2, int,mcOnly=True, help="PDF id of second parton"),
        NTupleVariable("pdf_scale",  lambda ev: ev.pdf_scale, float,mcOnly=True, help="PDF scale"),
    ]

#TrigAna.unrollbits=True

sequence = [LHEAna,FlagsAna, GenAna,VHGenAna,PUAna,TrigAna,VertexAna,LepAna,PhoAna,TauAna,JetAna,METAna,PdfAna,VHbb,TTHtoTauTau,TTHtoTauTauGen,treeProducer]#,sh]


from PhysicsTools.Heppy.utils.miniAodFiles import miniAodFiles
sample = cfg.MCComponent(
    files = [

#       'root://eoscms//eos/cms//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/14053C6D-AD09-E411-A748-00266CFFC9C4.root',
#       'root://eoscms//eos/cms//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/62ED6255-AE09-E411-97CB-00266CFFBF88.root',
#       'root://eoscms//eos/cms//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/72C26B45-AD09-E411-A77C-00266CFFBF80.root',
#       'root://eoscms//eos/cms//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/BAEE7255-AE09-E411-8F9F-00266CFFBF88.root',
#       'root://eoscms//eos/cms//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/D600138D-AD09-E411-917F-00266CFFBF88.root'
#'root://xrootd.ba.infn.it//store/mc/Phys14DR/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/12328AE8-796B-E411-9D32-002590A831B4.root'
#"32ABFE4A-916B-E411-B2FA-00266CFFBC60.root" #Hbb
#"04860BAA-B673-E411-8B20-002481E0D50C.root" #DY 600
"TTPU20-007B37D4-8B70-E411-BC2D-0025905A6066.root" # 
#"root://xrootd.ba.infn.it//store/mc/Phys14DR/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/32ABFE4A-916B-E411-B2FA-00266CFFBC60.root"
##"root://xrootd.ba.infn.it//store/mc/Phys14DR/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/32ABFE4A-916B-E411-B2FA-00266CFFBC60.root"
#"root://xrootd.ba.infn.it//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/#141029_PU40bx50_PLS170_V6AN2-v1/10000/80161D59-6665-E411-9B4F-C4346BB25698.root",
#"root://xrootd.ba.infn.it:1194//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/10000/8A345C56-6665-E411-9C25-1CC1DE04DF20.root",
#"root://stormgf1.pi.infn.it//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/10000/8A345C56-6665-E411-9C25-1CC1DE04DF20.root",
#"/gpfs/ddn/srm/cms/store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/10000/8A345C56-6665-E411-9C25-1CC1DE04DF20.root",
#"root://xrootd.ba.infn.it//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/10000/9C477248-6665-E411-A9A6-1CC1DE1D0600.root",
#"root://xrootd.ba.infn.it//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/10000/C6D4D875-6665-E411-9E35-00266CF91A18.root"
#"/home/joosep/mac-docs/tth/data/phys14/tth_hbb_phys14_08B36E8F-5E7F-E411-9D5A-002590200AE4.root"
#"/home/joosep/mac-docs/tth/data/phys14/ttjets_phys14_00C90EFC-3074-E411-A845-002590DB9262.root"
#'root://eoscms//eos/cms//store/user/veelken/Phys14/miniAODs/selEvents_TTH_3l1tauh_passingRecTTHtoTauTauDecayModeFilter_original_miniAOD.root'        
],

    #files = ["226BB247-A565-E411-91CF-00266CFF0AF4.root"],
    name="ZHLL125", isEmbed=False,
    splitFactor = 5
    )
sample.isMC=True

#"root://xrootd.ba.infn.it//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/00000/226BB247-A565-E411-91CF-00266CFF0AF4.root"
#/store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/00000/226BB247-A565-E411-91CF-00266CFF0AF4.root
#/store/mc/Spring14miniaod/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/10000/0004A557-C666-E411-8698-549F35AD8B61.root

# the following is declared in case this cfg is used in input to the heppy.py script
selectedComponents = [sample]
from PhysicsTools.HeppyCore.framework.eventsfwlite import Events
config = cfg.Config( components = selectedComponents,
                     sequence = sequence, 
		     services = [output_service],
                     events_class = Events)

class TestFilter(logging.Filter):
    def filter(self, record):
        print record

# and the following runs the process directly 
if __name__ == '__main__':
    from PhysicsTools.HeppyCore.framework.looper import Looper 
    looper = Looper( 'Loop', config, nPrint = 1, nEvents = 1000)

    import time
    import cProfile
    p = cProfile.Profile(time.clock)
    p.runcall(looper.loop)
    p.print_stats()
    looper.write()
