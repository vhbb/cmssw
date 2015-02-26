#! /usr/bin/env python
fastObjects=True

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
		 NTupleVariable("lheV_pt",  lambda ev: ev.lheV_pt, float,mcOnly=True, help="Vector pT at LHE level"),
		 NTupleVariable("lheHT",  lambda ev: ev.lheHT, float,mcOnly=True, help="HT at LHE level"),

 
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
	        "hJets"       : NTupleCollection("hJets",     jetTypeVHbb, 8, sortDescendingBy = lambda jet : jet.btag('combinedSecondaryVertexBJetTags'),help="Higgs jets"),
	        "aJets"       : NTupleCollection("aJets",     jetTypeVHbb, 8, sortDescendingBy = lambda jet : jet.btag('combinedSecondaryVertexBJetTags'),help="Additional jets"),

                "hjidx"       : NTupleCollection("hJidx",    objectInt, 2,help="Higgs jet indices"),
                "hjidxDiJetPtByCSV"       : NTupleCollection("hJidx_sortcsv",    objectInt, 2,help="Higgs jet indices within hJets with CSV sorting "),
                "ajidx"       : NTupleCollection("aJidx",    objectInt, 2,help="additional jet indices"),
                "hjidxCSV"       : NTupleCollection("hJCidx",    objectInt, 2,help="Higgs jet indices CSV"),
                "ajidxCSV"       : NTupleCollection("aJCidx",    objectInt, 2,help="additional jet indices CSV"),
                "hjidx3cj"       : NTupleCollection("hJ3Cidx",    objectInt, 2,help="Higgs jet indices 3 cen jets"),
                "ajidx3cj"       : NTupleCollection("aJ3Cidx",    objectInt, 2,help="additional jet indices 3 cen  jets"),
                "cleanJetsAll"       : NTupleCollection("Jet",     jetTypeVHbb, 15, help="Cental+fwd jets after full selection and cleaning, sorted by b-tag"),
                "selectedTaus"    : NTupleCollection("TauGood", tauType, 3, help="Taus after the preselection"),

		#dump of gen objects
                "genJets"    : NTupleCollection("GenJet",   genParticleType, 15, help="Generated top quarks from hard scattering",filter=lambda x: x.pt() > 20,mcOnly=True),
                "gentopquarks"    : NTupleCollection("GenTop",     genParticleType, 2, help="Generated top quarks from hard scattering"),
                "genbquarksFromH"      : NTupleCollection("GenBQuarkFromH",  genParticleType, 2, help="Generated bottom quarks from Higgs decays"),
                "genbquarksFromHafterISR"      : NTupleCollection("GenBQuarkFromHafterISR",  genParticleType, 2, help="Generated bottom quarks from Higgs decays"),
                "genwzquarks"     : NTupleCollection("GenWZQuark",   genParticleType, 6, help="Generated quarks from W/Z decays"),
                "genleps"         : NTupleCollection("GenLep",     genParticleType, 2, help="Generated leptons from W/Z decays"),
                "gentauleps"      : NTupleCollection("GenLepFromTau", genParticleType, 6, help="Generated leptons from decays of taus from W/Z decays"),
		"genHiggsBoson"   : NTupleCollection("GenHiggsBoson", genParticleType, 1, help="Generated Higgs boson "),
		#"genZbosonsToLL"  : NTupleCollection("GenZbosonsToLL", genParticleType, 6, help="Generated W or Z bosons decaying to LL"),
		#"genWbosonsToLL"  : NTupleCollection("GenWbosonsToLL", genParticleType, 6, help="Generated W or Z bosons decaying to LL"),
		"genvbosons"       : NTupleCollection("GenVbosons", genParticleType, 6, help="Generated W or Z bosons, mass > 30"),
	}
	)

# Lepton Analyzer, take its default config
from PhysicsTools.Heppy.analyzers.objects.LeptonAnalyzer import LeptonAnalyzer
LepAna = LeptonAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.objects.VertexAnalyzer import VertexAnalyzer
VertexAna = VertexAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.objects.PhotonAnalyzer import PhotonAnalyzer
PhoAna = PhotonAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.objects.TauAnalyzer import TauAnalyzer
TauAna = TauAnalyzer.defaultConfig

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


from VHbbAnalysis.Heppy.AdditionalBTag import AdditionalBTag
from VHbbAnalysis.Heppy.VHbbAnalyzer import VHbbAnalyzer
JetAna.jetPt = 15
btagana=cfg.Analyzer(
    verbose=False,
    class_object=AdditionalBTag,
)
VHbb= cfg.Analyzer(
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

    )

#from VHbbAnalysis.Heppy.HeppyShell import HeppyShell
#sh = cfg.Analyzer( class_object=HeppyShell)

from PhysicsTools.Heppy.analyzers.core.TriggerBitAnalyzer import TriggerBitAnalyzer
TrigAna= cfg.Analyzer(
   verbose=False,
   class_object=TriggerBitAnalyzer,
   triggerBits={
   "METBTAG":["HLT_PFMET120_NoiseCleaned_BTagCSV07_v*"],
   "MET":["HLT_PFMET170_NoiseCleaned_v*"],
   "DIELE" : ["HLT_Ele23_Ele12_CaloId_TrackId_Iso_v*"],
   "ELE": ["HLT_Ele32_eta2p1_WP85_Gsf_v*","HLT_Ele32_eta2p1_WP85_Gsf_v*"],
   "DIMU" : ["HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*"],
   "MU" : ["HLT_IsoTkMu24_eta2p1_IterTrk02_v*","HLT_IsoTkMu24_IterTrk02_v*"],
   "TAU": ["HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v*"],
   },
#   processName='HLT',
#   outprefix='HLT'
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

#TrigAna.unrollbits=True

sequence = [LHEAna,FlagsAna, GenAna,VHGenAna,PUAna,TrigAna,VertexAna,LepAna,TauAna,PhoAna,JetAna,METAna,btagana,VHbb,treeProducer]#,sh]


from PhysicsTools.Heppy.utils.miniAodFiles import miniAodFiles
sample = cfg.MCComponent(
    files = [

#       'root://eoscms//eos/cms//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/14053C6D-AD09-E411-A748-00266CFFC9C4.root',
#       'root://eoscms//eos/cms//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/62ED6255-AE09-E411-97CB-00266CFFBF88.root',
#       'root://eoscms//eos/cms//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/72C26B45-AD09-E411-A77C-00266CFFBF80.root',
#       'root://eoscms//eos/cms//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/BAEE7255-AE09-E411-8F9F-00266CFFBF88.root',
#       'root://eoscms//eos/cms//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/D600138D-AD09-E411-917F-00266CFFBF88.root'],
#'root://xrootd.ba.infn.it//store/mc/Phys14DR/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/12328AE8-796B-E411-9D32-002590A831B4.root'
#"32ABFE4A-916B-E411-B2FA-00266CFFBC60.root" #Hbb
#"04860BAA-B673-E411-8B20-002481E0D50C.root" #DY 600
#"ZLL-8A345C56-6665-E411-9C25-1CC1DE04DF20.root" #Hbb Zll
"TTPU20-007B37D4-8B70-E411-BC2D-0025905A6066.root" #Hbb Zll
#"root://xrootd.ba.infn.it//store/mc/Phys14DR/ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/32ABFE4A-916B-E411-B2FA-00266CFFBC60.root"
#"root://xrootd.ba.infn.it//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/#141029_PU40bx50_PLS170_V6AN2-v1/10000/80161D59-6665-E411-9B4F-C4346BB25698.root",
#"root://xrootd.ba.infn.it:1194//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/10000/8A345C56-6665-E411-9C25-1CC1DE04DF20.root",
#"root://stormgf1.pi.infn.it//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/10000/8A345C56-6665-E411-9C25-1CC1DE04DF20.root",
#"/gpfs/ddn/srm/cms/store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/10000/8A345C56-6665-E411-9C25-1CC1DE04DF20.root",
#"root://xrootd.ba.infn.it//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/10000/9C477248-6665-E411-A9A6-1CC1DE1D0600.root",
#"root://xrootd.ba.infn.it//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/10000/C6D4D875-6665-E411-9E35-00266CF91A18.root"
],

    #files = ["226BB247-A565-E411-91CF-00266CFF0AF4.root"],
    name="ZHLL125", isEmbed=False,
    splitFactor = 5
    )
sample.isMC=True

#"root://xrootd.ba.infn.it//store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/00000/226BB247-A565-E411-91CF-00266CFF0AF4.root"
#/store/mc/Spring14miniaod/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/00000/226BB247-A565-E411-91CF-00266CFF0AF4.root
#/store/mc/Spring14miniaod/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/10000/0004A557-C666-E411-8698-549F35AD8B61.root

from PhysicsTools.Heppy.utils.cmsswPreprocessor import CmsswPreprocessor
preprocessor = CmsswPreprocessor("newbtag.py")

# the following is declared in case this cfg is used in input to the heppy.py script
selectedComponents = [sample]
from PhysicsTools.HeppyCore.framework.eventsfwlite import Events
config = cfg.Config( components = selectedComponents,
                     sequence = sequence, 
		     services = [output_service],
                     preprocessor=preprocessor, #this would run cmsRun before running Heppy
                     events_class = Events)

class TestFilter(logging.Filter):
    def filter(self, record):
        print record



# and the following runs the process directly 
if __name__ == '__main__':
    from PhysicsTools.HeppyCore.framework.looper import Looper 
    looper = Looper( 'Loop', config, nPrint = 2, nEvents = 45000)
    import time
    import cProfile
    p = cProfile.Profile(time.clock)
    p.runcall(looper.loop)
    p.print_stats()
#    looper.loop()
    looper.write()
