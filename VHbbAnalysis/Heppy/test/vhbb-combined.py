#! /usr/bin/env python

from PhysicsTools.Heppy.utils.cmsswPreprocessor import CmsswPreprocessor

from vhbb import *
# from VHbbAnalysis.Heppy.AdditionalBTag import AdditionalBTag
from VHbbAnalysis.Heppy.AdditionalBoost import AdditionalBoost
from VHbbAnalysis.Heppy.GenHFHadronMatcher import GenHFHadronMatcher


# Add Boosted Information
boostana=cfg.Analyzer(
    verbose=False,
    class_object=AdditionalBoost,
)
sequence.insert(sequence.index(VHbb),boostana)


genhfana=cfg.Analyzer(
    verbose=False,
    class_object=GenHFHadronMatcher,
)
sequence.insert(sequence.index(VHbb),genhfana)


treeProducer.collections["ak08"] = NTupleCollection("FatjetAK08ungroomed",  ak8FatjetType,  10,
                                                    help = "AK, R=0.8, pT > 200 GeV, no grooming")

treeProducer.collections["ak08pruned"] = NTupleCollection("FatjetAK08pruned",
                                                            fourVectorType,
                                                            10,
                                                            help="AK, R=0.8, pT > 200 GeV, pruned zcut=0.1, rcut=0.5, n=2")

treeProducer.collections["ak08prunedsubjets"] = NTupleCollection("SubjetAK08pruned",
                                                                 subjetType,
                                                                 10,
                                                                 help="Subjets of AK, R=0.8, pT > 200 GeV, pruned zcut=0.1, rcut=0.5, n=2")

treeProducer.collections["ak0softdropsubjets"] = NTupleCollection("SubjetAK08softdrop",
                                                                 patSubjetType,
                                                                 10,
                                                                 help="Subjets of AK, R=0.8 softdrop")

if not AdditionalBoost.skip_ca15:
    treeProducer.collections["ca15ungroomed"] = NTupleCollection("FatjetCA15ungroomed",  fatjetType,  10,
                                                                 help = "CA, R=1.5, pT > 200 GeV, no grooming")

    treeProducer.collections["ca15softdrop"] = NTupleCollection("FatjetCA15softdrop",
                                                                fourVectorType,
                                                                10,
                                                                help="CA, R=1.5, pT > 200 GeV, softdrop zcut=0.1, beta=0")

    # four-vector + n-subjettiness
    treeProducer.collections["ca15softdropz2b1"] = NTupleCollection("FatjetCA15softdropz2b1",
                                                                    fatjetTauType,
                                                                    10,
                                                                    help="CA, R=1.5, pT > 200 GeV, softdrop zcut=0.2, beta=1")

    treeProducer.collections["ca15trimmed"] = NTupleCollection("FatjetCA15trimmed",
                                                                fourVectorType,
                                                                10,
                                                                help="CA, R=1.5, pT > 200 GeV, trimmed r=0.2, f=0.06")

    treeProducer.collections["ca15pruned"] = NTupleCollection("FatjetCA15pruned",
                                                                fourVectorType,
                                                                10,
                                                                help="CA, R=1.5, pT > 200 GeV, pruned zcut=0.1, rcut=0.5, n=2")

    treeProducer.collections["ca15prunedsubjets"] = NTupleCollection("SubjetCA15pruned",
                                                                     subjetType,
                                                                     10,
                                                                     help="Subjets of AK, R=1.5, pT > 200 GeV, pruned zcut=0.1, rcut=0.5, n=2")

    treeProducer.collections["httCandidates"] = NTupleCollection("httCandidates",
                                                                 httType,
                                                                 10,
                                                                 help="OptimalR HEPTopTagger Candidates")


# # Add b-Tagging Information
# 
# btagana=cfg.Analyzer(
#     verbose=False,
#     class_object=AdditionalBTag,
# )
# sequence.insert(sequence.index(VHbb),btagana)

# Add Information on generator level hadronic tau decays

from VHbbAnalysis.Heppy.TauGenJetAnalyzer import TauGenJetAnalyzer
TauGenJet = cfg.Analyzer(
    verbose = False,
    class_object = TauGenJetAnalyzer,
)
sequence.insert(sequence.index(VHbb),TauGenJet)

treeProducer.collections["tauGenJets"] = NTupleCollection("GenHadTaus", genTauJetType, 15, help="Generator level hadronic tau decays")
treeProducer.collections["genJetsHadronMatcher"] = NTupleCollection("GenJet",   genJetType, 15, help="Generated jets with hadron matching, sorted by pt descending",filter=lambda x: x.pt() > 20,mcOnly=True)

# Run Everything

preprocessor = CmsswPreprocessor("combined_cmssw.py")
config.preprocessor=preprocessor
if __name__ == '__main__':
    from PhysicsTools.HeppyCore.framework.looper import Looper 
    looper = Looper( 'Loop', config, nPrint = 1, nEvents = 1000)
    import time
    import cProfile
    p = cProfile.Profile(time.clock)
    p.runcall(looper.loop)
    p.print_stats()
    looper.write()
