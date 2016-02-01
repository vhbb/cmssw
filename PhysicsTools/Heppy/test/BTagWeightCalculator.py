#add per-jet b-tag systematic weight
import os
import numpy as np
from PhysicsTools.Heppy.physicsutils.BTagWeightCalculator import BTagWeightCalculator
import ROOT

csvpath = os.environ['CMSSW_BASE']+"/src/PhysicsTools/Heppy/data"
bweightcalc = BTagWeightCalculator(
    csvpath + "/csv_rwt_fit_hf_2015_12_14.root",
    csvpath + "/csv_rwt_fit_lf_2015_12_14.root"
)
bweightcalc.btag = "btag"

class Jet:
    pass

inf = open("jets.csv")

systs = ["nominal", "JESUp", "JESDown", "HFUp", "HFDown", "LFUp", "LFDown", "Stats1Up", "Stats1Down", "Stats2Up", "Stats2Down",]

hs = {
    s: ROOT.TH1D("weight_"+s, "weight", 100, 0, 2) for s in systs
}

ws = {
    s: 1.0 for s in systs
}
iev_prev = 0
for line in inf.readlines():
    iev, ij, pt, eta, csv, fl = map(float, line.split(","))
    iev = int(iev)
    ij = int(ij)
    fl = int(fl)
    if iev != iev_prev:
        print "tot", ws["nominal"]
        for s in systs:
            hs[s].Fill(ws[s])
            ws[s] = 1.0
    jet = Jet()
    jet.pt = pt
    jet.eta = eta
    jet.hadronFlavour = fl
    jet.btag = csv
    
    w = {
        s: bweightcalc.calcJetWeight(jet, "final", s) for s in systs
    }
    print line.strip(), w["nominal"]
    

    for s in systs:
        ws[s] *= w[s]
    
    iev_prev = iev

for s in systs:
    print s, hs[s].GetMean(), hs[s].GetRMS()
