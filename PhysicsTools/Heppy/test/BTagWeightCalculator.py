#add per-jet b-tag systematic weight
import os
import numpy as np
from PhysicsTools.Heppy.physicsutils.BTagWeightCalculator import BTagWeightCalculator
csvpath = os.environ['CMSSW_BASE']+"/src/PhysicsTools/Heppy/data"
bweightcalc = BTagWeightCalculator(
    csvpath + "/csv_rwt_fit_hf_2015_11_20.root",
    csvpath + "/csv_rwt_fit_lf_2015_11_20.root"
)
bweightcalc.btag = "btag"

class Jet:
    pass
for btag in np.linspace(0,1,11):
    jet = Jet()
    jet.pt = 40
    jet.eta = 0.3
    jet.mcFlavour = 5
    jet.btag = btag
    
    wup = bweightcalc.calcJetWeight(jet, "final", "JESUp")
    wdown = bweightcalc.calcJetWeight(jet, "final", "JESDown")
    w = bweightcalc.calcJetWeight(jet, "final", "nominal")
    print jet.__dict__, w, wup, wdown
