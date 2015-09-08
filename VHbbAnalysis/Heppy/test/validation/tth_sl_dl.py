import sys, os
sys.path.append(os.environ.get("CMSSW_BASE") + "/src/VHbbAnalysis/Heppy/test")

from vhbb_combined import *
components = [
    cfg.MCComponent(
        files = [
            "/shome/jpata/tth_spring15_miniaod.root"
        ],
        name = "tth_hbb",
        isMC = True
    ),
    cfg.MCComponent(
        files = [
            "root://xrootd.unl.edu//store/mc/RunIISpring15DR74/TT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/009A8083-EC2E-E511-ABB1-3417EBE64B91.root"
        ],
        name = "ttjets",
        isMC = True
    )
]

if __name__ == '__main__':
    from PhysicsTools.HeppyCore.framework.looper import Looper
    for comp in components:
        looper = Looper( 'Loop_validation_tth_sl_dl_' + comp.name, config, nPrint = 0, nEvents = 1000)
        config.components = [comp] 
        looper.loop()
        looper.write()
