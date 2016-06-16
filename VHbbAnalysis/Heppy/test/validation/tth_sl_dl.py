import sys, os
sys.path.append(os.environ.get("CMSSW_BASE") + "/src/VHbbAnalysis/Heppy/test")

from vhbb_combined import *
components = [
    #cfg.MCComponent(
    #    files = [
    #        "root://xrootd-cms.infn.it///store/mc/RunIISpring15DR74/ttHTobb_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/141B9915-1F08-E511-B9FF-001E675A6AB3.root"
    #    ],
    #    name = "tth_hbb",
    #    isMC = True
    #),
    cfg.MCComponent(
        files = [
            "root://xrootd-cms.infn.it///store/mc/RunIISpring16MiniAODv2/TT_TuneEE5C_13TeV-powheg-herwigpp/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/C2EC599F-7C1B-E611-BCA9-B083FECF837B.root"
        ],
        name = "ttjets",
        isMC = True
    )
]

if __name__ == '__main__':
    from PhysicsTools.HeppyCore.framework.looper import Looper
    for comp in components:
        print "processing",comp
        config.components = [comp] 
        looper = Looper( 'Loop_validation_tth_sl_dl_' + comp.name, config, nPrint = 0, nEvents = 100)
        looper.loop()
        looper.write()
