import sys, os
sys.path.append(os.environ.get("CMSSW_BASE") + "/src/VHbbAnalysis/Heppy/test")

from vhbb_combined import *
components = [
    cfg.MCComponent(
        files = [
            "root://xrootd-cms.infn.it///store/mc/RunIISpring16MiniAODv2/ttHTobb_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/40000/0089CC67-6338-E611-947D-0025904C4E2A.root"
        ],
        name = "tth_hbb",
        isMC = True
    ),
]

def validate_genjet(fn):
    tf = ROOT.TFile.Open(fn)
    tree = tf.Get("vhbb/tree")
    of = ROOT.TFile("out.root", "RECREATE")
    h = ROOT.TH1D("h", "h", 100, 0, 300)
    for ev in tree:
    #    for iGenJet in range(ev.nGenJet):
    #        print "gen", ev.GenJet_pt[iGenJet], iGenJet
        for iJet in range(ev.nJet):
            mcIdx = ev.Jet_mcIdx[iJet]
            if mcIdx >= ev.nGenJet:
                h.Fill(ev.Jet_pt[iJet])
    of.Write()
    of.Close()
    tf.Close()

if __name__ == '__main__':
#    from PhysicsTools.HeppyCore.framework.looper import Looper
#    for comp in components:
#        print "processing",comp.name
#        config.components = [comp]
#        looper = Looper( 'Loop_validation_genjet_' + comp.name, config, nPrint = 0, nEvents = 100) 
#        looper.loop()
#        looper.write()

    validate_genjet("root://storage01.lcg.cscs.ch/pnfs/lcg.cscs.ch/cms/trivcat/store/user/jpata/tth/Jan19_leptonic_nome/ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/Jan19_leptonic_nome/170119_145327/0000/tree_10.root")
