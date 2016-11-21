import ROOT, sys
from PhysicsTools.Heppy.physicsutils.JetReCalibrator import JetReCalibrator

files = ["root://storage01.lcg.cscs.ch/pnfs/lcg.cscs.ch/cms/trivcat/store/user/jpata/tth/Nov18_test_v1/ttHTobb_M125_13TeV_powheg_pythia8/Nov18_test_v1/161118_143453/0000/tree_{0}.root".format(ifile) for ifile in range(1, 21)]

of = ROOT.TFile("out.root", "RECREATE")
tt = ROOT.TChain("vhbb/tree")
for fi in files:
    tt.AddFile(fi)

of.cd()
tt.Draw("Jet_pt >>+ nominal(100,0,500)")
tt.Draw("Jet_pt*Jet_corr >>+ corr(100,0,500)")
for jet_corr in JetReCalibrator.factorizedJetCorrections + ["JEC"]:
    s1 = "(Jet_pt * Jet_corr_{0}Up) >>+ h__{0}__Up(100,0,500)".format(jet_corr)
    s2 = "(Jet_pt * Jet_corr_{0}Down) >>+ h__{0}__Down(100,0,500)".format(jet_corr)
    print s1
    print s2
    tt.Draw(s1)
    tt.Draw(s2)
of.Write()
of.Close()
