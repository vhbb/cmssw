import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)

debug = False

def hist1D(name, title, var, nBins, xLow, xHigh):
    h = ROOT.TH1D(name, title, nBins, xLow, xHigh)
    h.GetXaxis().SetTitle( var )
    return h

def hist2D(name, title, varx, xBins, xLow, xHigh, vary, yBins, yLow, yHigh):
    h = ROOT.TH2D(name, title, xBins, xLow, xHigh, yBins, yLow, yHigh)
    h.GetXaxis().SetTitle( varx )
    h.GetYaxis().SetTitle( vary )
    return h

def findbin_jetEta(eta):
    bin = "Bin0"
    if abs(eta)<1.5:
        bin = "Bin0"
    elif abs(eta)>1.5 and abs(eta)<2.5:
        bin = "Bin1"
    elif abs(eta)>2.5 and abs(eta)<4.5:
        bin = "Bin2"
    return bin

def findbin_jetPt(pt):
    bin = "Bin0"
    if pt>20 and pt<30:
        bin = "Bin0"
    elif pt>30 and pt<50:
        bin = "Bin1"
    elif pt>50 and pt<100:
        bin = "Bin2"
    else:
        bin = "Bin3"
    return bin

def findbin_genMetPt(met_genPt):
    genPtBin = ""
    if  met_genPt<50.:
         genPtBin = "_0to50"
    elif met_genPt<100.:
         genPtBin = "_50to100"
    else:
        genPtBin = "_100toinf"
    return genPtBin

def weight(ev):
    return ev.puWeight*abs(ev.genWeight)

def category(count_genjet_b, count_genjet_2b, count_genjet_c, count_genjet_2c):
    gen_cat_count = -1
    if (count_genjet_b + count_genjet_2b)>=2:
        gen_cat_count = 0
    elif count_genjet_b == 1:
        gen_cat_count = 1
    elif count_genjet_2b == 1:
        gen_cat_count = 2
    elif (count_genjet_c + count_genjet_2c)>=2:
        gen_cat_count = 3
    elif count_genjet_c ==1:
        gen_cat_count = 4        
    elif count_genjet_2c ==1:
        gen_cat_count = 5        
    return gen_cat_count


def pass_event(ev):
    return (ev.Vtype==0 or ev.Vtype==1 or ev.Vtype==2 or ev.Vtype==3)

##############################################################################################################

#path   = "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat///store/t3groups/ethz-higgs/run2/V12/"
path   = "dcap://t3se01.psi.ch:22125///pnfs/psi.ch/cms/trivcat///store/t3groups/ethz-higgs/run2/VHBBHeppyV13/"

#sample = "ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/VHBB_HEPPY_V13_ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/151002_083741/0000/" 
sample = "TT_TuneCUETP8M1_13TeV-powheg-pythia8/VHBB_HEPPY_V13_TT_TuneCUETP8M1_13TeV-powheg-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext3-v1/151002_060637/0000/"

file_names = [
    "tree_1.root",
    ]

for p in range(2, 10):
    file_names.append( "tree_"+str(p)+".root" )

chain = ROOT.TChain("tree")
for file_name in file_names:
    chain.AddFile( path + "/" + sample + "/" + file_name )
chain.SetBranchStatus("*",        False)
chain.SetBranchStatus("nJet",      True)
chain.SetBranchStatus("Jet_*",     True)
chain.SetBranchStatus("nGenJet",   True)
chain.SetBranchStatus("GenJet_*",  True)
chain.SetBranchStatus("hJidx",     True)
chain.SetBranchStatus("met*",      True)
chain.SetBranchStatus("tkMet*",    True)
chain.SetBranchStatus("Vtype",     True)
chain.SetBranchStatus("V_*",       True)
chain.SetBranchStatus("H_*",       True)
chain.SetBranchStatus("puWeight",  True)
chain.SetBranchStatus("genWeight", True)

out      = ROOT.TFile("JetMetValidation.root", "RECREATE")
hists = {

    "MET" : {
        
        "Pt_systematics" : {
            "JetEnUp_0to50"    : hist1D("MET_Pt_JetEnUp_type1_0",   "MET p_{T} resolution JetEnUp, 0<gen E_{T}^{miss}<50 GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            "JetEnUp_50to100"  : hist1D("MET_Pt_JetEnUp_type1_1",   "MET p_{T} resolution JetEnUp, 50<gen E_{T}^{miss}<100 GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            "JetEnUp_100toinf" : hist1D("MET_Pt_JetEnUp_type1_2",   "MET p_{T} resolution JetEnUp, 100<gen E_{T}^{miss} GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            "JetEnDown_0to50"    : hist1D("MET_Pt_JetEnDown_type1_0",   "MET p_{T} resolution JetEnDown, 0<gen E_{T}^{miss}<50 GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            "JetEnDown_50to100"  : hist1D("MET_Pt_JetEnDown_type1_1",   "MET p_{T} resolution JetEnDown, 50<gen E_{T}^{miss}<100 GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            "JetEnDown_100toinf" : hist1D("MET_Pt_JetEnDown_type1_2",   "MET p_{T} resolution JetEnDown, 100<gen E_{T}^{miss} GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,

            "JetResUp_0to50"    : hist1D("MET_Pt_JetResUp_type1_0",   "MET p_{T} resolution JetResUp, 0<gen E_{T}^{miss}<50 GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            "JetResUp_50to100"  : hist1D("MET_Pt_JetResUp_type1_1",   "MET p_{T} resolution JetResUp, 50<gen E_{T}^{miss}<100 GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            "JetResUp_100toinf" : hist1D("MET_Pt_JetResUp_type1_2",   "MET p_{T} resolution JetResUp, 100<gen E_{T}^{miss} GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            "JetResDown_0to50"    : hist1D("MET_Pt_JetResDown_type1_0",   "MET p_{T} resolution JetResDown, 0<gen E_{T}^{miss}<50 GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            "JetResDown_50to100"  : hist1D("MET_Pt_JetResDown_type1_1",   "MET p_{T} resolution JetResDown, 50<gen E_{T}^{miss}<100 GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            "JetResDown_100toinf" : hist1D("MET_Pt_JetResDown_type1_2",   "MET p_{T} resolution JetResDown, 100<gen E_{T}^{miss} GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,

            "UnclusteredEnUp_0to50"    : hist1D("MET_Pt_UnclusteredEnUp_type1_0",   "MET p_{T} resolution UnclusteredEnUp, 0<gen E_{T}^{miss}<50 GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            "UnclusteredEnUp_50to100"  : hist1D("MET_Pt_UnclusteredEnUp_type1_1",   "MET p_{T} resolution UnclusteredEnUp, 50<gen E_{T}^{miss}<100 GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            "UnclusteredEnUp_100toinf" : hist1D("MET_Pt_UnclusteredEnUp_type1_2",   "MET p_{T} resolution UnclusteredEnUp, 100<gen E_{T}^{miss} GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            "UnclusteredEnDown_0to50"    : hist1D("MET_Pt_UnclusteredEnDown_type1_0",   "MET p_{T} resolution UnclusteredEnDown, 0<gen E_{T}^{miss}<50 GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            "UnclusteredEnDown_50to100"  : hist1D("MET_Pt_UnclusteredEnDown_type1_1",   "MET p_{T} resolution UnclusteredEnDown, 50<gen E_{T}^{miss}<100 GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            "UnclusteredEnDown_100toinf" : hist1D("MET_Pt_UnclusteredEnDown_type1_2",   "MET p_{T} resolution UnclusteredEnDown, 100<gen E_{T}^{miss} GeV",  "Variation/Nominal", 40, 0.5, 1.5) ,
            },
        
        "Pt_response" : {
            "raw_0to50"     :  hist1D("MET_Pt_Resolution_raw_0",   "Raw MET p_{T} resolution, 0<gen E_{T}^{miss}<50 GeV",     "Reco/Gen", 40, 0, 2) ,
            "raw_50to100"   :  hist1D("MET_Pt_Resolution_raw_1",   "Raw MET p_{T} resolution, 50<gen E_{T}^{miss}<100 GeV",   "Reco/Gen", 40, 0, 2) ,
            "raw_100toinf"  :  hist1D("MET_Pt_Resolution_raw_2",   "Raw MET p_{T} resolution, 100<gen E_{T}^{miss} GeV",      "Reco/Gen", 40, 0, 2) ,
            "type1_0to50"   :  hist1D("MET_Pt_Resolution_type1_0", "Type1 MET p_{T} resolution, 0<gen E_{T}^{miss}<50 GeV",   "Reco/Gen", 40, 0, 2) ,
            "type1_50to100" :  hist1D("MET_Pt_Resolution_type1_1", "Type1 MET p_{T} resolution, 50<gen E_{T}^{miss}<100 GeV", "Reco/Gen", 40, 0, 2) ,
            "type1_100toinf":  hist1D("MET_Pt_Resolution_type1_2", "Type1 MET p_{T} resolution, 100<gen E_{T}^{miss} GeV",    "Reco/Gen", 40, 0, 2) ,
            "type1p2_0to50"   :  hist1D("MET_Pt_Resolution_type1p2_0", "Type1p2 MET p_{T} resolution, 0<gen E_{T}^{miss}<50 GeV",   "Reco/Gen", 40, 0, 2) ,
            "type1p2_50to100" :  hist1D("MET_Pt_Resolution_type1p2_1", "Type1p2 MET p_{T} resolution, 50<gen E_{T}^{miss}<100 GeV", "Reco/Gen", 40, 0, 2) ,
            "type1p2_100toinf":  hist1D("MET_Pt_Resolution_type1p2_2", "Type1p2 MET p_{T} resolution, 100<gen E_{T}^{miss} GeV",    "Reco/Gen", 40, 0, 2) ,
            "puppi_0to50"   :  hist1D("MET_Pt_Resolution_puppi_0", "Puppi MET p_{T} resolution, 0<gen E_{T}^{miss}<50 GeV",   "Reco/Gen", 40, 0, 2) ,
            "puppi_50to100" :  hist1D("MET_Pt_Resolution_puppi_1", "Puppi MET p_{T} resolution, 50<gen E_{T}^{miss}<100 GeV", "Reco/Gen", 40, 0, 2) ,
            "puppi_100toinf":  hist1D("MET_Pt_Resolution_puppi_2", "Puppi MET p_{T} resolution, 100<gen E_{T}^{miss} GeV",    "Reco/Gen", 40, 0, 2) ,
            "tkMet_0to50"   :  hist1D("MET_Pt_Resolution_tkMet_0", "TkMet MET p_{T} resolution, 0<gen E_{T}^{miss}<50 GeV",   "Reco/Gen", 40, 0, 2) ,
            "tkMet_50to100" :  hist1D("MET_Pt_Resolution_tkMet_1", "TkMet MET p_{T} resolution, 50<gen E_{T}^{miss}<100 GeV", "Reco/Gen", 40, 0, 2) ,
            "tkMet_100toinf":  hist1D("MET_Pt_Resolution_tkMet_2", "TkMet MET p_{T} resolution, 100<gen E_{T}^{miss} GeV",    "Reco/Gen", 40, 0, 2) ,
            },
        "PxPy_correlation" : {
            "type1_0to50"    :  hist2D("MET_Pxy_Correlation_type1_0", "Type1 MET #Deltap_{x.y} correlation,  0<gen E_{T}^{miss}<50  GeV", "(Reco-Gen)_{x}", 40, -100, 100, "(Reco-Gen)_{y}", 40, -100, 100) ,
            "type1_50to100"  :  hist2D("MET_Pxy_Correlation_type1_1", "Type1 MET #Deltap_{x.y} correlation, 50<gen E_{T}^{miss}<100 GeV", "(Reco-Gen)_{x}", 40, -100, 100, "(Reco-Gen)_{y}", 40, -100, 100) ,
            "type1_100toinf" :  hist2D("MET_Pxy_Correlation_type1_2", "Type1 MET #Deltap_{x.y} correlation, 100<gen E_{T}^{miss} GeV",    "(Reco-Gen)_{x}", 40, -100, 100, "(Reco-Gen)_{y}", 40, -100, 100) ,
            },
        "Phi_resolution_2D" : {
            "type1" :  hist2D("MET_Phi_Resolution_type1", "Type1 MET #phi resolution", "Gen p_{T}", 40, 0, 200, "(Reco-Gen)", 40, -2., 2.) ,
            },
        "Px_resolution_vsSumET_2D" : {
            "type1" :  hist2D("MET_Px_Resolution_type1_vsSumET", "Type1 MET p_{x} resolution vs SumET", "SumET", 20, 400, 2000, "(Reco-Gen)", 40, -100, 100) ,
            }

        },

    "Classification" : {
        "Gen" : {
            "ByGenJets"   : hist1D("GenJets",   "Gen-level classification based on jet counting",       "category", 7, -1, 6 ),
            "ByHiggsJets" : hist1D("HiggsJets", "Gen-level classification based on Higgs jet counting", "category", 7, -1, 6 ),
            "Compare_2D"  : hist2D("CompareClass",  "Gen-level classification", "ByGenJets" , 7, -1, 6, "ByHiggsJets",  7, -1, 6 ),
            "Compare_HVPt100_2D"  : hist2D("CompareClass_HVPt100",  "Gen-level classification, V,H p_{T}>100 GeV", "ByGenJets" , 7, -1, 6, "ByHiggsJets",  7, -1, 6 ),
            "Compare_H_btag_2D"   : hist2D("CompareClass_H_btag",  "Gen-level classification, H jets passing CSVM", "ByGenJets" , 7, -1, 6, "ByHiggsJets",  7, -1, 6 ),
            },
        },

    "Jet" : {

        "PuId_match_2D" : {
            "Bin0" : hist2D("Jet_puId_match_Bin0", "Gen-matched Jet pile up ID efficiency, 0.0<|#eta|<1.5", "p_{T} [GeV]", 20,  20, 400, "score", 2, -0.5, 1.5) ,
            "Bin1" : hist2D("Jet_puId_match_Bin1", "Gen-matched Jet pile up ID efficiency, 1.5<|#eta|<2.5", "p_{T} [GeV]", 20,  20, 400, "score", 2, -0.5, 1.5) ,
            "Bin2" : hist2D("Jet_puId_match_Bin2", "Gen-matched Jet pile up ID efficiency, 2.5<|#eta|<4.5", "p_{T} [GeV]", 20,  20, 400, "score", 2, -0.5, 1.5) ,
            },

        "PuId_unmatch_2D" : {
            "Bin0" : hist2D("Jet_puId_unmatch_Bin0", "Gen-unmatched Jet pile up ID efficiency, 0.0<|#eta|<1.5", "p_{T} [GeV]", 20,  20, 400, "score", 2, -0.5, 1.5) ,
            "Bin1" : hist2D("Jet_puId_unmatch_Bin1", "Gen-unmatched Jet pile up ID efficiency, 1.5<|#eta|<2.5", "p_{T} [GeV]", 20,  20, 400, "score", 2, -0.5, 1.5) ,
            "Bin2" : hist2D("Jet_puId_unmatch_Bin2", "Gen-unmatched Jet pile up ID efficiency, 2.5<|#eta|<4.5", "p_{T} [GeV]", 20,  20, 400, "score", 2, -0.5, 1.5) ,
            },

        "JES_2D" : {
            "Bin0" : hist2D("Jet_Pt_CorrJES_Bin0", "Jet p_{T} correction, 0.0<|#eta|<1.5", "raw p_{T} [GeV]", 20,  20, 400, "JES", 40, 0.6, 1.2) ,
            "Bin1" : hist2D("Jet_Pt_CorrJES_Bin1", "Jet p_{T} correction, 1.5<|#eta|<2.5", "raw p_{T} [GeV]", 20,  20, 400, "JES", 40, 0.6, 1.2) ,
            "Bin2" : hist2D("Jet_Pt_CorrJES_Bin2", "Jet p_{T} correction, 2.5<|#eta|<4.5", "raw p_{T} [GeV]", 20,  20, 400, "JES", 40, 0.6, 1.2) ,
            },

        "JESUp_2D" : {
            "Bin0" : hist2D("Jet_Pt_CorrJESUp_Bin0", "Jet p_{T} correction uncert. Up, 0.0<|#eta|<1.5", "raw p_{T} [GeV]", 20,  20, 400, "JES Up", 40, 0.9, 1.1) ,
            "Bin1" : hist2D("Jet_Pt_CorrJESUp_Bin1", "Jet p_{T} correction uncert. Up, 1.5<|#eta|<2.5", "raw p_{T} [GeV]", 20,  20, 400, "JES Up", 40, 0.9, 1.1) ,
            "Bin2" : hist2D("Jet_Pt_CorrJESUp_Bin2", "Jet p_{T} correction uncert. Up, 2.5<|#eta|<4.5", "raw p_{T} [GeV]", 20,  20, 400, "JES Up", 40, 0.9, 1.1) ,
            },

        "JESDown_2D" : {
            "Bin0" : hist2D("Jet_Pt_CorrJESDown_Bin0", "Jet p_{T} correction uncert. Down, 0.0<|#eta|<1.5", "raw p_{T} [GeV]", 20,  20, 400, "JES Down", 40, 0.9, 1.1) ,
            "Bin1" : hist2D("Jet_Pt_CorrJESDown_Bin1", "Jet p_{T} correction uncert. Down, 1.5<|#eta|<2.5", "raw p_{T} [GeV]", 20,  20, 400, "JES Down", 40, 0.9, 1.1) ,
            "Bin2" : hist2D("Jet_Pt_CorrJESDown_Bin2", "Jet p_{T} correction uncert. Down, 2.5<|#eta|<4.5", "raw p_{T} [GeV]", 20,  20, 400, "JES Down", 40, 0.9, 1.1) ,
            },

        "JER_2D" : {
            "Bin0" : hist2D("Jet_Pt_CorrJER_Bin0", "Jet p_{T} correction, 0.0<|#eta|<1.5", "p_{T} [GeV]", 20,  20, 400, "JER", 40, 0.6, 1.2) ,
            "Bin1" : hist2D("Jet_Pt_CorrJER_Bin1", "Jet p_{T} correction, 1.5<|#eta|<2.5", "p_{T} [GeV]", 20,  20, 400, "JER", 40, 0.6, 1.2) ,
            "Bin2" : hist2D("Jet_Pt_CorrJER_Bin2", "Jet p_{T} correction, 2.5<|#eta|<4.5", "p_{T} [GeV]", 20,  20, 400, "JER", 40, 0.6, 1.2) ,
            },

        "JERUp_2D" : {
            "Bin0" : hist2D("Jet_Pt_CorrJERUp_Bin0", "Jet p_{T} correction uncert. Up, 0.0<|#eta|<1.5", "p_{T} [GeV]", 20,  20, 400, "JER Up", 40, 0.9, 1.1) ,
            "Bin1" : hist2D("Jet_Pt_CorrJERUp_Bin1", "Jet p_{T} correction uncert. Up, 1.5<|#eta|<2.5", "p_{T} [GeV]", 20,  20, 400, "JER Up", 40, 0.9, 1.1) ,
            "Bin2" : hist2D("Jet_Pt_CorrJERUp_Bin2", "Jet p_{T} correction uncert. Up, 2.5<|#eta|<4.5", "p_{T} [GeV]", 20,  20, 400, "JER Up", 40, 0.9, 1.1) ,
            },

        "JERDown_2D" : {
            "Bin0" : hist2D("Jet_Pt_CorrJERDown_Bin0", "Jet p_{T} correction uncert. Down, 0.0<|#eta|<1.5", "p_{T} [GeV]", 20,  20, 400, "JER Down", 40, 0.9, 1.1) ,
            "Bin1" : hist2D("Jet_Pt_CorrJERDown_Bin1", "Jet p_{T} correction uncert. Down, 1.5<|#eta|<2.5", "p_{T} [GeV]", 20,  20, 400, "JER Down", 40, 0.9, 1.1) ,
            "Bin2" : hist2D("Jet_Pt_CorrJERDown_Bin2", "Jet p_{T} correction uncert. Down, 2.5<|#eta|<4.5", "p_{T} [GeV]", 20,  20, 400, "JER Down", 40, 0.9, 1.1) ,
            },

        "Pt_resolution_wElectron_2D"   : {
            "Bin0" : hist2D("Jet_Pt_Resolution_wElectron_Bin0", "Jet with leptons p_{T} response, 0.0<|#eta|<1.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin1" : hist2D("Jet_Pt_Resolution_wElectron_Bin1", "Jet with leptons p_{T} response, 1.5<|#eta|<2.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin2" : hist2D("Jet_Pt_Resolution_wElectron_Bin2", "Jet with leptons p_{T} response, 2.5<|#eta|<4.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            },

        "Pt_resolution_wMuon_2D"   : {
            "Bin0" : hist2D("Jet_Pt_Resolution_wMuon_Bin0", "Jet with leptons p_{T} response, 0.0<|#eta|<1.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin1" : hist2D("Jet_Pt_Resolution_wMuon_Bin1", "Jet with leptons p_{T} response, 1.5<|#eta|<2.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin2" : hist2D("Jet_Pt_Resolution_wMuon_Bin2", "Jet with leptons p_{T} response, 2.5<|#eta|<4.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            },

        "Pt_resolution_2D"   : {
            "Bin0" : hist2D("Jet_Pt_Resolution_Bin0", "Jet p_{T} response, 0.0<|#eta|<1.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin1" : hist2D("Jet_Pt_Resolution_Bin1", "Jet p_{T} response, 1.5<|#eta|<2.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin2" : hist2D("Jet_Pt_Resolution_Bin2", "Jet p_{T} response, 2.5<|#eta|<4.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            },

        "Eta_resolution_2D"  : {
            "Bin0" : hist2D("Jet_Eta_Resolution_Bin0", "Jet #eta response, 0.0<|#eta|<1.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)",     40, -0.2, 0.2) ,
            "Bin1" : hist2D("Jet_Eta_Resolution_Bin1", "Jet #eta response, 1.5<|#eta|<2.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)",     40, -0.2, 0.2) ,
            "Bin2" : hist2D("Jet_Eta_Resolution_Bin2", "Jet #eta response, 2.5<|#eta|<4.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)",     40, -0.2, 0.2) ,
            }, 

        "Phi_resolution_2D"  : {
            "Bin0" : hist2D("Jet_Phi_Resolution_Bin0", "Jet #phi response, 0.0<|#eta|<1.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)",     40, -0.2, 0.2) ,
            "Bin1" : hist2D("Jet_Phi_Resolution_Bin1", "Jet #phi response, 1.5<|#eta|<2.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)",     40, -0.2, 0.2) ,
            "Bin2" : hist2D("Jet_Phi_Resolution_Bin2", "Jet #phi response, 2.5<|#eta|<4.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)",     40, -0.2, 0.2) ,
            }, 


        },

    "RegressedJet" : {
        "Pt_resolution_2D"   : {
            "Bin0" : hist2D("RegressedJet_Pt_Resolution_Bin0", "Regressed Jet p_{T} response, 0.0<|#eta|<1.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin1" : hist2D("RegressedJet_Pt_Resolution_Bin1", "Regressed Jet p_{T} response, 1.5<|#eta|<2.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin2" : hist2D("RegressedJet_Pt_Resolution_Bin2", "Regressed Jet p_{T} response, 2.5<|#eta|<4.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            },       
        "Pt_resolution_wNu_2D"   : {
            "Bin0" : hist2D("RegressedJet_Pt_Resolution_wNu_Bin0", "Regressed Jet p_{T} response with #nu, 0.0<|#eta|<1.5", "p_{T} gen + #nu [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin1" : hist2D("RegressedJet_Pt_Resolution_wNu_Bin1", "Regressed Jet p_{T} response with #nu, 1.5<|#eta|<2.5", "p_{T} gen + #nu [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin2" : hist2D("RegressedJet_Pt_Resolution_wNu_Bin2", "Regressed Jet p_{T} response with #nu, 2.5<|#eta|<4.5", "p_{T} gen + #nu [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            },        
 
        },

    "BJet" : {
        "Pt_resolution_2D"   : {
            "Bin0" : hist2D("BJet_Pt_Resolution_Bin0", "Jet (b) p_{T} response, 0.0<|#eta|<1.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin1" : hist2D("BJet_Pt_Resolution_Bin1", "Jet (b) p_{T} response, 1.5<|#eta|<2.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin2" : hist2D("BJet_Pt_Resolution_Bin2", "Jet (b) p_{T} response, 2.5<|#eta|<4.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            },
        "NumBHadrons" : {
            "Bin0" : hist1D("BJet_numBHadrons_Bin0", "Number of b hadrons in b-jets, 0.0<|#eta|<1.5", "number of b hadrons", 5,-1,4),
            "Bin1" : hist1D("BJet_numBHadrons_Bin1", "Number of b hadrons in b-jets, 1.5<|#eta|<2.5", "number of b hadrons", 5,-1,4),
            "Bin2" : hist1D("BJet_numBHadrons_Bin2", "Number of b hadrons in b-jets, 2.5<|#eta|<4.5", "number of b hadrons", 5,-1,4),
            },

        "Btag_weight_2D" : {
            "Bin0_Bin0" : hist2D("BJet_BTagWeight_Bin0_Bin0", "BJet btag weight, 0.0<|#eta|<1.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            "Bin0_Bin1" : hist2D("BJet_BTagWeight_Bin0_Bin1", "BJet btag weight, 0.0<|#eta|<1.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            "Bin0_Bin2" : hist2D("BJet_BTagWeight_Bin0_Bin2", "BJet btag weight, 0.0<|#eta|<1.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            "Bin0_Bin3" : hist2D("BJet_BTagWeight_Bin0_Bin3", "BJet btag weight, 0.0<|#eta|<1.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            "Bin1_Bin0" : hist2D("BJet_BTagWeight_Bin1_Bin0", "BJet btag weight, 1.5<|#eta|<2.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            "Bin1_Bin1" : hist2D("BJet_BTagWeight_Bin1_Bin1", "BJet btag weight, 1.5<|#eta|<2.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            "Bin1_Bin2" : hist2D("BJet_BTagWeight_Bin1_Bin2", "BJet btag weight, 1.5<|#eta|<2.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            "Bin1_Bin3" : hist2D("BJet_BTagWeight_Bin1_Bin3", "BJet btag weight, 1.5<|#eta|<2.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            },

        "Btag_weightJES_2D" : {
            "Bin0_Bin0" : hist2D("BJet_BTagWeightJES_Bin0_Bin0", "BJet btag weight JES Up/Down, 0.0<|#eta|<1.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin1" : hist2D("BJet_BTagWeightJES_Bin0_Bin1", "BJet btag weight JES Up/Down, 0.0<|#eta|<1.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin2" : hist2D("BJet_BTagWeightJES_Bin0_Bin2", "BJet btag weight JES Up/Down, 0.0<|#eta|<1.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin3" : hist2D("BJet_BTagWeightJES_Bin0_Bin3", "BJet btag weight JES Up/Down, 0.0<|#eta|<1.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin0" : hist2D("BJet_BTagWeightJES_Bin1_Bin0", "BJet btag weight JES Up/Down, 1.5<|#eta|<2.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin1" : hist2D("BJet_BTagWeightJES_Bin1_Bin1", "BJet btag weight JES Up/Down, 1.5<|#eta|<2.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin2" : hist2D("BJet_BTagWeightJES_Bin1_Bin2", "BJet btag weight JES Up/Down, 1.5<|#eta|<2.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin3" : hist2D("BJet_BTagWeightJES_Bin1_Bin3", "BJet btag weight JES Up/Down, 1.5<|#eta|<2.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            },

        "Btag_weightLF_2D" : {
            "Bin0_Bin0" : hist2D("BJet_BTagWeightLF_Bin0_Bin0", "BJet btag weight LF Up/Down, 0.0<|#eta|<1.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin1" : hist2D("BJet_BTagWeightLF_Bin0_Bin1", "BJet btag weight LF Up/Down, 0.0<|#eta|<1.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin2" : hist2D("BJet_BTagWeightLF_Bin0_Bin2", "BJet btag weight LF Up/Down, 0.0<|#eta|<1.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin3" : hist2D("BJet_BTagWeightLF_Bin0_Bin3", "BJet btag weight LF Up/Down, 0.0<|#eta|<1.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin0" : hist2D("BJet_BTagWeightLF_Bin1_Bin0", "BJet btag weight LF Up/Down, 1.5<|#eta|<2.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin1" : hist2D("BJet_BTagWeightLF_Bin1_Bin1", "BJet btag weight LF Up/Down, 1.5<|#eta|<2.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin2" : hist2D("BJet_BTagWeightLF_Bin1_Bin2", "BJet btag weight LF Up/Down, 1.5<|#eta|<2.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin3" : hist2D("BJet_BTagWeightLF_Bin1_Bin3", "BJet btag weight LF Up/Down, 1.5<|#eta|<2.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            },

        "Btag_weightStats1_2D" : {
            "Bin0_Bin0" : hist2D("BJet_BTagWeightStats1_Bin0_Bin0", "BJet btag weight Stats1 Up/Down, 0.0<|#eta|<1.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin1" : hist2D("BJet_BTagWeightStats1_Bin0_Bin1", "BJet btag weight Stats1 Up/Down, 0.0<|#eta|<1.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin2" : hist2D("BJet_BTagWeightStats1_Bin0_Bin2", "BJet btag weight Stats1 Up/Down, 0.0<|#eta|<1.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin3" : hist2D("BJet_BTagWeightStats1_Bin0_Bin3", "BJet btag weight Stats1 Up/Down, 0.0<|#eta|<1.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin0" : hist2D("BJet_BTagWeightStats1_Bin1_Bin0", "BJet btag weight Stats1 Up/Down, 1.5<|#eta|<2.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin1" : hist2D("BJet_BTagWeightStats1_Bin1_Bin1", "BJet btag weight Stats1 Up/Down, 1.5<|#eta|<2.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin2" : hist2D("BJet_BTagWeightStats1_Bin1_Bin2", "BJet btag weight Stats1 Up/Down, 1.5<|#eta|<2.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin3" : hist2D("BJet_BTagWeightStats1_Bin1_Bin3", "BJet btag weight Stats1 Up/Down, 1.5<|#eta|<2.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            },

        "Btag_weightStats2_2D" : {
            "Bin0_Bin0" : hist2D("BJet_BTagWeightStats2_Bin0_Bin0", "BJet btag weight Stats2 Up/Down, 0.0<|#eta|<1.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin1" : hist2D("BJet_BTagWeightStats2_Bin0_Bin1", "BJet btag weight Stats2 Up/Down, 0.0<|#eta|<1.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin2" : hist2D("BJet_BTagWeightStats2_Bin0_Bin2", "BJet btag weight Stats2 Up/Down, 0.0<|#eta|<1.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin3" : hist2D("BJet_BTagWeightStats2_Bin0_Bin3", "BJet btag weight Stats2 Up/Down, 0.0<|#eta|<1.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin0" : hist2D("BJet_BTagWeightStats2_Bin1_Bin0", "BJet btag weight Stats2 Up/Down, 1.5<|#eta|<2.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin1" : hist2D("BJet_BTagWeightStats2_Bin1_Bin1", "BJet btag weight Stats2 Up/Down, 1.5<|#eta|<2.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin2" : hist2D("BJet_BTagWeightStats2_Bin1_Bin2", "BJet btag weight Stats2 Up/Down, 1.5<|#eta|<2.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin3" : hist2D("BJet_BTagWeightStats2_Bin1_Bin3", "BJet btag weight Stats2 Up/Down, 1.5<|#eta|<2.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            },

        },

    "LJet" : {
        "Pt_resolution_2D"   : {
            "Bin0" : hist2D("LJet_Pt_Resolution_Bin0", "Jet (udcs) p_{T} response, 0.0<|#eta|<1.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin1" : hist2D("LJet_Pt_Resolution_Bin1", "Jet (udcs) p_{T} response, 1.5<|#eta|<2.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin2" : hist2D("LJet_Pt_Resolution_Bin2", "Jet (udcs) p_{T} response, 2.5<|#eta|<4.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            },
        "NumCHadrons" : {
            "Bin0" : hist1D("LJet_numCHadrons_Bin0", "Number of c hadrons in l-jets, 0.0<|#eta|<1.5", "number of c hadrons", 5,-1,4),
            "Bin1" : hist1D("LJet_numCHadrons_Bin1", "Number of c hadrons in l-jets, 1.5<|#eta|<2.5", "number of c hadrons", 5,-1,4),
            "Bin2" : hist1D("LJet_numCHadrons_Bin2", "Number of c hadrons in l-jets, 2.5<|#eta|<4.5", "number of c hadrons", 5,-1,4),
            },

        "Btag_weight_2D" : {
            "Bin0_Bin0" : hist2D("LJet_BTagWeight_Bin0_Bin0", "LJet btag weight, 0.0<|#eta|<1.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            "Bin0_Bin1" : hist2D("LJet_BTagWeight_Bin0_Bin1", "LJet btag weight, 0.0<|#eta|<1.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            "Bin0_Bin2" : hist2D("LJet_BTagWeight_Bin0_Bin2", "LJet btag weight, 0.0<|#eta|<1.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            "Bin0_Bin3" : hist2D("LJet_BTagWeight_Bin0_Bin3", "LJet btag weight, 0.0<|#eta|<1.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            "Bin1_Bin0" : hist2D("LJet_BTagWeight_Bin1_Bin0", "LJet btag weight, 1.5<|#eta|<2.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            "Bin1_Bin1" : hist2D("LJet_BTagWeight_Bin1_Bin1", "LJet btag weight, 1.5<|#eta|<2.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            "Bin1_Bin2" : hist2D("LJet_BTagWeight_Bin1_Bin2", "LJet btag weight, 1.5<|#eta|<2.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            "Bin1_Bin3" : hist2D("LJet_BTagWeight_Bin1_Bin3", "LJet btag weight, 1.5<|#eta|<2.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "Weight",  50, 0.5, 1.8) ,
            },

        "Btag_weightJES_2D" : {
            "Bin0_Bin0" : hist2D("LJet_BTagWeightJES_Bin0_Bin0", "LJet btag weight JES Up/Down, 0.0<|#eta|<1.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin1" : hist2D("LJet_BTagWeightJES_Bin0_Bin1", "LJet btag weight JES Up/Down, 0.0<|#eta|<1.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin2" : hist2D("LJet_BTagWeightJES_Bin0_Bin2", "LJet btag weight JES Up/Down, 0.0<|#eta|<1.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin3" : hist2D("LJet_BTagWeightJES_Bin0_Bin3", "LJet btag weight JES Up/Down, 0.0<|#eta|<1.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin0" : hist2D("LJet_BTagWeightJES_Bin1_Bin0", "LJet btag weight JES Up/Down, 1.5<|#eta|<2.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin1" : hist2D("LJet_BTagWeightJES_Bin1_Bin1", "LJet btag weight JES Up/Down, 1.5<|#eta|<2.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin2" : hist2D("LJet_BTagWeightJES_Bin1_Bin2", "LJet btag weight JES Up/Down, 1.5<|#eta|<2.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin3" : hist2D("LJet_BTagWeightJES_Bin1_Bin3", "LJet btag weight JES Up/Down, 1.5<|#eta|<2.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            },

        "Btag_weightHF_2D" : {
            "Bin0_Bin0" : hist2D("LJet_BTagWeightHF_Bin0_Bin0", "LJet btag weight HF Up/Down, 0.0<|#eta|<1.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin1" : hist2D("LJet_BTagWeightHF_Bin0_Bin1", "LJet btag weight HF Up/Down, 0.0<|#eta|<1.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin2" : hist2D("LJet_BTagWeightHF_Bin0_Bin2", "LJet btag weight HF Up/Down, 0.0<|#eta|<1.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin3" : hist2D("LJet_BTagWeightHF_Bin0_Bin3", "LJet btag weight HF Up/Down, 0.0<|#eta|<1.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin0" : hist2D("LJet_BTagWeightHF_Bin1_Bin0", "LJet btag weight HF Up/Down, 1.5<|#eta|<2.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin1" : hist2D("LJet_BTagWeightHF_Bin1_Bin1", "LJet btag weight HF Up/Down, 1.5<|#eta|<2.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin2" : hist2D("LJet_BTagWeightHF_Bin1_Bin2", "LJet btag weight HF Up/Down, 1.5<|#eta|<2.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin3" : hist2D("LJet_BTagWeightHF_Bin1_Bin3", "LJet btag weight HF Up/Down, 1.5<|#eta|<2.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            },

        "Btag_weightStats1_2D" : {
            "Bin0_Bin0" : hist2D("LJet_BTagWeightStats1_Bin0_Bin0", "LJet btag weight Stats1 Up/Down, 0.0<|#eta|<1.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin1" : hist2D("LJet_BTagWeightStats1_Bin0_Bin1", "LJet btag weight Stats1 Up/Down, 0.0<|#eta|<1.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin2" : hist2D("LJet_BTagWeightStats1_Bin0_Bin2", "LJet btag weight Stats1 Up/Down, 0.0<|#eta|<1.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin3" : hist2D("LJet_BTagWeightStats1_Bin0_Bin3", "LJet btag weight Stats1 Up/Down, 0.0<|#eta|<1.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin0" : hist2D("LJet_BTagWeightStats1_Bin1_Bin0", "LJet btag weight Stats1 Up/Down, 1.5<|#eta|<2.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin1" : hist2D("LJet_BTagWeightStats1_Bin1_Bin1", "LJet btag weight Stats1 Up/Down, 1.5<|#eta|<2.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin2" : hist2D("LJet_BTagWeightStats1_Bin1_Bin2", "LJet btag weight Stats1 Up/Down, 1.5<|#eta|<2.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin3" : hist2D("LJet_BTagWeightStats1_Bin1_Bin3", "LJet btag weight Stats1 Up/Down, 1.5<|#eta|<2.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            },

        "Btag_weightStats2_2D" : {
            "Bin0_Bin0" : hist2D("LJet_BTagWeightStats2_Bin0_Bin0", "LJet btag weight Stats2 Up/Down, 0.0<|#eta|<1.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin1" : hist2D("LJet_BTagWeightStats2_Bin0_Bin1", "LJet btag weight Stats2 Up/Down, 0.0<|#eta|<1.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin2" : hist2D("LJet_BTagWeightStats2_Bin0_Bin2", "LJet btag weight Stats2 Up/Down, 0.0<|#eta|<1.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin0_Bin3" : hist2D("LJet_BTagWeightStats2_Bin0_Bin3", "LJet btag weight Stats2 Up/Down, 0.0<|#eta|<1.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin0" : hist2D("LJet_BTagWeightStats2_Bin1_Bin0", "LJet btag weight Stats2 Up/Down, 1.5<|#eta|<2.5 and 20<p_{T}<30",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin1" : hist2D("LJet_BTagWeightStats2_Bin1_Bin1", "LJet btag weight Stats2 Up/Down, 1.5<|#eta|<2.5 and 30<p_{T}<50",  "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin2" : hist2D("LJet_BTagWeightStats2_Bin1_Bin2", "LJet btag weight Stats2 Up/Down, 1.5<|#eta|<2.5 and 50<p_{T}<100", "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            "Bin1_Bin3" : hist2D("LJet_BTagWeightStats2_Bin1_Bin3", "LJet btag weight Stats2 Up/Down, 1.5<|#eta|<2.5 and 100<p_{T}",    "CSV",  20,  -0.1, 1.1, "max(Variation)/Nominal",  50, 0.5, 1.8) ,
            },

        },

    "GJet" : {
        "Pt_resolution_2D"   : {
            "Bin0" : hist2D("GJet_Pt_Resolution_Bin0", "Jet (g) p_{T} response, 0.0<|#eta|<1.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin1" : hist2D("GJet_Pt_Resolution_Bin1", "Jet (g) p_{T} response, 1.5<|#eta|<2.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            "Bin2" : hist2D("GJet_Pt_Resolution_Bin2", "Jet (g) p_{T} response, 2.5<|#eta|<4.5", "p_{T} gen [GeV]", 30,  20, 400, "(Reco-Gen)/Gen", 40, -1, 1) ,
            },
        },

}

for histo in ["ByGenJets", "ByHiggsJets", "Compare_2D"]:
    hists["Classification"]["Gen"][histo].GetXaxis().SetBinLabel(1, "ll")
    hists["Classification"]["Gen"][histo].GetXaxis().SetBinLabel(2, "bb, bB, BB")
    hists["Classification"]["Gen"][histo].GetXaxis().SetBinLabel(3, "bc, bl")
    hists["Classification"]["Gen"][histo].GetXaxis().SetBinLabel(4, "Bc, Bl")
    hists["Classification"]["Gen"][histo].GetXaxis().SetBinLabel(5, "cc, cC, CC")
    hists["Classification"]["Gen"][histo].GetXaxis().SetBinLabel(6, "cl")
    hists["Classification"]["Gen"][histo].GetXaxis().SetBinLabel(7, "Cl")
    if "Compare" in histo:
        hists["Classification"]["Gen"][histo].GetYaxis().SetBinLabel(1, "ll")
        hists["Classification"]["Gen"][histo].GetYaxis().SetBinLabel(2, "bb, bB, BB")
        hists["Classification"]["Gen"][histo].GetYaxis().SetBinLabel(3, "bc, bl")
        hists["Classification"]["Gen"][histo].GetYaxis().SetBinLabel(4, "Bc, Bl")
        hists["Classification"]["Gen"][histo].GetYaxis().SetBinLabel(5, "cc, cC, CC")
        hists["Classification"]["Gen"][histo].GetYaxis().SetBinLabel(6, "cl")
        hists["Classification"]["Gen"][histo].GetYaxis().SetBinLabel(7, "Cl")
        

out.cd()
for k in hists.keys(): 
    out.mkdir(k)
    for k2 in hists[k].keys():
        out.mkdir(k+"/"+k2)
out.cd()

print chain.GetEntries()
for iev in range( min(100000, chain.GetEntries()) ):

    chain.GetEntry(iev)
    ev = chain
    if iev%1000 == 0:
        print "Processing event ", iev
    
    if not pass_event(ev):
        continue

    count_genjet_b  = 0
    count_genjet_2b = 0
    count_genjet_c  = 0
    count_genjet_2c = 0
    for j in range(ev.nGenJet):
        if not( ev.GenJet_pt[j]>20. and abs(ev.GenJet_eta[j])<2.5 ):
            continue
        if ev.GenJet_numBHadrons[j]==1:
            count_genjet_b  += 1
        elif ev.GenJet_numBHadrons[j]>1:
            count_genjet_2b += 1
        elif ev.GenJet_numCHadrons[j]==1:
            count_genjet_c  += 1
        elif ev.GenJet_numCHadrons[j]>1:
            count_genjet_2c += 1

    gen_cat_count = category(count_genjet_b, count_genjet_2b, count_genjet_c, count_genjet_2c)
    hists["Classification"]["Gen"]["ByGenJets"].Fill(gen_cat_count, weight(ev) )

    count_genjet_b  = 0
    count_genjet_2b = 0
    count_genjet_c  = 0
    count_genjet_2c = 0
    for j in [ ev.hJidx[0],  ev.hJidx[1] ]:
        if j>=0:
            mcIdx = ev.Jet_mcIdx[j] 
            if mcIdx>=0 and ev.Jet_mcPt[j]>20.:
                if not( ev.GenJet_pt[mcIdx]>20. and abs(ev.GenJet_eta[mcIdx])<2.5 ):
                    continue
                if ev.GenJet_numBHadrons[mcIdx]==1:
                    count_genjet_b  += 1
                elif ev.GenJet_numBHadrons[mcIdx]>1:
                    count_genjet_2b += 1
                elif ev.GenJet_numCHadrons[mcIdx]==1:
                    count_genjet_c  += 1
                elif ev.GenJet_numCHadrons[mcIdx]>1:
                    count_genjet_2c += 1
    gen_cat_hjets = category(count_genjet_b, count_genjet_2b, count_genjet_c, count_genjet_2c)
    hists["Classification"]["Gen"]["ByHiggsJets"].Fill(gen_cat_hjets, weight(ev))

    hists["Classification"]["Gen"]["Compare_2D"].Fill(gen_cat_count, gen_cat_hjets, weight(ev))
    if ev.V_pt>100 and ev.H_pt>100:                
        hists["Classification"]["Gen"]["Compare_HVPt100_2D"].Fill(gen_cat_count, gen_cat_hjets, weight(ev))
    if ev.Jet_btagCSV[ev.hJidx[0]]>0.890 and ev.Jet_btagCSV[ev.hJidx[1]]>0.890:
        hists["Classification"]["Gen"]["Compare_H_btag_2D"].Fill(gen_cat_count, gen_cat_hjets, weight(ev))

    genPtBin = findbin_genMetPt(ev.met_genPt)

    hists["MET"]["Pt_systematics"]["JetEnUp"+genPtBin].Fill(ev.met_shifted_JetEnUp_pt/ev.met_pt if ev.met_pt>0 else 0., weight(ev))
    hists["MET"]["Pt_systematics"]["JetEnDown"+genPtBin].Fill(ev.met_shifted_JetEnDown_pt/ev.met_pt if ev.met_pt>0 else 0., weight(ev))
    hists["MET"]["Pt_systematics"]["JetResUp"+genPtBin].Fill(ev.met_shifted_JetResUp_pt/ev.met_pt if ev.met_pt>0 else 0., weight(ev))
    hists["MET"]["Pt_systematics"]["JetResDown"+genPtBin].Fill(ev.met_shifted_JetResDown_pt/ev.met_pt if ev.met_pt>0 else 0., weight(ev))
    hists["MET"]["Pt_systematics"]["UnclusteredEnUp"+genPtBin].Fill(ev.met_shifted_UnclusteredEnUp_pt/ev.met_pt if ev.met_pt>0 else 0., weight(ev))
    hists["MET"]["Pt_systematics"]["UnclusteredEnDown"+genPtBin].Fill(ev.met_shifted_UnclusteredEnDown_pt/ev.met_pt if ev.met_pt>0 else 0., weight(ev))

    hists["MET"]["Pt_response"]["raw"+genPtBin].Fill(ev.met_rawpt/ev.met_genPt if ev.met_genPt>0 else 0., weight(ev))
    hists["MET"]["Pt_response"]["type1"+genPtBin].Fill(ev.met_pt/ev.met_genPt if ev.met_genPt>0 else 0., weight(ev))
    hists["MET"]["Pt_response"]["type1p2"+genPtBin].Fill(ev.metType1p2_pt/ev.met_genPt if ev.met_genPt>0 else 0., weight(ev))
    hists["MET"]["Pt_response"]["puppi"+genPtBin].Fill(ev.metPuppi_pt/ev.met_genPt if ev.met_genPt>0 else 0., weight(ev))
    hists["MET"]["Pt_response"]["tkMet"+genPtBin].Fill(ev.tkMet_pt/ev.met_genPt if ev.met_genPt>0 else 0., weight(ev))

    hists["MET"]["PxPy_correlation"]["type1"+genPtBin].Fill( ev.met_pt*ROOT.TMath.Cos(ev.met_phi)  - ev.met_genPt*ROOT.TMath.Cos(ev.met_genPhi),
                                                            ev.met_pt*ROOT.TMath.Sin(ev.met_phi)  - ev.met_genPt*ROOT.TMath.Sin(ev.met_genPhi), weight(ev))
    hists["MET"]["Phi_resolution_2D"]["type1"].Fill( ev.met_genPt, ev.met_phi-ev.met_genPhi , weight(ev))
    hists["MET"]["Px_resolution_vsSumET_2D"]["type1"].Fill( ev.met_sumEt, ev.met_pt*ROOT.TMath.Cos(ev.met_phi) - ev.met_genPt*ROOT.TMath.Cos(ev.met_genPhi) )


    for j in range(ev.nJet):
        # reco jet
        pt   = ev.Jet_pt[j]
        eta  = ev.Jet_eta[j]
        phi  = ev.Jet_phi[j]
        mass = ev.Jet_mass[j]
        hadronFlavour = ev.Jet_hadronFlavour[j]
        partonFlavour = ev.Jet_mcFlavour[j]
        mcIdx = ev.Jet_mcIdx[j] 
        puId  = ev.Jet_puId[j]
        pt_reg = ev.Jet_pt_reg[j]
        pt_raw = ev.Jet_rawPt[j]
        pt_corrJES     = ev.Jet_corr[j]
        pt_corrJESUp   = ev.Jet_corr_JECUp[j]
        pt_corrJESDown = ev.Jet_corr_JECDown[j]
        pt_corrJER     = 1.
        pt_corrJERUp   = 1.
        pt_corrJERDown = 1.
        if hasattr(ev, "Jet_corr_JER"):
            pt_corrJER     = ev.Jet_corr_JER[j]
            pt_corrJERUp   = ev.Jet_corr_JERUp[j]
            pt_corrJERDown = ev.Jet_corr_JERDown[j]

        if debug:
            print ev.Jet_mcPt[j], " <=> " , ev.GenJet_pt[mcIdx]        

        # gen jet
        pt_mc   = ev.Jet_mcPt[j]
        eta_mc  = ev.Jet_mcEta[j]
        phi_mc  = ev.Jet_mcPhi[j]
        mass_mc = ev.Jet_mcM[j]

        numBHadrons   = -1
        numCHadrons   = -1
        
        pt_mc_wNu = -1
        if mcIdx >= 0 and mcIdx<ev.nGenJet and pt_mc>20.:
            numBHadrons = ev.GenJet_numBHadrons[mcIdx]
            numCHadrons = ev.GenJet_numCHadrons[mcIdx]
            pt_mc_wNu   = ev.GenJet_wNuPt[mcIdx] 

        if pt>20. and abs(eta)<4.5:
            binEta = findbin_jetEta(eta)
            binPt  = findbin_jetPt(pt)

            hists["Jet"]["JES_2D"][binEta].Fill( pt_raw, pt_corrJES, weight(ev))
            hists["Jet"]["JESUp_2D"][binEta].Fill( pt_raw, pt_corrJESUp/pt_corrJES, weight(ev))
            hists["Jet"]["JESDown_2D"][binEta].Fill( pt_raw, pt_corrJESDown/pt_corrJES, weight(ev))
            hists["Jet"]["JER_2D"][binEta].Fill( pt, pt_corrJER, weight(ev))
            hists["Jet"]["JERUp_2D"][binEta].Fill( pt, pt_corrJERUp/pt_corrJER, weight(ev))
            hists["Jet"]["JERDown_2D"][binEta].Fill( pt, pt_corrJERDown/pt_corrJER, weight(ev))

            if pt_mc>0.:
                hists["Jet"]["PuId_match_2D"][binEta].Fill( pt, max(puId, 0.), weight(ev))
            else:
                hists["Jet"]["PuId_unmatch_2D"][binEta].Fill( pt, max(puId, 0.), weight(ev))


            if pt_mc>0.:
                hists["Jet"]["Pt_resolution_2D"][binEta].Fill( pt_mc,  (pt-pt_mc)/pt_mc, weight(ev))
                if ev.Jet_leptonPt[j]>0.:
                    if abs(ev.Jet_leptonPdgId[j])==11:
                        hists["Jet"]["Pt_resolution_wElectron_2D"][binEta].Fill( pt_mc,  (pt-pt_mc)/pt_mc, weight(ev))
                    if abs(ev.Jet_leptonPdgId[j])==13:
                        hists["Jet"]["Pt_resolution_wMuon_2D"][binEta].Fill( pt_mc,  (pt-pt_mc)/pt_mc, weight(ev))

                hists["Jet"]["Eta_resolution_2D"][binEta].Fill( pt_mc, (eta-eta_mc), weight(ev))
                hists["Jet"]["Phi_resolution_2D"][binEta].Fill( pt_mc, (phi-phi_mc), weight(ev))
                if pt_reg>0.:
                    hists["RegressedJet"]["Pt_resolution_2D"][binEta].Fill( pt_mc,      (pt_reg-pt_mc)/pt_mc, weight(ev))
                    hists["RegressedJet"]["Pt_resolution_wNu_2D"][binEta].Fill( pt_mc,  (pt_reg-pt_mc_wNu)/pt_mc_wNu, weight(ev))
                    
            if abs(hadronFlavour)==5:
                if pt_mc>0.:
                    hists["BJet"]["Pt_resolution_2D"][binEta].Fill( pt_mc,  (pt-pt_mc)/pt_mc, weight(ev))
                hists["BJet"]["NumBHadrons"][binEta].Fill(numBHadrons)
                if binEta in ["Bin0", "Bin1"]:
                    hists["BJet"]["Btag_weight_2D"]   [binEta+"_"+binPt].Fill( max(ev.Jet_btagCSV[j], -0.1), ev.Jet_bTagWeight[j], weight(ev) )
                    for syst in ["JES", "LF", "Stats1", "Stats2"]:
                        hists["BJet"]["Btag_weight"+syst+"_2D"][binEta+"_"+binPt].Fill( max(ev.Jet_btagCSV[j], -0.1), 1+max( abs(1-getattr(ev,"Jet_bTagWeight"+syst+"Up")[j]  /ev.Jet_bTagWeight[j]), 
                                                                                                                             abs(1-getattr(ev,"Jet_bTagWeight"+syst+"Down")[j]/ev.Jet_bTagWeight[j])) , weight(ev))

            if abs(hadronFlavour)<5:
                if pt_mc>0.:
                    if partonFlavour != 21:
                        hists["LJet"]["Pt_resolution_2D"][binEta].Fill( pt_mc,  (pt-pt_mc)/pt_mc, weight(ev))
                    else:
                        hists["GJet"]["Pt_resolution_2D"][binEta].Fill( pt_mc,  (pt-pt_mc)/pt_mc, weight(ev))

                hists["LJet"]["NumCHadrons"][binEta].Fill(numBHadrons)
                if binEta in ["Bin0", "Bin1"]:
                    hists["LJet"]["Btag_weight_2D"]   [binEta+"_"+binPt].Fill( max(ev.Jet_btagCSV[j], -0.1), ev.Jet_bTagWeight[j], weight(ev) )
                    for syst in ["JES", "HF", "Stats1", "Stats2"]:
                        hists["LJet"]["Btag_weight"+syst+"_2D"][binEta+"_"+binPt].Fill( max(ev.Jet_btagCSV[j], -0.1), 1+max( abs(1-getattr(ev,"Jet_bTagWeight"+syst+"Up")[j]  /ev.Jet_bTagWeight[j]), 
                                                                                                                             abs(1-getattr(ev,"Jet_bTagWeight"+syst+"Down")[j]/ev.Jet_bTagWeight[j])) , weight(ev))

            #if mass_mc>0.:
            #    hists["Jet"]["mass_resolution"][binEta].Fill( pt_mc, (mass-mass_mc)/mass_mc)

        lv = ROOT.TLorentzVector()
        lv.SetPtEtaPhiM(pt,eta,phi,mass)
        if debug:
            print "\tJet %s: (%s, %s, %s, %s)" % (j, lv.Px(), lv.Py(), lv.Pz(), lv.E() )
            


print "Normalising histograms..."
for k in hists.keys():        
    for k2 in hists[k].keys():
        for k3 in hists[k][k2].keys():
            h = hists[k][k2][k3]
            out.cd(k+"/"+k2)
            if not ("2D" in k2 or "2D" in k3):
                h.Write("", ROOT.TObject.kOverwrite)
                continue
            hp = h.ProfileX("_pfx", 1, -1, "")   
            hp.GetYaxis().SetTitle( "mean of ("+h.GetYaxis().GetTitle()+")" )
            nbinsX = h.GetXaxis().GetNbins()
            nbinsY = h.GetYaxis().GetNbins()
            for binx in range(nbinsX):
                norm = 0.
                for biny in range(nbinsY):
                    norm += h.GetBinContent(binx+1, biny+1)
                if norm>0.:
                    for biny in range(nbinsY):
                        old = h.GetBinContent(binx+1, biny+1)
                        h.SetBinContent(binx+1, biny+1, old/norm)
            h.Write("", ROOT.TObject.kOverwrite)
            hp.Write("", ROOT.TObject.kOverwrite)

#out.Write("", ROOT.TObject.kOverwrite)
out.Close()  
