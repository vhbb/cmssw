{
RooDataHist *RDH = (RooDataHist*) ZeeLowPt_8TeV->data("ZjLFCMS_res_jUp");
RooRealVar BDT("CMS_vhbb_BDT_Zll_8TeV","CMS_vhbb_BDT_Zll_8TeV",-1,1) ;
RooPlot* xframe = BDT.frame() ;
RDH->plotOn(xframe) ;
xframe.Draw();
}
