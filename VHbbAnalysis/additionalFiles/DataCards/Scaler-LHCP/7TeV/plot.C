{
RooDataHist *RDH = (RooDataHist*) ZeeLoose_7TeV->data("ZjLFCMS_res_jUp");
RooRealVar BDT("CMS_vhbb_BDT_Zll","CMS_vhbb_BDT_Zll",-1,1) ;
RooPlot* xframe = BDT.frame() ;
RDH->plotOn(xframe) ;
xframe.Draw();
}
