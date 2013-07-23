{
RooDataHist *RDH = (RooDataHist*) ZnunuHighPt_8TeV->data("ZH");
RooRealVar BDT("CMS_vhbb_BDT_ZnunuHighPt_8TeV","CMS_vhbb_BDT_ZnunuHighPt_8TeV",-1,1) ;
RooPlot* xframe = BDT.frame() ;
RDH->plotOn(xframe) ;
xframe.Draw();
}
