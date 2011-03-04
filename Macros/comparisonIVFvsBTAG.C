void controlPlot(string strigger="Jet30"){
  gStyle->SetOptStat(0);
  gROOT->ProcessLine(".L ~/tdrstyle.C");
  setTDRStyle();
  //style modifications:
  gStyle->SetErrorX(0.5);
  gStyle->SetPadTopMargin(0.08);//.005
  gStyle->SetPadBottomMargin(0.10);//0.13
  gStyle->SetMarkerStyle(1);

  TCanvas *c = new TCanvas(("ivfVSbtag"+strigger).c_str(),("btag vs. ivf vertices "+strigger).c_str(),100,300,500,500);

  //string strigger="jet30"; 
  double rebin= 5.0;


  string plot1 = "Jet_dR_2b_CORR"; 
  string plot2 = "Vert_dR_2b_CORR"; 

  string file = "/scratch/leo/finalOpenBins36X_No6U_NewDPhiBins/histo-Data-HLT_" + strigger + "U-BBCorr_HadrJetMatchPt_OpenHLT-total.root"; 
  string label1 = "Jet based";
  string label2 = "IVF based";
  //string plot="Vert_dR_2b";  

  TFile *f1 = TFile::Open(file.c_str(),"READ");
  TFile *f2 = TFile::Open(file.c_str(),"READ");

  TH1F *hp1; f1->GetObject((plot1).c_str(),hp1);
  TH1F *hp2; f2->GetObject((plot2).c_str(),hp2);

  hp1->Rebin(rebin);
  hp2->Rebin(rebin);

  hp1->SetLineColor(4);
  hp1->SetMarkerColor(4);
  hp1->SetMarkerStyle(5);
  hp2->SetLineColor(1);
  hp2->SetMarkerColor(1);
  hp2->SetMarkerStyle(4);
  hp1->SetMarkerSize(1.13);

  bool hp1first = false;
  if(hp1->GetMaximum()>hp2->GetMaximum()) hp1first=true;
  double max = -1; 
  for(unsigned int b=2; b<11; b++){
    if(hp1->GetBinContent(b)>max){
      hp1first=true; max = hp1->GetBinContent(b); 
    }
    if(hp2->GetBinContent(b)>max){
      hp1first=false; max = hp2->GetBinContent(b); 
    }
  }
  std::cout << "max " << max << "\n";

//   double max;
  if(hp1first){
    hp1->GetXaxis()->SetTitle("#Delta R");
    hp1->GetYaxis()->SetTitle("d#sigma /d#Delta R");
    hp1->GetYaxis()->SetTitleOffset(1.9);
//     hp1->GetXaxis()->SetRangeUser(0,3.9);
    hp1->GetYaxis()->SetRangeUser(0,max*1.4);
    hp1->DrawCopy("PE");
    hp2->DrawCopy("PEsames");
//     hp1->DrawCopy("PEsames");
//     max = hp1->GetMaximum();
  }
  else{
    hp2->GetXaxis()->SetTitle("#Delta R");
    hp2->GetYaxis()->SetTitle("d#sigma /d#Delta R");
    hp2->GetYaxis()->SetTitleOffset(1.9);
//     hp2->GetXaxis()->SetRangeUser(0,3.9);
    hp2->GetYaxis()->SetRangeUser(0,max*1.4);
    hp2->DrawCopy("PE");
    hp1->DrawCopy("PEsames");
    hp2->DrawCopy("PEsames");
//     max = hp2->GetMaximum();
  }

  TLegend *leg = new TLegend(0.65,0.68,0.94,0.88);
  leg->SetHeader(("HLT_"+strigger+"U").c_str());
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(hp1,label1.c_str(),"PE");
  leg->AddEntry(hp2,label2.c_str(),"PE");
  leg->Draw();

  TLatex *   tex;
  if(strigger=="Jet15") tex = new TLatex(0,max*1.45,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
  else if(strigger=="Jet30") tex = new TLatex(0,max*1.45,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
  else if(strigger=="Jet50") tex = new TLatex(0,max*1.45,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
  tex->SetTextSize(0.040); //0.044
  tex->SetLineWidth(2);
  tex->Draw();

  gPad->RedrawAxis();

}
