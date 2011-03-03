

void controlPlots(string strigger="Jet30"){

  gStyle->SetOptStat(0);
  gROOT->ProcessLine(".L ~/tdrstyle.C");
  setTDRStyle();
  //style modifications:
  gStyle->SetErrorX(0.5);
  gStyle->SetPadTopMargin(0.08);//.005
  gStyle->SetPadBottomMargin(0.10);//0.13
  gStyle->SetMarkerStyle(1);

  TCanvas *c = new TCanvas(("ptCutComp"+strigger).c_str(),("pt cut comparison (8 vs. 10 GeV), "+strigger).c_str(),100,300,500,500);
  TCanvas *cr = new TCanvas(("ptCutCompRat"+strigger).c_str(),("pt cut comparison (8 vs. 10 GeV) Ratio, "+strigger).c_str(),500,300,500,500);
  c->cd(); 

  //string strigger="jet30"; 
  double rebin= 5.0; 
  string file1 = "PT8vsPT10/ptcut8/histo-Data-HLT_" + strigger + "U-BBCorr_HadrJetMatchPt_OpenHLT-total.root"; 
  string file2 = "PT8vsPT10/ptcut10/histo-Data-HLT_" + strigger + "U-BBCorr_HadrJetMatchPt_OpenHLT-total.root";
  string plot="Vert_dR_2b_CORR";
  string label1 = "p_{T} > 8 GeV";   
  string label2 = "p_{T} > 10 GeV";   
  //string plot="Vert_dR_2b";  

  TFile *f1 = TFile::Open(file1.c_str(),"READ");
  TFile *f2 = TFile::Open(file2.c_str(),"READ");
  
  TH1F *hp1; f1->GetObject((plot).c_str(),hp1);
  TH1F *hp2; f2->GetObject((plot).c_str(),hp2);

  hp1->Rebin(rebin);
  hp2->Rebin(rebin);
  
  hp1->SetLineColor(1);
  hp1->SetMarkerColor(1);
  hp1->SetMarkerStyle(4);
  hp2->SetLineColor(2);
  hp2->SetMarkerColor(2);
  hp2->SetMarkerStyle(20);
  hp2->SetMarkerSize(1.13);
  

  bool hp1first = false;
  if(hp1->GetMaximum()>hp2->GetMaximum()) hp1first=true;
  double max; 
  if(hp1first){
    hp1->GetXaxis()->SetTitle("#Delta R"); 
    hp1->GetYaxis()->SetTitle("d#sigma /d#Delta R"); 
    hp1->GetYaxis()->SetTitleOffset(1.9); 
    hp1->GetXaxis()->SetRangeUser(0,3.9); 
    hp1->DrawCopy("PE");
    hp2->DrawCopy("PEsames");
    hp1->DrawCopy("PEsames");
    max = hp1->GetMaximum(); 
  }
  else{
    hp2->DrawCopy("PE");
    hp2->GetXaxis()->SetTitle("#Delta R"); 
    hp2->GetYaxis()->SetTitle("d#sigma /d#Delta R"); 
    hp2->GetYaxis()->SetTitleOffset(1.9); 
    hp2->GetXaxis()->SetRangeUser(0,3.9); 
    hp1->DrawCopy("PEsames");
    max = hp2->GetMaximum(); 
  }

  TLegend *leg = new TLegend(0.6,0.65,0.89,0.84);
  leg->SetHeader(("HLT_"+strigger+"U").c_str());
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(hp1,label1.c_str(),"PE");
  leg->AddEntry(hp2,label2.c_str(),"PE");
  leg->Draw();

  TLatex *   tex; 
  if(strigger=="Jet15") tex = new TLatex(0,max*1.26,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
  else if(strigger=="Jet30") tex = new TLatex(0,max*1.17,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
  else if(strigger=="Jet50") tex = new TLatex(0,max*1.13,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
  tex->SetTextSize(0.040); //0.044
  tex->SetLineWidth(2);
  tex->Draw();

  gPad->RedrawAxis(); 

  cr->cd(); 
  hp2->Rebin(2); 
  hp1->Rebin(2); 
  TH1F *hrat = hp2->Clone(); 
  hrat->Divide(hp2,hp1,1,1);
  hrat->GetXaxis()->SetRangeUser(0,3.9); 
  hrat->GetXaxis()->SetTitle("#Delta R"); 
  hrat->GetYaxis()->SetTitle("ratio"); 
  hrat->GetYaxis()->SetTitleOffset(1.7); 
  hrat->SetLineColor(1);
  hrat->SetMarkerColor(1);  
  hrat->Draw(); 
  TLatex *tex2; 
  if(strigger=="Jet15") tex2 = new TLatex(0,1.53,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
  else if(strigger=="Jet30") tex2 = new TLatex(0,1.187,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
  else if(strigger=="Jet50") tex2 = new TLatex(0,1.143,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
  tex2->SetTextSize(0.040); //0.044
  tex2->SetLineWidth(2);
  tex2->Draw();

}
