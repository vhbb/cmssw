//use mybbbarcorr.C
// init() first!!!
//ivfDXbyflavor("dR","rootoutput/tempBDCOMPplots/defaultAddedBinsWITHOUTBDCleaning.root", 2, 1.0)
//ivfDXbyflavor("dR","defaultAddedBinsWITHOUTBDCleaning.root", 2, 1.0)
//ivfDXbyflavor("dR","rootoutput/tempBDCOMPplots/defaultAddedBinsWITHBDCleaning.root",2,1.0)
//ivfDXbyflavor("dR","defaultAddedBinsWITHBDCleaning.root",2,1.0)


string strigger="jet30"; 

void withBDcleaning(string trig="jet30"){
  strigger=trig; 
  ivfDXbyflavor("dR","BDcomparisonPlots/defaultAddedBinsWITHBDCleaning.root",2,1.0,"withBDcleaning"); 
}
void withoutBDcleaning(string trig="jet30"){
  strigger=trig; 
  ivfDXbyflavor("dR","BDcomparisonPlots/defaultAddedBinsWITHOUTBDCleaning.root", 2, 1.0,"withoutBDcleaning");
}

void ivfDXbyflavor(string todraw="dR", string fileN="rootoutput/demoJetVertStudy.root", int rebin=1, double rescale=1.0, string canv="DXivf"){

  gROOT->ProcessLine(".L ~/tdrstyle.C");
  setTDRStyle();
  gStyle->SetPadTopMargin(0.08);//.005
  gStyle->SetPadBottomMargin(0.10);//0.13

  TFile *f = TFile::Open(fileN.c_str(),"UPDATE");
//   f->cd(strigger.c_str());

  if(todraw=="dR"){
    TH1F *h2b; f->GetObject((strigger+"/Vert_dR_2bVert_2b").c_str(),h2b);
    TH1F *h1b1c; f->GetObject((strigger+"/Vert_dR_2bVert_1b1c").c_str(),h1b1c);
    TH1F *h1b1l; f->GetObject((strigger+"/Vert_dR_2bVert_1b1l").c_str(),h1b1l);
    TH1F *h1c1c; f->GetObject((strigger+"/Vert_dR_2bVert_1c1c").c_str(),h1c1c);
    TH1F *h1c1l; f->GetObject((strigger+"/Vert_dR_2bVert_1c1l").c_str(),h1c1l);
    TH1F *h1l1l; f->GetObject((strigger+"/Vert_dR_2bVert_1l1l").c_str(),h1l1l);
    h2b->GetXaxis()->SetTitle("#Delta R");
  }
  else if(todraw=="dPhi"){
    TH1F *h2b; f->GetObject((strigger+"/Vert_dPhi_2bVert_2b").c_str(),h2b);
    TH1F *h1b1c; f->GetObject((strigger+"/Vert_dPhi_2bVert_1b1c").c_str(),h1b1c);
    TH1F *h1b1l; f->GetObject((strigger+"/Vert_dPhi_2bVert_1b1l").c_str(),h1b1l);
    TH1F *h1c1c; f->GetObject((strigger+"/Vert_dPhi_2bVert_1c1c").c_str(),h1c1c);
    TH1F *h1c1l; f->GetObject((strigger+"/Vert_dPhi_2bVert_1c1l").c_str(),h1c1l);
    TH1F *h1l1l; f->GetObject((strigger+"/Vert_dPhi_2bVert_1l1l").c_str(),h1l1l);
    h2b->GetXaxis()->SetTitle("#Delta #phi");
  }
  else if(todraw=="dEta"){
    TH1F *h2b; f->GetObject((strigger+"/Vert_dEta_2bVert_2b").c_str(),h2b);
    TH1F *h1b1c; f->GetObject((strigger+"/Vert_dEta_2bVert_1b1c").c_str(),h1b1c);
    TH1F *h1b1l; f->GetObject((strigger+"/Vert_dEta_2bVert_1b1l").c_str(),h1b1l);
    TH1F *h1c1c; f->GetObject((strigger+"/Vert_dEta_2bVert_1c1c").c_str(),h1c1c);
    TH1F *h1c1l; f->GetObject((strigger+"/Vert_dEta_2bVert_1c1l").c_str(),h1c1l);
    TH1F *h1l1l; f->GetObject((strigger+"/Vert_dEta_2bVert_1l1l").c_str(),h1l1l);
    h2b->GetXaxis()->SetTitle("#Delta #eta");
  }

  else return;

  gStyle->SetOptStat(0);
  TCanvas *cDX = new TCanvas(canv.c_str(),("IVF DX by flavor " + canv).c_str(),100,100,700,500);
//   cDX->SetLogy();

  int nbb = h2b->Integral();
  int nbc = h1b1c->Integral();
  int nbl = h1b1l->Integral();
  int ncc = h1c1c->Integral();
  int ncl = h1c1l->Integral();
  int nll = h1l1l->Integral();

  std::stringstream ssnbb; ssnbb << nbb;
  std::stringstream ssnbc; ssnbc << nbc;
  std::stringstream ssnbl; ssnbl << nbl;
  std::stringstream ssncc; ssncc << ncc;
  std::stringstream ssncl; ssncl << ncl;
  std::stringstream ssnll; ssnll << nll;

  h1c1l->Add(h1l1l);
  h1c1c->Add(h1c1l);
  h1b1l->Add(h1c1c);
  h1b1c->Add(h1b1l);
  h2b->Add(h1b1c);


  h2b->SetFillColor(8);
  h1b1c->SetFillColor(7);
  h1b1l->SetFillColor(2);
  h1c1c->SetFillColor(6);
  h1c1l->SetFillColor(4);
  h1l1l->SetFillColor(1);

  h2b->SetLineColor(8);
  h1b1c->SetLineColor(7);
  h1b1l->SetLineColor(2);
  h1c1c->SetLineColor(6);
  h1c1l->SetLineColor(4);
  h1l1l->SetLineColor(1);

  if(rebin!=1){
    h2b->Rebin(rebin);
    h1b1c->Rebin(rebin);
    h1b1l->Rebin(rebin);
    h1c1c->Rebin(rebin);
    h1c1l->Rebin(rebin);
    h1l1l->Rebin(rebin);
  }

  h2b->Scale(rescale);
  h1b1c->Scale(rescale);
  h1b1l->Scale(rescale);
  h1c1c->Scale(rescale);
  h1c1l->Scale(rescale);
  h1l1l->Scale(rescale);

  h2b->GetYaxis()->SetRangeUser(0, h2b->GetMaximum()*1.1);
  h2b->SetTitle(""); 
  h2b->GetYaxis()->SetTitle("number of selected B candidate pairs"); 

  h2b->Draw("hist");
  h1b1c->Draw("histsame");
  h1b1l->Draw("histsame");
  h1c1c->Draw("histsame");
  h1c1l->Draw("histsame");
  h1l1l->Draw("histsame");

  TLegend *leg = new TLegend(0.7,0.5,0.89,0.89);
  leg->SetFillColor(0);
  leg->SetBorderSize(0); 
  leg->AddEntry(h2b,("bb: " + ssnbb.str()).c_str(),"F");
  leg->AddEntry(h1b1c,("bc: " + ssnbc.str()).c_str(),"F");
  leg->AddEntry(h1b1l,("bl: " + ssnbl.str()).c_str(),"F");
  leg->AddEntry(h1c1c,("cc: " + ssncc.str()).c_str(),"F");
  leg->AddEntry(h1c1l,("cl: " + ssncl.str()).c_str(),"F");
  leg->AddEntry(h1l1l,("ll: " + ssnll.str()).c_str(),"F");
  leg->Draw();

  double axismax;
  double axismax  = h2b->GetMaximum();
  
  TLatex *   tex = new TLatex(0.0,axismax*1.05,"CMS    #sqrt{s} = 7 TeV, Simulation");
  tex->SetTextSize(0.040); //0.044
  tex->SetLineWidth(2);
  tex->Draw();
  
  gPad->RedrawAxis();

}


