//for the moment just line to draw plot using mybbbarcorr.C: 

// ivfDXbyflavor("dR","rootoutput/tempBDCOMPplots/defaultAddedBinsWITHBDCleaning.root",1,1.0); 
// ivfDXbyflavor("dR","rootoutput/tempBDCOMPplots/defaultAddedBinsWITHOUTBDCleaning.root",1,1.0); 


//MADGRAPH

// ivfDXbyflavor("dR","rootoutput/tempMADGRAPH/defaultAddedBins.root",1,1.0); 

void controlPlots(){
  gStyle->SetOptStat(0);

  string dir="/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/BBC/CMSSW_3_7_0_patch2/crab/analyze/rootoutput/tempBDCOMPplots/"; 
  string strigger="jet30/";
  int rebin=4; 

  TFile *fmit = TFile::Open((dir+"defaultAddedBinsWITHBDCleaning.root").c_str()); 
  TFile *fohne= TFile::Open((dir+"defaultAddedBinsWITHOUTBDCleaning.root").c_str()); 

  TH1F *hmitAll; fmit->GetObject((strigger+"Vert_dR_2bVert").c_str(),hmitAll);
  TH1F *hmitBD; fmit->GetObject((strigger+"Vert_dR_2bVert_1b1c").c_str(),hmitBD);
  TH1F *hohneAll; fohne->GetObject((strigger+"Vert_dR_2bVert").c_str(),hohneAll);
  TH1F *hohneBD; fohne->GetObject((strigger+"Vert_dR_2bVert_1b1c").c_str(),hohneBD);

  hmitAll->SetTitle(""); 
  hmitAll->GetXaxis()->SetTitle("#Delta R"); 

  hmitAll->SetLineColor(1); 
  hohneAll->SetLineColor(4); 
  hmitBD->SetLineColor(1); 
  hohneBD->SetLineColor(4); 

  hmitAll->SetFillColor(0); 
  hohneAll->SetFillColor(0); 
  hmitBD->SetFillColor(0); 
  hohneBD->SetFillColor(0); 

  hmitAll->SetLineWidth(2); 
  hohneAll->SetLineWidth(2); 
  hmitBD->SetLineWidth(2); 
  hohneBD->SetLineWidth(2); 

  hmitBD->SetLineStyle(2); 
  hohneBD->SetLineStyle(2); 

  hmitAll->Rebin(rebin); 
  hohneAll->Rebin(rebin); 
  hmitBD->Rebin(rebin); 
  hohneBD->Rebin(rebin); 

  double max = hmitAll->GetMaximum(); 
  if(max<hohneAll->GetMaximum()) max = hohneAll->GetMaximum(); 
  hmitAll->GetYaxis()->SetRangeUser(0.01,max+max/10.0); 

  hmitAll->Draw("hist"); 
  hohneAll->Draw("histsames"); 
  hmitBD->Draw("histsames"); 
  hohneBD->Draw("histsames"); 
  
  //leg
  TLegend *leg = new TLegend(0.35,0.65,0.9,0.9); 
  leg->SetFillColor(0); 
  leg->AddEntry(hmitAll,"identification applied","L"); 
  leg->AddEntry(hohneAll,"before B->D identification","L");
  leg->AddEntry(hmitBD,"events with B and D vertex (ident. applied)","L"); 
  leg->AddEntry(hohneBD,"events with B and D vertex (no ident. applied)","L");
  leg->Draw(); 

}
