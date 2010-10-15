{
  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();

  //const int nPaths = 4; string HLTs[4] = {"HLT_Jet15U","HLT_Jet30U","HLT_Jet50U",""}

  const int nPaths = 1;  string HLTs[nPaths] = {""};
  //const int nMC = 2;  string MCs[nMC] = {"Pythia6","Madgraph"};

  TFile*files[nPaths];
  TH1F *Jet_dPhi_2b[nPaths],*Jet_dPhi[nPaths];
  TMultiGraph *mg_dPhi = new TMultiGraph(); TGraphAsymmErrors *pur_Jet_dPhi[nPaths];
  TH1F *Jet_dR_2b[nPaths],*Jet_dR[nPaths];
  TMultiGraph *mg_dR = new TMultiGraph(); TGraphAsymmErrors *pur_Jet_dR[nPaths];
  TH1F *Jet_pt_2b[nPaths],*Jet_pt[nPaths];
  TMultiGraph *mg_pt = new TMultiGraph(); TGraphAsymmErrors *pur_Jet_pt[nPaths];
  TH1F *Jet_dEta_2b[nPaths],*Jet_dEta[nPaths];
  TMultiGraph *mg_dEta = new TMultiGraph(); TGraphAsymmErrors *pur_Jet_dEta[nPaths];

  TCanvas *c_dEtadPhi[nPaths];
  TH2F *Jet_dEtadPhi_2b[nPaths],*Jet_dEtadPhi[nPaths];
  TCanvas *c_dEtadPhi_clone[nPaths];
  TH2F *Jet_dEtadPhi_2b_clone[nPaths],*Jet_dEtadPhi_clone[nPaths];

  TCanvas *c_pt1pt2[nPaths];
  TH2F *Jet_pt1pt2_2b[nPaths],*Jet_pt1pt2[nPaths];


  TLegend *leg = new TLegend(0.1,0.2,0.3,0.4);

  for(int hlt=0;hlt<nPaths;hlt++){
    string HLT = HLTs[hlt];
    //string HLT = MCs[hlt];

    if(HLT!="") HLT+="-";
    //files[hlt] = TFile::Open(("histo-Madgraph-"+HLT+"total.root").c_str());
    files[hlt] = TFile::Open(("histo-Pythia6-"+HLT+"total.root").c_str());

    files[hlt]->GetObject("Jet_dPhi_2b",Jet_dPhi[hlt]); Jet_dPhi[hlt]->Rebin(2);
    files[hlt]->GetObject("Jet_dPhi_2b_2b",Jet_dPhi_2b[hlt]); Jet_dPhi_2b[hlt]->Rebin(2);
    pur_Jet_dPhi[hlt] = new TGraphAsymmErrors(Jet_dPhi_2b[hlt],Jet_dPhi[hlt]);
    pur_Jet_dPhi[hlt]->SetLineColor(4-hlt);     pur_Jet_dPhi[hlt]->SetMarkerColor(4-hlt);
    mg_dPhi->Add(pur_Jet_dPhi[hlt]);
    string label = HLT;
    if(HLT=="") label = "Total";
    leg->AddEntry(pur_Jet_dPhi[hlt],label.c_str(),"pl");

    files[hlt]->GetObject("Jet_dR_2b",Jet_dR[hlt]); Jet_dR[hlt]->Rebin(2);
    files[hlt]->GetObject("Jet_dR_2b_2b",Jet_dR_2b[hlt]); Jet_dR_2b[hlt]->Rebin(2);
    pur_Jet_dR[hlt] = new TGraphAsymmErrors(Jet_dR_2b[hlt],Jet_dR[hlt]);
    pur_Jet_dR[hlt]->SetLineColor(4-hlt);     pur_Jet_dR[hlt]->SetMarkerColor(4-hlt);
    mg_dR->Add(pur_Jet_dR[hlt]);

    files[hlt]->GetObject("Jet_dEta_2b",Jet_dEta[hlt]); Jet_dEta[hlt]->Rebin(2);
    files[hlt]->GetObject("Jet_dEta_2b_2b",Jet_dEta_2b[hlt]); Jet_dEta_2b[hlt]->Rebin(2);
    pur_Jet_dEta[hlt] = new TGraphAsymmErrors(Jet_dEta_2b[hlt],Jet_dEta[hlt]);
    pur_Jet_dEta[hlt]->SetLineColor(4-hlt);     pur_Jet_dEta[hlt]->SetMarkerColor(4-hlt);
    mg_dEta->Add(pur_Jet_dEta[hlt]);

    files[hlt]->GetObject("Jet_dEtadPhi_2b",Jet_dEtadPhi[hlt]); Jet_dEtadPhi[hlt]->Rebin2D(2,2);
    files[hlt]->GetObject("Jet_dEtadPhi_2b_2b",Jet_dEtadPhi_2b[hlt]); Jet_dEtadPhi_2b[hlt]->Rebin2D(2,2);
    Jet_dEtadPhi_2b_clone[hlt] = (TH2F*)Jet_dEtadPhi_2b[hlt]->Clone();
    Jet_dEtadPhi_2b[hlt]->Divide(Jet_dEtadPhi[hlt]);
    files[hlt]->GetObject("Jet_pt1pt2_2b",Jet_pt1pt2[hlt]); //Jet_pt1pt2[hlt]->Rebin2D(2,2);
    files[hlt]->GetObject("Jet_pt1pt2_2b_2b",Jet_pt1pt2_2b[hlt]); //Jet_pt1pt2_2b[hlt]->Rebin2D(2,2);
    Jet_pt1pt2_2b[hlt]->Divide(Jet_pt1pt2[hlt]);


    files[hlt]->GetObject("Jet_pt_2b",Jet_pt[hlt]); Jet_pt[hlt]->Rebin(1);
    files[hlt]->GetObject("Jet_pt_2b_2b",Jet_pt_2b[hlt]); Jet_pt_2b[hlt]->Rebin(1);
    pur_Jet_pt[hlt] = new TGraphAsymmErrors(Jet_pt_2b[hlt],Jet_pt[hlt]);
    pur_Jet_pt[hlt]->SetLineColor(4-hlt);     pur_Jet_pt[hlt]->SetMarkerColor(4-hlt);
    mg_pt->Add(pur_Jet_pt[hlt]);


  }

  //HLT = "MadgraphVsPythia";
  TCanvas *c_dPhi = new TCanvas( ("BTagPurVsDphi"+HLT).c_str(), ("BTagPurVsDphi"+HLT).c_str());
  mg_dPhi->Draw("ap");
  leg->Draw();

  TCanvas *c_dEta = new TCanvas( ("BTagPurVsDeta"+HLT).c_str(), ("BTagPurVsDeta"+HLT).c_str());
  mg_dEta->Draw("ap");
  leg->Draw();

  TCanvas *c_dR = new TCanvas( ("BTagPurVsDR"+HLT).c_str(), ("BTagPurVsDR"+HLT).c_str());
  mg_dR->Draw("ap");
  leg->Draw();

  /*
  for(int hlt=0;hlt<nMC;hlt++){
    //string HLT = MCs[hlt];
    //if(HLT!="") HLT+="-";
    
    c_dEtadPhi[hlt] = new TCanvas( ("BTagPurVsDetaDphi"+HLT).c_str(), ("BTagPurVsDetadPhi"+HLT).c_str());
    Jet_dEtadPhi_2b[hlt]->Draw("");

    c_dEtadPhi_clone[hlt] = new TCanvas( ("DetaDphi"+HLT).c_str(), ("DetadPhi"+HLT).c_str());
    Jet_dEtadPhi_2b_clone[hlt]->Draw("");

    c_pt1pt2[hlt] = new TCanvas( ("BTagPurVspt1pt2"+HLT).c_str(), ("BTagPurVspt1pt2"+HLT).c_str());
    Jet_pt1pt2_2b[hlt]->Draw("");

  }
  */
  TCanvas *c_pt = new TCanvas( ("BTagPurVsPt"+HLT).c_str(), ("BTagPurVsPt"+HLT).c_str());
  mg_pt->Draw("ap");
  leg->Draw();
  pur_Jet_pt[0]->Fit("pol1","F","sames",30,250);
  TF1 *f_pt_l = pur_Jet_pt[0]->GetFunction("pol1");
  f_pt_l->SetName("f_purity_pt");
  f_pt_l->SetLineColor(kBlack);f_pt_l->SetLineWidth(2);
  

    TFile *_fileOut = TFile::Open("histo-Pythia6-flavourComp.root","RECREATE");
    //TFile *_fileOut = TFile::Open("histo-Pythia6vsMadgraph-flavourComp.root","RECREATE");
  c_dPhi->Write(); c_dEta->Write(); c_dR->Write(); 
  c_pt->Write();
  mg_pt->Write(); 
  f_pt_l->Write();
  
  for(int hlt=0;hlt<nPaths;hlt++){
    string HLT = HLTs[hlt];
    if(HLT!="") HLT+="-";
    Jet_dEtadPhi_2b[hlt]->Write(); Jet_dEtadPhi_2b_clone[hlt]->Write(); Jet_dEtadPhi_2b[hlt]->Write(); Jet_pt1pt2_2b[hlt]->Write();
  }

  _fileOut->Write();
  _fileOut->Close();



}
