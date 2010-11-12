{
  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();

  const int nPaths = 4; string HLTs[nPaths] = {"HLT_L1Jet6U","HLT_Jet15U","HLT_Jet30U","HLT_Jet50U"}

  TFile*files[nPaths];
  TH1F *Jet_dPhi_2b[nPaths],*Jet_dPhi[nPaths];
  //  TMultiGraph *mg_dPhi = new TMultiGraph(); TGraphAsymmErrors *pur_Jet_dPhi[nPaths];
  TH1F *Jet_dR_2b[nPaths],*Jet_dR[nPaths];
  //TMultiGraph *mg_dR = new TMultiGraph(); TGraphAsymmErrors *pur_Jet_dR[nPaths];
  TH1F *Jet_pt_2b[nPaths],*Jet_pt[nPaths];
  //TMultiGraph *mg_pt = new TMultiGraph(); TGraphAsymmErrors *pur_Jet_pt[nPaths];
  TH1F *Jet_dEta_2b[nPaths],*Jet_dEta[nPaths];
  //TMultiGraph *mg_dEta = new TMultiGraph(); TGraphAsymmErrors *pur_Jet_dEta[nPaths];

  TCanvas *c_dEtadPhi[nPaths];
  TH2F *Jet_dEtadPhi_2b[nPaths],*Jet_dEtadPhi[nPaths];
  TCanvas *c_dEtadPhi_clone[nPaths];
  TH2F *Jet_dEtadPhi_2b_clone[nPaths],*Jet_dEtadPhi_clone[nPaths];

  TCanvas *c_pt1pt2[nPaths];
  TH2F *Jet_pt1pt2_2b[nPaths],*Jet_pt1pt2[nPaths];

  TCanvas *c_BHadr_dR[nPaths];
  TH1F *h_BHadr_dR[nPaths], *h_BHadrEff_dR[nPaths];
  TH1F* den_BHadr_dR[nPaths];

  TLegend *leg = new TLegend(0.1,0.2,0.3,0.4);


  int hlt=0;
  string HLT = "";
  //files[hlt] = TFile::Open(("histo-Madgraph-"+HLT+"total.root").c_str());
  files[hlt] = TFile::Open("histo-Pythia6-total.root");
  
  files[hlt]->GetObject("Jet_dPhi_2b",Jet_dPhi[hlt]); Jet_dPhi[hlt]->Rebin(2);
  files[hlt]->GetObject("Jet_dPhi_2b_2b",Jet_dPhi_2b[hlt]); Jet_dPhi_2b[hlt]->Rebin(2);

  //pur_Jet_dPhi[hlt] = new TGraphAsymmErrors(Jet_dPhi_2b[hlt],Jet_dPhi[hlt]);
  Jet_dPhi_2b[hlt]->Divide(Jet_dPhi_2b[hlt],Jet_dPhi[hlt],1,1,"B");
  // pur_Jet_dPhi[hlt]->SetLineColor(4-hlt);     pur_Jet_dPhi[hlt]->SetMarkerColor(4-hlt);
  //mg_dPhi->Add(pur_Jet_dPhi[hlt]);
  string label = HLT;
  if(HLT=="") label = "Total";
  leg->AddEntry(Jet_dPhi_2b[hlt],label.c_str(),"pl");


  files[hlt]->GetObject("Jet_dR_2b",Jet_dR[hlt]); Jet_dR[hlt]->Rebin(2);
  files[hlt]->GetObject("Jet_dR_2b_2b",Jet_dR_2b[hlt]); Jet_dR_2b[hlt]->Rebin(2);
  //pur_Jet_dR[hlt] = new TGraphAsymmErrors(Jet_dR_2b[hlt],Jet_dR[hlt]);
  //pur_Jet_dR[hlt]->SetLineColor(4-hlt);     pur_Jet_dR[hlt]->SetMarkerColor(4-hlt);
  //mg_dR->Add(pur_Jet_dR[hlt]);
  Jet_dR_2b[hlt]->Divide(Jet_dR_2b[hlt],Jet_dR[hlt],1,1,"B");
  Jet_dR_2b[hlt]->SetName("Jet_BPurity_dR");

  files[hlt]->GetObject("Jet_dEta_2b",Jet_dEta[hlt]); Jet_dEta[hlt]->Rebin(2);
  files[hlt]->GetObject("Jet_dEta_2b_2b",Jet_dEta_2b[hlt]); Jet_dEta_2b[hlt]->Rebin(2);
  //pur_Jet_dEta[hlt] = new TGraphAsymmErrors(Jet_dEta_2b[hlt],Jet_dEta[hlt]);
  //pur_Jet_dEta[hlt]->SetLineColor(4-hlt);     pur_Jet_dEta[hlt]->SetMarkerColor(4-hlt);
  //mg_dEta->Add(pur_Jet_dEta[hlt]);
  Jet_dEta_2b[hlt]->Divide(Jet_dEta_2b[hlt],Jet_dEta[hlt],1,1,"B");
  Jet_dEta_2b[hlt]->SetName("Jet_BPurity_dEta");

  /*
  files[hlt]->GetObject("Jet_dEtadPhi_2b",Jet_dEtadPhi[hlt]); Jet_dEtadPhi[hlt]->Rebin2D(2,2);
  files[hlt]->GetObject("Jet_dEtadPhi_2b_2b",Jet_dEtadPhi_2b[hlt]); Jet_dEtadPhi_2b[hlt]->Rebin2D(2,2);
  Jet_dEtadPhi_2b_clone[hlt] = (TH2F*)Jet_dEtadPhi_2b[hlt]->Clone();
  Jet_dEtadPhi_2b[hlt]->Divide(Jet_dEtadPhi[hlt]);
  files[hlt]->GetObject("Jet_pt1pt2_2b",Jet_pt1pt2[hlt]); //Jet_pt1pt2[hlt]->Rebin2D(2,2);
  files[hlt]->GetObject("Jet_pt1pt2_2b_2b",Jet_pt1pt2_2b[hlt]); //Jet_pt1pt2_2b[hlt]->Rebin2D(2,2);
  Jet_pt1pt2_2b[hlt]->Divide(Jet_pt1pt2[hlt]);
  */
  
  files[hlt]->GetObject("Jet_pt_2b",Jet_pt[hlt]); Jet_pt[hlt]->Rebin(1);
  files[hlt]->GetObject("Jet_pt_2b_2b",Jet_pt_2b[hlt]); Jet_pt_2b[hlt]->Rebin(1);
  //pur_Jet_pt[hlt] = new TGraphAsymmErrors(Jet_pt_2b[hlt],Jet_pt[hlt]);
  //pur_Jet_pt[hlt]->SetLineColor(4-hlt);     pur_Jet_pt[hlt]->SetMarkerColor(4-hlt);
  //mg_pt->Add(pur_Jet_pt[hlt]);
  Jet_pt_2b[hlt]->Divide(Jet_pt_2b[hlt],Jet_pt[hlt],1,1,"B");
  Jet_pt_2b[hlt]->SetName("Jet_BPurity_pt");


  for(int hlt=0;hlt<nPaths;hlt++){
    string HLT = HLTs[hlt];
    
    TFile *_file0 = TFile::Open(("Pythia6BB/histo-"+HLT+"-BBCorr_scaledEff-anV3-InclusiveBB_Pt30_Spring10-V8b.root").c_str());
    
    c_BHadr_dR[hlt] = new TCanvas( ("BHadrEff_dR_"+HLT).c_str(), ("BHadrEff_dR_"+HLT).c_str());
    _file0->GetObject("Jet_BHadrJetMatched_dR",h_BHadr_dR[hlt]);  _file0->GetObject("BHadr_dR",den_BHadr_dR[hlt]);
    h_BHadrEff_dR[hlt] = (TH1F*)h_BHadr_dR[hlt]->Clone();
    h_BHadr_dR[hlt]->Rebin(2); h_BHadrEff_dR[hlt]->Rebin(2); den_BHadr_dR[hlt]->Rebin(2);
    h_BHadrEff_dR[hlt]->Divide(h_BHadr_dR[hlt],den_BHadr_dR[hlt],1,1,"B"); 
    h_BHadrEff_dR[hlt]->SetName(("BHadrEff_dR_"+HLT).c_str());
    h_BHadr_dR[hlt]->SetName(("BHadrJetMatched_dR_"+HLT).c_str());
    den_BHadr_dR[hlt]->SetName(("BHadr_dR_"+HLT).c_str());

  }

  //HLT = "MadgraphVsPythia";
  /*
  TCanvas *c_dPhi = new TCanvas( ("BTagPurVsDphi"+HLT).c_str(), ("BTagPurVsDphi"+HLT).c_str());
  mg_dPhi->Draw("ap");
  temp2->Draw("sames");
  leg->Draw();

  TCanvas *c_dEta = new TCanvas( ("BTagPurVsDeta"+HLT).c_str(), ("BTagPurVsDeta"+HLT).c_str());
  mg_dEta->Draw("ap");
  leg->Draw();

  TCanvas *c_dR = new TCanvas( ("BTagPurVsDR"+HLT).c_str(), ("BTagPurVsDR"+HLT).c_str());
  mg_dR->Draw("ap");
  leg->Draw();

  TCanvas *c_pt = new TCanvas( ("BTagPurVsPt"+HLT).c_str(), ("BTagPurVsPt"+HLT).c_str());
  mg_pt->Draw("ap");
  leg->Draw();
  pur_Jet_pt[0]->Fit("pol1","F","sames",30,250);
  TF1 *f_pt_l = pur_Jet_pt[0]->GetFunction("pol1");
  f_pt_l->SetName("f_purity_pt");
  f_pt_l->SetLineColor(kBlack);f_pt_l->SetLineWidth(2);
  */
  
  TFile *_fileOut = TFile::Open("histo-Pythia6-flavourComp.root","RECREATE");
  //TFile *_fileOut = TFile::Open("histo-Pythia6vsMadgraph-flavourComp.root","RECREATE");
  /*c_dPhi->Write(); c_dEta->Write(); c_dR->Write(); 
  c_pt->Write();
  mg_pt->Write(); 
  f_pt_l->Write();
  */
  Jet_dPhi_2b[0]->Write();
  Jet_dR_2b[0]->Write();
  Jet_dEta_2b[0]->Write();
  Jet_pt_2b[0]->Write();

  for(int hlt=0;hlt<nPaths;hlt++){
    string HLT = HLTs[hlt];
    //c_BHadr_dR[hlt]->Write();
    den_BHadr_dR[hlt]->Write();
    h_BHadr_dR[hlt]->Write();
    h_BHadrEff_dR[hlt]->Write();

    //    if(HLT!="") HLT+="-";
    //Jet_dEtadPhi_2b[hlt]->Write(); Jet_dEtadPhi_2b_clone[hlt]->Write(); Jet_dEtadPhi_2b[hlt]->Write(); Jet_pt1pt2_2b[hlt]->Write();
  }

  _fileOut->Write();
  _fileOut->Close();



}
