

void combine(bool write=false, bool drawContrib=false){
  
  string dir="mcatnloJET30/"; 
//   string dir="mcatnloJET15/"; 
//   string dir="mcatnloJET50/"; 
  double lumiToNorm=1.0; //pb
    
    // 	  TFile *fbin1 = TFile::Open("bbdist_mcatnlo_pthat0_20M_pthat.root");
    // 	  TFile *fbin2 = TFile::Open("bbdist_mcatnlo_pthat30_pthat.root");
    // 	  TFile *fbin3 = TFile::Open("bbdist_mcatnlo_pthat80_20M_pthat.root");
  
  //jet30
   TFile *fbin1 = TFile::Open((dir+"bbdist_mcatnlo_pthat0_lj84_eta2.root").c_str());
   TFile *fbin2 = TFile::Open((dir+"bbdist_mcatnlo_pthat30_lj84_eta2.root").c_str());
   TFile *fbin3 = TFile::Open((dir+"bbdist_mcatnlo_pthat80_lj84_eta2.root").c_str());
  //jet15
//   TFile *fbin1 = TFile::Open((dir+"bbdist_mcatnlo_pthat0_lj56_eta2.root").c_str());
//   TFile *fbin2 = TFile::Open((dir+"bbdist_mcatnlo_pthat30_lj56_eta2.root").c_str());
//   TFile *fbin3 = TFile::Open((dir+"bbdist_mcatnlo_pthat80_lj56_eta2.root").c_str());
  //jet50
//   TFile *fbin1 = TFile::Open((dir+"bbdist_mcatnlo_pthat0_lj120_eta2.root").c_str());
//   TFile *fbin2 = TFile::Open((dir+"bbdist_mcatnlo_pthat30_lj120_eta2.root").c_str());
//   TFile *fbin3 = TFile::Open((dir+"bbdist_mcatnlo_pthat80_lj120_eta2.root").c_str());
  
  if(write) TFile *fout = TFile::Open((dir+"mcatnlo_combined.root").c_str(),"UPDATE");

  TH1F *hdrbin1; fbin1->GetObject("bbdistmcatnlo/hdR_mcatnlo",hdrbin1);
  TH1F *hdrbin2; fbin2->GetObject("bbdistmcatnlo/hdR_mcatnlo",hdrbin2);
  TH1F *hdrbin3; fbin3->GetObject("bbdistmcatnlo/hdR_mcatnlo",hdrbin3);
  TH1F *hdphibin1; fbin1->GetObject("bbdistmcatnlo/hdPhi_mcatnlo",hdphibin1);
  TH1F *hdphibin2; fbin2->GetObject("bbdistmcatnlo/hdPhi_mcatnlo",hdphibin2);
  TH1F *hdphibin3; fbin3->GetObject("bbdistmcatnlo/hdPhi_mcatnlo",hdphibin3);
  TH1F *hdetabin1; fbin1->GetObject("bbdistmcatnlo/hdEta_mcatnlo",hdetabin1);
  TH1F *hdetabin2; fbin2->GetObject("bbdistmcatnlo/hdEta_mcatnlo",hdetabin2);
  TH1F *hdetabin3; fbin3->GetObject("bbdistmcatnlo/hdEta_mcatnlo",hdetabin3);
  
  TH1F *hpthatbin1; fbin1->GetObject("bbdistmcatnlo/hpthat",hpthatbin1);
  TH1F *hpthatbin2; fbin2->GetObject("bbdistmcatnlo/hpthat",hpthatbin2);
  TH1F *hpthatbin3; fbin3->GetObject("bbdistmcatnlo/hpthat",hpthatbin3);
  
  TH1F *hdR_mcatnlo = new TH1F("hdR_mcatnlo","#Delta R between two B hadrons",60,0,6);
  TH1F *hdPhi_mcatnlo = new TH1F("hdPhi_mcatnlo","#Delta #phi between two B hadrons",32,0,3.2);
  TH1F *hdEta_mcatnlo = new TH1F("hdEta_mcatnlo","#Delta #eta between two B hadrons",100,0,10);
  TH1F *hpthat = new TH1F("hpthat","pthat after cut",100,0,400);
  hdR_mcatnlo->GetXaxis()->SetTitle("#Delta R");
  hdPhi_mcatnlo->GetXaxis()->SetTitle("#Delta #phi");
  hdEta_mcatnlo->GetXaxis()->SetTitle("#Delta #eta");
  hpthat->GetXaxis()->SetTitle("p_{t} of hardest b quark");
  
  //double lf = 0.188921/0.04762; //pb/pb
//   double lf = 0.296705/0.04762; //pb/pb
  double lf = lumiToNorm/0.04762; //pb/pb
  hdrbin1->Scale(lf*20.);
  hdrbin2->Scale(lf);
  hdrbin3->Scale(lf/20.);
  hdphibin1->Scale(lf*20.);
  hdphibin2->Scale(lf);
  hdphibin3->Scale(lf/20.);
  hdetabin1->Scale(lf*20.);
  hdetabin2->Scale(lf);
  hdetabin3->Scale(lf/20.);
  hpthatbin1->Scale(lf*20.);
  hpthatbin2->Scale(lf);
  hpthatbin3->Scale(lf/20.);
  
  
  hdR_mcatnlo->Add(hdrbin1);
  hdR_mcatnlo->Add(hdrbin2);
  hdR_mcatnlo->Add(hdrbin3);
  hdPhi_mcatnlo->Add(hdphibin1);
  hdPhi_mcatnlo->Add(hdphibin2);
  hdPhi_mcatnlo->Add(hdphibin3);
  hdEta_mcatnlo->Add(hdetabin1);
  hdEta_mcatnlo->Add(hdetabin2);
  hdEta_mcatnlo->Add(hdetabin3);
  hpthat->Add(hpthatbin1);
  hpthat->Add(hpthatbin2);
  hpthat->Add(hpthatbin3);
  

  std::cout << "Factor 030 " << lf*20.0 << " " << hdrbin1->Integral() << std::endl;
  std::cout << "Factor 3080 " << lf << " " << hdrbin2->Integral() << std::endl;
  std::cout << "Factor 80inf " << lf/20.0 << " " << hdrbin3->Integral() << std::endl;
  
  
  
  TCanvas *cmcatnloDR = new TCanvas("cmcatnloDR","#Delta R",100,100,600,400);
  hdR_mcatnlo->Draw();
  if(drawContrib){
    hdrbin1->SetLineColor(2);
    hdrbin2->SetLineColor(4);
    hdrbin3->SetLineColor(8);
    hdrbin1->Draw("same");
    hdrbin2->Draw("same");
    hdrbin3->Draw("same");	  
  }
  
  TCanvas *cmcatnloDPHI = new TCanvas("cmcatnloDPHI","#Delta #phi",400,100,600,400);
  hdPhi_mcatnlo->Draw();
  if(drawContrib){
    hdphibin1->SetLineColor(2);
    hdphibin2->SetLineColor(4);
    hdphibin3->SetLineColor(8);
    hdphibin1->Draw("same");
    hdphibin2->Draw("same");
    hdphibin3->Draw("same");	  
  }
  
  TCanvas *cmcatnloDETA = new TCanvas("cmcatnloDETA","#Delta #eta",100,400,600,400);
  hdEta_mcatnlo->Draw();
  if(drawContrib){
    hdetabin1->SetLineColor(2);
    hdetabin2->SetLineColor(4);
    hdetabin3->SetLineColor(8);
    hdetabin1->Draw("same");
    hdetabin2->Draw("same");
    hdetabin3->Draw("same");	  
  }
  
  TCanvas *cmcatnloPTH = new TCanvas("cmcatnloPTH","PThat",400,400,600,400);
  hpthat->Draw();
  
  if(write){
    fout->Write();
    fout->Close();
  }
}
