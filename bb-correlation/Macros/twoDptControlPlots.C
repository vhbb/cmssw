#include "vector.h"
#include "TFile.h"
#include <string>
#include <sstream>
#include <iostream>



void controlPlots(){
  gStyle->SetOptStat(0);

//   string dir="/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/BBC/CMSSW_3_7_0_patch2/src/crab/analyze/rootoutput/tempPTplots3/"; 
  string dir="/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/BBC/CMSSW_3_7_0_patch2/crab/analyze/rootoutput/tempPTplotsALL/"; 
  string dirMC="mc/"; 
  string dirDATA="data/";
  string strigger="jet30/";

  //LUMI DATA FOR JET30: 0.313154 JET50U: 3.06903 ((IN /PB))
  //MC LUMI IS 1/PB
  bool luminorm=0; double lumiDATA=0.313154;
  bool log=1; 

  string xaxislabel ="#Delta R";

  
//   //pt hard reco (2V events)
//   string histoname = "Vert_ptrechardVSdR_2bVert"; 
//   string axis = "harder vertex pt vs dR"; 
//   int nbins = 60; double lower=0.0, upper = 6.0;
//   int nbins2 = 100; double lower2=0.0, upper2 = 200;
//   unsigned int rebinBy=2; 
//   string filename="harderptvsdr.root";

//   //pt soft reco (2V events)
//   string histoname = "Vert_ptrecsoftVSdR_2bVert"; 
//   string axis = "softer vertex pt vs dR"; 
//   int nbins = 60; double lower=0.0, upper = 6.0;
//   int nbins2 = 100; double lower2=0.0, upper2 = 200;
//   unsigned int rebinBy=2; 
//   string filename="softerptvsdr.root";

//   //average pt reco (2V events)
//   string histoname = "Vert_ptrecaverVSdR_2bVert"; 
//   string axis = "average vertex pt vs dR"; 
//   int nbins = 60; double lower=0.0, upper = 6.0;
//   int nbins2 = 100; double lower2=0.0, upper2 = 200;
//   unsigned int rebinBy=2; 
//   string filename="averageptvsdr.root";

//   //asymmetry pt reco (2V events)
//   string histoname = "Vert_ptrecasymVSdR_2bVert"; 
//   string axis = "asymmetry of vertex pt vs dR"; 
//   int nbins = 60; double lower=0.0, upper = 6.0;
//   int nbins2 = 50; double lower2=0.0, upper2 = 1.0;
//   unsigned int rebinBy=2; 
//   string filename="asymmetryptvsdr.root";

//PTRECO/PTSIM PLOTS
  //asymmetry pt reco (2V events)
  string histoname = "Vert_ptrecasymVSdR_2bVert"; 
  string axis = "asymmetry of vertex pt vs dR"; 
  int nbins = 60; double lower=0.0, upper = 6.0;
  int nbins2 = 50; double lower2=0.0, upper2 = 1.0;
  unsigned int rebinBy=2; 
  string filename="asymmetryptvsdr.root";





//   if(twoBonly) histoname = histoname + "_2b";

  TFile *fout = TFile::Open((filename).c_str(),"RECREATE"); 

  TFile *fMC1 = TFile::Open((dir+dirMC+"histos_qcdpt15_test_v2.root").c_str()); 
  TFile *fMC2 = TFile::Open((dir+dirMC+"histos_qcdpt30_test_v2_125of127.root").c_str()); 
  TFile *fMC3 = TFile::Open((dir+dirMC+"histos_qcdpt80_test_v3.root").c_str()); 
  TFile *fMC4 = TFile::Open((dir+dirMC+"histos_qcdpt170_test_v1_79of80.root").c_str()); 

  TFile *fDATA = TFile::Open((dir+dirDATA+"histo_allDATAupTo144114.root").c_str()); 
  

  TH2F *hMC1; fMC1->GetObject((strigger+histoname).c_str(),hMC1);
  TH2F *hMC2; fMC2->GetObject((strigger+histoname).c_str(),hMC2);
  TH2F *hMC3; fMC3->GetObject((strigger+histoname).c_str(),hMC3); 
  TH2F *hMC4; fMC4->GetObject((strigger+histoname).c_str(),hMC4);

  TH2F *hMCcombined= new TH2F("hMCcombined",axis.c_str(),nbins,lower,upper,nbins2,lower2,upper2);
  TH2F *hDATA; fDATA->GetObject((strigger+histoname).c_str(),hDATA);
  hDATA->SetName("hDATA");



  hMC1->Scale(143.866);
  hMC2->Scale(12.1072);
  hMC3->Scale(0.310862);
  hMC4->Scale(0.0082391);

  hMCcombined->Add(hMC1);
  hMCcombined->Add(hMC2);
  hMCcombined->Add(hMC3);
  hMCcombined->Add(hMC4);



//   //numbers, cuts
//   std::cout << "data lumi (/pb): " << lumiDATA << std::endl;
//   std::cout << "trigger used: " << strigger << std::endl;
//   std::cout << "Numbers\n";
//   std::cout << "MC:\t" << hMCcombined->Integral() << std::endl;
//   std::cout << "DATA:\t" << hDATA->Integral() << std::endl;

  TCanvas *c = new TCanvas("quant",axis.c_str(),100,300,500,500);

  
  hMCcombined->GetXaxis()->SetTitle(xaxislabel.c_str()); 
  hDATA->GetXaxis()->SetTitle(xaxislabel.c_str()); 

  //draw
  hMCcombined->Draw("colz");
  if(log) c->SetLogz();

  TCanvas *c2 = new TCanvas("quant2",axis.c_str(),300,300,500,500);
  hDATA->Draw("colz");
  if(log) c2->SetLogz();

  TCanvas *c3 = new TCanvas("quant3",axis.c_str(),500,300,500,500);
  //profiles
//   hMCcombined->ProfileX();
  hDATA->ProfileX("_pfx");
  hMCcombined->ProfileX();

  //rebin
  hMCcombined_pfx->Rebin(rebinBy);
  hDATA_pfx->Rebin(rebinBy);


  //format
  double dmax = hMCcombined_pfx->GetMaximum(); 
  if(dmax<hDATA_pfx->GetMaximum()) dmax = hDATA_pfx->GetMaximum(); 
  hMCcombined_pfx->GetYaxis()->SetRangeUser(0.01,dmax+dmax/10.0); 
  hDATA_pfx->SetLineColor(2); 

  hMCcombined_pfx->Draw();
  hDATA_pfx->Draw("sames"); 


  //leg
  TLegend *leg = new TLegend(0.65,0.75,0.9,0.9); 
  leg->SetFillColor(0); 
  leg->AddEntry(hDATA_pfx,"Data","L"); 
  leg->AddEntry(hMCcombined_pfx,"Pythia","L");
  leg->Draw(); 

  fout->cd();
  c->Write(); 
  c2->Write();
  c3->Write();
  fout->Close();

}
