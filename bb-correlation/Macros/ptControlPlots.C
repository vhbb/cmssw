#include "vector.h"
#include "TFile.h"
#include <string>
#include <sstream>
#include <iostream>



void controlPlots(){
  gStyle->SetOptStat(0);

  std::vector<bool> *bvec = new std::vector<bool>; 

  bvec->push_back(0); 

  std::cout << bvec->size() << std::endl;

//   string dir="/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/BBC/CMSSW_3_7_0_patch2/src/crab/analyze/rootoutput/tempPTplots2/"; 
  string dir="/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/BBC/CMSSW_3_7_0_patch2/crab/analyze/rootoutput/tempPTplotsALL/"; 
  string dirMC="mc/"; 
  string dirDATA="data/";
  string strigger="jet30/";

  //LUMI DATA FOR JET30: 0.313154 JET50U: 3.06903 ((IN /PB))
  //MC LUMI IS 1/PB
  bool luminorm=0; double lumiDATA=0.313154;
  bool log=0; 

//   bool twoBonly = true; 

//   //pt hard reco (2V events)
//   string histoname = "Vert_ptrechard_2bVert"; 
//   string axis = "harder vertex pt"; 
//   int nbins = 100; double lower=0.0, upper = 200;
//   unsigned int rebinBy=2; 

  //pt soft reco (2V events)
  string histoname = "Vert_ptrecsoft_2bVert"; 
  string axis = "softer vertex pt"; 
  int nbins = 100; double lower=0.0, upper = 200;
  unsigned int rebinBy=2; 

//   //average pt reco (2V events)
//   string histoname = "Vert_ptrecaver_2bVert"; 
//   string axis = "average vertex pt"; 
//   int nbins = 100; double lower=0.0, upper = 200;
//   unsigned int rebinBy=2; 

//   //asymmetry pt reco (2V events)
//   string histoname = "Vert_ptrecasym_2bVert"; 
//   string axis = "asymmetry of vertex pt"; 
//   int nbins = 50; double lower=0.0, upper = 1.0;
//   unsigned int rebinBy=2; 


//   if(twoBonly) histoname = histoname + "_2b";


  TFile *fMC1 = TFile::Open((dir+dirMC+"histos_qcdpt15_test_v2.root").c_str()); 
  TFile *fMC2 = TFile::Open((dir+dirMC+"histos_qcdpt30_test_v2_125of127.root").c_str()); 
  TFile *fMC3 = TFile::Open((dir+dirMC+"histos_qcdpt80_test_v3.root").c_str()); 
  TFile *fMC4 = TFile::Open((dir+dirMC+"histos_qcdpt170_test_v1_79of80.root").c_str()); 

  TFile *fDATA = TFile::Open((dir+dirDATA+"histo_allDATAupTo144114.root").c_str()); 
  

  TH1F *hMC1; fMC1->GetObject((strigger+histoname).c_str(),hMC1);
  TH1F *hMC2; fMC2->GetObject((strigger+histoname).c_str(),hMC2);
  TH1F *hMC3; fMC3->GetObject((strigger+histoname).c_str(),hMC3); 
  TH1F *hMC4; fMC4->GetObject((strigger+histoname).c_str(),hMC4);

  TH1F *hMCcombined= new TH1F("hMCcombined",axis.c_str(),nbins,lower,upper);
  TH1F *hDATA; fDATA->GetObject((strigger+histoname).c_str(),hDATA);

  TCanvas *c = new TCanvas("quant",axis.c_str(),100,300,500,500);



  hMC1->SetFillColor(0);
  hMC2->SetFillColor(0);
  hMC3->SetFillColor(0);
  hMC4->SetFillColor(0);
  hMCcombined->SetFillColor(0);
  hDATA->SetFillColor(0);

  hMC1->Scale(143.866);
  hMC2->Scale(12.1072);
  hMC3->Scale(0.310862);
  hMC4->Scale(0.0082391);

  hMCcombined->Add(hMC1);
  hMCcombined->Add(hMC2);
  hMCcombined->Add(hMC3);
  hMCcombined->Add(hMC4);


  //norm  
  if(luminorm){
    hMCcombined->Scale(lumiDATA); 
  }
  else{
    double intMC=hMCcombined->Integral(); 
    double intDATA=hDATA->Integral(); 
    hMCcombined->Scale(intDATA/intMC); 
  }

  //numbers, cuts
  std::cout << "data lumi (/pb): " << lumiDATA << std::endl;
  std::cout << "trigger used: " << strigger << std::endl;
  std::cout << "Numbers\n";
  std::cout << "MC:\t" << hMCcombined->Integral() << std::endl;
  std::cout << "DATA:\t" << hDATA->Integral() << std::endl;

  //rebin
  hMCcombined->Rebin(rebinBy);
  hDATA->Rebin(rebinBy);

  //format
  if(log) c->SetLogy();
  double dmax = hMCcombined->GetMaximum(); 
  if(dmax<hDATA->GetMaximum()) dmax = hDATA->GetMaximum(); 
  if(log)  hMCcombined->GetYaxis()->SetRangeUser(0.01,dmax*10.0); 
  else hMCcombined->GetYaxis()->SetRangeUser(0.01,dmax+dmax/10.0); 
  
  hMCcombined->GetXaxis()->SetTitle(axis.c_str()); 
  hMCcombined->SetFillColor(8); 
  hMCcombined->SetLineColor(8); 
  hDATA->SetMarkerStyle(21); 
  
  //draw
  hMCcombined->Draw("hist");
  hDATA->Draw("samesPE");

  //leg
  TLegend *leg = new TLegend(0.65,0.75,0.9,0.9); 
  leg->SetFillColor(0); 
  leg->AddEntry(hDATA,"Data","PE"); 
  leg->AddEntry(hMCcombined,"Pythia","F");
  leg->Draw(); 


}
