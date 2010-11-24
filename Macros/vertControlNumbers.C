#include "vector.h"
#include "TFile.h"
#include <string>
#include <sstream>
#include <iostream>



void controlPlots(){
  gStyle->SetOptStat(0);

  string dir="/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/BBC/CMSSW_3_7_0_patch2/crab/analyze/rootfiles/";
  string dirMC="mcBins/withIP/";
  string dirDATA="dataFiles/withIP/";
  string cut="BEvents.jet30==1 && BEvents.ptHardestPJ>84";

//   //LUMI DATA FOR JET30: 0.313154 JET50U: 3.06903 ((IN /PB))
//   //MC LUMI IS 1/PB
//   bool luminorm=0; double lumiDATA=0.313154;
//   bool log=0; 


  string axis = "nvert"; 
  string toDraw = "nV";
  double lower=0, upper=10; unsigned int nbins=10; 

  TFile *fMC1 = TFile::Open((dir+dirMC+"qcdpt15_test_v3.root").c_str());
  TFile *fMC2 = TFile::Open((dir+dirMC+"qcdpt30_test_v3_125of127.root").c_str());
  TFile *fMC3 = TFile::Open((dir+dirMC+"qcdpt80_test_v4.root").c_str());
  TFile *fMC4 = TFile::Open((dir+dirMC+"qcdpt170_test_v2_78of80.root").c_str());

  TFile *fDATA = TFile::Open((dir+dirDATA+"allDATAupTo144114.root").c_str());



  TTree *tMC1 = (TTree*)fMC1->Get("bcanalyzer/tEvents");
  TTree *tMC2 = (TTree*)fMC2->Get("bcanalyzer/tEvents");
  TTree *tMC3 = (TTree*)fMC3->Get("bcanalyzer/tEvents");
  TTree *tMC4 = (TTree*)fMC4->Get("bcanalyzer/tEvents");
  TTree *tDATA = (TTree*)fDATA->Get("bcanalyzer/tEvents");

  TH1F *hMC1= new TH1F("hMC1",axis.c_str(),nbins,lower,upper);
  TH1F *hMC2= new TH1F("hMC2",axis.c_str(),nbins,lower,upper);
  TH1F *hMC3= new TH1F("hMC3",axis.c_str(),nbins,lower,upper);
  TH1F *hMC4= new TH1F("hMC4",axis.c_str(),nbins,lower,upper);
  TH1F *hMCcombined= new TH1F("hMCcombined",axis.c_str(),nbins,lower,upper);
  TH1F *hDATA= new TH1F("hDATA",axis.c_str(),nbins,lower,upper);


  tMC1->Draw(("BEvents." + toDraw + ">>hMC1").c_str(),(cut).c_str());
  std::cout << "done 1...";
  tMC2->Draw(("BEvents." + toDraw + ">>hMC2").c_str(),(cut).c_str());
  std::cout << "2...";
  tMC3->Draw(("BEvents." + toDraw + ">>hMC3").c_str(),(cut).c_str());
  std::cout << "3...";
  tMC4->Draw(("BEvents." + toDraw + ">>hMC4").c_str(),(cut).c_str());
  std::cout << "4...";
  tDATA->Draw(("BEvents." + toDraw + ">>hDATA").c_str(),(cut).c_str());
  std::cout << "all";



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


//   //norm  
//   if(luminorm){
//     hMCcombined->Scale(lumiDATA); 
//   }
//   else{
//     double intMC=hMCcombined->Integral(); 
//     double intDATA=hDATA->Integral(); 
//     hMCcombined->Scale(intDATA/intMC); 
//   }

//   //numbers, cuts
//   std::cout << "data lumi (/pb): " << lumiDATA << std::endl;
//   std::cout << "trigger used: " << strigger << std::endl;
//   std::cout << "Numbers\n";
//   std::cout << "MC:\t" << hMCcombined->Integral() << std::endl;
//   std::cout << "DATA:\t" << hDATA->Integral() << std::endl;

//   //rebin
//   hMCcombined->Rebin(rebinBy);
//   hDATA->Rebin(rebinBy);

//   //format
//   if(log) c->SetLogy();
//   double dmax = hMCcombined->GetMaximum(); 
//   if(dmax<hDATA->GetMaximum()) dmax = hDATA->GetMaximum(); 
//   if(log)  hMCcombined->GetYaxis()->SetRangeUser(0.01,dmax*10.0); 
//   else hMCcombined->GetYaxis()->SetRangeUser(0.01,dmax+dmax/10.0); 
  
//   hMCcombined->GetXaxis()->SetTitle(axis.c_str()); 
//   hMCcombined->SetFillColor(8); 
//   hMCcombined->SetLineColor(8); 
//   hDATA->SetMarkerStyle(21); 
  
  //draw
  hMCcombined->Draw("hist");
  hDATA->Draw("sames");


  //bin 0: overflow: bin i: i-1 vertex
  double rat32MC = (hMCcombined->GetBinContent(4))/(hMCcombined->GetBinContent(3)); 
  double rat32DATA = (hDATA->GetBinContent(4))/(hDATA->GetBinContent(3)); 

  double rat42MC = (hMCcombined->GetBinContent(5))/(hMCcombined->GetBinContent(3)); 
  double rat42DATA = (hDATA->GetBinContent(5))/(hDATA->GetBinContent(3)); 

  std::cout << "CUTS USED " << cut << std::endl;
  std::cout << "RATIO 3vert/2vert\n";
  std::cout << "MC \t" << rat32MC << std::endl;
  std::cout << "DATA \t" << rat32DATA << std::endl;

  std::cout << "RATIO 4vert/2vert\n";
  std::cout << "MC \t" << rat42MC << std::endl;
  std::cout << "DATA \t" << rat42DATA << std::endl;

//   //leg
//   TLegend *leg = new TLegend(0.65,0.75,0.9,0.9); 
//   leg->SetFillColor(0); 
//   leg->AddEntry(hDATA,"Data","PE"); 
//   leg->AddEntry(hMCcombined,"Pythia","F");
//   leg->Draw(); 


}
