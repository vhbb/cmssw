#include "vector.h"
#include "TFile.h"
#include <string>
#include <sstream>
#include <iostream>



void controlPlots(){
  gStyle->SetOptStat(0);

//   string dir="/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/BBC/CMSSW_3_7_0_patch2/src/crab/analyze/rootfiles/"; 
//   string dirMC="mcBins/withIP/"; 
//   string dirDATA="dataFiles/withIP/";
//   string cut="BCands.selected==1 && BCands.jet30==1 && BCands.ptHardestPJ>84";

//   string cut="BCands.selected==1 && BCands.ptHardestPJ>84";

  string dirMC="/scratch/leo/";
  string dirDATA="/scratch/wehrlilu/BBCORR/MACROS/TREES/DATA/"; 

  string cut="BCands.selected==1 && BCands.jet30==1 && BCands.ptHardestPJ>84";


  //LUMI DATA FOR JET30: 0.313154 JET50U: 3.06903 ((IN /PB))
  //MC LUMI IS 1/PB
  bool luminorm=0; double lumiDATA=0.313154;
  bool log=1; 

//   //mass
//   string toDraw="massBcand";
//   string axis="vertex mass";
//   int nbins=40; double lower=0.0, upper=8.0; 

//   //dist3D
//   string toDraw="dist3D";
//   string axis="flight distance (3D)";
//   int nbins=50; double lower=0.0, upper=10.0; 

  //distSig3D
  string toDraw="distSig3D";
  string axis="flight distance significance (3D)";
  int nbins=40; double lower=0.0, upper=300.0; 

  TFile *fMC1 = TFile::Open((dirMC+"anV3-QCD_Pt15_Spring10-V8b.root").c_str()); 
  TFile *fMC2 = TFile::Open((dirMC+"anV3-QCD_Pt30_Spring10-V8b.root").c_str()); 
  TFile *fMC3 = TFile::Open((dirMC+"anV3-QCD_Pt80_Spring10-V8b.root").c_str()); 
  TFile *fMC4 = TFile::Open((dirMC+"anV3-QCD_Pt170_Spring10-V8b.root").c_str()); 


  TFile *fDATA = TFile::Open((dirDATA+"mergedDATA.root").c_str()); 

  TTree *tMC1 = (TTree*)fMC1->Get("bcanalyzer/tBCands");
  TTree *tMC2 = (TTree*)fMC2->Get("bcanalyzer/tBCands");
  TTree *tMC3 = (TTree*)fMC3->Get("bcanalyzer/tBCands");
  TTree *tMC4 = (TTree*)fMC4->Get("bcanalyzer/tBCands");
  TTree *tDATA = (TTree*)fDATA->Get("bcanalyzer/tBCands");

  TH1F *hMC1= new TH1F("hMC1",axis.c_str(),nbins,lower,upper);
  TH1F *hMC2= new TH1F("hMC2",axis.c_str(),nbins,lower,upper);
  TH1F *hMC3= new TH1F("hMC3",axis.c_str(),nbins,lower,upper);
  TH1F *hMC4= new TH1F("hMC4",axis.c_str(),nbins,lower,upper);
  TH1F *hMCcombined= new TH1F("hMCcombined",axis.c_str(),nbins,lower,upper);
  TH1F *hDATA= new TH1F("hDATA",axis.c_str(),nbins,lower,upper);

  TCanvas *c = new TCanvas("quant",axis.c_str(),100,300,500,500);

  tMC1->Draw(("BCands." + toDraw + ">>hMC1").c_str(),(cut).c_str());
  tMC2->Draw(("BCands." + toDraw + ">>hMC2").c_str(),(cut).c_str());
  tMC3->Draw(("BCands." + toDraw + ">>hMC3").c_str(),(cut).c_str());
  tMC4->Draw(("BCands." + toDraw + ">>hMC4").c_str(),(cut).c_str());
  tDATA->Draw(("BCands." + toDraw + ">>hDATA").c_str(),(cut).c_str());


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

  std::cout << "numbers test tot " << hMCcombined->Integral() << std::endl;

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
  std::cout << "\ncut used: " << cut << std::endl;
  std::cout << "data lumi (/pb): " << lumiDATA << std::endl;
  std::cout << "Numbers\n";
  std::cout << "MC:\t" << hMCcombined->Integral() << std::endl;
  std::cout << "DATA:\t" << hDATA->Integral() << std::endl;

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
  hMCcombined->Draw();
  hDATA->Draw("samesPE");

  //leg
  TLegend *leg = new TLegend(0.65,0.75,0.9,0.9); 
  leg->SetFillColor(0); 
  leg->AddEntry(hDATA,"Data","PE"); 
  leg->AddEntry(hMCcombined,"Pythia","F");
  leg->Draw(); 


}
