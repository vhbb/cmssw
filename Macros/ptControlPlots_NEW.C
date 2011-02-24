#include "vector.h"
#include "TFile.h"
#include <string>
#include <sstream>
#include <iostream>



void controlPlots(string quant="ptsoft", string trigger="30"){
  gStyle->SetOptStat(0);
  gROOT->ProcessLine(".L ~/tdrstyle.C");
  setTDRStyle(); 

  //style modifications:
  gStyle->SetErrorX(0.5);
  gStyle->SetPadTopMargin(0.099);//.005 --> 0.08
  gStyle->SetPadBottomMargin(0.105);//0.13 --> 0.1
  gStyle->SetMarkerStyle(1); 

  string ljcut = "84";
  if(trigger=="50") ljcut="120";
  if(trigger=="15") ljcut="56";

//   string dir="/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/BBC/CMSSW_3_7_0_patch2/crab/analyze/rootoutput/tempPTplotsALL/"; 

  string dir="/scratch/leo/finalOpenBins36X_No6U/";
  string mcfile="histo-Pythia6-HLT_Jet" + trigger + "U-BBCorr_HadrJetMatchPt_OpenHLT-total.root"; 
  string datafile="histo-Data-HLT_Jet" + trigger + "U-BBCorr_HadrJetMatchPt_OpenHLT-total.root"; 

  //LUMI DATA FOR JET30: 0.313154 JET50U: 3.06903 ((IN /PB))
  //MC LUMI IS 1/PB
//   bool luminorm=0; double lumiDATA=0.313154;
  bool log=0; 


  string histoname; 
  string axis; 
  int nbins; 
  unsigned int rebinBy; 

  if(quant=="pthard"){
    //pt hard reco (2V events)
    histoname = "Vert_ptrechard_2b"; 
    axis = "p_{T}^{rec} harder B candidate (GeV)"; 
    nbins = 100; double lower=0.0, upper = 140;
    rebinBy=2; 
  }

  if(quant=="ptsoft"){
  //pt soft reco (2V events)
    histoname = "Vert_ptrecsoft_2b"; 
    axis = "p_{T}^{rec} softer B candidate (GeV)"; 
    nbins = 100; double lower=0.0, upper = 110;
    rebinBy=2; //2
  }

  if(quant=="ptaver"){
    //average pt reco (2V events)
    histoname = "Vert_ptrecaver_2b"; 
    axis = "average B candidate p_{T}^{rec} (GeV)"; 
    nbins = 100; double lower=0.0, upper = 120;
    rebinBy=2; 
  }

  if(quant=="ptasym"){
    //asymmetry pt reco (2V events)
    histoname = "Vert_ptrecasym_2b"; 
    axis = "asymmetry of B candidate p_{T}^{rec}"; 
    nbins = 50; double lower=0.0, upper = 1.0;
    rebinBy=2; 
  }

//   if(twoBonly) histoname = histoname + "_2b";


  TFile *fMC = TFile::Open((dir+mcfile).c_str()); 

  TFile *fDATA = TFile::Open((dir+datafile).c_str()); 
  

  TH1F *hMC; fMC->GetObject((histoname).c_str(),hMC);
  TH1F *hDATA; fDATA->GetObject((histoname).c_str(),hDATA);

  TCanvas *c = new TCanvas("quant",axis.c_str(),100,300,500,500);


  std::cout << "BINS " << hMC->GetNbinsX() << "\n";

  hMC->SetFillColor(0);
  hDATA->SetFillColor(0);

  //norm to data area
  double intMC=hMC->Integral(); 
  double intDATA=hDATA->Integral(); 
  hMC->Scale(intDATA/intMC); 
  

//   //numbers, cuts
//   std::cout << "data lumi (/pb): " << lumiDATA << std::endl;
//   std::cout << "trigger used: " << strigger << std::endl;
//   std::cout << "Numbers\n";

  std::cout << "MC:\t" << intMC << std::endl;
  std::cout << "DATA:\t" << intDATA << std::endl;

  //rebin
  hMC->Rebin(rebinBy);
  hDATA->Rebin(rebinBy);

  //format
  if(log) c->SetLogy();
  double dmax = hMC->GetMaximum(); 
  if(dmax<hDATA->GetMaximum()) dmax = hDATA->GetMaximum();
   //manual correction
  if(quant=="ptasym") dmax = 70; 
  if(log)  hMC->GetYaxis()->SetRangeUser(0.01,dmax*10.0); 
  else hMC->GetYaxis()->SetRangeUser(0.01,dmax+dmax/10.0); 
  

  hMC->GetXaxis()->SetTitle(axis.c_str()); 
  hMC->GetYaxis()->SetTitle("number of events"); 
  hMC->GetYaxis()->SetTitleOffset(1.5); 
  hMC->GetXaxis()->SetTitleOffset(1.25); 
  hMC->GetXaxis()->SetRangeUser(lower, upper); 

  hMC->SetFillColor(8); 
  hMC->SetLineColor(8); 
  hDATA->SetMarkerStyle(21); 
  
  //draw
  hMC->Draw("hist");
  hDATA->Draw("samesPE1");

  //leg
  TLegend *leg = new TLegend(0.5,0.72,0.89,0.86); 
  leg->SetHeader(("leading jet p_{T} > " + ljcut + " GeV").c_str()); 
  leg->SetFillColor(0); 
  leg->AddEntry(hDATA,"Data","PE1"); 
  leg->AddEntry(hMC,"PYTHIA","F");
  leg->SetBorderSize(0); 
  leg->Draw(); 

  //cms label 
  TLatex *   tex = new TLatex(0,dmax+dmax/6.3,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}"); //was 50 //was /8
  tex->SetTextSize(0.040); //0.044
  tex->SetLineWidth(2);
  tex->Draw();

//   //leading jet label
//   TLatex *   tex2 = new TLatex(0,dmax+dmax/8.0,""); //was 50
//   tex2->SetTextSize(0.040); //0.044
//   tex2->SetLineWidth(2);
//   tex2->Draw();

  c->RedrawAxis();
  c->SetTickx(); 
  c->SetTicky();
  //setTDRStyle(); 

}
