#include "vector.h"
#include "TFile.h"
#include <string>
#include <sstream>
#include <iostream>



void controlPlots(string todr="ntracks"){
  gStyle->SetOptStat(0);
  bool haveFracInLeg = 0; 

  string dir="/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/BBC/CMSSW_3_7_0_patch2/crab/analyze/rootfiles/"; 
  string dirMC="withNTR/mc/"; 
  string dirDATA="withNTR/data/";
  string cut="BCands.selected==1 && BCands.jet30==1 && BCands.ptHardestPJ>84";

  //same:
  //  string cut="BCands.massBcand>1.4 && BCands.distSig3D>5.0 && abs(BCands.eta)<2.0 && BCands.jet30==1 && BCands.ptHardestPJ>84";

//   string cut="BCands.selected==1 && BCands.ptHardestPJ>84";

  //LUMI DATA FOR JET30: 0.313154 JET50U: 3.06903 ((IN /PB))
  //MC LUMI IS 1/PB
  bool luminorm=0; double lumiDATA=0.313154;
  bool log=1; 


  string toDraw, axis; 
  int nbins; double lower, upper; 

  string yaxis="# of B candidates";

  if(todr=="ntracks"){
    //mass
    toDraw="ntracks";
    axis="B candidate # of tracks";
    nbins=15; lower=2.0, upper=17.0; 
    cut="BCands.distSig3D>5.0 && abs(BCands.eta)<2.0 && BCands.jet30==1 && BCands.ptHardestPJ>84";
  }



  TFile *fMC1 = TFile::Open((dir+dirMC+"qcdpt15_withNTR.root").c_str()); 
  TFile *fMC2 = TFile::Open((dir+dirMC+"qcdpt30_withNTR_125of127.root").c_str()); 
  TFile *fMC3 = TFile::Open((dir+dirMC+"qcdpt80_withNTR.root").c_str()); 
  TFile *fMC4 = TFile::Open((dir+dirMC+"qcdpt170_withNTR_79of80.root").c_str()); 

  TFile *fDATA = TFile::Open((dir+dirDATA+"data_ALL.root").c_str()); 
  
  TTree *tMC1 = (TTree*)fMC1->Get("bcanalyzer/tBCands");
  TTree *tMC2 = (TTree*)fMC2->Get("bcanalyzer/tBCands");
  TTree *tMC3 = (TTree*)fMC3->Get("bcanalyzer/tBCands");
  TTree *tMC4 = (TTree*)fMC4->Get("bcanalyzer/tBCands");
  TTree *tDATA = (TTree*)fDATA->Get("bcanalyzer/tBCands");

  TH1F *hMC1_b= new TH1F("hMC1_b","",nbins,lower,upper);
  TH1F *hMC2_b= new TH1F("hMC2_b","",nbins,lower,upper);
  TH1F *hMC3_b= new TH1F("hMC3_b","",nbins,lower,upper);
  TH1F *hMC4_b= new TH1F("hMC4_b","",nbins,lower,upper);
  TH1F *hMC1_c= new TH1F("hMC1_c","",nbins,lower,upper);
  TH1F *hMC2_c= new TH1F("hMC2_c","",nbins,lower,upper);
  TH1F *hMC3_c= new TH1F("hMC3_c","",nbins,lower,upper);
  TH1F *hMC4_c= new TH1F("hMC4_c","",nbins,lower,upper);
  TH1F *hMC1_l= new TH1F("hMC1_l","",nbins,lower,upper);
  TH1F *hMC2_l= new TH1F("hMC2_l","",nbins,lower,upper);
  TH1F *hMC3_l= new TH1F("hMC3_l","",nbins,lower,upper);
  TH1F *hMC4_l= new TH1F("hMC4_l","",nbins,lower,upper);

  TH1F *hMCcombined_b= new TH1F("hMCcombined_b","",nbins,lower,upper);
  TH1F *hMCcombined_c= new TH1F("hMCcombined_c","",nbins,lower,upper);
  TH1F *hMCcombined_l= new TH1F("hMCcombined_l","",nbins,lower,upper);
  TH1F *hDATA= new TH1F("hDATA","",nbins,lower,upper);

  TCanvas *c = new TCanvas("quantFLAV",axis.c_str(),100,300,500,500);
  c->SetLeftMargin(0.1229839);

  tMC1->Draw(("BCands." + toDraw + ">>hMC1_b").c_str(),(cut+" && BCands.flavor==5").c_str());
  tMC2->Draw(("BCands." + toDraw + ">>hMC2_b").c_str(),(cut+" && BCands.flavor==5").c_str());
  tMC3->Draw(("BCands." + toDraw + ">>hMC3_b").c_str(),(cut+" && BCands.flavor==5").c_str());
  tMC4->Draw(("BCands." + toDraw + ">>hMC4_b").c_str(),(cut+" && BCands.flavor==5").c_str());
  tMC1->Draw(("BCands." + toDraw + ">>hMC1_c").c_str(),(cut+" && BCands.flavor==4").c_str());
  tMC2->Draw(("BCands." + toDraw + ">>hMC2_c").c_str(),(cut+" && BCands.flavor==4").c_str());
  tMC3->Draw(("BCands." + toDraw + ">>hMC3_c").c_str(),(cut+" && BCands.flavor==4").c_str());
  tMC4->Draw(("BCands." + toDraw + ">>hMC4_c").c_str(),(cut+" && BCands.flavor==4").c_str());
  tMC1->Draw(("BCands." + toDraw + ">>hMC1_l").c_str(),(cut+" && BCands.flavor==1").c_str());
  tMC2->Draw(("BCands." + toDraw + ">>hMC2_l").c_str(),(cut+" && BCands.flavor==1").c_str());
  tMC3->Draw(("BCands." + toDraw + ">>hMC3_l").c_str(),(cut+" && BCands.flavor==1").c_str());
  tMC4->Draw(("BCands." + toDraw + ">>hMC4_l").c_str(),(cut+" && BCands.flavor==1").c_str());
  tDATA->Draw(("BCands." + toDraw + ">>hDATA").c_str(),(cut).c_str());


//   hMC1_b->SetFillColor(0);
//   hMC2_b->SetFillColor(0);
//   hMC3_b->SetFillColor(0);
//   hMC4_b->SetFillColor(0);
//   hMCcombined_b->SetFillColor(0);
//   hDATA->SetFillColor(0);

  hMC1_b->Scale(143.866);
  hMC2_b->Scale(12.1072);
  hMC3_b->Scale(0.310862);
  hMC4_b->Scale(0.0082391);
  hMC1_c->Scale(143.866);
  hMC2_c->Scale(12.1072);
  hMC3_c->Scale(0.310862);
  hMC4_c->Scale(0.0082391);
  hMC1_l->Scale(143.866);
  hMC2_l->Scale(12.1072);
  hMC3_l->Scale(0.310862);
  hMC4_l->Scale(0.0082391);

  hMCcombined_b->Add(hMC1_b);
  hMCcombined_b->Add(hMC2_b);
  hMCcombined_b->Add(hMC3_b);
  hMCcombined_b->Add(hMC4_b);
  hMCcombined_c->Add(hMC1_c);
  hMCcombined_c->Add(hMC2_c);
  hMCcombined_c->Add(hMC3_c);
  hMCcombined_c->Add(hMC4_c);
  hMCcombined_l->Add(hMC1_l);
  hMCcombined_l->Add(hMC2_l);
  hMCcombined_l->Add(hMC3_l);
  hMCcombined_l->Add(hMC4_l);


  double totB = hMCcombined_b->Integral(); 
  double totC = hMCcombined_c->Integral(); 
  double totL = hMCcombined_l->Integral(); 
  double totAll = totB + totC + totL; 

  double fracB = ((int)(totB/totAll*1000.0))/10.0; 
  double fracC = ((int)(totC/totAll*1000.0))/10.0; 
  double fracL = ((int)(totL/totAll*1000.0))/10.0; 
  
  std::stringstream sdb; sdb << fracB;
  std::stringstream sdc; sdc << fracC;
  std::stringstream sdl; sdl << fracL;


  std::cout << "numbers test bcl " << totB << " " << totC << " " << totL << "\n"; 
  std::cout << "numbers test bcl " << fracB << " " << fracC << " " << fracL << "\n"; 
    

  //adding up the mc files: 
  hMCcombined_c->Add(hMCcombined_l); 
  hMCcombined_b->Add(hMCcombined_c); 

  std::cout << "numbers test tot " << hMCcombined_b->Integral() << std::endl;

  //norm  
  if(luminorm){
    hMCcombined_b->Scale(lumiDATA); 
    hMCcombined_c->Scale(lumiDATA); 
    hMCcombined_l->Scale(lumiDATA); 

  }
  else{
    double intMC=hMCcombined_b->Integral(); 
    double intDATA=hDATA->Integral(); 
    hMCcombined_b->Scale(intDATA/intMC); 
    hMCcombined_c->Scale(intDATA/intMC); 
    hMCcombined_l->Scale(intDATA/intMC); 

  }

  //numbers, cuts
  std::cout << "\ncut used: " << cut << std::endl;
  std::cout << "data lumi (/pb): " << lumiDATA << std::endl;
  std::cout << "Numbers\n";
  std::cout << "MC:\t" << hMCcombined_b->Integral() << std::endl;
  std::cout << "DATA:\t" << hDATA->Integral() << std::endl;

  //format
  if(log) c->SetLogy();
  double dmax = hMCcombined_b->GetMaximum(); 
  if(dmax<hDATA->GetMaximum()) dmax = hDATA->GetMaximum(); 
  if(log)  hMCcombined_b->GetYaxis()->SetRangeUser(0.05,dmax*10.0); 
  else hMCcombined_b->GetYaxis()->SetRangeUser(0.05,dmax+dmax/10.0); 
  
  hMCcombined_b->GetXaxis()->SetTitle(axis.c_str()); 
  hMCcombined_b->GetYaxis()->SetTitle(yaxis.c_str()); 
  hMCcombined_b->GetYaxis()->SetTitleOffset(1.3);
  hMCcombined_b->SetFillColor(kRed); 
  hMCcombined_b->SetLineColor(kRed); 
  hMCcombined_c->SetFillColor(kGreen); 
  hMCcombined_c->SetLineColor(kGreen); 
  hMCcombined_l->SetFillColor(kBlue); 
  hMCcombined_l->SetLineColor(kBlue); 
  hDATA->SetMarkerStyle(21); 
  

//   histoQb->SetLineColor(kRed+1);
//   histoQc->SetLineColor(kGreen+2);
//   histoQl->SetLineColor(kBlue+1);


  //draw
  hMCcombined_b->Draw();
  hMCcombined_c->Draw("sames");
  hMCcombined_l->Draw("sames");
  hDATA->Draw("samesPE");

  //leg
  TLegend *leg = new TLegend(0.70,0.75,0.9,0.9); 
  leg->SetFillColor(0); 
  leg->AddEntry(hDATA,"Data","PE");
  if(haveFracInLeg){
    leg->AddEntry(hMCcombined_b,("Beauty " + sdb.str() +" %").c_str(),"F");
    leg->AddEntry(hMCcombined_c,("Charm " + sdc.str() +" %").c_str(),"F");
    leg->AddEntry(hMCcombined_l,("Light " + sdl.str() +" %").c_str(),"F");
  }
  else{
    leg->AddEntry(hMCcombined_b,"B","F");
    leg->AddEntry(hMCcombined_c,"D","F");
    leg->AddEntry(hMCcombined_l,"Light","F");
  }
  leg->Draw(); 


  //cms label 
  double axismax  = hMCcombined_b->GetMaximum();
  TLatex *   tex = new TLatex(2,axismax+axismax/50.0,"CMS #sqrt{s} = 7 TeV, L = 3 pb^{-1}");
  tex->SetTextSize(0.040); //0.044
  tex->SetLineWidth(2);
  tex->Draw();


}
