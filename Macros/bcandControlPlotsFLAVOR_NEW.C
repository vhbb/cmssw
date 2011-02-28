#include "vector.h"
#include "TFile.h"
#include <string>
#include <sstream>
#include <iostream>



void controlPlots(string todr="dist3D"){
  gStyle->SetOptStat(0);
  gROOT->ProcessLine(".L ~/tdrstyle.C");
  setTDRStyle();
  //style modifications:
  gStyle->SetErrorX(0.5);
  gStyle->SetPadTopMargin(0.08);//.005
  gStyle->SetPadBottomMargin(0.10);//0.13
  gStyle->SetMarkerStyle(1); 

  bool haveFracInLeg = 0; 

//   string dir="/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/BBC/CMSSW_3_7_0_patch2/crab/analyze/rootfiles/"; 
//   string dirMC="mcBins/withIP/"; 
//   string dirDATA="dataFiles/withIP/";

  string dirMC="/scratch/leo/";
//   string dirDATA="/scratch/wehrlilu/BBCORR/MACROS/TREES/DATA/"; 
  string dirDATA="dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/wehrlilu/FILES/DATA/"; 

  string cut="BCands.selected==1 && BCands.jet30==1 && BCands.ptHardestPJ>84";


  //LUMI DATA FOR JET30: 0.313154 JET50U: 3.06903 ((IN /PB))
  //MC LUMI IS 1/PB
  bool luminorm=0; double lumiDATA=0.313154;
  bool log=1; 


  string toDraw, axis; 
  int nbins; double lower, upper; 

  string yaxis="number of B candidates";

  if(todr=="mass"){
    //mass
    toDraw="massBcand";
    axis="vertex mass (GeV)";
    nbins=60; lower=0.0, upper=8.0; 
    cut="BCands.distSig3D>5.0 && abs(BCands.eta)<2.0 && BCands.pt>8 && BCands.jet30==1 && BCands.ptHardestPJ>84";
  }

  if(todr=="dist3D"){
  //dist3D
  string toDraw="dist3D";
  string axis="3D flight distance (cm)";
//   int nbins=50; double lower=0.0, upper=10.0; 
  int nbins=40; double lower=0.0, upper=4.0; 
  }

  if(todr=="distSig3D"){
  //distSig3D
  string toDraw="distSig3D";
  string axis="3D flight distance significance";
  int nbins=400; double lower=0.0, upper=300.0; 
  cut="BCands.massBcand>1.4 && abs(BCands.eta)<2.0 && BCands.pt>8 && BCands.jet30==1 && BCands.ptHardestPJ>84";
  }

  if(todr=="pt"){
  //distSig3D
  string toDraw="pt";
  string axis="transverse momentum (GeV)";
  int nbins=40; double lower=0.0, upper=300.0; 
  }



  TFile *fMC1 = TFile::Open((dirMC+"anV3-QCD_Pt15_Spring10-V8b.root").c_str()); 
  TFile *fMC2 = TFile::Open((dirMC+"anV3-QCD_Pt30_Spring10-V8b.root").c_str()); 
  //TFile *fMC2 = TFile::Open((dirMC+"anV3-InclusiveBB_Pt30_Spring10-V8b.root").c_str()); 
  TFile *fMC3 = TFile::Open((dirMC+"anV3-QCD_Pt80_Spring10-V8b.root").c_str()); 
  TFile *fMC4 = TFile::Open((dirMC+"anV3-QCD_Pt170_Spring10-V8b.root").c_str()); 
  TFile *fMC5 = TFile::Open((dirMC+"anV3-QCD_Pt300_Spring10-V8b.root").c_str()); 


  TFile *fDATA = TFile::Open((dirDATA+"mergedDATA.root").c_str()); 
  
  TTree *tMC1 = (TTree*)fMC1->Get("bcanalyzer/tBCands");
  TTree *tMC2 = (TTree*)fMC2->Get("bcanalyzer/tBCands");
  TTree *tMC3 = (TTree*)fMC3->Get("bcanalyzer/tBCands");
  TTree *tMC4 = (TTree*)fMC4->Get("bcanalyzer/tBCands");
  TTree *tMC5 = (TTree*)fMC5->Get("bcanalyzer/tBCands");
  TTree *tDATA = (TTree*)fDATA->Get("bcanalyzer/tBCands");

  TH1F *hMC1_b= new TH1F("hMC1_b","",nbins,lower,upper);
  TH1F *hMC2_b= new TH1F("hMC2_b","",nbins,lower,upper);
  TH1F *hMC3_b= new TH1F("hMC3_b","",nbins,lower,upper);
  TH1F *hMC4_b= new TH1F("hMC4_b","",nbins,lower,upper);
  TH1F *hMC5_b= new TH1F("hMC5_b","",nbins,lower,upper);
  TH1F *hMC1_c= new TH1F("hMC1_c","",nbins,lower,upper);
  TH1F *hMC2_c= new TH1F("hMC2_c","",nbins,lower,upper);
  TH1F *hMC3_c= new TH1F("hMC3_c","",nbins,lower,upper);
  TH1F *hMC4_c= new TH1F("hMC4_c","",nbins,lower,upper);
  TH1F *hMC5_c= new TH1F("hMC5_c","",nbins,lower,upper);
  TH1F *hMC1_l= new TH1F("hMC1_l","",nbins,lower,upper);
  TH1F *hMC2_l= new TH1F("hMC2_l","",nbins,lower,upper);
  TH1F *hMC3_l= new TH1F("hMC3_l","",nbins,lower,upper);
  TH1F *hMC4_l= new TH1F("hMC4_l","",nbins,lower,upper);
  TH1F *hMC5_l= new TH1F("hMC5_l","",nbins,lower,upper);

  TH1F *hMCcombined_b= new TH1F("hMCcombined_b","",nbins,lower,upper);
  TH1F *hMCcombined_c= new TH1F("hMCcombined_c","",nbins,lower,upper);
  TH1F *hMCcombined_l= new TH1F("hMCcombined_l","",nbins,lower,upper);
  TH1F *hDATA= new TH1F("hDATA","",nbins,lower,upper);

  TCanvas *c = new TCanvas("quantFLAV",axis.c_str(),100,300,500,500);
  c->SetLeftMargin(0.1229839);
  //c->SetBottomMargin(0.1); 
  //c->Range(-1.266217,-1.997143,9.02958,4.963989);

  tMC1->Draw(("BCands." + toDraw + ">>hMC1_b").c_str(),(cut+" && BCands.flavor==5 && BCands.pthat<30").c_str());
  tMC2->Draw(("BCands." + toDraw + ">>hMC2_b").c_str(),(cut+" && BCands.flavor==5 && BCands.pthat<80").c_str());
  tMC3->Draw(("BCands." + toDraw + ">>hMC3_b").c_str(),(cut+" && BCands.flavor==5 && BCands.pthat<170").c_str());
  tMC4->Draw(("BCands." + toDraw + ">>hMC4_b").c_str(),(cut+" && BCands.flavor==5").c_str());
//   tMC5->Draw(("BCands." + toDraw + ">>hMC5_b").c_str(),(cut+" && BCands.flavor==5").c_str());
  tMC1->Draw(("BCands." + toDraw + ">>hMC1_c").c_str(),(cut+" && BCands.flavor==4 && BCands.pthat<30").c_str());
  tMC2->Draw(("BCands." + toDraw + ">>hMC2_c").c_str(),(cut+" && BCands.flavor==4 && BCands.pthat<80").c_str());
  tMC3->Draw(("BCands." + toDraw + ">>hMC3_c").c_str(),(cut+" && BCands.flavor==4 && BCands.pthat<170").c_str());
  tMC4->Draw(("BCands." + toDraw + ">>hMC4_c").c_str(),(cut+" && BCands.flavor==4").c_str());
//   tMC5->Draw(("BCands." + toDraw + ">>hMC5_c").c_str(),(cut+" && BCands.flavor==4").c_str());
  tMC1->Draw(("BCands." + toDraw + ">>hMC1_l").c_str(),(cut+" && BCands.flavor==1 && BCands.pthat<30").c_str());
  tMC2->Draw(("BCands." + toDraw + ">>hMC2_l").c_str(),(cut+" && BCands.flavor==1 && BCands.pthat<80").c_str());
  tMC3->Draw(("BCands." + toDraw + ">>hMC3_l").c_str(),(cut+" && BCands.flavor==1 && BCands.pthat<170").c_str());
  tMC4->Draw(("BCands." + toDraw + ">>hMC4_l").c_str(),(cut+" && BCands.flavor==1").c_str());
//   tMC5->Draw(("BCands." + toDraw + ">>hMC5_l").c_str(),(cut+" && BCands.flavor==1").c_str());
  tDATA->Draw(("BCands." + toDraw + ">>hDATA").c_str(),(cut).c_str());


//   hMC1_b->SetFillColor(0);
//   hMC2_b->SetFillColor(0);
//   hMC3_b->SetFillColor(0);
//   hMC4_b->SetFillColor(0);
//   hMCcombined_b->SetFillColor(0);
//   hDATA->SetFillColor(0);

  hMC1_b->Scale(143.866);
  hMC2_b->Scale(12.1072);
  //hMC2_b->Scale(4.06863); //incl. BB
  hMC3_b->Scale(0.310862);
  hMC4_b->Scale(0.0082391);
  hMC1_c->Scale(143.866);
  hMC2_c->Scale(12.1072);
  //hMC2_c->Scale(4.06863); //incl. BB
  hMC3_c->Scale(0.310862);
  hMC4_c->Scale(0.0082391);
  hMC1_l->Scale(143.866);
  hMC2_l->Scale(12.1072);
  //hMC2_l->Scale(4.06863); //incl. BB
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
//   hMCcombined_b->GetYaxis()->SetTitleOffset(1.3);
  hMCcombined_b->GetYaxis()->SetTitleOffset(1.1);
  //new
  hMCcombined_b->GetYaxis()->SetTitleSize(0.05);
  hMCcombined_b->GetYaxis()->SetLabelSize(0.04);
  hMCcombined_b->GetXaxis()->SetTitleSize(0.0435);
  hMCcombined_b->GetXaxis()->SetLabelSize(0.04);


  Int_t ci; 
  ci = TColor::GetColor("#ff7777"); 
  hMCcombined_b->SetFillColor(2); //kRed or ci
  hMCcombined_b->SetLineColor(2); //kRed or ci (bright red)
  ci = TColor::GetColor("#00cc00");
  hMCcombined_c->SetFillColor(ci); //kGreen
  hMCcombined_c->SetLineColor(ci); //kGreen
  ci = TColor::GetColor("#AABBff");
  hMCcombined_l->SetFillColor(ci); //kBlue 
  hMCcombined_l->SetLineColor(ci); //kBlue
  hDATA->SetMarkerStyle(21); 
  

//   histoQb->SetLineColor(kRed+1);
//   histoQc->SetLineColor(kGreen+2);
//   histoQl->SetLineColor(kBlue+1);

  TLegend *leg; 


  //changes because of npad (to be made before drawing)
  if(todr=="distSig3D"){
    hMCcombined_b->SetMinimum(0.05);
    hMCcombined_b->SetMaximum(7070);
    TH1F *hMCcombined_bi =  hMCcombined_b->Clone(); 
    TH1F *hMCcombined_ci =  hMCcombined_c->Clone(); 
    TH1F *hMCcombined_li =  hMCcombined_l->Clone(); 
    TH1F *hDATAi = hDATA->Clone(); 
    hMCcombined_b->Rebin(10); 
    hMCcombined_c->Rebin(10); 
    hMCcombined_l->Rebin(10); 
    hDATA->Rebin(10); 
    hMCcombined_b->GetXaxis()->SetRange(1,39);
    leg = new TLegend(0.328629,0.720339,0.5302419,0.8898305,NULL,"brNDC");

  }
  else{
    leg = new TLegend(0.6,0.7,0.8,0.85,NULL,"brNDC");
  }
  //draw
  hMCcombined_b->DrawCopy();
  hMCcombined_c->DrawCopy("sames");
  hMCcombined_l->DrawCopy("sames");
  hDATA->DrawCopy("samesPE");

  //leg
  leg->SetBorderSize(0);
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
//   double axismax  = hMCcombined_b->GetMaximum();
  if(todr=="mass") TLatex *   tex = new TLatex(0,35566.42,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
  else if(todr=="dist3D"){
    TLatex *   tex = new TLatex(0,36247.64,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
  }
  else if(todr=="distSig3D"){
    //TLatex *   tex = new TLatex(0,36247.64,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
    TLatex *   tex = new TLatex(0,13333.58,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
    TPad *npad = new TPad("npad", "",0.5604839,0.5063559,0.9475806,0.8855932);
    npad->Draw();
    npad->cd();
    npad->Range(-2.161763,-106.7171,25.79117,954.5241);
    npad->SetFillColor(10);
    npad->SetBorderMode(0);
    npad->SetBorderSize(2);
    npad->SetTickx(1);
    npad->SetTicky(1);
    npad->SetRightMargin(0.01041655);
    npad->SetTopMargin(0.01117313);
    npad->SetFrameFillColor(0);
    npad->SetFrameFillStyle(0);
    npad->SetFrameBorderMode(0);
    npad->SetFrameFillStyle(0);
    npad->SetFrameBorderMode(0);

    hMCcombined_bi->SetMinimum(0);
    hMCcombined_bi->SetMaximum(942.6667); 
    hMCcombined_bi->GetXaxis()->SetRange(2,34);
    hMCcombined_bi->GetXaxis()->SetNdivisions(3);
    hMCcombined_bi->GetXaxis()->SetLabelFont(42);
    hMCcombined_bi->GetXaxis()->SetLabelOffset(0.007);
    hMCcombined_bi->GetXaxis()->SetLabelSize(0.06);
    hMCcombined_bi->GetXaxis()->SetTitleSize(0.06);
    hMCcombined_bi->GetXaxis()->SetTitleOffset(0.9);
    hMCcombined_bi->GetXaxis()->SetTitleFont(42);
    hMCcombined_bi->GetYaxis()->SetNdivisions(5);
    hMCcombined_bi->GetYaxis()->SetLabelFont(42);
    hMCcombined_bi->GetYaxis()->SetLabelOffset(0.007);
    hMCcombined_bi->GetYaxis()->SetLabelSize(0.06);
    hMCcombined_bi->GetYaxis()->SetTitleSize(0.05);
    hMCcombined_bi->GetYaxis()->SetTitleOffset(1.1);
    hMCcombined_bi->GetYaxis()->SetTitleFont(42);
    hMCcombined_bi->GetZaxis()->SetLabelFont(42);
    hMCcombined_bi->GetZaxis()->SetLabelOffset(0.007);
    hMCcombined_bi->GetZaxis()->SetLabelSize(0.05);
    hMCcombined_bi->GetZaxis()->SetTitleSize(0.06);
    hMCcombined_bi->GetZaxis()->SetTitleFont(42);
    hMCcombined_bi->GetXaxis()->SetTitle("");
    hMCcombined_bi->GetYaxis()->SetTitle("");

    hDATAi->SetMarkerSize(0.6);

    hMCcombined_bi->Draw("");
    hMCcombined_ci->Draw("sames");
    hMCcombined_li->Draw("sames");
    hDATAi->Draw("samesPE");
    npad->Modified();
   
    c->cd(); 
  }
  else{
    TLatex *   tex = new TLatex(0,1000000,"CMS #sqrt{s} = 7 TeV, L = 3 pb^{-1}");
  }
  tex->SetTextSize(0.040); //0.044
  tex->SetLineWidth(2);
  tex->Draw();
  gPad->RedrawAxis();


}
