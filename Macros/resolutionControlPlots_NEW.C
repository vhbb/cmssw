#include "vector.h"
#include "TFile.h"
#include <string>
#include <sstream>
#include <iostream>


//give as arguments 
//twoD (bool): to plot the two dim. resolution (true) or the projection on diagonal axis (false)
//drawAll: set to false (true only if you want to draw contribution from different bins)
//strigger (string): jet15, jet30,jet50 works...

void controlPlots(bool twoD = true, bool drawAll = false, string strigger="jet30"){

//   string strigger = "jet30";
  string ljptcut = ""; 
  if(strigger=="jet15") ljptcut = "56";
  if(strigger=="jet30") ljptcut = "84";
  if(strigger=="jet50") ljptcut = "120";

  std::cout << "strigger " << strigger << " lj pt cut " << ljptcut << std::endl; 

  if(twoD) gStyle->SetOptStat(0);
  else {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
  }
  

  string dirMC="/scratch/leo/";
  string cut="BEvents.selected && BEvents.nV==2 && BEvents.dRvv>-1 && BEvents.dRbb>-1 && BEvents.nB==2 && BEvents." + strigger + "==1 && BEvents.ptHardestPJ>" + ljptcut + " && (BEvents.massV1+BEvents.massV2)>4.5 && BEvents.ptV1>8 && BEvents.ptV2>8";

//   string cut="BCands.selected==1 && BCands.ptHardestPJ>84";

  bool log=1; 


  string toDraw; 
  string xaxis, yaxis;
  int nbins; double lower, upper; 
  double alimit = 0.0; 
  if(twoD){
    toDraw="dRvv:BEvents.dRbb";
    xaxis="#Delta R between B and Bbar";
    yaxis="#Delta R between vertices"; 
//     nbins=60; lower=0.0, upper=6.0; 
    nbins=12; lower=0.0, upper=4.8; 
  }
  else{
    toDraw="dRvv-BEvents.dRbb";
    xaxis="#Delta R_{vv}-#Delta R_{bb}";
    yaxis=""; 
    nbins=1500; lower=-6.0, upper=6.0; 
    alimit = 0.4; 
    log = 0; 
  }
  TFile *fMC1 = TFile::Open((dirMC+"anV3-QCD_Pt15_Spring10-V8b.root").c_str()); 
  TFile *fMC2 = TFile::Open((dirMC+"anV3-QCD_Pt30_Spring10-V8b.root").c_str()); 
  TFile *fMC3 = TFile::Open((dirMC+"anV3-QCD_Pt80_Spring10-V8b.root").c_str()); 
  TFile *fMC4 = TFile::Open((dirMC+"anV3-QCD_Pt170_Spring10-V8b.root").c_str()); 

  
  TTree *tMC1 = (TTree*)fMC1->Get("bcanalyzer/tEvents");
  TTree *tMC2 = (TTree*)fMC2->Get("bcanalyzer/tEvents");
  TTree *tMC3 = (TTree*)fMC3->Get("bcanalyzer/tEvents");
  TTree *tMC4 = (TTree*)fMC4->Get("bcanalyzer/tEvents");


  if(twoD){
    TH2F *hMC1= new TH2F("hMC1","",nbins,lower,upper,nbins,lower,upper);
    TH2F *hMC2= new TH2F("hMC2","",nbins,lower,upper,nbins,lower,upper);
    TH2F *hMC3= new TH2F("hMC3","",nbins,lower,upper,nbins,lower,upper);
    TH2F *hMC4= new TH2F("hMC4","",nbins,lower,upper,nbins,lower,upper);
    TH2F *hMCcombined= new TH2F("hMCcombined","",nbins,lower,upper,nbins,lower,upper);
  }
  else{
    TH1F *hMC1= new TH1F("hMC1","",nbins,lower,upper);
    TH1F *hMC2= new TH1F("hMC2","",nbins,lower,upper);
    TH1F *hMC3= new TH1F("hMC3","",nbins,lower,upper);
    TH1F *hMC4= new TH1F("hMC4","",nbins,lower,upper);
    TH1F *hMCcombined= new TH1F("hMCcombined","",nbins,lower,upper);
  }

  TCanvas *c = new TCanvas("quant","resolution",100,300,500,500);

  tMC1->Draw(("BEvents." + toDraw + ">>hMC1").c_str(),(cut).c_str());
  tMC2->Draw(("BEvents." + toDraw + ">>hMC2").c_str(),(cut).c_str());
  tMC3->Draw(("BEvents." + toDraw + ">>hMC3").c_str(),(cut).c_str());
  tMC4->Draw(("BEvents." + toDraw + ">>hMC4").c_str(),(cut).c_str());

  std::cout << "int qcd15 " << hMC1->Integral() << std::endl;
  std::cout << "int qcd30 " << hMC2->Integral() << std::endl;
  std::cout << "int qcd80 " << hMC3->Integral() << std::endl;
  std::cout << "int qcd170 " << hMC4->Integral() << std::endl;

//   hMC1->Scale(143.866);
//   hMC2->Scale(12.1072);
//   hMC3->Scale(0.310862);
//   hMC4->Scale(0.0082391);
//LEO TREES:
  hMC1->Scale(143.866);
  hMC2->Scale(12.2051);
  hMC3->Scale(0.310862);
  hMC4->Scale(0.00813388);

  std::cout << "int qcd15 " << hMC1->Integral() << std::endl;
  std::cout << "int qcd30 " << hMC2->Integral() << std::endl;
  std::cout << "int qcd80 " << hMC3->Integral() << std::endl;
  std::cout << "int qcd170 " << hMC4->Integral() << std::endl;


  if(drawAll){
    if(twoD){
      TCanvas *c1 = new TCanvas("quant15","resolution",100,300,500,500);
      c1->SetLogz();
      hMC1->GetZaxis()->SetRangeUser(0.1,700); 
      hMC1->Draw("colz");
      hMC1->Draw("textsame");
      TCanvas *c2 = new TCanvas("quant30","resolution",100,300,500,500);
      c2->SetLogz();
      hMC2->GetZaxis()->SetRangeUser(0.1,700); 
      hMC2->Draw("colz");
      hMC2->Draw("textsame");
      TCanvas *c3 = new TCanvas("quant80","resolution",100,300,500,500);
      c3->SetLogz();
      hMC3->GetZaxis()->SetRangeUser(0.1,700); 
      hMC3->Draw("colz");
      hMC3->Draw("textsame");
      TCanvas *c4 = new TCanvas("quant170","resolution",100,300,500,500);
      c4->SetLogz();
      hMC4->GetZaxis()->SetRangeUser(0.1,700); 
      hMC4->Draw("colz");
      hMC4->Draw("textsame");
    }
    else {
      TCanvas *c1 = new TCanvas("quant15","resolution",100,300,500,500);
      c1->SetLogy();
      hMC1->GetYaxis()->SetRangeUser(0.005,4000); 
      hMC1->Draw("");
      TCanvas *c2 = new TCanvas("quant30","resolution",100,300,500,500);
      c2->SetLogy();
      hMC2->GetYaxis()->SetRangeUser(0.005,4000); 
      hMC2->Draw("");
      TCanvas *c3 = new TCanvas("quant80","resolution",100,300,500,500);
      c3->SetLogy();
      hMC3->GetYaxis()->SetRangeUser(0.005,4000); 
      hMC3->Draw("");
      TCanvas *c4 = new TCanvas("quant170","resolution",100,300,500,500);
      c4->SetLogy();
      hMC4->GetYaxis()->SetRangeUser(0.005,4000); 
      hMC4->Draw("");
    }
  }

  c->cd();
  hMCcombined->Add(hMC1);
  hMCcombined->Add(hMC2);
  hMCcombined->Add(hMC3);
  hMCcombined->Add(hMC4);



  //numbers, cuts
  std::cout << "\ncut used: " << cut << std::endl;
  std::cout << "Numbers\n";
  std::cout << "MC:\t" << hMCcombined->Integral() << std::endl;

  //format
  if(log && twoD) c->SetLogz();
  if(log && !twoD) c->SetLogy(); 

  hMCcombined->GetXaxis()->SetTitle(xaxis.c_str()); 
  if(twoD) hMCcombined->GetYaxis()->SetTitle(yaxis.c_str()); 
  
  ///////
  if(!twoD){
    //count amount of events being outside diagonal region...
    std::cout << "nbins " << nbins << "\n";
    double tot=0.0, totDRl04 =0.0; 
    for(unsigned int i=0; i<nbins; i++){
      tot += hMCcombined->GetBinContent(i); 
      if(fabs(hMCcombined->GetBinCenter(i))>0.4) totDRl04+=hMCcombined->GetBinContent(i); 
    }
    std::cout << "TOT " << tot << "\nLE 0.4 " << totDRl04 << " " << totDRl04/tot*100.0 << "\n";
  }
  ////////


  //draw
  if(twoD) {
    quant->SetRightMargin(0.12);
    hMCcombined->GetYaxis()->SetTitleOffset(1.2);
    hMCcombined->GetXaxis()->SetTitleOffset(1.2);
    hMCcombined->SetMarkerSize(1.4);
    hMCcombined->DrawCopy("colz");
    for(unsigned int x=0; x<=nbins; x++){
      for(unsigned int y=0; y<=nbins; y++){
// 	std::cout << x << " " << y << " " << hMCcombined->GetBinContent(x,y) << std::endl;
	hMCcombined->SetBinContent(x,y,(int)(hMCcombined->GetBinContent(x,y)));
// 	std::cout << x << " " << y << " " << hMCcombined->GetBinContent(x,y) << std::endl; 
      }
    }
    hMCcombined->Draw("textsame");
  }
  else {

//     hMCcombined->GetXaxis()->SetRangeUser(-alimit,alimit);
    hMCcombined->Draw("");
    hMCcombined->Fit("gaus");
    if(gaus->GetMaximum()>hMCcombined->GetMaximum()){
      hMCcombined->GetYaxis()->SetRangeUser(0,gaus->GetMaximum()+gaus->GetMaximum()/20.0); 
    }
//     std::cout << "fit max " << gaus->GetMaximum() << std::endl;
//     TPaveStats *st = (TPaveStats*)hMCcombined->FindObject("stats"); 
//     st->SetX1NDC(0.1);
//     st->Draw();
  }



  //cms label 
  double axismax; 
  if(twoD){
    double axismax  = 5;
    TLatex *   tex = new TLatex(0.0,axismax,"CMS #sqrt{s} = 7 TeV, Simulation");
  }
  else{
    double axismax  = hMCcombined->GetMaximum(); 
    if(log) TLatex *   tex = new TLatex(-alimit-alimit/10.0,axismax*3,"CMS #sqrt{s} = 7 TeV, Simulation");
    else{
      TLatex *   tex = new TLatex(-alimit-alimit/20.,axismax+axismax/50.0,"CMS #sqrt{s} = 7 TeV, Simulation");
    }
  }
  tex->SetTextSize(0.040); //0.044
  tex->SetLineWidth(2);
  tex->Draw();


}
