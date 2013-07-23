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
  gROOT->ProcessLine(".L ~/tdrstyle.C");
  setTDRStyle();

  //style modifications:
  gStyle->SetErrorX(0.5);
  gStyle->SetPadTopMargin(0.08);//.005
  gStyle->SetPadBottomMargin(0.10);//0.13
  gStyle->SetMarkerStyle(1); 
  gStyle->SetPalette(1);
  gStyle->SetFuncWidth(2);

//   string strigger = "jet30";
  string ljptcut = ""; 
  if(strigger=="jet6") ljptcut = "27";
  if(strigger=="jet15") ljptcut = "56";
  if(strigger=="jet30") ljptcut = "84";
  if(strigger=="jet50") ljptcut = "120";

  std::cout << "strigger " << strigger << " lj pt cut " << ljptcut << std::endl; 

  if(twoD) gStyle->SetOptStat(0);
  else {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
  }
  

  string dirMC="/scratch/leo/";
  string dirMCSE="dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/b-physics/bbcorr/mc/";

  string cut="BEvents.selected && BEvents.nV==2 && BEvents.dRvv>-1 && BEvents.dRbb>-1 && BEvents.nB==2 && BEvents." + strigger + "==1 && BEvents.ptHardestPJ>" + ljptcut + " && (BEvents.massV1+BEvents.massV2)>4.5 && BEvents.ptV1>8 && BEvents.ptV2>8";

//   string cut="BCands.selected==1 && BCands.ptHardestPJ>84";

  bool log=1; 


  string toDraw; 
  string xaxis, yaxis;
  int nbins; double lower, upper; 
  double alimit = 0.0; 
  if(twoD){
    toDraw="dRvv:BEvents.dRbb";
    xaxis="#Delta R_{BB}";
    yaxis="#Delta R_{VV}"; 
//     nbins=60; lower=0.0, upper=6.0; 
    nbins=12; lower=0.0, upper=4.8; 
//    nbins=24; lower=0.0, upper=4.8; 
  }
  else{
    toDraw="dRvv-BEvents.dRbb";
    xaxis="#Delta R_{VV} - #Delta R_{BB}";
    yaxis="number of events"; 
    nbins=3000; lower=-6.0, upper=6.0; 
    alimit = 0.2; 
    log = 0; 
  }

TFile *fMC1 = TFile::Open((dirMC+"anV3-QCD_Pt15_Spring10-V8b.root").c_str()); 
//TFile *fMC1 = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/b-physics/bbcorr/mc/anV5-QCD_Pt-15to30_bEnriched_TuneZ2_7TeV-pythia6-evtgen-V11.root");
//nfs/psi.ch/cms/trivcat//store/user/b-physics/bbcorr/anV5-QCD_Pt-15to30_bEnriched_TuneZ2_7TeV-pythia6-evtgen-V11.root");
//TFile *fMC1 = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/b-physics/bbcorr/mc/anV5-QCD_Pt_15_TuneZ2_7TeV_pythia6-Summer10-V8b.root");


//  TFile *fMC2 = TFile::Open((dirMC+"anV3-QCD_Pt30_Spring10-V8b.root").c_str()); 
//srm://t3se01.psi.ch:8443/pnfs/psi.ch/cms/trivcat//store/user/b-physics/bbcorr//anV5-QCD_Pt-30to50_bEnriched_TuneZ2_7TeV-pythia6-evtgen-V11.root
  TFile *fMC2 = TFile::Open((dirMCSE+"anV3-InclusiveBB_Pt30_Spring10-V8b.root").c_str());
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

  tMC1->Draw(("BEvents." + toDraw + ">>hMC1").c_str(),(cut+" && BEvents.pthat<30").c_str());
  tMC2->Draw(("BEvents." + toDraw + ">>hMC2").c_str(),(cut+" && BEvents.pthat<80").c_str());
  tMC3->Draw(("BEvents." + toDraw + ">>hMC3").c_str(),(cut+" && BEvents.pthat<170").c_str());
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
//  hMC1->Scale(7.6);
//  hMC1->Scale(0);
//  hMC2->Scale(12.2051);
  hMC2->Scale(4.1);
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
  hMCcombined->GetYaxis()->SetTitle(yaxis.c_str()); 
  
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
    hMCcombined->GetYaxis()->SetTitleOffset(1.1);
    hMCcombined->GetXaxis()->SetTitleOffset(1.);
    //new
    hMCcombined->GetYaxis()->SetTitleSize(0.05); 
    hMCcombined->GetYaxis()->SetLabelSize(0.04);
    hMCcombined->GetXaxis()->SetTitleSize(0.0435); 
    hMCcombined->GetXaxis()->SetLabelSize(0.04);
    hMCcombined->GetZaxis()->SetLabelSize(0.04);

    hMCcombined->SetMarkerSize(1.4);
    hMCcombined->GetYaxis()->SetRangeUser(0.001, hMCcombined->GetMaximum()); 
    hMCcombined->DrawCopy("colz");
    double diag24 = 0;
    double diag08 = 0;
    double migrated08to24 = hMCcombined->Integral(1,2,7,12); 
    double migrated24to08 = hMCcombined->Integral(7,12,1,2); 
    double tot = hMCcombined->Integral(0,12,0,12); 
    diag = 0;
    double tot2 = 0;
    for(unsigned int x=0; x<=nbins; x++){
      for(unsigned int y=0; y<=nbins; y++){
//        if(x==y && x <=2 && x >=1) diag08+=hMCcombined->GetBinContent(x,y);
//        if(x==y && x >=7 && x <=12) diag24+=hMCcombined->GetBinContent(x,y);
        if(y<=2 && y>=1 && x <=2 && x >=1) diag08+=hMCcombined->GetBinContent(x,y);
        if(y<=12 && y>=7 && x >=7 && x <=12) diag24+=hMCcombined->GetBinContent(x,y);
        if(x==y || abs(x-y) <=1) diag+=hMCcombined->GetBinContent(x,y);
        tot2+=hMCcombined->GetBinContent(x,y);

// 	std::cout << x << " " << y << " " << hMCcombined->GetBinContent(x,y) << std::endl;
//         if(hMCcombined->GetBinContent(x,y) < 0.1) //0.1
// 	hMCcombined->SetBinContent(x,y,0);

	hMCcombined->SetBinContent(x,y,(int)(hMCcombined->GetBinContent(x,y)));


//TMath::Nint(hMCcombined->GetBinContent(x,y)));
// 	std::cout << x << " " << y << " " << hMCcombined->GetBinContent(x,y) << std::endl; 
      }
    }
//     gStyle->SetPaintTextFormat("1.1f"); //1.1f
    hMCcombined->Draw("textsame");
    cout << "Migration Summary " << strigger << " (integral: " << tot << ")" <<  endl;
    cout << " < 0.8 to > 2.4 : " << migrated08to24/diag24 << "   N: " << migrated08to24 << " " << diag24  << endl;
    cout << " > 2.4 to < 0.8 : " << migrated24to08/diag08 << "   N: " << migrated24to08 << " " << diag08 << endl;
    cout << "on-diagonal    : " << diag/tot << "   N:" << diag << " / " << tot <<  " " << tot2 << endl;

  }
  else {

//     hMCcombined->GetXaxis()->SetRangeUser(-alimit,alimit);
    hMCcombined->GetYaxis()->SetTitleSize(0.05); 
    hMCcombined->GetYaxis()->SetLabelSize(0.04);
    hMCcombined->GetXaxis()->SetTitleSize(0.0435); 
    hMCcombined->GetXaxis()->SetLabelSize(0.04);
    hMCcombined->Draw("");
    /*hMCcombined->Fit("gaus");
    if(gaus->GetMaximum()>hMCcombined->GetMaximum()){
      hMCcombined->GetYaxis()->SetRangeUser(0,gaus->GetMaximum()+gaus->GetMaximum()/20.0); 
      }*/
//     std::cout << "fit max " << gaus->GetMaximum() << std::endl;
//     TPaveStats *st = (TPaveStats*)hMCcombined->FindObject("stats"); 
//     st->SetX1NDC(0.1);
//     st->Draw();
  }



  //cms label 
  double axismax; 
  if(twoD){
    double axismax  = 5;
    TLatex *   tex = new TLatex(0.0,axismax,"CMS    #sqrt{s} = 7 TeV, Simulation");
  }
  else{
    double axismax  = hMCcombined->GetMaximum(); 
    if(log) TLatex *   tex = new TLatex(-alimit-alimit/10.0,axismax*3,"CMS    #sqrt{s} = 7 TeV, Simulation");
    else{
      std::cout << "axislimit " << alimit << " " << axismax << std::endl;
      TLatex *   tex = new TLatex(-alimit,axismax+axismax/11.0,"CMS    #sqrt{s} = 7 TeV, Simulation");
      /*TLatex *   tex2 = new TLatex(0.02, 250,"#splitline{Mean (0.8 #pm 2.7) #times 10^{-4}}{Sigma 0.0132 #pm 0.0003}");
      tex2->SetTextSize(0.035); //0.044
      tex2->SetLineWidth(2);
      tex2->Draw(); */
      hMCcombined->GetXaxis()->SetRangeUser(-0.2,0.2); 
      hMCcombined->GetYaxis()->SetTitleOffset(1.3); 
    }
  }
  tex->SetTextSize(0.040); //0.044
  tex->SetLineWidth(2);
  tex->Draw();

  c->RedrawAxis();
  c->SetTickx(); 
  c->SetTicky();

  c->SaveAs(("resolutions_"+strigger+((twoD)?"_n":"_1d")+".png").c_str());
  c->SaveAs(("resolutions_"+strigger+((twoD)?"_n":"_1d")+".root").c_str());

}
