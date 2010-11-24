#include "vector.h"
#include "TFile.h"
#include <string>
#include <sstream>
#include <iostream>

void resultPlots(string todraw="dR", bool corrected=true, string output="finalPlots.root"){

  ////
  //non ratio plots
  double scaleFactorNR = 2.0;

  //ratio plots
  //log
  bool ratioPlotsInLog = false;  
  double ratioFactorLog = 100.0; 
  //lin
  //lumi
  double DratioOffsetLin_lumi = 2.0; 
  if(todraw=="dR"){
    TH1F *ratioOffsetLin_lumi = new TH1F("offsetlumi","offset",15,0,6); 
    for(unsigned int i=0; i<=15; i++) {
      ratioOffsetLin_lumi->SetBinContent(i,DratioOffsetLin_lumi);  
      ratioOffsetLin_lumi->SetBinError(i,0.0);  
    }
    //others
    double DratioOffsetLin = 6.0; 
    TH1F *ratioOffsetLin = new TH1F("offset","offset",15,0,6); 
    for(unsigned int i=0; i<=15; i++) {
      ratioOffsetLin->SetBinContent(i,DratioOffsetLin);  
      ratioOffsetLin->SetBinError(i,0.0);  
    }
  }
  else if(todraw=="dPhi"){
    TH1F *ratioOffsetLin_lumi = new TH1F("offsetlumi","offset",8,0,3.2); 
    for(unsigned int i=0; i<=8; i++) {
      ratioOffsetLin_lumi->SetBinContent(i,DratioOffsetLin_lumi);  
      ratioOffsetLin_lumi->SetBinError(i,0.0);  
    }
    //others
    double DratioOffsetLin = 4.0; 
    TH1F *ratioOffsetLin = new TH1F("offset","offset",8,0,3.2); 
    for(unsigned int i=0; i<=8; i++) {
      ratioOffsetLin->SetBinContent(i,DratioOffsetLin);  
      ratioOffsetLin->SetBinError(i,0.0);  
    }
  }
  ////

  ////
  //style
  gStyle->SetEndErrorSize(3); 
//   gStyle->SetOptStat(0);
//   gStyle->SetPadTickX(1); 
//   gStyle->SetPadTickY(1);
  ////

  ////
  //mc event generators to include
  bool bmadgraph=true; int madgCol = kBlue-7; int madgColF = kBlue-7; //4
  bool bmcatnlo=true; int mcatCol = kRed-7; int mcatColF = kRed-7; //2
  bool bcascade=true; int cascCol = kMagenta-7; int cascColF = kMagenta-7; //6
  string mc_draw_option = "E3same"; //"histLsame","E2same","histCsame"
  int boxcolor = kGray; //kRed-10

  ////

  ////
  //NOT SO NICE...
  //Aequivalent lumi fro different triggers in inv. pb
  double lumi50=3069.0/1000.0; //jet50
  double lumi30=313.0/1000.0; //jet30
  double lumi15=31.0/1000.0; //jet15
  ////

  ////
  //rebinning
  double rebinData=4; 
  double rebin=4;
  ////

  ////
  //how to handle error of correction: 
  bool addFlatErrorForCorrectionSyst = true;
  ////


  //scaled to 1/pb
//   string fileMC="/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/BBC/CMSSW_3_7_0_patch2/crab/analyze/rootoutput/temp8215_finalBinEffic/defaultAddedBins.root";
//   string fileDATA="/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/BBC/CMSSW_3_7_0_patch2/crab/analyze/rootoutput/temp8215_finalBinEffic/defaultAddedDATABins.root";
  string fileMC15="/scratch/wehrlilu/BBCORR/MACROS/SUMFILES/MC/histo-Pythia6-HLT_Jet15U-BBCorr_scaledEff-total_ETACUTOPEN.root";
  string fileMC30="/scratch/wehrlilu/BBCORR/MACROS/SUMFILES/MC/histo-Pythia6-HLT_Jet30U-BBCorr_scaledEff-total_ETACUTOPEN.root";
  string fileMC50="/scratch/wehrlilu/BBCORR/MACROS/SUMFILES/MC/histo-Pythia6-HLT_Jet50U-BBCorr_scaledEff-total_ETACUTOPEN.root";

  string fileDATA15="/scratch/wehrlilu/BBCORR/MACROS/SUMFILES/DATA/histo-Data-HLT_Jet15U-BBCorr_scaledEff-total_ETACUTOPEN.root";
  string fileDATA30="/scratch/wehrlilu/BBCORR/MACROS/SUMFILES/DATA/histo-Data-HLT_Jet30U-BBCorr_scaledEff-total_ETACUTOPEN.root";
  string fileDATA50="/scratch/wehrlilu/BBCORR/MACROS/SUMFILES/DATA/histo-Data-HLT_Jet50U-BBCorr_scaledEff-total_ETACUTOPEN.root";



  ////
  //open files
  TFile *ffinal; 
  ffinal = TFile::Open(output.c_str(),"RECREATE");
  TFile *fDATA15 = TFile::Open(fileDATA15.c_str(),"READ");
  TFile *fDATA30 = TFile::Open(fileDATA30.c_str(),"READ");
  TFile *fDATA50 = TFile::Open(fileDATA50.c_str(),"READ");
  TFile *fMC15 = TFile::Open(fileMC15.c_str(),"READ");
  TFile *fMC30 = TFile::Open(fileMC30.c_str(),"READ");
  TFile *fMC50 = TFile::Open(fileMC50.c_str(),"READ");
  if(fDATA15->IsZombie() || fMC15->IsZombie() || fDATA30->IsZombie() || fMC30->IsZombie() || fDATA50->IsZombie() || fMC50->IsZombie()){
    std::cout << "At least one DATA and/or MC file can not be opened\n";
    return;
  }
  ////

  ////
  //get histos
  TH1F *data15, *data30, *data50; 
  TH1F *data15_err, *data30_err, *data50_err; 
  TH1F *data15_toterr, *data30_toterr, *data50_toterr; 
  TH1F *data15_toterrLUMI, *data30_toterrLUMI, *data50_toterrLUMI; 
  TH1F *mc15, *mc30, *mc50; 
  TH1F *mc15_a, *mc30_a, *mc50_a; 
  TH1F *mc15_b, *mc30_b, *mc50_b; 
  if(!corrected){
    fDATA15->GetObject(("Vert_" + todraw + "_2b").c_str(),data15); data15->Sumw2();
    fDATA30->GetObject(("Vert_" + todraw + "_2b").c_str(),data30); data30->Sumw2();
    fDATA50->GetObject(("Vert_" + todraw + "_2b").c_str(),data50); data50->Sumw2();
    fMC15->GetObject(("Vert_" + todraw + "_2b").c_str(),mc15); mc15->Sumw2();
    fMC30->GetObject(("Vert_" + todraw + "_2b").c_str(),mc30); mc30->Sumw2();
    fMC50->GetObject(("Vert_" + todraw + "_2b").c_str(),mc50); mc50->Sumw2();
  }
  else{
    fDATA15->GetObject(("Vert_" + todraw + "_2b_CORR").c_str(),data15); data15->Sumw2();
    fDATA30->GetObject(("Vert_" + todraw + "_2b_CORR").c_str(),data30); data30->Sumw2();
    fDATA50->GetObject(("Vert_" + todraw + "_2b_CORR").c_str(),data50); data50->Sumw2();
    fDATA15->GetObject(("Vert_" + todraw + "_2b_ErrCORR").c_str(),data15_err); data15_err->Sumw2();
    fDATA30->GetObject(("Vert_" + todraw + "_2b_ErrCORR").c_str(),data30_err); data30_err->Sumw2();
    fDATA50->GetObject(("Vert_" + todraw + "_2b_ErrCORR").c_str(),data50_err); data50_err->Sumw2();
    fMC15->GetObject(("Vert_" + todraw + "_efficden").c_str(),mc15); mc15->Sumw2();
    fMC30->GetObject(("Vert_" + todraw + "_efficden").c_str(),mc30); mc30->Sumw2();
    fMC50->GetObject(("Vert_" + todraw + "_efficden").c_str(),mc50); mc50->Sumw2();
  }
  //histos for total errors (not only stat ones)
  data15_toterr = (TH1F*)data15->Clone();
  data30_toterr = (TH1F*)data30->Clone(); 
  data50_toterr = (TH1F*)data50->Clone(); 
  data15_toterrLUMI = (TH1F*)data15->Clone();
  data30_toterrLUMI = (TH1F*)data30->Clone(); 
  data50_toterrLUMI = (TH1F*)data50->Clone(); 
  ////


  ////
  //rebin data/mc
  data15->Rebin(rebinData);
  data30->Rebin(rebinData);
  data50->Rebin(rebinData);
  data15_toterr->Rebin(rebinData);
  data30_toterr->Rebin(rebinData);
  data50_toterr->Rebin(rebinData);
  data15_toterrLUMI->Rebin(rebinData);
  data30_toterrLUMI->Rebin(rebinData);
  data50_toterrLUMI->Rebin(rebinData);
  mc15->Rebin(rebin);
  mc30->Rebin(rebin);
  mc50->Rebin(rebin);
  ////


  ////
  //first rebin, then add errors from error histo and other syst error and then rescale: 
  double systErr = 0.113; 
  double systErrCORRECTION = 0.13; 
  double systErrLUMI = 0.4; 
  bool addErrCorrectionQuadratic = true; 
  for(unsigned int i=0; i<data15_toterr->GetNbinsX(); i++){
    //LATER
    if(addFlatErrorForCorrectionSyst){
	double toterr = sqrt(data15_toterr->GetBinError(i)*data15_toterr->GetBinError(i) +  (systErrCORRECTION*data15_toterr->GetBinContent(i))*(systErrCORRECTION*data15_toterr->GetBinContent(i)) +  (systErr*data15_toterr->GetBinContent(i))*(systErr*data15_toterr->GetBinContent(i)));
	data15_toterr->SetBinError(i,toterr); 
	double toterr = sqrt(data30_toterr->GetBinError(i)*data30_toterr->GetBinError(i) +  (systErrCORRECTION*data30_toterr->GetBinContent(i))*(systErrCORRECTION*data30_toterr->GetBinContent(i)) +  (systErr*data30_toterr->GetBinContent(i))*(systErr*data30_toterr->GetBinContent(i)));
	data30_toterr->SetBinError(i,toterr); 
	double toterr = sqrt(data50_toterr->GetBinError(i)*data50_toterr->GetBinError(i) +  (systErrCORRECTION*data50_toterr->GetBinContent(i))*(systErrCORRECTION*data50_toterr->GetBinContent(i)) +  (systErr*data50_toterr->GetBinContent(i))*(systErr*data50_toterr->GetBinContent(i)));
	data50_toterr->SetBinError(i,toterr); 

    }
    else{
      if(addErrCorrectionQuadratic){ //add squares of error of correction
	double toterr = sqrt(data15_toterr->GetBinError(i)*data15_toterr->GetBinError(i) + data15_err->GetBinContent(i)*data15_err->GetBinContent(i) +  (systErr*data15_toterr->GetBinContent(i))*(systErr*data15_toterr->GetBinContent(i)));
	data15_toterr->SetBinError(i,toterr); 
	
	double toterr = sqrt(data30_toterr->GetBinError(i)*data30_toterr->GetBinError(i) + data30_err->GetBinContent(i)*data30_err->GetBinContent(i) +  (systErr*data30_toterr->GetBinContent(i))*(systErr*data30_toterr->GetBinContent(i))); 
	data30_toterr->SetBinError(i,toterr); 
	
	double toterr = sqrt(data50_toterr->GetBinError(i)*data50_toterr->GetBinError(i) + data50_err->GetBinContent(i)*data50_err->GetBinContent(i) +  (systErr*data50_toterr->GetBinContent(i))*(systErr*data50_toterr->GetBinContent(i))); 
    data50_toterr->SetBinError(i,toterr); 
      }
      else{ //add error linear
	double toterr = sqrt(data15_toterr->GetBinError(i)*data15_toterr->GetBinError(i) + data15_err->GetBinContent(i) +  (systErr*data15_toterr->GetBinContent(i))*(systErr*data15_toterr->GetBinContent(i))); 
	data15_toterr->SetBinError(i,toterr); 
	
	double toterr = sqrt(data30_toterr->GetBinError(i)*data30_toterr->GetBinError(i) + data30_err->GetBinContent(i) +  (systErr*data30_toterr->GetBinContent(i))*(systErr*data30_toterr->GetBinContent(i))); 
	data30_toterr->SetBinError(i,toterr); 
	
	double toterr = sqrt(data50_toterr->GetBinError(i)*data50_toterr->GetBinError(i) + data50_err->GetBinContent(i) +  (systErr*data50_toterr->GetBinContent(i))*(systErr*data50_toterr->GetBinContent(i))); 
      data50_toterr->SetBinError(i,toterr); 
      }
    }

    //add error for absolute (LUMI) normalization
    double toterrLUMI = systErrLUMI*data15_toterrLUMI->GetBinContent(i);  
    toterrLUMI = sqrt(toterrLUMI*toterrLUMI+data15_toterr->GetBinError(i)*data15_toterr->GetBinError(i)); 
    data15_toterrLUMI->SetBinError(i,toterrLUMI); 

    double toterrLUMI = systErrLUMI*data30_toterrLUMI->GetBinContent(i);  
    toterrLUMI = sqrt(toterrLUMI*toterrLUMI+data30_toterr->GetBinError(i)*data30_toterr->GetBinError(i)); 
    data30_toterrLUMI->SetBinError(i,toterrLUMI); 

    double toterrLUMI = systErrLUMI*data50_toterrLUMI->GetBinContent(i); 
    toterrLUMI = sqrt(toterrLUMI*toterrLUMI+data50_toterr->GetBinError(i)*data50_toterr->GetBinError(i)); 
    data50_toterrLUMI->SetBinError(i,toterrLUMI); 

  }
  ////

//   std::cout << "LUMITEST " << data15_toterrLUMI->GetBinError(1) << " " << data15_toterrLUMI->GetBinContent(1) << std::endl;


  ////
  //rescale to 1/pb and divide by bin width
  double d = data15->GetBinWidth(1); 
  unsigned int nbins = data50->GetNbinsX(); 
  unsigned int nbinsMC = mc50->GetNbinsX(); 
  data15->Scale(1./lumi15/d); 
  data30->Scale(1./lumi30/d); 
  data50->Scale(1./lumi50/d);
  data15_toterr->Scale(1./lumi15/d); 
  data30_toterr->Scale(1./lumi30/d); 
  data50_toterr->Scale(1./lumi50/d); 
  data15_toterrLUMI->Scale(1./lumi15/d); 
  data30_toterrLUMI->Scale(1./lumi30/d); 
  data50_toterrLUMI->Scale(1./lumi50/d); 
  //mc files are already scaled to 1/pb
  mc15->Scale(1./d); 
  mc30->Scale(1./d); 
  mc50->Scale(1./d); 
  ////

  ////
  //add additional scale factor
  data15->Scale(scaleFactorNR*scaleFactorNR); 
  data15_toterr->Scale(scaleFactorNR*scaleFactorNR); 
  data15_toterrLUMI->Scale(scaleFactorNR*scaleFactorNR); 
  mc15->Scale(scaleFactorNR*scaleFactorNR); 
  data30->Scale(scaleFactorNR);
  data30_toterr->Scale(scaleFactorNR);
  data30_toterrLUMI->Scale(scaleFactorNR);
  mc30->Scale(scaleFactorNR);  
  ////


  ////
  //set style   
  data15->SetStats(kFALSE);
  data15->SetMarkerStyle(21); 
  data30->SetMarkerStyle(25); 
  data50->SetMarkerStyle(20); 
  data15_toterr->SetMarkerStyle(21); 
  data30_toterr->SetMarkerStyle(25); 
  data50_toterr->SetMarkerStyle(20); 
//   data15_toterrLUMI->SetMarkerStyle(21); 
//   data30_toterrLUMI->SetMarkerStyle(25); 
//   data50_toterrLUMI->SetMarkerStyle(20); 
  data15_toterrLUMI->SetFillColor(kYellow-7); 
  data30_toterrLUMI->SetFillColor(kYellow-7); 
  data50_toterrLUMI->SetFillColor(kYellow-7); 
  mc15->SetFillColor(8);
  mc15->SetLineWidth(3);
  mc15->SetLineColor(8);
  mc15->SetMarkerColor(8);
  mc30->SetFillColor(8);
  mc30->SetLineWidth(3);
  mc30->SetLineColor(8);
  mc30->SetMarkerColor(8);
  mc50->SetFillColor(8);
  mc50->SetLineWidth(3);
  mc50->SetLineColor(8);
  mc50->SetMarkerColor(8);
  ////

  ////
  //set range user (min, max)
  double max = -1; 
  if(data15->GetMaximum()>max) max = data15->GetMaximum();
  if(mc15->GetMaximum()>max) max = mc15->GetMaximum();
  double min=1000;  
  for(unsigned int i=0; i<nbins; i++)
    if(data50->GetBinContent(i)>0 && data50->GetBinContent(i)<min) 
      min=data50->GetBinContent(i); 
  
  for(unsigned int i=0; i<nbinsMC; i++)
    if(mc50->GetBinContent(i)>0 && mc50->GetBinCenter(i)<4.5 && mc50->GetBinContent(i)<min) 
      min=mc50->GetBinContent(i); 

  data15->GetXaxis()->SetRangeUser(0,4.39);
  //data15->GetYaxis()->SetRangeUser(min*0.2,20*scaleFactorNR*max); //ohne scale factor: min*0.5, 20*max | mit scale factor 10: min*0.2, 200*max
  double axismax = 2*max; 
  std::cout << "AXISMAX " << axismax << std::endl;
  data15->GetYaxis()->SetRangeUser(min*0.2,axismax); 

  ////

  ////
  //set axis properties
  data15->GetYaxis()->SetTitleOffset(1.5);
  data15->SetTitle("");
  if(todraw=="dR") {
    if(!corrected) data15->GetYaxis()->SetTitle("#frac{d #sigma}{d #Delta R} uncorrected [pb]");
    else data15->GetYaxis()->SetTitle("#frac{d #sigma}{d #Delta R} [pb]");
    data15->GetXaxis()->SetTitle("#Delta R");
  }
  else if (todraw=="dPhi") {
    if(!corrected) data15->GetYaxis()->SetTitle("#frac{d #sigma}{d #Delta #phi} uncorrected [pb]");
    else data15->GetYaxis()->SetTitle("#frac{d #sigma}{d #Delta #phi} [pb]");
    data15->GetXaxis()->SetTitle("#Delta #phi");
  }
  else if (todraw=="dEta") {
    if(!corrected) data15->GetYaxis()->SetTitle("#frac{d #sigma}{d #Delta #eta} uncorrected [pb]");
    else data15->GetYaxis()->SetTitle("#frac{d #sigma}{d #Delta #eta} [pb]");
    data15->GetXaxis()->SetTitle("#Delta #eta");
  }
  ////

  ////
  //Cloning for other normalizations
  mc15_a = (TH1F)mc15->Clone();
  mc30_a = (TH1F)mc30->Clone();
  mc50_a = (TH1F)mc50->Clone();
  mc15_b = (TH1F)mc15->Clone();
  mc30_b = (TH1F)mc30->Clone();
  mc50_b = (TH1F)mc50->Clone();
  ////
       

  /////////////////
  //LUMINORM
  TCanvas *cucorrcomp1 = new TCanvas(("lumi_"+todraw).c_str(),"DATA vs. MC before correction (lumi)",100,100,500,500);
  cucorrcomp1->SetLogy();
  cucorrcomp1->SetLeftMargin(0.16);

//   std::cout << "LUMITEST " << data15_toterrLUMI->GetBinError(1) << " " << data15_toterrLUMI->GetBinContent(1) << std::endl;

  data15->DrawCopy("E1"); 
  data15_toterrLUMI->Draw("E3same"); 
  data30_toterrLUMI->Draw("E3same"); 
  data50_toterrLUMI->Draw("E3same"); 
  mc15->DrawCopy("E2same"); 
  mc30->DrawCopy("E2same"); 
  mc50->DrawCopy("E2same"); 
  data15_toterr->DrawCopy("E1same"); 
  data30_toterr->DrawCopy("E1same"); 
  data50_toterr->DrawCopy("E1same");
  data15->DrawCopy("E1same"); 
  data30->DrawCopy("E1same"); 
  data50->DrawCopy("E1same"); 

  ////
  //Legend and CMS label
  //TLegend *leg = new TLegend(0.5,0.72,0.9,0.9);
  TLegend *leg = new TLegend(0.17,0.15,0.55,0.35);
  leg->SetFillColor(0);
  leg->SetShadowColor(0);
  leg->SetBorderSize(0); 
//   leg->AddEntry(data15,"leading jet p_{t}>56 GeV","PE");
//   leg->AddEntry(data30,"leading jet p_{t}>84 GeV","PE");
//   leg->AddEntry(data50,"leading jet p_{t}>120 GeV","PE");
  std::stringstream ssSFNR; ssSFNR << scaleFactorNR;
  std::stringstream ssSFNR2; ssSFNR2 << scaleFactorNR*scaleFactorNR; 

  leg->AddEntry(data15,("leading jet p_{t}>56 GeV (#times " + ssSFNR2.str() + ")").c_str(),"PE");
  leg->AddEntry(data30,("leading jet p_{t}>84 GeV (#times " + ssSFNR.str() + ")").c_str(),"PE");
  leg->AddEntry(data50,"leading jet p_{t}>120 GeV","PE");
  leg->AddEntry(mc15,"Pythia","l");
  leg->Draw();  
  TLatex *   tex = new TLatex(0.0,axismax*1.5,"CMS #sqrt{s} = 7 TeV, L = 3 pb^{-1}"); //without scale factor 21 | with scale factor 210
  tex->SetTextSize(0.040); //0.044
  tex->SetLineWidth(2);
  tex->Draw();
  gPad->RedrawAxis();
  ////
  /////////////////

  ////
  //area normalizatin
  double factor = (double)nbinsMC/(double)nbins; 
  double intd15 = data15->Integral();
  double intm15 = mc15_a->Integral()/factor;
  double intd30 = data30->Integral();
  double intm30 = mc30_a->Integral()/factor;
  double intd50 = data50->Integral();
  double intm50 = mc50_a->Integral()/factor;
  mc15_a->Scale(intd15/intm15);
  mc30_a->Scale(intd30/intm30);
  mc50_a->Scale(intd50/intm50);
  std::cout << "int dat " << data15->Integral() << std::endl;
  ////

  /////////////////
  //AREANORM
  TCanvas *cucorrcomp1_1 = new TCanvas(("final_area_"+todraw).c_str(),"DATA vs. MC after correction (area)",100,300,500,500);
  cucorrcomp1_1->SetLogy();
  cucorrcomp1_1->SetLeftMargin(0.16);

  data15->DrawCopy("E1"); 
  mc15_a->DrawCopy("E2same"); 
  mc30_a->DrawCopy("E2same"); 
  mc50_a->DrawCopy("E2same"); 
  data15_toterr->DrawCopy("E1same"); 
  data30_toterr->DrawCopy("E1same"); 
  data50_toterr->DrawCopy("E1same"); 
  data15->DrawCopy("E1same"); 
  data30->DrawCopy("E1same"); 
  data50->DrawCopy("E1same");
  leg->Draw();
  tex->Draw();
  gPad->RedrawAxis();
  /////////////////
  
//------------------------------------------------------------------------
//other generators, init and make area normalized plots

  /////////////////////////////
  //MADGRAPH
  if(bmadgraph && corrected){
    TFile *fmadgraph = TFile::Open("/scratch/wehrlilu/BBCORR/MACROS/OTHERMC/MADGRAPH/MADGRAPHaddedBins.root","");

    TH1F *mg15_a, *mg30_a, *mg50_a;
    TH1F *mg15_b, *mg30_b, *mg50_b;

//     fmadgraph15->GetObject("sum4",mg15_a); mg15_a->Sumw2();
//     fmadgraph30->GetObject("sum4",mg30_a); mg30_a->Sumw2();
//     fmadgraph50->GetObject("sum4",mg50_a); mg50_a->Sumw2();
    fmadgraph->GetObject(("jet15/Vert_" + todraw + "_efficden").c_str(),mg15_a); mg15_a->Sumw2();
    fmadgraph->GetObject(("jet30/Vert_" + todraw + "_efficden").c_str(),mg30_a); mg30_a->Sumw2();
    fmadgraph->GetObject(("jet50/Vert_" + todraw + "_efficden").c_str(),mg50_a); mg50_a->Sumw2();
    mg15_a->Rebin(rebin); 
    mg30_a->Rebin(rebin); 
    mg50_a->Rebin(rebin); 
    mg15_a->SetLineWidth(2); 
    mg30_a->SetLineWidth(2); 
    mg50_a->SetLineWidth(2); 
    mg15_a->SetFillColor(0); 
    mg30_a->SetFillColor(0); 
    mg50_a->SetFillColor(0); 


    mg15_b = (TH1F)mg15_a->Clone();
    mg30_b = (TH1F)mg30_a->Clone(); 
    mg50_b = (TH1F)mg50_a->Clone(); 

    mg15_a->SetLineColor(madgCol); 
    mg30_a->SetLineColor(madgCol); 
    mg50_a->SetLineColor(madgCol); 
    mg15_b->SetLineColor(madgCol); 
    mg30_b->SetLineColor(madgCol); 
    mg50_b->SetLineColor(madgCol); 
    mg15_a->SetMarkerColor(madgCol); 
    mg30_a->SetMarkerColor(madgCol); 
    mg50_a->SetMarkerColor(madgCol); 
    mg15_b->SetMarkerColor(madgCol); 
    mg30_b->SetMarkerColor(madgCol); 
    mg50_b->SetMarkerColor(madgCol); 
    mg15_a->SetFillColor(madgColF); 
    mg30_a->SetFillColor(madgColF); 
    mg50_a->SetFillColor(madgColF); 
    mg15_b->SetFillColor(madgColF); 
    mg30_b->SetFillColor(madgColF); 
    mg50_b->SetFillColor(madgColF); 


    ////
    //normalize to area
    double intmg15_a = mg15_a->Integral(); 
    double intmg30_a = mg30_a->Integral(); 
    double intmg50_a = mg50_a->Integral(); 
    mg15_a->Scale(intd15/intmg15_a);
    mg30_a->Scale(intd30/intmg30_a);
    mg50_a->Scale(intd50/intmg50_a);
    
    //no mg in non-ratio plots
//     mg15_a->DrawCopy("E2same");
//     mg30_a->DrawCopy("E2same");
//     mg50_a->DrawCopy("E2same");
//     leg->AddEntry(mg15_a,"Madgraph","l");
    /////
  }
  /////////////////////////////


  /////////////////////////////
  //MC@NLO
  if(bmcatnlo && corrected){
    TFile *fmcatnlo15 = TFile::Open("/scratch/wehrlilu/BBCORR/MACROS/OTHERMC/MCATNLO/MCATNLO_cobinedEta2_jet15.root","");
    TFile *fmcatnlo30 = TFile::Open("/scratch/wehrlilu/BBCORR/MACROS/OTHERMC/MCATNLO/MCATNLO_cobinedEta2_jet30.root","");
    TFile *fmcatnlo50 = TFile::Open("/scratch/wehrlilu/BBCORR/MACROS/OTHERMC/MCATNLO/MCATNLO_cobinedEta2_jet50.root","");

    TH1F *mcat15_a, *mcat30_a, *mcat50_a;
    TH1F *mcat15_b, *mcat30_b, *mcat50_b;

    fmcatnlo15->GetObject(("h" + todraw + "_mcatnlo").c_str(),mcat15_a); mcat15_a->Sumw2();
    fmcatnlo30->GetObject(("h" + todraw + "_mcatnlo").c_str(),mcat30_a); mcat30_a->Sumw2();
    fmcatnlo50->GetObject(("h" + todraw + "_mcatnlo").c_str(),mcat50_a); mcat50_a->Sumw2();

    mcat15_a->Rebin(4); 
    mcat30_a->Rebin(4); 
    mcat50_a->Rebin(4); 
    mcat15_a->SetLineWidth(2); 
    mcat30_a->SetLineWidth(2); 
    mcat50_a->SetLineWidth(2); 
    mcat15_a->SetFillColor(0); 
    mcat30_a->SetFillColor(0); 
    mcat50_a->SetFillColor(0); 

    mcat15_b = (TH1F)mcat15_a->Clone();
    mcat30_b = (TH1F)mcat30_a->Clone(); 
    mcat50_b = (TH1F)mcat50_a->Clone(); 

    mcat15_a->SetLineColor(mcatCol); 
    mcat30_a->SetLineColor(mcatCol); 
    mcat50_a->SetLineColor(mcatCol); 
    mcat15_b->SetLineColor(mcatCol); 
    mcat30_b->SetLineColor(mcatCol); 
    mcat50_b->SetLineColor(mcatCol); 
    mcat15_a->SetMarkerColor(mcatCol); 
    mcat30_a->SetMarkerColor(mcatCol); 
    mcat50_a->SetMarkerColor(mcatCol); 
    mcat15_b->SetMarkerColor(mcatCol); 
    mcat30_b->SetMarkerColor(mcatCol); 
    mcat50_b->SetMarkerColor(mcatCol); 
    mcat15_a->SetFillColor(mcatColF); 
    mcat30_a->SetFillColor(mcatColF); 
    mcat50_a->SetFillColor(mcatColF); 
    mcat15_b->SetFillColor(mcatColF); 
    mcat30_b->SetFillColor(mcatColF); 
    mcat50_b->SetFillColor(mcatColF); 

    ////
    //normalize to area
    double intmcat15_a = mcat15_a->Integral(); 
    double intmcat30_a = mcat30_a->Integral(); 
    double intmcat50_a = mcat50_a->Integral(); 
    mcat15_a->Scale(intd15/intmcat15_a);
    mcat30_a->Scale(intd30/intmcat30_a);
    mcat50_a->Scale(intd50/intmcat50_a);
    //no mcat in non-ratio plots
//     std::cout << "int mcat " << mcat15_a->Integral() << std::endl;
//     mcat15_a->DrawCopy("E2same");
//     mcat30_a->DrawCopy("E2same");
//     mcat50_a->DrawCopy("E2same");
//     leg->AddEntry(mcat15_a,"MC@NLO","l");
    ////
  }
  /////////////////////////////

  /////////////////////////////
  //CASCADE
  if(bcascade && corrected){
    //beforeAddingGG 
    TFile *fcascade1 = TFile::Open("/scratch/wehrlilu/BBCORR/MACROS/OTHERMC/CASCADE/cascade-bottom-ipro=11-ifps=2-ppNEW.root","");
    //added gg->gg with g->bbbar
    TFile *fcascade2 = TFile::Open("/scratch/wehrlilu/BBCORR/MACROS/OTHERMC/CASCADE/cascade-bottom-ipro=10-irpa=1-irpb=1-irpc=1-ifps=2-scal=0.25-pp.root","");

    TH1F *mcasc15_a, *mcasc30_a, *mcasc50_a;
    TH1F *mcasc15_a1, *mcasc30_a1, *mcasc50_a1; //direct bb
    TH1F *mcasc15_a2, *mcasc30_a2, *mcasc50_a2; //gg->gg with g-> bbbar
    TH1F *mcasc15_b, *mcasc30_b, *mcasc50_b;

    //ETA2 
    if(todraw=="dR"){
      fcascade1->GetObject("BPH10_010/h2101",mcasc15_a1); mcasc15_a1->Sumw2();
      fcascade1->GetObject("BPH10_010/h2102",mcasc30_a1); mcasc30_a1->Sumw2();
      fcascade1->GetObject("BPH10_010/h2103",mcasc50_a1); mcasc50_a1->Sumw2();
      fcascade2->GetObject("BPH10_010/h2101",mcasc15_a2); mcasc15_a2->Sumw2();
      fcascade2->GetObject("BPH10_010/h2102",mcasc30_a2); mcasc30_a2->Sumw2();
      fcascade2->GetObject("BPH10_010/h2103",mcasc50_a2); mcasc50_a2->Sumw2();
      //     //ETA2.4 IS EXACTLY THE SAME!!! (FOR ALL BINS)
      //     fcascade->GetObject("BPH10_010/h2001",mcasc15_a); mcasc15_a->Sumw2();
      //     fcascade->GetObject("BPH10_010/h2002",mcasc30_a); mcasc30_a->Sumw2();
      //     fcascade->GetObject("BPH10_010/h2003",mcasc50_a); mcasc50_a->Sumw2();
    }
    else if(todraw=="dPhi"){
      fcascade1->GetObject("BPH10_010/h3101",mcasc15_a1); mcasc15_a1->Sumw2();
      fcascade1->GetObject("BPH10_010/h3102",mcasc30_a1); mcasc30_a1->Sumw2();
      fcascade1->GetObject("BPH10_010/h3103",mcasc50_a1); mcasc50_a1->Sumw2();
      fcascade2->GetObject("BPH10_010/h3101",mcasc15_a2); mcasc15_a2->Sumw2();
      fcascade2->GetObject("BPH10_010/h3102",mcasc30_a2); mcasc30_a2->Sumw2();
      fcascade2->GetObject("BPH10_010/h3103",mcasc50_a2); mcasc50_a2->Sumw2();
    }

    ////////
    //add direct bb and gsp:
    mcasc15_a = (TH1F)mcasc15_a1->Clone();
    mcasc30_a = (TH1F)mcasc30_a1->Clone();
    mcasc50_a = (TH1F)mcasc50_a1->Clone();
    mcasc15_a->Add(mcasc15_a2); 
    mcasc30_a->Add(mcasc30_a2); 
    mcasc50_a->Add(mcasc50_a2); 
    ////////


    mcasc15_a->Rebin(4); 
    mcasc30_a->Rebin(4); 
    mcasc50_a->Rebin(4); 
    mcasc15_a->SetLineWidth(2); 
    mcasc30_a->SetLineWidth(2); 
    mcasc50_a->SetLineWidth(2); 
    mcasc15_a->SetFillColor(0); 
    mcasc30_a->SetFillColor(0); 
    mcasc50_a->SetFillColor(0); 

    mcasc15_b = (TH1F)mcasc15_a->Clone();
    mcasc30_b = (TH1F)mcasc30_a->Clone(); 
    mcasc50_b = (TH1F)mcasc50_a->Clone(); 

    mcasc15_a->SetLineColor(cascCol); 
    mcasc30_a->SetLineColor(cascCol); 
    mcasc50_a->SetLineColor(cascCol); 
    mcasc15_b->SetLineColor(cascCol); 
    mcasc30_b->SetLineColor(cascCol); 
    mcasc50_b->SetLineColor(cascCol); 
    mcasc15_a->SetFillColor(cascColF); 
    mcasc30_a->SetFillColor(cascColF); 
    mcasc50_a->SetFillColor(cascColF); 
    mcasc15_a->SetMarkerColor(cascCol); 
    mcasc30_a->SetMarkerColor(cascCol); 
    mcasc50_a->SetMarkerColor(cascCol); 
    mcasc15_b->SetFillColor(cascColF); 
    mcasc30_b->SetFillColor(cascColF); 
    mcasc50_b->SetFillColor(cascColF); 
    mcasc15_b->SetMarkerColor(cascCol); 
    mcasc30_b->SetMarkerColor(cascCol); 
    mcasc50_b->SetMarkerColor(cascCol); 

    ////
    //normalize to area
    double intmcasc15_a = mcasc15_a->Integral(); 
    double intmcasc30_a = mcasc30_a->Integral(); 
    double intmcasc50_a = mcasc50_a->Integral(); 
    mcasc15_a->Scale(intd15/intmcasc15_a);
    mcasc30_a->Scale(intd30/intmcasc30_a);
    mcasc50_a->Scale(intd50/intmcasc50_a);
    //no mcasc in non-ratio plots
//     std::cout << "int casc" << mcasc15_a->Integral() << std::endl;
//     mcasc15_a->DrawCopy("E2same");
//     mcasc30_a->DrawCopy("E2same");
//     mcasc50_a->DrawCopy("E2same");
//     leg->AddEntry(mcasc15_a,"Cascade","l");
    ////
  }
  /////////////////////////////
//------------------------------------------------------------------------



  ////
  //back-to-back norm (fcr norm)
  unsigned int firstBin = 100, firstBinMC = 100;
  for(unsigned int i=0; i<nbins; i++)
    if(data15->GetBinCenter(i)>2.4 && firstBin==100)
      firstBin = i; 
  
  for(unsigned int i=0; i<nbinsMC; i++)
    if(mc15_b->GetBinCenter(i)>2.4 && firstBinMC==100)
      firstBinMC = i; 
  
  intd15 = data15->Integral(firstBin,nbins);
  intm15 = mc15_b->Integral(firstBinMC,nbinsMC)/factor;
  intd30 = data30->Integral(firstBin,nbins);
  intm30 = mc30_b->Integral(firstBinMC,nbinsMC)/factor;
  intd50 = data50->Integral(firstBin,nbins);
  intm50 = mc50_b->Integral(firstBinMC,nbinsMC)/factor;

  double scale15, scale30, scale50; 
  std::cout << "firstBin/Center (data   mc) " << firstBin << "/" << data15->GetBinCenter(firstBin) << " " << firstBinMC << "/" << mc15_b->GetBinCenter(firstBinMC) << std::endl;
    scale15 = intd15/intm15; 
    scale30 = intd30/intm30; 
    scale50 = intd50/intm50; 

  std::cout << "SCALEFACTORS B2B\n";
  std::cout << scale15 << std::endl;
  std::cout << scale30 << std::endl;
  std::cout << scale50 << std::endl;


  mc15_b->Scale(scale15);
  mc30_b->Scale(scale30);
  mc50_b->Scale(scale50);
  ////

  /////////////////
  //B2BNORM
  TCanvas *cucorrcomp1_2 = new TCanvas(("final_btb_"+todraw).c_str(),"DATA vs. MC after correction (back-to-back)",100,500,500,500);
  cucorrcomp1_2->SetLogy();
  cucorrcomp1_2->SetLeftMargin(0.16);

  if(todraw=="dR"){
    TH1F *line = new TH1F("line","line",15,0,6); 
    TH1F *sepline = new TH1F("sepline","sepline",15,0,6); 
    for(unsigned int i=0; i<=15; i++) {
      line->SetBinContent(i,1);  
      line->SetBinError(i,0);  
      sepline->SetBinContent(i,0);  
      sepline->SetBinError(i,0);  
    }
    TBox *box = new TBox(2.4,0,4.4,axismax-axismax/15.0);
    box->SetFillColor(boxcolor); 
    box->SetLineColor(boxcolor); 
  }
  else if(todraw=="dPhi"){
    TH1F *line = new TH1F("line","line",8,0,3.2); 
    TH1F *sepline = new TH1F("sepline","sepline",8,0,3.2); 
    for(unsigned int i=0; i<=8; i++) {
      line->SetBinContent(i,1);  
      line->SetBinError(i,0);  
      sepline->SetBinContent(i,0);  
      sepline->SetBinError(i,0);  
    }
    TBox *box = new TBox(2.4,0,3.2,axismax-axismax/15.0);
    box->SetFillColor(boxcolor); 
    box->SetLineColor(boxcolor); 
  }
  data15->DrawCopy("E1"); 
  box->Draw("same"); 
  mc15_b->DrawCopy("E2same"); 
  mc30_b->DrawCopy("E2same"); 
  mc50_b->DrawCopy("E2same");
  data15_toterr->DrawCopy("E1same"); 
  data30_toterr->DrawCopy("E1same"); 
  data50_toterr->DrawCopy("E1same"); 
  data15->DrawCopy("E1same"); 
  data30->DrawCopy("E1same"); 
  data50->DrawCopy("E1same"); 

  TLine *l=new TLine(2.4,0,2.4,100000000);
  l->SetLineWidth(2);
  l->SetLineStyle(2);  
  //l->Draw("same");
  TLegend *legB = leg->Clone(); 
  TH1F *dummy = new TH1F(); dummy->SetFillColor(boxcolor);
  dummy->SetLineColor(boxcolor); 
  legB->AddEntry(dummy, "Normalization region","F");
  legB->Draw();
  tex->Draw();
  gPad->RedrawAxis();
  /////////////////

//------------------------------------------------------------------------  
//other generators, init and make back-to-back normalized plots

  /////////////////////////////
  //MADGRAPH
  if(bmadgraph && corrected){
    double intmg15_b = mg15_b->Integral(firstBin,nbins); 
    double intmg30_b = mg30_b->Integral(firstBin,nbins); 
    double intmg50_b = mg50_b->Integral(firstBin,nbins); 
    mg15_b->Scale(intd15/intmg15_b);
    mg30_b->Scale(intd30/intmg30_b);
    mg50_b->Scale(intd50/intmg50_b);
    //no mg in non-ratio plots
//     mg15_b->DrawCopy("E2same");
//     mg30_b->DrawCopy("E2same");
//     mg50_b->DrawCopy("E2same");
  }
  /////////////////////////////

  /////////////////////////////
  //MC@NLO
  if(bmcatnlo && corrected){
    double intmcat15_b = mcat15_b->Integral(firstBin,nbins); 
    double intmcat30_b = mcat30_b->Integral(firstBin,nbins); 
    double intmcat50_b = mcat50_b->Integral(firstBin,nbins); 
    mcat15_b->Scale(intd15/intmcat15_b);
    mcat30_b->Scale(intd30/intmcat30_b);
    mcat50_b->Scale(intd50/intmcat50_b);
    //no mcat in non-ratio plot
//     mcat15_b->DrawCopy("E2same");
//     mcat30_b->DrawCopy("E2same");
//     mcat50_b->DrawCopy("E2same");
  }
  /////////////////////////////

  /////////////////////////////
  //CASCADE
  if(bcascade && corrected){
    double intmcasc15_b = mcasc15_b->Integral(firstBin,nbins); 
    double intmcasc30_b = mcasc30_b->Integral(firstBin,nbins); 
    double intmcasc50_b = mcasc50_b->Integral(firstBin,nbins); 
    mcasc15_b->Scale(intd15/intmcasc15_b);
    mcasc30_b->Scale(intd30/intmcasc30_b);
    mcasc50_b->Scale(intd50/intmcasc50_b);
    //no mcasc in non-ratio plot
//     mcasc15_b->DrawCopy("E2same");
//     mcasc30_b->DrawCopy("E2same");
//     mcasc50_b->DrawCopy("E2same");
  }
  /////////////////////////////
//------------------------------------------------------------------------

  //***************************
  //***RATIOS
  //***************************

  ////
  //histos
  TH1F *data15_r, *data30_r, *data50_r; 
  data15_r = (TH1F)data15->Clone();
  data30_r = (TH1F)data30->Clone();
  data50_r = (TH1F)data50->Clone();
  TH1F *data15_toterr_r, *data30_toterr_r, *data50_toterr_r; 
  data15_toterr_r = (TH1F)data15_toterr->Clone();
  data30_toterr_r = (TH1F)data30_toterr->Clone();
  data50_toterr_r = (TH1F)data50_toterr->Clone();
  TH1F *data15_toterrLUMI_r, *data30_toterrLUMI_r, *data50_toterrLUMI_r; 
  data15_toterrLUMI_r = (TH1F)data15_toterrLUMI->Clone();
  data30_toterrLUMI_r = (TH1F)data30_toterrLUMI->Clone();
  data50_toterrLUMI_r = (TH1F)data50_toterrLUMI->Clone();
  ////

  ////
  //must have the same number of bins for data/mc histograms to divide
  if(nbins!=nbinsMC){
    std::cout << "factor " << factor << std::endl;
    mc15->Rebin((int)factor);
    mc15->Scale(1./factor);
    mc30->Rebin((int)factor);
    mc30->Scale(1./factor);
    mc50->Rebin((int)factor);
    mc50->Scale(1./factor);
    
    mc15_a->Rebin((int)factor);
    mc15_a->Scale(1./factor);
    mc30_a->Rebin((int)factor);
    mc30_a->Scale(1./factor);
    mc50_a->Rebin((int)factor);
    mc50_a->Scale(1./factor);
    
    mc15_b->Rebin((int)factor);
    mc15_b->Scale(1./factor);
    mc30_b->Rebin((int)factor);
    mc30_b->Scale(1./factor);
    mc50_b->Rebin((int)factor);
    mc50_b->Scale(1./factor);    
  }
  ////
  
  ////
  //divide histos, scale
  if(!corrected) data15_r->GetYaxis()->SetTitle("Ratio to Pythia (uncorrected)");
  else  data15_r->GetYaxis()->SetTitle("Ratio to Pythia"); 
//   //%%%
//   std::cout << "ERROR DIVIDE TEST \nDAT: " << data15_r->GetBinContent(1) << " +- " << data15_r->GetBinError(1) << std::endl;
//   std::cout << "MC : " << mc15->GetBinContent(1) << " +- " << mc15->GetBinError(1) << std::endl;
//   //%%%
  data15_r->Divide(mc15); 
//   //%%%
//   std::cout << "DIV : " << data15_r->GetBinContent(1) << " +- " << data15_r->GetBinError(1) << std::endl;

  data30_r->Divide(mc30); 
  data50_r->Divide(mc50); 
  data15_toterr_r->Divide(mc15); 
  data30_toterr_r->Divide(mc30); 
  data50_toterr_r->Divide(mc50); 
  data15_toterrLUMI_r->Divide(mc15); 
  data30_toterrLUMI_r->Divide(mc30); 
  data50_toterrLUMI_r->Divide(mc50); 

  if(ratioPlotsInLog){
     data15_r->Scale(ratioFactorLog*ratioFactorLog); 
     data30_r->Scale(ratioFactorLog); 
     data15_toterr_r->Scale(ratioFactorLog*ratioFactorLog); 
     data30_toterr_r->Scale(ratioFactorLog); 
     data15_toterrLUMI_r->Scale(ratioFactorLog*ratioFactorLog); 
     data30_toterrLUMI_r->Scale(ratioFactorLog); 
     data15_r->GetYaxis()->SetRangeUser(0.1,10000000);
  }
  else{
    //add offset
    std::cout << "errortest " << data15_r->GetBinError(1)
	      << " " << data15_r->GetBinError(3) << std::endl;
    data15_r->Add(ratioOffsetLin_lumi,2.0); 
    std::cout << "errortest " << data15_r->GetBinError(1)
	      << " " << data15_r->GetBinError(3) << std::endl;
    data30_r->Add(ratioOffsetLin_lumi); 
    data15_toterr_r->Add(ratioOffsetLin_lumi,2.0); 
    data30_toterr_r->Add(ratioOffsetLin_lumi); 
    data15_toterrLUMI_r->Add(ratioOffsetLin_lumi,2.0); 
    data30_toterrLUMI_r->Add(ratioOffsetLin_lumi); 
    data15_r->GetYaxis()->SetRangeUser(0,3*DratioOffsetLin_lumi+2);
  }
  ////

  ////
  //format line
  line->SetLineWidth(2);
  line->SetFillColor(0);
  line->SetLineColor(8);
  TH1F *line2, *line3, *line2_lumi, *line3_lumi, *sepline_lumi, *sepline2, *sepline2_lumi; 
  line2 = (TH1F)line->Clone();
  line3 = (TH1F)line2->Clone();
  line2_lumi = (TH1F)line->Clone();
  line3_lumi = (TH1F)line2->Clone();
  sepline->SetFillColor(0);
  sepline_lumi = (TH1F)sepline->Clone();
  sepline2 = (TH1F)sepline->Clone();
  sepline2_lumi = (TH1F)sepline->Clone();
  if(ratioPlotsInLog){
     line2->Scale(ratioFactorLog); 
     line3->Scale(ratioFactorLog*ratioFactorLog); 
     line2_lumi->Scale(ratioFactorLog); 
     line3_lumi->Scale(ratioFactorLog*ratioFactorLog); 
  }
  else{
    line2_lumi->Add(ratioOffsetLin_lumi); 
    line3_lumi->Add(ratioOffsetLin_lumi,2.0); 
    line2->Add(ratioOffsetLin); 
    line3->Add(ratioOffsetLin,2.0); 
    sepline->Add(ratioOffsetLin);
    sepline2->Add(ratioOffsetLin,2.0);
    sepline_lumi->Add(ratioOffsetLin_lumi);
    sepline2_lumi->Add(ratioOffsetLin_lumi,2.0);
  }
  ////

  std::cout << "LUMINORM RATIO\n";
  /////////////////
  //LUMINORM
  TCanvas *cucorrcomp2 = new TCanvas(("final_ratio_lumi_"+todraw).c_str(),"Ratio plots",700,100,500,500);
  if(ratioPlotsInLog) cucorrcomp2->SetLogy();
  cucorrcomp2->SetLeftMargin(0.14);
  data15_r->GetYaxis()->SetLabelSize(0.0); 
  data15_r->Draw("E1");
  data15_toterrLUMI_r->Draw("E3same");
  data30_toterrLUMI_r->Draw("E3same"); 
  data50_toterrLUMI_r->Draw("E3same"); 
  //lines
  line->DrawCopy("Csame");
  line2_lumi->DrawCopy("Csame");
  line3_lumi->DrawCopy("Csame");
  data15_toterr_r->Draw("E1same");
  data30_toterr_r->Draw("E1same"); 
  data50_toterr_r->Draw("E1same"); 
  data15_r->Draw("E1same");
  data30_r->Draw("E1same"); 
  data50_r->Draw("E1same"); 
  //legend,cms label

  //ratio area & b2b norm 
  std::stringstream ssfact1, ssfact2, ssfact1_lumi, ssfact2_lumi; 
//   if(ratioPlotsInLog) {
//     ssfact1 << " (#times " << ratioFactorLog << ")"; 
//     ssfact2 << " (#times " << ratioFactorLog*ratioFactorLog << ")"; 
//     ssfact1_lumi << " (#times " << ratioFactorLog << ")"; 
//     ssfact2_lumi << " (#times " << ratioFactorLog*ratioFactorLog << ")"; 
//   }
//   else {
//     ssfact1 << " (+ " << DratioOffsetLin << ")"; 
//     ssfact2 << " (+ " << DratioOffsetLin*2.0 << ")"; 
//     ssfact1_lumi << " (+ " << DratioOffsetLin_lumi << ")"; 
//     ssfact2_lumi << " (+ " << DratioOffsetLin_lumi*2.0 << ")"; 
//   }

  TLegend *leg2 = new TLegend(0.5,0.8,0.899,0.89);
  leg2->SetFillColor(0);
  leg2->SetShadowColor(0);
  leg2->SetLineColor(0); 
  leg2->AddEntry(data15,("leading jet p_{t}>56 GeV"+ ssfact2.str()).c_str(),"PE");
  leg2->AddEntry(data30,("leading jet p_{t}>84 GeV" + ssfact1.str()).c_str(),"PE");
  leg2->AddEntry(data50,"leading jet p_{t}>120 GeV","PE");
  //ratio lumi plot
  TLegend *leg2_1 = new TLegend(0.5,0.7,0.89,0.89);
  leg2_1->SetFillColor(0);
  leg2_1->SetShadowColor(0);
  leg2_1->SetBorderSize(0); 
  leg2_1->AddEntry(data15,("leading jet p_{t}>56 GeV" + ssfact2_lumi.str()).c_str(),"PE");
  leg2_1->AddEntry(data30,("leading jet p_{t}>84 GeV" + ssfact1_lumi.str()).c_str(),"PE");
  leg2_1->AddEntry(data50,"leading jet p_{t}>120 GeV","PE");
  leg2_1->AddEntry(mc15,"Pythia","l");
  //ratio area & b2b norm
  TLegend *leg3 = new TLegend(0.2,0.75,0.48,0.89);
  leg3->SetFillColor(0);
  leg3->SetShadowColor(0);
  leg3->SetBorderSize(0);  
  leg3->AddEntry(mc15,"Pythia","l");	
  if(bmadgraph) leg3->AddEntry(mg15_a,"Madgraph","l");
  if(bmcatnlo) leg3->AddEntry(mcat15_a,"MC@NLO","l");
  if(bcascade) leg3->AddEntry(mcasc15_a,"Cascade","l");
  leg2_1->Draw();
  if(ratioPlotsInLog) {
    TLatex *   tex2 = new TLatex(0.15,15000000,"CMS #sqrt{s} = 7 TeV, L = 3 pb^{-1}");
    TLatex *   tex2_lumi = new TLatex(0.15,15000000,"CMS #sqrt{s} = 7 TeV, L = 3 pb^{-1}");
  }
  else {
    TLatex *   tex2_lumi = new TLatex(0.0,DratioOffsetLin_lumi*3.0+2.2,"CMS #sqrt{s} = 7 TeV, L = 3 pb^{-1}");
    TLatex *   tex2 = new TLatex(0.0,DratioOffsetLin*3.0+2.3,"CMS #sqrt{s} = 7 TeV, L = 3 pb^{-1}");
  }
  tex2->SetTextSize(0.040);
  tex2->SetLineWidth(2);
  tex2_lumi->SetTextSize(0.040);
  tex2_lumi->SetLineWidth(2);
  tex2_lumi->Draw();

  //seplines:
  sepline_lumi->Draw("same");
  sepline2_lumi->Draw("same");

  //axis
  //ylabels: (look for data15_r->...->SetLabelSize(0.0); to have def. lab)
  TLatex *t = new TLatex();
  t->SetTextSize(0.042); 
//   t->SetTextFont(42); 
  t->SetTextAlign(31);
  const int nylab = 7; 
  const double ylabels[nylab] = {0,1,0,1,0,1,2}; 
  for(int ii=0; ii<nylab; ii++){
    t->DrawTextNDC(0.13,0.083+0.1*ii, Form("%1.0f",ylabels[ii]));
  }
  
  /////////////////
  

  //area norm
  TH1F *data15_r_a, *data30_r_a, *data50_r_a; 
  data15_r_a = (TH1F)data15->Clone();
  data30_r_a = (TH1F)data30->Clone();
  data50_r_a = (TH1F)data50->Clone();
  TH1F *data15_toterr_r_a, *data30_toterr_r_a, *data50_toterr_r_a; 
  data15_toterr_r_a = (TH1F)data15_toterr->Clone();
  data30_toterr_r_a = (TH1F)data30_toterr->Clone();
  data50_toterr_r_a = (TH1F)data50_toterr->Clone();
  
  ////
  //divide and scale, axis
  if(!corrected) data15_r_a->GetYaxis()->SetTitle("Ratio to Pythia (uncorrected)"); 
  else data15_r_a->GetYaxis()->SetTitle("Ratio to Pythia"); 
  data15_r_a->Divide(mc15_a); 
  data30_r_a->Divide(mc30_a); 
  data50_r_a->Divide(mc50_a); 
  data15_toterr_r_a->Divide(mc15_a); 
  data30_toterr_r_a->Divide(mc30_a); 
  data50_toterr_r_a->Divide(mc50_a); 

  if(ratioPlotsInLog){
     data15_r_a->Scale(ratioFactorLog*ratioFactorLog); 
     data30_r_a->Scale(ratioFactorLog); 
     data15_toterr_r_a->Scale(ratioFactorLog*ratioFactorLog); 
     data30_toterr_r_a->Scale(ratioFactorLog); 
     data15_r_a->GetYaxis()->SetRangeUser(0.1,10000000);
  }
  else{
    //add offset
    data15_r_a->Add(ratioOffsetLin,2.0); 
    data30_r_a->Add(ratioOffsetLin); 
    data15_toterr_r_a->Add(ratioOffsetLin,2.0); 
    data30_toterr_r_a->Add(ratioOffsetLin); 
    data15_r_a->GetYaxis()->SetRangeUser(0,3*DratioOffsetLin+2);
  }
  ////

  std::cout << "AREANORM RATIO\n";
  /////////////////
  //AREANORM
  TCanvas *cucorrcomp2_1 = new TCanvas(("final_ratio_area_"+todraw).c_str(),"Ratio plots",700,300,500,500);
  if(ratioPlotsInLog) cucorrcomp2_1->SetLogy();
  cucorrcomp2_1->SetLeftMargin(0.14);
  data15_r_a->GetYaxis()->SetLabelSize(0.0);
  data15_r_a->Draw("E1");

//------------------------------------------------------------------------
//other generators
  //MADGRAPH
  if(bmadgraph && corrected){
    std::cout << "mg bins " << mg15_a->GetNbinsX() << std::endl;
    mg15_a->Divide(mc15_a);
    mg30_a->Divide(mc30_a);
    mg50_a->Divide(mc50_a);
    if(ratioPlotsInLog){
      mg15_a->Scale(ratioFactorLog*ratioFactorLog);
      mg30_a->Scale(ratioFactorLog);
    }
    else{
      mg15_a->Add(ratioOffsetLin,2.0); 
      mg30_a->Add(ratioOffsetLin); 
    }
    mg15_a->DrawCopy(mc_draw_option.c_str());
    mg30_a->DrawCopy(mc_draw_option.c_str());
    mg50_a->DrawCopy(mc_draw_option.c_str());
    mg15_a->SetFillColor(0);
    mg30_a->SetFillColor(0);
    mg50_a->SetFillColor(0);
    mg15_a->DrawCopy("Chistsame");
    mg30_a->DrawCopy("Chistsame");
    mg50_a->DrawCopy("Chistsame");
  }
  //MCATNLO
  if(bmcatnlo && corrected){
    mcat15_a->Divide(mc15_a);
    mcat30_a->Divide(mc30_a);
    mcat50_a->Divide(mc50_a);
    if(ratioPlotsInLog){
      mcat15_a->Scale(ratioFactorLog*ratioFactorLog);
      mcat30_a->Scale(ratioFactorLog);
    }
    else{
      mcat15_a->Add(ratioOffsetLin,2.0); 
      mcat30_a->Add(ratioOffsetLin); 
    }
    mcat15_a->DrawCopy(mc_draw_option.c_str());
    mcat30_a->DrawCopy(mc_draw_option.c_str());
    mcat50_a->DrawCopy(mc_draw_option.c_str());
    mcat15_a->SetFillColor(0);
    mcat30_a->SetFillColor(0);
    mcat50_a->SetFillColor(0);
    mcat15_a->DrawCopy("Chistsame");
    mcat30_a->DrawCopy("Chistsame");
    mcat50_a->DrawCopy("Chistsame");

  }
  //CASCADE
  if(bcascade && corrected){
//     std::cout << "CASC " << mcasc15_a->GetXaxis()->GetXmin() << " " << mcasc15_a->GetXaxis()->GetXmax() << " " << mcasc15_a->GetNbinsX() << std::endl;
//     std::cout << "PYTH " << mc15_a->GetXaxis()->GetXmin() << " " << mc15_a->GetXaxis()->GetXmax() << " " << mc15_a->GetNbinsX() << std::endl;
    mcasc15_a->Divide(mc15_a);
    mcasc30_a->Divide(mc30_a);
    mcasc50_a->Divide(mc50_a);
//     std::cout << "CASC " << mcasc15_a->GetXaxis()->GetXmin() << " " << mcasc15_a->GetXaxis()->GetXmax() << " " << mcasc15_a->GetNbinsX() << std::endl;
//     std::cout << "PYTH " << mc15_a->GetXaxis()->GetXmin() << " " << mc15_a->GetXaxis()->GetXmax() << " " << mc15_a->GetNbinsX() << std::endl;
    if(ratioPlotsInLog){
      mcasc15_a->Scale(ratioFactorLog*ratioFactorLog);
      mcasc30_a->Scale(ratioFactorLog);
    }
    else{
      mcasc15_a->Add(ratioOffsetLin,2.0); 
      mcasc30_a->Add(ratioOffsetLin); 
    }
    mcasc15_a->DrawCopy(mc_draw_option.c_str());
    mcasc30_a->DrawCopy(mc_draw_option.c_str());
    mcasc50_a->DrawCopy(mc_draw_option.c_str());
    mcasc15_a->SetFillColor(0);
    mcasc30_a->SetFillColor(0);
    mcasc50_a->SetFillColor(0);
    mcasc15_a->DrawCopy("Chistsame");
    mcasc30_a->DrawCopy("Chistsame");
    mcasc50_a->DrawCopy("Chistsame");

  }
//------------------------------------------------------------------------

  line->DrawCopy("Csame");
  line2->DrawCopy("Csame");
  line3->DrawCopy("Csame");
  data15_toterr_r_a->Draw("E1same");
  data30_toterr_r_a->Draw("E1same"); 
  data50_toterr_r_a->Draw("E1same"); 
  data15_r_a->Draw("E1same");
  data30_r_a->Draw("E1same"); 
  data50_r_a->Draw("E1same"); 
  leg2->Draw();
  leg3->Draw();
  tex2->Draw("same");

  //seplines:
  sepline->Draw("same");
  sepline2->Draw("same");


  //axis
  //ylabels: (look for data15_r->...->SetLabelSize(0.0); to have def. lab)
  TLatex *t = new TLatex();
  t->SetTextSize(0.035);//0.042 
//   t->SetTextFont(42); 
  t->SetTextAlign(31);
  const int nylab2 = 18; 
  const double ylabels2[nylab2] = {0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5}; 
  for(int ii=0; ii<nylab2; ii++){
    t->DrawTextNDC(0.13,0.086+0.04*ii, Form("%1.0f",ylabels2[ii]));
  }


  /////////////////


  //back-to-back norm
  TH1F *data15_r_b, *data30_r_b, *data50_r_b; 
  data15_r_b = (TH1F)data15->Clone();
  data30_r_b = (TH1F)data30->Clone();
  data50_r_b = (TH1F)data50->Clone();
  TH1F *data15_toterr_r_b, *data30_toterr_r_b, *data50_toterr_r_b; 
  data15_toterr_r_b = (TH1F)data15_toterr->Clone();
  data30_toterr_r_b = (TH1F)data30_toterr->Clone();
  data50_toterr_r_b = (TH1F)data50_toterr->Clone();
  
  ////
  //divide and scale, axis
  if(!corrected) data15_r_b->GetYaxis()->SetTitle("Ratio to Pythia (uncorrected)"); 
  else data15_r_b->GetYaxis()->SetTitle("Ratio to Pythia"); 
  std::cout << "## " << data15_r_b->GetBinCenter(firstBin) << " " << data15_r_b->GetBinCenter(nbins) << " " << data15_r_b->Integral(firstBin, nbins) << std::endl;
  data15_r_b->Divide(mc15_b); 
  data30_r_b->Divide(mc30_b); 
  data50_r_b->Divide(mc50_b); 
  data15_toterr_r_b->Divide(mc15_b); 
  data30_toterr_r_b->Divide(mc30_b); 
  data50_toterr_r_b->Divide(mc50_b); 

  if(ratioPlotsInLog){
    data15_r_b->Scale(ratioFactorLog*ratioFactorLog); 
    data30_r_b->Scale(ratioFactorLog); 
    data15_toterr_r_b->Scale(ratioFactorLog*ratioFactorLog); 
    data30_toterr_r_b->Scale(ratioFactorLog); 
    data15_r_b->GetYaxis()->SetRangeUser(0.1,10000000);
  }
  else{
    //add offset
    data15_r_b->Add(ratioOffsetLin,2.0); 
    data30_r_b->Add(ratioOffsetLin); 
    data15_toterr_r_b->Add(ratioOffsetLin,2.0); 
    data30_toterr_r_b->Add(ratioOffsetLin); 
    data15_r_b->GetYaxis()->SetRangeUser(0,3*DratioOffsetLin+2);
  }
  ////

  std::cout << "B2BNORM RATIO\n";
  /////////////////
  //B2BNORM
  TCanvas *cucorrcomp2_2 = new TCanvas(("final_ratio_btb_"+todraw).c_str(),"Ratio plots",700,500,500,500);
  if(ratioPlotsInLog) cucorrcomp2_2->SetLogy();
  cucorrcomp2_2->SetLeftMargin(0.14);
  data15_r_b->GetYaxis()->SetLabelSize(0.0);

  data15_r_b->Draw("E1");
  if(todraw=="dR"){
    if(ratioPlotsInLog) TBox *box2 = new TBox(2.4,0,4.4,1000000);
    else TBox *box2 = new TBox(2.4,0,4.4,3*DratioOffsetLin);
  }
  else if(todraw=="dPhi"){
    if(ratioPlotsInLog) TBox *box2 = new TBox(2.4,0,3.2,1000000);
    else TBox *box2 = new TBox(2.4,0,3.2,3*DratioOffsetLin);
  }
  box2->SetFillColor(boxcolor); 
  box2->SetLineColor(boxcolor); 
  box2->Draw();
//------------------------------------------------------------------------
//other generators
  //MADGRAPH
  if(bmadgraph && corrected){
    std::cout << "## " << mg15_b->GetBinCenter(firstBin) << " " << mg15_b->GetBinCenter(nbins) << " " << mg15_b->Integral(firstBin, nbins) << std::endl;
    mg15_b->Divide(mc15_b);
    mg30_b->Divide(mc30_b);
    mg50_b->Divide(mc50_b);
    if(ratioPlotsInLog){
      mg15_b->Scale(ratioFactorLog*ratioFactorLog);
      mg30_b->Scale(ratioFactorLog);
    }
    else{
      mg15_b->Add(ratioOffsetLin,2.0); 
      mg30_b->Add(ratioOffsetLin); 
    }

    mg15_b->DrawCopy(mc_draw_option.c_str());
    mg30_b->DrawCopy(mc_draw_option.c_str());
    mg50_b->DrawCopy(mc_draw_option.c_str());
    mg15_b->SetFillColor(0);
    mg30_b->SetFillColor(0);
    mg50_b->SetFillColor(0);
    mg15_b->DrawCopy("Chistsame");
    mg30_b->DrawCopy("Chistsame");
    mg50_b->DrawCopy("Chistsame");

  }
  //MCATNLO
  if(bmcatnlo && corrected){    
    std::cout << "## " << mcat15_b->GetBinCenter(firstBin) << " " << mcat15_b->GetBinCenter(nbins) << " " << mcat15_b->Integral(firstBin, nbins) << std::endl;
    mcat15_b->Divide(mc15_b);
    mcat30_b->Divide(mc30_b);
    mcat50_b->Divide(mc50_b);
    if(ratioPlotsInLog){
      mcat15_b->Scale(ratioFactorLog*ratioFactorLog);
      mcat30_b->Scale(ratioFactorLog);
    }
    else{
      mcat15_b->Add(ratioOffsetLin,2.0); 
      mcat30_b->Add(ratioOffsetLin); 
    }
    mcat15_b->DrawCopy(mc_draw_option.c_str());
    mcat30_b->DrawCopy(mc_draw_option.c_str());
    mcat50_b->DrawCopy(mc_draw_option.c_str());
    mcat15_b->SetFillColor(0);
    mcat30_b->SetFillColor(0);
    mcat50_b->SetFillColor(0);
    mcat15_b->DrawCopy("Chistsame");
    mcat30_b->DrawCopy("Chistsame");
    mcat50_b->DrawCopy("Chistsame");

  }
  //CASCADE
  if(bcascade && corrected){
    std::cout << "## " << mcasc15_b->GetBinCenter(firstBin) << " " << mcasc15_b->GetBinCenter(nbins) << " " << mcasc15_b->Integral(firstBin, nbins) << std::endl;
    mcasc15_b->Divide(mc15_b);
    mcasc30_b->Divide(mc30_b);
    mcasc50_b->Divide(mc50_b);
    if(ratioPlotsInLog){
      mcasc15_b->Scale(ratioFactorLog*ratioFactorLog);
      mcasc30_b->Scale(ratioFactorLog);
    }
    else{
      mcasc15_b->Add(ratioOffsetLin,2.0); 
      mcasc30_b->Add(ratioOffsetLin); 
    }
   
    mcasc15_b->DrawCopy(mc_draw_option.c_str());
    mcasc30_b->DrawCopy(mc_draw_option.c_str());
    mcasc50_b->DrawCopy(mc_draw_option.c_str());
    mcasc15_b->SetFillColor(0);
    mcasc30_b->SetFillColor(0);
    mcasc50_b->SetFillColor(0);
    mcasc15_b->DrawCopy("Chistsame");
    mcasc30_b->DrawCopy("Chistsame");
    mcasc50_b->DrawCopy("Chistsame");

  }
//------------------------------------------------------------------------

  line->DrawCopy("Csame");
  line2->DrawCopy("Csame");
  line3->DrawCopy("Csame");
  data15_toterr_r_b->Draw("E1same");
  data30_toterr_r_b->Draw("E1same"); 
  data50_toterr_r_b->Draw("E1same"); 
  data15_r_b->Draw("E1same");
  data30_r_b->Draw("E1same"); 
  data50_r_b->Draw("E1same"); 

  ////
  //vertical line
//   l->Draw("same"); 
  leg2->Draw();
  leg3->Draw();
  tex2->Draw();

  //seplines:
  sepline->Draw("same");
  sepline2->Draw("same");

  //axis
  for(int ii=0; ii<nylab2; ii++){
    t->DrawTextNDC(0.13,0.086+0.04*ii, Form("%1.0f",ylabels2[ii]));
  }


  //
  /////////////////


  
  
  //writing
  ffinal->cd();
  if(!corrected) {
    ffinal->mkdir("uncorrected");
    ffinal->cd("uncorrected");
  }
  else {
    ffinal->mkdir("corrected");
    ffinal->cd("corrected");
  }
  

  cucorrcomp1->Write();
  cucorrcomp1_1->Write();
  cucorrcomp1_2->Write();
  cucorrcomp2->Write();
  cucorrcomp2_1->Write();
  cucorrcomp2_2->Write();
 
  //printing
//     cucorrcomp1->Print(("./plotsFINAL/final_lumi_"+todraw+".pdf").c_str());
//     cucorrcomp1_1->Print(("./plotsFINAL/final_area_"+todraw+".pdf").c_str());
//     cucorrcomp1_2->Print(("./plotsFINAL/final_btb_"+todraw+".pdf").c_str());
//     cucorrcomp2->Print(("./plotsFINAL/final_ratio_lumi_"+todraw+".pdf").c_str());
//     cucorrcomp2_1->Print(("./plotsFINAL/final_ratio_area_"+todraw+".pdf").c_str());
//     cucorrcomp2_2->Print(("./plotsFINAL/final_ratio_btb_"+todraw+".pdf").c_str());

  ffinal->Close();
    
}


