#include "vector.h"
#include "TFile.h"
#include <string>
#include <sstream>
#include <iostream>


double a[15] = {0.0,0.2,0.4,0.8,2.0,2.8,3.2,3.6,4.4,6.0}; 
int nbins = 9; 

TH1F *makeNonEqualBinning(TH1F *histo){

  std::cout << histo->GetName() << "\n";
//   double a[nbins] = {0.0,0.2,0.4,0.8,1.2,2.0,2.4,2.8,3.2,3.6,4.0,6.0}; 
  TH1F *temp = new TH1F("temp","temp",nbins,a); temp->Sumw2(); 
  temp->SetTitle(histo->GetTitle()); 
  
  for(unsigned int i=0; i<histo->GetNbinsX(); i++){
    int bin = temp->FindBin(histo->GetBinCenter(i));
    temp->SetBinContent(bin,temp->GetBinContent(bin)+histo->GetBinContent(i)); 
    temp->SetBinError(bin,sqrt(temp->GetBinError(bin)*temp->GetBinError(bin)+histo->GetBinError(i)*histo->GetBinError(i))); 

  }
  return temp; 
}

void addQCDandInclBBALLTP(){
  //string DIR = "/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/BBC/FIXLEOMACRO/CMSSW_3_7_0_patch2/src/UserCode/BbCorrelation/Macros/SUMFILES/MC/"; 

  string DIR = "/scratch/wehrlilu/BBCORR/MACROS/SUMFILES/MC/";


  string file15qcd = DIR+"histo-Pythia6-HLT_Jet15U-BBCorr_scaledEff-total_ETACUTOPEN.root"; 
  string file30qcd = DIR+"histo-Pythia6-HLT_Jet30U-BBCorr_scaledEff-total_ETACUTOPEN.root"; 
  string file50qcd = DIR+"histo-Pythia6-HLT_Jet50U-BBCorr_scaledEff-total_ETACUTOPEN.root"; 
  string file15inclbb = DIR+"histo-HLT_Jet15U-BBCorr_scaledEff-anV3-InclusiveBB_Pt30_Spring10-V8b_ETACUTOPEN.root"; 
  string file30inclbb = DIR+"histo-HLT_Jet30U-BBCorr_scaledEff-anV3-InclusiveBB_Pt30_Spring10-V8b_ETACUTOPEN.root"; 
  string file50inclbb = DIR+"histo-HLT_Jet50U-BBCorr_scaledEff-anV3-InclusiveBB_Pt30_Spring10-V8b_ETACUTOPEN.root"; 

  string strigger = "jet15"; addQCDandInclBB(file15qcd, file15inclbb, strigger, true); 
  strigger = "jet30"; addQCDandInclBB(file30qcd,file30inclbb,strigger, false); 
  strigger = "jet50"; addQCDandInclBB(file50qcd,file50inclbb,strigger, false); 
}

void addQCDandInclBB(string fileQCD="SUMFILES/MC/histo-Pythia6-HLT_Jet30U-BBCorr_scaledEff-total_ETACUTOPEN.root", string fileinclBB="SUMFILES/MC/histo-HLT_Jet30U-BBCorr_scaledEff-anV3-InclusiveBB_Pt30_Spring10-V8b_ETACUTOPEN.root", string strigger="jet30", bool recreate=true){

//   string strigger="jet30"; 

  TFile *fQCD = TFile::Open(fileQCD.c_str(),"READ");
  TFile *finclBB = TFile::Open(fileinclBB.c_str(),"READ");

  TH1F *t_efficnumBB, *t_efficdenBB, *t_puritynumBB, *t_puritydenBB; 
  finclBB->GetObject("Vert_dR_efficnum",t_efficnumBB);
  finclBB->GetObject("Vert_dR_efficden",t_efficdenBB);
  finclBB->GetObject("Vert_dR_puritynum",t_puritynumBB);
  finclBB->GetObject("Vert_dR_puritydenBB",t_puritydenBB);
//   //incl BB, norm to 1/pb: 
//   double neventsBB = 1017541/0.069; 
//   double xsecBB = 60411000; 
//   double lumiBB = neventsBB/xsecBB; 
//   t_efficnumBB->Scale(1./lumiBB); 
//   t_efficdenBB->Scale(1./lumiBB); 
//   t_puritynumBB->Scale(1./lumiBB); 
//   t_puritydenBB->Scale(1./lumiBB); 

  //make non equal binning
  TH1F *efficnumBB = makeNonEqualBinning(t_efficnumBB); 
  TH1F *efficdenBB = makeNonEqualBinning(t_efficdenBB); 
  TH1F *puritynumBB = makeNonEqualBinning(t_puritynumBB); 
  TH1F *puritydenBB = makeNonEqualBinning(t_puritydenBB); 

//   TCanvas *temp = new TCanvas("temp","tempBB",100,100,300,300); 
//   t_efficnumBB->Draw(); 
//   TCanvas *temp2 = new TCanvas("temp2","tempQCD",200,100,300,300); 
//   efficnumBB->SetLineColor(4);
//   efficnumBB->Draw();
//   t_efficnumBB->Draw("sames");

//   //test: 
//   std::cout << "B1/e B2/e " << t_efficnumBB->GetBinContent(1) << "/" << t_efficnumBB->GetBinError(1) << "    " << t_efficnumBB->GetBinContent(2) << "/" 
// 	    << t_efficnumBB->GetBinError(2) << std::endl; 
//   t_efficnumBB->Rebin(2); 
//   std::cout << "Rebin2 -> " << t_efficnumBB->GetBinContent(1) << "/" << t_efficnumBB->GetBinError(1) << std::endl; 

  

  TH1F *t_efficnumQCD, *t_efficdenQCD, *t_puritynumQCD, *t_puritydenQCDBB, *t_puritydenQCDNONBB; 
  fQCD->GetObject("Vert_dR_efficnum",t_efficnumQCD);
  fQCD->GetObject("Vert_dR_efficden",t_efficdenQCD);
  fQCD->GetObject("Vert_dR_puritynum",t_puritynumQCD);
  fQCD->GetObject("Vert_dR_puritydenBB",t_puritydenQCDBB);
  fQCD->GetObject("Vert_dR_puritydenNONBB",t_puritydenQCDNONBB);

  //make non equal binning
  TH1F *efficnumQCD = makeNonEqualBinning(t_efficnumQCD); 
  TH1F *efficdenQCD = makeNonEqualBinning(t_efficdenQCD); 
  TH1F *puritynumQCD = makeNonEqualBinning(t_puritynumQCD); 
  TH1F *puritydenQCDBB = makeNonEqualBinning(t_puritydenQCDBB); 
  TH1F *puritydenQCDNONBB = makeNonEqualBinning(t_puritydenQCDNONBB); 



  //output file
  TFile *fout; 
  if(recreate) fout = (TFile*)TFile::Open("efficCorrection.root","RECREATE"); 
  else fout = (TFile*)TFile::Open("efficCorrection.root","UPDATE"); 
  fout->cd();
  fout->mkdir(strigger.c_str());
  fout->cd(strigger.c_str());


//   TCanvas *temp = new TCanvas("temp","tempBB",100,100,300,300); 
//   efficnumBB->Draw(); 
//   TCanvas *temp2 = new TCanvas("temp2","tempQCD",200,100,300,300); 
//   efficnumQCD->Draw();

  //adding the histos for all apart purityden
  TH1F *efficnumEND = efficnumQCD->Clone(); 
  TH1F *efficdenEND = efficdenQCD->Clone(); 
  TH1F *puritynumEND = puritynumQCD->Clone(); 
  TH1F *puritydenEND = puritydenQCDBB->Clone(); 


  double inv_wBB, inv_wQCD, wBB, wQCD, binRes, errRes;

  //efficnum 
  for(unsigned int i=1; i<efficnumQCD->GetNbinsX(); i++){
    inv_wBB = efficnumBB->GetBinError(i)*efficnumBB->GetBinError(i); 
    inv_wQCD = efficnumQCD->GetBinError(i)*efficnumQCD->GetBinError(i); 
    if(inv_wBB==0) wBB = 0.0; 
    else wBB = 1./inv_wBB; 
    if(inv_wQCD==0) wQCD = 0.0; 
    else wQCD = 1./inv_wQCD; 
    if(inv_wBB+inv_wQCD==0){
      binRes = 0; 
      if(efficnumBB->GetBinContent(i)!=0 || efficnumQCD->GetBinContent(i)!=0) std::cout << "STRANGE\n";
    }
    else binRes = (efficnumBB->GetBinContent(i)*wBB + efficnumQCD->GetBinContent(i)*wQCD)/(wBB+wQCD);
    if(inv_wBB+inv_wQCD==0) errRes = 0; 
    else errRes = 1./sqrt(wBB+wQCD); 
    efficnumEND->SetBinContent(i,binRes); 
    efficnumEND->SetBinError(i,errRes); 
  }

  //efficden 
  for(unsigned int i=1; i<efficdenQCD->GetNbinsX(); i++){
    inv_wBB = efficdenBB->GetBinError(i)*efficdenBB->GetBinError(i); 
    inv_wQCD = efficdenQCD->GetBinError(i)*efficdenQCD->GetBinError(i); 
    if(inv_wBB==0) wBB = 0.0; 
    else wBB = 1./inv_wBB; 
    if(inv_wQCD==0) wQCD = 0.0; 
    else wQCD = 1./inv_wQCD; 
    if(inv_wBB+inv_wQCD==0){
      binRes = 0; 
      if(efficdenBB->GetBinContent(i)!=0 || efficdenQCD->GetBinContent(i)!=0) std::cout << "STRANGE\n";
    }
    else binRes = (efficdenBB->GetBinContent(i)*wBB + efficdenQCD->GetBinContent(i)*wQCD)/(wBB+wQCD);
    if(inv_wBB+inv_wQCD==0) errRes = 0; 
    else errRes = 1./sqrt(wBB+wQCD); 
    efficdenEND->SetBinContent(i,binRes); 
    efficdenEND->SetBinError(i,errRes); 
  }

  //puritynum 
  for(unsigned int i=1; i<puritynumQCD->GetNbinsX(); i++){
    inv_wBB = puritynumBB->GetBinError(i)*puritynumBB->GetBinError(i); 
    inv_wQCD = puritynumQCD->GetBinError(i)*puritynumQCD->GetBinError(i); 
    if(inv_wBB==0) wBB = 0.0; 
    else wBB = 1./inv_wBB; 
    if(inv_wQCD==0) wQCD = 0.0; 
    else wQCD = 1./inv_wQCD; 
    if(inv_wBB+inv_wQCD==0){
      binRes = 0; 
      if(puritynumBB->GetBinContent(i)!=0 || puritynumQCD->GetBinContent(i)!=0) std::cout << "STRANGE\n";
    }
    else binRes = (puritynumBB->GetBinContent(i)*wBB + puritynumQCD->GetBinContent(i)*wQCD)/(wBB+wQCD);
    if(inv_wBB+inv_wQCD==0) errRes = 0; 
    else errRes = 1./sqrt(wBB+wQCD); 
    puritynumEND->SetBinContent(i,binRes); 
    puritynumEND->SetBinError(i,errRes); 
  }

  //purityden 
  for(unsigned int i=1; i<puritydenQCDBB->GetNbinsX(); i++){
    inv_wBB = puritydenBB->GetBinError(i)*puritydenBB->GetBinError(i); 
    inv_wQCD = puritydenQCDBB->GetBinError(i)*puritydenQCDBB->GetBinError(i); 
    if(inv_wBB==0) wBB = 0.0; 
    else wBB = 1./inv_wBB; 
    if(inv_wQCD==0) wQCD = 0.0; 
    else wQCD = 1./inv_wQCD; 
    if(inv_wBB+inv_wQCD==0){
      binRes = 0; 
      if(puritydenBB->GetBinContent(i)!=0 || puritydenQCDBB->GetBinContent(i)!=0) std::cout << "STRANGE\n";
    }
    else binRes = (puritydenBB->GetBinContent(i)*wBB + puritydenQCDBB->GetBinContent(i)*wQCD)/(wBB+wQCD);
    if(inv_wBB+inv_wQCD==0) errRes = 0; 
    else errRes = 1./sqrt(wBB+wQCD); 
    puritydenEND->SetBinContent(i,binRes); 
    puritydenEND->SetBinError(i,errRes); 
  }


  puritydenEND->Add(puritydenQCDNONBB); 

  //from andrea mail
// float inv_w1 = h1->GetBinError(i)*h1->GetBinError(i);
// float inv_w2 = h2->GetBinError(i)*h2->GetBinError(i);
// float w1 = 1./inv_w1;
// float w2 = 1./inv_w2;
// float bin_i = (h1->GetBinContent(i)*w1 +
// h1->GetBinContent(i)*w2) / (w1+w2);
// float err_i = 1./sqrt(w1+w2)
// hres->SetBinContent(bin_i);
// hres->SetBinError(err_i);


  //compute actual corrections

  TCanvas *ce = new TCanvas("efficiency_corr","efficiency correction curve",100,100,600,400);
  TCanvas *cp = new TCanvas("purity_corr","purity correction curve",100,100,600,400);

//   double a[12] = {0.0,0.2,0.4,0.8,1.2,2.0,2.4,2.8,3.2,3.6,4.0,6.0}; 

  TH1F *effic = new TH1F("effic","effic",nbins,a); effic->Sumw2();
  TH1F *purity = new TH1F("purity","purity",nbins,a); purity->Sumw2();
  TH1F *corrfunc = new TH1F("corrfunc","corrfunc",nbins,a); corrfunc->Sumw2();

  TH1F * effic1 = new TH1F("efficbb","effic",nbins,a);
  TH1F * purity1 = new TH1F("puritybb","purity",nbins,a);
  TH1F * corrfunc1 = new TH1F("corrfuncbb","corrfunc",nbins,a);
  TH1F * effic2 = new TH1F("efficqcd","effic",nbins,a);
  TH1F * purity2 = new TH1F("purityqcd","purity",nbins,a);
  TH1F * corrfunc2 = new TH1F("corrfuncqcd","corrfunc",nbins,a);
  
  effic1->Sumw2();
  purity1->Sumw2();
  corrfunc1->Sumw2();
  
  effic2->Sumw2();
  purity2->Sumw2();
  corrfunc2->Sumw2();




  effic->Divide(efficnumEND,efficdenEND,1,1,"B");

  //new computation of purity...
  purity->Divide(puritynumEND,puritydenEND,1,1,"B");
  corrfunc->Divide(purity,effic,1,1);

  //other effic/purity of single samples...
  effic1->Divide(efficnumBB,efficdenBB,1,1,"B");
  purity1->Divide(puritynumBB,puritydenBB,1,1,"B");
  corrfunc1->Divide(purity1,effic1,1,1);

  effic2->Divide(efficnumQCD,efficdenQCD,1,1,"B");
  TH1F *puritydenQCD = (TH1F*)puritydenQCDBB->Clone(); 
  puritydenQCD->Add(puritydenQCDNONBB); 
  purity2->Divide(puritynumQCD,puritydenQCD,1,1,"B");
  corrfunc2->Divide(purity2,effic2,1,1);




  double eff = 0, pur = 0;
  if(efficdenEND->Integral()>0) eff = efficnumEND->Integral()/efficdenEND->Integral(); 
  if(puritydenEND->Integral()>0) pur = puritynumEND->Integral()/puritydenEND->Integral(); 
  std::cout << "Effic " << efficnumEND->Integral() << "/" << efficdenEND->Integral() << " = " 
	    << eff << std::endl;
  std::cout << "Purity " << puritynumEND->Integral() << "/" << puritydenEND->Integral() << " = " 
	    << pur << std::endl;

  std::stringstream spur; spur << (int)(pur*100.0); 
  std::stringstream seff; seff << (double)((int)(eff*1000.0))/10.0; 

  effic->GetYaxis()->SetRangeUser(0,0.2);
  effic->SetLineColor(4);
  effic->SetMarkerColor(4);
  effic->SetMarkerStyle(21);
  purity->GetYaxis()->SetRangeUser(0,1.1);
  purity->SetLineColor(8);
  purity->SetMarkerColor(8);
  purity->SetMarkerStyle(21);
  ce->cd();
  effic->Draw("PE");
  double maxe = effic->GetMaximum();
  TPaveLabel *ple = new TPaveLabel(0.5,maxe-0.05,1.5,maxe-0.01,(seff.str()+"%").c_str());
  ple->SetFillColor(0);
  ple->SetTextColor(4);
  ple->SetBorderSize(0); 
  ple->Draw();
  cp->cd();
  purity->Draw("PE");
  TPaveLabel *plp = new TPaveLabel(1,0.3,3,0.5,(spur.str() + "%").c_str());
  plp->SetFillColor(0);
  plp->SetTextColor(8);
  plp->SetBorderSize(0); 
  plp->Draw();



  fout->cd(strigger.c_str());
  effic->Write();
  purity->Write();
  corrfunc->Write();
  effic1->Write();
  purity1->Write();
  corrfunc1->Write();
  effic2->Write();
  purity2->Write();
  corrfunc2->Write();
  ce->Write();
  cp->Write();
  fout->Write();
  fout->Close();


}

void plotBBvsQCDvsCombined(string histo="corrfunc", string strigger="jet30"){
  gStyle->SetOptStat(0); 
  fout = (TFile*)TFile::Open(("efficCorrection" + strigger + ".root").c_str(),"READ"); 
  TH1F *hcomb, *hbb, *hqcd; 
  fout->GetObject((strigger+"/"+histo).c_str(),hcomb);
  fout->GetObject((strigger+"/"+histo+"bb").c_str(),hbb);
  fout->GetObject((strigger+"/"+histo+"qcd").c_str(),hqcd);
  string striggerFull = "HLT_Jet30U";
  if(strigger == "jet15") striggerFull = "HLT_Jet15U";
  if(strigger == "jet50") striggerFull = "HLT_Jet50U";

  hcomb->SetLineColor(1); 
  hcomb->SetMarkerColor(1); 
  hcomb->SetMarkerStyle(1); 

  hcomb->SetLineWidth(2); 
  hcomb->GetXaxis()->SetTitle("#Delta R");
  hcomb->SetTitle("");
  hqcd->SetFillColor(2); 
  hqcd->SetFillStyle(3004); 
  hbb->SetFillColor(4); 
  hbb->SetFillStyle(3005); 
  hbb->SetLineColor(0); 
  hqcd->SetLineColor(0); 

  TH1F *hcomb_err13 = hcomb->Clone(); 
  int color = kBlue-7; 
  if(histo=="corrfunc") color = kYellow-7; 
  if(histo=="purity") color = kGreen-7;

  hcomb_err13->SetFillColor(color); 
  hcomb_err13->SetLineColor(color); 
  for(unsigned int i=1; i<10; i++){
    hcomb_err13->SetBinError(i,0.13*hcomb_err13->GetBinContent(i)); 
  }
  std::cout << "test " << hcomb_err13->GetBinError(1)/hcomb_err13->GetBinContent(1) << "\n";
  TCanvas *cp = new TCanvas("corr","e/p/correction curve",100,100,600,600);
  hcomb->GetXaxis()->SetRangeUser(0,3.99); 
		
  hcomb->Draw(); 
  if(histo=="corrfunc") hcomb_err13->Draw("E3same"); 
  hcomb->Draw("same"); 

//   hqcd->Draw("E2same");
//   hbb->Draw("E2same"); 

  TLegend *leg = new TLegend(0.15,0.7,0.4,0.9);
  string header; int lineC = 1;
  if(histo=="corrfunc"){
    header="Correction function";
  }
  if(histo=="effic") {
    header="Efficiency";
    lineC = 4; 
  }
  if(histo=="purity") {
    header="Purity"; 
    lineC = kGreen+2; 
  }
  hcomb->SetLineColor(lineC); 
  leg->SetHeader(header.c_str()); 
  leg->SetFillColor(0);
  leg->AddEntry(hqcd,"QCD","F"); 
  leg->AddEntry(hbb,"Incl. BB","F");
  leg->AddEntry(hcomb,"Combined"); 
//   leg->Draw();


    //cms label 
  double axismax; 

  double axismax  = hcomb->GetMaximum(); 
  if(histo!="corrfunc") 
    {
      TLatex * tex = new TLatex(0.2,axismax-axismax/20.0,header.c_str());
      TLatex * tex2 = new TLatex(2.6,axismax-axismax/20.0,striggerFull.c_str());
    }
  else {
    TLatex * tex = new TLatex(2.6,axismax+axismax/20.0,header.c_str());
    TLatex * tex2 = new TLatex(2.6,axismax-axismax/20.0,striggerFull.c_str());
  }
  tex->SetTextSize(0.040);
  tex->SetLineWidth(2);
  tex->SetTextColor(lineC); 
  //std::cout << "linec "<< lineC << std::endl;
  tex->Draw();
  tex2->SetTextSize(0.040);
  tex2->SetLineWidth(2);
  tex2->Draw();



}
