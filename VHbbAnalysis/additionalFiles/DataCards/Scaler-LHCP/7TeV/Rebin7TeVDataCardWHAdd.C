#include "TH1.h"
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "math.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "THStack.h"
#include "TList.h"
#include "xsec7TeV.h"
#include "xsec8TeV.h"
#include "filesWH.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "tdrstyle.h"

const int oldBins = 50;
const int rebin = 2;

int kount;
int kount2;
TString IFILE;
int maxCount;


void makeSystPlot( TFile * f, TString oldFolder, RooWorkspace *WS,  string channel, string syst, int rightBinNo, int addRightBin)  
{

  RooArgList  * hobs = new RooArgList("hobs");
  RooRealVar BDT("CMS_vhbb_BDT_Wln", "BDT", -1, 1);///OLD VARIABLE NAME HERE
  hobs->add(*WS->var("CMS_vhbb_BDT_Wln"));  ///NEW VARIABLE NAME HERE
  RooWorkspace *tempWS =  (RooWorkspace*) f->Get(oldFolder.Data());
  TString systT(syst);
  TString chanT(channel);

  if((kount < 3) && (channel=="data_obs"))
  {
    kount++;
    std::string namen  = channel;
    
    std::cout << oldFolder.Data() << std::endl;
    std::cout << namen << std::endl;
    RooDataHist* tempRooDataHistNom = (RooDataHist*)  tempWS->data(namen.c_str());
    TH1 *tempHistNom = tempRooDataHistNom->createHistogram(namen.c_str(),BDT,RooFit::Binning(oldBins));

    tempHistNom->Rebin(rebin);

    if(addRightBin == 1)
    {
     float  err0 = tempHistNom->GetBinError(rightBinNo);
     float  con0 = tempHistNom->GetBinContent(rightBinNo);
 
     float  err1 = tempHistNom->GetBinError(rightBinNo-1);
     float  con1 = tempHistNom->GetBinContent(rightBinNo-1);
       
     tempHistNom->SetBinContent(rightBinNo,0);
     tempHistNom->SetBinError(rightBinNo,0);
     tempHistNom->SetBinContent(rightBinNo-1,con0+con1);
     tempHistNom->SetBinError(rightBinNo-1,sqrt(err0*err0+err1*err1));
    }



    RooDataHist *DHnom = new RooDataHist(channel.c_str(),"",*hobs,tempHistNom);  
    WS->import(*(new RooHistPdf(channel.c_str(),"",*hobs,*DHnom)));


 
 }

 if (channel!="data_obs")
{
  std::string nameUp; 
  std::string namen; 
  std::string nameDown;


  if(syst == "stat")
  {
   if(oldFolder.Contains("Wenu"))
   { 
     nameUp  = channel + "_CMS_vhbb_stat" + channel + "_WenuUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_stat" + channel + "_WenuDown";

     if(channel == "s_Top")
     {
     nameUp  = channel + "_CMS_vhbb_statsTop_WenuUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_statsTop_WenuDown";
     }

   }
   else
   {
     nameUp  = channel + "_CMS_vhbb_stat" + channel + "_WmunuUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_stat" + channel + "_WmunuDown";

     if(channel == "s_Top")
     {
     nameUp  = channel + "_CMS_vhbb_statsTop_WmunuUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_statsTop_WmunuDown";
     }


   }
  }
  else
  {
  nameUp  = channel + "_CMS_" + syst + "Up";
  namen  = channel;
  nameDown = channel + "_CMS_" + syst + "Down";
  }
 
  if((syst == "stat") && (oldFolder.Contains("2")))
  {
   if(oldFolder.Contains("Wenu"))
   { 
     nameUp  = channel + "_CMS_vhbb_stat" + channel + "_Wenu2Up";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_stat" + channel + "_Wenu2Down";

     if(channel == "s_Top")
     {
     nameUp  = channel + "_CMS_vhbb_statsTop_Wenu2Up";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_statsTop_Wenu2Down";
     }

   }
   else
   {
     nameUp  = channel + "_CMS_vhbb_stat" + channel +  "_Wmunu2Up";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_stat" + channel + "_Wmunu2Down";

     if(channel == "s_Top")
     {
     nameUp  = channel + "_CMS_vhbb_statsTop_Wmunu2Up";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_statsTop_Wmunu2Down";
     }


   }
  }

  if(systT.Contains("Model"))
  {
     nameUp  = channel + "_CMS_vhbb_WModelUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_WModelDown";
  }

  if( systT.Contains("stat") && (oldFolder.Contains("Wenu")) && IFILE.Contains("8TeV") && !(oldFolder.Contains("2")))
  { 
     nameUp  = channel + "_CMS_vhbb_stat" + channel + "_Wenu_8TeVUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_stat" + channel + "_Wenu_8TeVDown";

   if(channel == "s_Top")
     {
     nameUp  = channel + "_CMS_vhbb_statsTop_Wenu_8TeVUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_statsTop_Wenu_8TeVDown";
     }


  }

  if( systT.Contains("stat") && (oldFolder.Contains("Wmunu")) && IFILE.Contains("8TeV") && !(oldFolder.Contains("2")))
  { 
     nameUp  = channel + "_CMS_vhbb_stat" + channel + "_Wmnu_8TeVUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_stat" + channel + "_Wmnu_8TeVDown";

   if(channel == "s_Top")
     {
     nameUp  = channel + "_CMS_vhbb_statsTop_Wmnu_8TeVUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_statsTop_Wmnu_8TeVDown";
     }


  }


  if( systT.Contains("stat") && (oldFolder.Contains("Wenu")) && IFILE.Contains("8TeV") && (oldFolder.Contains("2")))
  { 
     nameUp  = channel + "_CMS_vhbb_stat" + channel + "_Wenu2_8TeVUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_stat" + channel + "_Wenu2_8TeVDown";

   if(channel == "s_Top")
     {
     nameUp  = channel + "_CMS_vhbb_statsTop_Wenu2_8TeVUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_statsTop_Wenu2_8TeVDown";
     }


  }

  if( systT.Contains("stat") && (oldFolder.Contains("Wmunu")) && IFILE.Contains("8TeV") && (oldFolder.Contains("2")))
  { 
     nameUp  = channel + "_CMS_vhbb_stat" + channel + "_Wmnu2_8TeVUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_stat" + channel + "_Wmnu2_8TeVDown";

   if(channel == "s_Top")
     {
     nameUp  = channel + "_CMS_vhbb_statsTop_Wmnu2_8TeVUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_statsTop_Wmnu2_8TeVDown";
     }


  }







  RooDataHist* tempRooDataHistUp = (RooDataHist*)  tempWS->data(nameUp.c_str());
  RooDataHist* tempRooDataHistDown = (RooDataHist*)  tempWS->data(nameDown.c_str());
  RooDataHist* tempRooDataHistNom = (RooDataHist*)  tempWS->data(namen.c_str());


  std::cout << oldFolder.Data() << std::endl; 
  std::cout << nameUp.c_str() << std::endl; 
  

  


  TH1 *tempHistUp = tempRooDataHistUp->createHistogram(nameUp.c_str(),BDT,RooFit::Binning(oldBins));
  TH1 *tempHistDown = tempRooDataHistDown->createHistogram(nameDown.c_str(),BDT,RooFit::Binning(oldBins));
  TH1 *tempHistNom = tempRooDataHistNom->createHistogram(namen.c_str(),BDT,RooFit::Binning(oldBins));
  
  
  std::cout<< "channel--> " << channel << std::endl;


    
  tempHistUp->SetLineColor(kRed);
  tempHistUp->SetLineWidth(3);
  tempHistUp->SetFillColor(0);
  
  tempHistDown->SetLineColor(kBlue);
  tempHistDown->SetFillColor(0);
  tempHistDown->SetLineWidth(3);
  
  tempHistNom->SetFillColor(0);
  tempHistNom->SetMarkerStyle(20);
    
  tempHistUp->SetTitle((channel + syst).c_str());


  tempHistNom->Rebin(rebin);
  tempHistUp->Rebin(rebin);
  tempHistDown->Rebin(rebin);
 
  if(addRightBin == 1)
    {
     float  err0 = tempHistNom->GetBinError(rightBinNo);
     float  con0 = tempHistNom->GetBinContent(rightBinNo);
 
     float  err1 = tempHistNom->GetBinError(rightBinNo-1);
     float  con1 = tempHistNom->GetBinContent(rightBinNo-1);
       
     tempHistNom->SetBinContent(rightBinNo,0);
     tempHistNom->SetBinError(rightBinNo,0);
     tempHistNom->SetBinContent(rightBinNo-1,con0+con1);
     tempHistNom->SetBinError(rightBinNo-1,sqrt(err0*err0+err1*err1));


     err0 = tempHistUp->GetBinError(rightBinNo);
     con0 = tempHistUp->GetBinContent(rightBinNo);
 
     err1 = tempHistUp->GetBinError(rightBinNo-1);
     con1 = tempHistUp->GetBinContent(rightBinNo-1);
       
     tempHistUp->SetBinContent(rightBinNo,0);
     tempHistUp->SetBinError(rightBinNo,0);
     tempHistUp->SetBinContent(rightBinNo-1,con0+con1);
     tempHistUp->SetBinError(rightBinNo-1,sqrt(err0*err0+err1*err1));



     err0 = tempHistDown->GetBinError(rightBinNo);
     con0 = tempHistDown->GetBinContent(rightBinNo);
 
     err1 = tempHistDown->GetBinError(rightBinNo-1);
     con1 = tempHistDown->GetBinContent(rightBinNo-1);
       
     tempHistDown->SetBinContent(rightBinNo,0);
     tempHistDown->SetBinError(rightBinNo,0);
     tempHistDown->SetBinContent(rightBinNo-1,con0+con1);
     tempHistDown->SetBinError(rightBinNo-1,sqrt(err0*err0+err1*err1));

    }

  RooDataHist *DHnom;
  RooDataHist *DHup = new RooDataHist(nameUp.c_str(),"",*hobs,tempHistUp);  
  if(kount2 < 3) DHnom = new RooDataHist(namen.c_str(),"",*hobs,tempHistNom);
  RooDataHist *DHdown = new RooDataHist(nameDown.c_str(),"",*hobs,tempHistDown);  

  WS->import(*(new RooHistPdf(nameUp.c_str(),"",*hobs,*DHup)));
  WS->import(*(new RooHistPdf(nameDown.c_str(),"",*hobs,*DHdown)));
  if(kount2 < 3){ WS->import(*(new RooHistPdf(namen.c_str(),"",*hobs,*DHnom))); kount2++;}

 }


}






void Process(TString fname, TString oldFolder)
{


  std::string channels[] = {"data_obs", "WH", "TT",  "WjLF", "WjHF", "ZjLF", "ZjHF" , "VV" , "s_Top"};
  std::string systs[] = {"eff_b", "fake_b", "res_j", "scale_j" , "stat" };
 
  TFile * file = new TFile(fname.Data(), "READ");
  std::cout << "reading " << fname.Data() << std::endl;
  fname.ReplaceAll(".root","_Rebin.root");


///BEGIN CHECK RBIN
  int addRightBin, rightBinNo;
  RooWorkspace *mytempWS =  (RooWorkspace*) file->Get(oldFolder.Data());
  RooRealVar BDT("CMS_vhbb_BDT_Wln", "BDT", -1, 1);

  float Signal = 0;
  float Background = 0;
  RooDataHist* tempRooDataHistNomS = (RooDataHist*)  mytempWS->data(channels[1].c_str());
  TH1 *tempHistNomS = tempRooDataHistNomS->createHistogram(channels[1].c_str(),BDT,RooFit::Binning(oldBins));
  tempHistNomS->Rebin(rebin);

  for(int i = 1; i < tempHistNomS->GetNbinsX(); i++)
  { 
    if(tempHistNomS->GetBinContent(i) > 0)
     { rightBinNo = i; Signal = tempHistNomS->GetBinContent(i);}
  }

  for(int i = 0; i < 9; i++)
  {
   RooDataHist* tempRooDataHistNom = (RooDataHist*)  mytempWS->data(channels[i].c_str());
   TH1 *tempHistNom = tempRooDataHistNom->createHistogram(channels[i].c_str(),BDT,RooFit::Binning(oldBins));
   tempHistNom->Rebin(rebin);
   Background += tempHistNom->GetBinContent(rightBinNo);
   
   float a;
   if(tempHistNom->GetBinContent(rightBinNo+1) > 0)
   {
   std::cout << "ARGHGHGHGHGH" << std::endl;
   std::cin >> a;
   }
  }

if(Background < 0.1) addRightBin = 1;
else  addRightBin = 0;
std::cout << "right bin no " << rightBinNo << " signal: " << Signal << " bkg: " << Background << std::endl;
std::cout << "REBINNING: " << addRightBin << std::endl;  
///END CHECK RBIN

  kount = 0;  

  gROOT->SetStyle("Plain");
  setTDRStyle();

  TFile * outfile = new TFile(fname.Data(), "RECREATE");

  using namespace RooFit;
  RooWorkspace *myWS = new RooWorkspace(oldFolder.Data(),oldFolder.Data());
  myWS->factory("CMS_vhbb_BDT_Wln[-1.,1.]"); ///NEW VARIABLE NAME HERE 

  TString oldFolder2(oldFolder.Data());
  oldFolder2.Append("2");
  RooWorkspace *myWS2 = new RooWorkspace(oldFolder2.Data(),oldFolder2.Data());
  myWS2->factory("CMS_vhbb_BDT_Wln[-1.,1.]"); ///NEW VARIABLE NAME HERE 


  
  for (int c =0; c<9; c++)
  {
     kount2 = 0;  
    for (int s =0; s<5 ; s++ ){
      makeSystPlot( file, oldFolder, myWS,  channels[c], systs[s], rightBinNo, addRightBin  );
      makeSystPlot( file, oldFolder2, myWS2,  channels[c], systs[s], rightBinNo, addRightBin );
    }
  }


  if(!(IFILE.Contains("8TeV")))
  {
  makeSystPlot(file, oldFolder, myWS, "WjLF", "WModel", rightBinNo, addRightBin);
  makeSystPlot(file, oldFolder, myWS, "WjHF", "WModel", rightBinNo, addRightBin);


  makeSystPlot(file, oldFolder2, myWS2, "WjLF", "WModel", rightBinNo, addRightBin);
  makeSystPlot(file, oldFolder2, myWS2, "WjHF", "WModel", rightBinNo, addRightBin);
  }

  myWS->writeToFile(fname.Data());  
  std::cout << std::endl << std::endl << std::endl << std::endl << "///////////////////////////" << std::endl;
  std::cout << fname.Data() << " written" << std::endl;
  std::cout << "///////////////////////////" << std::endl << std::endl << std::endl;


  outfile->Write();
  outfile->Close();
  fname.ReplaceAll("Rebin","2_Rebin");
  TFile * outfile2 = new TFile(fname.Data(), "RECREATE");
  myWS2->writeToFile(fname.Data());  
  std::cout << std::endl << std::endl << std::endl << std::endl << "///////////////////////////" << std::endl;
  std::cout << fname.Data() << " written" << std::endl;
  std::cout << "///////////////////////////" << std::endl << std::endl << std::endl;


}

  



void Rebin7TeVDataCardWHAdd()
{

maxCount=0;



int n = 134;

for(int i = 0; i < n; i++)
{
  TString oldFolder;
  IFILE = files[i];
  if(IFILE.Contains("Wenu")) oldFolder = "Wenu";
  if(IFILE.Contains("Wm"))oldFolder = "Wmunu";
  Process(IFILE, oldFolder);
  

  
  
}
}


