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
#include "filesZll.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "tdrstyle.h"

const int bins = 10;
int kount;
int kount2;
TString IFILE;
int maxCount;
std::string massS[51]={
"110","110_5","111","111_5","112","112_5","113","113_5","114","114_5",
"115","115_5","116","116_5","117","117_5","118","118_5","119","119_5",
"120","120_5","121","121_5","122","122_5","123","123_5","124","124_5",
"125","125_5","126","126_5","127","127_5","128","128_5","129","129_5",
"130","130_5","131","131_5","132","132_5","133","133_5","134","134_5",
"135"};


void makeSystPlot( TFile * f, TString oldFolder, RooWorkspace *WS,  string channel, string syst, int toMassNo, int fromMassNo ,int rightBinNo, int addRightBin) //massNo 0-51, see xSec7TeV.h 
{

  RooArgList  * hobs = new RooArgList("hobs");
  RooRealVar BDT("CMS_vhbb_BDT_Zll", "CMS_vhbb_BDT_Zll", -1, 1);///OLD VARIABLE NAME HERE
  hobs->add(*WS->var("CMS_vhbb_BDT_Zll"));  ///NEW VARIABLE NAME HERE
  RooWorkspace *tempWS =  (RooWorkspace*) f->Get(oldFolder.Data());
  TString systT(syst);
  TString chanT(channel);
  bool writeIt = 1;
 
  if(chanT.Contains("QCD") || chanT.Contains("Wj"))
  if(!(systT.Contains("stat"))) writeIt = 0;
 

  if((kount < 3) && (channel=="data_obs"))
  {
    kount++;
    std::string namen  = channel;
    
    std::cout << "reading WS "<< oldFolder.Data() << std::endl;
    std::cout << namen << std::endl;
    RooDataHist* tempRooDataHistNom = (RooDataHist*)  tempWS->data(namen.c_str());
    TH1 *tempHistNom = tempRooDataHistNom->createHistogram(namen.c_str(),BDT,RooFit::Binning(bins));
    std::cout << namen << std::endl;


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


  if((syst == "stat"))
  {
    if(IFILE.Contains("7TeV"))
    {
      nameUp  = channel + "CMS_vhbb_stats_" + channel + "_" + oldFolder.Data() + "Up";
      namen  = channel;
      nameDown  = channel + "CMS_vhbb_stats_" + channel + "_" + oldFolder.Data() + "Down";
    }
    if(IFILE.Contains("8TeV"))
    { 
      nameUp  = channel + "CMS_vhbb_stats_" + channel + "_" + oldFolder.Data() + "Up";
      namen  = channel;
      nameDown  = channel + "CMS_vhbb_stats_" + channel + "_" + oldFolder.Data() + "Down";
    }

  }
  else
  {
  nameUp  = channel + "CMS_" + syst + "Up";
  namen  = channel;
  nameDown = channel + "CMS_" + syst + "Down";
  }
  if((syst == "ZJModel"))
  {
    if(IFILE.Contains("7TeV"))
    { 
      nameUp  = channel + "CMS_vhbb_ZJModel_" + oldFolder.Data() + "_7TeVUp";
      namen  = channel;
      nameDown  = channel + "CMS_vhbb_ZJModel_" + oldFolder.Data() + "_7TeVDown";
    }
    if(IFILE.Contains("8TeV"))
    { 
      nameUp  = channel + "CMS_vhbb_ZJModel_" + oldFolder.Data() + "_8TeVUp";
      namen  = channel;
      nameDown  = channel + "CMS_vhbb_ZJModel_" + oldFolder.Data() + "_8TeVDown";
    }

  }
  
  if(writeIt)
   {
    RooDataHist* tempRooDataHistUp = (RooDataHist*)  tempWS->data(nameUp.c_str());
    RooDataHist* tempRooDataHistDown = (RooDataHist*)  tempWS->data(nameDown.c_str());
    RooDataHist* tempRooDataHistNom = (RooDataHist*)  tempWS->data(namen.c_str());


   
    std::cout << oldFolder.Data() << std::endl; 
    std::cout << nameUp.c_str() << std::endl; 
 	 
    TH1 *tempHistUp = tempRooDataHistUp->createHistogram(nameUp.c_str(),BDT,RooFit::Binning(bins));
    TH1 *tempHistDown = tempRooDataHistDown->createHistogram(nameDown.c_str(),BDT,RooFit::Binning(bins));
    std::cout << namen.c_str() << std::endl;
    TH1 *tempHistNom = tempRooDataHistNom->createHistogram(namen.c_str(),BDT,RooFit::Binning(bins));

    if(chanT.Contains("ZH") && IFILE.Contains("7TeV"))
    {
     tempHistUp->Scale(xSec7ZH[toMassNo]/xSec7ZH[fromMassNo]);
     tempHistDown->Scale(xSec7ZH[toMassNo]/xSec7ZH[fromMassNo]);
     tempHistNom->Scale(xSec7ZH[toMassNo]/xSec7ZH[fromMassNo]);
    }  
 
   if(chanT.Contains("ZH") && IFILE.Contains("8TeV"))
   {
     tempHistUp->Scale(xSec8ZH[toMassNo]/xSec8ZH[fromMassNo]);
     tempHistDown->Scale(xSec8ZH[toMassNo]/xSec8ZH[fromMassNo]);
     tempHistNom->Scale(xSec8ZH[toMassNo]/xSec8ZH[fromMassNo]);
   }  
   
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


}






void Process(TString fname, TString oldFolder, int toMass, int fromMass)
{




  std::string channels[] = {"data_obs", "ZH", "TT",  "WjLF", "WjHF", "ZjLF", "ZjHF" , "VV" , "s_Top", "QCD"};
  std::string systs[] = {"eff_b", "fake_b", "res_j", "scale_j" , "stat" };

  kount = 0;  

  gROOT->SetStyle("Plain");
  setTDRStyle();

  TFile * file = new TFile(fname.Data(), "READ");
  std::cout << "reading " << fname.Data() << std::endl;
  TString outname(massS[toMass]);
  outname.Append(".root");
  fname.ReplaceAll(".root",outname.Data());
  

///BEGIN CHECK RBIN
  int addRightBin, rightBinNo;
  float a,overflow;
  RooWorkspace *mytempWS =  (RooWorkspace*) file->Get(oldFolder.Data());
  RooRealVar BDT("CMS_vhbb_BDT_Zll", "BDT", -1, 1);

  float Signal = 0;
  float Background = 0;
  RooDataHist* tempRooDataHistNomS = (RooDataHist*)  mytempWS->data(channels[1].c_str());
  TH1 *tempHistNomS = tempRooDataHistNomS->createHistogram(channels[1].c_str(),BDT,RooFit::Binning(bins));

  for(int i = 1; i <= tempHistNomS->GetNbinsX(); i++)
  { 
    if(tempHistNomS->GetBinContent(i) > 0)
     { rightBinNo = i; Signal = tempHistNomS->GetBinContent(i);}
  }

  for(int i = 2; i < 9; i++)
  {
   RooDataHist* tempRooDataHistNom = (RooDataHist*)  mytempWS->data(channels[i].c_str());
   TH1 *tempHistNom = tempRooDataHistNom->createHistogram(channels[i].c_str(),BDT,RooFit::Binning(bins));
   Background += tempHistNom->GetBinContent(rightBinNo);
   
   overflow = tempHistNom->GetBinContent(rightBinNo+1) ;
   
   if(tempHistNom->GetBinContent(rightBinNo+1) > 0)
   {
   std::cout << "ARGHGHGHGHGH" << std::endl;
   //   std::cin >> ;
   }
   a+= overflow;
  }

  if( (Background ==0) || ( (Signal/Background) > 3. )  ) addRightBin = 1;
else  addRightBin = 0;
std::cout << "################# folder" << oldFolder << std::endl;
 std::cout << "################# CHECK RBIN:: right bin n  " << rightBinNo << " signal: " << Signal << " bkg: " << Background <<  " at right there is an overflow of: "<< a << std::endl;
std::cout << "########################### CHECK RBIN:: REBINNING: " << addRightBin << std::endl;  
///END CHECK RBIN




  TFile * outfile = new TFile(fname.Data(), "RECREATE");

  using namespace RooFit;
  RooWorkspace *myWS = new RooWorkspace(oldFolder.Data(),oldFolder.Data());
  myWS->factory("CMS_vhbb_BDT_Zll[-1.,1.]"); ///NEW VARIABLE NAME HERE 

  
  
  for (int c =0; c<10; c++)
  {
    kount2 = 0;  
    for (int s =0; s<5 ; s++ ){
      makeSystPlot( file, oldFolder, myWS,  channels[c], systs[s], toMass, fromMass ,rightBinNo, addRightBin );
    }
  }


  



  myWS->writeToFile(fname.Data());  
  std::cout << std::endl << std::endl << std::endl << std::endl << "///////////////////////////" << std::endl;
  std::cout << fname.Data() << " written" << std::endl;
  std::cout << "///////////////////////////" << std::endl << std::endl << std::endl;


  outfile->Write();
  outfile->Close();


}

  



void Rebin7TeVDataCardZllAdd()
{

maxCount=0;




for(int i = 0; i < n; i++)
{
  TString oldFolder;
  IFILE = files[i];

  if(IFILE.Contains ("Zee1") && IFILE.Contains("7TeV"))  oldFolder = "ZeeLoose_7TeV";
  if(IFILE.Contains ("Zmm1") && IFILE.Contains("7TeV"))  oldFolder = "ZmmLoose7TeV";
  if(IFILE.Contains ("Zee2") && IFILE.Contains("7TeV"))  oldFolder = "ZeeTight_7TeV";
  if(IFILE.Contains ("Zmm2") && IFILE.Contains("7TeV"))  oldFolder = "ZmmTight_7TeV";

  if(IFILE.Contains ("Zee") && IFILE.Contains("1") && IFILE.Contains("8TeV"))  oldFolder = "ZeeLoose_8TeV";
  if(IFILE.Contains ("Zmm") && IFILE.Contains("1") && IFILE.Contains("8TeV"))  oldFolder = "ZmmLoose_8TeV";
  if(IFILE.Contains ("Zee") && IFILE.Contains("2") && IFILE.Contains("8TeV"))  oldFolder = "ZeeTight_8TeV";
  if(IFILE.Contains ("Zmm") && IFILE.Contains("2") && IFILE.Contains("8TeV"))  oldFolder = "ZmmTight_8TeV";


  if((IFILE.Contains("110")))
  {
     Process(IFILE, oldFolder, 0,0);
     Process(IFILE, oldFolder, 1,0);
     Process(IFILE, oldFolder, 2,0);
     Process(IFILE, oldFolder, 3,0);
     Process(IFILE, oldFolder, 4,0);
  }
  if((IFILE.Contains("115")))
  {
     for(int to = 5; to < 15; to++)   Process(IFILE, oldFolder, to , 10);
  }

  if((IFILE.Contains("120")))
  {
     for(int to = 15; to < 25; to++)   Process(IFILE, oldFolder,to , 20);
  }

  if((IFILE.Contains("125")))
  {
     for(int to = 25; to < 35; to++)   Process(IFILE, oldFolder,to , 30);
  }

  if((IFILE.Contains("130")))
  {
     for(int to = 35; to < 45; to++)   Process(IFILE, oldFolder,to , 40);
  }


  if((IFILE.Contains("135")))
  {
     for(int to = 45; to < 51; to++)   Process(IFILE, oldFolder,to , 50);
  }

  
  
}
}


