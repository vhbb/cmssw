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

const int bins = 25;
int kount;
int kount2;
TString IFILE;
int maxCount;

std::string massS[67]={
"110","110.5","111","111.5","112","112.5","113","113.5","114","114.5",
"115","115.5","116","116.5","117","117.5","118","118.5","119","119.5",
"120","120.5","121","121.5","122","122.5","123","123.5","124","124.5",
"125","125.5","126","126.5","127","127.5","128","128.5","129","129.5",
"130","130.5","131","131.5","132","132.5","133","133.5","134","134.5",
"135","124.6","124.7","124.8","124.9","125.1","125.2","125.3","125.4",
"125.6","125.7","125.8","125.9","126.1","126.2","126.3","126.4"
};



void makeSystPlot( TFile * f, TString myRooWS, RooWorkspace *WS,  string channel, string syst, int toMassNo, int fromMassNo) //massNo 0-51, see xSec7TeV.h 
{

  RooArgList  * hobs = new RooArgList("hobs");
  RooRealVar BDT("CMS_vhbb_BDT_Wln_8TeV", "CMS_vhbb_BDT_Wln_8TeV", -1, 1);///OLD VARIABLE NAME HERE
  hobs->add(*WS->var("CMS_vhbb_BDT_Wln_8TeV"));  ///NEW VARIABLE NAME HERE
  RooWorkspace *tempWS =  (RooWorkspace*) f->Get(myRooWS.Data());
  TString systT(syst);
  TString chanT(channel);

  bool writeIt = 1;
 
  if(chanT.Contains("QCD"))
  if(!(systT.Contains("stat"))) writeIt = 0;
  

  if((kount < 3) && (channel=="data_obs"))
  {
    kount++;
    std::string namen  = channel;
    
    std::cout << myRooWS.Data() << std::endl;
    std::cout << namen << std::endl;
    RooDataHist* tempRooDataHistNom = (RooDataHist*)  tempWS->data(namen.c_str());
    TH1 *tempHistNom = tempRooDataHistNom->createHistogram(namen.c_str(),BDT,Binning(bins));

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
   if(myRooWS.Contains("Wen"))
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
 
  if((syst == "stat") && (myRooWS.Contains("High")))
  {
   if(myRooWS.Contains("Wen"))
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

  if( systT.Contains("stat") && (myRooWS.Contains("Wen")) && IFILE.Contains("8TeV") && !(myRooWS.Contains("High")))
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

  if( systT.Contains("stat") && (myRooWS.Contains("Wmn")) && IFILE.Contains("8TeV") && !(myRooWS.Contains("High")))
  { 
     nameUp  = channel + "_CMS_vhbb_stat" + channel + "_Wmunu_8TeVUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_stat" + channel + "_Wmunu_8TeVDown";

   if(channel == "s_Top")
     {
     nameUp  = channel + "_CMS_vhbb_statsTop_Wmunu_8TeVUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_statsTop_Wmunu_8TeVDown";
     }


  }


  if( systT.Contains("stat") && (myRooWS.Contains("Wen")) && IFILE.Contains("8TeV") && (myRooWS.Contains("High")))
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

  if( systT.Contains("stat") && (myRooWS.Contains("Wmn")) && IFILE.Contains("8TeV") && (myRooWS.Contains("High")))
  { 
     nameUp  = channel + "_CMS_vhbb_stat" + channel + "_Wmunu2_8TeVUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_stat" + channel + "_Wmunu2_8TeVDown";

   if(channel == "s_Top")
     {
     nameUp  = channel + "_CMS_vhbb_statsTop_Wmunu2_8TeVUp";
     namen  = channel;
     nameDown  = channel + "_CMS_vhbb_statsTop_Wmunu2_8TeVDown";
     }


  }
if(writeIt)
{






  RooDataHist* tempRooDataHistUp = (RooDataHist*)  tempWS->data(nameUp.c_str());
  RooDataHist* tempRooDataHistDown = (RooDataHist*)  tempWS->data(nameDown.c_str());
  RooDataHist* tempRooDataHistNom = (RooDataHist*)  tempWS->data(namen.c_str());


  std::cout << myRooWS.Data() << std::endl; 
  std::cout << nameUp.c_str() << std::endl; 
  

  


  TH1 *tempHistUp = tempRooDataHistUp->createHistogram(nameUp.c_str(),BDT,Binning(bins));
  TH1 *tempHistDown = tempRooDataHistDown->createHistogram(nameDown.c_str(),BDT,Binning(bins));
  TH1 *tempHistNom = tempRooDataHistNom->createHistogram(namen.c_str(),BDT,Binning(bins));
  
  if(chanT.Contains("VH") && IFILE.Contains("7TeV"))
  {
     tempHistUp->Scale(xSec7WH[toMassNo]/xSec7WH[fromMassNo]);
     tempHistDown->Scale(xSec7WH[toMassNo]/xSec7WH[fromMassNo]);
     tempHistNom->Scale(xSec7WH[toMassNo]/xSec7WH[fromMassNo]);
  }  
 
  if(chanT.Contains("VH") && IFILE.Contains("8TeV"))
  {
     tempHistUp->Scale(xSec8WH[toMassNo]/xSec8WH[fromMassNo]);
     tempHistDown->Scale(xSec8WH[toMassNo]/xSec8WH[fromMassNo]);
     tempHistNom->Scale(xSec8WH[toMassNo]/xSec8WH[fromMassNo]);
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






void Process(TString fname, TString myRooWS, int toMass, int fromMass)
{




  std::string channels[] = {"data_obs", "VH", "TT",  "WjLF", "WjHF", "ZjLF", "ZjHF" , "VV" , "s_Top" };
  std::string systs[] = {"eff_b", "fake_b_8TeV", "res_j", "scale_j" , "stat" };

  kount = 0;  

  gROOT->SetStyle("Plain");
  setTDRStyle();

  TFile * file = new TFile(fname.Data(), "READ");
  std::cout << "reading " << fname.Data() << std::endl;
  TString outname(massS[toMass]);
  outname.Append(".root");
  fname.ReplaceAll(".root",outname.Data());
  TFile * outfile = new TFile(fname.Data(), "RECREATE");

  using namespace RooFit;
  RooWorkspace *myWS = new RooWorkspace(myRooWS.Data(),myRooWS.Data());
  myWS->factory("CMS_vhbb_BDT_Wln_8TeV[-1.,1.]"); ///NEW VARIABLE NAME HERE 

  
  for (int c =0; c<9; c++)
  {
     kount2 = 0;  
    for (int s =0; s<5 ; s++ ){
      makeSystPlot( file, myRooWS, myWS,  channels[c], systs[s], toMass, fromMass );
    }
  }


  if(!(IFILE.Contains("8TeV")))
  {
  makeSystPlot(file, myRooWS, myWS, "WjLF", "WModel",toMass, fromMass);
  makeSystPlot(file, myRooWS, myWS, "WjHF", "WModel",toMass, fromMass);
  }

  myWS->writeToFile(fname.Data());  
  std::cout << std::endl << std::endl << std::endl << std::endl << "///////////////////////////" << std::endl;
  std::cout << fname.Data() << " written" << std::endl;
  std::cout << "///////////////////////////" << std::endl << std::endl << std::endl;
  outfile->Write();
  outfile->Close();


}

  



void IntermediateMassMakerWH()
{

maxCount=0;




for(int i = 0; i < n; i++)
{
  TString myRooWS;
  IFILE = files[i];

if(IFILE.Contains("7TeV"))
{
  if(IFILE.Contains("Wenu")) myRooWS = "Wenu";
  if(IFILE.Contains("Wmn"))myRooWS = "Wmunu";
}

  if(IFILE.Contains ("Wmn") && IFILE.Contains("Low") && IFILE.Contains("8TeV"))  myRooWS = "WmnLowPt_8TeV";
  if(IFILE.Contains ("Wen") && IFILE.Contains("Low") && IFILE.Contains("8TeV"))  myRooWS = "WenLowPt_8TeV";
  if(IFILE.Contains ("Wmn") && IFILE.Contains("High") && IFILE.Contains("8TeV"))  myRooWS = "WmnHighPt_8TeV";
  if(IFILE.Contains ("Wen") && IFILE.Contains("High") && IFILE.Contains("8TeV"))  myRooWS = "WenHighPt_8TeV";

  if((IFILE.Contains("110")))
  {
     Process(IFILE, myRooWS, 0,0);
     Process(IFILE, myRooWS, 1,0);
     Process(IFILE, myRooWS, 2,0);
     Process(IFILE, myRooWS, 3,0);
     Process(IFILE, myRooWS, 4,0);
  }
  if((IFILE.Contains("115")))
  {
     for(int to = 5; to < 15; to++)   Process(IFILE, myRooWS, to , 10);
  }

  if((IFILE.Contains("120")))
  {
     for(int to = 15; to < 25; to++)   Process(IFILE, myRooWS,to , 20);
  }

  if((IFILE.Contains("125")))
  {
     for(int to = 25; to < 35; to++)   Process(IFILE, myRooWS,to , 30);
     for(int to = 51; to < 67; to++)   Process(IFILE, myRooWS,to , 30);
  }

  if((IFILE.Contains("130")))
  {
     for(int to = 35; to < 45; to++)   Process(IFILE, myRooWS,to , 40);
  }


  if((IFILE.Contains("135")))
  {
     for(int to = 45; to < 51; to++)   Process(IFILE, myRooWS,to , 50);
  }

  
  
}
}


