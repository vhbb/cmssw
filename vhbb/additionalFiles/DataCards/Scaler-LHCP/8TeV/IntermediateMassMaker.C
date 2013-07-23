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
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "SYSTS.h"

int bins = 40; //mBDT
TString IFILE;
RooDataHist* tempRooDataHistNom;
TH1 *tempHistNom;

std::string massS[97]={
"110","110.5","111","111.5","112","112.5","113","113.5","114","114.5",
"115","115.5","116","116.5","117","117.5","118","118.5","119","119.5",
"120","120.5","121","121.5","122","122.5","123","123.5","124","124.5",
"125","125.5","126","126.5","127","127.5","128","128.5","129","129.5",
"130","130.5","131","131.5","132","132.5","133","133.5","134","134.5",
"135","135.5","136","136.5","137","137.5","138","138.5","139","139.5",
"140","140.5","141","141.5","142","142.5","143","143.5","144","144.5",
"145","145.5","146","146.5","147","147.5","148","148.5","149","149.5",
"150",

"124.6","124.7","124.8","124.9","125.1","125.2","125.3","125.4",
"125.6","125.7","125.8","125.9","126.1","126.2","126.3","126.4"
};



void makeSystPlot(RooWorkspace *tempWS, RooWorkspace *WS,  TString systT, int Bin, int toMassNo, int fromMassNo)  
{

  std::cout<< "SystName--> " << systT.Data() << std::endl;
  RooArgList  * hobs = new RooArgList("hobs");

  TString xname;

  if(IFILE.Contains ("Zmm") | IFILE.Contains ("Zee"))
  { 
    xname = "CMS_vhbb_BDT_Zll_8TeV";
  }
  if(IFILE.Contains("We") | IFILE.Contains("Wm")) 
  {
    xname = "CMS_vhbb_BDT_Wln_8TeV";
  }
  if(IFILE.Contains ("Zn") && (Bin == 0))
  {
    xname = "CMS_vhbb_BDT_ZnunuLowPt_8TeV";
  }

  if(IFILE.Contains ("Zn") && (Bin == 1)) 
  {
    xname = "CMS_vhbb_BDT_ZnunuMedPt_8TeV";
  }

  if(IFILE.Contains ("Zn") && (Bin == 2)) 
  {
    xname = "CMS_vhbb_BDT_ZnunuHighPt_8TeV";
  }
  if(IFILE.Contains("Wt")) 
  {
    xname = "BDT";
  }

  std::cout<< "BDT Name--> " << xname.Data() << std::endl;
  hobs->add(*WS->var(xname.Data()));
  RooRealVar BDT(xname.Data(),xname.Data(),-1,1); 

  tempRooDataHistNom = (RooDataHist*)  tempWS->data(systT.Data());
  tempHistNom = tempRooDataHistNom->createHistogram(systT.Data(),BDT,RooFit::Binning(bins));
  
  if(systT.Contains("WH"))
  {
     tempHistNom->Scale(xSec8WH[toMassNo]/xSec8WH[fromMassNo]);
  }  
 
  if(systT.Contains("ZH"))
  {
     tempHistNom->Scale(xSec8ZH[toMassNo]/xSec8ZH[fromMassNo]);
  }  
 
  RooDataHist *DHnom;
  DHnom = new RooDataHist(systT.Data(),"",*hobs,tempHistNom);
  WS->import(*(new RooHistPdf(systT.Data(),"",*hobs,*DHnom))); 
}









void Process(TString fname, TString myRooWS, int toMass, int fromMass, int Bin)
{

  gROOT->SetStyle("Plain");
  TFile * file = new TFile(fname.Data(), "READ");
  std::cout << "reading " << fname.Data() << std::endl;
  TString outname("newcards/");
  outname.Append(massS[toMass]);
  outname.Append("/");
  fname.ReplaceAll("110","");
  fname.ReplaceAll("115","");
  fname.ReplaceAll("120","");
  fname.ReplaceAll("125","");
  fname.ReplaceAll("130","");
  fname.ReplaceAll("135","");
  fname.ReplaceAll("140","");
  fname.ReplaceAll("145","");
  fname.ReplaceAll("150","");
  fname.ReplaceAll("/","");

  if(fname.Contains("Zn") && (Bin == 0)) fname.ReplaceAll(".root","Low.root");
  if(fname.Contains("Zn") && (Bin == 1)) fname.ReplaceAll(".root","Med.root");
  if(fname.Contains("Zn") && (Bin == 2)) fname.ReplaceAll(".root","High.root");
  
  outname.Append(fname.Data());
  outname.ReplaceAll("Wtn_BDT_newBinning_","vhbb_Wtn_8TeV");
  outname.ReplaceAll("WS_BDT_H_","");
  outname.ReplaceAll("WS_BDT_M_","");
  std::cout << "FILENAME: " << outname.Data() <<std::endl;

  TFile * outfile = new TFile(outname.Data(), "RECREATE");

  using namespace RooFit;
  RooWorkspace *myWS = new RooWorkspace(myRooWS.Data(),myRooWS.Data());

  if(fname.Contains ("Zmm") | fname.Contains ("Zee") )  {myWS->factory("CMS_vhbb_BDT_Zll_8TeV[-1.,1.]"); bins = 15;}
  else if(fname.Contains ("Zn") && (Bin == 0))  {myWS->factory("CMS_vhbb_BDT_ZnunuLowPt_8TeV[-1.,1.]"); bins = 24;}
  else if(fname.Contains ("Zn") && (Bin == 1))  {myWS->factory("CMS_vhbb_BDT_ZnunuMedPt_8TeV[-1.,1.]"); bins = 32;}
  else if(fname.Contains ("Zn") && (Bin == 2))  {myWS->factory("CMS_vhbb_BDT_ZnunuHighPt_8TeV[-1.,1.]"); bins = 40;}
  else if(fname.Contains("We") | fname.Contains("Wm")) {myWS->factory("CMS_vhbb_BDT_Wln_8TeV[-1.,1.]"); bins = 48;}
  else if(fname.Contains("Wt")) {myWS->factory("BDT[-1.,1.]"); bins = 18;}


  RooWorkspace *tempWS =  (RooWorkspace*) file->Get(myRooWS.Data());
  
  for(int s=0; s<NS; s++ ){
      makeSystPlot(tempWS, myWS,  systs[s], Bin, toMass, fromMass );
    }
  


  myWS->writeToFile(outname.Data());  
  std::cout << std::endl << std::endl << std::endl << std::endl << "///////////////////////////" << std::endl;
  std::cout << outname.Data() << " written" << std::endl;
  std::cout << "///////////////////////////" << std::endl << std::endl << std::endl;
  outfile->Write();
  outfile->Close();


}

  



void IntermediateMassMaker(int Bin = 0)
{





  TString myRooWS;
  IFILE = "XXX";


  if(IFILE.Contains ("Wmn") && IFILE.Contains("Low") && IFILE.Contains("8TeV"))  myRooWS = "WmnLowPt_8TeV";
  if(IFILE.Contains ("Wen") && IFILE.Contains("Low") && IFILE.Contains("8TeV"))  myRooWS = "WenLowPt_8TeV";
  if(IFILE.Contains ("Wmn") && IFILE.Contains("Mid") && IFILE.Contains("8TeV"))  myRooWS = "WmnMidPt_8TeV";
  if(IFILE.Contains ("Wen") && IFILE.Contains("Mid") && IFILE.Contains("8TeV"))  myRooWS = "WenMidPt_8TeV";
  if(IFILE.Contains ("Wmn") && IFILE.Contains("High") && IFILE.Contains("8TeV"))  myRooWS = "WmnHighPt_8TeV";
  if(IFILE.Contains ("Wen") && IFILE.Contains("High") && IFILE.Contains("8TeV"))  myRooWS = "WenHighPt_8TeV";
  if(IFILE.Contains ("Zn") && (Bin == 0))  myRooWS = "ZnunuLowPt_8TeV";
  if(IFILE.Contains ("Zn") && (Bin == 1))  myRooWS = "ZnunuMedPt_8TeV";
  if(IFILE.Contains ("Zn") && (Bin == 2))  myRooWS = "ZnunuHighPt_8TeV";
  if(IFILE.Contains ("ZmmL")) myRooWS = "ZmmLowPt_8TeV";
  if(IFILE.Contains ("ZmmH")) myRooWS = "ZmmHighPt_8TeV";
  if(IFILE.Contains ("ZeeL")) myRooWS = "ZeeLowPt_8TeV";
  if(IFILE.Contains ("ZeeH")) myRooWS = "ZeeHighPt_8TeV";
  if(IFILE.Contains ("Wt")) myRooWS = "Wtn";


  if((IFILE.Contains("110")))
  {
       Process(IFILE, myRooWS, 0,0, Bin);
       Process(IFILE, myRooWS, 1,0, Bin);
       Process(IFILE, myRooWS, 2,0, Bin);
       Process(IFILE, myRooWS, 3,0, Bin);
     Process(IFILE, myRooWS, 4,0, Bin);
  }
  if((IFILE.Contains("115")))
  {
     for(int to = 5; to < 15; to++)   Process(IFILE, myRooWS, to , 10, Bin);
  }

  if((IFILE.Contains("120")))
  {
     for(int to = 15; to < 25; to++)   Process(IFILE, myRooWS,to , 20, Bin);
  }

  if((IFILE.Contains("125")))
  {
     for(int to = 25; to < 35; to++)   Process(IFILE, myRooWS,to , 30, Bin);
     for(int to = 81; to < 97; to++)   Process(IFILE, myRooWS,to , 30, Bin);
  }

  if((IFILE.Contains("130")))
  {
     for(int to = 35; to < 45; to++)   Process(IFILE, myRooWS,to , 40, Bin);
  }


  if((IFILE.Contains("135")))
  {
     for(int to = 45; to < 55; to++)   Process(IFILE, myRooWS,to , 50, Bin);
  }

   if((IFILE.Contains("140")))
  {
     for(int to = 55; to < 65; to++)   Process(IFILE, myRooWS,to , 60, Bin);
  }

    if((IFILE.Contains("145")))
  {
     for(int to = 65; to < 75; to++)   Process(IFILE, myRooWS,to , 70, Bin);
  }

    if((IFILE.Contains("150")))
  {
     for(int to = 75; to < 81; to++)   Process(IFILE, myRooWS,to , 80, Bin);
  }


  

}




