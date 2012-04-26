#include <iostream>
#include <boost/shared_ptr.hpp>
#include "Math/GenVector/Boost.h"

#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"

#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuSummaryHandler.h"
#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuPhysicsEvent.h"
#include "CMGTools/HtoZZ2l2nu/interface/METUtils.h"
#include "CMGTools/HtoZZ2l2nu/interface/GammaEventHandler.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/plotter.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"
#include "CMGTools/HtoZZ2l2nu/interface/SmartSelectionMonitor.h"
#include "CMGTools/HtoZZ2l2nu/interface/TMVAUtils.h"
#include "CMGTools/HtoZZ2l2nu/interface/MacroUtils.h"
#include "CMGTools/HtoZZ2l2nu/interface/EventCategory.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TEventList.h"
#include "TROOT.h"
 
using namespace std;


int main(int argc, char* argv[])
{
  //##############################################
  //########    GLOBAL INITIALIZATION     ########
  //##############################################

  // check arguments
  if(argc<2){ std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl; exit(0); }
  
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  // configure the process
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  TString url=runProcess.getParameter<std::string>("input");
  TString outdir=runProcess.getParameter<std::string>("outdir");
  TString outUrl( outdir );
  gSystem->Exec("mkdir -p " + outUrl);

  bool isMC = runProcess.getParameter<bool>("isMC");
  int mctruthmode=runProcess.getParameter<int>("mctruthmode");
  
  TString outTxtUrl= outUrl + "/" + gSystem->BaseName(url) + ".txt";
  FILE* outTxtFile = NULL;
  if(!isMC)outTxtFile = fopen(outTxtUrl.Data(), "w");
  printf("TextFile URL = %s\n",outTxtUrl.Data());

  //handler for gamma processes
  GammaEventHandler *gammaEvHandler=0;
  if(mctruthmode==22){
     isMC=false;
     gammaEvHandler = new GammaEventHandler(runProcess);
  }

  //tree info
  int evStart=runProcess.getParameter<int>("evStart");
  int evEnd=runProcess.getParameter<int>("evEnd");
  TString dirname = runProcess.getParameter<std::string>("dirName");

  //systematic uncertainties
  TString uncFile =  runProcess.getParameter<std::string>("jesUncFileName"); gSystem->ExpandPathName(uncFile);
  JetCorrectionUncertainty jecUnc(uncFile.Data());
  bool runSystematics                        = runProcess.getParameter<bool>("runSystematics");
  TString varNames[]={"","_jesup","_jesdown","_jerup","_jerdown","_puup","_pudown","_renup","_rendown","_factup","_factdown","_btagup","_btagdown"};
  size_t nvarsToInclude(1);
  if(runSystematics) 
    { 
      cout << "Systematics will be computed for this analysis" << endl;
      nvarsToInclude=sizeof(varNames)/sizeof(TString);
    }
  
  //isMC_VBF sample ?
  double HiggsMass=0; string VBFString = ""; string GGString("");
  bool isMC_GG  = isMC && ( string(url.Data()).find("GG" )  != string::npos);
  bool isMC_VBF = isMC && ( string(url.Data()).find("VBF")  != string::npos);
  std::vector<TGraph *> hWeightsGrVec;
  if(isMC_GG){
    size_t GGStringpos =  string(url.Data()).find("GG");
    string StringMass = string(url.Data()).substr(GGStringpos+5,3);  sscanf(StringMass.c_str(),"%lf",&HiggsMass);
    GGString = string(url.Data()).substr(GGStringpos);
    TString hqtWeightsFileURL = runProcess.getParameter<std::string>("hqtWeightsFile"); gSystem->ExpandPathName(hqtWeightsFileURL);
    TFile *fin=TFile::Open(hqtWeightsFileURL);
    if(fin){
      cout << "HpT weights (and uncertainties) will be applied from: " << hqtWeightsFileURL << endl;
        for(int ivar=0; ivar<5; ivar++){
          double ren=HiggsMass; if(ivar==ZZ2l2nuSummary_t::hKfactor_renDown)  ren/=2; if(ivar==ZZ2l2nuSummary_t::hKfactor_renUp)  ren *= 2;
          double fac=HiggsMass; if(ivar==ZZ2l2nuSummary_t::hKfactor_factDown) fac/=2; if(ivar==ZZ2l2nuSummary_t::hKfactor_factUp) fac *= 2;
          char buf[100]; sprintf(buf,"kfactors_mh%3.0f_ren%3.0f_fac%3.0f",HiggsMass,ren,fac);
          TGraph *gr= (TGraph *) fin->Get(buf);
          if(gr) hWeightsGrVec.push_back((TGraph *)gr->Clone());
        }
        fin->Close();
        delete fin;
    }
  }else if(isMC_VBF){
    size_t VBFStringpos =  string(url.Data()).find("VBF");
    string StringMass = string(url.Data()).substr(VBFStringpos+6,3);  sscanf(StringMass.c_str(),"%lf",&HiggsMass);
    VBFString = string(url.Data()).substr(VBFStringpos);
  }

  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;
  TH1F* Hcutflow  = (TH1F*) mon.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,5,0,5) ) ;
  TH1F *h=(TH1F*) mon.addHistogram( new TH1F ("eventflow", ";Step;Events", 9,0,9) );
  h->GetXaxis()->SetBinLabel(1,"Preselected");
  h->GetXaxis()->SetBinLabel(2,"|M-M_{Z}|<15");
  h->GetXaxis()->SetBinLabel(3,"Lepton id");
  h->GetXaxis()->SetBinLabel(4,"Lepton iso");
  h->GetXaxis()->SetBinLabel(5,"Z_{pT}>25");
  h->GetXaxis()->SetBinLabel(6,"3^{rd}-lepton veto");
  h->GetXaxis()->SetBinLabel(7,"b-veto");
  h->GetXaxis()->SetBinLabel(8,"#delta #phi(jet,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(9,"E_{T}^{miss}>70");

  mon.addHistogram( new TH1F( "zeta", ";#eta^{ll};Events", 50,-10,10) );
  mon.addHistogram( new TH1F( "zpt", ";p_{T}^{ll};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "zmass", ";M^{ll};Events", 100,40,250) );

  //lepton control
  h=(TH1F *)mon.addHistogram( new TH1F( "zeeselafterid", ";Id;Events", 7,0,7) );
  h->GetXaxis()->SetBinLabel(1,"Veto");
  h->GetXaxis()->SetBinLabel(2,"Loose");
  h->GetXaxis()->SetBinLabel(3,"Medium");
  h->GetXaxis()->SetBinLabel(4,"Tight");
  h->GetXaxis()->SetBinLabel(5,"2011");
  h->GetXaxis()->SetBinLabel(6,"Conversion");
  h->GetXaxis()->SetBinLabel(7,"Trigger+Veto");
  TString lepTypes[]={"e","mu"};
  for(size_t ilep=0; ilep<2; ilep++)
    {
      TString pre=lepTypes[ilep];

      if(pre=="e")
	{
	  mon.addHistogram(new TH1F(pre+"detain",   ";#Delta#eta_{in};Leptons",20,0,0.01) );
	  mon.addHistogram(new TH1F(pre+"dphiin",   ";#Delta#phi_{in};Leptons",20,0,0.1) );
	  mon.addHistogram(new TH1F(pre+"sihih",    ";#sigma_{i#eta i#eta};Leptons",20,0,0.05) );
	  mon.addHistogram(new TH1F(pre+"sipip",    ";#sigma_{i#phi i#phi};Leptons",20,0,0.05) );
	  mon.addHistogram(new TH1F(pre+"r9",       ";R_{9};Leptons",20,0,1.) );
	  mon.addHistogram(new TH1F(pre+"hoe",      ";h/e;Leptons",20,0,0.2) );
	  mon.addHistogram(new TH1F(pre+"ooemoop",  ";1/E-1/p;Leptons",50,0,0.05) );
	  mon.addHistogram(new TH1F(pre+"eopin",    ";E/p;Leptons",50,0,2) );
	  mon.addHistogram(new TH1F(pre+"fbrem",    ";f_{brem};Leptons",50,0,2) );
	}
      if(pre=="mu")
	{
	  mon.addHistogram(new TH1F(pre+"nmatches", ";Muon matches;Leptons",15,0,15) );
	  mon.addHistogram(new TH1F(pre+"nmuonhits", ";Muon hits;Leptons",30,0,30) );
	}
      mon.addHistogram(new TH1F(pre+"d0",            ";d_{0};Leptons",50,0,0.1) );
      mon.addHistogram(new TH1F(pre+"dZ",            ";d_{Z};Leptons",50,0,0.1) );
      mon.addHistogram(new TH1F(pre+"trkchi2",       ";#chi^{2};Leptons",50,0,10) );
      mon.addHistogram(new TH1F(pre+"trkvalidpixel",  ";Valid pixel hits;Leptons",20,0,20) );
      mon.addHistogram(new TH1F(pre+"trkvalidtracker",  ";Valid tracker hits;Leptons",50,0,50) );
      mon.addHistogram(new TH1F(pre+"losthits", ";Lost hits;Leptons",4,0,4) );
      mon.addHistogram(new TH1F(pre+"reliso",        ";RelIso;Leptons",50,0,2) );
      mon.addHistogram(new TH1F(pre+"relisorho",     ";RelIso(#rho);Leptons",50,0,2) );
      mon.addHistogram(new TH1F(pre+"pfreliso",      ";PFRelIso;Leptons",50,0,2) );
      mon.addHistogram(new TH1F(pre+"pfrelisodbeta", ";PFRelIso(#Delta #beta);Leptons",50,0,2) );
      mon.addHistogram(new TH1F(pre+"ecaliso",       ";ECALIso;Leptons",50,0,5) );
      mon.addHistogram(new TH1F(pre+"hcaliso",       ";HCALIso;Leptons",50,0,5) );
      mon.addHistogram(new TH1F(pre+"trkiso",        ";TrkIso;Leptons",50,0,5) );
      mon.addHistogram(new TH1F(pre+"giso",          ";GammaIso;Leptons",50,0,5) );
      mon.addHistogram(new TH1F(pre+"chiso",         ";ChHadronIso;Leptons",50,0,5) );
      mon.addHistogram(new TH1F(pre+"nhiso",         ";NeutralHadron;Leptons",50,0,5) );
      mon.addHistogram(new TH2F(pre+"relisovspu",         ";RelIso;Pileup;Leptons",50,0,50,50,0,2) );
      mon.addHistogram(new TH2F(pre+"relisorhovspu",      ";RelIso(#rho);Pileup;Leptons",50,0,50,50,0,2) ); 
      mon.addHistogram(new TH2F(pre+"pfrelisovspu",       ";PFRelIso;Pileup;Leptons",50,0,50,50,0,2) ); 
      mon.addHistogram(new TH2F(pre+"pfrelisodbetavspu",  ";PFRelIso(#Delta #beta);Pileup;Leptons",50,0,50,50,0,2) );  
      mon.addHistogram(new TH2F(pre+"ecalisovspu",        ";ECALIso;Pileup;Leptons",50,0,50,50,0,5) ); 
      mon.addHistogram(new TH2F(pre+"hcalisovspu",        ";HCALIso;Pileup;Leptons",50,0,50,50,0,5) ); 
      mon.addHistogram(new TH2F(pre+"trkisovspu",         ";TrkIso;Pileup;Leptons",50,0,50,50,0,5) ); 
      mon.addHistogram(new TH2F(pre+"gisovspu",           ";GammaIso;Pileup;Leptons",50,0,50,50,0,5) ); 
      mon.addHistogram(new TH2F(pre+"chisovspu",          ";ChHadronIso;Pileup;Leptons",50,0,50,50,0,5) ); 
      mon.addHistogram(new TH2F(pre+"nhisovspu",          ";NeutralHadronIso;Pileup;Leptons",50,0,50,50,0,5) ); 
    }


  mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "rho",";#rho;Events",50,0,25) ); 
  mon.addHistogram( new TH1F( "rho25",";#rho(#eta<2.5);Events",50,0,25) ); 

  TString jetRegs[]={"TK","HEin","HEout","HF"};
  for(size_t ireg=0; ireg<4; ireg++)
    {
      mon.addHistogram( new TH1F(jetRegs[ireg]+"jetbeta"       , ";#beta;Events",50,0,1) );
      mon.addHistogram( new TH1F(jetRegs[ireg]+"jetbetastar"       , ";#beta *;Events",50,0,1) );
      mon.addHistogram( new TH1F(jetRegs[ireg]+"jetdrmean"       , ";<#Delta R>;Events",50,0,0.5) );
      mon.addHistogram( new TH1F(jetRegs[ireg]+"jetptd"       , ";p_{T}D [GeV/c];Events",50,0,1) );
      mon.addHistogram( new TH1F(jetRegs[ireg]+"jetptrms"       , ";RMS p_{T} [GeV/c];Events",50,0,0.5) );
      mon.addHistogram( new TH1F(jetRegs[ireg]+"jetmva"      , ";MVA;Events",50,0,1) );
    }
  mon.addHistogram( new TH1F("jetpt"       , ";p_{T} [GeV/c];Events",50,0,250) );
  mon.addHistogram( new TH1F("jeteta"       , ";|#eta|;Events",25,0,5) );
  mon.addHistogram( new TH1F("njets"       , ";Jet multiplicity (p_{T}>30 GeV/c);Events",5,0,5) );
  int jetids[]  ={JETID_LOOSE, JETID_TIGHT,JETID_PHILV1_LOOSE, JETID_PHILV1_MEDIUM, JETID_MIN_LOOSE, JETID_MIN_MEDIUM, JETID_OPT_LOOSE, JETID_OPT_MEDIUM};
  for(size_t ijetid=0; ijetid<sizeof(jetids)/sizeof(int); ijetid++)
    {
      char buf[15];
      sprintf(buf,"njets%d",jetids[ijetid]);
      mon.addHistogram( new TH1F(TString(buf), ";Jet multiplicity (p_{T}>30 GeV/c);Events",5,0,5) );
      mon.addHistogram( new TH1F(TString(buf)+"lowpu", ";Jet multiplicity (p_{T}>30 GeV/c);Events",5,0,5) );
      mon.addHistogram( new TH1F(TString(buf)+"cenpu", ";Jet multiplicity (p_{T}>30 GeV/c);Events",5,0,5) );
      mon.addHistogram( new TH1F(TString(buf)+"highpu", ";Jet multiplicity (p_{T}>30 GeV/c);Events",5,0,5) );

      char buf2[100];
      sprintf(buf2,"met_clusteredmet%d",jetids[ijetid]);
      mon.addHistogram( new TH1F( buf2  , ";clustered E_{T}^{miss};Events", 50,0,500) );
      mon.addHistogram( new TH2F( TString(buf2)+"_vspu" , ";Pileup events; clustered E_{T}^{miss};Events",50,0,50,50,0,500) );
    }

  mon.addHistogram( new TH1F ("nbtags", ";b-tag multiplicity; Events", 5,0,5) );  
  for(size_t ibin=1; ibin<=5; ibin++){
    TString label("");
    if(ibin==5) label +="#geq";
    else        label +="=";
    label += (ibin-1);
    mon.getHisto("njets")->GetXaxis()->SetBinLabel(ibin,label);
    mon.getHisto("nbtags")->GetXaxis()->SetBinLabel(ibin,label);
  }
  mon.addHistogram( new TH1F( "met_met"  , ";E_{T}^{miss};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "met_min3Met"  , ";min(E_{T}^{miss},assoc-E_{T}^{miss},clustered-E_{T}^{miss});Events", 50,0,500) );
  mon.addHistogram( new TH1F( "met_redMet"  , ";red(E_{T}^{miss},clustered-E_{T}^{miss});Events", 50,0,500) );
  mon.addHistogram( new TH1F( "met_redMetL"  , ";red(E_{T}^{miss},clustered-E_{T}^{miss}) - longi.;Events", 50,-250,250) );
  mon.addHistogram( new TH1F( "met_redMetT"  , ";red(E_{T}^{miss},clustered-E_{T}^{miss}) - perp.;Events", 50,-250,250) );
  mon.addHistogram( new TH1F( "mt"  , ";M_{T};Events", 100,0,1000) );

  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################

  //open the file and get events tree
  ZZ2l2nuSummaryHandler evSummaryHandler;
  TFile *file = TFile::Open(url);
  printf("Looping on %s\n",url.Data());
  if(file==0) return -1;
  if(file->IsZombie()) return -1;
  if( !evSummaryHandler.attachToTree( (TTree *)file->Get(dirname) ) ){
      file->Close();
      return -1;
  }


  //check run range to compute scale factor (if not all entries are used)
  const Int_t totalEntries= evSummaryHandler.getEntries();
  float rescaleFactor( evEnd>0 ?  float(totalEntries)/float(evEnd-evStart) : -1 );
  if(evEnd<0 || evEnd>evSummaryHandler.getEntries() ) evEnd=totalEntries;
  if(evStart > evEnd ){
      file->Close();
      return -1;
  }

  //MC normalization (to 1/pb)
  float cnorm=1.0;
  if(isMC){
      TH1F* cutflowH = (TH1F *) file->Get("evAnalyzer/h2zz/cutflow");
      if(cutflowH) cnorm=cutflowH->GetBinContent(1);
      if(rescaleFactor>0) cnorm /= rescaleFactor;
      printf("cnorm = %f\n",cnorm);
  }
  Hcutflow->SetBinContent(1,cnorm);


  //pileup weighting: based on vtx for now...
  std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
  std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
  std::vector<float> mcPileupDistribution(dataPileupDistribution.size(),0);
  if(isMC)
    {
      for( int iev=evStart; iev<evEnd; iev++)
	{
	  evSummaryHandler.getEntry(iev);
	  int nvtx=evSummaryHandler.evSummary_.nvtx;
	  if(nvtx>int(mcPileupDistribution.size()-1)) continue;
	  mcPileupDistribution[nvtx] += 1;
      }
    }
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
  edm::LumiReWeighting *LumiWeights=0;
  reweight::PoissonMeanShifter *PShiftUp=0, *PShiftDown=0;
  if(isMC)
    {
      LumiWeights= new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
      PShiftUp = new reweight::PoissonMeanShifter(+0.8);
      PShiftDown = new reweight::PoissonMeanShifter(-0.8);
    }

  //event Categorizer
  EventCategory eventCategoryInst(0); //inclusive analysis

  //##############################################
  //########           EVENT LOOP         ########
  //##############################################
  //loop on all the events
  printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Scanning the ntuple :");
  int treeStep = (evEnd-evStart)/50;if(treeStep==0)treeStep=1;
  for( int iev=evStart; iev<evEnd; iev++){
      if((iev-evStart)%treeStep==0){printf(".");fflush(stdout);}

      //##############################################   EVENT LOOP STARTS   ##############################################
   
      //load the event content from tree
      evSummaryHandler.getEntry(iev);
      ZZ2l2nuSummary_t &ev=evSummaryHandler.getEvent();
      PhysicsEvent_t phys=getPhysicsEventFrom(ev);
      
      
      //categorize events
      TString tag_cat;
      switch(ev.cat){
      case EMU : tag_cat = "emu";   break;
      case MUMU: tag_cat = "mumu";  break;
      case EE  : tag_cat = "ee";    break;
      default  : continue;
      }
      //      if(isMC && mctruthmode==1 && !isDYToLL(ev.mccat) && !isZZ2l2nu(ev.mccat) ) continue;
      if(isMC && mctruthmode==1 && !isDYToLL(ev.mccat) ) continue;
      if(isMC && mctruthmode==2 && !isDYToTauTau(ev.mccat) ) continue;
      

      bool isGammaEvent = false;
      if(gammaEvHandler)
	{
          isGammaEvent=gammaEvHandler->isGood(phys);
          if(mctruthmode==22 && !isGammaEvent) continue;
          tag_cat = "gamma";
	}
     
      int eventSubCat  = eventCategoryInst.Get(phys);
      TString tag_subcat = eventCategoryInst.GetLabel(eventSubCat);

      //prepare the tag's vectors for histo filling
      std::vector<TString> tags_full;
      tags_full.push_back("all");
      tags_full.push_back(tag_cat);
//      if(tag_subcat=="vbf") tags_full.push_back(tag_cat+"_"+tag_subcat);
      if(tag_subcat!="")tags_full.push_back(tag_cat + tag_subcat);
      if(tag_subcat=="geq2jets" || tag_subcat=="vbf")tags_full.push_back(tag_cat + "geq2jetsInc");
      if(tag_subcat!="vbf")tags_full.push_back(tag_cat + "novbf");

      //pileup and Higgs pT weight
      //float weight=ev.puWeight;
      float weight = 1.0;
      double TotalWeight_plus = 1.0;
      double TotalWeight_minus = 1.0;
      if(isMC){
        weight = LumiWeights->weight(ev.nvtx);
        TotalWeight_plus = PShiftUp->ShiftWeight(ev.nvtx);
        TotalWeight_minus = PShiftDown->ShiftWeight(ev.nvtx);
        if(isMC_VBF) weight *= weightVBF(VBFString,HiggsMass, phys.genhiggs[0].mass() );         
        if(isMC_GG)  {
          for(size_t iwgt=0; iwgt<hWeightsGrVec.size(); iwgt++) ev.hptWeights[iwgt] = hWeightsGrVec[iwgt]->Eval(phys.genhiggs[0].pt());
          weight *= ev.hptWeights[0];
        }
      }
      Hcutflow->Fill(1,1);
      Hcutflow->Fill(2,weight);
      Hcutflow->Fill(3,weight*TotalWeight_minus);
      Hcutflow->Fill(4,weight*TotalWeight_plus);


      //analyze the leptons
      LorentzVector lep1=phys.leptons[0];
      LorentzVector lep2=phys.leptons[1];
      std::vector<bool> eidPass(7,false);
      if(ev.cat==EE)
	{
	  int lpid[]={phys.leptons[0].pid,phys.leptons[1].pid};
	  int lid[]={ev.en_idbits[lpid[0]],ev.en_idbits[lpid[1]]};
	  bool lpassConversion[]={hasObjectId(lid[0], EID_CONVERSIONVETO),hasObjectId(lid[1],EID_CONVERSIONVETO)};
	  int cutBasedIdsToTest[]={EgammaCutBasedEleId::VETO, EgammaCutBasedEleId::LOOSE, EgammaCutBasedEleId::MEDIUM, EgammaCutBasedEleId::TIGHT };
	  for(int iid=0; iid<4; iid++)
	    {
	      bool result=true;
	      for(size_t ilep=0; ilep<2; ilep++)
		{
		  result &= EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::WorkingPoint(cutBasedIdsToTest[iid]),
							fabs(phys.leptons[0].eta()<1.4442),
							phys.leptons[0].pt(), phys.leptons[0].eta(),
							ev.en_detain[lpid[ilep]],  ev.en_dphiin[lpid[ilep]], ev.en_sihih[lpid[ilep]], ev.en_hoe[lpid[ilep]],
							ev.en_ooemoop[lpid[ilep]], phys.leptons[0].d0, phys.leptons[0].dZ,
							0., 0., 0.,
							!lpassConversion[ilep],0,ev.rho);
		}
	      eidPass[iid]=result;
	    }
	  eidPass[4] = (hasObjectId(lid[0], EID_VBTF2011)       && hasObjectId(lid[1],EID_VBTF2011));
	}

      LorentzVector zll=lep1+lep2;

      //analyze JET/MET
      LorentzVectorCollection jetsP4;
      std::vector<double> genJetsPt;
      for(size_t ijet=0; ijet<phys.ajets.size(); ijet++)
	{
	  jetsP4.push_back( phys.ajets[ijet] );
	  genJetsPt.push_back( phys.ajets[ijet].genPt);
	}
      //base raw METs
      LorentzVector assocMetP4(phys.met[1]);
      LorentzVector zvvRaw(phys.met[0]);
      LorentzVector rawClusteredMet(zll);            rawClusteredMet *= -1;
      for(size_t ijet=0; ijet<jetsP4.size(); ijet++) rawClusteredMet -= jetsP4[ijet];
      LorentzVector redMetRaw(METUtils::redMET(METUtils::INDEPENDENTLYMINIMIZED, lep1, 0, lep2, 0, rawClusteredMet, zvvRaw,false));
      Float_t mtRaw( METUtils::transverseMass(zll,zvvRaw,true) );
      
      //prepare variations (first variation is the baseline, corrected for JER) 
      LorentzVectorCollection zvvs,redMets, min3Mets;
      std::vector<Float_t>  mts,mt3s,redMetLs,redMetTs;
      std::vector<LorentzVectorCollection> jets;
      METUtils::computeVariation(jetsP4, genJetsPt, zvvRaw, jets, zvvs, &jecUnc);
      for(size_t ivars=0; ivars<zvvs.size(); ivars++)
	{
	  LorentzVector clusteredMetP4(zll); clusteredMetP4 *= -1;
	  for(size_t ijet=0; ijet<jets[ivars].size(); ijet++) clusteredMetP4 -= jets[ivars][ijet];
	  METUtils::stRedMET redMetOut; 
	  redMets.push_back(   METUtils::redMET(METUtils::INDEPENDENTLYMINIMIZED, lep1, 0, lep2, 0, clusteredMetP4, zvvs[ivars],false,&redMetOut) );
	  redMetLs.push_back( redMetOut.redMET_l );
	  redMetTs.push_back( redMetOut.redMET_t );
	  mts.push_back(METUtils::transverseMass(zll,zvvs[ivars],true));
	  mt3s.push_back(phys.leptons.size()>2 ? METUtils::transverseMass(phys.leptons[2],zvvs[ivars],false) : 0. );
	  min3Mets.push_back( min(zvvs[ivars], min(assocMetP4,clusteredMetP4)) );
	}
      
      //
      //run the variations
      //
      for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
	float iweight = weight;                                               //nominal
	if(ivar==5)                        iweight *=TotalWeight_plus;        //pu up
	if(ivar==6)                        iweight *=TotalWeight_minus;       //pu down
	if(ivar<=10 && ivar>=7 && isMC_GG) iweight *=ev.hptWeights[ivar-6];   //ren/fact scales
	
	Float_t zmass=zll.mass();
	Float_t zpt=zll.pt();
	Float_t zeta=zll.eta();
	LorentzVectorCollection &origJetsP4=jets[ivar>4?0:ivar];            
	LorentzVector zvv    = zvvs[ivar>4?0:ivar];
	LorentzVector min3Met = min3Mets[ivar>4?0:ivar];
	LorentzVector redMet = redMets[ivar>4?0:ivar];
	Float_t redMetL      = redMetLs[ivar>4?0:ivar];
	Float_t redMetT      = redMetTs[ivar>4?0:ivar];
	Float_t mt           = mts[ivar>4?0:ivar];
	Float_t mt3          = mt3s[ivar>4?0:ivar]; 
	int nBtaggedVsDisc[6]={0,0,0,0,0,0};
	int njets(0),nbtags(0);

	int nidjets[]                    = { 0,           0,          0,                  0,                   0,               0,                0,               0};
	LorentzVector idclusteredmets[]  = { -zll,        -zll,       -zll,               -zll,                -zll,            -zll,             -zll,            -zll};
	Float_t mindphijmet(9999.);
	Float_t btagcut(2.0); if(ivar==10) btagcut=2.03; else if(ivar==11) btagcut=1.97;
	for(size_t ijet=0; ijet<phys.ajets.size(); ijet++)
	  {
	    
	    Float_t idphijmet=fabs(deltaPhi(zvv.phi(),phys.ajets[ijet].phi()));
	    if(idphijmet<mindphijmet) mindphijmet=idphijmet;
	    for(size_t ijetid=0; ijetid<sizeof(jetids)/sizeof(int); ijetid++)
	      {
		if( hasObjectId( phys.ajets[ijet].pid, jetids[ijetid] ) )
		  idclusteredmets[ijetid] -= phys.ajets[ijet];
	      }
	    
	    if(origJetsP4[ijet].pt()>30)
	      {
		njets++;
		for(size_t ijetid=0; ijetid<sizeof(jetids)/sizeof(int); ijetid++)  nidjets[ijetid] += hasObjectId( phys.ajets[ijet].pid, jetids[ijetid] );
		
		nbtags += (phys.ajets[ijet].btag1>btagcut);
		nBtaggedVsDisc[0] +=(phys.ajets[ijet].btag1>1.7);
		nBtaggedVsDisc[1] +=(phys.ajets[ijet].btag1>2.0);
		nBtaggedVsDisc[2] +=(phys.ajets[ijet].btag1>3.3);
		nBtaggedVsDisc[3] +=(phys.ajets[ijet].btag2>0.244);
		nBtaggedVsDisc[4] +=(phys.ajets[ijet].btag2>0.679);
		nBtaggedVsDisc[5] +=(phys.ajets[ijet].btag2>0.898);
	      }
	  }
	
	//##############################################
	//########     PRESELECTION             ########
	//##############################################
	if(zmass<40) continue; // this is required by the HZZ skim anyhow
	bool passZmass(fabs(zmass-91)<15);
	bool passId(true); if(ev.cat==EE && !eidPass[1]) passId=false;
	bool passIso(phys.leptons[0].relIsoRho(ev.rho)<0.15 && phys.leptons[1].relIsoRho(ev.rho)<0.15);
	bool isZsideBand( (zmass>40 && zmass<70) || (zmass>110 && zmass<200));
	bool isZsideBandPlus( (zmass>110 && zmass<200));
	bool passDphijmet(mindphijmet>0.5); 
	bool passZpt(zpt>55);
	bool pass3dLeptonVeto(true); int nExtraLep(0); for(unsigned int i=2;i<phys.leptons.size();i++) { if(phys.leptons[i].pt()>10){ nExtraLep++; pass3dLeptonVeto=false;} }
	bool passBveto(nbtags==0);
	bool passBaseMet(zvv.pt()>70);
     

// 	LorentzVector genRes(0,0,0,0);
// 	for(size_t igl=0;igl<phys.genleptons.size(); igl++) genRes+= phys.genleptons[igl];
// 	genRes += phys.genmet[0];
// 	cout << genRes.mass() << " " << genRes.px() << " " << genRes.py() << " " << genRes.pz() << endl;
    
	//##############################################  
	//########         GENERAL PLOTS        ########                                                                                                                  
	//##############################################  

	//N-1 PLOTS FOR PRESELECTION
	if(ivar==0){
	  mon.fillHisto  ("eventflow",tags_full,0,iweight);
	  if(passId && passIso) mon.fillHisto  ("zmass"    ,tags_full,zmass,iweight);
	  
	  if(passZmass){

	    if(ev.cat==EE)
	      {
		for(size_t eidPassCtr=0; eidPassCtr<7; eidPassCtr++)
		  if(eidPass[eidPassCtr]) 
		    mon.fillHisto("zeeselafterid",tags_full,eidPassCtr,weight);
	      }
	    for(size_t ilep=0; ilep<2; ilep++)
	      {
		TString pre("mu");
		if(fabs(phys.leptons[ilep].id)==ELECTRON) pre ="e";
		if(pre=="e")
		  {
		    mon.fillHisto(pre+"detain",  tags_full,fabs(ev.en_detain[phys.leptons[ilep].pid]),weight);
		    mon.fillHisto(pre+"dphiin",  tags_full,fabs(ev.en_dphiin[phys.leptons[ilep].pid]),weight);
		    mon.fillHisto(pre+"sihih",   tags_full,fabs(ev.en_sihih[phys.leptons[ilep].pid]),weight);
		    mon.fillHisto(pre+"sipip",   tags_full,fabs(ev.en_sipip[phys.leptons[ilep].pid]),weight);
		    mon.fillHisto(pre+"r9",      tags_full,fabs(ev.en_r9[phys.leptons[ilep].pid]),weight);
		    mon.fillHisto(pre+"hoe",     tags_full,fabs(ev.en_hoe[phys.leptons[ilep].pid]),weight);
		    mon.fillHisto(pre+"eopin",   tags_full,fabs(ev.en_eopin[phys.leptons[ilep].pid]),weight);
		    mon.fillHisto(pre+"fbrem",   tags_full,fabs(ev.en_fbrem[phys.leptons[ilep].pid]),weight);
		    mon.fillHisto(pre+"ooemoop", tags_full,fabs(ev.en_ooemoop[phys.leptons[ilep].pid]),weight);
		  }
		else
		  {
		    mon.fillHisto(pre+"nmuonhits",  tags_full,fabs(ev.mn_validMuonHits[phys.leptons[ilep].pid]),weight);
		    mon.fillHisto(pre+"nmatches",   tags_full,fabs(ev.mn_nMatches[phys.leptons[ilep].pid]),weight);
		  }
		mon.fillHisto(pre+"d0",              tags_full,fabs(phys.leptons[ilep].d0),weight);
		mon.fillHisto(pre+"dZ",              tags_full,fabs(phys.leptons[ilep].dZ),weight);
		mon.fillHisto(pre+"trkchi2",         tags_full,fabs(phys.leptons[ilep].trkchi2),weight);
		mon.fillHisto(pre+"trkvalidpixel",   tags_full,fabs(phys.leptons[ilep].trkValidPixelHits),weight);
		mon.fillHisto(pre+"trkvalidtracker", tags_full,fabs(phys.leptons[ilep].trkLostInnerHits),weight);
		mon.fillHisto(pre+"losthits",        tags_full,fabs(phys.leptons[ilep].trkLostInnerHits),weight);
	      }		
	    
	    mon.fillHisto("eventflow",tags_full,1,iweight);
	    mon.fillHisto("zeta"     ,tags_full,zeta   ,iweight);
	    mon.fillHisto("zpt"      ,tags_full,zpt     ,iweight);
	    mon.fillHisto("nvtx"     ,tags_full,ev.nvtx,iweight);
	    mon.fillHisto("rho"      ,tags_full,ev.rho,iweight);
	    mon.fillHisto("rho25"    ,tags_full,ev.rho25,iweight);
	    
	    if(passId)
	      {
		mon.fillHisto("eventflow",tags_full,2,iweight);
		for(size_t ilep=0; ilep<2; ilep++)
		  {
		    TString pre("mu");
		    if(fabs(phys.leptons[ilep].id)==ELECTRON) pre ="e";
		    mon.fillHisto(pre+"reliso",        tags_full,phys.leptons[ilep].relIso(),          iweight);
		    mon.fillHisto(pre+"relisorho",     tags_full,phys.leptons[ilep].relIsoRho(ev.rho), iweight);
		    mon.fillHisto(pre+"pfreliso",      tags_full,phys.leptons[ilep].pfRelIso(),        iweight);
		    mon.fillHisto(pre+"pfrelisodbeta", tags_full,phys.leptons[ilep].pfRelIsoDbeta(),   iweight);
		    mon.fillHisto(pre+"ecaliso",       tags_full,phys.leptons[ilep].ecalIso,           iweight);
		    mon.fillHisto(pre+"hcaliso",       tags_full,phys.leptons[ilep].hcalIso,           iweight);
		    mon.fillHisto(pre+"trkiso",        tags_full,phys.leptons[ilep].trkIso,            iweight);
		    mon.fillHisto(pre+"giso",          tags_full,phys.leptons[ilep].gIso,              iweight);
		    mon.fillHisto(pre+"chiso",         tags_full,phys.leptons[ilep].chIso,             iweight);
		    mon.fillHisto(pre+"nhiso",         tags_full,phys.leptons[ilep].nhIso,             iweight);
		    if(isMC)
		      {
			mon.fill2DHisto(pre+"relisovspu",        tags_full, ev.ngenITpu, phys.leptons[ilep].relIso(),          iweight);
			mon.fill2DHisto(pre+"relisorhovspu",     tags_full, ev.ngenITpu, phys.leptons[ilep].relIsoRho(ev.rho), iweight);
			mon.fill2DHisto(pre+"pfrelisovspu",      tags_full, ev.ngenITpu, phys.leptons[ilep].pfRelIso(),        iweight);
			mon.fill2DHisto(pre+"pfrelisodbetavspu", tags_full, ev.ngenITpu, phys.leptons[ilep].pfRelIsoDbeta(),   iweight);
			mon.fill2DHisto(pre+"ecalisovspu",       tags_full, ev.ngenITpu, phys.leptons[ilep].ecalIso,           iweight);
			mon.fill2DHisto(pre+"hcalisovspu",       tags_full, ev.ngenITpu, phys.leptons[ilep].hcalIso,           iweight);
			mon.fill2DHisto(pre+"trkisovspu",        tags_full, ev.ngenITpu, phys.leptons[ilep].trkIso,            iweight);
			mon.fill2DHisto(pre+"gisovspu",          tags_full, ev.ngenITpu, phys.leptons[ilep].gIso,              iweight);
			mon.fill2DHisto(pre+"chisovspu",         tags_full, ev.ngenITpu, phys.leptons[ilep].chIso,             iweight);
			mon.fill2DHisto(pre+"nhisovspu",         tags_full, ev.ngenITpu, phys.leptons[ilep].nhIso,             iweight);
		      }
		  }

		if(passIso)
		  {
		    mon.fillHisto("eventflow",tags_full,3,iweight);

		    if(passZpt){
		      mon.fillHisto("eventflow",tags_full,4,iweight);
                      
		      if(pass3dLeptonVeto){
			mon.fillHisto("eventflow",tags_full,6,iweight);
			mon.fillHisto("nbtags",tags_full, nbtags,iweight);
			for(int ibtag=0; ibtag<6; ibtag++)
			  {
			    if(nBtaggedVsDisc[ibtag]==0)
			      mon.fillHisto("btagvetosel",tags_full,ibtag,iweight);
			  }
			
			if(passBveto){
			  mon.fillHisto("eventflow",tags_full,7,iweight);
			  mon.fillHisto("njets",tags_full,njets,iweight);
			  mon.fillHisto("mindphijmet",tags_full,mindphijmet,iweight);
			  for(size_t ijetid=0; ijetid<sizeof(jetids)/sizeof(int); ijetid++)  
			    {
			      char buf[15];
			      sprintf(buf,"njets%d",jetids[ijetid]);
			      mon.fillHisto(buf,tags_full,nidjets[ijetid],iweight);

			      TString puReg("lowpu");
			      if(ev.ngenITpu>15) puReg="cenpu";
			      if(ev.ngenITpu>25) puReg="highpu";
			      mon.fillHisto(buf+puReg,tags_full,nidjets[ijetid],iweight);

			      char buf2[100];
			      sprintf(buf2,"met_clusteredmet%d",jetids[ijetid]);
			      mon.fillHisto(buf2,tags_full  , idclusteredmets[ijetid].pt(),iweight);
			      mon.fill2DHisto(TString(buf2)+"_vspu" ,tags_full, ev.ngenITpu,idclusteredmets[ijetid].pt(),iweight);
			    }

			  for(size_t ijet=0; ijet<phys.ajets.size(); ijet++)
			    {
			      TString reg("TK");
			      if(fabs(phys.ajets[ijet].eta())>2.5)  reg="HEin";
			      if(fabs(phys.ajets[ijet].eta())>2.75) reg="HEout";
			      if(fabs(phys.ajets[ijet].eta())>3)    reg="HF";

			      mon.fillHisto(reg+"jetbeta",tags_full,phys.ajets[ijet].beta,iweight);
			      mon.fillHisto(reg+"jetbetastar",tags_full,phys.ajets[ijet].betaStar,iweight);
			      mon.fillHisto(reg+"jetdrmean",tags_full,phys.ajets[ijet].dRMean,iweight);
			      mon.fillHisto(reg+"jetptd",tags_full,phys.ajets[ijet].ptD,iweight);
			      mon.fillHisto(reg+"jetptrms",tags_full,phys.ajets[ijet].ptRMS,iweight);
			      
			      mon.fillHisto(reg+"jetmva",tags_full,phys.ajets[ijet].pumva,iweight);
			      
			      mon.fillHisto("jetpt",tags_full,phys.ajets[ijet].pt(),iweight);
			      mon.fillHisto("jeteta",tags_full,fabs(phys.ajets[ijet].eta()),iweight);
			    }
			  
			  
			  if(passDphijmet){
			    mon.fillHisto("eventflow",tags_full,8,iweight);
			    mon.fillHisto("met_met",tags_full,zvv.pt(),iweight);
			    mon.fillHisto("met_min3Met",tags_full,min3Met.pt(),iweight);
			    mon.fillHisto("met_met_vspu",tags_full,ev.ngenITpu,zvv.pt(),iweight);
			    mon.fillHisto("met_metRaw",tags_full,zvvRaw.pt(),iweight);
			    mon.fillHisto("met_redMet",tags_full,redMet.pt(),iweight);
			    mon.fillHisto("met_redMet_vspu",tags_full,ev.ngenITpu,redMet.pt(),iweight);
			    mon.fillHisto("met_min3Met_vspu",tags_full,ev.ngenITpu,min3Met.pt(),iweight);
			    mon.fillHisto("met_redMetL",tags_full,redMetT,iweight);
			    mon.fillHisto("met_redMetT",tags_full,redMetL,iweight);
			    mon.fillHisto("met_redMetRaw",tags_full,redMetRaw.pt(),iweight);
			    mon.fillHisto("mt",tags_full,mt,iweight);
			    mon.fillHisto("mtRaw",tags_full,mtRaw,iweight);
			    
			    if(passBaseMet){
			      mon.fillHisto  ("eventflow",tags_full,9,iweight);
			    }
			  }
			}
		      }
		    }
		  }
	      }
	  }
	}
      }
      
      //##############################################   EVENT LOOP ENDS   ##############################################
  }
  
  printf("\n"); 
  file->Close();

  //##############################################
  //########     SAVING HISTO TO FILE     ########
  //##############################################
  //save control plots to file
  outUrl += "/";
  outUrl += gSystem->BaseName(url);
  printf("Results save in %s\n", outUrl.Data());

  //save all to the file
  TFile *ofile=TFile::Open(outUrl, "recreate");
  mon.Write();
  ofile->Close();

  if(outTxtFile)fclose(outTxtFile);
}  





