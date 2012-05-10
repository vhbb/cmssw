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

int muonIdLevel(int id, float chi2,float muonHits, float muonMatches, float d0, float dZ, float pixelHits, float trkLayers)
{
  if(!hasObjectId(id,MID_PF) || !hasObjectId(id,MID_GLOBAL))    return -1;
  if(chi2>10)                     return 1;
  if(muonHits==0)                 return 2;
  if(muonMatches<=1)              return 3;
  if(fabs(d0)>0.2)                return 4;
  if(fabs(dZ)>0.5)                return 5;
  if(pixelHits==0)                return 6;
  if(trkLayers<=5)                return 7;
  return 8;
}

TString getJetRegion(float eta)
{
  TString reg("TK");
  if(fabs(eta)>2.5)  reg="HEin";
  if(fabs(eta)>2.75) reg="HEout";
  if(fabs(eta)>3)    reg="HF";
  return reg;
}


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
  bool runBlinded = runProcess.getParameter<bool>("runBlinded"); 
  int mctruthmode=runProcess.getParameter<int>("mctruthmode");
  
  TString outTxtUrl= outUrl + "/" + gSystem->BaseName(url) + ".txt";
  FILE* outTxtFile = NULL;
  if(!isMC)outTxtFile = fopen(outTxtUrl.Data(), "w");
  printf("TextFile URL = %s\n",outTxtUrl.Data());

  //tree info
  int evStart=runProcess.getParameter<int>("evStart");
  int evEnd=runProcess.getParameter<int>("evEnd");
  TString dirname = runProcess.getParameter<std::string>("dirName");

  //jet energy scale uncertainties
  TString uncFile =  runProcess.getParameter<std::string>("jesUncFileName"); gSystem->ExpandPathName(uncFile);
  JetCorrectionUncertainty jecUnc(uncFile.Data());


  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;
  TH1F* Hcutflow  = (TH1F*) mon.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,5,0,5) ) ;
  TH1F *h=(TH1F*) mon.addHistogram( new TH1F ("eventflow", ";Step;Events", 9,0,9) );
  h->GetXaxis()->SetBinLabel(1,"Preselected");
  h->GetXaxis()->SetBinLabel(2,"|M-M_{Z}|<15");
  h->GetXaxis()->SetBinLabel(3,"Lepton id+iso");
  h->GetXaxis()->SetBinLabel(4,"p_{T}^{ll}>55");
  h->GetXaxis()->SetBinLabel(5,"b-veto"); 
  h->GetXaxis()->SetBinLabel(6,"3^{rd}-lepton veto");
  h->GetXaxis()->SetBinLabel(7,"#delta #phi(jet,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(8,"E_{T}^{miss}>70");

  mon.addHistogram( new TH1F( "leadpt", ";p_{T}^{l};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "leadeta", ";#eta^{l};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "trailerpt", ";p_{T}^{l};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "trailereta", ";#eta^{l};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "zpt", ";p_{T}^{ll};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "zeta", ";#eta^{ll};Events", 50,-10,10) );
  mon.addHistogram( new TH1F( "zmass", ";M^{ll};Events", 100,40,250) );

  mon.addHistogram( new TH1F( "thirdleptonpt", ";p_{T}^{l};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "thirdleptoneta", ";#eta^{l};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "thirdleptonmt", ";M_{T} [GeV/c^{2}];Events", 50,0,200) );

  //lepton control
  for(size_t ilep=0; ilep<2; ilep++)
    {
      TString lepStr(ilep==0? "mu" :"e");
      mon.addHistogram(new TH1F(lepStr+"genpt",   ";p_{T} [GeV/c];Leptons",50,0,250) );
      mon.addHistogram(new TH1F(lepStr+"geneta",   ";#eta;Leptons",50,-5,5) );
      mon.addHistogram(new TH1F(lepStr+"genpu",   ";Pileup;Leptons",50,0,50) );
      for(int iid=0; iid<4; iid++)
	{
	  TString idctr(""); idctr+=iid;
	  mon.addHistogram(new TH1F(lepStr+idctr+"pt",   ";p_{T} [GeV/c];Leptons",50,0,250) );
	  mon.addHistogram(new TH1F(lepStr+idctr+"eta",   ";#eta;Leptons",50,-5,5) );
	  mon.addHistogram(new TH1F(lepStr+idctr+"pu",   ";Pileup;Leptons",50,0,50) );
	  mon.addHistogram(new TH1F(lepStr+idctr+"isopt",   ";p_{T} [GeV/c];Leptons",50,0,250) );
	  mon.addHistogram(new TH1F(lepStr+idctr+"isoeta",   ";#eta;Leptons",50,-5,5) );
	  mon.addHistogram(new TH1F(lepStr+idctr+"isopu",   ";Pileup;Leptons",50,0,50) );
	}
      if(ilep==1)
	{
	  for(size_t ireg=0; ireg<2; ireg++)
	    { 
	      TString reg(ireg==0?"eb":"ee");
	      mon.addHistogram(new TH1F(lepStr+reg+"detain",   ";#Delta#eta_{in};Leptons",50,0,0.01) );
	      mon.addHistogram(new TH1F(lepStr+reg+"dphiin",   ";#Delta#phi_{in};Leptons",50,0,0.1) );
	      mon.addHistogram(new TH1F(lepStr+reg+"sihih",    ";#sigma_{i#eta i#eta};Leptons",50,0,0.05) );
	      mon.addHistogram(new TH1F(lepStr+reg+"sipip",    ";#sigma_{i#phi i#phi};Leptons",50,0,0.05) );
	      mon.addHistogram(new TH1F(lepStr+reg+"r9",       ";R_{9};Leptons",50,0,1.) );
	      mon.addHistogram(new TH1F(lepStr+reg+"hoe",      ";h/e;Leptons",50,0,0.2) );
	      mon.addHistogram(new TH1F(lepStr+reg+"ooemoop",  ";1/E-1/p;Leptons",100,0,0.05) );
	      mon.addHistogram(new TH1F(lepStr+reg+"eopin",    ";E/p;Leptons",100,0,2) );
	      mon.addHistogram(new TH1F(lepStr+reg+"fbrem",    ";f_{brem};Leptons",100,0,2) );
	    }
	}
      else
	{
	  mon.addHistogram(new TH1F(lepStr+"nmatches", ";Muon matches;Leptons",15,0,15) );
	  mon.addHistogram(new TH1F(lepStr+"nmuonhits", ";Muon hits;Leptons",30,0,30) );
	}
      mon.addHistogram(new TH1F(lepStr+"d0",            ";d_{0};Leptons",50,0,0.05) );
      mon.addHistogram(new TH1F(lepStr+"dZ",            ";d_{Z};Leptons",50,0,0.1) );
      mon.addHistogram(new TH1F(lepStr+"trkchi2",       ";#chi^{2};Leptons",50,0,10) );
      mon.addHistogram(new TH1F(lepStr+"trkvalidpixel",  ";Valid pixel hits;Leptons",20,0,20) );
      mon.addHistogram(new TH1F(lepStr+"trkvalidtracker",  ";Valid tracker hits;Leptons",50,0,50) );
      mon.addHistogram(new TH1F(lepStr+"losthits",         ";Lost hits;Leptons",4,0,4) );
      mon.addHistogram(new TH1F(lepStr+"reliso",           ";RelIso;Leptons",50,0,2) );
      mon.addHistogram(new TH1F(lepStr+"reliso2011",        ";RelIso(#rho);Leptons",50,0,2) );
    }
  
  mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "rho",";#rho;Events",50,0,25) ); 
  mon.addHistogram( new TH1F( "rho25",";#rho(#eta<2.5);Events",50,0,25) ); 

  TString jetTypes[]={"pf","pfchs"};
  TString jetRegs[]={"TK","HEin","HEout","HF"};
  TString btagAlgos[]={"TCHE","CSV","JP"};
  Double_t btagAlgoMin[]={-5,-0.5,0};
  Double_t btagAlgoMax[]={15,1.5,5};
  for(size_t i=0; i<2; i++)
    {
      for(size_t ireg=0; ireg<4; ireg++)
	{
	  mon.addHistogram( new TH1F(jetRegs[ireg]+jetTypes[i]+"jetbeta"    , ";#beta;Events",50,0,1) );
	  mon.addHistogram( new TH1F(jetRegs[ireg]+jetTypes[i]+"jetbetastar", ";#beta *;Events",50,0,1) );
	  mon.addHistogram( new TH1F(jetRegs[ireg]+jetTypes[i]+"jetdrmean"  , ";<#Delta R>;Events",50,0,0.5) );
	  mon.addHistogram( new TH1F(jetRegs[ireg]+jetTypes[i]+"jetptd"     , ";p_{T}D [GeV/c];Events",50,0,1) );
	  mon.addHistogram( new TH1F(jetRegs[ireg]+jetTypes[i]+"jetptrms"   , ";RMS p_{T} [GeV/c];Events",50,0,0.5) );
	  mon.addHistogram( new TH1F(jetRegs[ireg]+jetTypes[i]+"jetmva"     , ";MVA;Events",50,-1,1) );
	}
      for(size_t ibtag=0; ibtag<3; ibtag++)
	{
	  mon.addHistogram( new TH1F(btagAlgos[ibtag]+"b"+jetTypes[i]+"jetstags",     ";b tags;Events",100,btagAlgoMin[ibtag],btagAlgoMax[ibtag]) );
	  mon.addHistogram( new TH1F(btagAlgos[ibtag]+"other"+jetTypes[i]+"jetstags", ";udscg tags;Events",100,btagAlgoMin[ibtag],btagAlgoMax[ibtag]) );
	  mon.addHistogram( new TH1F(btagAlgos[ibtag]+jetTypes[i]+"jetstags",         ";"+btagAlgos[ibtag]+";Events",100,btagAlgoMin[ibtag],btagAlgoMax[ibtag]) );
	  mon.addHistogram( new TH1F("n"+jetTypes[i]+"jetsbtags"+btagAlgos[ibtag],    ";b-tag multiplicity ("+btagAlgos[ibtag] +");Events",5,0,5) );
	}
      
      mon.addHistogram( new TH1F(jetTypes[i]+"jetpt"       , ";p_{T} [GeV/c];Events",50,0,250) );
      mon.addHistogram( new TH1F(jetTypes[i]+"jeteta"       , ";|#eta|;Events",25,0,5) );
      mon.addHistogram( new TH2F("n"+jetTypes[i]+"jetsvspu",          ";Pileup interactions;Jet multiplicity (p_{T}>30 GeV/c);Events",50,0,50,5,0,5) );
      mon.addHistogram( new TH2F("n"+jetTypes[i]+"jetstightvspu",     ";Pileup interactions;Jet multiplicity (p_{T}>30 GeV/c);Events",50,0,50,5,0,5) );
      mon.addHistogram( new TH2F("n"+jetTypes[i]+"jetspuidloosevspu", ";Pileup interactions;Jet multiplicity (p_{T}>30 GeV/c);Events",50,0,50,5,0,5) );
      mon.addHistogram( new TH2F("n"+jetTypes[i]+"jetspuidmediumvspu",";Pileup interactions;Jet multiplicity (p_{T}>30 GeV/c);Events",50,0,50,5,0,5) );  
    }


  //   mon.addHistogram( new TH1F("njets"       , ";Jet multiplicity (p_{T}>30 GeV/c);Events",5,0,5) );
  //   int jetids[]  ={JETID_LOOSE, JETID_TIGHT,JETID_CUTBASED_LOOSE, JETID_CUTBASED_MEDIUM, JETID_MIN_LOOSE, JETID_MIN_MEDIUM, JETID_OPT_LOOSE, JETID_OPT_MEDIUM};
  //   for(size_t ijetid=0; ijetid<sizeof(jetids)/sizeof(int); ijetid++)
  //     {
  //       char buf[15];
  //       sprintf(buf,"njets%d",jetids[ijetid]);
  //       mon.addHistogram( new TH1F(TString(buf), ";Jet multiplicity (p_{T}>30 GeV/c);Events",5,0,5) );
  //       mon.addHistogram( new TH1F(TString(buf)+"lowpu", ";Jet multiplicity (p_{T}>30 GeV/c);Events",5,0,5) );
  //       mon.addHistogram( new TH1F(TString(buf)+"cenpu", ";Jet multiplicity (p_{T}>30 GeV/c);Events",5,0,5) );
  //       mon.addHistogram( new TH1F(TString(buf)+"highpu", ";Jet multiplicity (p_{T}>30 GeV/c);Events",5,0,5) );
  
  //       char buf2[100];
  //       sprintf(buf2,"met_clusteredmet%d",jetids[ijetid]);
  //       mon.addHistogram( new TH1F( buf2  , ";clustered E_{T}^{miss};Events", 50,0,500) );
  //       mon.addHistogram( new TH2F( TString(buf2)+"_vspu" , ";Pileup events; clustered E_{T}^{miss};Events",50,0,50,50,0,500) );
  //     }
  
  //   mon.addHistogram( new TH1F ("nbtags", ";b-tag multiplicity; Events", 5,0,5) );  
  //   for(size_t ibin=1; ibin<=5; ibin++){
  //     TString label("");
  //     if(ibin==5) label +="#geq";
  //     else        label +="=";
  //     label += (ibin-1);
  //     mon.getHisto("njets")->GetXaxis()->SetBinLabel(ibin,label);
  //     mon.getHisto("nbtags")->GetXaxis()->SetBinLabel(ibin,label);
  //   }
  //   mon.addHistogram( new TH1F( "met_met"  , ";E_{T}^{miss};Events", 50,0,500) );
  //   mon.addHistogram( new TH1F( "met_min3Met"  , ";min(E_{T}^{miss},assoc-E_{T}^{miss},clustered-E_{T}^{miss});Events", 50,0,500) );
  //   mon.addHistogram( new TH1F( "met_redMet"  , ";red(E_{T}^{miss},clustered-E_{T}^{miss});Events", 50,0,500) );
  //   mon.addHistogram( new TH1F( "met_redMetL"  , ";red(E_{T}^{miss},clustered-E_{T}^{miss}) - longi.;Events", 50,-250,250) );
  //   mon.addHistogram( new TH1F( "met_redMetT"  , ";red(E_{T}^{miss},clustered-E_{T}^{miss}) - perp.;Events", 50,-250,250) );
  //   mon.addHistogram( new TH1F( "mt"  , ";M_{T};Events", 100,0,1000) );
  
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
  std::vector<float> mcPileupDistribution;
  if(isMC){
    TH1F* histo = (TH1F *) file->Get("evAnalyzer/h2zz/pileup");
    if(!histo)std::cout<<"pileup histogram is null!!!\n";
    for(int i=1;i<=histo->GetNbinsX();i++){mcPileupDistribution.push_back(histo->GetBinContent(i));}
    delete histo;
  }
  while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
  while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);
  /*
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
  */
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
      
      //load the event content from tree
      evSummaryHandler.getEntry(iev);
      ZZ2l2nuSummary_t &ev=evSummaryHandler.getEvent();
      PhysicsEvent_t phys=getPhysicsEventFrom(ev);

      //event category
      TString tag_cat;
      switch(ev.cat){
      case MUMU: tag_cat = "mumu";  break;
      case EE  : tag_cat = "ee";    break;
      default  : continue;
      }
      if(isMC && mctruthmode==1 && !isDYToLL(ev.mccat) ) continue;
      if(isMC && mctruthmode==2 && !isDYToTauTau(ev.mccat) ) continue;
      int eventSubCat  = eventCategoryInst.Get(phys);
      TString tag_subcat = eventCategoryInst.GetLabel(eventSubCat);
      
      //prepare the tag's vectors for histo filling
      std::vector<TString> tags_full(1,"all");
      // if(tag_subcat!="")tags_full.push_back(tag_cat + tag_subcat);
      
      //pileup weight
      float weight = 1.0;
      double TotalWeight_plus = 1.0;
      double TotalWeight_minus = 1.0;
      if(isMC){
//         weight = LumiWeights->weight(ev.nvtx);
//         TotalWeight_plus = PShiftUp->ShiftWeight(ev.nvtx);
//         TotalWeight_minus = PShiftDown->ShiftWeight(ev.nvtx);
        weight = LumiWeights->weight(ev.ngenITpu);
        TotalWeight_plus = PShiftUp->ShiftWeight(ev.ngenITpu);
        TotalWeight_minus = PShiftDown->ShiftWeight(ev.ngenITpu);
      }
      Hcutflow->Fill(1,1);
      Hcutflow->Fill(2,weight);
      Hcutflow->Fill(3,weight*TotalWeight_minus);
      Hcutflow->Fill(4,weight*TotalWeight_plus);
      
      //
      // LEPTON ANALYSIS
      //
      LorentzVector lep1=phys.leptons[0];
      LorentzVector lep2=phys.leptons[1];
      bool passIdAndIso(true);
      LorentzVector zll(lep1+lep2);
      bool passZmass(fabs(zll.mass()-91)<15);
      bool passZpt(zll.pt()>55);

      //check alternative selections for the dilepton
      for(int ilep=0; ilep<2; ilep++)
	{
	  //lepton type
	  TString lepStr((fabs(phys.leptons[ilep].id)==13?"mu":"e"));
	  
	  //generator level matching
	  int matchid(0);
	  LorentzVector genP4(0,0,0,0);
	  for(size_t igl=0;igl<phys.genleptons.size(); igl++) 
	    {
	      if(deltaR(phys.genleptons[igl],phys.leptons[ilep])>0.1) continue;
	      genP4=phys.genleptons[igl];
	      matchid=phys.genleptons[igl].id;
	    }
	  
	  //id and isolation
	  int lpid=phys.leptons[ilep].pid;
	  float relIso2011    = phys.leptons[ilep].relIsoRho(ev.rho);
	  float relIso = (lepStr=="mu") ? 
	    phys.leptons[ilep].muPFRelIsoCorrected2012(ev.rho25Neut):
	    phys.leptons[ilep].ePFRelIsoCorrected2012(ev.rho);
	  std::vector<int> passIds;
	  std::map<int,bool> passIsos;
	  bool hasGoodId(false);
	  if(fabs(phys.leptons[ilep].id)==13)
	    {
	      if( hasObjectId(ev.mn_idbits[lpid], MID_LOOSE) )    { passIds.push_back(0); passIsos[0]=(relIso<0.2); }
	      if( hasObjectId(ev.mn_idbits[lpid], MID_TIGHT) )    { passIds.push_back(1); passIsos[1]=(relIso<0.2); hasGoodId=true; }
	      if( hasObjectId(ev.mn_idbits[lpid], MID_VBTF2011) ) { passIds.push_back(2); passIsos[2]=(relIso2011<0.15); }
	      if( hasObjectId(ev.mn_idbits[lpid], MID_SOFT) )     { passIds.push_back(3); passIsos[3]=true;}
	    }
	  else
	    {
	      if(hasObjectId(ev.en_idbits[lpid], EID_LOOSE))      { passIds.push_back(0); passIsos[0]=(relIso<0.15); }
	      if(hasObjectId(ev.en_idbits[lpid], EID_MEDIUM))     { passIds.push_back(1); passIsos[1]=(relIso<0.15); hasGoodId=true; }
	      if(hasObjectId(ev.en_idbits[lpid], EID_VBTF2011))   { passIds.push_back(2); passIsos[2]=(relIso2011<0.10);}
	      if(hasObjectId(ev.en_idbits[lpid], EID_VETO))       { passIds.push_back(3); passIsos[3]=(relIso<0.15); }
	    }
	  if(!hasGoodId)        passIdAndIso=false;
	  else if(!passIsos[1]) passIdAndIso=false;     

	  //fill control histograms (constrained to the Z mass)
	  if(passZmass)
	    {
	      if(matchid!=0)
		{
		  mon.fillHisto(lepStr+"genpt",tags_full, genP4.pt(), weight);
		  mon.fillHisto(lepStr+"geneta",tags_full,genP4.eta(), weight);
		  mon.fillHisto(lepStr+"genpu",tags_full,ev.ngenITpu, weight);
		  for(size_t iid=0; iid<passIds.size(); iid++)
		    {
		      TString idStr(lepStr);  idStr += passIds[iid];
		      mon.fillHisto(idStr+"pt",tags_full, genP4.pt(), weight);
		      mon.fillHisto(idStr+"eta",tags_full,genP4.eta(), weight);
		      mon.fillHisto(idStr+"pu",tags_full,ev.ngenITpu, weight);
		      if(!passIsos[ passIds[iid] ]) continue;
		      mon.fillHisto(idStr+"isopt",tags_full, genP4.pt(), weight);
		      mon.fillHisto(idStr+"isoeta",tags_full,genP4.eta(), weight);
		      mon.fillHisto(idStr+"isopu",tags_full,ev.ngenITpu, weight);
		    }
		}
	  
	      mon.fillHisto(lepStr+"d0",              tags_full,fabs(phys.leptons[ilep].d0),weight);
	      mon.fillHisto(lepStr+"dZ",              tags_full,fabs(phys.leptons[ilep].dZ),weight);
	      mon.fillHisto(lepStr+"trkchi2",         tags_full,fabs(phys.leptons[ilep].trkchi2),weight);
	      mon.fillHisto(lepStr+"trkvalidpixel",   tags_full,fabs(phys.leptons[ilep].trkValidPixelHits),weight);
	      mon.fillHisto(lepStr+"trkvalidtracker", tags_full,fabs(phys.leptons[ilep].trkLostInnerHits),weight);
	      mon.fillHisto(lepStr+"losthits",        tags_full,fabs(phys.leptons[ilep].trkLostInnerHits),weight);

	      if(lepStr=="e")
		{
		  TString reg="ee";
		  if(fabs(phys.leptons[ilep].eta())<1.442) reg="eb";
		  mon.fillHisto(lepStr+reg+"detain",  tags_full,fabs(ev.en_detain[lpid]),weight);
		  mon.fillHisto(lepStr+reg+"dphiin",  tags_full,fabs(ev.en_dphiin[lpid]),weight);
		  mon.fillHisto(lepStr+reg+"sihih",   tags_full,fabs(ev.en_sihih[lpid]),weight);
		  mon.fillHisto(lepStr+reg+"sipip",   tags_full,fabs(ev.en_sipip[lpid]),weight);
		  mon.fillHisto(lepStr+reg+"r9",      tags_full,fabs(ev.en_r9[lpid]),weight);
		  mon.fillHisto(lepStr+reg+"hoe",     tags_full,fabs(ev.en_hoe[lpid]),weight);
		  mon.fillHisto(lepStr+reg+"eopin",   tags_full,fabs(ev.en_eopin[lpid]),weight);
		  mon.fillHisto(lepStr+reg+"fbrem",   tags_full,fabs(ev.en_fbrem[lpid]),weight);
		  mon.fillHisto(lepStr+reg+"ooemoop", tags_full,fabs(ev.en_ooemoop[lpid]),weight);
		}
	      else
		{
		  mon.fillHisto(lepStr+"nmuonhits",  tags_full,fabs(ev.mn_validMuonHits[lpid]),weight);
		  mon.fillHisto(lepStr+"nmatches",   tags_full,fabs(ev.mn_nMatches[lpid]),weight);
		}
	      
	      if(hasGoodId)
		{
		  mon.fillHisto(lepStr+"reliso",     tags_full, relIso,   weight);
		  mon.fillHisto(lepStr+"reliso2011", tags_full, relIso2011, weight);
		}
	    }
	}
      
      //select dilepton

      mon.fillHisto  ("eventflow",tags_full,0,weight);

      if(!passIdAndIso) continue;
      tags_full.push_back(tag_cat);
      mon.fillHisto  ("eventflow",tags_full,1,weight);
      mon.fillHisto("zmass",       tags_full, zll.mass(), weight);  

      if(!passZmass) continue;
      mon.fillHisto  ("eventflow",tags_full,2,weight);
      mon.fillHisto("zpt"      ,   tags_full,zll.pt()     ,weight);      
 
      if(!passZpt) continue;
      mon.fillHisto  ("eventflow",tags_full,3,weight);

      //analyze dilepton kinematics
      LorentzVector leadingLep(phys.leptons[0].pt()>phys.leptons[1].pt() ? phys.leptons[0]: phys.leptons[1]);
      LorentzVector trailerLep(phys.leptons[0].pt()>phys.leptons[1].pt() ? phys.leptons[1]: phys.leptons[0]);
      mon.fillHisto("leadeta"     ,   tags_full,leadingLep.eta()   ,weight);
      mon.fillHisto("leadpt"      ,   tags_full,leadingLep.pt()     ,weight);
      mon.fillHisto("trailereta"     ,   tags_full,trailerLep.eta()   ,weight);
      mon.fillHisto("trailerpt"      ,   tags_full,trailerLep.pt()     ,weight);
      mon.fillHisto("zeta"     ,   tags_full,zll.eta()   ,weight);
      mon.fillHisto("nvtx"        ,tags_full, ev.nvtx,  weight);
      mon.fillHisto("rho"         ,tags_full, ev.rho,   weight);
      mon.fillHisto("rho25"       ,tags_full, ev.rho25Neut, weight);
    
      //
      //JET ANALYSIS
      //

      //std PF
      PhysicsObjectJetCollection aJetsP4;
      std::vector<double> genAJetsPt;
      for(size_t ijet=0; ijet<phys.ajets.size(); ijet++)
	{
	  aJetsP4.push_back( phys.ajets[ijet] );
	  genAJetsPt.push_back( phys.ajets[ijet].genPt);
	}
      LorentzVectorCollection zvvs;
      std::vector<PhysicsObjectJetCollection> variedAJets;
      METUtils::computeVariation(aJetsP4, phys.met[0], variedAJets, zvvs, &jecUnc);
      PhysicsObjectJetCollection aJets=variedAJets[0];
      LorentzVector aClusteredMetP4(zll); aClusteredMetP4 *= -1;
      int nAJetsLoose(0), nAJetsTight(0), nAJetsPUIdLoose(0), nAJetsPUIdMedium(0);
      int nABtags[3]={0,0,0};
      for(size_t ijet=0; ijet<aJets.size(); ijet++) 
	{
	  aClusteredMetP4 -= aJets[ijet];
	  mon.fillHisto("pfjetpt",  tags_full, aJets[ijet].pt(),weight);
	  mon.fillHisto("pfjeteta",  tags_full, fabs(aJets[ijet].eta()),weight);

	  if(aJets[ijet].pt()<20) continue;
	  if(fabs(aJets[ijet].eta())<2.5) 
	    {
	      nABtags[0] += (phys.ajets[ijet].btag1>2.0);
	      nABtags[1] += (phys.ajets[ijet].btag2>0.244);
	      nABtags[2] += (phys.ajets[ijet].btag3>0.275);
	      Float_t ijetbtags[]={phys.ajets[ijet].btag1, phys.ajets[ijet].btag2, phys.ajets[ijet].btag3};
	      for(size_t ibtag=0; ibtag<3; ibtag++)
		{
		  if(isMC)
		    {
		      if(fabs(phys.ajets[ijet].flavid)==5) mon.fillHisto(btagAlgos[ibtag]+"bpfjetstags",     tags_full, ijetbtags[ibtag], weight);
		      else                                 mon.fillHisto(btagAlgos[ibtag]+"otherpfjetstags", tags_full, ijetbtags[ibtag], weight);
		    }
		  mon.fillHisto(btagAlgos[ibtag]+"pfjetstags",  tags_full, ijetbtags[ibtag],weight);
		}
	    }
	  if(aJets[ijet].pt()<30 ) continue;
	  nAJetsLoose      += hasObjectId(phys.ajets[ijet].pid,JETID_LOOSE);
	  nAJetsTight      += hasObjectId(phys.ajets[ijet].pid,JETID_TIGHT);
	  nAJetsPUIdLoose  += hasObjectId(phys.ajets[ijet].pid,JETID_MIN_LOOSE);
	  nAJetsPUIdMedium += hasObjectId(phys.ajets[ijet].pid,JETID_MIN_MEDIUM);

	  TString reg=getJetRegion(phys.ajets[ijet].eta());
	  mon.fillHisto(reg+"pfjetbeta",     tags_full,phys.ajets[ijet].beta,     weight);
	  mon.fillHisto(reg+"pfjetbetastar", tags_full,phys.ajets[ijet].betaStar, weight);
	  mon.fillHisto(reg+"pfjetdrmean",   tags_full,phys.ajets[ijet].dRMean,   weight);
	  mon.fillHisto(reg+"pfjetptd",      tags_full,phys.ajets[ijet].ptD,      weight);
	  mon.fillHisto(reg+"pfjetptrms",    tags_full,phys.ajets[ijet].ptRMS,    weight);
	  mon.fillHisto(reg+"pfjetmva",      tags_full,phys.ajets[ijet].pumva,    weight);
	}
      METUtils::stRedMET aRedMetOut; 
      LorentzVector aRedMet=METUtils::redMET(METUtils::INDEPENDENTLYMINIMIZED, lep1, 0, lep2, 0, aClusteredMetP4, zvvs[0],false,&aRedMetOut);
      double aRedMetL=aRedMetOut.redMET_l;
      double aRedMetT=aRedMetOut.redMET_t;
      double aMT=METUtils::transverseMass(zll,zvvs[0],true);
      
      //AK5 CHS
      PhysicsObjectJetCollection jetsP4;
      std::vector<double> genJetsPt;
      for(size_t ijet=0; ijet<phys.jets.size(); ijet++)
	{
	  jetsP4.push_back( phys.jets[ijet] );
	  genJetsPt.push_back( phys.jets[ijet].genPt);
	}
      LorentzVectorCollection zvvsCHS;
      std::vector<PhysicsObjectJetCollection> variedJets;
      METUtils::computeVariation(jetsP4, phys.met[0], variedJets, zvvsCHS, &jecUnc);
      PhysicsObjectJetCollection jets=variedJets[0];
      LorentzVector clusteredMetP4(zll); clusteredMetP4 *= -1;
      int nJetsLoose(0), nJetsTight(0);
      int nBtags[3]={0,0,0};
      for(size_t ijet=0; ijet<jets.size(); ijet++) 
	{
	  clusteredMetP4 -= jets[ijet];
	  mon.fillHisto("pfchsjetpt",  tags_full, jets[ijet].pt(),weight);
	  mon.fillHisto("pfchsjeteta",  tags_full, fabs(jets[ijet].eta()),weight);

	  if(jets[ijet].pt()<20) continue;
	  if(fabs(jets[ijet].eta())<2.5) 
	    {
	      nBtags[0] += (phys.jets[ijet].btag1>2.0);
	      nBtags[1] += (phys.jets[ijet].btag2>0.244);
	      nBtags[2] += (phys.jets[ijet].btag3>0.275);
	      Float_t ijetbtags[]={phys.jets[ijet].btag1, phys.jets[ijet].btag2, phys.jets[ijet].btag3};
	      for(size_t ibtag=0; ibtag<3; ibtag++)
		{
		  if(isMC)
		    {
		      if(fabs(phys.jets[ijet].flavid)==5) mon.fillHisto(btagAlgos[ibtag]+"bpfchsjetstags",     tags_full, ijetbtags[ibtag], weight);
		      else                               mon.fillHisto(btagAlgos[ibtag]+"otherpfchsjetstags", tags_full, ijetbtags[ibtag], weight);
		    }
		  mon.fillHisto(btagAlgos[ibtag]+"pfchsjetstags",  tags_full, ijetbtags[ibtag],weight);
		}
	    }
	  if(jets[ijet].pt()<30 ) continue;
	  nJetsLoose      += hasObjectId(phys.jets[ijet].pid,JETID_LOOSE);
	  nJetsTight      += hasObjectId(phys.jets[ijet].pid,JETID_TIGHT);

	  TString reg=getJetRegion(phys.jets[ijet].eta());
	  mon.fillHisto(reg+"pfchsjetbeta",     tags_full,phys.jets[ijet].beta,     weight);
	  mon.fillHisto(reg+"pfchsjetbetastar", tags_full,phys.jets[ijet].betaStar, weight);
	  mon.fillHisto(reg+"pfchsjetdrmean",   tags_full,phys.jets[ijet].dRMean,   weight);
	  mon.fillHisto(reg+"pfchsjetptd",      tags_full,phys.jets[ijet].ptD,      weight);
	  mon.fillHisto(reg+"pfchsjetptrms",    tags_full,phys.jets[ijet].ptRMS,    weight);
	  mon.fillHisto(reg+"pfchsjetmva",      tags_full,phys.jets[ijet].pumva,    weight);
	}
      METUtils::stRedMET redMetOut; 
      LorentzVector redMet=METUtils::redMET(METUtils::INDEPENDENTLYMINIMIZED, lep1, 0, lep2, 0, clusteredMetP4, zvvsCHS[0],false,&redMetOut);
      double redMetL=redMetOut.redMET_l;
      double redMetT=redMetOut.redMET_t;
      double mT=METUtils::transverseMass(zll,zvvsCHS[0],true);
      

      for(size_t ibtagalgo=0; ibtagalgo<3; ibtagalgo++)
	{
	  mon.fillHisto("npfjetsbtags"+btagAlgos[ibtagalgo],    tags_full, nABtags[ibtagalgo], weight);
	  mon.fillHisto("npfchsjetsbtags"+btagAlgos[ibtagalgo], tags_full, nBtags[ibtagalgo], weight);
	}
      if(nABtags[2]>0) continue;
      mon.fillHisto  ("eventflow",tags_full,4,weight);

      mon.fillHisto("npfjetsvspu",          tags_full, ev.ngenITpu, nAJetsLoose,weight);
      mon.fillHisto("npfjetstightvspu",     tags_full, ev.ngenITpu, nAJetsTight,weight);
      mon.fillHisto("npfjetspuidloosevspu", tags_full, ev.ngenITpu, nAJetsPUIdLoose,weight);
      mon.fillHisto("npfjetspuidmediumvspu",tags_full, ev.ngenITpu, nAJetsPUIdMedium,weight);

      mon.fillHisto("npfchsjetsvspu",          tags_full, ev.ngenITpu, nJetsLoose,weight);
      mon.fillHisto("npfchsjetstightvspu",     tags_full, ev.ngenITpu, nJetsTight,weight);


      //
      // 3rd LEPTON ANALYSIS
      //
      for(size_t ilep=2; ilep<phys.leptons.size(); ilep++)
	{
	  //lepton type
	  TString lepStr((fabs(phys.leptons[ilep].id)==13?"mu":"e"));
	  
	  //generator level matching
	  int matchid(0);
	  LorentzVector genP4(0,0,0,0);
	  for(size_t igl=0;igl<phys.genleptons.size(); igl++) 
	    {
	      if(deltaR(phys.genleptons[igl],phys.leptons[ilep])>0.1) continue;
	      genP4=phys.genleptons[igl];
	      matchid=phys.genleptons[igl].id;
	    }

	  mon.fillHisto("thirdleptoneta",   tags_full,phys.leptons[ilep].eta()   ,weight);
	  mon.fillHisto("thirdleptonpt" ,   tags_full,phys.leptons[ilep].pt()     ,weight);
	  mon.fillHisto("thirdleptonmt" ,   tags_full,METUtils::transverseMass(phys.leptons[ilep],zvvs[0],false),weight);
	  if(matchid!=0)
	    {
	      
	    }
	  

	}

      /*
      
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
	bool passId(idPass[0]);
	bool passIso(phys.leptons[0].relIsoRho(ev.rho)<0.15 && phys.leptons[1].relIsoRho(ev.rho)<0.15);
	bool isZsideBand( (zmass>40 && zmass<70) || (zmass>110 && zmass<200));
	bool isZsideBandPlus( (zmass>110 && zmass<200));
	bool passDphijmet(mindphijmet>0.5); 
	bool passZpt(zpt>55);
	bool pass3dLeptonVeto(true); int nExtraLep(0); for(unsigned int i=2;i<phys.leptons.size();i++) { if(phys.leptons[i].pt()>10){ nExtraLep++; pass3dLeptonVeto=false;} }
	bool passBveto(nbtags==0);
	bool passBaseMet(zvv.pt()>70);
     

	//match the leptons
	Float_t genPt[]={-1,-1};
	Float_t genEta[]={0,0,};
	for(size_t igl=0;igl<phys.genleptons.size(); igl++) 
	  {
	    for(size_t irec=0; irec<2; irec++) 
	      if(deltaR(phys.genleptons[igl],phys.leptons[irec])<0.1) {
		genPt[irec]=phys.genleptons[igl].pt(); 
		genEta[irec]=phys.genleptons[igl].eta(); 
	      }
	  }

    
	//##############################################  
	//########         GENERAL PLOTS        ########                                                                                                                  
	//##############################################  

	//N-1 PLOTS FOR PRESELECTION
	if(ivar==0){

	  if(passId && passIso) mon.fillHisto  ("zmass"    ,tags_full,zmass,iweight);
	  

	      {
		TString pre("ee");
		if(ev.cat==MUMU) pre="mumu";
		for(size_t idPassCtr=0; idPassCtr<4; idPassCtr++)
		  if(idPass[idPassCtr])  
		    mon.fillHisto("z"+pre+"selafterid",tags_full,idPassCtr,weight);

		pre="e";
		Int_t idsToTest[]={EID_LOOSE,EID_MEDIUM,EID_TIGHT,EID_VBTF2011};
		if(ev.cat==MUMU)
		  {
		    pre="mu";
		    idsToTest[0]=MID_LOOSE;
		    idsToTest[1]=MID_TIGHT;
		    idsToTest[2]=MID_HIGHPT;
		    idsToTest[3]=MID_VBTF2011;
		  }
		for(int ilep=0; ilep<2; ilep++)
		  {

		    if(genPt[ilep]>0)
		      {
			mon.fillHisto(pre+"genpt",tags_full, genPt[ilep], weight);
			mon.fillHisto(pre+"geneta",tags_full, genEta[ilep], weight);
			for(int iid=0; iid<4; iid++)
			  {
			    if(hasObjectId(lid[ilep], idsToTest[iid]))
			      { 
				TString idctr(""); idctr+=iid;
				mon.fillHisto(pre+idctr+"pt",tags_full, genPt[ilep], weight);
				mon.fillHisto(pre+idctr+"eta",tags_full, genEta[ilep], weight);
			      }
			  }
		      }
		  }
	      }
	    
	    for(size_t ilep=0; ilep<2; ilep++)
	      {
	      }		
	    
	    mon.fillHisto("eventflow",tags_full,1,iweight);
	    
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
		
		mon.fillHisto("zeta"     ,tags_full,zeta   ,iweight);
		mon.fillHisto("zpt"      ,tags_full,zpt     ,iweight);
		mon.fillHisto("nvtx"     ,tags_full,ev.nvtx,iweight);
		mon.fillHisto("rho"      ,tags_full,ev.rho,iweight);
		mon.fillHisto("rho25"    ,tags_full,ev.rho25,iweight);
		mon.fill2DHisto("rhovspu"      ,tags_full,ev.ngenITpu,ev.rho,iweight);
		mon.fill2DHisto("rho25vspu"    ,tags_full,ev.ngenITpu,ev.rho25,iweight);

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


			      
			      mon.fillHisto("jetpt",tags_full,phys.ajets[ijet].pt(),iweight);
			      mon.fillHisto("jeteta",tags_full,fabs(phys.ajets[ijet].eta()),iweight);
			    }
			  
			  
			  if(passDphijmet){
			    if(runBlinded && !isMC && evSummaryHandler.hasSpoilerAlert(!isMC,tags_full)) continue;
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
      */
  
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





