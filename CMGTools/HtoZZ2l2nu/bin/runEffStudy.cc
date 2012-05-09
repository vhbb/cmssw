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
  cout << hasObjectId(id,MID_PF) << " " << hasObjectId(id,MID_GLOBAL) << endl;
  if(!hasObjectId(id,MID_PF) )    return 0;
  if(!hasObjectId(id,MID_GLOBAL)) return 1;
  if(chi2>10)                     return 2;
  if(muonHits==0)                 return 3;
  if(muonMatches<=1)              return 4;
  if(fabs(d0)>0.2)                return 5;
  if(fabs(dZ)>0.5)                return 6;
  if(pixelHits==0)                return 7;
  if(trkLayers<=5)                return 8;
  return 9;
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

  mon.addHistogram( new TH1F( "pt", ";p_{T}^{l};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "eta", ";#eta^{l};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "zpt", ";p_{T}^{ll};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "zeta", ";#eta^{ll};Events", 50,-10,10) );
  mon.addHistogram( new TH1F( "zmass", ";M^{ll};Events", 100,40,250) );

  //lepton control
  for(size_t ilep=0; ilep<2; ilep++)
    {
      TString lepStr(ilep==0? "mu" :"e");
      mon.addHistogram(new TH1F(lepStr+"genpt",   ";p_{T} [GeV/c];Leptons",25,0,250) );
      mon.addHistogram(new TH1F(lepStr+"geneta",   ";#eta;Leptons",20,-5,5) );
      for(int iid=0; iid<10; iid++)
	{
	  TString idctr(""); idctr+=iid;
	  mon.addHistogram(new TH1F(lepStr+idctr+"pt",   ";p_{T} [GeV/c];Leptons",25,0,250) );
	  mon.addHistogram(new TH1F(lepStr+idctr+"eta",   ";#eta;Leptons",20,-5,5) );
	}
      if(ilep==1)
	{
	  for(size_t ireg=0; ireg<2; ireg++)
	    { 
	      TString reg(ireg==0?"eb":"ee");
	      mon.addHistogram(new TH1F(lepStr+reg+"detain",   ";#Delta#eta_{in};Leptons",20,0,0.01) );
	      mon.addHistogram(new TH1F(lepStr+reg+"dphiin",   ";#Delta#phi_{in};Leptons",20,0,0.1) );
	      mon.addHistogram(new TH1F(lepStr+reg+"sihih",    ";#sigma_{i#eta i#eta};Leptons",20,0,0.05) );
	      mon.addHistogram(new TH1F(lepStr+reg+"sipip",    ";#sigma_{i#phi i#phi};Leptons",20,0,0.05) );
	      mon.addHistogram(new TH1F(lepStr+reg+"r9",       ";R_{9};Leptons",20,0,1.) );
	      mon.addHistogram(new TH1F(lepStr+reg+"hoe",      ";h/e;Leptons",20,0,0.2) );
	      mon.addHistogram(new TH1F(lepStr+reg+"ooemoop",  ";1/E-1/p;Leptons",50,0,0.05) );
	      mon.addHistogram(new TH1F(lepStr+reg+"eopin",    ";E/p;Leptons",50,0,2) );
	      mon.addHistogram(new TH1F(lepStr+reg+"fbrem",    ";f_{brem};Leptons",50,0,2) );
	    }
	}
      else
	{
	  mon.addHistogram(new TH1F(lepStr+"nmatches", ";Muon matches;Leptons",15,0,15) );
	  mon.addHistogram(new TH1F(lepStr+"nmuonhits", ";Muon hits;Leptons",30,0,30) );
	}
      mon.addHistogram(new TH1F(lepStr+"d0",            ";d_{0};Leptons",50,0,0.05) );
      mon.addHistogram(new TH2F(lepStr+"d0vsreliso",    ";d_{0};Leptons",50,0,0.05,50,0,2) );
      mon.addHistogram(new TH1F(lepStr+"dZ",            ";d_{Z};Leptons",50,0,0.1) );
      mon.addHistogram(new TH1F(lepStr+"trkchi2",       ";#chi^{2};Leptons",50,0,10) );
      mon.addHistogram(new TH1F(lepStr+"trkvalidpixel",  ";Valid pixel hits;Leptons",20,0,20) );
      mon.addHistogram(new TH1F(lepStr+"trkvalidtracker",  ";Valid tracker hits;Leptons",50,0,50) );
      mon.addHistogram(new TH1F(lepStr+"losthits",         ";Lost hits;Leptons",4,0,4) );
    //   mon.addHistogram(new TH1F(lepStr"reliso",           ";RelIso;Leptons",50,0,2) );
//       mon.addHistogram(new TH1F(lepStr"relisorho",        ";RelIso(#rho);Leptons",50,0,2) );
//       mon.addHistogram(new TH1F(lepStr"pfreliso",         ";PFRelIso;Leptons",50,0,2) );
//       mon.addHistogram(new TH1F(lepStr"pfrelisodbeta",    ";PFRelIso(#Delta #beta);Leptons",50,0,2) );
//       mon.addHistogram(new TH1F(lepStr"ecaliso",          ";ECALIso;Leptons",50,0,5) );
//       mon.addHistogram(new TH1F(lepStr"hcaliso",          ";HCALIso;Leptons",50,0,5) );
//       mon.addHistogram(new TH1F(lepStr"trkiso",           ";TrkIso;Leptons",50,0,5) );
//       mon.addHistogram(new TH1F(lepStr"giso",             ";GammaIso;Leptons",50,0,5) );
//       mon.addHistogram(new TH1F(lepStr"chiso",            ";ChHadronIso;Leptons",50,0,5) );
//       mon.addHistogram(new TH1F(lepStr"nhiso",            ";NeutralHadron;Leptons",50,0,5) );
//       mon.addHistogram(new TH2F(lepStr"relisovspu",         ";RelIso;Pileup;Leptons",50,0,50,50,0,2) );
//       mon.addHistogram(new TH2F(lepStr"relisorhovspu",      ";RelIso(#rho);Pileup;Leptons",50,0,50,50,0,2) ); 
//       mon.addHistogram(new TH2F(lepStr"pfrelisovspu",       ";PFRelIso;Pileup;Leptons",50,0,50,50,0,2) ); 
//       mon.addHistogram(new TH2F(lepStr"pfrelisodbetavspu",  ";PFRelIso(#Delta #beta);Pileup;Leptons",50,0,50,50,0,2) );  
//       mon.addHistogram(new TH2F(lepStr"ecalisovspu",        ";ECALIso;Pileup;Leptons",50,0,50,50,0,5) ); 
//       mon.addHistogram(new TH2F(lepStr"hcalisovspu",        ";HCALIso;Pileup;Leptons",50,0,50,50,0,5) ); 
//       mon.addHistogram(new TH2F(lepStr"trkisovspu",         ";TrkIso;Pileup;Leptons",50,0,50,50,0,5) ); 
//       mon.addHistogram(new TH2F(lepStr"gisovspu",           ";GammaIso;Pileup;Leptons",50,0,50,50,0,5) ); 
//       mon.addHistogram(new TH2F(lepStr"chisovspu",          ";ChHadronIso;Pileup;Leptons",50,0,50,50,0,5) ); 
//       mon.addHistogram(new TH2F(lepStr"nhisovspu",          ";NeutralHadronIso;Pileup;Leptons",50,0,50,50,0,5) ); 
    }

  mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "rho",";#rho;Events",50,0,25) ); 
  mon.addHistogram( new TH1F( "rho25",";#rho(#eta<2.5);Events",50,0,25) ); 
  mon.addHistogram(new TH2F( "rhovspu",          ";#rho;Pileup;Leptons",50,0,50,50,0,10) ); 
  mon.addHistogram(new TH2F( "rho25vspu",          ";#rho(|#eta|<2.5);Pileup;Leptons",50,0,50,50,0,10) ); 

//   TString jetRegs[]={"TK","HEin","HEout","HF"};
//   for(size_t ireg=0; ireg<4; ireg++)
//     {
//       mon.addHistogram( new TH1F(jetRegs[ireg]+"jetbeta"       , ";#beta;Events",50,0,1) );
//       mon.addHistogram( new TH1F(jetRegs[ireg]+"jetbetastar"       , ";#beta *;Events",50,0,1) );
//       mon.addHistogram( new TH1F(jetRegs[ireg]+"jetdrmean"       , ";<#Delta R>;Events",50,0,0.5) );
//       mon.addHistogram( new TH1F(jetRegs[ireg]+"jetptd"       , ";p_{T}D [GeV/c];Events",50,0,1) );
//       mon.addHistogram( new TH1F(jetRegs[ireg]+"jetptrms"       , ";RMS p_{T} [GeV/c];Events",50,0,0.5) );
//       mon.addHistogram( new TH1F(jetRegs[ireg]+"jetmva"      , ";MVA;Events",50,-1,1) );
//     }
//   mon.addHistogram( new TH1F("jetpt"       , ";p_{T} [GeV/c];Events",50,0,250) );
//   mon.addHistogram( new TH1F("jeteta"       , ";|#eta|;Events",25,0,5) );
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

  //   std::vector<float> mcPileupDistribution(dataPileupDistribution.size(),0);
  //   if(isMC)
  //     {
  //       for( int iev=evStart; iev<evEnd; iev++)
  // 	{
  // 	  evSummaryHandler.getEntry(iev);
  // 	  int nvtx=evSummaryHandler.evSummary_.nvtx;
  // 	  if(nvtx>int(mcPileupDistribution.size()-1)) continue;
  // 	  mcPileupDistribution[nvtx] += 1;
  //       }
  //     }

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
        weight = LumiWeights->weight(ev.nvtx);
        TotalWeight_plus = PShiftUp->ShiftWeight(ev.nvtx);
        TotalWeight_minus = PShiftDown->ShiftWeight(ev.nvtx);
      }
      Hcutflow->Fill(1,1);
      Hcutflow->Fill(2,weight);
      Hcutflow->Fill(3,weight*TotalWeight_minus);
      Hcutflow->Fill(4,weight*TotalWeight_plus);
      
      //analyze the dilepton
      LorentzVector lep1=phys.leptons[0];
      LorentzVector lep2=phys.leptons[1];
      LorentzVector zll(lep1+lep2);
      bool passZmass(fabs(zll.mass()-91)<15);

      if(!passZmass) continue;
      mon.fillHisto  ("eventflow",tags_full,0,weight);
      
      //match the leptons
      std::vector<LorentzVector> lGenP4(2,LorentzVector(0,0,0,0));
      for(size_t igl=0;igl<phys.genleptons.size(); igl++) 
	{
	  for(size_t irec=0; irec<2; irec++) 
	    {
	      if(deltaR(phys.genleptons[igl],phys.leptons[irec])>0.1) continue;
	      lGenP4[irec]=phys.genleptons[igl];
	    }
	}

      //check id and isolation
      std::vector< std::vector<int> > lPassIds(2);

      //muon specific
      TString lepStr("mu");
      if (ev.cat==MUMU)
	{
	  int m1pid=phys.leptons[0].pid;
	  int m1Level = muonIdLevel(ev.mn_idbits[m1pid], ev.l1_trkchi2, ev.mn_validMuonHits[m1pid], ev.mn_nMatches[m1pid], 
				    ev.l1_d0, ev.l1_dZ, ev.l1_trkValidPixelHits, ev.mn_trkLayersWithMeasurement[m1pid]);
	  // cout << m1pid << " " << hex << ev.mn_idbits[m1pid] << dec << " " <<  m1Level << endl;
	  for(int iid=1; iid<=m1Level; iid++) lPassIds[0].push_back(iid);
	  if( hasObjectId(ev.mn_idbits[m1pid], MID_VBTF2011) ) lPassIds[0].push_back(0);
	  
	  int m2pid=phys.leptons[1].pid;
	  int m2Level = muonIdLevel(ev.mn_idbits[m2pid], ev.l2_trkchi2, ev.mn_validMuonHits[m2pid], ev.mn_nMatches[m2pid], 
				    ev.l2_d0, ev.l2_dZ, ev.l2_trkValidPixelHits, ev.mn_trkLayersWithMeasurement[m2pid]);
	  for(int iid=1; iid<=m2Level; iid++) lPassIds[1].push_back(iid);
	  if( hasObjectId(ev.mn_idbits[m2pid], MID_VBTF2011) ) lPassIds[1].push_back(0);	  
	}

      //electron specific
      else if(ev.cat==EE)
	{
	  lepStr="e";
	  int e1pid=phys.leptons[0].pid;
	  if(hasObjectId(e1pid, EID_VBTF2011)) lPassIds[0].push_back(0);
	  if(hasObjectId(e1pid, EID_MEDIUM))   lPassIds[0].push_back(1);
	  if(hasObjectId(e1pid, EID_TIGHT))    lPassIds[0].push_back(2);

	  int e2pid=phys.leptons[1].pid;
	  if(hasObjectId(e2pid, EID_VBTF2011)) lPassIds[1].push_back(0);
	  if(hasObjectId(e2pid, EID_MEDIUM))   lPassIds[1].push_back(1);
	  if(hasObjectId(e2pid, EID_TIGHT))    lPassIds[1].push_back(2);
	}

      //fill histograms
      for(int ilep=0; ilep<2; ilep++)
	{
	  if(lGenP4[ilep].pt()>0)
	    {
	      mon.fillHisto(lepStr+"genpt",tags_full, lGenP4[ilep].pt(), weight);
	      mon.fillHisto(lepStr+"geneta",tags_full,lGenP4[ilep].eta(), weight);
	      for(size_t iid=0; iid<lPassIds[ilep].size(); iid++)
		{
		  TString idStr(""); idStr += lPassIds[ilep][iid];
		  mon.fillHisto(lepStr+idStr+"pt",tags_full, lGenP4[ilep].pt(), weight);
		  mon.fillHisto(lepStr+idStr+"eta",tags_full,lGenP4[ilep].eta(), weight);
		}
	    }
	  mon.fillHisto(lepStr+"d0",              tags_full,fabs(phys.leptons[ilep].d0),weight);
	  mon.fillHisto(lepStr+"dZ",              tags_full,fabs(phys.leptons[ilep].dZ),weight);
	  mon.fillHisto(lepStr+"trkchi2",         tags_full,fabs(phys.leptons[ilep].trkchi2),weight);
	  mon.fillHisto(lepStr+"trkvalidpixel",   tags_full,fabs(phys.leptons[ilep].trkValidPixelHits),weight);
	  mon.fillHisto(lepStr+"trkvalidtracker", tags_full,fabs(phys.leptons[ilep].trkLostInnerHits),weight);
	  mon.fillHisto(lepStr+"losthits",        tags_full,fabs(phys.leptons[ilep].trkLostInnerHits),weight);

	  TString reg="ee";
	  if(fabs(phys.leptons[ilep].eta())<1.442) reg="eb";
	  mon.fillHisto(lepStr+reg+"detain",  tags_full,fabs(ev.en_detain[ilep]),weight);
	  mon.fillHisto(lepStr+reg+"dphiin",  tags_full,fabs(ev.en_dphiin[ilep]),weight);
	  mon.fillHisto(lepStr+reg+"sihih",   tags_full,fabs(ev.en_sihih[ilep]),weight);
	  mon.fillHisto(lepStr+reg+"sipip",   tags_full,fabs(ev.en_sipip[ilep]),weight);
	  mon.fillHisto(lepStr+reg+"r9",      tags_full,fabs(ev.en_r9[ilep]),weight);
	  mon.fillHisto(lepStr+reg+"hoe",     tags_full,fabs(ev.en_hoe[ilep]),weight);
	  mon.fillHisto(lepStr+reg+"eopin",   tags_full,fabs(ev.en_eopin[ilep]),weight);
	  mon.fillHisto(lepStr+reg+"fbrem",   tags_full,fabs(ev.en_fbrem[ilep]),weight);
	  mon.fillHisto(lepStr+reg+"ooemoop", tags_full,fabs(ev.en_ooemoop[ilep]),weight);

	  mon.fillHisto(lepStr+"nmuonhits",  tags_full,fabs(ev.mn_validMuonHits[ilep]),weight);
	  mon.fillHisto(lepStr+"nmatches",   tags_full,fabs(ev.mn_nMatches[ilep]),weight);


	  // mon.fillHisto(lepStr+"pt",  tags_full, phys.leptons[ilep].pt(),weight);
	  // mon.fillHisto(lepStr+"eta", tags_full, phys.leptons[ilep].eta(),weight);
	}
		
      //mon.fillHisto("zeta"     ,tags_full,zeta   ,iweight);
      //	mon.fillHisto("zpt"      ,tags_full,zpt     ,iweight);
      // mon.fillHisto("zmass",      tags_full, zmass, iweight);  
      mon.fillHisto("nvtx"        ,tags_full, ev.nvtx,  weight);
      mon.fillHisto("rho"         ,tags_full, ev.rho,   weight);
      mon.fillHisto("rho25"       ,tags_full, ev.rho25, weight);
      mon.fill2DHisto("rhovspu"   ,tags_full, ev.ngenITpu, ev.rho,   weight);
      mon.fill2DHisto("rho25vspu" ,tags_full, ev.ngenITpu, ev.rho25, weight);
    


      /*
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





