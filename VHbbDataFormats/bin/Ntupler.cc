#include <TH1F.h>
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "VHbbAnalysis/VHbbDataFormats/interface/HbbCandidateFinderAlgo.h" 
#include "VHbbAnalysis/VHbbDataFormats/src/HbbCandidateFinderAlgo.cc"

#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEventAuxInfo.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/TriggerReader.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/TopMassReco.h"

#include <sstream>
#include <string>

#define MAXJ 30
#define MAXL 10


TTree * tscaleHLTmu = 0;
TTree * tscaleIDmu = 0;

float ScaleCSV(float CSV)
{
if(CSV < 0.68) return 1.0;
if(CSV < 0.9) return  0.96;
else  return  0.94;
}


float ScaleIsoHLT(float pt1, float eta1)
{
float ptMin,ptMax,etaMin,etaMax,scale,error;
float s1 = 0;
int count = 0;
tscaleHLTmu->SetBranchAddress("ptMin",&ptMin);
tscaleHLTmu->SetBranchAddress("ptMax",&ptMax);
tscaleHLTmu->SetBranchAddress("etaMin",&etaMin);
tscaleHLTmu->SetBranchAddress("etaMax",&etaMax);
tscaleHLTmu->SetBranchAddress("scale",&scale);
tscaleHLTmu->SetBranchAddress("error",&error);

for(int jentry = 0; jentry < tscaleHLTmu->GetEntries(); jentry++)
  {
   tscaleHLTmu->GetEntry(jentry);
   if((pt1 > ptMin) && (pt1 < ptMax) && (eta1 > etaMin) && (eta1 < etaMax))
    {
    s1 = scale;
    count++;   
    }
  }

if(count == 0 || s1 == 0) 
{
//caleFile->Close();
 return 1;
}


//aleFile->Close();
return (s1);
}



float ScaleID(float pt1, float eta1)
{

float ptMin,ptMax,etaMin,etaMax,scale,error;
float s1 = 0;
int count = 0;
tscaleIDmu->SetBranchAddress("ptMin",&ptMin);
tscaleIDmu->SetBranchAddress("ptMax",&ptMax);
tscaleIDmu->SetBranchAddress("etaMin",&etaMin);
tscaleIDmu->SetBranchAddress("etaMax",&etaMax);
tscaleIDmu->SetBranchAddress("scale",&scale);
tscaleIDmu->SetBranchAddress("error",&error);

for(int jentry = 0; jentry < tscaleIDmu->GetEntries(); jentry++)
  {

   tscaleIDmu->GetEntry(jentry);
   if((pt1 > ptMin) && (pt1 < ptMax) && (eta1 > etaMin) && (eta1 < etaMax))
    {
    s1 = scale;
    count++;   
    }
   }

if(count == 0 || s1 == 0) 
{
//caleFile->Close();
 return 1;
}

//aleFile->Close();
return (s1);

}

bool jsonContainsEvent (const std::vector< edm::LuminosityBlockRange > &jsonVec,
                        const edm::EventBase &event)
{
   // if the jsonVec is empty, then no JSON file was provided so all
   // events should pass
   if (jsonVec.empty())
   {
      return true;
   }
   bool (* funcPtr) (edm::LuminosityBlockRange const &,
                     edm::LuminosityBlockID const &) = &edm::contains;
   edm::LuminosityBlockID lumiID (event.id().run(), 
                                  event.id().luminosityBlock());
   std::vector< edm::LuminosityBlockRange >::const_iterator iter = 
      std::find_if (jsonVec.begin(), jsonVec.end(),
                    boost::bind(funcPtr, _1, lumiID) );
   return jsonVec.end() != iter;

}



  typedef struct 
  {
    float mass;  //MT in case of W
    float pt;
    float eta;
    float phi;
  } TrackInfo;
  

struct  LeptonInfo
  {
    void reset()
    {
    for(int i =0; i < MAXL;i++){      mass[i]=-99; pt[i]=-99; eta[i]=-99; phi[i]=-99; aodCombRelIso[i]=-99; pfCombRelIso[i]=-99; photonIso[i]=-99; neutralHadIso[i]=-99; chargedHadIso[i]=-99; particleIso[i]=-99; dxy[i]=-99; dz[i]=-99; type[i]=-99; }
    }

    template <class Input> void set(const Input & i, int j,int t)
    {
      type[j]=t;
      pt[j]=i.p4.Pt(); 
      mass[j]=i.p4.M();
      eta[j]=i.p4.Eta();
      phi[j]=i.p4.Phi();
      aodCombRelIso[j]=(i.hIso+i.eIso+i.tIso)/i.p4.Pt();
      pfCombRelIso[j]=(i.pfChaIso+i.pfPhoIso+i.pfNeuIso)/i.p4.Pt();
      photonIso[j]=i.pfPhoIso;
      neutralHadIso[j]=i.pfNeuIso;
      chargedHadIso[j]=i.pfChaIso;
      setID(i,j);
//FIXME: whats this?      particleIso;
    }
     template <class Input> void setID(const Input & i, int j)
     {
      id[j]=-99;
     }      
 
      
     

    float mass[MAXL];  //MT in case of W
    float pt[MAXL];
    float eta[MAXL];
    float phi[MAXL];
    float aodCombRelIso[MAXL];
    float pfCombRelIso[MAXL];
    float photonIso[MAXL];
    float neutralHadIso[MAXL];
    float chargedHadIso[MAXL];
    float particleIso[MAXL];
    float dxy[MAXL];
    float dz[MAXL];
    int type[MAXL];
    float id[MAXL];
  };
  
 template <> void LeptonInfo::setID<VHbbEvent::ElectronInfo>(const VHbbEvent::ElectronInfo & i, int j){
     id[j]=i.id80r;
  }

  typedef struct 
  {
    float et; 
    float sumet;   
    float sig;
    float phi;
  } METInfo;
  
  typedef struct 
  {
    float mht;
    float ht;  
    float sig;
    float phi;
  } MHTInfo;
  
  typedef struct 
  {
    float mass;
    float pt;
    float wMass;
  } TopInfo;

  typedef struct 
  {
    int run;
    int lumi;
    int event;
    int json;
  } EventInfo;
  
  typedef struct 
  {
    void set(const VHbbEvent::SimpleJet & j, int i) 
    {
     pt[i]=j.p4.Pt();
     eta[i]=j.p4.Eta();
     phi[i]=j.p4.Phi();
     csv[i]=j.csv;
     //FIXME:     numTracksSV  (NEED EDM FIX)
     //FIXME:    chf;    float nhf;    float cef;    float nef;    float nch; nconstituents;  (NEED EDM FIX)

     flavour[i]=j.flavour;
     //FIXME: genPt parton or genjet?  (NEED EDM FIX)

     if(j.bestMCp4.Pt() > 0)
     {
      genPt[i]=j.bestMCp4.Pt();
      genEta[i]=j.bestMCp4.Eta();
      genPhi[i]=j.bestMCp4.Phi();
     }
     //FIXME JECUnc  (NEED EDM FIX)

    }
   void reset()
   {
   for(int i=0;i<MAXJ;i++) {
    pt[i]=-99; eta[i]=-99; phi[i]=-99; csv[i]=-99; cosTheta[i]=-99; numTracksSV[i]=-99; chf[i]=-99; nhf[i]=-99; cef[i]=-99; nef[i]=-99; nch[i]=-99; nconstituents[i]=-99; flavour[i]=-99; genPt[i]=-99; genEta[i]=-99; genPhi[i]=-99; JECUnc[i]=-99;
   }
   }
    float pt[MAXJ];
    float eta[MAXJ];
    float phi[MAXJ];
    float csv[MAXJ];
    float cosTheta[MAXJ];
    int numTracksSV[MAXJ];
    float chf[MAXJ];
    float nhf[MAXJ];
    float cef[MAXJ];
    float nef[MAXJ];
    float nch[MAXJ];
    float nconstituents[MAXJ];
    float flavour[MAXJ];
    float genPt[MAXJ];
    float genEta[MAXJ];
    float genPhi[MAXJ];
    float JECUnc[MAXJ];

  } JetInfo;
  
int main(int argc, char* argv[]) 
{
  gROOT->Reset();

  TTree *_outTree;
  METInfo MET;
  MHTInfo MHT;
  TopInfo top;
  EventInfo EVENT;
//  JetInfo jet1,jet2, addJet1, addJet2;
  // lepton1,lepton2;
  JetInfo hJets, aJets;
  LeptonInfo vLeptons, aLeptons;
  int naJets=0, nhJets=0;
  TrackInfo H;
  TrackInfo V;
  int nvlep=0,nalep=0; 
  
 float jjdr,jjdPhi,jjdEta,HVdPhi,VMt,deltaPullAngle,deltaPullAngleAK7,gendrcc,gendrbb, genZpt, genWpt, weightTrig,addJet3Pt, minDeltaPhijetMET,  jetPt_minDeltaPhijetMET , PUweight;
   int nofLeptons15,nofLeptons20, Vtype,numJets,numBJets,eventFlav;
//   bool isMET80_CJ80, ispfMHT150, isMET80_2CJ20,isMET65_2CJ20, isJETID,isIsoMu17;
   bool triggerFlags[500];
  // ----------------------------------------------------------------------
  // First Part: 
  //
  //  * enable the AutoLibraryLoader 
  //  * book the histograms of interest 
  //  * open the input file
  // ----------------------------------------------------------------------

  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  gSystem->Load("libDataFormatsFWLite");
  AutoLibraryLoader::enable();
  
  // parse arguments
  if ( argc < 2 ) {
    return 0;
  }

  std::vector<VHbbCandidate> * candZlocal = new std::vector<VHbbCandidate>;
  std::vector<VHbbCandidate> * candWlocal = new std::vector<VHbbCandidate>;

  // get the python configuration
  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in  = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput" );
  const edm::ParameterSet& out = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteOutput");
  const edm::ParameterSet& ana = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("Analyzer");
  std::vector<edm::LuminosityBlockRange> jsonVector;
  if ( in.exists("lumisToProcess") ) 
  {
       std::vector<edm::LuminosityBlockRange> const & lumisTemp =
         in.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> > ("lumisToProcess");
      jsonVector.resize( lumisTemp.size() );
      copy( lumisTemp.begin(), lumisTemp.end(), jsonVector.begin() );
  }
  
  // now get each parameter
  int maxEvents_( in.getParameter<int>("maxEvents") );
  unsigned int outputEvery_( in.getParameter<unsigned int>("outputEvery") );
  std::string outputFile_( out.getParameter<std::string>("fileName" ) );
  
  
  std::vector<std::string> triggers( ana.getParameter<std::vector<std::string> >("triggers") );
  double btagThr =  ana.getParameter<double>("bJetCountThreshold" );
  bool fromCandidate = ana.getParameter<bool>("readFromCandidates");
  HbbCandidateFinderAlgo * algoZ = new HbbCandidateFinderAlgo(ana.getParameter<bool>("verbose"), ana.getParameter<double>("jetPtThresholdZ"),
                                     ana.getParameter<bool>("useHighestPtHiggsZ")                         );
  HbbCandidateFinderAlgo * algoW = new HbbCandidateFinderAlgo(ana.getParameter<bool>("verbose"), ana.getParameter<double>("jetPtThresholdW"),
                                     ana.getParameter<bool>("useHighestPtHiggsW")                         );

  TFile *hltMuFile = new TFile (ana.getParameter<std::string> ("hltMuFileName").c_str(),"read");
  tscaleHLTmu = (TTree*) hltMuFile->Get("tree");
//  hltMuFile->Close();
  TFile *idMuFile = new TFile (ana.getParameter<std::string> ("idMuFileName").c_str(),"read"); //"ScaleEffs42.root","read");
  tscaleIDmu = (TTree*) idMuFile->Get("tree");
//  idMuFile->Close();

  std::vector<std::string> inputFiles_( in.getParameter<std::vector<std::string> >("fileNames") );
//  std::string inputFile( in.getParameter<std::string> ("fileName") );

  
  std::string PUmcfileName_ = in.getParameter<std::string> ("PUmcfileName") ;
  std::string PUdatafileName_ = in.getParameter<std::string> ("PUdatafileName") ;
  bool isMC_( ana.getParameter<bool>("isMC") );  
  TriggerReader trigger(isMC_);
 
//   TFile *_outPUFile	= new TFile((outputFile_+"_PU").c_str(), "recreate");	
//   TH1F * pu = new TH1F("pileup","",-0.5,24.5,25);
   TFile *_outFile	= new TFile(outputFile_.c_str(), "recreate");	
  _outTree = new TTree("tree", "myTree");
  
  _outTree->Branch("H"		,  &H	            ,  "mass/F:pt/F:eta:phi/F");
  _outTree->Branch("V"		,  &V	            ,  "mass/F:pt/F:eta:phi/F");
  _outTree->Branch("nhJets"		,  &nhJets	            ,  "nhJets/I");
  _outTree->Branch("naJets"		,  &naJets	            ,  "naJets/I");

  _outTree->Branch("hJet_pt",hJets.pt ,"pt[nhJets]/F");
  _outTree->Branch("hJet_eta",hJets.eta ,"eta[nhJets]/F");
  _outTree->Branch("hJet_phi",hJets.phi ,"phi[nhJets]/F");
  _outTree->Branch("hJet_csv",hJets.csv ,"csv[nhJets]/F");
  _outTree->Branch("hJet_cosTheta",hJets.cosTheta ,"cosTheta[nhJets]/F");
  _outTree->Branch("hJet_numTracksSV",hJets.numTracksSV ,"numTracksSV[nhJets]/I");
  _outTree->Branch("hJet_chf",hJets.chf ,"chf[nhJets]/F");
  _outTree->Branch("hJet_nhf",hJets.nhf ,"nhf[nhJets]/F");
  _outTree->Branch("hJet_cef",hJets.cef ,"cef[nhJets]/F");
  _outTree->Branch("hJet_nef",hJets.nef ,"nef[nhJets]/F");
  _outTree->Branch("hJet_nch",hJets.nch ,"nch[nhJets]/F");
  _outTree->Branch("hJet_nconstituents",hJets.nconstituents ,"nconstituents[nhJets]");
  _outTree->Branch("hJet_flavour",hJets.flavour ,"flavour[nhJets]/F");
  _outTree->Branch("hJet_genPt",hJets.genPt ,"genPt[nhJets]/F");
  _outTree->Branch("hJet_genEta",hJets.genEta ,"genEta[nhJets]/F");
  _outTree->Branch("hJet_genPhi",hJets.genPhi ,"genPhi[nhJets]/F");
  _outTree->Branch("hJet_JECUnc",hJets.JECUnc ,"JECUnc[nhJets]/F");

  _outTree->Branch("aJet_pt",aJets.pt ,"pt[naJets]/F");
  _outTree->Branch("aJet_eta",aJets.eta ,"eta[naJets]/F");
  _outTree->Branch("aJet_phi",aJets.phi ,"phi[naJets]/F");
  _outTree->Branch("aJet_csv",aJets.csv ,"csv[naJets]/F");
  _outTree->Branch("aJet_cosTheta",aJets.cosTheta ,"cosTheta[naJets]/F");
  _outTree->Branch("aJet_numTracksSV",aJets.numTracksSV ,"numTracksSV[naJets]/I");
  _outTree->Branch("aJet_chf",aJets.chf ,"chf[naJets]/F");
  _outTree->Branch("aJet_nhf",aJets.nhf ,"nhf[naJets]/F");
  _outTree->Branch("aJet_cef",aJets.cef ,"cef[naJets]/F");
  _outTree->Branch("aJet_nef",aJets.nef ,"nef[naJets]/F");
  _outTree->Branch("aJet_nch",aJets.nch ,"nch[naJets]/F");
  _outTree->Branch("aJet_nconstituents",aJets.nconstituents ,"nconstituents[naJets]");
  _outTree->Branch("aJet_flavour",aJets.flavour ,"flavour[naJets]/F");
  _outTree->Branch("aJet_genPt",aJets.genPt ,"genPt[naJets]/F");
  _outTree->Branch("aJet_genEta",aJets.genEta ,"genEta[naJets]/F");
  _outTree->Branch("aJet_genPhi",aJets.genPhi ,"genPhi[naJets]/F");
  _outTree->Branch("aJet_JECUnc",aJets.JECUnc ,"JECUnc[naJets]/F");


  _outTree->Branch("addJet3Pt", &addJet3Pt  ,  "addJet3Pt/F");
  _outTree->Branch("jjdr" 	,  &jjdr            ,  "jjdr/F"         );         	
  _outTree->Branch("jjdPhi"  	,  &jjdPhi          ,  "jjdPhi/F"       );            	
  _outTree->Branch("jjdEta"  	,  &jjdEta          ,  "jjdEta/F"       );            	
  _outTree->Branch("numJets"      ,  &numJets         ,  "numJets/I"       );                
  _outTree->Branch("numBJets"      ,  &numBJets         ,  "numBJets/I"       );                
  _outTree->Branch("nofLeptons15"   ,  &nofLeptons15      ,  "nofLeptons15/I"    );                
  _outTree->Branch("nofLeptons20"   ,  &nofLeptons20      ,  "nofLeptons20/I"    );                
  _outTree->Branch("deltaPullAngle", &deltaPullAngle  ,  "deltaPullAngle/F");
  _outTree->Branch("gendrcc"    , &gendrcc      ,  "gendrcc/F");
  _outTree->Branch("gendrbb"    , &gendrbb      ,  "gendrbb/F");
  _outTree->Branch("genZpt"    , &genZpt      ,  "genZpt/F");
  _outTree->Branch("genWpt"    , &genWpt      ,  "genWpt/F");
  _outTree->Branch("weightTrig"        , &weightTrig          ,  "weightTrig/F");
  _outTree->Branch("deltaPullAngleAK7", &deltaPullAngleAK7  ,  "deltaPullAngleAK7/F");
  _outTree->Branch("PUweight",       &PUweight  ,  "PUweight/F");
  _outTree->Branch("eventFlav",       &eventFlav  ,  "eventFlav/I");
 

   
    
  _outTree->Branch("Vtype"     ,  &Vtype   ,   "Vtype/I" );                
  _outTree->Branch("HVdPhi"     ,  &HVdPhi   ,   "HVdPhi/F" );                
  _outTree->Branch("VMt"  	,  &VMt      ,   "VMt/F"    );             	

  _outTree->Branch("nvlep"	,  &nvlep    ,   "nvlep/I");
  _outTree->Branch("nalep"	,  &nalep    ,   "nalep/I");

  _outTree->Branch("vLepton_mass",vLeptons.mass ,"mass[nvlep]/F");
  _outTree->Branch("vLepton_pt",vLeptons.pt ,"pt[nvlep]/F");
  _outTree->Branch("vLepton_eta",vLeptons.eta ,"eta[nvlep]");
  _outTree->Branch("vLepton_phi",vLeptons.phi ,"phi[nvlep]/F");
  _outTree->Branch("vLepton_aodCombRelIso",vLeptons.aodCombRelIso ,"aodCombRelIso[nvlep]/F");
  _outTree->Branch("vLepton_pfCombRelIso",vLeptons.pfCombRelIso ,"pfCombRelIso[nvlep]/F");
  _outTree->Branch("vLepton_photonIso",vLeptons.photonIso ,"photonIso[nvlep]/F");
  _outTree->Branch("vLepton_neutralHadIso",vLeptons.neutralHadIso ,"neutralHadIso[nvlep]/F");
  _outTree->Branch("vLepton_chargedHadIso",vLeptons.chargedHadIso ,"chargedHadIso[nvlep]/F");
  _outTree->Branch("vLepton_particleIso",vLeptons.particleIso ,"particleIso[nvlep]/F");
  _outTree->Branch("vLepton_dxy",vLeptons.dxy ,"dxy[nvlep]/F");
  _outTree->Branch("vLepton_dz",vLeptons.dz ,"dz[nvlep]/F");
  _outTree->Branch("vLepton_type",vLeptons.type ,"type[nvlep]/I");
  _outTree->Branch("vLepton_id",vLeptons.id ,"id[nvlep]/F");
 
  _outTree->Branch("aLepton_mass",aLeptons.mass ,"mass[nalep]/F");
  _outTree->Branch("aLepton_pt",aLeptons.pt ,"pt[nalep]/F");
  _outTree->Branch("aLepton_eta",aLeptons.eta ,"eta[nalep]");
  _outTree->Branch("aLepton_phi",aLeptons.phi ,"phi[nalep]/F");
  _outTree->Branch("aLepton_aodCombRelIso",aLeptons.aodCombRelIso ,"aodCombRelIso[nalep]/F");
  _outTree->Branch("aLepton_pfCombRelIso",aLeptons.pfCombRelIso ,"pfCombRelIso[nalep]/F");
  _outTree->Branch("aLepton_photonIso",aLeptons.photonIso ,"photonIso[nalep]/F");
  _outTree->Branch("aLepton_neutralHadIso",aLeptons.neutralHadIso ,"neutralHadIso[nalep]/F");
  _outTree->Branch("aLepton_chargedHadIso",aLeptons.chargedHadIso ,"chargedHadIso[nalep]/F");
  _outTree->Branch("aLepton_particleIso",aLeptons.particleIso ,"particleIso[nalep]/F");
  _outTree->Branch("aLepton_dxy",aLeptons.dxy ,"dxy[nalep]/F");
  _outTree->Branch("aLepton_dz",aLeptons.dz ,"dz[nalep]/F");
  _outTree->Branch("aLepton_type",aLeptons.type ,"type[nalep]/I");
  _outTree->Branch("aLepton_id",aLeptons.id ,"id[nalep]/F");
 
  _outTree->Branch("top"		,  &top	         ,   "mass/F:pt/F:wMass/F");

  _outTree->Branch("MET"		,  &MET	         ,   "et/F:sumet:sig/F:phi/F");
  _outTree->Branch("MHT"		,  &MHT	         ,   "mht/F:ht:sig/F:phi/F");
  _outTree->Branch("minDeltaPhijetMET"		,  &minDeltaPhijetMET	         ,   "minDeltaPhijetMET/F");
  _outTree->Branch("jetPt_minDeltaPhijetMET"		,  &jetPt_minDeltaPhijetMET	         ,   "jetPt_minDeltaPhijetMET/F");

   std::stringstream s;
   s << "triggerFlags[" << triggers.size() << "]/b";
   _outTree->Branch("triggerFlags", triggerFlags, s.str().c_str()); 
  
   _outTree->Branch("EVENT"		,  &EVENT	         ,   "run/I:lumi/I:event/I:json/I");
 
   /*
      FIXME - btag SF
      FIXME - HBHE filter (NEED EDM)
    */

    int ievt=0;  
    int totalcount=0;

//  TFile* inFile = new TFile(inputFile.c_str(), "read");
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile) {
    std::cout << iFile << std::endl;
    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
    if(inFile==0) continue;

  // loop the events
      
      fwlite::Event ev(inFile);
      for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt)
        {

 	  if(isMC_){
 	  // PU weights
          
// 	  edm::LumiReWeighting   LumiWeights_ = edm::LumiReWeighting(PUmcfileName_,PUdatafileName_ , "pileup", "pileup");
// 	  double avg=0;
//FIXME:  PU (NEED EDM FIX)
// 	   if( PUintimeSizes.isValid() && PUouttime1minusSizes.isValid() && PUouttime1plusSizes.isValid()){
// 	     avg = (double)( *PUintimeSizes );
// 	   }
 	   PUweight = 1.0; // FIXME: LumiWeights_.weight3BX( avg /3.);  (NEED EDM FIX)
 	  }
      
      //Write event info 
      EVENT.run = ev.id().run();
      EVENT.lumi = ev.id().luminosityBlock();
      EVENT.event = ev.id().event();
      EVENT.json = jsonContainsEvent (jsonVector, ev);

 
 const std::vector<VHbbCandidate> * candZ ;
 const std::vector<VHbbCandidate> * candW ;
 const VHbbEvent *  iEvent =0;
 if(fromCandidate)
  {
	fwlite::Handle< std::vector<VHbbCandidate> > vhbbCandHandleZ;
    	vhbbCandHandleZ.getByLabel(ev,"hbbBestCSVPt20Candidates");
    	candZ = vhbbCandHandleZ.product();

   	fwlite::Handle< std::vector<VHbbCandidate> > vhbbCandHandle;
    	vhbbCandHandle.getByLabel(ev,"hbbHighestPtHiggsPt30Candidates");
    	candW = vhbbCandHandle.product();
   }
 else
   {
      candZlocal->clear();
      candWlocal->clear();
      fwlite::Handle< VHbbEvent > vhbbHandle; 
      vhbbHandle.getByLabel(ev,"HbbAnalyzerNew");
      iEvent = vhbbHandle.product();
      algoZ->run(vhbbHandle.product(),*candZlocal);
      algoW->run(vhbbHandle.product(),*candWlocal);
      candZ= candZlocal; 
      candW= candWlocal; 


 
    }

        const std::vector<VHbbCandidate> * cand = candZ;


      fwlite::Handle< VHbbEventAuxInfo > vhbbAuxHandle; 
      vhbbAuxHandle.getByLabel(ev,"HbbAnalyzerNew");
      const VHbbEventAuxInfo & aux = *vhbbAuxHandle.product();
      
      /*  fwlite::Handle< VHbbEvent > vhbbHandle; 
          vhbbHandle.getByLabel(ev,"HbbAnalyzerNew");
          const VHbbEvent iEvent = *vhbbHandle.product();
      */

      //      std::clog << "Filling tree "<< std::endl;
      
     if(cand->size() == 0 or cand->at(0).H.jets.size() < 2) continue;
          if(cand->size() > 1 ) 
          {
           std::cout << "MULTIPLE CANDIDATES: " << cand->size() << std::endl;
          }
          if(cand->at(0).candidateType == VHbbCandidate::Wmun || cand->at(0).candidateType == VHbbCandidate::Wen ) cand=candW;
          if(cand->size() == 0) 
          {
//            std::cout << "W event loss due to tigther cuts" << std::endl;
            continue;
          }
          const VHbbCandidate & vhCand =  cand->at(0);
          trigger.setEvent(&ev);
          for(size_t j=0;j < triggers.size();j++)
          triggerFlags[j]=trigger.accept(triggers[j]);
     
          eventFlav=0;
          if(aux.mcBbar.size() > 0 || aux.mcB.size() > 0) eventFlav=5;
          else if(aux.mcC.size() > 0) eventFlav=4;
         

          H.mass = vhCand.H.p4.M();
          H.pt = vhCand.H.p4.Pt();
          H.eta = vhCand.H.p4.Eta();
          H.phi = vhCand.H.p4.Phi();
          V.mass = vhCand.V.p4.M();
          V.pt = vhCand.V.p4.Pt();
          V.eta = vhCand.V.p4.Eta();
          V.phi = vhCand.V.p4.Phi();
          nhJets=2;
          hJets.set(vhCand.H.jets[0],0);
          hJets.set(vhCand.H.jets[1],1);
          aJets.reset();
          naJets=vhCand.additionalJets.size();
          numBJets=0;
          if(vhCand.H.jets[0].csv> btagThr) numBJets++;
          if(vhCand.H.jets[1].csv> btagThr) numBJets++;
          for( int j=0; j < naJets && j < MAXJ; j++ ) 
          {
              aJets.set(vhCand.additionalJets[j],j);
              if(vhCand.additionalJets[j].csv> btagThr) numBJets++;
          }   
          numJets = vhCand.additionalJets.size()+2;
          jjdr = deltaR(vhCand.H.jets[0].p4.Eta(),vhCand.H.jets[0].p4.Phi(),vhCand.H.jets[1].p4.Eta(),vhCand.H.jets[1].p4.Phi());
          jjdPhi = deltaPhi(vhCand.H.jets[0].p4.Phi(),vhCand.H.jets[1].p4.Phi());
          jjdEta= TMath::Abs( vhCand.H.jets[0].p4.Eta() - vhCand.H.jets[1].p4.Eta() );
          HVdPhi = deltaPhi(vhCand.H.p4.Phi(),vhCand.V.p4.Phi()) ;
          VMt = vhCand.Mt() ;
          deltaPullAngle = vhCand.H.deltaTheta;
          float deltaPhipfMETjet1 = deltaPhi( vhCand.V.mets.at(0).p4.Phi(), vhCand.H.jets[0].p4.Phi() );
          float deltaPhipfMETjet2 = deltaPhi( vhCand.V.mets.at(0).p4.Phi(), vhCand.H.jets[1].p4.Phi() );
          if(deltaPhipfMETjet1 <= deltaPhipfMETjet2) 
          {
	     minDeltaPhijetMET=deltaPhipfMETjet1;
             jetPt_minDeltaPhijetMET=vhCand.H.jets[0].p4.Pt(); 
          }
	  else
          {
	     minDeltaPhijetMET=deltaPhipfMETjet2;
             jetPt_minDeltaPhijetMET=vhCand.H.jets[1].p4.Pt(); 
          }
  
          hJets.cosTheta[0]=  vhCand.H.helicities[0];
          hJets.cosTheta[1]=  vhCand.H.helicities[1];

          MET.et = vhCand.V.mets.at(0).p4.Pt();
          MET.phi = vhCand.V.mets.at(0).p4.Phi();
          MET.sumet = vhCand.V.mets.at(0).sumEt;
          MET.sig = vhCand.V.mets.at(0).metSig;
//FIXME  add MHT     _outTree->Branch("MHT"            ,  &MHT          ,   "mht/F:ht:sig/F:phi/F");  (NEED EDM FIX)
          Vtype = vhCand.candidateType;
          vLeptons.reset();
          weightTrig = 0.; 
	  TLorentzVector leptonForTop;
          size_t firstAddMu=0;
          size_t firstAddEle=0;
          if(Vtype == VHbbCandidate::Zmumu ){
                  vLeptons.set(vhCand.V.muons[0],0,13); 
                  vLeptons.set(vhCand.V.muons[1],1,13);
                  float cweightID = ScaleID(vLeptons.pt[0],vLeptons.eta[0]) * ScaleID(vLeptons.pt[1],vLeptons.eta[1]) ;
                  float weightTrig1 = ScaleIsoHLT(vLeptons.pt[0],vLeptons.eta[0]);
                  float weightTrig2 = ScaleIsoHLT(vLeptons.pt[1],vLeptons.eta[1]);
                  float cweightTrig = weightTrig1 + weightTrig2 - weightTrig1*weightTrig2;
                  weightTrig = cweightID * cweightTrig;
                  nvlep=2;
                  firstAddMu=2;
          }
          if( Vtype == VHbbCandidate::Zee ){
	          vLeptons.set(vhCand.V.electrons[0],0,11);
            	  vLeptons.set(vhCand.V.electrons[1],1,11);
                  nvlep=2;
		  firstAddEle=2;
                  //FIXME: trigger weights for electrons
           }
          if(Vtype == VHbbCandidate::Wmun ){
                  leptonForTop=vhCand.V.muons[0].p4;
                  vLeptons.set(vhCand.V.muons[0],0,13); 
                  float cweightID = ScaleID(vLeptons.pt[0],vLeptons.eta[0]);
                  float weightTrig1 = ScaleIsoHLT(vLeptons.pt[0],vLeptons.eta[0]);
                  float cweightTrig = weightTrig1;
                  weightTrig = cweightID * cweightTrig;
                  nvlep=1;
		  firstAddMu=1;
          }
          if( Vtype == VHbbCandidate::Wen ){
                  leptonForTop=vhCand.V.electrons[0].p4;
            	  vLeptons.set(vhCand.V.electrons[0],0,11);
		  nvlep=1;
		  firstAddEle=1;
           }
          if( Vtype == VHbbCandidate::Znn ){
                  nvlep=0;
                  //FIXME: trigger weights for Znn

          }
          
          aLeptons.reset();
          nalep=0;
          for(size_t j=firstAddMu;j< vhCand.V.muons.size();j++) aLeptons.set(vhCand.V.muons[j],nalep++,13);
          for(size_t j=firstAddEle;j< vhCand.V.electrons.size();j++) aLeptons.set(vhCand.V.electrons[j],nalep++,11);


          double maxBtag=-99999;
          TLorentzVector bJet;
          for(unsigned int j=0; j < vhCand.H.jets.size(); j++ ){
                if (vhCand.H.jets[j].csv > maxBtag) { bJet=vhCand.H.jets[j].p4 ; maxBtag =vhCand.H.jets[j].csv; }
          }
          for(unsigned int j=0; j < vhCand.additionalJets.size(); j++ ){
                if (vhCand.additionalJets[j].csv > maxBtag) { bJet=vhCand.additionalJets[j].p4 ; maxBtag =vhCand.additionalJets[j].csv; }
          }
          if(maxBtag > -99999)
          { 
           TopHypo topQuark = TopMassReco::topMass(leptonForTop,bJet,vhCand.V.mets.at(0).p4);
           top.mass = topQuark.p4.M();
           top.pt = topQuark.p4.Pt();
           top.wMass = topQuark.p4W.M();
          } else {
              top.mass = -99;
              top.pt = -99;
              top.wMass = -99;
          }
  

       
//FIXME  _outTree->Branch("nofLeptons15"   ,  &nofLeptons15      ,  "nofLeptons15/I"    );
          nofLeptons20= vhCand.additionalLeptons();
// if(aux.mcC.size() >=2)
// std::cout << "C Must not be zero and it is ... " << aux.mcC[1].p4.Pt() << std::endl;
// if(aux.mcB.size() >=1)
// std::cout << "B Must not be zero and it is ... " << aux.mcB[0].p4.Pt() << std::endl;

// FIXME         gendrcc=aux.genCCDeltaR();  (NEED EDM FIX)

// FIXME    gendrbb=aux.genBBDeltaR();  (NEED EDM FIX)
          genZpt=aux.mcZ.size() > 0 ? aux.mcZ[0].p4.Pt():-99;
          genWpt=aux.mcW.size() > 0 ? aux.mcW[0].p4.Pt():-99;

//FIXME:  _outTree->Branch("deltaPullAngleAK7", &deltaPullAngleAK7  ,  "deltaPullAngleAK7/F");



	  _outTree->Fill();

	}// closed event loop

      std::cout << "closing the file: " << inputFiles_[iFile] << std::endl;
      inFile->Close();
	   // close input file
   } // loop on files
     
  
    std::cout << "Events: " << ievt <<std::endl;
    std::cout << "TotalCount: " << totalcount <<std::endl;

    
    
    _outFile->cd();
    
    _outTree->Write();
    _outFile->Write();
    _outFile->Close();
    return 0;
}


