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


#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEventAuxInfo.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/TriggerReader.h"

#include <sstream>
#include <string>

#define MAXJ 30
#define MAXL 10

float ScaleCSV(float CSV)
{
if(CSV < 0.68) return 1.0;
if(CSV < 0.9) return  0.96;
else  return  0.94;
}


float ScaleIsoHLT(float pt1, float eta1)
{
//FIXME: get the files for HLT
return 1;
float ptMin,ptMax,etaMin,etaMax,scale,error;
float s1 = 0;

TFile *scaleFile = new TFile ("IsoToHLT42.root","read");
TTree *tscale = (TTree*) scaleFile->Get("tree");
int count = 0;
tscale->SetBranchAddress("ptMin",&ptMin);
tscale->SetBranchAddress("ptMax",&ptMax);
tscale->SetBranchAddress("etaMin",&etaMin);
tscale->SetBranchAddress("etaMax",&etaMax);
tscale->SetBranchAddress("scale",&scale);
tscale->SetBranchAddress("error",&error);

for(int jentry = 0; jentry < tscale->GetEntries(); jentry++)
  {
   tscale->GetEntry(jentry);
   if((pt1 > ptMin) && (pt1 < ptMax) && (eta1 > etaMin) && (eta1 < etaMax))
    {
    s1 = scale;
    count++;   
    }
  }

if(count == 0 || s1 == 0) 
{
 scaleFile->Close();
 return 1;
}


scaleFile->Close();
return (s1);
}



float ScaleID(float pt1, float eta1)
{
//FIXME: get the files for ID
return 1;

float ptMin,ptMax,etaMin,etaMax,scale,error;
float s1 = 0;

TFile *scaleFile = new TFile ("ScaleEffs42.root","read");
TTree *tscale = (TTree*) scaleFile->Get("tree");
int count = 0;
tscale->SetBranchAddress("ptMin",&ptMin);
tscale->SetBranchAddress("ptMax",&ptMax);
tscale->SetBranchAddress("etaMin",&etaMin);
tscale->SetBranchAddress("etaMax",&etaMax);
tscale->SetBranchAddress("scale",&scale);
tscale->SetBranchAddress("error",&error);

for(int jentry = 0; jentry < tscale->GetEntries(); jentry++)
  {

   tscale->GetEntry(jentry);
   if((pt1 > ptMin) && (pt1 < ptMax) && (eta1 > etaMin) && (eta1 < etaMax))
    {
    s1 = scale;
    count++;   
    }
   }

if(count == 0 || s1 == 0) 
{
 scaleFile->Close();
 return 1;
}

scaleFile->Close();
return (s1);

}



  typedef struct 
  {
    float mass;  //MT in case of W
    float pt;
    float eta;
    float phi;
  } _TrackInfo;
  

struct  _LeptonInfo
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
//FIXME: whats this?      particleIso;
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
  };
  
  typedef struct 
  {
    float et; 
    float sumet;   
    float sig;
    float phi;
  } _METInfo;
  
  typedef struct 
  {
    float mht;
    float ht;  
    float sig;
    float phi;
  } _MHTInfo;

  struct 
  {
    int run;
    int lumi;
    int event;
  } _EventInfo;
  

  typedef struct 
  {
    void set(const VHbbEvent::SimpleJet & j, int i) 
    {
     pt[i]=j.p4.Pt();
     eta[i]=j.p4.Eta();
     phi[i]=j.p4.Phi();
     csv[i]=j.csv;
     //FIXME: cosTheta=tVector;
     //FIXME:     numTracksSV
     //FIXME:    chf;    float nhf;    float cef;    float nef;    float nch; nconstituents;
     flavour[i]=j.flavour;
     //FIXME: genPt parton or genjet?
/*     genPt=j.bestMCp4.Pt();
     genEta=j.bestMCp4.Eta();
     genPhi=j.bestMCp4.Phi();*/
     //FIXME JECUnc
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

  } _JetInfo;
  
int main(int argc, char* argv[]) 
{
  gROOT->Reset();

  TTree *_outTree;
  _METInfo MET;
  _MHTInfo MHT;
//  _JetInfo jet1,jet2, addJet1, addJet2;
  _JetInfo hJets, aJets;
  int naJets=0, nhJets=0;
  _TrackInfo H;
  _TrackInfo V;
  _LeptonInfo leptons; // lepton1,lepton2;
  int nlep=0; 
  
 float jjdr,jjdPhi,HVdPhi,VMt,deltaPullAngle,deltaPullAngleAK7,gendrcc,gendrbb, genZpt, genWpt, weightTrig,addJet3Pt, minDeltaPhijetMET,  jetPt_minDeltaPhijetMET , PUweight;
   int nofLeptons15,nofLeptons20, Vtype,numJets,numBJets;
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

  // get the python configuration
  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in  = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput" );
  const edm::ParameterSet& out = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteOutput");
  const edm::ParameterSet& ana = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("Analyzer");
  
  // now get each parameter
  int maxEvents_( in.getParameter<int>("maxEvents") );
  unsigned int outputEvery_( in.getParameter<unsigned int>("outputEvery") );
  std::string outputFile_( out.getParameter<std::string>("fileName" ) );
  
  
  std::vector<std::string> triggers( ana.getParameter<std::vector<std::string> >("triggers") );

  std::vector<std::string> inputFiles_( in.getParameter<std::vector<std::string> >("fileNames") );
//  std::string inputFile( in.getParameter<std::string> ("fileName") );

  
  std::string PUmcfileName_;//FIXME ( in.getParameter<std::string> ("PUmcfileName") );
  std::string PUdatafileName_;//FIXME ( in.getParameter<std::string> ("PUdatafileName") );
  bool isMC_( ana.getParameter<bool>("isMC") );  
  TriggerReader trigger(isMC_);
 
   TFile *_outFile	= new TFile(outputFile_.c_str(), "recreate");	
  _outTree = new TTree("tree", "myTree");
  
  _outTree->Branch("H"		,  &H	            ,  "mass/F:pt/F:eta:phi/F");
  _outTree->Branch("V"		,  &V	            ,  "mass/F:pt/F:eta:phi/F");
  _outTree->Branch("nhJets"		,  &nhJets	            ,  "nhJets/I");
  _outTree->Branch("naJets"		,  &naJets	            ,  "naJets/I");
//  _outTree->Branch("hJet_"		,  &hJets	    ,  "pt[nhJets]/F:eta[nhJets]/F:phi[nhJets]/F:csv[nhJets]/F:cosTheta[nhJets]/F:numTracksSV[nhJets]/I:chf[nhJets]/F:nhf[nhJets]:cef[nhJets]:nef[nhJets]:nch[nhJets]:nconstituents[nhJets]:flavour[nhJets]:genPt[nhJets]:genEta[nhJets]:genPhi[nhJets]:JECUnc[nhJets]/F");
//  _outTree->Branch("aJet_"		,  &aJets	    ,  "pt[naJets]/F:eta[naJets]/F:phi[naJets]/F:csv[naJets]/F:cosTheta[naJets]/F:numTracksSV[naJets]/I:chf[naJets]/F:nhf[naJets]:cef[naJets]:nef[naJets]:nch[naJets]:nconstituents[naJets]:flavour[naJets]:genPt[naJets]:genEta[naJets]:genPhi[naJets]:JECUnc[naJets]/F");

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
 

   
    
  _outTree->Branch("Vtype"     ,  &Vtype   ,   "Vtype/I" );                
  _outTree->Branch("HVdPhi"     ,  &HVdPhi   ,   "HVdPhi/F" );                
  _outTree->Branch("VMt"  	,  &VMt      ,   "VMt/F"    );             	

  _outTree->Branch("nlep"	,  &nlep    ,   "nlep/I");
//_outTree->Branch("leptons"	,  &leptons    ,   "mass[nlep]/F:pt[nlep]/F:eta[nlep]:phi[nlep]/F:aodCombRelIso[nlep]/F:pfCombRelIso[nlep]/F:photonIso[nlep]/F:neutralHadIso[nlep]/F:chargedHadIso[nlep]/F:particleIso[nlep]/F:dxy[nlep]/F:dz[nlep]/F:type[nlep]/I");
  _outTree->Branch("lepton_mass",leptons.mass ,"mass[nlep]/F");
  _outTree->Branch("lepton_pt",leptons.pt ,"pt[nlep]/F");
  _outTree->Branch("lepton_eta",leptons.eta ,"eta[nlep]");
  _outTree->Branch("lepton_phi",leptons.phi ,"phi[nlep]/F");
  _outTree->Branch("lepton_aodCombRelIso",leptons.aodCombRelIso ,"aodCombRelIso[nlep]/F");
  _outTree->Branch("lepton_pfCombRelIso",leptons.pfCombRelIso ,"pfCombRelIso[nlep]/F");
  _outTree->Branch("lepton_photonIso",leptons.photonIso ,"photonIso[nlep]/F");
  _outTree->Branch("lepton_neutralHadIso",leptons.neutralHadIso ,"neutralHadIso[nlep]/F");
  _outTree->Branch("lepton_chargedHadIso",leptons.chargedHadIso ,"chargedHadIso[nlep]/F");
  _outTree->Branch("lepton_particleIso",leptons.particleIso ,"particleIso[nlep]/F");
  _outTree->Branch("lepton_dxy",leptons.dxy ,"dxy[nlep]/F");
  _outTree->Branch("lepton_dz",leptons.dz ,"dz[nlep]/F");
  _outTree->Branch("lepton_type",leptons.type ,"type[nlep]/I");

//_outTree->Branch("lepton1"		,  &lepton1    ,   "mass/F:pt/F:eta:phi/F:aodCombRelIso/F:pfCombRelIso/F:photonIso/F:neutralHadIso/F:chargedHadIso/F:particleIso/F:dxy/F:dz/F:type/I");
//_outTree->Branch("lepton2"		,  &lepton2    ,   "mass/F:pt/F:eta:phi/F:aodCombRelIso/F:pfCombRelIso/F:photonIso/F:neutralHadIso/F:chargedHadIso/F:particleIso/F:dxy/F:dz/F:type/I");

//FIXME: add something about ELE id
  _outTree->Branch("MET"		,  &MET	         ,   "et/F:sumet:sig/F:phi/F");
  _outTree->Branch("MHT"		,  &MHT	         ,   "mht/F:ht:sig/F:phi/F");
  _outTree->Branch("minDeltaPhijetMET"		,  &minDeltaPhijetMET	         ,   "minDeltaPhijetMET/F");
  _outTree->Branch("jetPt_minDeltaPhijetMET"		,  &jetPt_minDeltaPhijetMET	         ,   "jetPt_minDeltaPhijetMET/F");



/*for(size_t j=0;j<triggers.size();j++)
   {
    _outTree->Branch(triggers[j].c_str(), &triggerFlags+j, (triggers[j]+"/b".c_str())); 
   }
*/
   std::stringstream s;
   s << "triggerFlags[" << triggers.size() << "]/b";
   _outTree->Branch("triggerFlags", triggerFlags, s.str().c_str()); 
 

   /* FIXME:        - trigger bits
      FIXME - top event reco
      FIXME - electrons id ?
      FIXME - add B/C/Light classification
    */

    int ievt=0;  
    int totalcount=0;

//  TFile* inFile = new TFile(inputFile.c_str(), "read");
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile) {
    std::cout << iFile << std::endl;
    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
    if(inFile==0) continue;

  // loop the events
      
      // ----------------------------------------------------------------------
      // Second Part: 
      //
      //  * loop the events in the input file 
      //  * receive the collections of interest via fwlite::Handle
      //  * fill the histograms
      //  * after the loop close the input file
      // ----------------------------------------------------------------------
      fwlite::Event ev(inFile);
      for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt)
        {

 	  if(isMC_){
 	  // PU weights
 	  edm::LumiReWeighting   LumiWeights_ = edm::LumiReWeighting(PUmcfileName_,PUdatafileName_ , "pileup", "pileup");
 	  double avg=0;
//FIXME
// 	   if( PUintimeSizes.isValid() && PUouttime1minusSizes.isValid() && PUouttime1plusSizes.isValid()){
// 	     avg = (double)( *PUintimeSizes );
// 	   }
 	   PUweight = 1.0; // FIXME: LumiWeights_.weight3BX( avg /3.);
 	  }
 

      fwlite::Handle< std::vector<VHbbCandidate> > vhbbCandHandle; 
      vhbbCandHandle.getByLabel(ev,"hbbBestCSVPt20Candidates");
      const std::vector<VHbbCandidate> & cand = *vhbbCandHandle.product();


      fwlite::Handle< VHbbEventAuxInfo > vhbbAuxHandle; 
      vhbbAuxHandle.getByLabel(ev,"HbbAnalyzerNew");
      const VHbbEventAuxInfo & aux = *vhbbAuxHandle.product();
      
      /*  fwlite::Handle< VHbbEvent > vhbbHandle; 
          vhbbHandle.getByLabel(ev,"HbbAnalyzerNew");
          const VHbbEvent iEvent = *vhbbHandle.product();
      */
      trigger.setEvent(&ev);
      for(size_t j=0;j < triggers.size();j++)
        triggerFlags[j]=trigger.accept(triggers[j]);

      //      std::clog << "Filling tree "<< std::endl;
      
     if(cand.size() == 0 or cand.at(0).H.jets.size() < 2) continue;
          const VHbbCandidate & vhCand =  cand.at(0);
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
          for( int j=0; j < naJets && j < MAXJ; j++ ) aJets.set(vhCand.additionalJets[j],j);
          numJets = vhCand.additionalJets.size()+2;
//FIXME:  _outTree->Branch("numBJets"      ,  &numBJets         ,  "numBJets/I"       );
          jjdr = deltaR(vhCand.H.jets[0].p4.Eta(),vhCand.H.jets[0].p4.Phi(),vhCand.H.jets[1].p4.Eta(),vhCand.H.jets[1].p4.Phi());
          jjdPhi = deltaPhi(vhCand.H.jets[0].p4.Phi(),vhCand.H.jets[1].p4.Phi());
          HVdPhi = deltaPhi(vhCand.H.p4.Phi(),vhCand.V.p4.Phi()) ;
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
  
//FIXME: should we add?          DeltaEtabb = TMath::Abs( vhCand.H.jets[0].p4.Eta() - vhCand.H.jets[1].p4.Eta() );
//FIXME: should we add?          helicity = vhCand.H.helicities[0];

          MET.et = vhCand.V.mets.at(0).p4.Pt();
          MET.phi = vhCand.V.mets.at(0).p4.Phi();
          MET.sumet = vhCand.V.mets.at(0).sumEt;
          MET.sig = vhCand.V.mets.at(0).metSig;
//FIXME  add MHT     _outTree->Branch("MHT"            ,  &MHT          ,   "mht/F:ht:sig/F:phi/F");
          Vtype = vhCand.candidateType;
          leptons.reset();
          
          if(Vtype == VHbbCandidate::Zmumu ){
            leptons.set(vhCand.V.muons[0],0,12); //FIXME set the type
            leptons.set(vhCand.V.muons[1],1,12);
                  float cweightID = ScaleID(leptons.pt[0],leptons.eta[0]) * ScaleID(leptons.pt[1],leptons.eta[1]) ;
                  float weightTrig1 = ScaleIsoHLT(leptons.pt[0],leptons.eta[0]);
                  float weightTrig2 = ScaleIsoHLT(leptons.pt[1],leptons.eta[1]);
                  float cweightTrig = weightTrig1 + weightTrig2 - weightTrig1*weightTrig2;
                  weightTrig = cweightID * cweightTrig;
                  nlep=2;
          }
          if( Vtype == VHbbCandidate::Zee ){
	          leptons.set(vhCand.V.electrons[0],0,11);
            	  leptons.set(vhCand.V.electrons[1],1,11);
                  nlep=2;
           }
          if(Vtype == VHbbCandidate::Wmun ){
            leptons.set(vhCand.V.muons[0],0,12); //FIXME set the type
                  float cweightID = ScaleID(leptons.pt[0],leptons.eta[0]);
                  float weightTrig1 = ScaleIsoHLT(leptons.pt[0],leptons.eta[0]);
                  float cweightTrig = weightTrig1;
                  weightTrig = cweightID * cweightTrig;
                  nlep=1;
          }
          if( Vtype == VHbbCandidate::Wen ){
            	  leptons.set(vhCand.V.electrons[0],0,11);
		  nlep=1;
           }
//FIXME  _outTree->Branch("nofLeptons15"   ,  &nofLeptons15      ,  "nofLeptons15/I"    );
          nofLeptons20= vhCand.additionalLeptons();
// if(aux.mcC.size() >=2)
// std::cout << "C Must not be zero and it is ... " << aux.mcC[1].p4.Pt() << std::endl;
// if(aux.mcB.size() >=1)
// std::cout << "B Must not be zero and it is ... " << aux.mcB[0].p4.Pt() << std::endl;

// FIXME         gendrcc=aux.genCCDeltaR();
// FIXME    gendrbb=aux.genBBDeltaR();
          genZpt=aux.mcZ.size() > 0 ? aux.mcZ[0].p4.Pt():-99;
          genWpt=aux.mcW.size() > 0 ? aux.mcW[0].p4.Pt():-99;
//FIXME:  _outTree->Branch("weightTrig"        , &weightTrig          ,  "weightTrig/F");
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


