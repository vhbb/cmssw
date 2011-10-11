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

//for IVF
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include <DataFormats/GeometrySurface/interface/Surface.h>
#include "Math/SMatrix.h"

//Move class definition to Ntupler.h ?
//#include "VHbbAnalysis/VHbbDataFormats/interface/Ntupler.h"

//for SimBhad: FIXME
//#include "VHbbAnalysis/BAnalysis/interface/SimBHadron.h"
//btagging
#include "VHbbAnalysis/VHbbDataFormats/interface/BTagWeight.h"
//trigger
#include "VHbbAnalysis/VHbbDataFormats/interface/TriggerWeight.h"

#include <sstream>
#include <string>

#define MAXJ 30
#define MAXL 10
#define MAXB 10

struct CompareDeltaR {
  CompareDeltaR(TLorentzVector p4dir_): p4dir(p4dir_) {}
  bool operator()( const VHbbEvent::SimpleJet& j1, const  VHbbEvent::SimpleJet& j2 ) const {
    return j1.p4.DeltaR(p4dir) > j2.p4.DeltaR(p4dir);
  }
  TLorentzVector p4dir;
};

const GlobalVector flightDirection(const TVector3 pv, const reco::Vertex &sv){
  GlobalVector fdir(sv.position().X() - pv.X(),
                    sv.position().Y() - pv.Y(),
                    sv.position().Z() - pv.Z());
  return fdir;
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


//FIXME : need to update EDM ntuple with SimBhadron infos
// typedef struct
// {
//   void set(const SimBHadron & sbhc, int i){
//     mass[i] = sbhc.mass();
//     pt[i] = sbhc.pt();
//     eta[i] = sbhc.eta();
//     phi[i] = sbhc.phi();
//     vtx_x[i] = sbhc.decPosition.x();
//     vtx_y[i] = sbhc.decPosition.y();
//     vtx_z[i] = sbhc.decPosition.z();
//     pdgId[i] = sbhc.pdgId();
//     status[i] = sbhc.status();
//   };
//   void reset(){
//     for(int i=0; i < MAXB; ++i){
//       mass[i] = -99; pt[i] = -99; eta[i] = -99; phi[i] = -99; vtx_x[i] = -99; vtx_y[i] = -99; vtx_z[i] = -99; pdgId[i] = -99; status[i] = -99;
//     }
//   };
//   float mass[MAXB];
//   float pt[MAXB];
//   float eta[MAXB];
//   float phi[MAXB];
//   float vtx_x[MAXB];
//   float vtx_y[MAXB];
//   float vtx_z[MAXB];
//   int pdgId[MAXB];
//   int status[MAXB];
// //   int quarkStatus[MAXB];
// //   int brotherStatus[MAXB];
// //   int otherId[MAXB];
// //   bool etaOk[MAXB];
// //   bool simOk[MAXB];
// //   bool trackOk[MAXB];
// //   bool cutOk[MAXB];
// //   bool cutNewOk[MAXB];
// //   bool mcMatchOk[MAXB];
// //   bool matchOk[MAXB];
// } SimBHadronInfo;


typedef struct
{
  void set( const reco::SecondaryVertex & recoSv, const TVector3 recoPv, int isv){
    pt[isv]   = recoSv.p4().Pt();
    eta[isv]  = flightDirection(recoPv,recoSv).eta();
    phi[isv]  = flightDirection(recoPv,recoSv).phi();
    massBcand[isv] = recoSv.p4().M();
    massSv[isv] = recoSv.p4().M();
    dist3D[isv] = recoSv.dist3d().value();  
    distSig3D[isv] = recoSv.dist3d().significance(); 
    dist2D[isv] = recoSv.dist2d().value();  
    distSig2D[isv] = recoSv.dist2d().significance(); 
    dist3D_norm[isv] = recoSv.dist3d().value()/recoSv.p4().Gamma(); 
  };
  void reset(){
    for(int i = 0; i < MAXB; ++i){
      massBcand[i] = -99; massSv[i]= -99; pt[i] = -99; eta[i] = -99; phi[i] = -99; dist3D[i] = -99; distSig3D[i] = -99; dist2D[i] = -99; distSig2D[i] = -99; dist3D_norm[i] = -99;
    }
  };
  float massBcand[MAXB];
  float massSv[MAXB];
  float pt[MAXB]; 
  float eta[MAXB]; 
  float phi[MAXB];
  float dist3D[MAXB]; 
  float distSig3D[MAXB];
  float dist2D[MAXB]; 
  float distSig2D[MAXB];
  float dist3D_norm[MAXB];
} IVFInfo;


typedef struct 
{
  float mass; 
  float pt;
  float eta;
  float phi;
  float dR;
  float dPhi;
  float dEta;
} HiggsInfo;
  
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
    setSpecific(i,j);
  }
  template <class Input> void setSpecific(const Input & i, int j)
  {
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
  float id80[MAXL];
  float id95[MAXL];
};
  
template <> void LeptonInfo::setSpecific<VHbbEvent::ElectronInfo>(const VHbbEvent::ElectronInfo & i, int j){
  id80[j]=i.id80;
  id95[j]=i.id95;
}
template <> void LeptonInfo::setSpecific<VHbbEvent::MuonInfo>(const VHbbEvent::MuonInfo & i, int j){
  dxy[j]=i.ipDb;
  dz[j]=i.zPVPt;
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
    numTracksSV[i] = j.vtxMass;
    vtxMass[i]= j.vtxMass;
    vtx3dL[i] = j.vtx3dL;
    vtx3deL[i] = j.vtx3deL;
    chf[i]=j.chargedHadronEFraction;
    nhf[i]  =j.neutralHadronEFraction;
    cef[i]  =j.chargedEmEFraction;
    nef[i]  =j.neutralEmEFraction;
    nconstituents[i]  =j.nConstituents;
    nch[i]=j.ntracks;

    flavour[i]=j.flavour;
    if(j.bestMCp4.Pt() > 0)
      {
	genPt[i]=j.bestMCp4.Pt();
	genEta[i]=j.bestMCp4.Eta();
	genPhi[i]=j.bestMCp4.Phi();
      }
    JECUnc[i]=j.jecunc;
    id[i]=jetId(i);
  }
  bool jetId(int i)
  {
    if(nhf[i] > 0.99) return false;
    if(nef[i] > 0.99) return false;
    if(nconstituents[i]  <= 1) return false;
    if(fabs(eta[i])<2.4) {
    if(cef[i] > 0.99) return false;
    if(chf[i] == 0) return false;
    if(nch[i]== 0) return false;
    }
    return true;
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
  float vtxMass[MAXJ];
  float vtx3dL [MAXJ];
  float vtx3deL[MAXJ];
  bool id[MAXJ];
} JetInfo;
  
int main(int argc, char* argv[]) 
{
  gROOT->Reset();

  TTree *_outTree;
  IVFInfo IVF;
  //FIXME
  //  SimBHadronInfo SimBs;
  float rho,rho25;
  int nPVs;
  METInfo MET;
  METInfo METnoPU;
  MHTInfo MHT;
  TopInfo top;
  EventInfo EVENT;
  //  JetInfo jet1,jet2, addJet1, addJet2;
  // lepton1,lepton2;
  JetInfo hJets, aJets;
  LeptonInfo vLeptons, aLeptons;
  int naJets=0, nhJets=0;
  HiggsInfo H,SVH; //add here the fatjet higgs
  TrackInfo V;
  int nvlep=0,nalep=0; 
  
  float HVdPhi,HVMass,HMETdPhi,VMt,deltaPullAngle,deltaPullAngleAK7,gendrcc,gendrbb, genZpt, genWpt, weightTrig, weightTrigMay,weightTrigV4, weightTrigMET, minDeltaPhijetMET,  jetPt_minDeltaPhijetMET , PUweight;
  float weightEleRecoAndId,weightEleTrigJetMETPart, weightEleTrigElePart;

  //FIXME
  //  float SimBs_dr = -99, SimBs_dEta= -99, SimBs_dPhi = -99, SimBs_Hmass = -99, SimBs_Hpt = -99, SimBs_Heta = -99, SimBs_Hphi = -99;

  int Vtype,nSvs,nSimBs,numJets,numBJets,eventFlav;
  //   bool isMET80_CJ80, ispfMHT150, isMET80_2CJ20,isMET65_2CJ20, isJETID,isIsoMu17;
  bool triggerFlags[500],hbhe;
  float btag1T2CSF=1.,btag2TSF=1.,btag1TSF=1.,btagA0CSF=1., btag2CSF=1.;
  // ----------------------------------------------------------------------
  // First Part: 
  //
  //  * enable the AutoLibraryLoader 
  //  * book the histograms of interest 
  //  * open the input file
  // ----------------------------------------------------------------------

  // load framework libraries
  gSystem->Load("libFWCoreFWLite");
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
  int skipEvents_( in.getParameter<int>("skipEvents") );
  unsigned int outputEvery_( in.getParameter<unsigned int>("outputEvery") );
  std::string outputFile_( out.getParameter<std::string>("fileName" ) );
  std::vector<std::string> triggers( ana.getParameter<std::vector<std::string> >("triggers") );
  double btagThr =  ana.getParameter<double>("bJetCountThreshold" );
  bool fromCandidate = ana.getParameter<bool>("readFromCandidates");
  bool useHighestPtHiggsZ = ana.getParameter<bool>("useHighestPtHiggsZ");
  bool useHighestPtHiggsW = ana.getParameter<bool>("useHighestPtHiggsW");
  HbbCandidateFinderAlgo * algoZ = new HbbCandidateFinderAlgo(ana.getParameter<bool>("verbose"), ana.getParameter<double>("jetPtThresholdZ"),useHighestPtHiggsZ);
  HbbCandidateFinderAlgo * algoW = new HbbCandidateFinderAlgo(ana.getParameter<bool>("verbose"), ana.getParameter<double>("jetPtThresholdW"),useHighestPtHiggsW );

  TriggerWeight triggerWeight(ana);
  BTagWeight btag(2); // 2 operating points "Custom" = 0.5 and "Tight = 0.898"
  BTagSampleEfficiency btagEff( ana.getParameter<std::string>("btagEffFileName" ).c_str() ); 

  std::vector<std::string> inputFiles_( in.getParameter<std::vector<std::string> >("fileNames") );
  //  std::string inputFile( in.getParameter<std::string> ("fileName") );

  std::string PUmcfileName_ = in.getParameter<std::string> ("PUmcfileName") ;
  std::string PUdatafileName_ = in.getParameter<std::string> ("PUdatafileName") ;

  bool isMC_( ana.getParameter<bool>("isMC") );  
  TriggerReader trigger(isMC_);
  TriggerReader patFilters(false);

  edm::LumiReWeighting   lumiWeights;
  if(isMC_)
    {
        	   lumiWeights = edm::LumiReWeighting(PUmcfileName_,PUdatafileName_ , "pileup", "pileup");
    }
 
  //   TFile *_outPUFile	= new TFile((outputFile_+"_PU").c_str(), "recreate");	
  TFile *_outFile	= new TFile(outputFile_.c_str(), "recreate");	
  TH1F *  count = new TH1F("Count","Count", 1,0,2 );
  TH1F *  countWithPU = new TH1F("CountWithPU","CountWithPU", 1,0,2 );
  TH1F * pu = new TH1F("pileup","",51,-0.5,50.5);
  _outTree = new TTree("tree", "myTree");
  
  _outTree->Branch("H"		,  &H	            ,  "mass/F:pt/F:eta:phi/F:dR/F:dPhi/F:dEta/F");
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
  _outTree->Branch("hJet_vtxMass",hJets.vtxMass ,"vtxMass[nhJets]/F");
  _outTree->Branch("hJet_vtx3dL",hJets.vtx3dL ,"vtx3dL[nhJets]/F");
  _outTree->Branch("hJet_vtx3deL",hJets.vtx3deL ,"vtx3deL[nhJets]/F");
  _outTree->Branch("hJet_id",hJets.id ,"id[nhJets]/b");

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
  _outTree->Branch("aJet_vtxMass",aJets.vtxMass ,"vtxMass[naJets]/F");
  _outTree->Branch("aJet_vtx3dL",aJets.vtx3dL ,"vtx3dL[naJets]/F");
  _outTree->Branch("aJet_vtx3deL",aJets.vtx3deL ,"vtx3deL[naJets]/F");
  _outTree->Branch("aJet_id",aJets.id ,"id[naJets]/b");

  _outTree->Branch("numJets"      ,  &numJets         ,  "numJets/I"       );                
  _outTree->Branch("numBJets"      ,  &numBJets         ,  "numBJets/I"       );                
  _outTree->Branch("deltaPullAngle", &deltaPullAngle  ,  "deltaPullAngle/F");
  _outTree->Branch("gendrcc"    , &gendrcc      ,  "gendrcc/F");
  _outTree->Branch("gendrbb"    , &gendrbb      ,  "gendrbb/F");
  _outTree->Branch("genZpt"    , &genZpt      ,  "genZpt/F");
  _outTree->Branch("genWpt"    , &genWpt      ,  "genWpt/F");
  _outTree->Branch("weightTrig"        , &weightTrig          ,  "weightTrig/F");
  _outTree->Branch("weightTrigMay"        , &weightTrigMay          ,  "weightTrigMay/F");
  _outTree->Branch("weightTrigV4"        , &weightTrigV4          ,  "weightTrigV4/F");
  _outTree->Branch("weightTrigMET"        , &weightTrigMET          ,  "weightTrigMET/F");
  _outTree->Branch("weightEleRecoAndId"        , &weightEleRecoAndId     ,  "weightEleRecoAndId/F");
  _outTree->Branch("weightEleTrigJetMETPart"        , &weightEleTrigJetMETPart          ,  "weightEleTrigJetMETPart/F");
  _outTree->Branch("weightEleTrigElePart"        , &weightEleTrigElePart          ,  "weightEleTrigElePart/F");


  _outTree->Branch("deltaPullAngleAK7", &deltaPullAngleAK7  ,  "deltaPullAngleAK7/F");
  _outTree->Branch("PUweight",       &PUweight  ,  "PUweight/F");
  _outTree->Branch("eventFlav",       &eventFlav  ,  "eventFlav/I");
 

   
    
  _outTree->Branch("Vtype"     ,  &Vtype   ,   "Vtype/I" );                
  _outTree->Branch("HVdPhi"     ,  &HVdPhi   ,   "HVdPhi/F" );                
  _outTree->Branch("HVMass"     ,  &HVMass   ,   "HVMass/F" );                
  _outTree->Branch("HMETdPhi"     ,  &HMETdPhi   ,   "HMETdPhi/F" );                
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
  _outTree->Branch("vLepton_id80",vLeptons.id80 ,"id80[nvlep]/F");
  _outTree->Branch("vLepton_id95",vLeptons.id95 ,"id95[nvlep]/F");
 
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
  _outTree->Branch("aLepton_id80",aLeptons.id80 ,"id[nalep]/F");
  _outTree->Branch("aLepton_id95",aLeptons.id95 ,"id[nalep]/F");
 
  _outTree->Branch("top"		,  &top	         ,   "mass/F:pt/F:wMass/F");

  //IVF
  _outTree->Branch("nSvs",&nSvs ,"nSvs/I");
  _outTree->Branch("Sv_massBCand", &IVF.massBcand,"massBcand[nSvs]/F");
  _outTree->Branch("Sv_massSv", &IVF.massSv,"massSv[nSvs]/F");
  _outTree->Branch("Sv_pt", &IVF.pt,"pt[nSvs]/F");
  _outTree->Branch("Sv_eta", &IVF.eta,"eta[nSvs]/F");
  _outTree->Branch("Sv_phi", &IVF.phi,"phi[nSvs]/F");
  _outTree->Branch("Sv_dist3D", &IVF.dist3D,"dist3D[nSvs]/F");
  _outTree->Branch("Sv_dist2D", &IVF.dist2D,"dist2D[nSvs]/F");
  _outTree->Branch("Sv_distSim2D", &IVF.distSig2D,"distSig2D[nSvs]/F");
  _outTree->Branch("Sv_distSig3D", &IVF.distSig3D,"distSig3D[nSvs]/F");
  _outTree->Branch("Sv_dist3D_norm", &IVF.dist3D_norm,"dist3D_norm[nSvs]/F");
  //IVF higgs candidate
  _outTree->Branch("SVH"          ,  &SVH               ,  "mass/F:pt/F:eta:phi/F:dR/F:dPhi/F:dEta/F");


  //FIXME : need to update the EDMNtuple with BHadron infos
//   //SimBHadron
//   _outTree->Branch("nSimBs",&nSimBs ,"nSimBs/I");
//   _outTree->Branch("SimBs_mass", &SimBs.mass,"mass[nSvs]/F");
//   _outTree->Branch("SimBs_pt", &SimBs.pt,"pt[nSvs]/F");
//   _outTree->Branch("SimBs_eta", &SimBs.eta,"eta[nSvs]/F");
//   _outTree->Branch("SimBs_phi", &SimBs.phi,"phi[nSvs]/F");
//   _outTree->Branch("SimBs_vtx_x", &SimBs.vtx_x,"vtx_x[nSvs]/F");
//   _outTree->Branch("SimBs_vtx_y", &SimBs.vtx_y,"vtx_y[nSvs]/F");
//   _outTree->Branch("SimBs_vtx_z", &SimBs.vtx_z,"vtx_z[nSvs]/F");
//   _outTree->Branch("SimBs_pdgId", &SimBs.pdgId,"pdgId[nSvs]/F");
//   _outTree->Branch("SimBs_status", &SimBs.status,"status[nSvs]/F");
//   //SimBHadron Higgs Candidate
//   _outTree->Branch("SimBs_dr", &SimBs_dr,"SimBs_dr/F");
//   _outTree->Branch("SimBs_dPhi", &SimBs_dPhi,"SimBs_dPhi/F");
//   _outTree->Branch("SimBs_dEta", &SimBs_dEta,"SimBs_dEta/F");
//   _outTree->Branch("SimBs_Hmass", &SimBs_Hmass,"SimBs_Hmass/F");
//   _outTree->Branch("SimBs_Hpt", &SimBs_Hpt,"SimBs_Hpt/F");
//   _outTree->Branch("SimBs_Heta", &SimBs_Heta,"SimBs_Heta/F");
//   _outTree->Branch("SimBs_Hphi", &SimBs_Hphi,"SimBs_Hphi/F");


  _outTree->Branch("rho"		,  &rho	         ,   "rho/F");
  _outTree->Branch("rho25"		,  &rho25	         ,   "rho25/F");
  _outTree->Branch("nPVs"		,  &nPVs	         ,   "nPVs/I");
  _outTree->Branch("METnoPU"		,  &METnoPU	         ,   "et/F:sumet:sig/F:phi/F");
  _outTree->Branch("MET"		,  &MET	         ,   "et/F:sumet:sig/F:phi/F");
  _outTree->Branch("MHT"		,  &MHT	         ,   "mht/F:ht:sig/F:phi/F");
  _outTree->Branch("minDeltaPhijetMET"		,  &minDeltaPhijetMET	         ,   "minDeltaPhijetMET/F");
  _outTree->Branch("jetPt_minDeltaPhijetMET"		,  &jetPt_minDeltaPhijetMET	         ,   "jetPt_minDeltaPhijetMET/F");

  std::stringstream s;
  s << "triggerFlags[" << triggers.size() << "]/b";
  _outTree->Branch("triggerFlags", triggerFlags, s.str().c_str()); 
  
 
  _outTree->Branch("EVENT"		,  &EVENT	         ,   "run/I:lumi/I:event/I:json/I");
  _outTree->Branch("hbhe"		,  &hbhe	         ,   "hbhe/b");
  _outTree->Branch("btag1TSF"		,  &btag1TSF	         ,   "btag1TSF/F");
  _outTree->Branch("btag2TSF"		,  &btag2TSF	         ,   "btag2TSF/F");
  _outTree->Branch("btag1T2CSF"	,  &btag1T2CSF	         ,   "btag1T2CSF/F");
  _outTree->Branch("btag2CSF"	,  &btag2CSF	         ,   "btag2CSF/F");
  _outTree->Branch("btagA0CSF"	,  &btagA0CSF	         ,   "btagA0CSF/F");

  int ievt=0;  
  int totalcount=0;

  //  TFile* inFile = new TFile(inputFile.c_str(), "read");
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile) {
    std::cout << iFile << std::endl;
    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
    if(inFile==0) continue;

    // loop the events
      
      fwlite::Event ev(inFile);
      for(ev.toBegin(); !ev.atEnd() ; ++ev, ++ievt)
        {
          if (ievt <= skipEvents_) continue;
          if (maxEvents_ >= 0){
              if (ievt > maxEvents_ + skipEvents_) break;
          };
          count->Fill(1.);

      fwlite::Handle< VHbbEventAuxInfo > vhbbAuxHandle; 
      vhbbAuxHandle.getByLabel(ev,"HbbAnalyzerNew");
      const VHbbEventAuxInfo & aux = *vhbbAuxHandle.product();
      PUweight=1.;
 	  if(isMC_){
 
 	  // PU weights
	  std::map<int, unsigned int>::const_iterator puit = aux.puInfo.pus.find(0);
	  PUweight =  lumiWeights.weight( puit->second );        
	  pu->Fill(puit->second);
	}
	countWithPU->Fill(PUweight);
      
	//Write event info 
	EVENT.run = ev.id().run();
	EVENT.lumi = ev.id().luminosityBlock();
	EVENT.event = ev.id().event();
	EVENT.json = jsonContainsEvent (jsonVector, ev);

	//FIXME : need to update EDM ntuple with BHadron infos
// 	// simBHadrons
// 	fwlite::Handle<SimBHadronCollection> SBHC;
// 	SBHC.getByLabel(ev, "bhadrons");
// 	const SimBHadronCollection sbhc = *(SBHC.product());
	
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

	    /*     for(size_t m=0;m<iEvent->muInfo.size();m++)
		   { 

		   if( fabs(iEvent->muInfo[m].p4.Pt()-28.118684) < 0.0001 ||
		   fabs(iEvent->muInfo[m].p4.Pt()-34.853199) < 0.0001 )  
		   {
		   std::cout << "FOUND " << iEvent->muInfo[m].p4.Pt() <<  " " << EVENT.event << " " << candW->size() << " " << candZ->size() << std::endl;
		   }
		   }
	    */
 
	  }

        const std::vector<VHbbCandidate> * cand = candZ;


      
	/*  fwlite::Handle< VHbbEvent > vhbbHandle; 
	    vhbbHandle.getByLabel(ev,"HbbAnalyzerNew");
	    const VHbbEvent iEvent = *vhbbHandle.product();
	*/

	//      std::clog << "Filling tree "<< std::endl;
	bool isW=false;

 
	if(cand->size() == 0 or cand->at(0).H.jets.size() < 2) continue;
	if(cand->size() > 1 ) 
          {
	    std::cout << "MULTIPLE CANDIDATES: " << cand->size() << std::endl;
          }
	if(cand->at(0).candidateType == VHbbCandidate::Wmun || cand->at(0).candidateType == VHbbCandidate::Wen ) { cand=candW; isW=true; }
	if(cand->size() == 0) 
          {
	    //            std::cout << "W event loss due to tigther cuts" << std::endl;
            continue;
          }

	// secondary vtxs
	fwlite::Handle<std::vector<reco::Vertex> > SVC;
	SVC.getByLabel(ev,"selectedVertices");
	const std::vector<reco::Vertex> svc = *(SVC.product());

	const VHbbCandidate & vhCand =  cand->at(0);
	patFilters.setEvent(&ev,"VH");
	hbhe = patFilters.accept("hbhe");

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
        if(isW) V.mass = vhCand.Mt(); 
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
	H.dR = deltaR(vhCand.H.jets[0].p4.Eta(),vhCand.H.jets[0].p4.Phi(),vhCand.H.jets[1].p4.Eta(),vhCand.H.jets[1].p4.Phi());
	H.dPhi = deltaPhi(vhCand.H.jets[0].p4.Phi(),vhCand.H.jets[1].p4.Phi());
	H.dEta= TMath::Abs( vhCand.H.jets[0].p4.Eta() - vhCand.H.jets[1].p4.Eta() );
	HVdPhi = fabs( deltaPhi(vhCand.H.p4.Phi(),vhCand.V.p4.Phi()) ) ;
	HVMass = (vhCand.H.p4 + vhCand.V.p4).M() ;
	HMETdPhi = fabs( deltaPhi(vhCand.H.p4.Phi(),vhCand.V.mets.at(0).p4.Phi()) ) ;
	VMt = vhCand.Mt() ;
	deltaPullAngle = vhCand.H.deltaTheta;
  
	hJets.cosTheta[0]=  vhCand.H.helicities[0];
	hJets.cosTheta[1]=  vhCand.H.helicities[1];
// METInfo calomet;  METInfo tcmet;  METInfo pfmet;  METInfo mht;  METInfo metNoPU

	MET.et = vhCand.V.mets.at(0).p4.Pt();
	MET.phi = vhCand.V.mets.at(0).p4.Phi();
	MET.sumet = vhCand.V.mets.at(0).sumEt;
	MET.sig = vhCand.V.mets.at(0).metSig;

	METnoPU.et = iEvent->metNoPU.p4.Pt();
	METnoPU.phi = iEvent->metNoPU.p4.Phi();
	METnoPU.sumet = iEvent->metNoPU.sumEt;
	METnoPU.sig = iEvent->metNoPU.metSig;

        rho = aux.puInfo.rho;
        rho25 = aux.puInfo.rho25;
        nPVs=aux.pvInfo.nVertices; 
 
        if(!fromCandidate) {
	  MHT.mht = iEvent->mht.p4.Pt(); 
	  MHT.phi = iEvent->mht.p4.Phi(); 
	  MHT.ht = iEvent->mht.sumEt; 
	  MHT.sig = iEvent->mht.metSig; 
	}
 
	Vtype = vhCand.candidateType;


	//Secondary Vertices
	IVF.reset();
	nSvs = svc.size();
	const TVector3 recoPv = aux.pvInfo.firstPVInPT2;
	const math::XYZPoint myPv(recoPv);
	//look here for Matrix filling info http://project-mathlibs.web.cern.ch/project-mathlibs/sw/html/SMatrixDoc.html
	std::vector<double> fillMatrix(6);
	for (int i = 0; i<6; ++i) fillMatrix[i] = 0.;
	fillMatrix[0] = TMath::Power(0.002,2);
	fillMatrix[2] = TMath::Power(0.002,2);
	fillMatrix[5] = TMath::Power(0.002,2);
	const ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > myFakeMatrixError(fillMatrix.begin(),fillMatrix.end());
	const reco::Vertex recoVtxPv(myPv, myFakeMatrixError);
	for( int j=0; j < nSvs && j < MAXB; ++j ) {
	  const GlobalVector flightDir = flightDirection(recoPv,svc[j]);
	  reco::SecondaryVertex recoSv(recoVtxPv, svc[j], flightDir ,true);
	  IVF.set( recoSv, recoPv ,j);
	}
	if(nSvs > 2){
	  TLorentzVector BCands_H1, BCands_H2, BCands_H;
	  BCands_H1.SetPtEtaPhiM(IVF.pt[0], IVF.eta[0], IVF.phi[0], IVF.massBcand[0]);
	  BCands_H2.SetPtEtaPhiM(IVF.pt[1], IVF.eta[1], IVF.phi[1], IVF.massBcand[1]);
	  BCands_H = BCands_H1 + BCands_H2;
	  SVH.dR = deltaR(IVF.eta[0], IVF.phi[0], IVF.eta[1], IVF.phi[1] );
	  SVH.dPhi = deltaPhi(IVF.phi[0], IVF.phi[1] );
	  SVH.dEta = TMath::Abs(IVF.eta[0] - IVF.eta[1] );
	  SVH.mass = BCands_H.M();
	  SVH.pt = BCands_H.Pt();
	  SVH.eta = BCands_H.Eta();
	  SVH.phi = BCands_H.Phi();
	}

	//FIXME : need to update EDM ntuple with simBhadron info
// 	//SimBhadron
// 	SimBs.reset();
// 	nSimBs = sbhc.size();
// 	for( int j=0; j < nSimBs && j < MAXB; ++j )
// 	  SimBs.set( sbhc.at(j), j);
// 	if(nSimBs > 2){
// 	  TLorentzVector SimBs_H1, SimBs_H2, SimBs_H;
// 	  SimBs_H1.SetPtEtaPhiM(SimBs.pt[0], SimBs.eta[0], SimBs.phi[0], SimBs.mass[0]);
// 	  SimBs_H2.SetPtEtaPhiM(SimBs.pt[1], SimBs.eta[1], SimBs.phi[1], SimBs.mass[1]);
// 	  SimBs_H = SimBs_H1 + SimBs_H2;
// 	  SimBs_dr = deltaR(SimBs.eta[0], SimBs.phi[0], SimBs.eta[1], SimBs.phi[1] );
// 	  SimBs_dPhi = deltaPhi(SimBs.phi[0], SimBs.phi[1] );
// 	  SimBs_dEta = TMath::Abs(SimBs.eta[0] - SimBs.eta[1] );
// 	  SimBs_Hmass = SimBs_H.M();
// 	  SimBs_Hpt = SimBs_H.Pt();
// 	  SimBs_Heta = SimBs_H.Eta();
// 	  SimBs_Hphi = SimBs_H.Phi();
// 	}


	//Loop on jets
	double maxBtag=-99999;
	minDeltaPhijetMET = 999;
	TLorentzVector bJet;
	std::vector<std::vector<BTagWeight::JetInfo> > btagJetInfos;
	std::vector<float> jet30eta;
	std::vector<float> jet30pt;
	if(fromCandidate)
	  {
	    //Loop on Higgs Jets
	    for(unsigned int j=0; j < vhCand.H.jets.size(); j++ ){
	      if (vhCand.H.jets[j].csv > maxBtag) { bJet=vhCand.H.jets[j].p4 ; maxBtag =vhCand.H.jets[j].csv; }
	      if (fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), vhCand.H.jets[j].p4.Phi())) < minDeltaPhijetMET) 
		{
                  minDeltaPhijetMET=fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), vhCand.H.jets[j].p4.Phi())); 
                  jetPt_minDeltaPhijetMET=vhCand.H.jets[j].p4.Pt();
		}
	      btagJetInfos.push_back(btagEff.jetInfo(vhCand.H.jets[j]));
	    }
	    //Loop on Additional Jets
	    for(unsigned int j=0; j < vhCand.additionalJets.size(); j++ ){
	      if (vhCand.additionalJets[j].csv > maxBtag) { bJet=vhCand.additionalJets[j].p4 ; maxBtag =vhCand.additionalJets[j].csv; }
/*	      if (fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), vhCand.additionalJets[j].p4.Phi())) < minDeltaPhijetMET) 
		{
                  minDeltaPhijetMET=fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), vhCand.additionalJets[j].p4.Phi()));
                  jetPt_minDeltaPhijetMET=vhCand.additionalJets[j].p4.Pt();
		}*/
	      if( ( isW && ! useHighestPtHiggsW ) ||  ( ! isW && ! useHighestPtHiggsZ )  )  // btag SF computed using only H-jets if best-H made with dijetPt rather than best CSV
		{
		  if(vhCand.additionalJets[j].p4.Pt() > 20)
		    btagJetInfos.push_back(btagEff.jetInfo(vhCand.additionalJets[j]));
		}
	    }
	  } 
	else
	  {
            //Loop on all jets
	    for(unsigned int j=0; j < iEvent->simpleJets2.size(); j++ ){
	      if (iEvent->simpleJets2[j].csv > maxBtag) { bJet=iEvent->simpleJets2[j].p4 ; maxBtag =iEvent->simpleJets2[j].csv; }
	      if ( iEvent->simpleJets2[j].p4.Pt() > 20 &&  fabs(iEvent->simpleJets2[j].p4.Eta()) < 2.4&& fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), iEvent->simpleJets2[j].p4.Phi())) < minDeltaPhijetMET)
		{
                  minDeltaPhijetMET=fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), iEvent->simpleJets2[j].p4.Phi()));
                  jetPt_minDeltaPhijetMET=iEvent->simpleJets2[j].p4.Pt();
		}
	      if(iEvent->simpleJets2[j].p4.Pt() > 30)
		{
		  jet30eta.push_back(iEvent->simpleJets2[j].p4.Eta());
		  jet30pt.push_back(iEvent->simpleJets2[j].p4.Pt());
		}
   
	      //For events made with highest CSV, all jets in the event should be taken into account for "tagging" SF (anti tagging is a mess)
	      // because for example a light jet not used for the Higgs can have in reality a higher CSV due to SF > 1 and become a higgs jet
	      if( ( isW && ! useHighestPtHiggsW ) ||  ( ! isW && ! useHighestPtHiggsZ )  ) 
		{ 
		  if(iEvent->simpleJets2[j].p4.Pt() > 20)
		    btagJetInfos.push_back(btagEff.jetInfo(iEvent->simpleJets2[j]));
		}
	    }

            //Loop on Higgs jets

	    for(unsigned int j=0; j < vhCand.H.jets.size(); j++ ) {

	    //if we use the highest pt pair, only the two higgs jet should be used to compute the SF because the other jets are excluded 
	    // by a criteria (pt of the dijet) that is not btag SF dependent 
	    if(!( ( isW && ! useHighestPtHiggsW ) ||  ( ! isW && ! useHighestPtHiggsZ ) )) {
		   btagJetInfos.push_back(btagEff.jetInfo(vhCand.H.jets[j]));
	    }
           }

	  }
	vLeptons.reset();
	weightTrig = 1.; // better to default to 1 
	weightTrigMay = -1.;
	weightTrigV4 = -1.; 
	TLorentzVector leptonForTop;
	size_t firstAddMu=0;
	size_t firstAddEle=0;
	if(Vtype == VHbbCandidate::Zmumu ){
	  vLeptons.set(vhCand.V.muons[0],0,13); 
	  vLeptons.set(vhCand.V.muons[1],1,13);
	  float cweightID = triggerWeight.scaleMuID(vLeptons.pt[0],vLeptons.eta[0]) * triggerWeight.scaleMuID(vLeptons.pt[1],vLeptons.eta[1]) ;
	  float weightTrig1 = triggerWeight.scaleMuIsoHLT(vLeptons.pt[0],vLeptons.eta[0]);
	  float weightTrig2 = triggerWeight.scaleMuIsoHLT(vLeptons.pt[1],vLeptons.eta[1]);
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
	  std::vector<float> pt,eta;
	  pt.push_back(vLeptons.pt[0]); eta.push_back(vLeptons.eta[0]);
	  pt.push_back(vLeptons.pt[1]); eta.push_back(vLeptons.eta[1]);
	  weightEleRecoAndId=triggerWeight.scaleID95Ele(vLeptons.pt[0],vLeptons.eta[0]) * triggerWeight.scaleRecoEle(vLeptons.pt[0],vLeptons.eta[0]) *
	    triggerWeight.scaleID95Ele(vLeptons.pt[1],vLeptons.eta[1]) * triggerWeight.scaleRecoEle(vLeptons.pt[1],vLeptons.eta[1]);
	  weightEleTrigElePart = triggerWeight.scaleDoubleEle17Ele8(pt,eta); 
	  weightTrig = weightEleTrigElePart * weightEleRecoAndId;
 
                 

	}
	if(Vtype == VHbbCandidate::Wmun ){
	  leptonForTop=vhCand.V.muons[0].p4;
	  vLeptons.set(vhCand.V.muons[0],0,13); 
	  float cweightID = triggerWeight.scaleMuID(vLeptons.pt[0],vLeptons.eta[0]);
	  float weightTrig1 = triggerWeight.scaleMuIsoHLT(vLeptons.pt[0],vLeptons.eta[0]);
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
	  weightTrigMay = triggerWeight.scaleSingleEleMay(vLeptons.pt[0],vLeptons.eta[0]);
	  weightTrigV4 = triggerWeight.scaleSingleEleV4(vLeptons.pt[0],vLeptons.eta[0]);
	  weightEleRecoAndId=triggerWeight.scaleID80Ele(vLeptons.pt[0],vLeptons.eta[0]) * triggerWeight.scaleRecoEle(vLeptons.pt[0],vLeptons.eta[0]);
	  weightEleTrigJetMETPart=triggerWeight.scaleJet30Jet25(jet30eta,jet30pt)*triggerWeight.scalePFMHTEle(MET.et);
	  weightEleTrigElePart= weightTrigV4; //this is for debugging only, checking only the V4 part

	  weightTrigMay*=weightEleRecoAndId;
	  weightTrigV4*=weightEleRecoAndId;
	  weightTrigV4*=weightEleTrigJetMETPart;
	  weightTrig = weightTrigMay * 0.187 + weightTrigV4 * (1.-0.187); //FIXME: use proper lumi if we reload 2.fb
	}
	if( Vtype == VHbbCandidate::Znn ){
	  nvlep=0;
	  float weightTrig1 = triggerWeight.scaleMetHLT(vhCand.V.mets.at(0).p4.Pt());
	  weightTrig = weightTrig1;
	  weightTrigMET = weightTrig1;  
	}
      
        if(weightTrigMay < 0) weightTrigMay=weightTrig;
        if(weightTrigV4 < 0) weightTrigV4=weightTrig;
        if(!isMC_)
        {
         weightTrig = 1.; 
         weightTrigMay = 1.;
         weightTrigV4 = 1.;

        } 
	aLeptons.reset();
	nalep=0;
	if(fromCandidate)
          {
            for(size_t j=firstAddMu;j< vhCand.V.muons.size();j++) aLeptons.set(vhCand.V.muons[j],nalep++,13);
            for(size_t j=firstAddEle;j< vhCand.V.electrons.size();j++) aLeptons.set(vhCand.V.electrons[j],nalep++,11);
          }
	else
          {
	    for(size_t j=0;j< iEvent->muInfo.size();j++)
	      { 
                if((j!= vhCand.V.firstLepton && j!= vhCand.V.secondLepton) || ((Vtype != VHbbCandidate::Wmun ) && (Vtype != VHbbCandidate::Zmumu )) )
		  aLeptons.set(iEvent->muInfo[j],nalep++,13);
	      }
	    for(size_t j=0;j< iEvent->eleInfo.size();j++)
	      { 
                if((j!= vhCand.V.firstLepton && j!= vhCand.V.secondLepton) || ((Vtype != VHbbCandidate::Wen ) && (Vtype != VHbbCandidate::Zee )))
		  aLeptons.set(iEvent->eleInfo[j],nalep++,11);
	      }

          }


	if(isMC_)
	  {
	    //std::cout << "BTAGSF " <<  btagJetInfos.size() << " " << btag.weight<BTag1Tight2CustomFilter>(btagJetInfos) << std::endl;
	    if ( btagJetInfos.size()< 10) 
	      {
		btag1T2CSF = btag.weight<BTag1Tight2CustomFilter>(btagJetInfos);
		btag2TSF = btag.weight<BTag2TightFilter>(btagJetInfos);
		btag1TSF = btag.weight<BTag1TightFilter>(btagJetInfos);
		btagA0CSF = btag.weight<BTagAntiMax0CustomFilter>(btagJetInfos);
		btag2CSF = btag.weight<BTag2CustomFilter>(btagJetInfos);
	      }
	    else
	      {
		std::cout << "WARNING:  combinatorics for " << btagJetInfos.size() << " jets is too high (>=10). use SF=1 " << std::endl;
		//TODO: revert to random throw  for this cases
		btag1T2CSF = 1.;
		btag2TSF =  1.;
		btag1TSF =  1.;
	      }
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
  

       
	//FIXME: too much  warnings... figure out why 
	//         gendrcc=aux.genCCDeltaR(); 
	//         gendrbb=aux.genBBDeltaR(); 
	genZpt=aux.mcZ.size() > 0 ? aux.mcZ[0].p4.Pt():-99;
	genWpt=aux.mcW.size() > 0 ? aux.mcW[0].p4.Pt():-99;


        /// Compute pull angle from AK7
        if(!fromCandidate){
          std::vector<VHbbEvent::SimpleJet> ak7wrt1(iEvent->simpleJets3);
          std::vector<VHbbEvent::SimpleJet> ak7wrt2(iEvent->simpleJets3);
          if(ak7wrt1.size() > 1){
            CompareDeltaR deltaRComparatorJ1(vhCand.H.jets[0].p4);
            CompareDeltaR deltaRComparatorJ2(vhCand.H.jets[1].p4);
            std::sort( ak7wrt1.begin(),ak7wrt1.end(),deltaRComparatorJ1 );
            std::sort( ak7wrt2.begin(),ak7wrt2.end(),deltaRComparatorJ2 );
            std::vector<VHbbEvent::SimpleJet> ak7_matched;
            // if the matched are different save them
            if(ak7wrt1[0].p4.DeltaR(ak7wrt2[0].p4) > 0.1) {
              ak7_matched.push_back(ak7wrt1[0]);
              ak7_matched.push_back(ak7wrt2[0]);
            }
            // else look at the second best
            else{
              // ak7wrt1 is best
              if( ak7wrt1[1].p4.DeltaR(vhCand.H.jets[0].p4) < ak7wrt2[1].p4.DeltaR(vhCand.H.jets[1].p4))
                {
                  ak7_matched.push_back(ak7wrt1[1]);
                  ak7_matched.push_back(ak7wrt2[0]);
                }
              else
                {
                  ak7_matched.push_back(ak7wrt1[0]);
                  ak7_matched.push_back(ak7wrt2[1]);
                }
            }
            CompareJetPt ptComparator;
            std::sort( ak7_matched.begin(),ak7_matched.end(),ptComparator );
            if(ak7_matched[0].p4.DeltaR(vhCand.H.jets[0].p4) < 0.5
               and ak7_matched[1].p4.DeltaR(vhCand.H.jets[1].p4) < 0.5)
              {
                deltaPullAngleAK7 =  VHbbCandidateTools::getDeltaTheta(ak7_matched[0],ak7_matched[1]);
              }
          }
        }

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


