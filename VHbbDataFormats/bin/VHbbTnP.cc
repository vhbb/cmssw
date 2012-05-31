#include <TH1F.h>
#include <TH3F.h>
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
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

//for LHE info
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

//Move class definition to Ntupler.h ?
//#include "VHbbAnalysis/VHbbDataFormats/interface/Ntupler.h"

#include "ZSV/BAnalysis/interface/SimBHadron.h"
//btagging
#include "VHbbAnalysis/VHbbDataFormats/interface/BTagWeight.h"
//trigger
#include "VHbbAnalysis/VHbbDataFormats/interface/TriggerWeight.h"


//SCZ
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"
using namespace pat;
#define GENPTOLOR(a) TLorentzVector((a).px(), (a).py(), (a).pz(), (a).energy())
#define GENPTOLORP(a) TLorentzVector((a)->px(), (a)->py(), (a)->pz(), (a)->energy())
//#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
//#include "CLHEP/Random/RandFlat.h"
//#include "FWCore/Utilities/interface/Exception.h"
#include <TRandom3.h>

#include <sstream>
#include <string>

#define MAXJ 130
#define MAXL 110
#define MAXB 110
#define MAXT 160
#define nMetUnc 24 
//eleEnDown/Up, muEn, tauEn, JES, JER, Unclustered for type1 MET [0-11] and than type1p2 MET [12-23]


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

/*  
template <> void LeptonInfo::setSpecific<VHbbEvent::ElectronInfo>(const VHbbEvent::ElectronInfo & i, int j,const VHbbEventAuxInfo & aux){
  id80[j]=i.id80;
  id95[j]=i.id95;
  id80NoIso[j]=(i.innerHits ==0 && !(fabs(i.convDist)<0.02 && fabs(i.convDcot)<0.02) &&
((i.isEB && i.sihih<0.01 && fabs(i.Dphi)<0.06 && fabs(i.Deta)<0.004) || (i.isEE && i.sihih<0.03 && fabs(i.Dphi)<0.03  && fabs(i.Deta)<0.007)));

float mincor=0.0;
float minrho=0.0;
float rhoN = std::max(aux.puInfo.rhoNeutral,minrho);
float eta=i.p4.Eta();
float areagamma=0.5;
float areaNH=0.5;
float areaComb=0.5;
if(fabs(eta) <= 1.0 ) {areagamma=0.081; areaNH=0.024; areaComb=0.10;}
if(fabs(eta) > 1.0 &&  fabs(eta) <= 1.479 ) {areagamma=0.084; areaNH=0.037; areaComb=0.12;}
if(fabs(eta) > 1.479 &&  fabs(eta) <= 2.0 ) {areagamma=0.048; areaNH=0.037; areaComb=0.085;}
if(fabs(eta) > 2.0 &&  fabs(eta) <= 2.2 ) {areagamma=0.089; areaNH=0.023; areaComb=0.11;}
if(fabs(eta) > 2.2 &&  fabs(eta) <= 2.3 ) {areagamma=0.092; areaNH=0.023; areaComb=0.12;}
if(fabs(eta) > 2.3 &&  fabs(eta) <= 2.4 ) {areagamma=0.097; areaNH=0.021; areaComb=0.12;}
if(fabs(eta) > 2.4  ) {areagamma=0.11; areaNH=0.021; areaComb=0.13;}


pfCorrIso[j] = (i.pfChaIso+ std::max(i.pfPhoIso-rhoN*areagamma,mincor )+std::max(i.pfNeuIso-rhoN*areaNH,mincor))/i.p4.Pt();

id2012tight[j] = fabs(i.dxy) < 0.02  &&fabs(i.dz) < 0.1  &&(
(i.isEE  &&fabs(i.Deta) < 0.005 &&fabs(i.Dphi) < 0.02 &&i.sihih < 0.03  &&i.HoE < 0.10  &&fabs(i.fMVAVar_IoEmIoP) < 0.05
) ||
(i.isEB  &&fabs(i.Deta) < 0.004 &&fabs(i.Dphi) < 0.03 &&i.sihih < 0.01  &&i.HoE < 0.12  &&fabs(i.fMVAVar_IoEmIoP) < 0.05
 ));
float id=i.mvaOutTrig;
float iso=pfCorrIso[j];
wp70[j]=((fabs(eta) < 0.8 && id>0.977 && iso < 0.093) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.956 && iso < 0.095) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.966 && iso < 0.171));
wp80[j]=((fabs(eta) < 0.8 && id>0.913 && iso < 0.105) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.964 && iso < 0.178) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.899 && iso < 0.150));
wp85[j]=((fabs(eta) < 0.8 && id>0.929 && iso < 0.135) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.931 && iso < 0.159) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.805 && iso < 0.155));
wp90[j]=((fabs(eta) < 0.8 && id>0.877 && iso < 0.177) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.794 && iso < 0.180) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.846 && iso < 0.244));
wp95[j]=((fabs(eta) < 0.8 && id>0.858 && iso < 0.253) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.425 && iso < 0.225) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.759 && iso < 0.308));
wpHWW[j]=((fabs(eta) < 0.8 && id>0.94 && iso < 0.15) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.85 && iso < 0.15) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.92 && iso < 0.15));

 idMVAnotrig[j]=i.mvaOut;
 idMVAtrig[j]=i.mvaOutTrig;

}
template <> void LeptonInfo::setSpecific<VHbbEvent::MuonInfo>(const VHbbEvent::MuonInfo & i, int j,const VHbbEventAuxInfo & aux){
  dxy[j]=i.ipDb;
  dz[j]=i.zPVPt;
  vbtf[j]=( i.globChi2<10 && i.nPixelHits>= 1 && i.globNHits != 0 && i.nHits > 10 && i.cat & 0x1 && i.cat & 0x2 && i.nMatches >=2 && i.ipDb<.2 &&
        (i.pfChaIso+i.pfPhoIso+i.pfNeuIso)/i.p4.Pt()<.15  && fabs(i.p4.Eta())<2.4 && i.p4.Pt()>20 ) ;
float mincor=0.0;
float minrho=0.0;
float rhoN = std::max(aux.puInfo.rhoNeutral,minrho);
float eta=i.p4.Eta();
float area=0.5;
if(fabs(eta)>0.0 && fabs(eta) <= 1.0) {area=0.674;}
if(fabs(eta)>1.0 && fabs(eta) <= 1.5) {area=0.565;}
if(fabs(eta)>1.5 && fabs(eta) <= 2.0) {area=0.442;}
if(fabs(eta)>2.0 && fabs(eta) <= 2.2) {area=0.515;}
if(fabs(eta)>2.2 && fabs(eta) <= 2.3) {area=0.821;}
if(fabs(eta)>2.3 && fabs(eta) <= 2.4) {area=0.660;}
pfCorrIso[j] = (i.pfChaIso+ std::max(i.pfPhoIso+i.pfNeuIso-rhoN*area,mincor))/i.p4.Pt();
id2012tight[j]= i.isPF && i. globChi2<10 && i.nPixelHits>= 1 && i.globNHits != 0 && i.nValidLayers > 5 &&         (i.cat & 0x2) && i.nMatches >=2 && i.ipDb<.2;
}
*/

// Copied from electron selection and setSpecific
bool ElectronWP(const VHbbEvent::ElectronInfo &i,float rho,int wp) {
  float mincor=0.0;
  float minrho=0.0;
  float rhoN = std::max(rho,minrho);
  float eta=i.p4.Eta();
  float areagamma=0.5;
  float areaNH=0.5;
  float areaComb=0.5;
  if(fabs(eta) <= 1.0 ) {areagamma=0.081; areaNH=0.024; areaComb=0.10;}
  if(fabs(eta) > 1.0 &&  fabs(eta) <= 1.479 ) {areagamma=0.084; areaNH=0.037; areaComb=0.12;}
  if(fabs(eta) > 1.479 &&  fabs(eta) <= 2.0 ) {areagamma=0.048; areaNH=0.037; areaComb=0.085;}
  if(fabs(eta) > 2.0 &&  fabs(eta) <= 2.2 ) {areagamma=0.089; areaNH=0.023; areaComb=0.11;}
  if(fabs(eta) > 2.2 &&  fabs(eta) <= 2.3 ) {areagamma=0.092; areaNH=0.023; areaComb=0.12;}
  if(fabs(eta) > 2.3 &&  fabs(eta) <= 2.4 ) {areagamma=0.097; areaNH=0.021; areaComb=0.12;}
  if(fabs(eta) > 2.4  ) {areagamma=0.11; areaNH=0.021; areaComb=0.13;}


  float pfCorrIso = (i.pfChaIso+ std::max(i.pfPhoIso-rhoN*areagamma,mincor )+std::max(i.pfNeuIso-rhoN*areaNH,mincor))/i.p4.Pt();
  float iso=pfCorrIso;
  float id=i.mvaOutTrig;
  bool wp70=((fabs(eta) < 0.8 && id>0.977 && iso < 0.093) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.956 && iso < 0.095) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.966 && iso < 0.171));
  bool wp80=((fabs(eta) < 0.8 && id>0.913 && iso < 0.105) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.964 && iso < 0.178) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.899 && iso < 0.150));
  bool wp85=((fabs(eta) < 0.8 && id>0.929 && iso < 0.135) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.931 && iso < 0.159) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.805 && iso < 0.155));
  bool wp90=((fabs(eta) < 0.8 && id>0.877 && iso < 0.177) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.794 && iso < 0.180) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.846 && iso < 0.244));
  bool wp95=((fabs(eta) < 0.8 && id>0.858 && iso < 0.253) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.425 && iso < 0.225) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.759 && iso < 0.308));
  bool wpHWW=((fabs(eta) < 0.8 && id>0.94 && iso < 0.15) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.85 && iso < 0.15) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.92 && iso < 0.15));

  if (wp == 70) return wp70;
  if (wp == 80) return wp80;
  if (wp == 85) return wp85;
  if (wp == 90) return wp90;
  if (wp == 95) return wp95;

  return 0;
}

// Copied from muon selection and LeptonInfo::setSpecific
bool muonId2012Tight(const VHbbEvent::MuonInfo & i, float rho, bool requireiso = true) {
  float mincor=0.0;
  float minrho=0.0;
  float rhoN = std::max(rho,minrho);
  float eta=i.p4.Eta();
  float area=0.5;
  if(fabs(eta)>0.0 && fabs(eta) <= 1.0) {area=0.674;}
  if(fabs(eta)>1.0 && fabs(eta) <= 1.5) {area=0.565;}
  if(fabs(eta)>1.5 && fabs(eta) <= 2.0) {area=0.442;}
  if(fabs(eta)>2.0 && fabs(eta) <= 2.2) {area=0.515;}
  if(fabs(eta)>2.2 && fabs(eta) <= 2.3) {area=0.821;}
  if(fabs(eta)>2.3 && fabs(eta) <= 2.4) {area=0.660;}
  float pfCorrIso = (i.pfChaIso+ std::max(i.pfPhoIso+i.pfNeuIso-rhoN*area,mincor))/i.p4.Pt();
  return (i.isPF && i. globChi2<10 && i.nPixelHits>= 1 && i.globNHits != 0 && i.nValidLayers > 5 &&         (i.cat & 0x2) && i.nMatches >=2 && i.ipDb<.2 
	  && fabs(eta) < 2.4 && i.p4.Pt() > 20. &&
	  (i.cat & 0x1) && (!requireiso || pfCorrIso < 0.12)); // Added by SCZ from selection
}

typedef struct
{
  float et;
  float sumet;
  float sig;
  float phi;
} METInfo;

int main(int argc, char* argv[]) 
{
  gROOT->Reset();

  //  if (verbose_) cout << "Start of main" << endl;

  float rhoN;
  int nPVs;
  METInfo MET;

 float  PUweight, PUweight2011B;
  float PU0,PUp1,PUm1;

  // load framework libraries
  gSystem->Load("libFWCoreFWLite");
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
  bool verbose_ = ana.getParameter<bool>("verbose");

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
  int runMin_( in.getParameter<int>("runMin") );
  int runMax_( in.getParameter<int>("runMax") );
  unsigned int outputEvery_( in.getParameter<unsigned int>("outputEvery") );
  std::string outputFile_( out.getParameter<std::string>("fileName" ) );

  std::vector<std::string> inputFiles_( in.getParameter<std::vector<std::string> >("fileNames") );

  std::string PUmcfileName_ = in.getParameter<std::string> ("PUmcfileName") ;
  std::string PUmcfileName2011B_ = in.getParameter<std::string> ("PUmcfileName2011B") ;

  std::string PUdatafileName_ = in.getParameter<std::string> ("PUdatafileName") ;
  std::string PUdatafileName2011B_ = in.getParameter<std::string> ("PUdatafileName2011B") ;
  std::string Weight3DfileName_ = in.getParameter<std::string> ("Weight3DfileName") ;

  bool isMC_( ana.getParameter<bool>("isMC") );  

  float lumiWeight;

  edm::LumiReWeighting   lumiWeights;
  edm::Lumi3DReWeighting   lumiWeights2011B;
  if(isMC_)
    {
      lumiWeights = edm::LumiReWeighting(PUmcfileName_,PUdatafileName_ , "pileup", "pileup");

      lumiWeights2011B = edm::Lumi3DReWeighting(PUmcfileName2011B_,PUdatafileName2011B_ , "pileup", "pileup");
      if(Weight3DfileName_!="")
	{ lumiWeights2011B.weight3D_init(Weight3DfileName_.c_str()); }
      else
	{
	  lumiWeights2011B.weight3D_init(1.0); // generate the weights the fisrt time;
	}

    }


  TFile *_outFile	= new TFile(outputFile_.c_str(), "recreate");	
  TDirectory *muTrigDir = _outFile->mkdir("muTrigDir");
  TDirectory *muWCandDir =  _outFile->mkdir("muWCandDir");
  TDirectory *muRecoDir = _outFile->mkdir("muRecoDir");
  TDirectory * eleTrigDir = _outFile->mkdir("eleTrigDir");
  TDirectory *eleRecoDir = _outFile->mkdir("eleRecoDir");

  int eventNumber, runNumber, lb;

  float muTrigTree_tag_pt, muTrigTree_tag_eta, muTrigTree_probe_pt, muTrigTree_probe_eta, muTrigTree_mass, muTrigTree_tag_phi, muTrigTree_probe_phi;
  int muTrigTree_probe_passingIsoMu24L1, muTrigTree_probe_passingIsoMu24L2, muTrigTree_probe_passingIsoMu24L3, muTrigTree_probe_passingIsoMu24Iso;
  int muTrigTree_probe_passingMu40L1, muTrigTree_probe_passingMu40L2, muTrigTree_probe_passingMu40L3;
  int muTrigTree_probe_passingDiMuL1, muTrigTree_probe_passingDiMuL20, muTrigTree_probe_passingDiMuL210, muTrigTree_probe_passingDiMu8, muTrigTree_probe_passingDiMu17, muTrigTree_probe_passingDiMuDz;
  int muTrigTree_probe_passingDiMuTkL1, muTrigTree_probe_passingDiMuTkL2, muTrigTree_probe_passingDiMuTk17, muTrigTree_probe_passingDiMuTk8, muTrigTree_probe_passingDiMuTkDz;
  int muTrigTree_probe_passingDiMuL110, muTrigTree_probe_passingDiMuTkL110, muTrigTree_probe_passingDiMu17Dz, muTrigTree_probe_passingDiMuTk17Dz;
  int muTrigTree_probe_passingIsoMu20Iso;
  float muTrigTree_probe_WCandPt; int muTrigTree_probe_passingWCandPt; 
  int muTrigTree_event_Mu17_Mu8, muTrigTree_event_Mu17_TkMu8;
  int muTrigTree_probe_passingIsoMu24ORMu40, muTrigTree_probe_passingMu40ANDNOTIsoMu24, muTrigTree_probe_passingIsoMu20ANDNOTIsoMu24;
  int muTrigTree_probe_passingIsoMu20ANDNOTIsoMu24ANDNOTMu40, muTrigTree_probe_passingIsoMu24ORMu40ORIsoMu20;
  int muTrigTree_probe_passingDiMuTk17ANDNOTeventMu17Mu8, muTrigTree_probe_passingDiMuTk8ANDNOTeventMu17Mu8, muTrigTree_probe_passingDiMuTkDzANDNOTeventMu17Mu8;
  TTree *muTrigTree = new TTree("muTrigTree","muTrigTree"); muTrigTree->SetDirectory(muTrigDir);
  muTrigTree->Branch("tag_pt"            ,  &muTrigTree_tag_pt                 ,  "tag_pt/F");
  muTrigTree->Branch("tag_eta"            ,  &muTrigTree_tag_eta                 ,  "tag_eta/F");
  muTrigTree->Branch("tag_phi"            ,  &muTrigTree_tag_phi                 ,  "tag_phi/F");
  muTrigTree->Branch("probe_pt"            ,  &muTrigTree_probe_pt                 ,  "probe_pt/F");
  muTrigTree->Branch("probe_eta"            ,  &muTrigTree_probe_eta                 ,  "probe_eta/F");
  muTrigTree->Branch("probe_phi"            ,  &muTrigTree_probe_phi                 ,  "probe_phi/F");
  muTrigTree->Branch("probe_passingIsoMu24L1", &muTrigTree_probe_passingIsoMu24L1, "probe_passingIsoMu24L1/I");
  muTrigTree->Branch("probe_passingIsoMu24L2", &muTrigTree_probe_passingIsoMu24L2, "probe_passingIsoMu24L2/I");
  muTrigTree->Branch("probe_passingIsoMu24L3", &muTrigTree_probe_passingIsoMu24L3, "probe_passingIsoMu24L3/I");
  muTrigTree->Branch("probe_passingIsoMu24Iso", &muTrigTree_probe_passingIsoMu24Iso, "probe_passingIsoMu24Iso/I");
  muTrigTree->Branch("probe_passingMu40L1", &muTrigTree_probe_passingMu40L1, "probe_passingMu40L1/I");
  muTrigTree->Branch("probe_passingMu40L2", &muTrigTree_probe_passingMu40L2, "probe_passingMu40L2/I");
  muTrigTree->Branch("probe_passingMu40L3", &muTrigTree_probe_passingMu40L3, "probe_passingMu40L3/I");
  muTrigTree->Branch("probe_passingDiMuL1", &muTrigTree_probe_passingDiMuL1, "probe_passingDiMuL1/I");
  muTrigTree->Branch("probe_passingDiMuL20", &muTrigTree_probe_passingDiMuL20, "probe_passingDiMuL20/I");
  muTrigTree->Branch("probe_passingDiMuL210", &muTrigTree_probe_passingDiMuL210, "probe_passingDiMuL210/I");
  muTrigTree->Branch("probe_passingDiMu8", &muTrigTree_probe_passingDiMu8, "probe_passingDiMu8/I");
  muTrigTree->Branch("probe_passingDiMu17", &muTrigTree_probe_passingDiMu17, "probe_passingDiMu17/I");
  muTrigTree->Branch("probe_passingDiMuDz", &muTrigTree_probe_passingDiMuDz, "probe_passingDiMuDz/I");
  muTrigTree->Branch("probe_passingDiMuTkL1", &muTrigTree_probe_passingDiMuTkL1, "probe_passingDiMuTkL1/I");
  muTrigTree->Branch("probe_passingDiMuTkL2", &muTrigTree_probe_passingDiMuTkL2, "probe_passingDiMuTkL2/I");
  muTrigTree->Branch("probe_passingDiMuTk17", &muTrigTree_probe_passingDiMuTk17, "probe_passingDiMuTk17/I");
  muTrigTree->Branch("probe_passingDiMuTk8", &muTrigTree_probe_passingDiMuTk8, "probe_passingDiMuTk8/I");
  muTrigTree->Branch("probe_passingDiMuTkDz", &muTrigTree_probe_passingDiMuTkDz, "probe_passingDiMuTkDz/I");
  muTrigTree->Branch("probe_passingIsoMu20Iso", &muTrigTree_probe_passingIsoMu20Iso, "probe_passingIsoMu20Iso/I");
  muTrigTree->Branch("probe_passingDiMuL110", &muTrigTree_probe_passingDiMuL110, "probe_passingDiMuL110/I");
  muTrigTree->Branch("probe_passingDiMuTkL110", &muTrigTree_probe_passingDiMuTkL110, "probe_passingDiMuTkL110/I");
  muTrigTree->Branch("probe_passingDiMu17Dz", &muTrigTree_probe_passingDiMu17Dz, "probe_passingDiMu17Dz/I");
  muTrigTree->Branch("probe_passingDiMuTk17Dz", &muTrigTree_probe_passingDiMuTk17Dz, "probe_passingDiMuTk17Dz/I");



  muTrigTree->Branch("probe_passingIsoMu24ORMu40", &muTrigTree_probe_passingIsoMu24ORMu40, "probe_passingIsoMu24ORMu40/I");
  muTrigTree->Branch("probe_passingMu40ANDNOTIsoMu24", &muTrigTree_probe_passingMu40ANDNOTIsoMu24, "probe_passingMu40ANDNOTIsoMu24/I");
  muTrigTree->Branch("probe_passingIsoMu20ANDNOTIsoMu24", &muTrigTree_probe_passingIsoMu20ANDNOTIsoMu24, "probe_passingIsoMu20ANDNOTIsoMu24/I");
  muTrigTree->Branch("probe_passingIsoMu20ANDNOTIsoMu24ANDNOTMu40", &muTrigTree_probe_passingIsoMu20ANDNOTIsoMu24ANDNOTMu40, "probe_passingIsoMu20ANDNOTIsoMu24ANDNOTMu40/I");
  muTrigTree->Branch("probe_passingIsoMu24ORMu40ORIsoMu20", &muTrigTree_probe_passingIsoMu24ORMu40ORIsoMu20, "probe_passingIsoMu24ORMu40ORIsoMu20/I");

  muTrigTree->Branch("probe_WCandPt", &muTrigTree_probe_WCandPt, "probe_WCandPt/F");
  muTrigTree->Branch("probe_passingWCandPt", &muTrigTree_probe_passingWCandPt, "probe_passingWCandPt/I");

  muTrigTree->Branch("event_Mu17_Mu8", &muTrigTree_event_Mu17_Mu8, "event_Mu17_Mu8/I");
  muTrigTree->Branch("event_Mu17_TkMu8", &muTrigTree_event_Mu17_TkMu8, "event_Mu17_TkMu8/I");

  muTrigTree->Branch("probe_passingDiMuTk17ANDNOTeventMu17Mu8", &muTrigTree_probe_passingDiMuTk17ANDNOTeventMu17Mu8, "probe_passingDiMuTk17ANDNOTeventMu17Mu8/I");
 muTrigTree->Branch("probe_passingDiMuTk8ANDNOTeventMu17Mu8", &muTrigTree_probe_passingDiMuTk8ANDNOTeventMu17Mu8, "probe_passingDiMuTk8ANDNOTeventMu17Mu8/I");
 muTrigTree->Branch("probe_passingDiMuTkDzANDNOTeventMu17Mu8", &muTrigTree_probe_passingDiMuTkDzANDNOTeventMu17Mu8, "probe_passingDiMuTkDzANDNOTeventMu17Mu8/I");


  muTrigTree->Branch("mass", &muTrigTree_mass, "mass/F");
  muTrigTree->Branch("nPVs", &nPVs, "nPVs/I");
  muTrigTree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  muTrigTree->Branch("runNumber", &runNumber, "runNumber/I");
  muTrigTree->Branch("lb", &lb, "lb/I");
  muTrigTree->Branch("lumiWeight", &lumiWeight, "lumiWeight/F");


  float muWCandTree_tag_pt, muWCandTree_tag_eta, muWCandTree_probe_pt, muWCandTree_tag_phi, muWCandTree_probe_phi;
  float muWCandTree_probe_met, muWCandTree_probe_metPhi;
  int muWCandTree_probe_passingWCandPt;
  TTree *muWCandTree = new TTree("muWCandTree","muWCandTree"); muWCandTree->SetDirectory(muWCandDir);
  muWCandTree->Branch("tag_pt"            ,  &muWCandTree_tag_pt                 ,  "tag_pt/F");
  muWCandTree->Branch("tag_eta"            ,  &muWCandTree_tag_eta                 ,  "tag_eta/F");
  muWCandTree->Branch("tag_phi"            ,  &muWCandTree_tag_phi                 ,  "tag_phi/F");
  muWCandTree->Branch("probe_pt"            ,  &muWCandTree_probe_pt                 ,  "probe_pt/F");
  muWCandTree->Branch("probe_phi"            ,  &muWCandTree_probe_phi                 ,  "probe_phi/F");
  muWCandTree->Branch("probe_met"            ,  &muWCandTree_probe_met                 ,  "probe_met/F");
  muWCandTree->Branch("probe_metPhi"            ,  &muWCandTree_probe_metPhi                 ,  "probe_metPhi/F");
  muWCandTree->Branch("probe_passingWCandPt", &muWCandTree_probe_passingWCandPt, "probe_passingWCandPt/I");
  muWCandTree->Branch("nPVs", &nPVs, "nPVs/I");
  muWCandTree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  muWCandTree->Branch("runNumber", &runNumber, "runNumber/I");
  muWCandTree->Branch("lb", &lb, "lb/I");
  muWCandTree->Branch("lumiWeight", &lumiWeight, "lumiWeight/F");



  // (muons[it].cat & 0x2)

  float muRecoTree_tag_pt, muRecoTree_tag_eta, muRecoTree_probe_pt, muRecoTree_probe_eta, muRecoTree_mass, muRecoTree_tag_phi, muRecoTree_probe_phi;
  int muRecoTree_probe_tightNoIso, muRecoTree_probe_tightIso;
  TTree *muRecoTree = new TTree("muRecoTree","muRecoTree"); muRecoTree->SetDirectory(muRecoDir);
  muRecoTree->Branch("tag_pt"            ,  &muRecoTree_tag_pt                 ,  "tag_pt/F");
  muRecoTree->Branch("tag_eta"            ,  &muRecoTree_tag_eta                 ,  "tag_eta/F");
  muRecoTree->Branch("tag_phi"            ,  &muRecoTree_tag_phi                 ,  "tag_phi/F");
  muRecoTree->Branch("probe_pt"            ,  &muRecoTree_probe_pt                 ,  "probe_pt/F");
  muRecoTree->Branch("probe_eta"            ,  &muRecoTree_probe_eta                 ,  "probe_eta/F");
  muRecoTree->Branch("probe_phi"            ,  &muRecoTree_probe_phi                 ,  "probe_phi/F");
  muRecoTree->Branch("probe_tightNoIso"            ,  &muRecoTree_probe_tightNoIso                 ,  "probe_tightNoIso/I");
  muRecoTree->Branch("probe_tightIso"            ,  &muRecoTree_probe_tightIso                 ,  "probe_tightIso/I");
  muRecoTree->Branch("mass", &muRecoTree_mass, "mass/F");
  muRecoTree->Branch("nPVs", &nPVs, "nPVs/I");
  muRecoTree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  muRecoTree->Branch("runNumber", &runNumber, "runNumber/I");
  muRecoTree->Branch("lb", &lb, "lb/I");
  muRecoTree->Branch("lumiWeight", &lumiWeight, "lumiWeight/F");


  float eleTrigTree_tag_pt, eleTrigTree_tag_eta, eleTrigTree_probe_pt, eleTrigTree_probe_eta, eleTrigTree_mass, eleTrigTree_tag_phi, eleTrigTree_probe_phi;
  int eleTrigTree_probe_passingEle27L1, eleTrigTree_probe_passingEle27HLT;
  int eleTrigTree_probe_passingDiEle17, eleTrigTree_probe_passingDiEle8, eleTrigTree_probe_passingDiEleDz;
  int eleTrigTree_probe_wp70, eleTrigTree_probe_wp80, eleTrigTree_probe_wp85, eleTrigTree_probe_wp90, eleTrigTree_probe_wp95;
  TTree *eleTrigTree = new TTree("eleTrigTree","eleTrigTree"); eleTrigTree->SetDirectory(eleTrigDir);
  eleTrigTree->Branch("tag_pt"            ,  &eleTrigTree_tag_pt                 ,  "tag_pt/F");
  eleTrigTree->Branch("tag_eta"            ,  &eleTrigTree_tag_eta                 ,  "tag_eta/F");
  eleTrigTree->Branch("tag_phi"            ,  &eleTrigTree_tag_phi                 ,  "tag_phi/F");
  eleTrigTree->Branch("probe_pt"            ,  &eleTrigTree_probe_pt                 ,  "probe_pt/F");
  eleTrigTree->Branch("probe_eta"            ,  &eleTrigTree_probe_eta                 ,  "probe_eta/F");
  eleTrigTree->Branch("probe_phi"            ,  &eleTrigTree_probe_phi                 ,  "probe_phi/F");
  eleTrigTree->Branch("probe_passingEle27L1", &eleTrigTree_probe_passingEle27L1, "probe_passingEle27L1/I");
  eleTrigTree->Branch("probe_passingEle27HLT", &eleTrigTree_probe_passingEle27HLT, "probe_passingEle27HLT/I");
  eleTrigTree->Branch("probe_passingDiEle17", &eleTrigTree_probe_passingDiEle17, "probe_passingDiEle17/I");
  eleTrigTree->Branch("probe_passingDiEle8", &eleTrigTree_probe_passingDiEle8, "probe_passingDiEle8/I");
  eleTrigTree->Branch("probe_passingDiEleDz", &eleTrigTree_probe_passingDiEleDz, "probe_passingDiEleDz/I");
  eleTrigTree->Branch("probe_wp70"            ,  &eleTrigTree_probe_wp70                 ,  "probe_wp70/I");
  eleTrigTree->Branch("probe_wp80"            ,  &eleTrigTree_probe_wp80                 ,  "probe_wp80/I");
  eleTrigTree->Branch("probe_wp85"            ,  &eleTrigTree_probe_wp85                 ,  "probe_wp85/I");
  eleTrigTree->Branch("probe_wp90"            ,  &eleTrigTree_probe_wp90                 ,  "probe_wp90/I");
  eleTrigTree->Branch("probe_wp95"            ,  &eleTrigTree_probe_wp95                 ,  "probe_wp95/I");
  eleTrigTree->Branch("mass", &eleTrigTree_mass, "mass/F");
  eleTrigTree->Branch("nPVs", &nPVs, "nPVs/I");
  eleTrigTree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  eleTrigTree->Branch("runNumber", &runNumber, "runNumber/I");
  eleTrigTree->Branch("lb", &lb, "lb/I");
  eleTrigTree->Branch("lumiWeight", &lumiWeight, "lumiWeight/F");



  float eleRecoTree_tag_pt, eleRecoTree_tag_eta, eleRecoTree_probe_pt, eleRecoTree_probe_eta, eleRecoTree_mass, eleRecoTree_tag_phi, eleRecoTree_probe_phi;
  int eleRecoTree_probe_wp70, eleRecoTree_probe_wp80, eleRecoTree_probe_wp85, eleRecoTree_probe_wp90, eleRecoTree_probe_wp95;
  TTree *eleRecoTree = new TTree("eleRecoTree","eleRecoTree"); eleRecoTree->SetDirectory(eleRecoDir);
  eleRecoTree->Branch("tag_pt"            ,  &eleRecoTree_tag_pt                 ,  "tag_pt/F");
  eleRecoTree->Branch("tag_eta"            ,  &eleRecoTree_tag_eta                 ,  "tag_eta/F");
  eleRecoTree->Branch("tag_phi"            ,  &eleRecoTree_tag_phi                 ,  "tag_phi/F");
  eleRecoTree->Branch("probe_pt"            ,  &eleRecoTree_probe_pt                 ,  "probe_pt/F");
  eleRecoTree->Branch("probe_eta"            ,  &eleRecoTree_probe_eta                 ,  "probe_eta/F");
  eleRecoTree->Branch("probe_phi"            ,  &eleRecoTree_probe_phi                 ,  "probe_phi/F");
  eleRecoTree->Branch("probe_wp70"            ,  &eleRecoTree_probe_wp70                 ,  "probe_wp70/I");
  eleRecoTree->Branch("probe_wp80"            ,  &eleRecoTree_probe_wp80                 ,  "probe_wp80/I");
  eleRecoTree->Branch("probe_wp85"            ,  &eleRecoTree_probe_wp85                 ,  "probe_wp85/I");
  eleRecoTree->Branch("probe_wp90"            ,  &eleRecoTree_probe_wp90                 ,  "probe_wp90/I");
  eleRecoTree->Branch("probe_wp95"            ,  &eleRecoTree_probe_wp95                 ,  "probe_wp95/I");
  eleRecoTree->Branch("mass", &eleRecoTree_mass, "mass/F");
  eleRecoTree->Branch("nPVs", &nPVs, "nPVs/I");
  eleRecoTree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  eleRecoTree->Branch("runNumber", &runNumber, "runNumber/I");
  eleRecoTree->Branch("lb", &lb, "lb/I");
  eleRecoTree->Branch("lumiWeight", &lumiWeight, "lumiWeight/F");


  int ievt=0;
  int totalcount=0;

  TRandom3* rand = new TRandom3();


  //  TFile* inFile = new TFile(inputFile.c_str(), "read");
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile) {
    std::cout << iFile << std::endl;
    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
    if(inFile==0) { std::cout << "FAILED " << inputFiles_[iFile] << std::endl; continue; }


    // loop the events
      
      fwlite::Event ev(inFile);
      fwlite::Handle< VHbbEventAuxInfo > vhbbAuxHandle; 
      for(ev.toBegin(); !ev.atEnd() ; ++ev, ++ievt)
        {
          if (ievt <= skipEvents_) continue;
          if (maxEvents_ >= 0){
              if (ievt > maxEvents_ + skipEvents_) break;
          }
      const char * lab = "HbbAnalyzerNew";
      vhbbAuxHandle.getByLabel(ev,lab,0,0);
      const VHbbEventAuxInfo & aux = *vhbbAuxHandle.product();

      if(int(ev.id().run()) < runMin_ && runMin_ > 0) continue;
      if(int(ev.id().run()) > runMax_ && runMax_ > 0) continue;

      if (ievt%1000==0) cout << " Entry=" << ievt << " event=" << ev.id().event() << " run=" << ev.id().run() << " lb=" << ev.id().luminosityBlock() << endl;

      lumiWeight = 1.;
      if(isMC_){
 
	// PU weights // Run2011A
	std::map<int, unsigned int>::const_iterator puit = aux.puInfo.pus.find(0);
	int npu =puit->second ;
	//	PUweight =  lumiWeights.weight( npu );        
	//	pu->Fill(puit->second);
	std::map<int, unsigned int>::const_iterator puit0 =  aux.puInfo.pus.find(0);
	std::map<int, unsigned int>::const_iterator puitm1 = aux.puInfo.pus.find(-1);
	std::map<int, unsigned int>::const_iterator puitp1 = aux.puInfo.pus.find(+1);
	//	PU0=puit0->second;
	//	PUp1=puitp1->second;
	//	PUm1=puitm1->second;
	//	input3DPU->Fill(PUm1,PU0,PUp1);  
	lumiWeight = lumiWeights2011B.weight3D( puitm1->second, puit0->second,puitp1->second); 

      }


      runNumber = ev.id().run();
      eventNumber = ev.id().event();
      lb = ev.id().luminosityBlock();

      VHbbEvent modifiedEvent;;
      const VHbbEvent *  iEvent =0;

      fwlite::Handle< VHbbEvent > vhbbHandle;
      vhbbHandle.getByLabel(ev,"HbbAnalyzerNew");
      modifiedEvent = *vhbbHandle.product();

      iEvent = vhbbHandle.product();

	float trigMatchCut = 0.1;
	float L1MatchCut = 0.2;
        std::string mu24name = "HLT_IsoMu24_eta2p1_v";
	std::string mu24nameL1 = "hltL1sMu16Eta2p1";
	std::string mu24nameL2 = "hltL2fL1sMu16Eta2p1L1f0L2Filtered16Q";
	std::string mu24nameL3 = "hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered24Q";
	std::string mu24nameIso = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10";

	std::string mu40name = "HLT_Mu40_eta2p1_v";
        std::string mu40nameL1 = "hltL1sMu16Eta2p1";
        std::string mu40nameL2 = "hltL2fL1sMu16Eta2p1L1f0L2Filtered16Q";
	std::string mu40nameL3 = "hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q";

        std::string diMuName = "HLT_Mu17_Mu8_v";
        std::string diMuNameL1 = "hltL1sL1DoubleMu10MuOpen";
        std::string diMuNameL20 = "hltL2pfL1DoubleMu10MuOpenL1f0L2PreFiltered0";
        std::string diMuNameL210 = "hltL2fL1DoubleMu10MuOpenL1f0L2Filtered10";
        std::string diMuName8 = "hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8";
        std::string diMuName17 = "hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17";
        std::string diMuNameDz = "hltDiMuonMu17Mu8DzFiltered0p2";

        std::string diMuTkName = "HLT_Mu17_TkMu8_v";
        std::string diMuTkNameL1 = "hltL1sL1DoubleMu10MuOpen";
	std::string diMuTkNameL2 = "hltL2fL1sDoubleMu10MuOpenL1f0L2Filtered10";
        std::string diMuTkName17 = "hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17";
        std::string diMuTkName8 = "hltDiMuonGlbFiltered17TrkFiltered8";
        std::string diMuTkNameDz = "hltDiMuonGlb17Trk8DzFiltered0p2";

	std::string mu20name = "HLT_IsoMu20_eta2p1_WCandPt80_v";
        std::string mu20nameIso = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f20L3crIsoFiltered10";
	std::string mu20nameWCand = "hlt2IsoMu20eta2p1PFMHTPt80";

	std::string ele27name = "HLT_Ele27_WP80_v";
	std::string ele27nameL1 = "hltL1sL1SingleEG20ORL1SingleEG22";
	std::string ele27nameHLT = "hltEle27WP80TrackIsoFilter";

	std::string diEleName = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
	std::string diEleName17 = "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter";
	std::string diEleName8 = "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter";
	std::string diEleNameDz = "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ";

	std::vector<TLorentzVector> mu24L1;
	std::vector<TLorentzVector> mu24L2;
	std::vector<TLorentzVector> mu24L3;
	std::vector<TLorentzVector> mu24Iso;

        std::vector<TLorentzVector> mu40L1;
        std::vector<TLorentzVector> mu40L2;
        std::vector<TLorentzVector> mu40L3;

        std::vector<TLorentzVector> diMuL110;
	std::vector<TLorentzVector> diMuL1;
	std::vector<TLorentzVector> diMuL20;
	std::vector<TLorentzVector> diMuL210;
	std::vector<TLorentzVector> diMu8;
	std::vector<TLorentzVector> diMu17;
	std::vector<TLorentzVector> diMuDz;

        std::vector<TLorentzVector> diMuTkL110;
        std::vector<TLorentzVector> diMuTkL1;
        std::vector<TLorentzVector> diMuTkL2;
        std::vector<TLorentzVector> diMuTk17;
        std::vector<TLorentzVector> diMuTk8;
        std::vector<TLorentzVector> diMuTkDz;

        std::vector<TLorentzVector> mu20Iso;
        std::vector<TLorentzVector> mu20WCand;

	std::vector<TLorentzVector> ele27L1;
	std::vector<TLorentzVector> ele27HLT;
	
	std::vector<TLorentzVector> diEle17;
        std::vector<TLorentzVector> diEle8;
        std::vector<TLorentzVector> diEleDz;

        fwlite::Handle< TriggerEvent > triggerEvent;
        //      iEvent.getByLabel( "patTriggerEvent", triggerEvent );
        triggerEvent.getByLabel(ev,"patTriggerEvent");

        const TriggerPathRefVector trigPathRefs (triggerEvent->pathRefs() );

	for ( TriggerPathRefVector::const_iterator trigPathRef  = trigPathRefs.begin(); trigPathRef != trigPathRefs.end(); ++trigPathRef ) {
	  std::string name = (*trigPathRef)->name();
	  if (name.size() >= mu24name.size() && name.compare(0,mu24name.size(),mu24name) == 0) {
	    const TriggerFilterRefVector trigFilterRefs( triggerEvent->pathFilters(name) );
	    for ( TriggerFilterRefVector::const_iterator trigFilterRef  = trigFilterRefs.begin(); trigFilterRef != trigFilterRefs.end(); ++trigFilterRef ) {
	      std::string label = (*trigFilterRef)->label();
	      if (label == mu24nameL1) {
		if (verbose_) cout << label << endl;
		const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
		for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
		  mu24L1.push_back(GENPTOLORP(*trigObjRef));
		}
	      }
	      if (label == mu24nameL2) {
		if (verbose_) cout << label << endl;
		const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
		for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
		  mu24L2.push_back(GENPTOLORP(*trigObjRef));
		}
	      }
	      if (label == mu24nameL3) {
		if (verbose_) cout << label << endl;
		const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
		for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
		  mu24L3.push_back(GENPTOLORP(*trigObjRef));
		}
	      }
	      if (label == mu24nameIso) {
		if (verbose_) cout << label << endl;
		const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
		for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
		  mu24Iso.push_back(GENPTOLORP(*trigObjRef));
		}
	      }
	    }
	  }
          if (name.size() >= mu40name.size() && name.compare(0,mu40name.size(),mu40name) == 0) {
            const TriggerFilterRefVector trigFilterRefs( triggerEvent->pathFilters(name) );
            for ( TriggerFilterRefVector::const_iterator trigFilterRef  = trigFilterRefs.begin(); trigFilterRef != trigFilterRefs.end(); ++trigFilterRef ) {
	      std::string label = (*trigFilterRef)->label();
              if (label == mu40nameL1) {
                if (verbose_) cout << label << endl;
		const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
		for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
		  mu40L1.push_back(GENPTOLORP(*trigObjRef));
                }
              }
              if (label == mu40nameL2) {
		if (verbose_) cout << label << endl;
                const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
                for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  mu40L2.push_back(GENPTOLORP(*trigObjRef));
                }
              }
              if (label == mu40nameL3) {
		if (verbose_) cout << label << endl;
                const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
		for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  mu40L3.push_back(GENPTOLORP(*trigObjRef));
                }
              }
            }
	  }
          if (name.size() >= diMuName.size() && name.compare(0,diMuName.size(),diMuName) == 0) {
            muTrigTree_event_Mu17_Mu8 =(*trigPathRef)->wasAccept();
            const TriggerFilterRefVector trigFilterRefs( triggerEvent->pathFilters(name) );
            for ( TriggerFilterRefVector::const_iterator trigFilterRef  = trigFilterRefs.begin(); trigFilterRef != trigFilterRefs.end(); ++trigFilterRef ) {
	      std::string label = (*trigFilterRef)->label();
              if (label == diMuNameL1) {
		if (verbose_) cout << label << endl;
		const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
                for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  diMuL1.push_back(GENPTOLORP(*trigObjRef));
		  if ((*trigObjRef)->pt() >= 10.) diMuL110.push_back(GENPTOLORP(*trigObjRef));
		}
              }
	      if (label == diMuNameL20) {
		if (verbose_) cout << label << endl;
                const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
		for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
		  diMuL20.push_back(GENPTOLORP(*trigObjRef));
                }
	      }
              if (label == diMuNameL210) {
                if (verbose_) cout << label << endl;
                const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
                for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  diMuL210.push_back(GENPTOLORP(*trigObjRef));
                }
              }
              if (label == diMuName8) {
                if (verbose_) cout << label << endl;
		const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
                for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  diMu8.push_back(GENPTOLORP(*trigObjRef));
                }
              }
              if (label == diMuName17) {
                if (verbose_) cout << label << endl;
		const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
                for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  diMu17.push_back(GENPTOLORP(*trigObjRef));
		}
	      }
              if (label == diMuNameDz) {
                if (verbose_) cout << label << endl;
		const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
                for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  diMuDz.push_back(GENPTOLORP(*trigObjRef));
		}
	      }
            }
          }
          if (name.size() >= diMuTkName.size() && name.compare(0,diMuTkName.size(),diMuTkName) == 0) {
            muTrigTree_event_Mu17_TkMu8 =(*trigPathRef)->wasAccept();
            const TriggerFilterRefVector trigFilterRefs( triggerEvent->pathFilters(name) );
            for ( TriggerFilterRefVector::const_iterator trigFilterRef  = trigFilterRefs.begin(); trigFilterRef != trigFilterRefs.end(); ++trigFilterRef ) {
	      std::string label = (*trigFilterRef)->label();
              if (label == diMuTkNameL1) {
                if (verbose_) cout << label << endl;
                const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
                for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  diMuTkL1.push_back(GENPTOLORP(*trigObjRef));
                  if ((*trigObjRef)->pt() >= 10.) diMuTkL110.push_back(GENPTOLORP(*trigObjRef));
                }
              }
              if (label == diMuTkNameL2) {
		if (verbose_) cout << label << endl;
                const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
		for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  diMuTkL2.push_back(GENPTOLORP(*trigObjRef));
                }
              }
              if (label == diMuTkName8) {
                if (verbose_) cout << label << endl;
                const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
                for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  diMuTk8.push_back(GENPTOLORP(*trigObjRef));
                }
              }
              if (label == diMuTkName17) {
                if (verbose_) cout << label << endl;
                const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
                for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  diMuTk17.push_back(GENPTOLORP(*trigObjRef));
                }
              }
              if (label == diMuTkNameDz) {
                if (verbose_) cout << label << endl;
                const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
                for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  diMuTkDz.push_back(GENPTOLORP(*trigObjRef));
                }
              }
            }
          }
          if (name.size() >= mu20name.size() && name.compare(0,mu20name.size(),mu20name) == 0) {
            const TriggerFilterRefVector trigFilterRefs( triggerEvent->pathFilters(name) );
            for ( TriggerFilterRefVector::const_iterator trigFilterRef  = trigFilterRefs.begin(); trigFilterRef != trigFilterRefs.end(); ++trigFilterRef ) {
	      std::string label = (*trigFilterRef)->label();
              if (label == mu20nameIso) {
		if (verbose_) cout << label << endl;
                const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
                for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  mu20Iso.push_back(GENPTOLORP(*trigObjRef));
		}
              }
              if (label == mu20nameWCand) {
                if (verbose_) cout << label << endl;
                const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
                for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  mu20WCand.push_back(GENPTOLORP(*trigObjRef));
                }
              }
            }
          }
          if (name.size() >= ele27name.size() && name.compare(0,ele27name.size(),ele27name) == 0) {
            const TriggerFilterRefVector trigFilterRefs( triggerEvent->pathFilters(name) );
            for ( TriggerFilterRefVector::const_iterator trigFilterRef  = trigFilterRefs.begin(); trigFilterRef != trigFilterRefs.end(); ++trigFilterRef ) {
	      std::string label = (*trigFilterRef)->label();
              if (label == ele27nameL1) {
                if (verbose_) cout << label << endl;
                const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
                for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  ele27L1.push_back(GENPTOLORP(*trigObjRef));
                }
              }
              if (label == ele27nameHLT) {
		if (verbose_) cout << label << endl;
                const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
		for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
                  ele27HLT.push_back(GENPTOLORP(*trigObjRef));
                }
              }
	    }
	  }
          if (name.size() >= diEleName.size() && name.compare(0,diEleName.size(),diEleName) == 0) {
            const TriggerFilterRefVector trigFilterRefs( triggerEvent->pathFilters(name) );
            for ( TriggerFilterRefVector::const_iterator trigFilterRef  = trigFilterRefs.begin(); trigFilterRef != trigFilterRefs.end(); ++trigFilterRef ) {
	      std::string label = (*trigFilterRef)->label();
	      if (label == diEleName8) {
		if (verbose_) cout << label << endl;
		const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
		for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
		  diEle8.push_back(GENPTOLORP(*trigObjRef));
		}
	      }
	      if (label == diEleName17) {
		if (verbose_) cout << label << endl;
		const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
		for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
		  diEle17.push_back(GENPTOLORP(*trigObjRef));
		}
	      }
	      if (label == diEleNameDz) {
		if (verbose_) cout << label << endl;
		const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
		for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  if (verbose_) std::cout << "  pushing back " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
		  diEleDz.push_back(GENPTOLORP(*trigObjRef));
		}
	      }
	    }
	  }
	}

        nPVs=aux.pvInfo.nVertices;
	rhoN = aux.puInfo.rhoNeutral;

	if (verbose_) cout << iEvent->muInfo.size() << " " << iEvent->eleInfo.size() << endl;

	std::vector<unsigned int> potentialMuTags;

	for(size_t m=0;m<iEvent->muInfo.size();m++) {
	  if (muonId2012Tight(iEvent->muInfo[m],rhoN)) {
	    if (verbose_) std::cout << "FOUND identified Mu " << iEvent->muInfo[m].p4.Pt() << " " << iEvent->muInfo[m].p4.Eta() << " " << iEvent->muInfo[m].p4.Phi() << endl;
	    potentialMuTags.push_back(m);
	  }
	}

	if (verbose_) cout << " potentialMuTags.size()=" << potentialMuTags.size() << endl;

	if (potentialMuTags.size()>0) {
	  unsigned int tagNumber = potentialMuTags[rand->Integer(potentialMuTags.size())];
	  if (verbose_) std::cout << " Picked Tag Mu " << iEvent->muInfo[tagNumber].p4.Pt() << " " << iEvent->muInfo[tagNumber].p4.Eta() << " " << iEvent->muInfo[tagNumber].p4.Phi() << endl;
	  muTrigTree_tag_pt = iEvent->muInfo[tagNumber].p4.Pt();
          muTrigTree_tag_eta = iEvent->muInfo[tagNumber].p4.Eta();
          muTrigTree_tag_phi = iEvent->muInfo[tagNumber].p4.Phi();
          muRecoTree_tag_pt = iEvent->muInfo[tagNumber].p4.Pt();
          muRecoTree_tag_eta = iEvent->muInfo[tagNumber].p4.Eta();
          muRecoTree_tag_phi = iEvent->muInfo[tagNumber].p4.Phi();
          muWCandTree_tag_pt = iEvent->muInfo[tagNumber].p4.Pt();
          muWCandTree_tag_eta = iEvent->muInfo[tagNumber].p4.Eta();
          muWCandTree_tag_phi = iEvent->muInfo[tagNumber].p4.Phi();

	  bool matched = false;
	  for (size_t i=0 ; i < mu24Iso.size(); i++) {
	    if (iEvent->muInfo[tagNumber].p4.DeltaR(mu24Iso[i]) < trigMatchCut) {
	      matched = true;
	      if (verbose_) cout << "   Matches trigger " << mu24Iso[i].Pt() << " " << mu24Iso[i].Eta() << " " << mu24Iso[i].Phi() << endl;
	    }
	  }
	  if (!matched) if (verbose_) cout << "   Has no match!" << endl;
	  if (matched) {
	    std::vector<unsigned int> potentialHLTProbes;
	    for(size_t m=0;m<iEvent->muInfo.size();m++) {
	      if (m == tagNumber) continue;
	      float mass = (iEvent->muInfo[m].p4 + iEvent->muInfo[tagNumber].p4).M();
	      if (muonId2012Tight(iEvent->muInfo[m],rhoN) && (iEvent->muInfo[m].charge + iEvent->muInfo[tagNumber].charge) == 0 && (mass > 50. && mass < 200.)) {
		potentialHLTProbes.push_back(m);
	      }
	    }
	    if (potentialHLTProbes.size()>0) {
	      unsigned int probeNumber = potentialHLTProbes[rand->Integer(potentialHLTProbes.size())];
	      if (verbose_) std::cout << "   Picked Probe Mu " << iEvent->muInfo[probeNumber].p4.Pt() << " " << iEvent->muInfo[probeNumber].p4.Eta() << " " << iEvent->muInfo[probeNumber].p4.Phi() << endl;
	      muTrigTree_probe_pt= iEvent->muInfo[probeNumber].p4.Pt();
	      muTrigTree_probe_eta = iEvent->muInfo[probeNumber].p4.Eta();
	      muTrigTree_probe_phi = iEvent->muInfo[probeNumber].p4.Phi();
	      muTrigTree_probe_passingIsoMu24L1 = 0;
	      muTrigTree_probe_passingIsoMu24L2 = 0;
	      muTrigTree_probe_passingIsoMu24L3 = 0;
	      muTrigTree_probe_passingIsoMu24Iso = 0;
	      muTrigTree_probe_passingMu40L1 = 0;
	      muTrigTree_probe_passingMu40L2 = 0;
	      muTrigTree_probe_passingMu40L3 = 0;
	      muTrigTree_probe_passingDiMuL1 = 0;
	      muTrigTree_probe_passingDiMuL20 = 0;
	      muTrigTree_probe_passingDiMuL210 = 0;
	      muTrigTree_probe_passingDiMu8 = 0;
	      muTrigTree_probe_passingDiMu17 = 0;
	      muTrigTree_probe_passingDiMuDz = 0;
              muTrigTree_probe_passingDiMuTkL1 = 0;
              muTrigTree_probe_passingDiMuTkL2 = 0;
              muTrigTree_probe_passingDiMuTk8 = 0;
              muTrigTree_probe_passingDiMuTk17 = 0;
              muTrigTree_probe_passingDiMuTkDz = 0;
              muTrigTree_probe_passingIsoMu20Iso = 0;
	      muTrigTree_probe_passingWCandPt = 0;
	      muTrigTree_probe_passingDiMuTkL110 = 0;
              muTrigTree_probe_passingDiMuL110 = 0;
	      muTrigTree_mass = (iEvent->muInfo[probeNumber].p4 + iEvent->muInfo[tagNumber].p4).M();

              for (size_t i=0 ; i < mu24L1.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(mu24L1[i]) < L1MatchCut) muTrigTree_probe_passingIsoMu24L1 = 1;
              }
              for (size_t i=0 ; i < mu24L2.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(mu24L2[i]) < trigMatchCut) muTrigTree_probe_passingIsoMu24L2 = 1;
              }
              for (size_t i=0 ; i < mu24L3.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(mu24L3[i]) < trigMatchCut) muTrigTree_probe_passingIsoMu24L3 = 1;
              }
	      for (size_t i=0 ; i < mu24Iso.size(); i++) {
		if (iEvent->muInfo[probeNumber].p4.DeltaR(mu24Iso[i]) < trigMatchCut) muTrigTree_probe_passingIsoMu24Iso = 1;
	      }
              for (size_t i=0 ; i < mu40L1.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(mu40L1[i]) < L1MatchCut) muTrigTree_probe_passingMu40L1 = 1;
              }
              for (size_t i=0 ; i < mu40L2.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(mu40L2[i]) < trigMatchCut) muTrigTree_probe_passingMu40L2 = 1;
              }
              for (size_t i=0 ; i < mu40L3.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(mu40L3[i]) < trigMatchCut) muTrigTree_probe_passingMu40L3 = 1;
              }
	      for (size_t i=0 ; i < diMuL1.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuL1[i]) < L1MatchCut) muTrigTree_probe_passingDiMuL1 = 1;
              }
              for (size_t i=0 ; i < diMuL110.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuL110[i]) < L1MatchCut) muTrigTree_probe_passingDiMuL110 = 1;
              }
              for (size_t i=0 ; i < diMuL20.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuL20[i]) < trigMatchCut) muTrigTree_probe_passingDiMuL20 = 1;
	      }
              for (size_t i=0 ; i < diMuL210.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuL210[i]) < trigMatchCut) muTrigTree_probe_passingDiMuL210 = 1;
	      }
              for (size_t i=0 ; i < diMu8.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMu8[i]) < trigMatchCut) muTrigTree_probe_passingDiMu8 = 1;
	      }
              for (size_t i=0 ; i < diMu17.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMu17[i]) < trigMatchCut) muTrigTree_probe_passingDiMu17 = 1;
	      }
              for (size_t i=0 ; i < diMuDz.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuDz[i]) < trigMatchCut) muTrigTree_probe_passingDiMuDz = 1;
	      }
              for (size_t i=0 ; i < diMuTkL1.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuTkL1[i]) < L1MatchCut) muTrigTree_probe_passingDiMuTkL1 = 1;
              }
              for (size_t i=0 ; i < diMuTkL110.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuTkL110[i]) < L1MatchCut) muTrigTree_probe_passingDiMuTkL110 = 1;
              }
              for (size_t i=0 ; i < diMuTkL2.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuTkL2[i]) < trigMatchCut) muTrigTree_probe_passingDiMuTkL2 = 1;
              }
              for (size_t i=0 ; i < diMuTk8.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuTk8[i]) < trigMatchCut) muTrigTree_probe_passingDiMuTk8 = 1;
              }
              for (size_t i=0 ; i < diMuTk17.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuTk17[i]) < trigMatchCut) muTrigTree_probe_passingDiMuTk17 = 1;
              }
	      for (size_t i=0 ; i < diMuTkDz.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuTkDz[i]) < trigMatchCut) muTrigTree_probe_passingDiMuTkDz = 1;
              }
              for (size_t i=0 ; i < mu20Iso.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(mu20Iso[i]) < trigMatchCut) muTrigTree_probe_passingIsoMu20Iso = 1;
              }

	      muTrigTree_probe_passingIsoMu24ORMu40 = (int)(muTrigTree_probe_passingIsoMu24Iso||muTrigTree_probe_passingMu40L3);
              muTrigTree_probe_passingMu40ANDNOTIsoMu24 = (int)(muTrigTree_probe_passingMu40L3&&!muTrigTree_probe_passingIsoMu24Iso);
              muTrigTree_probe_passingIsoMu20ANDNOTIsoMu24 = (int)(muTrigTree_probe_passingIsoMu20Iso&&!muTrigTree_probe_passingIsoMu24Iso);
              muTrigTree_probe_passingIsoMu20ANDNOTIsoMu24ANDNOTMu40 = (int)(muTrigTree_probe_passingIsoMu20Iso&&(!muTrigTree_probe_passingIsoMu24Iso)&&(!muTrigTree_probe_passingMu40L3));
              muTrigTree_probe_passingIsoMu24ORMu40ORIsoMu20 = (int)(muTrigTree_probe_passingIsoMu24Iso||muTrigTree_probe_passingMu40L3||muTrigTree_probe_passingIsoMu20Iso);
	      muTrigTree_probe_passingDiMuTk17ANDNOTeventMu17Mu8 = (int)(muTrigTree_probe_passingDiMuTk17&&!muTrigTree_event_Mu17_Mu8);
              muTrigTree_probe_passingDiMuTk8ANDNOTeventMu17Mu8 = (int)(muTrigTree_probe_passingDiMuTk8&&!muTrigTree_event_Mu17_Mu8);
              muTrigTree_probe_passingDiMuTkDzANDNOTeventMu17Mu8 = (int)(muTrigTree_probe_passingDiMuTkDz&&!muTrigTree_event_Mu17_Mu8);

	      muTrigTree_probe_passingDiMu17Dz = muTrigTree_probe_passingDiMu17 && muTrigTree_probe_passingDiMuDz;
	      muTrigTree_probe_passingDiMuTk17Dz = muTrigTree_probe_passingDiMuTk17&& muTrigTree_probe_passingDiMuTkDz;

              TLorentzVector WCand = iEvent->pfmet.p4 + iEvent->muInfo[probeNumber].p4;

              muTrigTree_probe_WCandPt= WCand.Pt();
              for (size_t i=0 ; i < mu20WCand.size(); i++) {
		if (iEvent->muInfo[probeNumber].p4.DeltaR(mu20WCand[i]) < trigMatchCut) muTrigTree_probe_passingWCandPt = 1;
              }

	      muTrigTree->Fill(); // Only if trigger-matched tag and a probe

	    } // potential HLT probe size > 0

	    if (iEvent->pfmet.p4.Et() > 0.) {

	      TLorentzVector WCand = iEvent->pfmet.p4 + iEvent->muInfo[tagNumber].p4;

              muWCandTree_probe_pt= WCand.Pt();
              muWCandTree_probe_phi = WCand.Phi();
	      muWCandTree_probe_met = iEvent->pfmet.p4.Et();
	      muWCandTree_probe_metPhi = iEvent->pfmet.p4.Phi();
	      muWCandTree_probe_passingWCandPt = (mu20WCand.size() > 0);

	      muWCandTree->Fill();

	    }

	    std::vector<unsigned int> potentialRecoProbes;
            for(size_t m=0;m<iEvent->muInfo.size();m++) {
              if (m == tagNumber) continue;
              float mass = (iEvent->muInfo[m].p4 + iEvent->muInfo[tagNumber].p4).M();
              if ((iEvent->muInfo[m].cat & 0x2) && iEvent->muInfo[m].p4.Pt() > 20. && fabs(iEvent->muInfo[m].p4.Eta()) < 2.4
		  && (iEvent->muInfo[m].charge + iEvent->muInfo[tagNumber].charge) == 0 && (mass > 50. && mass < 200.)) {
		potentialRecoProbes.push_back(m);
              }
            }
            if (potentialRecoProbes.size()>0) {
              unsigned int probeNumber = potentialRecoProbes[rand->Integer(potentialRecoProbes.size())];
              muRecoTree_probe_pt= iEvent->muInfo[probeNumber].p4.Pt();
              muRecoTree_probe_eta = iEvent->muInfo[probeNumber].p4.Eta();
              muRecoTree_probe_phi = iEvent->muInfo[probeNumber].p4.Phi();
	      muRecoTree_probe_tightNoIso = muonId2012Tight(iEvent->muInfo[probeNumber],rhoN,false);
              muRecoTree_probe_tightIso = muonId2012Tight(iEvent->muInfo[probeNumber],rhoN);
              muRecoTree_mass = (iEvent->muInfo[probeNumber].p4 + iEvent->muInfo[tagNumber].p4).M();
	      
	      muRecoTree->Fill(); // Only if trigger-matched tag and a probe

	    } // potential Reco probe size > 0

	  } // trigger-matched tag
	} // potential tag size > 0

        std::vector<unsigned int> potentialEleTags;

        for(size_t m=0;m<iEvent->eleInfo.size();m++) {
          if (ElectronWP(iEvent->eleInfo[m],rhoN,80) && iEvent->eleInfo[m].p4.Pt() > 20.) {
	    if (verbose_) std::cout << "FOUND identified Ele " << iEvent->eleInfo[m].p4.Pt() << " " << iEvent->eleInfo[m].p4.Eta() << " " << iEvent->eleInfo[m].p4.Phi() << endl;
            potentialEleTags.push_back(m);
          }
        }

        if (verbose_) cout << " potentialEleTags.size()=" << potentialEleTags.size() << endl;

        if (potentialEleTags.size()>0) {
          unsigned int tagNumber = potentialEleTags[rand->Integer(potentialEleTags.size())];
	  if (verbose_) std::cout << " Picked Tag Ele " << iEvent->eleInfo[tagNumber].p4.Pt() << " " << iEvent->eleInfo[tagNumber].p4.Eta() << " " << iEvent->eleInfo[tagNumber].p4.Phi() << endl;
          eleTrigTree_tag_pt = iEvent->eleInfo[tagNumber].p4.Pt();
          eleTrigTree_tag_eta = iEvent->eleInfo[tagNumber].p4.Eta();
          eleTrigTree_tag_phi = iEvent->eleInfo[tagNumber].p4.Phi();
          eleRecoTree_tag_pt = iEvent->eleInfo[tagNumber].p4.Pt();
          eleRecoTree_tag_eta = iEvent->eleInfo[tagNumber].p4.Eta();
          eleRecoTree_tag_phi = iEvent->eleInfo[tagNumber].p4.Phi();

          bool matched = false;
          for (size_t i=0 ; i < ele27HLT.size(); i++) {
            if (iEvent->eleInfo[tagNumber].p4.DeltaR(ele27HLT[i]) < trigMatchCut) {
              matched = true;
              if (verbose_) cout << "   Matches trigger " << ele27HLT[i].Pt() << " " << ele27HLT[i].Eta() << " " << ele27HLT[i].Phi() << endl;
            }
          }
          if (!matched) if (verbose_) cout << "   Has no match!" << endl;
          if (matched) {
	    std::vector<unsigned int> potentialHLTProbes;
            for(size_t m=0;m<iEvent->eleInfo.size();m++) {
              if (m == tagNumber) continue;
              float mass = (iEvent->eleInfo[m].p4 + iEvent->eleInfo[tagNumber].p4).M();
              if (ElectronWP(iEvent->eleInfo[m],rhoN,95) && (iEvent->eleInfo[m].charge + iEvent->eleInfo[tagNumber].charge) == 0 && (mass > 50. && mass < 200.)) {
                potentialHLTProbes.push_back(m);
              }
            }
            if (potentialHLTProbes.size()>0) {
              unsigned int probeNumber = potentialHLTProbes[rand->Integer(potentialHLTProbes.size())];
	      if (verbose_) std::cout << "   Picked Probe Ele " << iEvent->eleInfo[probeNumber].p4.Pt() << " " << iEvent->eleInfo[probeNumber].p4.Eta() << " " << iEvent->eleInfo[probeNumber].p4.Phi() << endl;
              eleTrigTree_probe_pt= iEvent->eleInfo[probeNumber].p4.Pt();
              eleTrigTree_probe_eta = iEvent->eleInfo[probeNumber].p4.Eta();
              eleTrigTree_probe_phi = iEvent->eleInfo[probeNumber].p4.Phi();
              eleTrigTree_probe_wp70 = ElectronWP(iEvent->eleInfo[probeNumber],rhoN,70);
              eleTrigTree_probe_wp80 = ElectronWP(iEvent->eleInfo[probeNumber],rhoN,80);
              eleTrigTree_probe_wp85 = ElectronWP(iEvent->eleInfo[probeNumber],rhoN,85);
              eleTrigTree_probe_wp90 = ElectronWP(iEvent->eleInfo[probeNumber],rhoN,90);
              eleTrigTree_probe_wp95 = ElectronWP(iEvent->eleInfo[probeNumber],rhoN,95);
	      eleTrigTree_probe_passingEle27L1 = 0;
	      eleTrigTree_probe_passingEle27HLT = 0;
	      eleTrigTree_probe_passingDiEle17 = 0;
	      eleTrigTree_probe_passingDiEle8 = 0;
	      eleTrigTree_probe_passingDiEleDz = 0;
	      eleTrigTree_mass = (iEvent->eleInfo[probeNumber].p4 + iEvent->eleInfo[tagNumber].p4).M();
	      
	      for (size_t i=0 ; i < ele27L1.size(); i++) {
		if (iEvent->eleInfo[probeNumber].p4.DeltaR(ele27L1[i]) < L1MatchCut) eleTrigTree_probe_passingEle27L1 = 1;
	      }
              for (size_t i=0 ; i < ele27HLT.size(); i++) {
                if (iEvent->eleInfo[probeNumber].p4.DeltaR(ele27HLT[i]) < trigMatchCut) eleTrigTree_probe_passingEle27HLT = 1;
              }
              for (size_t i=0 ; i < diEle17.size(); i++) {
                if (iEvent->eleInfo[probeNumber].p4.DeltaR(diEle17[i]) < trigMatchCut) eleTrigTree_probe_passingDiEle17 = 1;
              }
              for (size_t i=0 ; i < diEle8.size(); i++) {
                if (iEvent->eleInfo[probeNumber].p4.DeltaR(diEle8[i]) < trigMatchCut) eleTrigTree_probe_passingDiEle8 = 1;
              }
              for (size_t i=0 ; i < diEleDz.size(); i++) {
                if (iEvent->eleInfo[probeNumber].p4.DeltaR(diEleDz[i]) < trigMatchCut) eleTrigTree_probe_passingDiEleDz = 1;
              }

	      eleTrigTree->Fill();

	    } // HLT probes > 0

	    std::vector<unsigned int> potentialRecoProbes;
            for(size_t m=0;m<iEvent->eleInfo.size();m++) {
              if (m == tagNumber) continue;
              float mass = (iEvent->eleInfo[m].p4 + iEvent->eleInfo[tagNumber].p4).M();
              if ((iEvent->eleInfo[m].p4.Pt() > 20. && fabs(iEvent->eleInfo[m].p4.Eta()) < 2.5)
                  && (iEvent->eleInfo[m].charge + iEvent->eleInfo[tagNumber].charge) == 0 && (mass > 50. && mass < 200.)) {
                potentialRecoProbes.push_back(m);
              }
            }
            if (potentialRecoProbes.size()>0) {
              unsigned int probeNumber = potentialRecoProbes[rand->Integer(potentialRecoProbes.size())];
              eleRecoTree_probe_pt= iEvent->eleInfo[probeNumber].p4.Pt();
              eleRecoTree_probe_eta = iEvent->eleInfo[probeNumber].p4.Eta();
              eleRecoTree_probe_phi = iEvent->eleInfo[probeNumber].p4.Phi();
	      eleRecoTree_probe_wp70 = ElectronWP(iEvent->eleInfo[probeNumber],rhoN,70);
	      eleRecoTree_probe_wp80 = ElectronWP(iEvent->eleInfo[probeNumber],rhoN,80);
	      eleRecoTree_probe_wp85 = ElectronWP(iEvent->eleInfo[probeNumber],rhoN,85);
	      eleRecoTree_probe_wp90 = ElectronWP(iEvent->eleInfo[probeNumber],rhoN,90);
	      eleRecoTree_probe_wp95 = ElectronWP(iEvent->eleInfo[probeNumber],rhoN,95);
              eleRecoTree_mass = (iEvent->eleInfo[probeNumber].p4 + iEvent->eleInfo[tagNumber].p4).M();

              eleRecoTree->Fill(); // Only if trigger-matched tag and a probe

            } // potential Reco probe size > 0

	  } // matched ele tag
	} // potential ele tag

	/*
        for(size_t m=0;m<iEvent->eleInfo.size();m++) {
	  if (verbose_) std::cout << "FOUND Ele " << iEvent->eleInfo[m].p4.Pt() <<  " " << EVENT.event << " " << candW->size() << " " << candZ->size() << std::endl;
        }
	cout << "Before PFMET" << endl;
	cout << "PFMET: " << iEvent->pfmet.p4.Et() << endl;
	*/

	// Trigger name debugging
	std::vector<std::string> interestingNames;
	//	interestingNames.push_back("HLT_IsoMu20_eta2p1_WCandPt80");
	//	interestingNames.push_back("HLT_IsoMu24_eta2p1_CentralPFJet30_CentralPFJet25_PFMET20");
	//	interestingNames.push_back("HLT_Ele27_WP80_WCandPt80");
	//	interestingNames.push_back("HLT_Ele27_WP80_CentralPFJet30_CentralPFJet25_PFMET20");
	interestingNames.push_back("HLT_IsoMu24_eta2p1_v");
	interestingNames.push_back("HLT_Mu40_eta2p1_v");
	interestingNames.push_back("HLT_Mu17_Mu8_v");
	interestingNames.push_back("HLT_Mu17_TkMu8_v");
	interestingNames.push_back("HLT_IsoMu20_eta2p1_WCandPt80_v");
	interestingNames.push_back("HLT_Ele27_WP80_v");
	interestingNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");

	//	const TriggerPathRefVector trigPathRefs (triggerEvent->pathRefs() );

	for ( TriggerPathRefVector::const_iterator trigPathRef  = trigPathRefs.begin(); trigPathRef != trigPathRefs.end(); ++trigPathRef ) {
	std::string name = (*trigPathRef)->name();
	bool interesting = false;
	std::string interestingName;
	for (std::vector<std::string>::iterator intName = interestingNames.begin() ; intName != interestingNames.end() ; intName++) {
	if (name.compare(0,intName->size(),*intName) == 0) {
        interesting = true;
        interestingName = *intName;
	}
	}
	if (interesting) {
	if (verbose_) std::cout << name << " " << (*trigPathRef)->wasRun() << " " << (*trigPathRef)->wasAccept() << std::endl;
	const TriggerFilterRefVector trigFilterRefs( triggerEvent->pathFilters(name) );
	for ( TriggerFilterRefVector::const_iterator trigFilterRef  = trigFilterRefs.begin(); trigFilterRef != trigFilterRefs.end(); ++trigFilterRef ) {
        std::string label = (*trigFilterRef)->label();
        if (verbose_) std::cout << "  " << label << " " << (*trigFilterRef)->status() << std::endl;
	//      if (verbose_) std::cout << "Looking for muon label " << muonLabels[interestingName] << " (and label is " << label << ")" << std::endl;
	//      if (muonLabels.count(name) && label == muonLabels[interestingName]) {
	//        if (verbose_) std::cout << "Found label " << muonLabels[interestingName] << "(and label is "<< label << ")" << std::endl;
	//        const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
	//        for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
	//          trigMuons[interestingName].push_back(GENPTOLORP(*trigObjRef));
	//          if (verbose_) std::cout << "   Pushing muon back to " << interestingName << std::endl;
	//        }
	//      }
        const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
        for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
	if (verbose_) std::cout << "    " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
        }
	} 
	}
	}


	} // closed event loop

    std::cout << "closing the file: " << inputFiles_[iFile] << std::endl;
    inFile->Close();
    // close input file
  } // loop on files
     
  
  std::cout << "Events: " << ievt <<std::endl;
  std::cout << "TotalCount: " << totalcount <<std::endl;

    
    
  _outFile->cd();
  muTrigDir->cd();  muTrigTree->Write();
  muWCandDir->cd();  muWCandTree->Write();
  muRecoDir->cd();  muRecoTree->Write();
  eleTrigDir->cd();  eleTrigTree->Write();
  eleRecoDir->cd();  eleRecoTree->Write();
  _outFile->Write();
  _outFile->Close();
  return 0;
}


