#include <TH1F.h>
#include <TH3F.h>
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
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
float ElectronIso(const VHbbEvent::ElectronInfo &i,float rho) {
  float mincor=0.0;
  float minrho=0.0;
  float rhoN = std::max(rho,minrho);
  float eta=i.p4.Eta();
  float areagamma=0.5;
  float areaNH=0.5;
  float areaComb=0.5;

  if(fabs(eta) <= 1.0 ) {areagamma=0.14; areaNH=0.044; areaComb=0.18;}
  if(fabs(eta) > 1.0 &&  fabs(eta) <= 1.479 ) {areagamma=0.13; areaNH=0.065; areaComb=0.20;}
  if(fabs(eta) > 1.479 &&  fabs(eta) <= 2.0 ) {areagamma=0.079; areaNH=0.068; areaComb=0.15;}
  if(fabs(eta) > 2.0 &&  fabs(eta) <= 2.2 ) {areagamma=0.13; areaNH=0.057; areaComb=0.19;}
  if(fabs(eta) > 2.2 &&  fabs(eta) <= 2.3 ) {areagamma=0.15; areaNH=0.058; areaComb=0.21;}
  if(fabs(eta) > 2.3 &&  fabs(eta) <= 2.4 ) {areagamma=0.16; areaNH=0.061; areaComb=0.22;}
  if(fabs(eta) > 2.4  ) {areagamma=0.18; areaNH=0.11; areaComb=0.29;}


  /*
  if(fabs(eta) <= 1.0 ) {areagamma=0.081; areaNH=0.024; areaComb=0.10;}
  if(fabs(eta) > 1.0 &&  fabs(eta) <= 1.479 ) {areagamma=0.084; areaNH=0.037; areaComb=0.12;}
  if(fabs(eta) > 1.479 &&  fabs(eta) <= 2.0 ) {areagamma=0.048; areaNH=0.037; areaComb=0.085;}
  if(fabs(eta) > 2.0 &&  fabs(eta) <= 2.2 ) {areagamma=0.089; areaNH=0.023; areaComb=0.11;}
  if(fabs(eta) > 2.2 &&  fabs(eta) <= 2.3 ) {areagamma=0.092; areaNH=0.023; areaComb=0.12;}
  if(fabs(eta) > 2.3 &&  fabs(eta) <= 2.4 ) {areagamma=0.097; areaNH=0.021; areaComb=0.12;}
  if(fabs(eta) > 2.4  ) {areagamma=0.11; areaNH=0.021; areaComb=0.13;}
  */

  //Correct electron photon double count
  float pho=i.pfPhoIso;
  if(i.innerHits>0) 
    { 
      pho-=i.pfPhoIsoDoubleCounted;
    }
  
  float pfCorrIso = (i.pfChaIso+ std::max(pho-rhoN*areagamma,mincor )+std::max(i.pfNeuIso-rhoN*areaNH,mincor))/i.p4.Pt();
  return pfCorrIso;
}

bool ElectronPresel(const VHbbEvent::ElectronInfo &i) {
  bool presel =   ((i.isEE  &&
                    //fabs(i.Deta) < 0.009 &&
                    //fabs(i.Dphi) < 0.1 &&
                    i.sihih < 0.03  &&
                    i.HoE < 0.10  &&
                    i.innerHits == 0  &&
                    (i.tIso/i.p4.Pt()) < 0.2 &&
                    (i.eIso/i.p4.Pt()) < 0.2 &&
                    (i.hIso/i.p4.Pt()) < 0.2)
                   ||
                   (i.isEB &&
                    //fabs(i.Deta) < 0.007 &&
                    //fabs(i.Dphi) < 0.015 &&
                    i.sihih < 0.01  &&
                    i.HoE < 0.12  &&
                    i.innerHits == 0  &&
                    (i.tIso/i.p4.Pt()) < 0.2 &&
                    (i.eIso/i.p4.Pt()) < 0.2 &&
                    (i.hIso/i.p4.Pt()) < 0.2)
                   );
  return presel;
}

bool ElectronWP(const VHbbEvent::ElectronInfo &i,float rho,int wp, bool requireIso = true,bool requireID = true, bool requirePresel=true) {

  float iso=ElectronIso(i,rho);
  float id=i.mvaOutTrig;
  float eta=i.p4.Eta();
  bool wp70=((fabs(eta) < 0.8 && id>0.977 && iso < 0.093) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.956 && iso < 0.095) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.966 && iso < 0.171));
  bool wp80=((fabs(eta) < 0.8 && id>0.913 && iso < 0.105) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.964 && iso < 0.178) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.899 && iso < 0.150));
  bool wp85=((fabs(eta) < 0.8 && id>0.929 && iso < 0.135) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.931 && iso < 0.159) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.805 && iso < 0.155));
  bool wp90=((fabs(eta) < 0.8 && id>0.877 && iso < 0.177) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.794 && iso < 0.180) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.846 && iso < 0.244));
  bool wp95=((fabs(eta) < 0.8 && id>0.858 && iso < 0.253) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.425 && iso < 0.225) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.759 && iso < 0.308));
  bool wpHWW=((fabs(eta) < 0.8 && id>0.94 && iso < 0.15) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.85 && iso < 0.15) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.92 && iso < 0.15));

  bool presel = ElectronPresel(i);

  if (!requireIso && !requireID) return presel;

  if (!requireIso) {
    //    cout << "Not requiring Iso; id=" << id << " iso=" <<iso << endl;
    wp70=((fabs(eta) < 0.8 && id>0.977) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.956) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.966));
    wp80=((fabs(eta) < 0.8 && id>0.913) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.964) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.899));
    wp85=((fabs(eta) < 0.8 && id>0.929) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.931) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.805));
    wp90=((fabs(eta) < 0.8 && id>0.877) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.794) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.846));
    wp95=((fabs(eta) < 0.8 && id>0.858) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.425) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.759));
    wpHWW=((fabs(eta) < 0.8 && id>0.94) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.85) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.92));
  }

  if (!requireID) {
    //    cout << "Not requiring ID; id=" << id << " iso=" << iso << endl;
    wp70=((fabs(eta) < 0.8 && iso < 0.093) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && iso < 0.095) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && iso < 0.171));
    wp80=((fabs(eta) < 0.8 && iso < 0.105) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && iso < 0.178) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && iso < 0.150));
    wp85=((fabs(eta) < 0.8 && iso < 0.135) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && iso < 0.159) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && iso < 0.155));
    wp90=((fabs(eta) < 0.8 && iso < 0.177) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && iso < 0.180) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && iso < 0.244));
    wp95=((fabs(eta) < 0.8 && iso < 0.253) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && iso < 0.225) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && iso < 0.308));
    wpHWW=((fabs(eta) < 0.8 && iso < 0.15) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && iso < 0.15) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && iso < 0.15));
  }

  if (requirePresel && ! presel) {
    return false;
  }

  if (wp == 70) return (wp70);
  if (wp == 80) return (wp80);
  if (wp == 85) return (wp85);
  if (wp == 90) return (wp90);
  if (wp == 95) return (wp95);
  return (wpHWW);  // use HWW as default
}

// Copied from muon selection and LeptonInfo::setSpecific
float muon2012PfCorrIso(const VHbbEvent::MuonInfo & i, float rho) {
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
  return pfCorrIso;
}

bool muonId2012Tight(const VHbbEvent::MuonInfo & i, float rho, bool requireiso = true) {
  float pfCorrIso = muon2012PfCorrIso(i,rho);
  return (i.isPF && i. globChi2<10 && i.nPixelHits>= 1 && i.globNHits != 0 && i.nValidLayers > 5 &&         (i.cat & 0x1) && i.nMatches >=2 && i.ipDb<.2 
	  && fabs(i.p4.Eta()) < 2.4 && i.p4.Pt() > 20. &&
	  i.zPVPt < 0.5 && (!requireiso || pfCorrIso < 0.12)); // Added by SCZ from selection
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
  float rho;
  float rhoForEleIso;
  int nPVs;
  float met;
  int nJets, nBJets;

  float probe_dPhiDiJet, probe_diJetPt, probe_dPhiJetMin, probe_closestJetPt;

  // float  PUweight, PUweight2011B;
  //  float PU0,PUp1,PUm1;

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
  float PUweight, PUweight2011B,PUweight1DObs;
  float PU0,PUp1,PUm1;


  /*
  edm::LumiReWeighting   lumiWeights;
///  edm::Lumi3DReWeighting   lumiWeights2011B;
  if(isMC_)
    {
      lumiWeights = edm::LumiReWeighting(PUmcfileName_,PUdatafileName_ , "pileup", "pileup");

   //   lumiWeights2011B = edm::Lumi3DReWeighting(PUmcfileName2011B_,PUdatafileName2011B_ , "pileup", "pileup");
    //  if(Weight3DfileName_!="")
//	{ lumiWeights2011B.weight3D_init(Weight3DfileName_.c_str()); }
 //     else
//	{
//	  lumiWeights2011B.weight3D_init(1.0); // generate the weights the fisrt time;
//	}

    }
  */

  edm::LumiReWeighting   lumiWeights;
  edm::LumiReWeighting   lumiWeights1DObs;
//  edm::Lumi3DReWeighting   lumiWeights2011B;
  if(isMC_)
    {
      lumiWeights = edm::LumiReWeighting(PUmcfileName_,PUdatafileName_ , "pileup", "pileup");
      lumiWeights1DObs = edm::LumiReWeighting(PUmcfileName2011B_,PUdatafileName2011B_ , "pileup", "pileup");

/*      lumiWeights2011B = edm::Lumi3DReWeighting(PUmcfileName2011B_,PUdatafileName2011B_ , "pileup", "pileup");
      if(Weight3DfileName_!="")
	{ lumiWeights2011B.weight3D_init(Weight3DfileName_.c_str()); }
      else
	{
	  lumiWeights2011B.weight3D_init(1.0); // generate the weights the fisrt time;
	}
*/
    }

  TFile *_outFile	= new TFile(outputFile_.c_str(), "recreate");	
  TDirectory *muTrigDir = _outFile->mkdir("muTrigDir");
  TDirectory *muWCandDir =  _outFile->mkdir("muWCandDir");
  TDirectory *muRecoDir = _outFile->mkdir("muRecoDir");
  TDirectory *muRecoIsoDir = _outFile->mkdir("muRecoIsoDir");
  TDirectory * eleTrigDir = _outFile->mkdir("eleTrigDir");
  TDirectory *eleRecoDir = _outFile->mkdir("eleRecoDir");
  TDirectory *eleRecoIsoDir = _outFile->mkdir("eleRecoIsoDir");

  int eventNumber, runNumber, lb;

  float muTrigTree_tag_pt, muTrigTree_tag_eta, muTrigTree_probe_pt, muTrigTree_probe_eta, muTrigTree_mass, muTrigTree_tag_phi, muTrigTree_probe_phi;
  int muTrigTree_probe_passingIsoMu24L1, muTrigTree_probe_passingIsoMu24L2, muTrigTree_probe_passingIsoMu24L3, muTrigTree_probe_passingIsoMu24Iso;
  int muTrigTree_probe_passingMu40L1, muTrigTree_probe_passingMu40L2, muTrigTree_probe_passingMu40L3;
  int muTrigTree_probe_passingDiMuL1, muTrigTree_probe_passingDiMuL20, muTrigTree_probe_passingDiMuL210, muTrigTree_probe_passingDiMu8, muTrigTree_probe_passingDiMu17, muTrigTree_probe_passingDiMuDz;
  int muTrigTree_probe_passingDiMuTkL1, muTrigTree_probe_passingDiMuTkL2, muTrigTree_probe_passingDiMuTk17, muTrigTree_probe_passingDiMuTk8, muTrigTree_probe_passingDiMuTkDz;
  int muTrigTree_probe_passingDiMuL110, muTrigTree_probe_passingDiMuTkL110, muTrigTree_probe_passingDiMu17Dz, muTrigTree_probe_passingDiMuTk17Dz;
  int muTrigTree_probe_passingIsoMu20Iso;
  int muTrigTree_probe_passingDiMuL13p5;
  int muTrigTree_probe_passingDiMuL110_match3, muTrigTree_probe_passingDiMuL13p5_match3, muTrigTree_probe_passingDiMuL1_match3;
  int muTrigTree_probe_passingDiMuL110_match5, muTrigTree_probe_passingDiMuL13p5_match5, muTrigTree_probe_passingDiMuL1_match5;
  int muTrigTree_probe_passingDiMuL110_matchinf, muTrigTree_probe_passingDiMuL13p5_matchinf, muTrigTree_probe_passingDiMuL1_matchinf;
  float muTrigTree_probe_WCandPt; int muTrigTree_probe_passingWCandPt; 
  int muTrigTree_event_Mu17_Mu8, muTrigTree_event_Mu17_TkMu8;
  int muTrigTree_probe_passingIsoMu24ORMu40, muTrigTree_probe_passingMu40ANDNOTIsoMu24, muTrigTree_probe_passingIsoMu20ANDNOTIsoMu24;
  int muTrigTree_probe_passingIsoMu20ANDNOTIsoMu24ANDNOTMu40, muTrigTree_probe_passingIsoMu24ORMu40ORIsoMu20;
  int muTrigTree_probe_passingDiMuTk17ANDNOTeventMu17Mu8, muTrigTree_probe_passingDiMuTk8ANDNOTeventMu17Mu8, muTrigTree_probe_passingDiMuTkDzANDNOTeventMu17Mu8;
  float muTrigTree_dR, muTrigTree_ZPt, muTrigTree_probe_abseta;
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
  muTrigTree->Branch("probe_passingDiMuL13p5", &muTrigTree_probe_passingDiMuL13p5, "probe_passingDiMuL13p5/I");
  muTrigTree->Branch("probe_passingDiMuL1", &muTrigTree_probe_passingDiMuL1, "probe_passingDiMuL1/I");

  muTrigTree->Branch("probe_passingDiMuL110_match3", &muTrigTree_probe_passingDiMuL110_match3, "probe_passingDiMuL110_match3/I");
  muTrigTree->Branch("probe_passingDiMuL13p5_match3", &muTrigTree_probe_passingDiMuL13p5_match3, "probe_passingDiMuL13p5_match3/I");
  muTrigTree->Branch("probe_passingDiMuL1_match3", &muTrigTree_probe_passingDiMuL1_match3, "probe_passingDiMuL1_match3/I");

  muTrigTree->Branch("probe_passingDiMuL110_match5", &muTrigTree_probe_passingDiMuL110_match5, "probe_passingDiMuL110_match5/I");
  muTrigTree->Branch("probe_passingDiMuL13p5_match5", &muTrigTree_probe_passingDiMuL13p5_match5, "probe_passingDiMuL13p5_match5/I");
  muTrigTree->Branch("probe_passingDiMuL1_match5", &muTrigTree_probe_passingDiMuL1_match5, "probe_passingDiMuL1_match5/I");

  muTrigTree->Branch("probe_passingDiMuL110_matchinf", &muTrigTree_probe_passingDiMuL110_matchinf, "probe_passingDiMuL110_matchinf/I");
  muTrigTree->Branch("probe_passingDiMuL13p5_matchinf", &muTrigTree_probe_passingDiMuL13p5_matchinf, "probe_passingDiMuL13p5_matchinf/I");
  muTrigTree->Branch("probe_passingDiMuL1_matchinf", &muTrigTree_probe_passingDiMuL1_matchinf, "probe_passingDiMuL1_matchinf/I");


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
  
  muTrigTree->Branch("dR", &muTrigTree_dR, "dR/F");
  muTrigTree->Branch("mass", &muTrigTree_mass, "mass/F");
  muTrigTree->Branch("ZPt", &muTrigTree_ZPt, "ZPt/F");
  muTrigTree->Branch("probe_abseta",&muTrigTree_probe_abseta, "probe_abseta/F");
  muTrigTree->Branch("nPVs", &nPVs, "nPVs/I");
  muTrigTree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  muTrigTree->Branch("runNumber", &runNumber, "runNumber/I");
  muTrigTree->Branch("lb", &lb, "lb/I");
  muTrigTree->Branch("lumiWeight", &lumiWeight, "lumiWeight/F");
  muTrigTree->Branch("PUweight", &PUweight, "PUweight/F");
  muTrigTree->Branch("PUweight2011B", &PUweight2011B, "PUweight2011B/F");
  muTrigTree->Branch("PUweight1DObs", &PUweight1DObs, "PUweight1DObs/F");

  muTrigTree->Branch("nJets", &nJets, "nJets/I");
  muTrigTree->Branch("nBJets", &nBJets, "nBJets/I");
  muTrigTree->Branch("met", &met, "met/F");
  muTrigTree->Branch("probe_dPhiDiJet", &probe_dPhiDiJet, "probe_dPhiDiJet/F");
  muTrigTree->Branch("probe_diJetPt", &probe_diJetPt, "probe_diJetPt/F");
  muTrigTree->Branch("probe_dPhiJetMin", &probe_dPhiJetMin, "probe_dPhiJetMin/F");
  muTrigTree->Branch("probe_closestJetPt", &probe_closestJetPt, "probe_closestJetPt/F");


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
  muWCandTree->Branch("PUweight", &PUweight, "PUweight/F");
  muWCandTree->Branch("PUweight2011B", &PUweight2011B, "PUweight2011B/F");
  muWCandTree->Branch("PUweight1DObs", &PUweight1DObs, "PUweight1DObs/F");



  //  return (i.isPF && i. globChi2<10 && i.nPixelHits>= 1 && i.globNHits != 0 && i.nValidLayers > 5 &&         (i.cat & 0x1) && i.nMatches >=2 && i.ipDb<.2
  //          && fabs(eta) < 2.4 && i.p4.Pt() > 20. &&
  //          /* (i.cat & 0x1) && */ (!requireiso || pfCorrIso < 0.12)); // Added by SCZ from selection
  // muon2012PfCorrIso

  float muRecoTree_tag_pt, muRecoTree_tag_eta, muRecoTree_probe_pt, muRecoTree_probe_eta, muRecoTree_mass, muRecoTree_tag_phi, muRecoTree_probe_phi, muRecoTree_dR;
  int muRecoTree_probe_tightNoIso, muRecoTree_probe_tightIso, muRecoTree_probe_nearTightNoIso, muRecoTree_probe_nearTightIso;
  float muRecoTree_probe_Iso, muRecoTree_probe_ipDb, muRecoTree_probe_globChi2, muRecoTree_probe_zPVPt;
  int muRecoTree_probe_isPF, muRecoTree_probe_tracker, muRecoTree_probe_nPixelHits, muRecoTree_probe_nValidLayers, muRecoTree_probe_nMatches;
  int muRecoTree_probe_passingIsoMu24L3, muRecoTree_probe_passingIsoMu24Iso;
  float muRecoTree_probe_abseta, muRecoTree_ZPt;
  TTree *muRecoTree = new TTree("muRecoTree","muRecoTree"); muRecoTree->SetDirectory(muRecoDir);
  muRecoTree->Branch("tag_pt"            ,  &muRecoTree_tag_pt                 ,  "tag_pt/F");
  muRecoTree->Branch("tag_eta"            ,  &muRecoTree_tag_eta                 ,  "tag_eta/F");
  muRecoTree->Branch("tag_phi"            ,  &muRecoTree_tag_phi                 ,  "tag_phi/F");
  muRecoTree->Branch("probe_pt"            ,  &muRecoTree_probe_pt                 ,  "probe_pt/F");
  muRecoTree->Branch("probe_eta"            ,  &muRecoTree_probe_eta                 ,  "probe_eta/F");
  muRecoTree->Branch("probe_phi"            ,  &muRecoTree_probe_phi                 ,  "probe_phi/F");
  muRecoTree->Branch("probe_tightNoIso"            ,  &muRecoTree_probe_tightNoIso                 ,  "probe_tightNoIso/I");
  muRecoTree->Branch("probe_tightIso"            ,  &muRecoTree_probe_tightIso                 ,  "probe_tightIso/I");
  muRecoTree->Branch("probe_nearTightNoIso"            ,  &muRecoTree_probe_nearTightNoIso                 ,  "probe_nearTightNoIso/I");
  muRecoTree->Branch("probe_nearTightIso"            ,  &muRecoTree_probe_nearTightIso                 ,  "probe_nearTightIso/I");
  muRecoTree->Branch("probe_Iso"            ,  &muRecoTree_probe_Iso                 ,  "probe_Iso/F");
  muRecoTree->Branch("probe_ipDb"            ,  &muRecoTree_probe_ipDb                 ,  "probe_ipDb/F");
  muRecoTree->Branch("probe_zPVPt"            ,  &muRecoTree_probe_zPVPt                 ,  "probe_zPVPt/F");
  muRecoTree->Branch("probe_globChi2"            ,  &muRecoTree_probe_globChi2                 ,  "probe_globChi2/F");
  muRecoTree->Branch("probe_isPF"            ,  &muRecoTree_probe_isPF                 ,  "probe_isPF/I");
  muRecoTree->Branch("probe_tracker"            ,  &muRecoTree_probe_tracker                 ,  "probe_tracker/I");
  muRecoTree->Branch("probe_nPixelHits"            ,  &muRecoTree_probe_nPixelHits                 ,  "probe_nPixelHits/I");
  muRecoTree->Branch("probe_nValidLayers"            ,  &muRecoTree_probe_nValidLayers                 ,  "probe_nValidLayers/I");
  muRecoTree->Branch("probe_nMatches"            ,  &muRecoTree_probe_nMatches                 ,  "probe_nMatches/I");
  muRecoTree->Branch("probe_passingIsoMu24L3", &muRecoTree_probe_passingIsoMu24L3, "probe_passingIsoMu24L3/I");
  muRecoTree->Branch("probe_passingIsoMu24Iso", &muRecoTree_probe_passingIsoMu24Iso, "probe_passingIsoMu24Iso/I");
  muRecoTree->Branch("dR", &muRecoTree_dR, "dR/F");
  muRecoTree->Branch("ZPt", &muRecoTree_ZPt, "ZPt/F");
  muRecoTree->Branch("probe_abseta",&muRecoTree_probe_abseta, "probe_abseta/F");
  muRecoTree->Branch("mass", &muRecoTree_mass, "mass/F");
  muRecoTree->Branch("nPVs", &nPVs, "nPVs/I");
  muRecoTree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  muRecoTree->Branch("runNumber", &runNumber, "runNumber/I");
  muRecoTree->Branch("lb", &lb, "lb/I");
  muRecoTree->Branch("lumiWeight", &lumiWeight, "lumiWeight/F");
  muRecoTree->Branch("PUweight", &PUweight, "PUweight/F");
  muRecoTree->Branch("PUweight2011B", &PUweight2011B, "PUweight2011B/F");
  muRecoTree->Branch("PUweight1DObs", &PUweight1DObs, "PUweight1DObs/F");
  muRecoTree->Branch("nJets", &nJets, "nJets/I");
  muRecoTree->Branch("nBJets", &nBJets, "nBJets/I");
  muRecoTree->Branch("met", &met, "met/F");
  muRecoTree->Branch("rho"            ,  &rho                 ,  "rho/F");
  muRecoTree->Branch("rhoN"            ,  &rhoN                 ,  "rhoN/F");
  muRecoTree->Branch("probe_dPhiDiJet", &probe_dPhiDiJet, "probe_dPhiDiJet/F");
  muRecoTree->Branch("probe_diJetPt", &probe_diJetPt, "probe_diJetPt/F");
  muRecoTree->Branch("probe_dPhiJetMin", &probe_dPhiJetMin, "probe_dPhiJetMin/F");
  muRecoTree->Branch("probe_closestJetPt", &probe_closestJetPt, "probe_closestJetPt/F");


  float muRecoIsoTree_tag_pt, muRecoIsoTree_tag_eta, muRecoIsoTree_probe_pt, muRecoIsoTree_probe_eta, muRecoIsoTree_mass, muRecoIsoTree_tag_phi, muRecoIsoTree_probe_phi, muRecoIsoTree_dR;
  int muRecoIsoTree_probe_tightNoIso, muRecoIsoTree_probe_tightIso, muRecoIsoTree_probe_nearTightNoIso, muRecoIsoTree_probe_nearTightIso;
  float muRecoIsoTree_probe_Iso, muRecoIsoTree_probe_ipDb, muRecoIsoTree_probe_globChi2, muRecoIsoTree_probe_zPVPt;
  int muRecoIsoTree_probe_isPF, muRecoIsoTree_probe_tracker, muRecoIsoTree_probe_nPixelHits, muRecoIsoTree_probe_nValidLayers, muRecoIsoTree_probe_nMatches;
  int muRecoIsoTree_probe_passingIsoMu24L3, muRecoIsoTree_probe_passingIsoMu24Iso;
  float muRecoIsoTree_probe_abseta, muRecoIsoTree_ZPt;
  TTree *muRecoIsoTree = new TTree("muRecoIsoTree","muRecoIsoTree"); muRecoIsoTree->SetDirectory(muRecoIsoDir);
  muRecoIsoTree->Branch("tag_pt"            ,  &muRecoIsoTree_tag_pt                 ,  "tag_pt/F");
  muRecoIsoTree->Branch("tag_eta"            ,  &muRecoIsoTree_tag_eta                 ,  "tag_eta/F");
  muRecoIsoTree->Branch("tag_phi"            ,  &muRecoIsoTree_tag_phi                 ,  "tag_phi/F");
  muRecoIsoTree->Branch("probe_pt"            ,  &muRecoIsoTree_probe_pt                 ,  "probe_pt/F");
  muRecoIsoTree->Branch("probe_eta"            ,  &muRecoIsoTree_probe_eta                 ,  "probe_eta/F");
  muRecoIsoTree->Branch("probe_phi"            ,  &muRecoIsoTree_probe_phi                 ,  "probe_phi/F");
  muRecoIsoTree->Branch("probe_tightNoIso"            ,  &muRecoIsoTree_probe_tightNoIso                 ,  "probe_tightNoIso/I");
  muRecoIsoTree->Branch("probe_tightIso"            ,  &muRecoIsoTree_probe_tightIso                 ,  "probe_tightIso/I");
  muRecoIsoTree->Branch("probe_nearTightNoIso"            ,  &muRecoIsoTree_probe_nearTightNoIso                 ,  "probe_nearTightNoIso/I");
  muRecoIsoTree->Branch("probe_nearTightIso"            ,  &muRecoIsoTree_probe_nearTightIso                 ,  "probe_nearTightIso/I");
  muRecoIsoTree->Branch("probe_Iso"            ,  &muRecoIsoTree_probe_Iso                 ,  "probe_Iso/F");
  muRecoIsoTree->Branch("probe_ipDb"            ,  &muRecoIsoTree_probe_ipDb                 ,  "probe_ipDb/F");
  muRecoIsoTree->Branch("probe_zPVPt"            ,  &muRecoIsoTree_probe_zPVPt                 ,  "probe_zPVPt/F");
  muRecoIsoTree->Branch("probe_globChi2"            ,  &muRecoIsoTree_probe_globChi2                 ,  "probe_globChi2/F");
  muRecoIsoTree->Branch("probe_isPF"            ,  &muRecoIsoTree_probe_isPF                 ,  "probe_isPF/I");
  muRecoIsoTree->Branch("probe_tracker"            ,  &muRecoIsoTree_probe_tracker                 ,  "probe_tracker/I");
  muRecoIsoTree->Branch("probe_nPixelHits"            ,  &muRecoIsoTree_probe_nPixelHits                 ,  "probe_nPixelHits/I");
  muRecoIsoTree->Branch("probe_nValidLayers"            ,  &muRecoIsoTree_probe_nValidLayers                 ,  "probe_nValidLayers/I");
  muRecoIsoTree->Branch("probe_nMatches"            ,  &muRecoIsoTree_probe_nMatches                 ,  "probe_nMatches/I");
  muRecoIsoTree->Branch("probe_passingIsoMu24L3", &muRecoIsoTree_probe_passingIsoMu24L3, "probe_passingIsoMu24L3/I");
  muRecoIsoTree->Branch("probe_passingIsoMu24Iso", &muRecoIsoTree_probe_passingIsoMu24Iso, "probe_passingIsoMu24Iso/I");
  muRecoIsoTree->Branch("dR", &muRecoIsoTree_dR, "dR/F");
  muRecoIsoTree->Branch("ZPt", &muRecoIsoTree_ZPt, "ZPt/F");
  muRecoIsoTree->Branch("probe_abseta",&muRecoIsoTree_probe_abseta, "probe_abseta/F");
  muRecoIsoTree->Branch("mass", &muRecoIsoTree_mass, "mass/F");
  muRecoIsoTree->Branch("nPVs", &nPVs, "nPVs/I");
  muRecoIsoTree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  muRecoIsoTree->Branch("runNumber", &runNumber, "runNumber/I");
  muRecoIsoTree->Branch("lb", &lb, "lb/I");
  muRecoIsoTree->Branch("lumiWeight", &lumiWeight, "lumiWeight/F");
  muRecoIsoTree->Branch("nJets", &nJets, "nJets/I");
  muRecoIsoTree->Branch("nBJets", &nBJets, "nBJets/I");
  muRecoIsoTree->Branch("met", &met, "met/F");
  muRecoIsoTree->Branch("PUweight", &PUweight, "PUweight/F");
  muRecoIsoTree->Branch("PUweight2011B", &PUweight2011B, "PUweight2011B/F");
  muRecoIsoTree->Branch("PUweight1DObs", &PUweight1DObs, "PUweight1DObs/F");
  muRecoIsoTree->Branch("rho"            ,  &rho                 ,  "rho/F");
  muRecoIsoTree->Branch("rhoN"            ,  &rhoN                 ,  "rhoN/F");
  muRecoIsoTree->Branch("probe_dPhiDiJet", &probe_dPhiDiJet, "probe_dPhiDiJet/F");
  muRecoIsoTree->Branch("probe_diJetPt", &probe_diJetPt, "probe_diJetPt/F");
  muRecoIsoTree->Branch("probe_dPhiJetMin", &probe_dPhiJetMin, "probe_dPhiJetMin/F");
  muRecoIsoTree->Branch("probe_closestJetPt", &probe_closestJetPt, "probe_closestJetPt/F");


  float eleTrigTree_tag_pt, eleTrigTree_tag_eta, eleTrigTree_probe_pt, eleTrigTree_probe_eta, eleTrigTree_mass, eleTrigTree_tag_phi, eleTrigTree_probe_phi, eleTrigTree_dR;
  int eleTrigTree_probe_passingEle27L1, eleTrigTree_probe_passingEle27HLT;
  int eleTrigTree_probe_passingDiEle17, eleTrigTree_probe_passingDiEle8, eleTrigTree_probe_passingDiEleDz;
  int eleTrigTree_probe_wp70, eleTrigTree_probe_wp80, eleTrigTree_probe_wp85, eleTrigTree_probe_wp90, eleTrigTree_probe_wp95;
  int eleTrigTree_probe_wpHWW;
  float eleTrigTree_probe_abseta, eleTrigTree_ZPt;
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
  eleTrigTree->Branch("probe_wpHWW"            ,  &eleTrigTree_probe_wpHWW                 ,  "probe_wpHWW/I");
  eleTrigTree->Branch("dR", &eleTrigTree_dR, "dR/F");
  eleTrigTree->Branch("ZPt", &eleTrigTree_ZPt, "ZPt/F");
  eleTrigTree->Branch("probe_abseta",&eleTrigTree_probe_abseta, "probe_abseta/F");
  eleTrigTree->Branch("mass", &eleTrigTree_mass, "mass/F");
  eleTrigTree->Branch("nPVs", &nPVs, "nPVs/I");
  eleTrigTree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  eleTrigTree->Branch("runNumber", &runNumber, "runNumber/I");
  eleTrigTree->Branch("lb", &lb, "lb/I");
  eleTrigTree->Branch("lumiWeight", &lumiWeight, "lumiWeight/F");
  eleTrigTree->Branch("nJets", &nJets, "nJets/I");
  eleTrigTree->Branch("nBJets", &nBJets, "nBJets/I");
  eleTrigTree->Branch("met", &met, "met/F");
  eleTrigTree->Branch("PUweight", &PUweight, "PUweight/F");
  eleTrigTree->Branch("PUweight2011B", &PUweight2011B, "PUweight2011B/F");
  eleTrigTree->Branch("PUweight1DObs", &PUweight1DObs, "PUweight1DObs/F");
  eleTrigTree->Branch("probe_dPhiDiJet", &probe_dPhiDiJet, "probe_dPhiDiJet/F");
  eleTrigTree->Branch("probe_diJetPt", &probe_diJetPt, "probe_diJetPt/F");
  eleTrigTree->Branch("probe_dPhiJetMin", &probe_dPhiJetMin, "probe_dPhiJetMin/F");
  eleTrigTree->Branch("probe_closestJetPt", &probe_closestJetPt, "probe_closestJetPt/F");



  float eleRecoTree_tag_pt, eleRecoTree_tag_eta, eleRecoTree_probe_pt, eleRecoTree_probe_eta, eleRecoTree_mass, eleRecoTree_tag_phi, eleRecoTree_probe_phi, eleRecoTree_dR;
  int eleRecoTree_probe_wp70, eleRecoTree_probe_wp80, eleRecoTree_probe_wp85, eleRecoTree_probe_wp90, eleRecoTree_probe_wp95;
  int eleRecoTree_probe_nearWp70, eleRecoTree_probe_nearWp80, eleRecoTree_probe_nearWp85, eleRecoTree_probe_nearWp90, eleRecoTree_probe_nearWp95;
  int eleRecoTree_probe_wp70noIso, eleRecoTree_probe_wp80noIso, eleRecoTree_probe_wp85noIso, eleRecoTree_probe_wp90noIso, eleRecoTree_probe_wp95noIso;
  int eleRecoTree_probe_wp70noId, eleRecoTree_probe_wp80noId, eleRecoTree_probe_wp85noId, eleRecoTree_probe_wp90noId, eleRecoTree_probe_wp95noId, eleRecoTree_probe_wpHWWnoId;
  int eleRecoTree_probe_wp70noPresel, eleRecoTree_probe_wp80noPresel, eleRecoTree_probe_wp85noPresel, eleRecoTree_probe_wp90noPresel, eleRecoTree_probe_wp95noPresel, eleRecoTree_probe_wpHWWnoPresel;
  int eleRecoTree_probe_passingEle27HLT;
  int eleRecoTree_probe_wpHWW, eleRecoTree_probe_nearWpHWW, eleRecoTree_probe_wpHWWnoIso;
  float eleRecoTree_probe_iso, eleRecoTree_probe_id;
  int eleRecoTree_probe_presel;
  float eleRecoTree_probe_abseta, eleRecoTree_ZPt;
  TTree *eleRecoTree = new TTree("eleRecoTree","eleRecoTree"); eleRecoTree->SetDirectory(eleRecoDir);
  eleRecoTree->Branch("tag_pt"            ,  &eleRecoTree_tag_pt                 ,  "tag_pt/F");
  eleRecoTree->Branch("tag_eta"            ,  &eleRecoTree_tag_eta                 ,  "tag_eta/F");
  eleRecoTree->Branch("tag_phi"            ,  &eleRecoTree_tag_phi                 ,  "tag_phi/F");
  eleRecoTree->Branch("probe_pt"            ,  &eleRecoTree_probe_pt                 ,  "probe_pt/F");
  eleRecoTree->Branch("probe_eta"            ,  &eleRecoTree_probe_eta                 ,  "probe_eta/F");
  eleRecoTree->Branch("probe_phi"            ,  &eleRecoTree_probe_phi                 ,  "probe_phi/F");
  eleRecoTree->Branch("probe_presel"            ,  &eleRecoTree_probe_presel                 ,  "probe_presel/I");
  eleRecoTree->Branch("probe_wp70"            ,  &eleRecoTree_probe_wp70                 ,  "probe_wp70/I");
  eleRecoTree->Branch("probe_wp80"            ,  &eleRecoTree_probe_wp80                 ,  "probe_wp80/I");
  eleRecoTree->Branch("probe_wp85"            ,  &eleRecoTree_probe_wp85                 ,  "probe_wp85/I");
  eleRecoTree->Branch("probe_wp90"            ,  &eleRecoTree_probe_wp90                 ,  "probe_wp90/I");
  eleRecoTree->Branch("probe_wp95"            ,  &eleRecoTree_probe_wp95                 ,  "probe_wp95/I");
  eleRecoTree->Branch("probe_wpHWW"            ,  &eleRecoTree_probe_wpHWW                 ,  "probe_wpHWW/I");
  eleRecoTree->Branch("probe_nearWp70"            ,  &eleRecoTree_probe_nearWp70                 ,  "probe_nearWp70/I");
  eleRecoTree->Branch("probe_nearWp80"            ,  &eleRecoTree_probe_nearWp80                 ,  "probe_nearWp80/I");
  eleRecoTree->Branch("probe_nearWp85"            ,  &eleRecoTree_probe_nearWp85                 ,  "probe_nearWp85/I");
  eleRecoTree->Branch("probe_nearWp90"            ,  &eleRecoTree_probe_nearWp90                 ,  "probe_nearWp90/I");
  eleRecoTree->Branch("probe_nearWp95"            ,  &eleRecoTree_probe_nearWp95                 ,  "probe_nearWp95/I");
  eleRecoTree->Branch("probe_nearWpHWW"            ,  &eleRecoTree_probe_nearWpHWW                 ,  "probe_nearWpHWW/I");
  eleRecoTree->Branch("probe_wp70noIso"            ,  &eleRecoTree_probe_wp70noIso                 ,  "probe_wp70noIso/I");
  eleRecoTree->Branch("probe_wp80noIso"            ,  &eleRecoTree_probe_wp80noIso                 ,  "probe_wp80noIso/I");
  eleRecoTree->Branch("probe_wp85noIso"            ,  &eleRecoTree_probe_wp85noIso                 ,  "probe_wp85noIso/I");
  eleRecoTree->Branch("probe_wp90noIso"            ,  &eleRecoTree_probe_wp90noIso                 ,  "probe_wp90noIso/I");
  eleRecoTree->Branch("probe_wp95noIso"            ,  &eleRecoTree_probe_wp95noIso                 ,  "probe_wp95noIso/I");
  eleRecoTree->Branch("probe_wpHWWnoIso"            ,  &eleRecoTree_probe_wpHWWnoIso                 ,  "probe_wpHWWnoIso/I");
  eleRecoTree->Branch("probe_wp70noId"            ,  &eleRecoTree_probe_wp70noId                 ,  "probe_wp70noId/I");
  eleRecoTree->Branch("probe_wp80noId"            ,  &eleRecoTree_probe_wp80noId                 ,  "probe_wp80noId/I");
  eleRecoTree->Branch("probe_wp85noId"            ,  &eleRecoTree_probe_wp85noId                 ,  "probe_wp85noId/I");
  eleRecoTree->Branch("probe_wp90noId"            ,  &eleRecoTree_probe_wp90noId                 ,  "probe_wp90noId/I");
  eleRecoTree->Branch("probe_wp95noId"            ,  &eleRecoTree_probe_wp95noId                 ,  "probe_wp95noId/I");
  eleRecoTree->Branch("probe_wpHWWnoId"            ,  &eleRecoTree_probe_wpHWWnoId                 ,  "probe_wpHWWnoId/I");
  eleRecoTree->Branch("probe_wp70noPresel"            ,  &eleRecoTree_probe_wp70noPresel                 ,  "probe_wp70noPresel/I");
  eleRecoTree->Branch("probe_wp80noPresel"            ,  &eleRecoTree_probe_wp80noPresel                 ,  "probe_wp80noPresel/I");
  eleRecoTree->Branch("probe_wp85noPresel"            ,  &eleRecoTree_probe_wp85noPresel                 ,  "probe_wp85noPresel/I");
  eleRecoTree->Branch("probe_wp90noPresel"            ,  &eleRecoTree_probe_wp90noPresel                 ,  "probe_wp90noPresel/I");
  eleRecoTree->Branch("probe_wp95noPresel"            ,  &eleRecoTree_probe_wp95noPresel                 ,  "probe_wp95noPresel/I");
  eleRecoTree->Branch("probe_wpHWWnoPresel"            ,  &eleRecoTree_probe_wpHWWnoPresel                 ,  "probe_wpHWWnoPresel/I");
  eleRecoTree->Branch("rho"            ,  &rho                 ,  "rho/F");
  eleRecoTree->Branch("rhoN"            ,  &rhoN                 ,  "rhoN/F");
  eleRecoTree->Branch("rhoForEleIso"            ,  &rhoForEleIso                 ,  "rhoForEleIso/F");
  eleRecoTree->Branch("probe_passingEle27HLT", &eleRecoTree_probe_passingEle27HLT, "probe_passingEle27HLT/I");
  eleRecoTree->Branch("probe_iso"            ,  &eleRecoTree_probe_iso                 ,  "probe_iso/F");
  eleRecoTree->Branch("probe_id"            ,  &eleRecoTree_probe_id                 ,  "probe_id/F");
  eleRecoTree->Branch("dR", &eleRecoTree_dR, "dR/F");
  eleRecoTree->Branch("ZPt", &eleRecoTree_ZPt, "ZPt/F");
  eleRecoTree->Branch("probe_abseta",&eleRecoTree_probe_abseta, "probe_abseta/F");
  eleRecoTree->Branch("mass", &eleRecoTree_mass, "mass/F");
  eleRecoTree->Branch("nPVs", &nPVs, "nPVs/I");
  eleRecoTree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  eleRecoTree->Branch("runNumber", &runNumber, "runNumber/I");
  eleRecoTree->Branch("lb", &lb, "lb/I");
  eleRecoTree->Branch("lumiWeight", &lumiWeight, "lumiWeight/F");
  eleRecoTree->Branch("nBJets", &nBJets, "nBJets/I");
  eleRecoTree->Branch("nJets", &nJets, "nJets/I");
  eleRecoTree->Branch("met", &met, "met/F");
  eleRecoTree->Branch("PUweight", &PUweight, "PUweight/F");
  eleRecoTree->Branch("PUweight2011B", &PUweight2011B, "PUweight2011B/F");
  eleRecoTree->Branch("PUweight1DObs", &PUweight1DObs, "PUweight1DObs/F");
  eleRecoTree->Branch("probe_dPhiDiJet", &probe_dPhiDiJet, "probe_dPhiDiJet/F");
  eleRecoTree->Branch("probe_diJetPt", &probe_diJetPt, "probe_diJetPt/F");
  eleRecoTree->Branch("probe_dPhiJetMin", &probe_dPhiJetMin, "probe_dPhiJetMin/F");
  eleRecoTree->Branch("probe_closestJetPt", &probe_closestJetPt, "probe_closestJetPt/F");



  float eleRecoIsoTree_tag_pt, eleRecoIsoTree_tag_eta, eleRecoIsoTree_probe_pt, eleRecoIsoTree_probe_eta, eleRecoIsoTree_mass, eleRecoIsoTree_tag_phi, eleRecoIsoTree_probe_phi, eleRecoIsoTree_dR;
  int eleRecoIsoTree_probe_wp70, eleRecoIsoTree_probe_wp80, eleRecoIsoTree_probe_wp85, eleRecoIsoTree_probe_wp90, eleRecoIsoTree_probe_wp95;
  int eleRecoIsoTree_probe_nearWp70, eleRecoIsoTree_probe_nearWp80, eleRecoIsoTree_probe_nearWp85, eleRecoIsoTree_probe_nearWp90, eleRecoIsoTree_probe_nearWp95;
  int eleRecoIsoTree_probe_wp70noIso, eleRecoIsoTree_probe_wp80noIso, eleRecoIsoTree_probe_wp85noIso, eleRecoIsoTree_probe_wp90noIso, eleRecoIsoTree_probe_wp95noIso;
  int eleRecoIsoTree_probe_wp70noId, eleRecoIsoTree_probe_wp80noId, eleRecoIsoTree_probe_wp85noId, eleRecoIsoTree_probe_wp90noId, eleRecoIsoTree_probe_wp95noId, eleRecoIsoTree_probe_wpHWWnoId;
  int eleRecoIsoTree_probe_passingEle27HLT;
  int eleRecoIsoTree_probe_wpHWW, eleRecoIsoTree_probe_nearWpHWW, eleRecoIsoTree_probe_wpHWWnoIso;
  float eleRecoIsoTree_probe_iso, eleRecoIsoTree_probe_id;
  float eleRecoIsoTree_probe_abseta, eleRecoIsoTree_ZPt;
  TTree *eleRecoIsoTree = new TTree("eleRecoIsoTree","eleRecoIsoTree"); eleRecoIsoTree->SetDirectory(eleRecoIsoDir);
  eleRecoIsoTree->Branch("tag_pt"            ,  &eleRecoIsoTree_tag_pt                 ,  "tag_pt/F");
  eleRecoIsoTree->Branch("tag_eta"            ,  &eleRecoIsoTree_tag_eta                 ,  "tag_eta/F");
  eleRecoIsoTree->Branch("tag_phi"            ,  &eleRecoIsoTree_tag_phi                 ,  "tag_phi/F");
  eleRecoIsoTree->Branch("probe_pt"            ,  &eleRecoIsoTree_probe_pt                 ,  "probe_pt/F");
  eleRecoIsoTree->Branch("probe_eta"            ,  &eleRecoIsoTree_probe_eta                 ,  "probe_eta/F");
  eleRecoIsoTree->Branch("probe_phi"            ,  &eleRecoIsoTree_probe_phi                 ,  "probe_phi/F");
  eleRecoIsoTree->Branch("probe_wp70"            ,  &eleRecoIsoTree_probe_wp70                 ,  "probe_wp70/I");
  eleRecoIsoTree->Branch("probe_wp80"            ,  &eleRecoIsoTree_probe_wp80                 ,  "probe_wp80/I");
  eleRecoIsoTree->Branch("probe_wp85"            ,  &eleRecoIsoTree_probe_wp85                 ,  "probe_wp85/I");
  eleRecoIsoTree->Branch("probe_wp90"            ,  &eleRecoIsoTree_probe_wp90                 ,  "probe_wp90/I");
  eleRecoIsoTree->Branch("probe_wp95"            ,  &eleRecoIsoTree_probe_wp95                 ,  "probe_wp95/I");
  eleRecoIsoTree->Branch("probe_wpHWW"            ,  &eleRecoIsoTree_probe_wpHWW                 ,  "probe_wpHWW/I");
  eleRecoIsoTree->Branch("probe_nearWp70"            ,  &eleRecoIsoTree_probe_nearWp70                 ,  "probe_nearWp70/I");
  eleRecoIsoTree->Branch("probe_nearWp80"            ,  &eleRecoIsoTree_probe_nearWp80                 ,  "probe_nearWp80/I");
  eleRecoIsoTree->Branch("probe_nearWp85"            ,  &eleRecoIsoTree_probe_nearWp85                 ,  "probe_nearWp85/I");
  eleRecoIsoTree->Branch("probe_nearWp90"            ,  &eleRecoIsoTree_probe_nearWp90                 ,  "probe_nearWp90/I");
  eleRecoIsoTree->Branch("probe_nearWp95"            ,  &eleRecoIsoTree_probe_nearWp95                 ,  "probe_nearWp95/I");
  eleRecoIsoTree->Branch("probe_nearWpHWW"            ,  &eleRecoIsoTree_probe_nearWpHWW                 ,  "probe_nearWpHWW/I");
  eleRecoIsoTree->Branch("probe_wp70noIso"            ,  &eleRecoIsoTree_probe_wp70noIso                 ,  "probe_wp70noIso/I");
  eleRecoIsoTree->Branch("probe_wp80noIso"            ,  &eleRecoIsoTree_probe_wp80noIso                 ,  "probe_wp80noIso/I");
  eleRecoIsoTree->Branch("probe_wp85noIso"            ,  &eleRecoIsoTree_probe_wp85noIso                 ,  "probe_wp85noIso/I");
  eleRecoIsoTree->Branch("probe_wp90noIso"            ,  &eleRecoIsoTree_probe_wp90noIso                 ,  "probe_wp90noIso/I");
  eleRecoIsoTree->Branch("probe_wp95noIso"            ,  &eleRecoIsoTree_probe_wp95noIso                 ,  "probe_wp95noIso/I");
  eleRecoIsoTree->Branch("probe_wpHWWnoIso"            ,  &eleRecoIsoTree_probe_wpHWWnoIso                 ,  "probe_wpHWWnoIso/I");
  eleRecoIsoTree->Branch("probe_wp70noId"            ,  &eleRecoIsoTree_probe_wp70noId                 ,  "probe_wp70noId/I");
  eleRecoIsoTree->Branch("probe_wp80noId"            ,  &eleRecoIsoTree_probe_wp80noId                 ,  "probe_wp80noId/I");
  eleRecoIsoTree->Branch("probe_wp85noId"            ,  &eleRecoIsoTree_probe_wp85noId                 ,  "probe_wp85noId/I");
  eleRecoIsoTree->Branch("probe_wp90noId"            ,  &eleRecoIsoTree_probe_wp90noId                 ,  "probe_wp90noId/I");
  eleRecoIsoTree->Branch("probe_wp95noId"            ,  &eleRecoIsoTree_probe_wp95noId                 ,  "probe_wp95noId/I");
  eleRecoIsoTree->Branch("probe_wpHWWnoId"            ,  &eleRecoIsoTree_probe_wpHWWnoId                 ,  "probe_wpHWWnoId/I");
  eleRecoIsoTree->Branch("rho"            ,  &rho                 ,  "rho/F");
  eleRecoIsoTree->Branch("rhoN"            ,  &rhoN                 ,  "rhoN/F");
  eleRecoIsoTree->Branch("rhoForEleIso"            ,  &rhoForEleIso                 ,  "rhoForEleIso/F");
  eleRecoIsoTree->Branch("probe_passingEle27HLT", &eleRecoIsoTree_probe_passingEle27HLT, "probe_passingEle27HLT/I");
  eleRecoIsoTree->Branch("probe_iso"            ,  &eleRecoIsoTree_probe_iso                 ,  "probe_iso/F");
  eleRecoIsoTree->Branch("probe_id"            ,  &eleRecoIsoTree_probe_id                 ,  "probe_id/F");
  eleRecoIsoTree->Branch("dR", &eleRecoIsoTree_dR, "dR/F");
  eleRecoIsoTree->Branch("ZPt", &eleRecoIsoTree_ZPt, "ZPt/F");
  eleRecoIsoTree->Branch("probe_abseta",&eleRecoIsoTree_probe_abseta, "probe_abseta/F");
  eleRecoIsoTree->Branch("mass", &eleRecoIsoTree_mass, "mass/F");
  eleRecoIsoTree->Branch("nPVs", &nPVs, "nPVs/I");
  eleRecoIsoTree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  eleRecoIsoTree->Branch("runNumber", &runNumber, "runNumber/I");
  eleRecoIsoTree->Branch("lb", &lb, "lb/I");
  eleRecoIsoTree->Branch("lumiWeight", &lumiWeight, "lumiWeight/F");
  eleRecoIsoTree->Branch("nJets", &nJets, "nJets/I");
  eleRecoIsoTree->Branch("nBJets", &nBJets, "nBJets/I");
  eleRecoIsoTree->Branch("met", &met, "met/F");
  eleRecoIsoTree->Branch("PUweight", &PUweight, "PUweight/F");
  eleRecoIsoTree->Branch("PUweight2011B", &PUweight2011B, "PUweight2011B/F");
  eleRecoIsoTree->Branch("PUweight1DObs", &PUweight1DObs, "PUweight1DObs/F");
  eleRecoIsoTree->Branch("probe_dPhiDiJet", &probe_dPhiDiJet, "probe_dPhiDiJet/F");
  eleRecoIsoTree->Branch("probe_diJetPt", &probe_diJetPt, "probe_diJetPt/F");
  eleRecoIsoTree->Branch("probe_dPhiJetMin", &probe_dPhiJetMin, "probe_dPhiJetMin/F");
  eleRecoIsoTree->Branch("probe_closestJetPt", &probe_closestJetPt, "probe_closestJetPt/F");



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

      /*
      if(isMC_){
 
	// PU weights // Run2011A
	std::map<int, unsigned int>::const_iterator puit = aux.puInfo.pus.find(0);
	int npu =puit->second ;
	//PUweight =  lumiWeights.weight( npu );        
	//pu->Fill(puit->second);
	std::map<int, unsigned int>::const_iterator puit0 =  aux.puInfo.pus.find(0);
	std::map<int, unsigned int>::const_iterator puitm1 = aux.puInfo.pus.find(-1);
	std::map<int, unsigned int>::const_iterator puitp1 = aux.puInfo.pus.find(+1);
	//PU0=puit0->second;
	//PUp1=puitp1->second;
	//PUm1=puitm1->second;
	//input3DPU->Fill(PUm1,PU0,PUp1);  
	lumiWeight = lumiWeights2011B.weight3D( puitm1->second, puit0->second,puitp1->second); 

      }
      */

      PUweight=1.;
      PUweight2011B=1.;
      PUweight1DObs=1.;
      if(isMC_){
	cout << "truePU = " << aux.puInfo.truePU << endl;
	// PU weights // Run2011A
	std::map<int, unsigned int>::const_iterator puit = aux.puInfo.pus.find(0);
	int npu =puit->second ;
	PUweight =  lumiWeights.weight( (int) aux.puInfo.truePU ); //use new method with "true PU"        
	//	pu->Fill(puit->second);
	// PU weight Run2011B
	// PU weight Run2011B
	std::map<int, unsigned int>::const_iterator puit0 =  aux.puInfo.pus.find(0);
	std::map<int, unsigned int>::const_iterator puitm1 = aux.puInfo.pus.find(-1);
	std::map<int, unsigned int>::const_iterator puitp1 = aux.puInfo.pus.find(+1);
	PU0=puit0->second;
	PUp1=puitp1->second;
	PUm1=puitm1->second;
	//	input3DPU->Fill(PUm1,PU0,PUp1);  
//	lumiWeight = lumiWeights2011B.weight3D( puitm1->second, puit0->second,puitp1->second); 
//	PUweight2011B = lumiWeights2011B.weight3D( puitm1->second, puit0->second,puitp1->second); 
	PUweight1DObs = lumiWeights1DObs.weight( npu); 
      }

      if (verbose_) cout << "Weights: " << PUweight2011B << " " << PUweight1DObs << " " << PUweight << endl;

      runNumber = ev.id().run();
      eventNumber = ev.id().event();
      lb = ev.id().luminosityBlock();

      VHbbEvent modifiedEvent;;
      const VHbbEvent *  iEvent =0;

      fwlite::Handle< VHbbEvent > vhbbHandle;
      vhbbHandle.getByLabel(ev,"HbbAnalyzerNew");
      modifiedEvent = *vhbbHandle.product();

      iEvent = vhbbHandle.product();

      met = iEvent->pfmet.p4.Et();


	float trigMatchCut = 0.1;
	float L1MatchCut = 0.2;
	float recoNearMatchCut = 0.1;

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

        if (runNumber>= 193834) {
          mu40name = "HLT_Mu40_v";
          mu40nameL1 = "hltL1sMu16";
          mu40nameL2 = "hltL2fL1sMu16L1f0L2Filtered16Q";
          mu40nameL3 = "hltL3fL1sMu16L1f0L2f16QL3Filtered40Q";

          mu24name = "HLT_IsoMu24_v";
          mu24nameL1 = "hltL1sMu16";
          mu24nameL2 = "hltL2fL1sMu16L1f0L2Filtered16Q";
          mu24nameL3 = "hltL3fL1sMu16L1f0L2f16QL3Filtered24Q";
          mu24nameIso ="hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";

	  mu20nameIso = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f20L3crIsoRhoFiltered0p15"; // n.b. eta restriction not lifted until later, see run 194270 below
        }


	if (runNumber >= 194270) {
	  mu20name = "HLT_IsoMu20_WCandPt80_v";
	  mu20nameIso = "hltL3crIsoL1sMu16L1f0L2f16QL3f20QL3crIsoRhoFiltered0p15";
	  mu20nameWCand = "hlt2IsoMu20PFMHTPt80";
	}


	std::vector<TLorentzVector> mu24L1;
	std::vector<TLorentzVector> mu24L2;
	std::vector<TLorentzVector> mu24L3;
	std::vector<TLorentzVector> mu24Iso;

        std::vector<TLorentzVector> mu40L1;
        std::vector<TLorentzVector> mu40L2;
        std::vector<TLorentzVector> mu40L3;

        std::vector<TLorentzVector> diMuL110;
	std::vector<TLorentzVector> diMuL13p5;
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
		  if ((*trigObjRef)->pt() > 9.999999) diMuL110.push_back(GENPTOLORP(*trigObjRef));
		  if ((*trigObjRef)->pt() > 3.499999) diMuL13p5.push_back(GENPTOLORP(*trigObjRef));
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
	rho = aux.puInfo.rho25;
	rhoForEleIso = aux.puInfo.rho25Iso;

	if (verbose_) cout << " rhoN=" << rhoN << " rho=" << rho << " rhoForEleIso=" << rhoForEleIso << endl;
	

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
	      muTrigTree_probe_abseta = fabs(muTrigTree_probe_eta);
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
              muTrigTree_probe_passingDiMuL13p5 = 0;
              muTrigTree_probe_passingDiMuL1 = 0;

              muTrigTree_probe_passingDiMuL110_match3 = 0;
              muTrigTree_probe_passingDiMuL13p5_match3 = 0;
              muTrigTree_probe_passingDiMuL1_match3 = 0;

              muTrigTree_probe_passingDiMuL110_match5 = 0;
              muTrigTree_probe_passingDiMuL13p5_match5 = 0;
              muTrigTree_probe_passingDiMuL1_match5 = 0;

              muTrigTree_probe_passingDiMuL110_matchinf = 0;
              muTrigTree_probe_passingDiMuL13p5_matchinf = 0;
              muTrigTree_probe_passingDiMuL1_matchinf = 0;


	      muTrigTree_mass = (iEvent->muInfo[probeNumber].p4 + iEvent->muInfo[tagNumber].p4).M();
              muTrigTree_dR = iEvent->muInfo[probeNumber].p4.DeltaR(iEvent->muInfo[tagNumber].p4);
              muTrigTree_ZPt = (iEvent->muInfo[probeNumber].p4 + iEvent->muInfo[tagNumber].p4).Pt();



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
              for (size_t i=0 ; i < diMuL13p5.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuL13p5[i]) < L1MatchCut) muTrigTree_probe_passingDiMuL13p5 = 1;
              }

              for (size_t i=0 ; i < diMuL1.size(); i++) {
		if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuL1[i]) < 0.3) muTrigTree_probe_passingDiMuL1_match3 = 1;
              }
              for (size_t i=0 ; i < diMuL110.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuL110[i]) < 0.3) muTrigTree_probe_passingDiMuL110_match3 = 1;
              }
              for (size_t i=0 ; i < diMuL13p5.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuL13p5[i]) < 0.3) muTrigTree_probe_passingDiMuL13p5_match3 = 1;
              }
              for (size_t i=0 ; i < diMuL1.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuL1[i]) < 0.5) muTrigTree_probe_passingDiMuL1_match5 = 1;
              }
              for (size_t i=0 ; i < diMuL110.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuL110[i]) < 0.5) muTrigTree_probe_passingDiMuL110_match5 = 1;
              }
              for (size_t i=0 ; i < diMuL13p5.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuL13p5[i]) < 0.5) muTrigTree_probe_passingDiMuL13p5_match5 = 1;
              }

              for (size_t i=0 ; i < diMuL1.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuL1[i]) < 999999.) muTrigTree_probe_passingDiMuL1_matchinf = 1;
              }
              for (size_t i=0 ; i < diMuL110.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuL110[i]) < 999999.) muTrigTree_probe_passingDiMuL110_matchinf = 1;
              }
              for (size_t i=0 ; i < diMuL13p5.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(diMuL13p5[i]) < 999999.) muTrigTree_probe_passingDiMuL13p5_matchinf = 1;
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

	      nJets = 0;
	      nBJets = 0;
	      probe_dPhiJetMin = 99999.;
	      probe_closestJetPt = 0.;
	      std::vector<TLorentzVector> dijets;
	      for(size_t j=0; j< iEvent->simpleJets2.size(); j++) {
		if (iEvent->simpleJets2[j].p4.Pt() > 20. && fabs(iEvent->simpleJets2[j].p4.Eta()) < 2.5 
		    && iEvent->muInfo[probeNumber].p4.DeltaR(iEvent->simpleJets2[j].p4) > 0.5 && iEvent->muInfo[tagNumber].p4.DeltaR(iEvent->simpleJets2[j].p4) > 0.5) {
		  nJets++;
		  if (iEvent->simpleJets2[j].csv>0.4) nBJets++;
		  if (fabs(iEvent->muInfo[probeNumber].p4.DeltaPhi(iEvent->simpleJets2[j].p4)) < fabs(probe_dPhiJetMin)) {
		    probe_dPhiJetMin = iEvent->muInfo[probeNumber].p4.DeltaPhi(iEvent->simpleJets2[j].p4);
		    probe_closestJetPt = iEvent->simpleJets2[j].p4.Pt();
		  }
		  for(size_t k=j+1; k < iEvent->simpleJets2.size(); k++) {
		    if (iEvent->simpleJets2[k].p4.Pt() > 20. && fabs(iEvent->simpleJets2[k].p4.Eta()) < 2.5
			&& iEvent->muInfo[probeNumber].p4.DeltaR(iEvent->simpleJets2[k].p4) > 0.5 && iEvent->muInfo[tagNumber].p4.DeltaR(iEvent->simpleJets2[k].p4) > 0.5) {
		      dijets.push_back(iEvent->simpleJets2[j].p4+iEvent->simpleJets2[k].p4);
		    }
		  }
		}
	      }
	      probe_dPhiDiJet = 99999.;
	      probe_diJetPt = 0.;
	      for (size_t j=0; j < dijets.size() ; j++) {
		if (dijets[j].Pt() > probe_diJetPt) {
		  probe_dPhiDiJet = iEvent->muInfo[probeNumber].p4.DeltaPhi(dijets[j]);
		  probe_diJetPt = dijets[j].Pt();
		}
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
	  }
	} // potential MuTags;
            
	std::vector<unsigned int> potentialMuRecoIsoTags;
	
	for(size_t m=0;m<iEvent->muInfo.size();m++) {
	  if (muonId2012Tight(iEvent->muInfo[m],rhoN,false) && iEvent->muInfo[m].p4.Pt() > 20. && fabs(iEvent->muInfo[m].p4.Eta()) < 2.4) {
	    if (verbose_) std::cout << "FOUND RecoIso tag Mu " << iEvent->muInfo[m].p4.Pt() << " " << iEvent->muInfo[m].p4.Eta() << " " << iEvent->muInfo[m].p4.Phi() << endl;
	    potentialMuRecoIsoTags.push_back(m);
	  }
	}
        
	if (verbose_) cout << " potentialMuRecoIsoTags.size()=" << potentialMuRecoIsoTags.size() << endl;
        
	if (potentialMuRecoIsoTags.size()>0) {
	  unsigned int tagNumber = potentialMuRecoIsoTags[rand->Integer(potentialMuRecoIsoTags.size())];
	  if (verbose_) std::cout << " Picked RecoIso Tag Mu " << iEvent->muInfo[tagNumber].p4.Pt() << " " << iEvent->muInfo[tagNumber].p4.Eta() << " " << iEvent->muInfo[tagNumber].p4.Phi() << endl;
	  muRecoIsoTree_tag_pt = iEvent->muInfo[tagNumber].p4.Pt();
	  muRecoIsoTree_tag_eta = iEvent->muInfo[tagNumber].p4.Eta();
	  muRecoIsoTree_tag_phi = iEvent->muInfo[tagNumber].p4.Phi();
          
	  bool matched = false;
	  for (size_t i=0 ; i < mu24Iso.size(); i++) {
	    if (iEvent->muInfo[tagNumber].p4.DeltaR(mu24Iso[i]) < trigMatchCut) {
	      matched = true;
	      if (verbose_) cout << "   Matches trigger " << mu24Iso[i].Pt() << " " << mu24Iso[i].Eta() << " " << mu24Iso[i].Phi() << endl;
	    }
	  }
	  if (!matched) if (verbose_) cout << "   Has no match!" << endl;
	  if (!muonId2012Tight(iEvent->muInfo[tagNumber],rhoN,false) && verbose_) {
	    cout << "Does not pass tightNoIso cuts";
	  }
	  if (matched&&muonId2012Tight(iEvent->muInfo[tagNumber],rhoN)) {
	    
	    std::vector<unsigned int> potentialRecoIsoProbes;
	    for(size_t m=0;m<iEvent->muInfo.size();m++) {
	      if (m == tagNumber) continue;
	      float mass = (iEvent->muInfo[m].p4 + iEvent->muInfo[tagNumber].p4).M();
	      if (muonId2012Tight(iEvent->muInfo[m],rhoN,false) && iEvent->muInfo[m].p4.Pt() > 20. && fabs(iEvent->muInfo[m].p4.Eta()) < 2.4
		  && (iEvent->muInfo[m].charge + iEvent->muInfo[tagNumber].charge) == 0 && (mass > 50. && mass < 200.)) {
		potentialRecoIsoProbes.push_back(m);
		if (verbose_) {
		  cout << "  Found a potential RecoIso probe: " <<  iEvent->muInfo[m].p4.Pt() << " " << iEvent->muInfo[m].p4.Eta() << " " << iEvent->muInfo[m].p4.Phi() << endl;
		  cout << "  Inv mass: " << mass << " Tight: " << muonId2012Tight(iEvent->muInfo[m],rhoN,false) << " TightIso: " << muonId2012Tight(iEvent->muInfo[m],rhoN) 
		       << " Iso: " << muon2012PfCorrIso(iEvent->muInfo[m],rhoN) << endl;
		}
	      }
	    }
	    if (potentialRecoIsoProbes.size()>0) {
	      unsigned int probeNumber = potentialRecoIsoProbes[rand->Integer(potentialRecoIsoProbes.size())];
	      muRecoIsoTree_probe_pt= iEvent->muInfo[probeNumber].p4.Pt();
	      muRecoIsoTree_probe_eta = iEvent->muInfo[probeNumber].p4.Eta();
	      muRecoIsoTree_probe_abseta = fabs(muRecoIsoTree_probe_eta);
	      muRecoIsoTree_probe_phi = iEvent->muInfo[probeNumber].p4.Phi();
	      muRecoIsoTree_probe_tightNoIso = muonId2012Tight(iEvent->muInfo[probeNumber],rhoN,false);
	      muRecoIsoTree_probe_tightIso = muonId2012Tight(iEvent->muInfo[probeNumber],rhoN);
	      muRecoIsoTree_mass = (iEvent->muInfo[probeNumber].p4 + iEvent->muInfo[tagNumber].p4).M();
              muRecoIsoTree_ZPt = (iEvent->muInfo[probeNumber].p4 + iEvent->muInfo[tagNumber].p4).Pt();
	      muRecoIsoTree_dR = iEvent->muInfo[probeNumber].p4.DeltaR(iEvent->muInfo[tagNumber].p4);
              
              
	      muRecoIsoTree_probe_nearTightNoIso = 0;
	      muRecoIsoTree_probe_nearTightIso = 0;
	      for (size_t m=0;m<iEvent->muInfo.size();m++) {
		if (m == tagNumber || iEvent->muInfo[probeNumber].p4.DeltaR(iEvent->muInfo[m].p4) > recoNearMatchCut) continue;
		if (muonId2012Tight(iEvent->muInfo[m],rhoN,false)) muRecoIsoTree_probe_nearTightNoIso = 1;
		if (muonId2012Tight(iEvent->muInfo[m],rhoN)) muRecoIsoTree_probe_nearTightIso = 1;
	      }
              
	      muRecoIsoTree_probe_Iso = muon2012PfCorrIso(iEvent->muInfo[probeNumber],rhoN);
	      muRecoIsoTree_probe_ipDb = iEvent->muInfo[probeNumber].ipDb;
	      muRecoIsoTree_probe_zPVPt = iEvent->muInfo[probeNumber].zPVPt;
	      muRecoIsoTree_probe_globChi2 = iEvent->muInfo[probeNumber].globChi2;
	      muRecoIsoTree_probe_isPF= iEvent->muInfo[probeNumber].isPF;
	      muRecoIsoTree_probe_tracker= (iEvent->muInfo[probeNumber].cat & 0x2);
	      muRecoIsoTree_probe_nPixelHits= iEvent->muInfo[probeNumber].nPixelHits;
	      muRecoIsoTree_probe_nValidLayers= iEvent->muInfo[probeNumber].nValidLayers;
	      muRecoIsoTree_probe_nMatches= iEvent->muInfo[probeNumber].nMatches;
              
	      muRecoIsoTree_probe_passingIsoMu24Iso = 0;
	      muRecoIsoTree_probe_passingIsoMu24L3 = 0;
	      for (size_t i=0 ; i < mu24L3.size(); i++) {
		if (iEvent->muInfo[probeNumber].p4.DeltaR(mu24L3[i]) < trigMatchCut) muRecoIsoTree_probe_passingIsoMu24L3 = 1;
	      }
	      for (size_t i=0 ; i < mu24Iso.size(); i++) {
		if (iEvent->muInfo[probeNumber].p4.DeltaR(mu24Iso[i]) < trigMatchCut) muRecoIsoTree_probe_passingIsoMu24Iso = 1;
	      }

              nJets = 0;
              nBJets = 0;
              probe_dPhiJetMin = 99999.;
              probe_closestJetPt = 0.;
	      std::vector<TLorentzVector> dijets;
              for(size_t j=0; j< iEvent->simpleJets2.size(); j++) {
                if (iEvent->simpleJets2[j].p4.Pt() > 20. && fabs(iEvent->simpleJets2[j].p4.Eta()) < 2.5
                    && iEvent->muInfo[probeNumber].p4.DeltaR(iEvent->simpleJets2[j].p4) > 0.5 && iEvent->muInfo[tagNumber].p4.DeltaR(iEvent->simpleJets2[j].p4) > 0.5) {
                  nJets++;
                  if (iEvent->simpleJets2[j].csv>0.4) nBJets++;
                  if (fabs(iEvent->muInfo[probeNumber].p4.DeltaPhi(iEvent->simpleJets2[j].p4)) < fabs(probe_dPhiJetMin)) {
                    probe_dPhiJetMin = iEvent->muInfo[probeNumber].p4.DeltaPhi(iEvent->simpleJets2[j].p4);
                    probe_closestJetPt = iEvent->simpleJets2[j].p4.Pt();
                  }
                  for(size_t k=j+1; k < iEvent->simpleJets2.size(); k++) {
                    if (iEvent->simpleJets2[k].p4.Pt() > 20. && fabs(iEvent->simpleJets2[k].p4.Eta()) < 2.5
                        && iEvent->muInfo[probeNumber].p4.DeltaR(iEvent->simpleJets2[k].p4) > 0.5 && iEvent->muInfo[tagNumber].p4.DeltaR(iEvent->simpleJets2[k].p4) > 0.5) {
                      dijets.push_back(iEvent->simpleJets2[j].p4+iEvent->simpleJets2[k].p4);
                    }
                  }
                }
              }
              probe_dPhiDiJet = 99999.;
              probe_diJetPt = 0.;
              for (size_t j=0; j < dijets.size() ; j++) {
                if (dijets[j].Pt() > probe_diJetPt) {
                  probe_dPhiDiJet = iEvent->muInfo[probeNumber].p4.DeltaPhi(dijets[j]);
                  probe_diJetPt = dijets[j].Pt();
                }
              }
             
	      muRecoIsoTree->Fill(); // Only if trigger-matched tag and a probe
              
	    } // potential RecoIso probe size > 0
            
	  } // trigger-matched tag
	} // potential tag size > 0
        
	
        std::vector<unsigned int> potentialMuRecoTags;
	
        for(size_t m=0;m<iEvent->muInfo.size();m++) {
          if ((iEvent->muInfo[m].cat & 0x1) && iEvent->muInfo[m].p4.Pt() > 20. && fabs(iEvent->muInfo[m].p4.Eta()) < 2.4) {
            if (verbose_) std::cout << "FOUND reco tag Mu " << iEvent->muInfo[m].p4.Pt() << " " << iEvent->muInfo[m].p4.Eta() << " " << iEvent->muInfo[m].p4.Phi() << endl;
            potentialMuRecoTags.push_back(m);
          }
        }

        if (verbose_) cout << " potentialMuRecoTags.size()=" << potentialMuRecoTags.size() << endl;

        if (potentialMuRecoTags.size()>0) {
          unsigned int tagNumber = potentialMuRecoTags[rand->Integer(potentialMuRecoTags.size())];
          if (verbose_) std::cout << " Picked reco Tag Mu " << iEvent->muInfo[tagNumber].p4.Pt() << " " << iEvent->muInfo[tagNumber].p4.Eta() << " " << iEvent->muInfo[tagNumber].p4.Phi() << endl;
          muRecoTree_tag_pt = iEvent->muInfo[tagNumber].p4.Pt();
          muRecoTree_tag_eta = iEvent->muInfo[tagNumber].p4.Eta();
          muRecoTree_tag_phi = iEvent->muInfo[tagNumber].p4.Phi();

          bool matched = false;
          for (size_t i=0 ; i < mu24Iso.size(); i++) {
            if (iEvent->muInfo[tagNumber].p4.DeltaR(mu24Iso[i]) < trigMatchCut) {
              matched = true;
              if (verbose_) cout << "   Matches trigger " << mu24Iso[i].Pt() << " " << mu24Iso[i].Eta() << " " << mu24Iso[i].Phi() << endl;
            }
          }
          if (!matched) if (verbose_) cout << "   Has no match!" << endl;
	  if (!muonId2012Tight(iEvent->muInfo[tagNumber],rhoN) && verbose_) {
	    cout << "Does not pass tight cuts";
	  }
          if (matched&&muonId2012Tight(iEvent->muInfo[tagNumber],rhoN,false)) {

	    std::vector<unsigned int> potentialRecoProbes;
            for(size_t m=0;m<iEvent->muInfo.size();m++) {
              if (m == tagNumber) continue;
              float mass = (iEvent->muInfo[m].p4 + iEvent->muInfo[tagNumber].p4).M();
              if ((iEvent->muInfo[m].cat & 0x1) && iEvent->muInfo[m].p4.Pt() > 20. && fabs(iEvent->muInfo[m].p4.Eta()) < 2.4
		  && (iEvent->muInfo[m].charge + iEvent->muInfo[tagNumber].charge) == 0 && (mass > 50. && mass < 200.)) {
		potentialRecoProbes.push_back(m);
		if (verbose_) {
		  cout << "  Found a potential reco probe: " <<  iEvent->muInfo[m].p4.Pt() << " " << iEvent->muInfo[m].p4.Eta() << " " << iEvent->muInfo[m].p4.Phi() << endl;
		  cout << "  Inv mass: " << mass << " Tight: " << muonId2012Tight(iEvent->muInfo[m],rhoN,false) << " TightIso: " << muonId2012Tight(iEvent->muInfo[m],rhoN) 
		       << " Iso: " << muon2012PfCorrIso(iEvent->muInfo[m],rhoN) << endl;
		}
              }
            }
            if (potentialRecoProbes.size()>0) {
              unsigned int probeNumber = potentialRecoProbes[rand->Integer(potentialRecoProbes.size())];
              muRecoTree_probe_pt= iEvent->muInfo[probeNumber].p4.Pt();
              muRecoTree_probe_eta = iEvent->muInfo[probeNumber].p4.Eta();
	      muRecoTree_probe_abseta = fabs(muRecoTree_probe_eta);
              muRecoTree_probe_phi = iEvent->muInfo[probeNumber].p4.Phi();
	      muRecoTree_probe_tightNoIso = muonId2012Tight(iEvent->muInfo[probeNumber],rhoN,false);
              muRecoTree_probe_tightIso = muonId2012Tight(iEvent->muInfo[probeNumber],rhoN);
              muRecoTree_mass = (iEvent->muInfo[probeNumber].p4 + iEvent->muInfo[tagNumber].p4).M();
	      muRecoTree_ZPt = (iEvent->muInfo[probeNumber].p4 + iEvent->muInfo[tagNumber].p4).Pt();
              muRecoTree_dR = iEvent->muInfo[probeNumber].p4.DeltaR(iEvent->muInfo[tagNumber].p4);


	      muRecoTree_probe_nearTightNoIso = 0;
	      muRecoTree_probe_nearTightIso = 0;
	      for (size_t m=0;m<iEvent->muInfo.size();m++) {
		if (m == tagNumber || iEvent->muInfo[probeNumber].p4.DeltaR(iEvent->muInfo[m].p4) > recoNearMatchCut) continue;
		if (muonId2012Tight(iEvent->muInfo[m],rhoN,false)) muRecoTree_probe_nearTightNoIso = 1;
                if (muonId2012Tight(iEvent->muInfo[m],rhoN)) muRecoTree_probe_nearTightIso = 1;
	      }

	      muRecoTree_probe_Iso = muon2012PfCorrIso(iEvent->muInfo[probeNumber],rhoN);
	      muRecoTree_probe_ipDb = iEvent->muInfo[probeNumber].ipDb;
              muRecoTree_probe_zPVPt = iEvent->muInfo[probeNumber].zPVPt;
	      muRecoTree_probe_globChi2 = iEvent->muInfo[probeNumber].globChi2;
	      muRecoTree_probe_isPF= iEvent->muInfo[probeNumber].isPF;
              muRecoTree_probe_tracker= (iEvent->muInfo[probeNumber].cat & 0x2);
              muRecoTree_probe_nPixelHits= iEvent->muInfo[probeNumber].nPixelHits;
              muRecoTree_probe_nValidLayers= iEvent->muInfo[probeNumber].nValidLayers;
              muRecoTree_probe_nMatches= iEvent->muInfo[probeNumber].nMatches;

	      muRecoTree_probe_passingIsoMu24Iso = 0;
	      muRecoTree_probe_passingIsoMu24L3 = 0;
              for (size_t i=0 ; i < mu24L3.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(mu24L3[i]) < trigMatchCut) muRecoTree_probe_passingIsoMu24L3 = 1;
              }
              for (size_t i=0 ; i < mu24Iso.size(); i++) {
                if (iEvent->muInfo[probeNumber].p4.DeltaR(mu24Iso[i]) < trigMatchCut) muRecoTree_probe_passingIsoMu24Iso = 1;
              }

              nJets = 0;
              nBJets = 0;
              probe_dPhiJetMin = 99999.;
              probe_closestJetPt = 0.;
	      std::vector<TLorentzVector> dijets;
              for(size_t j=0; j< iEvent->simpleJets2.size(); j++) {
                if (iEvent->simpleJets2[j].p4.Pt() > 20. && fabs(iEvent->simpleJets2[j].p4.Eta()) < 2.5
                    && iEvent->muInfo[probeNumber].p4.DeltaR(iEvent->simpleJets2[j].p4) > 0.5 && iEvent->muInfo[tagNumber].p4.DeltaR(iEvent->simpleJets2[j].p4) > 0.5) {
                  nJets++;
                  if (iEvent->simpleJets2[j].csv>0.4) nBJets++;
                  if (fabs(iEvent->muInfo[probeNumber].p4.DeltaPhi(iEvent->simpleJets2[j].p4)) < fabs(probe_dPhiJetMin)) {
                    probe_dPhiJetMin = iEvent->muInfo[probeNumber].p4.DeltaPhi(iEvent->simpleJets2[j].p4);
                    probe_closestJetPt = iEvent->simpleJets2[j].p4.Pt();
                  }
                  for(size_t k=j+1; k < iEvent->simpleJets2.size(); k++) {
                    if (iEvent->simpleJets2[k].p4.Pt() > 20. && fabs(iEvent->simpleJets2[k].p4.Eta()) < 2.5
                        && iEvent->muInfo[probeNumber].p4.DeltaR(iEvent->simpleJets2[k].p4) > 0.5 && iEvent->muInfo[tagNumber].p4.DeltaR(iEvent->simpleJets2[k].p4) > 0.5) {
                      dijets.push_back(iEvent->simpleJets2[j].p4+iEvent->simpleJets2[k].p4);
                    }
                  }
                }
              }
              probe_dPhiDiJet = 99999.;
              probe_diJetPt = 0.;
              for (size_t j=0; j < dijets.size() ; j++) {
                if (dijets[j].Pt() > probe_diJetPt) {
                  probe_dPhiDiJet = iEvent->muInfo[probeNumber].p4.DeltaPhi(dijets[j]);
                  probe_diJetPt = dijets[j].Pt();
                }
              }


	      muRecoTree->Fill(); // Only if trigger-matched tag and a probe

	    } // potential Reco probe size > 0

	  } // trigger-matched tag
	} // potential tag size > 0
            
            

        std::vector<unsigned int> potentialEleTags;

        for(size_t m=0;m<iEvent->eleInfo.size();m++) {
          if (ElectronWP(iEvent->eleInfo[m],rhoForEleIso,80) && iEvent->eleInfo[m].p4.Pt() > 20.) {
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
              if (ElectronWP(iEvent->eleInfo[m],rhoForEleIso,95) && (iEvent->eleInfo[m].charge + iEvent->eleInfo[tagNumber].charge) == 0 && (mass > 50. && mass < 200.)) {
                potentialHLTProbes.push_back(m);
              }
            }
            if (potentialHLTProbes.size()>0) {
              unsigned int probeNumber = potentialHLTProbes[rand->Integer(potentialHLTProbes.size())];
	      if (verbose_) std::cout << "   Picked Probe Ele " << iEvent->eleInfo[probeNumber].p4.Pt() << " " << iEvent->eleInfo[probeNumber].p4.Eta() << " " << iEvent->eleInfo[probeNumber].p4.Phi() << endl;
              eleTrigTree_probe_pt= iEvent->eleInfo[probeNumber].p4.Pt();
              eleTrigTree_probe_eta = iEvent->eleInfo[probeNumber].p4.Eta();
	      eleTrigTree_probe_abseta = fabs(eleTrigTree_probe_eta);
              eleTrigTree_probe_phi = iEvent->eleInfo[probeNumber].p4.Phi();
              eleTrigTree_probe_wpHWW = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,0);
              eleTrigTree_probe_wp70 = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,70);
              eleTrigTree_probe_wp80 = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,80);
              eleTrigTree_probe_wp85 = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,85);
              eleTrigTree_probe_wp90 = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,90);
              eleTrigTree_probe_wp95 = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,95);
	      eleTrigTree_probe_passingEle27L1 = 0;
	      eleTrigTree_probe_passingEle27HLT = 0;
	      eleTrigTree_probe_passingDiEle17 = 0;
	      eleTrigTree_probe_passingDiEle8 = 0;
	      eleTrigTree_probe_passingDiEleDz = 0;
	      eleTrigTree_mass = (iEvent->eleInfo[probeNumber].p4 + iEvent->eleInfo[tagNumber].p4).M();
	      eleTrigTree_ZPt = (iEvent->eleInfo[probeNumber].p4 + iEvent->eleInfo[tagNumber].p4).Pt();
              eleTrigTree_dR = iEvent->eleInfo[probeNumber].p4.DeltaR(iEvent->eleInfo[tagNumber].p4);

	      
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

              nJets = 0;
              nBJets = 0;
              probe_dPhiJetMin = 99999.;
              probe_closestJetPt = 0.;
	      std::vector<TLorentzVector> dijets;
              for(size_t j=0; j< iEvent->simpleJets2.size(); j++) {
                if (iEvent->simpleJets2[j].p4.Pt() > 20. && fabs(iEvent->simpleJets2[j].p4.Eta()) < 2.5
                    && iEvent->eleInfo[probeNumber].p4.DeltaR(iEvent->simpleJets2[j].p4) > 0.5 && iEvent->eleInfo[tagNumber].p4.DeltaR(iEvent->simpleJets2[j].p4) > 0.5) {
                  nJets++;
                  if (iEvent->simpleJets2[j].csv>0.4) nBJets++;
                  if (iEvent->eleInfo[probeNumber].p4.DeltaPhi(iEvent->simpleJets2[j].p4) < probe_dPhiJetMin) {
                    probe_dPhiJetMin = iEvent->eleInfo[probeNumber].p4.DeltaPhi(iEvent->simpleJets2[j].p4);
                    probe_closestJetPt = iEvent->simpleJets2[j].p4.Pt();
                  }
                  for(size_t k=j+1; k < iEvent->simpleJets2.size(); k++) {
                    if (iEvent->simpleJets2[k].p4.Pt() > 20. && fabs(iEvent->simpleJets2[k].p4.Eta()) < 2.5
                        && iEvent->eleInfo[probeNumber].p4.DeltaR(iEvent->simpleJets2[k].p4) > 0.5 && iEvent->eleInfo[tagNumber].p4.DeltaR(iEvent->simpleJets2[k].p4) > 0.5) {
                      dijets.push_back(iEvent->simpleJets2[j].p4+iEvent->simpleJets2[k].p4);
                    }
                  }
                }
              }
              probe_dPhiDiJet = 99999.;
              probe_diJetPt = 0.;
              for (size_t j=0; j < dijets.size() ; j++) {
                if (dijets[j].Pt() > probe_diJetPt) {
                  probe_dPhiDiJet = iEvent->eleInfo[probeNumber].p4.DeltaPhi(dijets[j]);
                  probe_diJetPt = dijets[j].Pt();
                }
              }

	      eleTrigTree->Fill();

	    } // HLT probes > 0
	  }
	} // potential ele tags 

        std::vector<unsigned int> potentialEleRecoTags;

        for(size_t m=0;m<iEvent->eleInfo.size();m++) {
          if (iEvent->eleInfo[m].p4.Pt() > 20. && fabs(iEvent->eleInfo[m].p4.Eta()) < 2.5) {
            if (verbose_) std::cout << "FOUND reco tag Ele " << iEvent->eleInfo[m].p4.Pt() << " " << iEvent->eleInfo[m].p4.Eta() << " " << iEvent->eleInfo[m].p4.Phi() << endl;
            potentialEleRecoTags.push_back(m);
          }
        }

        if (verbose_) cout << " potentialEleRecoTags.size()=" << potentialEleRecoTags.size() << endl;

        if (potentialEleRecoTags.size()>0) {
          unsigned int tagNumber = potentialEleRecoTags[rand->Integer(potentialEleRecoTags.size())];
          if (verbose_) std::cout << " Picked reco Tag Ele " << iEvent->eleInfo[tagNumber].p4.Pt() << " " << iEvent->eleInfo[tagNumber].p4.Eta() << " " << iEvent->eleInfo[tagNumber].p4.Phi() << endl;
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
	  if (!ElectronWP(iEvent->eleInfo[tagNumber],rhoForEleIso,80,false)&&verbose_) {
	    cout << "  Does not pass WP80 (no Iso) cuts" << endl;
	  }
          if (matched&&ElectronWP(iEvent->eleInfo[tagNumber],rhoForEleIso,80,false)) {
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
	      eleRecoTree_probe_presel = ElectronPresel(iEvent->eleInfo[probeNumber]);
              eleRecoTree_probe_pt= iEvent->eleInfo[probeNumber].p4.Pt();
              eleRecoTree_probe_eta = iEvent->eleInfo[probeNumber].p4.Eta();
	      eleRecoTree_probe_abseta = fabs(eleRecoTree_probe_eta);
              eleRecoTree_probe_phi = iEvent->eleInfo[probeNumber].p4.Phi();
              eleRecoTree_probe_wpHWW = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,0);
	      eleRecoTree_probe_wp70 = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,70);
	      eleRecoTree_probe_wp80 = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,80);
	      eleRecoTree_probe_wp85 = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,85);
	      eleRecoTree_probe_wp90 = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,90);
	      eleRecoTree_probe_wp95 = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,95);
              eleRecoTree_probe_wpHWWnoIso = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,0,false);
              eleRecoTree_probe_wp70noIso = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,70,false);
              eleRecoTree_probe_wp80noIso = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,80,false);
              eleRecoTree_probe_wp85noIso = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,85,false);
              eleRecoTree_probe_wp90noIso = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,90,false);
              eleRecoTree_probe_wp95noIso = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,95,false);
              eleRecoTree_probe_wpHWWnoId = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,0,true,false);
              eleRecoTree_probe_wp70noId = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,70,true,false);
              eleRecoTree_probe_wp80noId = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,80,true,false);
	      eleRecoTree_probe_wp85noId = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,85,true,false);
              eleRecoTree_probe_wp90noId = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,90,true,false);
	      eleRecoTree_probe_wp95noId = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,95,true,false);
              eleRecoTree_probe_wpHWWnoPresel = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,0,true,true,false);
              eleRecoTree_probe_wp70noPresel = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,70,true,true,false);
              eleRecoTree_probe_wp80noPresel = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,80,true,true,false);
              eleRecoTree_probe_wp85noPresel = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,85,true,true,false);
              eleRecoTree_probe_wp90noPresel = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,90,true,true,false);
              eleRecoTree_probe_wp95noPresel = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,95,true,true,false);
              eleRecoTree_mass = (iEvent->eleInfo[probeNumber].p4 + iEvent->eleInfo[tagNumber].p4).M();
	      eleRecoTree_ZPt = (iEvent->eleInfo[probeNumber].p4 + iEvent->eleInfo[tagNumber].p4).Pt();
              eleRecoTree_dR = iEvent->eleInfo[probeNumber].p4.DeltaR(iEvent->eleInfo[tagNumber].p4);
	      eleRecoTree_probe_iso = ElectronIso(iEvent->eleInfo[probeNumber],rhoForEleIso);
	      eleRecoTree_probe_id = iEvent->eleInfo[probeNumber].mvaOutTrig;

	      eleRecoTree_probe_nearWpHWW = 0;
	      eleRecoTree_probe_nearWp70 = 0;
              eleRecoTree_probe_nearWp80 = 0;
              eleRecoTree_probe_nearWp85 = 0;
              eleRecoTree_probe_nearWp90 = 0;
              eleRecoTree_probe_nearWp95 = 0;
	      for (size_t m=0;m<iEvent->eleInfo.size();m++) {
                if (m == tagNumber || iEvent->eleInfo[probeNumber].p4.DeltaR(iEvent->eleInfo[m].p4) > recoNearMatchCut) continue;
                if (ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,0)) eleRecoTree_probe_nearWpHWW = 1;
		if (ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,70)) eleRecoTree_probe_nearWp70 = 1;
		if (ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,80)) eleRecoTree_probe_nearWp80 = 1;
		if (ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,85)) eleRecoTree_probe_nearWp85 = 1;
		if (ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,90)) eleRecoTree_probe_nearWp90 = 1;
		if (ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,95)) eleRecoTree_probe_nearWp95 = 1;
              }

	      eleRecoTree_probe_passingEle27HLT = 0;
              for (size_t i=0 ; i < ele27HLT.size(); i++) {
                if (iEvent->eleInfo[probeNumber].p4.DeltaR(ele27HLT[i]) < trigMatchCut) eleRecoTree_probe_passingEle27HLT = 1;
              }

              nJets = 0;
              nBJets = 0;
	      probe_dPhiJetMin = 99999.;
              probe_closestJetPt = 0.;
	      std::vector<TLorentzVector> dijets;
              for(size_t j=0; j< iEvent->simpleJets2.size(); j++) {
		if (iEvent->simpleJets2[j].p4.Pt() > 20. && fabs(iEvent->simpleJets2[j].p4.Eta()) < 2.5
                    && iEvent->eleInfo[probeNumber].p4.DeltaR(iEvent->simpleJets2[j].p4) > 0.5 && iEvent->eleInfo[tagNumber].p4.DeltaR(iEvent->simpleJets2[j].p4) > 0.5) {
		  nJets++;
                  if (iEvent->simpleJets2[j].csv>0.4) nBJets++;
                  if (iEvent->eleInfo[probeNumber].p4.DeltaPhi(iEvent->simpleJets2[j].p4) < probe_dPhiJetMin) {
		    probe_dPhiJetMin = iEvent->eleInfo[probeNumber].p4.DeltaPhi(iEvent->simpleJets2[j].p4);
                    probe_closestJetPt = iEvent->simpleJets2[j].p4.Pt();
                  }
                  for(size_t k=j+1; k < iEvent->simpleJets2.size(); k++) {
                    if (iEvent->simpleJets2[k].p4.Pt() > 20. && fabs(iEvent->simpleJets2[k].p4.Eta()) < 2.5
                        && iEvent->eleInfo[probeNumber].p4.DeltaR(iEvent->simpleJets2[k].p4) > 0.5 && iEvent->eleInfo[tagNumber].p4.DeltaR(iEvent->simpleJets2[k].p4) > 0.5) {
                      dijets.push_back(iEvent->simpleJets2[j].p4+iEvent->simpleJets2[k].p4);
                    }
                  }
                }
              }
              probe_dPhiDiJet = 99999.;
              probe_diJetPt = 0.;
              for (size_t j=0; j < dijets.size() ; j++) {
                if (dijets[j].Pt() > probe_diJetPt) {
                  probe_dPhiDiJet = iEvent->eleInfo[probeNumber].p4.DeltaPhi(dijets[j]);
                  probe_diJetPt = dijets[j].Pt();
                }
              }


              eleRecoTree->Fill(); // Only if trigger-matched tag and a probe

            } // potential Reco probe size > 0

	  } // matched ele tag
	} // potential ele tag

            
	std::vector<unsigned int> potentialEleRecoIsoTags;
	
	for(size_t m=0;m<iEvent->eleInfo.size();m++) {
	  if (ElectronPresel(iEvent->eleInfo[m]) && iEvent->eleInfo[m].p4.Pt() > 20. && fabs(iEvent->eleInfo[m].p4.Eta()) < 2.5) {
	    if (verbose_) std::cout << "FOUND RecoIso tag Ele " << iEvent->eleInfo[m].p4.Pt() << " " << iEvent->eleInfo[m].p4.Eta() << " " << iEvent->eleInfo[m].p4.Phi() << endl;
	    potentialEleRecoIsoTags.push_back(m);
	  }
	}
        
	if (verbose_) cout << " potentialEleRecoIsoTags.size()=" << potentialEleRecoIsoTags.size() << endl;
        
	if (potentialEleRecoIsoTags.size()>0) {
	  unsigned int tagNumber = potentialEleRecoIsoTags[rand->Integer(potentialEleRecoIsoTags.size())];
	  if (verbose_) std::cout << " Picked RecoIso Tag Ele " << iEvent->eleInfo[tagNumber].p4.Pt() << " " << iEvent->eleInfo[tagNumber].p4.Eta() << " " << iEvent->eleInfo[tagNumber].p4.Phi() << endl;
	  eleRecoIsoTree_tag_pt = iEvent->eleInfo[tagNumber].p4.Pt();
	  eleRecoIsoTree_tag_eta = iEvent->eleInfo[tagNumber].p4.Eta();
	  eleRecoIsoTree_tag_phi = iEvent->eleInfo[tagNumber].p4.Phi();
          
	  bool matched = false;
	  for (size_t i=0 ; i < ele27HLT.size(); i++) {
	    if (iEvent->eleInfo[tagNumber].p4.DeltaR(ele27HLT[i]) < trigMatchCut) {
	      matched = true;
	      if (verbose_) cout << "   Matches trigger " << ele27HLT[i].Pt() << " " << ele27HLT[i].Eta() << " " << ele27HLT[i].Phi() << endl;
	    }
	  }
	  if (!matched) if (verbose_) cout << "   Has no match!" << endl;
	  if (!ElectronWP(iEvent->eleInfo[tagNumber],rhoForEleIso,80)&&verbose_) {
	    cout << "  Does not pass WP80 cuts" << endl;
	  }
	  if (matched&&ElectronWP(iEvent->eleInfo[tagNumber],rhoForEleIso,80)) {
	    std::vector<unsigned int> potentialRecoIsoProbes;
	    for(size_t m=0;m<iEvent->eleInfo.size();m++) {
	      if (m == tagNumber) continue;
	      float mass = (iEvent->eleInfo[m].p4 + iEvent->eleInfo[tagNumber].p4).M();
	      if ((ElectronPresel(iEvent->eleInfo[m]) && iEvent->eleInfo[m].p4.Pt() > 20. && fabs(iEvent->eleInfo[m].p4.Eta()) < 2.5)
		  && (iEvent->eleInfo[m].charge + iEvent->eleInfo[tagNumber].charge) == 0 && (mass > 50. && mass < 200.)) {
		potentialRecoIsoProbes.push_back(m);
	      }
	    }
	    if (potentialRecoIsoProbes.size()>0) {
	      unsigned int probeNumber = potentialRecoIsoProbes[rand->Integer(potentialRecoIsoProbes.size())];
	      eleRecoIsoTree_probe_pt= iEvent->eleInfo[probeNumber].p4.Pt();
	      eleRecoIsoTree_probe_eta = iEvent->eleInfo[probeNumber].p4.Eta();
	      eleRecoIsoTree_probe_abseta = fabs(eleRecoIsoTree_probe_eta);
	      eleRecoIsoTree_probe_phi = iEvent->eleInfo[probeNumber].p4.Phi();
	      eleRecoIsoTree_probe_wpHWW = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,0);
	      eleRecoIsoTree_probe_wp70 = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,70);
	      eleRecoIsoTree_probe_wp80 = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,80);
	      eleRecoIsoTree_probe_wp85 = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,85);
	      eleRecoIsoTree_probe_wp90 = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,90);
	      eleRecoIsoTree_probe_wp95 = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,95);
	      eleRecoIsoTree_probe_wpHWWnoIso = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,0,false);
	      eleRecoIsoTree_probe_wp70noIso = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,70,false);
	      eleRecoIsoTree_probe_wp80noIso = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,80,false);
	      eleRecoIsoTree_probe_wp85noIso = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,85,false);
	      eleRecoIsoTree_probe_wp90noIso = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,90,false);
	      eleRecoIsoTree_probe_wp95noIso = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,95,false);
              eleRecoIsoTree_probe_wpHWWnoId = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,0,true,false);
              eleRecoIsoTree_probe_wp70noId = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,70,true,false);
              eleRecoIsoTree_probe_wp80noId = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,80,true,false);
	      eleRecoIsoTree_probe_wp85noId = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,85,true,false);
              eleRecoIsoTree_probe_wp90noId = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,90,true,false);
	      eleRecoIsoTree_probe_wp95noId = ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,95,true,false);
	      eleRecoIsoTree_mass = (iEvent->eleInfo[probeNumber].p4 + iEvent->eleInfo[tagNumber].p4).M();
	      eleRecoIsoTree_ZPt = (iEvent->eleInfo[probeNumber].p4 + iEvent->eleInfo[tagNumber].p4).Pt();
	      eleRecoIsoTree_dR = iEvent->eleInfo[probeNumber].p4.DeltaR(iEvent->eleInfo[tagNumber].p4);
	      eleRecoIsoTree_probe_iso = ElectronIso(iEvent->eleInfo[tagNumber],rhoForEleIso);
	      eleRecoIsoTree_probe_id = iEvent->eleInfo[tagNumber].mvaOutTrig;

              
	      eleRecoIsoTree_probe_nearWpHWW = 0;
	      eleRecoIsoTree_probe_nearWp70 = 0;
	      eleRecoIsoTree_probe_nearWp80 = 0;
	      eleRecoIsoTree_probe_nearWp85 = 0;
	      eleRecoIsoTree_probe_nearWp90 = 0;
	      eleRecoIsoTree_probe_nearWp95 = 0;
	      for (size_t m=0;m<iEvent->eleInfo.size();m++) {
		if (m == tagNumber || iEvent->eleInfo[probeNumber].p4.DeltaR(iEvent->eleInfo[m].p4) > recoNearMatchCut) continue;
		if (ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,0)) eleRecoIsoTree_probe_nearWpHWW = 1;
		if (ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,70)) eleRecoIsoTree_probe_nearWp70 = 1;
		if (ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,80)) eleRecoIsoTree_probe_nearWp80 = 1;
		if (ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,85)) eleRecoIsoTree_probe_nearWp85 = 1;
		if (ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,90)) eleRecoIsoTree_probe_nearWp90 = 1;
		if (ElectronWP(iEvent->eleInfo[probeNumber],rhoForEleIso,95)) eleRecoIsoTree_probe_nearWp95 = 1;
	      }
              
	      eleRecoIsoTree_probe_passingEle27HLT = 0;
	      for (size_t i=0 ; i < ele27HLT.size(); i++) {
		if (iEvent->eleInfo[probeNumber].p4.DeltaR(ele27HLT[i]) < trigMatchCut) eleRecoIsoTree_probe_passingEle27HLT = 1;
	      }

              nJets = 0;
              nBJets = 0;
	      probe_dPhiJetMin = 99999.;
              probe_closestJetPt = 0.;
	      std::vector<TLorentzVector> dijets;
              for(size_t j=0; j< iEvent->simpleJets2.size(); j++) {
		if (iEvent->simpleJets2[j].p4.Pt() > 20. && fabs(iEvent->simpleJets2[j].p4.Eta()) < 2.5
                    && iEvent->eleInfo[probeNumber].p4.DeltaR(iEvent->simpleJets2[j].p4) > 0.5 && iEvent->eleInfo[tagNumber].p4.DeltaR(iEvent->simpleJets2[j].p4) > 0.5) {
		  nJets++;
                  if (iEvent->simpleJets2[j].csv>0.4) nBJets++;
                  if (iEvent->eleInfo[probeNumber].p4.DeltaPhi(iEvent->simpleJets2[j].p4) < probe_dPhiJetMin) {
		    probe_dPhiJetMin = iEvent->eleInfo[probeNumber].p4.DeltaPhi(iEvent->simpleJets2[j].p4);
                    probe_closestJetPt = iEvent->simpleJets2[j].p4.Pt();
                  }
                  for(size_t k=j+1; k < iEvent->simpleJets2.size(); k++) {
                    if (iEvent->simpleJets2[k].p4.Pt() > 20. && fabs(iEvent->simpleJets2[k].p4.Eta()) < 2.5
                        && iEvent->eleInfo[probeNumber].p4.DeltaR(iEvent->simpleJets2[k].p4) > 0.5 && iEvent->eleInfo[tagNumber].p4.DeltaR(iEvent->simpleJets2[k].p4) > 0.5) {
                      dijets.push_back(iEvent->simpleJets2[j].p4+iEvent->simpleJets2[k].p4);
                    }
                  }
                }
              }
              probe_dPhiDiJet = 99999.;
              probe_diJetPt = 0.;
              for (size_t j=0; j < dijets.size() ; j++) {
                if (dijets[j].Pt() > probe_diJetPt) {
                  probe_dPhiDiJet = iEvent->eleInfo[probeNumber].p4.DeltaPhi(dijets[j]);
                  probe_diJetPt = dijets[j].Pt();
                }
              }
              
	      eleRecoIsoTree->Fill(); // Only if trigger-matched tag and a probe
              
	    } // potential RecoIso probe size > 0
            
	  } // matched ele tag
	} // potential ele tag
        

	/*
        for(size_t m=0;m<iEvent->eleInfo.size();m++) {
	  if (verbose_) std::cout << "FOUND Ele " << iEvent->eleInfo[m].p4.Pt() <<  " " << EVENT.event << " " << candW->size() << " " << candZ->size() << std::endl;
        }
	cout << "Before PFMET" << endl;
	cout << "PFMET: " << iEvent->pfmet.p4.Et() << endl;
	*/

	if (verbose_) {
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
          interestingNames.push_back("HLT_IsoMu24_v");
          interestingNames.push_back("HLT_Mu40_v");
          interestingNames.push_back("HLT_IsoMu20_WCandPt80_v");

	  
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
	      std::cout << name << " " << (*trigPathRef)->wasRun() << " " << (*trigPathRef)->wasAccept() << std::endl;
	      const TriggerFilterRefVector trigFilterRefs( triggerEvent->pathFilters(name) );
	      for ( TriggerFilterRefVector::const_iterator trigFilterRef  = trigFilterRefs.begin(); trigFilterRef != trigFilterRefs.end(); ++trigFilterRef ) {
		std::string label = (*trigFilterRef)->label();
		std::cout << "  " << label << " " << (*trigFilterRef)->status() << std::endl;
		//      std::cout << "Looking for muon label " << muonLabels[interestingName] << " (and label is " << label << ")" << std::endl;
		//      if (muonLabels.count(name) && label == muonLabels[interestingName]) {
		//        std::cout << "Found label " << muonLabels[interestingName] << "(and label is "<< label << ")" << std::endl;
		//        const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
		//        for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		//          trigMuons[interestingName].push_back(GENPTOLORP(*trigObjRef));
		//          std::cout << "   Pushing muon back to " << interestingName << std::endl;
		//        }
		//      }
		const TriggerObjectRefVector trigObjRefs(triggerEvent->filterObjects( label ));
		for ( TriggerObjectRefVector::const_iterator trigObjRef = trigObjRefs.begin() ; trigObjRef != trigObjRefs.end() ; trigObjRef++) {
		  std::cout << "    " << (*trigObjRef)->pt() << " " << (*trigObjRef)->eta() << " " << (*trigObjRef)->phi() << std::endl;
		}
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
  muRecoIsoDir->cd() ; muRecoIsoTree->Write();  
  eleTrigDir->cd();  eleTrigTree->Write();
  eleRecoDir->cd();  eleRecoTree->Write();
    eleRecoIsoDir->cd() ; eleRecoIsoTree->Write();  
  _outFile->Write();
  _outFile->Close();
  return 0;
}


