// system include files
#include <memory>
#include <iostream>
using namespace std;

// user include files
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TH1.h"
#include "TTree.h"
#include "TMath.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include <TString.h>

#include "DataFormats/GeometryVector/interface/Phi.h"

#include<TVector2.h>

#include "TArrayD.h"
#include "TLorentzVector.h"

#include "DataFormats/METReco/interface/PFMET.h"

//from .cc
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEventAuxInfo.h"

#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "RecoBTag/Records/interface/BTagPerformanceRecord.h"
#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"
#include "RecoBTag/PerformanceDB/interface/BtagPerformance.h"

//Include files needed for CSV Variables
#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputer.h"


// class decleration
//
class JetCorrectionUncertainty;


class HbbAnalyzerAlgo {
 
 struct BTagSFContainer{
   const BtagPerformance * BTAGSF_CSVL;
   const BtagPerformance * BTAGSF_CSVM;
   const BtagPerformance * BTAGSF_CSVT;
   const BtagPerformance * MISTAGSF_CSVL;
   const BtagPerformance * MISTAGSF_CSVM;
   const BtagPerformance * MISTAGSF_CSVT;
};

 public:
  explicit HbbAnalyzerAlgo();
  ~HbbAnalyzerAlgo();
  virtual void produce( fwlite::Event &,VHbbEvent&, VHbbEventAuxInfo& );
  
 protected:
  TVector2 getTvect( const pat::Jet* patJet );

  TLorentzVector getChargedTracksMomentum(const pat::Jet* patJet );
  
 private:
  virtual void beginJob() ;
  virtual void endJob() ;
  virtual void fillMuBlock(std::vector<pat::Muon>::const_iterator mu, int *muInfo);
  virtual void fillScaleFactors(VHbbEvent::SimpleJet&, BTagSFContainer);
  
  // ----------member data ---------------------------
  
  edm::InputTag eleLabel_;
  edm::InputTag muoLabel_;
  double lep_ptCutForBjets_;
  //  edm::InputTag elenoCutsLabel_;
  //  edm::InputTag muonoCutsLabel_;
  edm::InputTag jetLabel_;
  //  edm::InputTag subjetLabel_;
  //  edm::InputTag filterjetLabel_;
  //  edm::InputTag simplejet1Label_;
  edm::InputTag simplejet2Label_;
  //  edm::InputTag simplejet3Label_;
  //  edm::InputTag simplejet4Label_;
  edm::InputTag tauLabel_;
  edm::InputTag metLabel_;
  edm::InputTag phoLabel_;
  edm::InputTag hltResults_;

  bool runOnMC_;
  
  //  TMatrixD *pointerPt;
  TMatrixD *pointerEta;
  TMatrixD *pointerPhi;
  
  //The computer for the CSV variables
  const GenericMVAJetTagComputer *computer;

  bool verbose_;
 protected:
  void fillSimpleJet (VHbbEvent::SimpleJet& sj, std::vector<pat::Jet>::const_iterator iter);
  void setJecUnc(VHbbEvent::SimpleJet& sj,JetCorrectionUncertainty* jecunc);
  float metSignificance(const reco::MET * met);
};


