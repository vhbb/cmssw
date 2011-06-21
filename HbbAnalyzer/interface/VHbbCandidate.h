#ifndef VHbbCandidate__H 
#define VHbbCandidate__H 

#include <TLorentzVector.h>
#include <TVector2.h>
#include <vector>

#include "VHbbAnalysis/HbbAnalyzer/interface/VHbbEvent.h"

class VHbbCandidate {
 public:
  class VectorCandidate {
  public:

    TLorentzVector fourMomentum;
    std::vector<VHbbEvent::MuonInfo> muons;
    std::vector<VHbbEvent::ElectronInfo> electrons;
    std::vector<VHbbEvent::TauInfo> taus;
    std::vector<VHbbEvent::METInfo> mets;
    
  };
  
  class HiggsCandidate {
  public:
    TLorentzVector fourMomentum;
    std::vector<VHbbEvent::SimpleJet> jets;
    float deltaTheta;
  public:
    VHbbEvent::SimpleJet& firstJet(){return jets[0];}
    VHbbEvent::SimpleJet& secondJet(){return jets[1];}
 };
  
  
 public:
  HiggsCandidate H;
  VectorCandidate V;
  std::vector<VHbbEvent::SimpleJet> additionalJets;
};


#endif
