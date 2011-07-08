#ifndef HBBCANDIDATEFINDERALGO_HH
#define HBBCANDIDATEFINDERALGO_HH

#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"

class HbbCandidateFinderAlgo {
 public:

  explicit HbbCandidateFinderAlgo(bool verbose, float jetPt): verbose_(verbose), jetPtThreshold(jetPt){}


  void run (const VHbbEvent*, std::vector<VHbbCandidate>  &);
  float getDeltaTheta( const VHbbEvent::SimpleJet & j1, const VHbbEvent::SimpleJet & j2 ) const ;

  
 protected:
  
 
  bool  findDiJets (const std::vector<VHbbEvent::SimpleJet>& , VHbbEvent::SimpleJet& , VHbbEvent::SimpleJet& ,std::vector<VHbbEvent::SimpleJet>& );
  
  void findMuons (const std::vector<VHbbEvent::MuonInfo>& muons, std::vector<VHbbEvent::MuonInfo>& out);
  void findElectrons(const std::vector<VHbbEvent::ElectronInfo>& electrons, std::vector<VHbbEvent::ElectronInfo>& out);
  void findMET(const VHbbEvent::METInfo& met, std::vector<VHbbEvent::METInfo>& out);
  
 private:
  bool verbose_;
  float jetPtThreshold;
  


};

#endif
