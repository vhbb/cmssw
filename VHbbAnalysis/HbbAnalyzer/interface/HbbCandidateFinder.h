// system include files
#include <memory>
#include <iostream>
using namespace std;

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TH1.h"
#include "TTree.h"
#include "TMath.h"

#include <TString.h>


#include "VHbbAnalysis/VHbbDataFormats/interface/HbbCandidateFinderAlgo.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"

//
// class declaration
//



class HbbCandidateFinder : public edm::EDFilter {
  
 public:
  explicit HbbCandidateFinder(const edm::ParameterSet&);
  ~HbbCandidateFinder();
  bool filter( edm::Event&, const edm::EventSetup& );

  float getDeltaTheta( const VHbbEvent::SimpleJet & j1, const VHbbEvent::SimpleJet & j2 ) const ;


protected:

  void run (const VHbbEvent*, std::auto_ptr<std::vector<VHbbCandidate> > &);
  
 bool  findDiJets (const std::vector<VHbbEvent::SimpleJet>& , VHbbEvent::SimpleJet& , VHbbEvent::SimpleJet& ,std::vector<VHbbEvent::SimpleJet>& );

 // voif findVectorCandidate()


 void findMuons (const std::vector<VHbbEvent::MuonInfo>& muons, std::vector<VHbbEvent::MuonInfo>& out, std::vector<unsigned int>&);
 void findElectrons(const std::vector<VHbbEvent::ElectronInfo>& electrons, std::vector<VHbbEvent::ElectronInfo>& out, std::vector<unsigned int>&);
 void findMET(const VHbbEvent::METInfo& met, std::vector<VHbbEvent::METInfo>& out);

 private:
  virtual void beginJob() ;
  virtual void endJob() ;
  
  edm::InputTag vhbbevent_;
 HbbCandidateFinderAlgo *algo_;
 bool verbose_;
 bool useHighestHiggs;
 bool applyFilter;
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//

