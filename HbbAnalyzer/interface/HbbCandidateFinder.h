// system include files
#include <memory>
#include <iostream>
using namespace std;

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

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

#include "VHbbAnalysis/HbbAnalyzer/interface/VHbbCandidate.h"
#include "VHbbAnalysis/HbbAnalyzer/interface/VHbbEvent.h"




//
// class declaration
//



class HbbCandidateFinder : public edm::EDProducer {
  
 public:
  explicit HbbCandidateFinder(const edm::ParameterSet&);
  ~HbbCandidateFinder();
  void produce( edm::Event&, const edm::EventSetup& );

  float getDeltaTheta( VHbbEvent::SimpleJet * j1, VHbbEvent::SimpleJet * j2 );


protected:

  void run (const VHbbEvent*, std::auto_ptr<std::vector<VHbbCandidate> > &);
  
  std::pair <int, int>  findDiJets (const std::vector<VHbbEvent::SimpleJet>& jets);


  int findDiMuon (const std::vector<VHbbEvent::DiMuonInfo>&);
  int findDiElectron (const std::vector<VHbbEvent::DiElectronInfo>&);


 private:
  virtual void beginJob() ;
  virtual void endJob() ;
  
  edm::InputTag vhbbevent_;

  
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//

