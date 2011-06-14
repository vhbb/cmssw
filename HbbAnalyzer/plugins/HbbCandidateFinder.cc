#include "VHbbAnalysis/HbbAnalyzer/interface/HbbCandidateFinder.h"




HbbCandidateFinder::HbbCandidateFinder(const edm::ParameterSet& iConfig):   vhbbevent_(iConfig.getParameter<edm::InputTag>("VHbbEventLabel")) {
  produces<std::vector<VHbbCandidate > >();
}

HbbCandidateFinder::~HbbCandidateFinder(){}

void HbbCandidateFinder::beginJob(){}
void HbbCandidateFinder::endJob(){}


float HbbCandidateFinder::getDeltaTheta( VHbbEvent::SimpleJet * j1, VHbbEvent::SimpleJet * j2 ){return -1.;}



void HbbCandidateFinder::produce( edm::Event& iEvent, const edm::EventSetup& iEventSetup){
  
  std::auto_ptr<std::vector<VHbbCandidate> >  vHbbCandidates( new std::vector<VHbbCandidate>  );

  edm::Handle<VHbbEvent>  vHbbEvent;
  iEvent.getByLabel(vhbbevent_, vHbbEvent);
  

  //
  // start searching for candidates
  //

  //  hbbCandidateFinderAlgo(vHbbCandidates, vHbbEvent-> result());
  // do nothing for a test
  
  iEvent.put(vHbbCandidates);  

}


//define this as a plug-in
DEFINE_FWK_MODULE(HbbCandidateFinder);
