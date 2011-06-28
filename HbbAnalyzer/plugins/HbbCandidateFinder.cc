#include "VHbbAnalysis/HbbAnalyzer/interface/HbbCandidateFinder.h"

HbbCandidateFinder::HbbCandidateFinder(const edm::ParameterSet& iConfig): verbose_(iConfig.getParameter<bool>("verbose") ){
  algo_ = new HbbCandidateFinderAlgo(iConfig.getParameter<bool>("verbose"), iConfig.getParameter<double>("jetPtThreshold"));
  produces<std::vector<VHbbCandidate > >();
}

HbbCandidateFinder::~HbbCandidateFinder(){delete algo_;}

void HbbCandidateFinder::beginJob(){}
void HbbCandidateFinder::endJob(){}

void HbbCandidateFinder::produce( edm::Event& iEvent, const edm::EventSetup& iEventSetup){
  
  std::auto_ptr<std::vector<VHbbCandidate> >  vHbbCandidates( new std::vector<VHbbCandidate>  );

  edm::Handle<VHbbEvent>  vHbbEvent;
  //  iEvent.getByLabel(vhbbevent_, vHbbEvent);
  iEvent.getByType(vHbbEvent);
  

  //
  // start searching for candidates
  //

  //  hbbCandidateFinderAlgo(vHbbCandidates, vHbbEvent-> result());
  // do nothing for a test
  
  algo_->run(vHbbEvent.product(),(*vHbbCandidates));
  
  
  if (verbose_)
    std::cout <<" Pushing VHbb candidates: "<<vHbbCandidates->size()<<std::endl;
  
  iEvent.put(vHbbCandidates);  
  
}


//define this as a plug-in
DEFINE_FWK_MODULE(HbbCandidateFinder);
