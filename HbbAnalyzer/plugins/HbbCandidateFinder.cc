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
  
  run(vHbbEvent.product(), vHbbCandidates);


  iEvent.put(vHbbCandidates);  

}

void HbbCandidateFinder::run (const VHbbEvent* event, std::auto_ptr<std::vector<VHbbCandidate> > & candidates){
  //
  // first find the jets
  //
  std::pair<int,int> jets = findDiJets(event->simpleJets2);
  
  if (jets.first<0 || jets.second<0) return ;

  //
  // search for a dilepton - just 
  //
  int ele =    findDiElectron(event->diElectronInfo);
  int mu =  findDiMuon(event->diMuonInfo);

  if (ele<0 && mu < 0 ) return;

  //
  // fill!
  //
  VHbbCandidate temp;
  temp.H.jets.push_back(event->simpleJets2[jets.first]);
  temp.H.jets.push_back(event->simpleJets2[jets.second]);
  temp.H.fourMomentum = (event->simpleJets2[jets.first]).fourMomentum+(event->simpleJets2[jets.second]).fourMomentum;
  if (mu>-1){
    temp.V.muons.push_back((event->diMuonInfo[mu]).daughter1);
    temp.V.muons.push_back((event->diMuonInfo[mu]).daughter2);
    temp.V.fourMomentum  = (temp.V.muons[0]).fourMomentum+(temp.V.muons[1]).fourMomentum;
  }else{
    temp.V.electrons.push_back((event->diElectronInfo[ele]).daughter1);
    temp.V.electrons.push_back((event->diElectronInfo[ele]).daughter2);
    temp.V.fourMomentum  = (temp.V.electrons[0]).fourMomentum+(temp.V.electrons[1]).fourMomentum;
  }

  candidates->push_back(temp);

}

std::pair <int, int>  HbbCandidateFinder::findDiJets (const std::vector<VHbbEvent::SimpleJet>& jets){
  
  //
  // select jets
  //
  //

  unsigned int maxJets = jets.size();
  int i1(-99),i2(-99);
  float btag_max(-100);
  //
  // do the combinatorics

  if (maxJets<2) return pair<int,int>(-1,-1);

  for (unsigned int j1=0; j1<maxJets-1; j1++) {
    for (unsigned int j2=j1+1; j2<maxJets; j2++) {
      float j1_btag = (jets)[j1].csv;
      float j2_btag = (jets)[j2].csv;
      
      if (j1_btag<=0.0 || j2_btag<=0.0) continue;

      if ((jets)[j1].fourMomentum.Pt()< 30. ||
	  (jets)[j2].fourMomentum.Pt()< 30.) continue;

      if (j1_btag+j2_btag>btag_max) { 
	btag_max = j1_btag+j2_btag; 
	i1 = j1; 
	i2 = j2; 
      }
    }
  }
  if (i1>=0 && i2>=0) return pair<int,int>(i1,i2);
  else return pair<int,int>(-1,-1);
}

int HbbCandidateFinder::findDiMuon (const std::vector<VHbbEvent::DiMuonInfo>& dimu){
  int res = -99;
  if (dimu.size()>0) res= 1;
  return res;

}
int HbbCandidateFinder::findDiElectron (const std::vector<VHbbEvent::DiElectronInfo>& die){
  int res = -99;
  if (die.size()>0) res= 1;
  return res;
}



//define this as a plug-in
DEFINE_FWK_MODULE(HbbCandidateFinder);
