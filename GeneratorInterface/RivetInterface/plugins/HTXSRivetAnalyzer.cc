#include "GeneratorInterface/RivetInterface/interface/HTXSRivetAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

using namespace Rivet;
using namespace edm;

HTXSRivetAnalyzer::HTXSRivetAnalyzer(const edm::ParameterSet& pset) : 
_HTXS()
{
  _hepmcCollection = consumes<HepMCProduct>(pset.getParameter<edm::InputTag>("HepMCCollection"));
}

HTXSRivetAnalyzer::~HTXSRivetAnalyzer(){
}

void HTXSRivetAnalyzer::beginJob(){
}

void HTXSRivetAnalyzer::beginRun(const edm::Run& iRun,const edm::EventSetup& iSetup){
  return;
}

void HTXSRivetAnalyzer::analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup){
  
  //get the hepmc product from the event
  edm::Handle<HepMCProduct> evt;
  iEvent.getByToken(_hepmcCollection, evt);

  // get HepMC GenEvent
  const HepMC::GenEvent *myGenEvent = evt->GetEvent();
  
  std::cout<<"HTXSRivetAnalyzer begin classify"<<endl;

  // get the classification                                                                                                   
  //HiggsClassification cat = classifyEvent(*myGenEvent,m_HiggsProdMode);
  HiggsClassification cat = _HTXS.classifyEvent(*myGenEvent,HTXS::QQ2ZH);

  std::cout<<"HTXSRivetAnalyzer cat.prodMode "<<cat.prodMode<<std::endl;

}


void HTXSRivetAnalyzer::endRun(const edm::Run& iRun,const edm::EventSetup& iSetup){
  return;
}

void HTXSRivetAnalyzer::endJob(){
}

DEFINE_FWK_MODULE(HTXSRivetAnalyzer);
