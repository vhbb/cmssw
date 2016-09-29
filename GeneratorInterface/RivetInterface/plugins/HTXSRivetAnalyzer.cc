#include "GeneratorInterface/RivetInterface/interface/HTXSRivetAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

using namespace Rivet;
using namespace edm;

HTXSRivetAnalyzer::HTXSRivetAnalyzer(const edm::ParameterSet& pset) : 
_analysisHandler(),
_isFirstEvent(true),
_HTXS()
{
  _hepmcCollection = consumes<HepMCProduct>(pset.getParameter<edm::InputTag>("HepMCCollection"));
  std::vector<std::string> analysisNames = pset.getParameter<std::vector<std::string> >("AnalysisNames");
}

HTXSRivetAnalyzer::~HTXSRivetAnalyzer(){
    delete &_HTXS;
}

void HTXSRivetAnalyzer::beginJob(){
}

void HTXSRivetAnalyzer::beginRun(const edm::Run& iRun,const edm::EventSetup& iSetup){
  _analysisHandler.addAnalysis(&_HTXS);
  return;
}

void HTXSRivetAnalyzer::analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup){
  
  //get the hepmc product from the event
  edm::Handle<HepMCProduct> evt;
  iEvent.getByToken(_hepmcCollection, evt);

  // get HepMC GenEvent
  const HepMC::GenEvent *myGenEvent = evt->GetEvent();

  if (_isFirstEvent){
      _HTXS.setHiggsProdMode(HTXS::QQ2ZH);
     _analysisHandler.init(*myGenEvent);
     _isFirstEvent = false;
  }

  // Run the analyses
  HiggsClassification cat = _HTXS.classifyEvent(*myGenEvent,HTXS::QQ2ZH);

  std::cout<<"HTXSRivetAnalyzer cat.prodMode "<<cat.prodMode <<std::endl;
  std::cout<<"HTXSRivetAnalyzer cat.errorCode "<< cat.errorCode <<std::endl;
  std::cout<<"HTXSRivetAnalyzer cat.stage0_cat "<< cat.stage0_cat <<std::endl; 
  std::cout<<"HTXSRivetAnalyzer cat.stage1_cat_pTjet30GeV "<< cat.stage1_cat_pTjet30GeV <<std::endl;
  std::cout<<"HTXSRivetAnalyzer cat.jets30 "<< cat.jets30 <<std::endl;
  
}


void HTXSRivetAnalyzer::endRun(const edm::Run& iRun,const edm::EventSetup& iSetup){
    _HTXS.printClassificationSummary();
  return;
}

void HTXSRivetAnalyzer::endJob(){
}

DEFINE_FWK_MODULE(HTXSRivetAnalyzer);
