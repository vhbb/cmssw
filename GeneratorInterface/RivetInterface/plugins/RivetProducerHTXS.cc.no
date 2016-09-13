#include "GeneratorInterface/RivetInterface/interface/RivetProducerHTXS.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Analysis.hh"

using namespace Rivet;
using namespace edm;

RivetProducerHTXS::RivetProducerHTXS(const edm::ParameterSet& pset) : 
_analysisHandler(),
_isFirstEvent(true)
{
  //retrive the analysis name from parameter set
  std::vector<std::string> analysisNames = pset.getParameter<std::vector<std::string> >("AnalysisNames");
  
  _hepmcCollection = consumes<HepMCProduct>(pset.getParameter<edm::InputTag>("HepMCCollection"));

  _useExternalWeight = pset.getParameter<bool>("UseExternalWeight");
  if (_useExternalWeight) {
    if (!pset.exists("GenEventInfoCollection")){
      throw cms::Exception("RivetProducerHTXS") << "when using an external event weight you have to specify the GenEventInfoProduct collection from which the weight has to be taken " ; 
    }
    _genEventInfoCollection = consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("GenEventInfoCollection"));
    _LHECollection          = consumes<LHEEventProduct>(pset.getParameter<edm::InputTag>("LHECollection"));
    _useLHEweights          = pset.getParameter<bool>("useLHEweights");
    _LHEweightNumber        = pset.getParameter<int>("LHEweightNumber");    
    
  }

  //get the analyses
  _analysisHandler.addAnalyses(analysisNames);

  //go through the analyses and check those that need the cross section
  const std::set< AnaHandle, CmpAnaHandle > & analyses = _analysisHandler.analyses();

  std::set< AnaHandle, CmpAnaHandle >::const_iterator ibeg = analyses.begin();
  std::set< AnaHandle, CmpAnaHandle >::const_iterator iend = analyses.end();
  std::set< AnaHandle, CmpAnaHandle >::const_iterator iana; 
  double xsection = -1.;
  xsection = pset.getParameter<double>("CrossSection");
  for (iana = ibeg; iana != iend; ++iana){
    if ((*iana)->needsCrossSection())
    (*iana)->setCrossSection(xsection);
  }

  produces<int>();

}

RivetProducerHTXS::~RivetProducerHTXS(){
}

void RivetProducerHTXS::beginJob(){
  //set the environment, very ugly but rivet is monolithic when it comes to paths
  char * cmsswbase    = getenv("CMSSW_BASE");
  char * cmsswrelease = getenv("CMSSW_RELEASE_BASE");
  std::string rivetref, rivetinfo;
  rivetref = "RIVET_REF_PATH=" + string(cmsswbase) + "/src/GeneratorInterface/RivetInterface/data:" + string(cmsswrelease) + "/src/GeneratorInterface/RivetInterface/data";
  rivetinfo = "RIVET_INFO_PATH=" + string(cmsswbase) + "/src/GeneratorInterface/RivetInterface/data:" + string(cmsswrelease) + "/src/GeneratorInterface/RivetInterface/data";
  putenv(strdup(rivetref.c_str()));
  putenv(strdup(rivetinfo.c_str()));
}

void RivetProducerHTXS::beginRun(const edm::Run& iRun,const edm::EventSetup& iSetup){
  return;
}

void RivetProducerHTXS::analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup){
  
  //get the hepmc product from the event
  edm::Handle<HepMCProduct> evt;
  iEvent.getByToken(_hepmcCollection, evt);

  // get HepMC GenEvent
  const HepMC::GenEvent *myGenEvent = evt->GetEvent();
  
  //if you want to use an external weight we have to clone the GenEvent and change the weight  
  if ( _useExternalWeight ){
    
    HepMC::GenEvent * tmpGenEvtPtr = new HepMC::GenEvent( *(evt->GetEvent()) );
    if (tmpGenEvtPtr->weights().size() == 0) {
      throw cms::Exception("RivetProducerHTXS") << "Original weight container has 0 size ";
    }
    if (tmpGenEvtPtr->weights().size() > 1) {
      edm::LogWarning("RivetProducerHTXS") << "Original event weight size is " << tmpGenEvtPtr->weights().size() << ". Will change only the first one ";  
    }
    
    if(!_useLHEweights){
      edm::Handle<GenEventInfoProduct> genEventInfoProduct;
      iEvent.getByToken(_genEventInfoCollection, genEventInfoProduct);
      tmpGenEvtPtr->weights()[0] = genEventInfoProduct->weight();
    }else{
      edm::Handle<LHEEventProduct> lheEventHandle;
      iEvent.getByToken(_LHECollection,lheEventHandle);
      const LHEEventProduct::WGT& wgt = lheEventHandle->weights().at(_LHEweightNumber);
      tmpGenEvtPtr->weights()[0] = wgt.wgt;
    }
    myGenEvent = tmpGenEvtPtr;

  }
  

  //apply the beams initialization on the first event
  if (_isFirstEvent){
    _analysisHandler.init(*myGenEvent);
    _isFirstEvent = false;
  }

  //run the analysis
  _analysisHandler.analyze(*myGenEvent);

  //if we have cloned the GenEvent, we delete it
  if ( _useExternalWeight ) 
  delete myGenEvent;
}


void RivetProducerHTXS::endRun(const edm::Run& iRun,const edm::EventSetup& iSetup){

  return;
}



//from Rivet 2.X: Analysis.hh (cls 18Feb2014)
/// List of registered analysis data objects
//const vector<AnalysisObjectPtr>& analysisObjects() const {
//return _analysisobjects;
//}



void RivetProducerHTXS::endJob(){
}


void RivetProducerHTXS::normalizeTree()    {
  using namespace YODA;
  std::vector<string> analyses = _analysisHandler.analysisNames();
  
}



DEFINE_FWK_MODULE(RivetProducerHTXS);
