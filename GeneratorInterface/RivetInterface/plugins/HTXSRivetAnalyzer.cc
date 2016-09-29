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

  std::cout<<"HTXSRivetAnalyzer before declaring myGenEvent"<<endl;
  // get HepMC GenEvent
  const HepMC::GenEvent *myGenEvent = evt->GetEvent();
  std::cout<<"HTXSRivetAnalyzer after declaring myGenEvent"<<endl;

  /* This works 
  if (_isFirstEvent){
     _analysisHandler.init(*myGenEvent);
     _isFirstEvent = false;
  }
  _HTXS.setHiggsProdMode(HTXS::QQ2ZH);
  _analysisHandler.analyze(*myGenEvent);
  */


  /// Get the collection of currently registered analyses.
  const std::set< AnaHandle, CmpAnaHandle > & analyses = _analysisHandler.analyses();
  std::set< AnaHandle, CmpAnaHandle >::const_iterator ibeg = analyses.begin();
  std::set< AnaHandle, CmpAnaHandle >::const_iterator iend = analyses.end();
  std::set< AnaHandle, CmpAnaHandle >::const_iterator iana;
  for (iana = ibeg; iana != iend; ++iana){
      std::cout<<(*iana)->name()<<std::endl;
      //(*iana)->_allowProjReg = true;
  }
  std::cout<<"HTXSRivetAnalyzer before first event"<<endl;

  if (_isFirstEvent){
      //_analysisHandler.addAnalysis(&_HTXS);
      _HTXS.setHiggsProdMode(HTXS::QQ2ZH);
      std::cout<<"HTXSRivetAnalyzer after setHiggsProdMode"<<endl;
     _analysisHandler.init(*myGenEvent);
      std::cout<<"HTXSRivetAnalyzer after _analysisHandler.init"<<endl;
     _isFirstEvent = false;
  }
  std::cout<<"HTXSRivetAnalyzer after first event"<<endl;

  // Run the analyses
  std::cout<<"HTXSRivetAnalyzer begin classify"<<endl;
  
  HiggsClassification cat = _HTXS.classifyEvent(*myGenEvent,HTXS::QQ2ZH);


  //const AnaHandle analysis(const std::string& analysisname) const;
  //const AnaHandle & htxsana = _analysisHandler.analysis("HiggsTemplateCrossSections");
  //htxsana.analyze(*myGenEvent);

  std::cout<<"HTXSRivetAnalyzer cat.prodMode "<<cat.prodMode <<std::endl;
  std::cout<<"HTXSRivetAnalyzer cat.errorCode "<< cat.errorCode <<std::endl;
  std::cout<<"HTXSRivetAnalyzer cat.stage0_cat "<< cat.stage0_cat <<std::endl; 
  std::cout<<"HTXSRivetAnalyzer cat.stage1_cat_pTjet25GeV "<< cat.stage1_cat_pTjet25GeV <<std::endl;
  std::cout<<"HTXSRivetAnalyzer cat.stage1_cat_pTjet30GeV "<< cat.stage1_cat_pTjet30GeV <<std::endl;
  
}


void HTXSRivetAnalyzer::endRun(const edm::Run& iRun,const edm::EventSetup& iSetup){
    _HTXS.printClassificationSummary();
  return;
}

void HTXSRivetAnalyzer::endJob(){
}

DEFINE_FWK_MODULE(HTXSRivetAnalyzer);
