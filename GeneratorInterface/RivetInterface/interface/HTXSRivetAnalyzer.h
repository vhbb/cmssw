#ifndef GeneratorInterface_RivetInterface_HTXSRivetAnalyzer
#define GeneratorInterface_RivetInterface_HTXSRivetAnalyzer

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "Rivet/AnalysisHandler.hh"

#include "GeneratorInterface/RivetInterface/src/HiggsTemplateCrossSections.cc"

#include <vector>
#include <string>

class HTXSRivetAnalyzer : public edm::EDAnalyzer
{

  public:

    explicit HTXSRivetAnalyzer(const edm::ParameterSet&);
    ~HTXSRivetAnalyzer();
  
  private:
    
    virtual void beginJob();
    virtual void endJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void beginRun(const edm::Run&, const edm::EventSetup&);
    virtual void endRun(const edm::Run&, const edm::EventSetup&);

    edm::EDGetTokenT<edm::HepMCProduct> _hepmcCollection;
    Rivet::AnalysisHandler* _analysisHandler;
    bool _isFirstEvent;
    Rivet::HiggsTemplateCrossSections* _HTXS;

};

#endif
