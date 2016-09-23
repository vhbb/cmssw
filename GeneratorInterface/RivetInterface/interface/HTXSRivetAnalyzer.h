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
  HTXSRivetAnalyzer(const edm::ParameterSet&);

  virtual ~HTXSRivetAnalyzer();

  virtual void beginJob() override;

  virtual void endJob() override;

  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;

  virtual void endRun(const edm::Run&, const edm::EventSetup&) override;
  
  private:

  edm::EDGetTokenT<edm::HepMCProduct> _hepmcCollection;
  Rivet::AnalysisHandler   _analysisHandler;
  bool                     _isFirstEvent;
  Rivet::HiggsTemplateCrossSections _HTXS;

};

#endif
