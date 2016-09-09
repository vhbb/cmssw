#ifndef GeneratorInterface_RivetInterface_RivetProducerHTXS
#define GeneratorInterface_RivetInterface_RivetProducerHTXS

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "Rivet/AnalysisHandler.hh"

//DQM services
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "Rivet/Tools/RivetYODA.hh"
//#include "YODA/ROOTCnv.h"

#include <iostream>
#include <vector>
#include <string>
#include <map>



#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"


class RivetProducerHTXS : public edm::EDProducer
{
  public:
  RivetProducerHTXS(const edm::ParameterSet&);

  virtual ~RivetProducerHTXS();

  virtual void beginJob();

  virtual void endJob();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void beginRun(const edm::Run&, const edm::EventSetup&);

  virtual void endRun(const edm::Run&, const edm::EventSetup&);
  
  private:

  void normalizeTree();

  edm::EDGetTokenT<edm::HepMCProduct> _hepmcCollection;
  bool                     _useExternalWeight;
  bool                     _useLHEweights;
  int                      _LHEweightNumber;
  edm::EDGetTokenT<LHEEventProduct> _LHECollection;
  edm::EDGetTokenT<GenEventInfoProduct> _genEventInfoCollection;
  Rivet::AnalysisHandler   _analysisHandler;   
  bool                     _isFirstEvent;
  std::string              _outFileName;
  bool                     _doFinalize;
  bool                     _produceDQM;

  void produce( edm::Event &, const edm::EventSetup & );

  DQMStore *dbe;
  std::vector<MonitorElement *> _mes;
};


/*
template<class YODATYPE, class ROOTTYPE> 
ROOTTYPE* 
  RivetProducerHTXS::prebook(const YODATYPE* yodah, const std::string& name){
  ROOTTYPE* h = 0;
  if (yodah->axis().isFixedBinning() ) {//equidistant binning (easier case)
    int nbins = yodah->axis().bins();
    h = new ROOTTYPE(name.c_str(), name.c_str(), nbins, yodah->axis().lowerEdge(), yodah->axis().upperEdge());
  } else {
    int nbins = yodah->axis().bins();
    const YODA::Axis1D* vax = dynamic_cast<const YODA::Axis1D*>(&yodah->axis());
    if (! vax ){
      throw cms::Exception("RivetProducerHTXS") << "Cannot dynamix cast an YODA axis to VariAxis ";
    }
    double* bins = new double[nbins+1];
    for (int i=0; i<nbins; ++i) {
      bins[i] = vax->binEdges(i).first;
    }
    bins[nbins] = vax->binEdges(nbins-1).second; //take last bin right border
    h = new ROOTTYPE(name.c_str(), name.c_str(), nbins, bins);
    delete bins;
  }
  return h; 
}

template<> 
TH1F* RivetProducerHTXS::yoda2root(const YODA::IHistogram1D* yodah, const std::string& name){
  TH1F* h = prebook<YODA::Histo1D, TH1F>(yodah, name);
  for (int i = 0; i < yodah->axis().bins(); ++i){
    h->SetBinContent(i+1, yodah->binHeight(i));
    h->SetBinError(i+1, yodah->binError(i));
  }  
  return h;
}

template<>
TProfile* RivetProducerHTXS::yoda2root(const YODA::IProfile1D* yodah, const std::string& name){
  TProfile* h = prebook<YODA::IProfile1D, TProfile>(yodah, name);
  for (int i = 0; i < yodah->axis().bins(); ++i){
    h->SetBinContent(i+1, yodah->binMean(i));
    h->SetBinError(i+1, yodah->binRms(i));
  }
  return h;
}
*/

#endif
