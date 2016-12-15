/* \class PDFWeightsProducer
 *
 * \author Luca Lista, INFN
 *
 * Convert PDF replicas in Hessian representation
 * Originally developed by Josh BEndavid, adapted by Luca Perrozzi
 *
 *
 */

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "PhysicsTools/HepMCCandAlgos/interface/PDFWeightsHelper.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//#include "SimDataFormats/HiGenData/interface/SubEventMap.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/transform.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <fstream>
#include <algorithm>

#include <vector>
#include <string>
#include <map>
#include <set>


using namespace edm;
// using namespace reco;
using namespace std;
using namespace HepMC;

namespace edm { class ParameterSet; }
namespace HepMC { class GenParticle; class GenEvent; }

class PDFWeightsProducer : public edm::EDProducer {
 public:
  /// constructor
  PDFWeightsProducer( const edm::ParameterSet & );
  /// destructor
  ~PDFWeightsProducer();

  /// process one event
  virtual void produce( edm::Event& e, const edm::EventSetup&) override;

 private:
  PDFWeightsHelper pdfweightshelper_;
  EDGetTokenT<LHEEventProduct> srcToken_;
  EDGetTokenT<LHEEventProduct> srcTokenAlt_;
  EDGetTokenT<GenEventInfoProduct> srcTokenGen_;

  unsigned int pdfWeightOffset_;
  unsigned int nPdfWeights_;
  unsigned int nPdfEigWeights_;
  
  std::vector<float> pdfweights_;
  // std::vector<float> pdfeigweights_;

  float weight_;

};



PDFWeightsProducer::PDFWeightsProducer( const ParameterSet & cfg ) :
  srcToken_( consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"))),
  srcTokenAlt_( consumes<LHEEventProduct>(edm::InputTag("source"))),
  srcTokenGen_( consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  pdfWeightOffset_(cfg.getParameter<unsigned int>("pdfWeightOffset")),
  nPdfWeights_(cfg.getParameter<unsigned int>("nPdfWeights")),
  nPdfEigWeights_(cfg.getParameter<unsigned int>("nPdfEigWeights")),
  pdfweights_(nPdfWeights_)
  // pdfeigweights_(nPdfEigWeights_)
{
  
  produces<std::vector<float>>();

  edm::FileInPath mc2hessianCSV = cfg.getParameter<edm::FileInPath>("mc2hessianCSV");
  pdfweightshelper_.Init(nPdfWeights_,nPdfEigWeights_,mc2hessianCSV);
}

PDFWeightsProducer::~PDFWeightsProducer() {
}

void PDFWeightsProducer::produce( Event& evt, const EventSetup& es ) {

    auto_ptr<std::vector<float>> pdfeigweights_ptr( new std::vector<float>() );
    // pdfeigweights_ = * pdfeigweights_ptr;
    Handle<LHEEventProduct> lheInfo;
    evt.getByToken( srcToken_, lheInfo );
    
    if (!lheInfo.isValid()) {
      evt.getByToken( srcTokenAlt_, lheInfo );
    }

    double nomlheweight = lheInfo->weights()[0].wgt;

    Handle<GenEventInfoProduct> genInfo;
    evt.getByToken( srcTokenGen_, genInfo );
    
    weight_ = genInfo->weight();
    
    //get the original mc replica weights
    std::vector<double> inpdfweights(nPdfWeights_);
    for (unsigned int ipdf=0; ipdf<nPdfWeights_; ++ipdf) {
      unsigned int iwgt = ipdf + pdfWeightOffset_;
      
      //this is the weight to be used for evaluating uncertainties with mc replica weights
      pdfweights_[ipdf] = lheInfo->weights()[iwgt].wgt*weight_/nomlheweight;
      
      //this is the raw weight to be fed to the mc2hessian convertor
      inpdfweights[ipdf] = lheInfo->weights()[iwgt].wgt;
      
    }
    
    std::vector<double> outpdfweights(nPdfEigWeights_);
    //do the actual conversion, where the nominal lhe weight is needed as the reference point for the linearization
    pdfweightshelper_.DoMC2Hessian(nomlheweight,inpdfweights.data(),outpdfweights.data());
    
    for (unsigned int iwgt=0; iwgt<nPdfEigWeights_; ++iwgt) {
      double wgtval = outpdfweights[iwgt];
      
      //the is the weight to be used for evaluating uncertainties with hessian weights
      // pdfeigweights_[iwgt] = wgtval*weight_/nomlheweight;
      pdfeigweights_ptr->push_back(wgtval*weight_/nomlheweight);
    }    
  
  evt.put( pdfeigweights_ptr );

}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( PDFWeightsProducer );

