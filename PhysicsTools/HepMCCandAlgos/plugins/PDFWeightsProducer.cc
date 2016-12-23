/* \class PDFWeightsProducer
*
* \author Luca Perrozzi, ETHZ
*
* Convert PDF replicas in Hessian representation
* Originally developed by Josh BEndavid, adapted by Luca Perrozzi
*
*
*/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Run.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "PhysicsTools/HepMCCandAlgos/interface/PDFWeightsHelper.h"

#include <vector>
#include <stdio.h>
#include <string.h>
#include <TString.h>

using namespace edm;
using namespace std;
using namespace lhef;

class PDFWeightsProducer : public edm::EDProducer {
public:
  
  explicit PDFWeightsProducer( const ParameterSet & cfg ) :
  _lheRunInfo(consumes<LHERunInfoProduct,edm::InRun>(edm::InputTag("externalLHEProducer"))),
  srcTokenlheEvent_( consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"))),
  srcTokenGen_( consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  pdfWeightLHAnumber_(cfg.getParameter<unsigned int>("pdfWeightLHAnumber")),
  nPdfWeights_(cfg.getParameter<unsigned int>("nPdfWeights")),
  nPdfEigWeights_(cfg.getParameter<unsigned int>("nPdfEigWeights"))
  {
    produces<std::vector<float>>("outputHessianWeights");
    produces<unsigned int>("inputPDFset");
    mc2hessianCSV = cfg.getParameter<edm::FileInPath>("mc2hessianCSV");
  }

  ~PDFWeightsProducer();
  
private:
  
  virtual void beginRun(edm::Run const& run, edm::EventSetup const& es);
  virtual void endRun(edm::Run const& run, edm::EventSetup const& es);
  virtual void produce( edm::Event&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();
  
  edm::FileInPath mc2hessianCSV;
  PDFWeightsHelper pdfweightshelper_;
  EDGetTokenT<LHERunInfoProduct> _lheRunInfo;
  EDGetTokenT<LHEEventProduct> srcTokenlheEvent_;
  EDGetTokenT<GenEventInfoProduct> srcTokenGen_;
  
  unsigned int pdfWeightLHAnumber_;
  unsigned int nPdfWeights_;
  unsigned int nPdfEigWeights_;
  
  // std::vector<float> pdfweights_;
  Handle<LHERunInfoProduct> run;

  float weight_;
  std::vector<unsigned int> PDFweightsLHEorder_;
  std::vector<unsigned int> PDFweightsLHEsets_;
  
};

PDFWeightsProducer::~PDFWeightsProducer(){
}

void PDFWeightsProducer::beginJob(){
}

void PDFWeightsProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& es) {

  // Access LHERunInfoProduct (i.e. LHE header) to retrieve the list of weights in the sample
  iRun.getByLabel( edm::InputTag("externalLHEProducer"), run );

  typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
  LHERunInfoProduct myLHERunInfoProduct = *(run.product());
  for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
    // cout << "iter->tag()" << iter->tag() << endl;
    std::vector<std::string> lines = iter->lines();
    for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
      std::string line=lines.at(iLine);
      TString line_tstr = line.c_str();
      TString buffer;
      if (line_tstr.Contains("<weight id=")) {
        // cout<<iLine<< " "<<line<<endl;
        line_tstr = line_tstr.ReplaceAll("\n","").ReplaceAll("<weight id=\"","").ReplaceAll("</weight>","");
        TObjArray *tx = line_tstr.Tokenize("> ");
        line_tstr = ((TObjString *)(tx->At(tx->GetEntries()-1)))->String();
        unsigned int pdfnumber = line_tstr.ReplaceAll("PDF set = ","").Atoi();
        PDFweightsLHEorder_.push_back(pdfnumber);
        if(pdfnumber%100==1){
          cout << "PDF SET CANDIDATE FOR REWEIGHTING IN THIS SAMPLE= " << (pdfnumber-1) << endl;
          PDFweightsLHEsets_.push_back(pdfnumber-1);
        }
      }
    }
  }

  // In case this is required, check whether the sample contains the weights 
  // for the appropriate PDF set to be converted to Hessian.
  // A check is performed on the NLO NNPDF samples to catch the _pdfas PDF set in case is used
  if(pdfWeightLHAnumber_<1 || nPdfWeights_<1 || nPdfEigWeights_<1){
    TString mc2hessianCSV_str = mc2hessianCSV.fullPath();
    cout << "Either pdfWeightLHAnumber_<1 or nPdfWeights_<1 or nPdfEigWeights_<1: using predefined settings for" << endl;
    cout << "mc2hessianCSV " << mc2hessianCSV_str << endl;
    if(mc2hessianCSV_str.Contains("NNPDF30_nlo_as_0118_hessian_60.csv")){
      if(std::find(PDFweightsLHEsets_.begin(), PDFweightsLHEsets_.end(), 260000) != PDFweightsLHEsets_.end()){
        // 260000	NNPDF30_nlo_as_0118	      101
        pdfWeightLHAnumber_ = 260000;
      }else{
        // 292200	NNPDF30_nlo_nf_5_pdfas	103	
        pdfWeightLHAnumber_ = 292200;
      }
      nPdfWeights_        = 100;
      nPdfEigWeights_     = 60;
      
    }else if(mc2hessianCSV_str.Contains("NNPDF30_nlo_as_0118_nf_4_hessian_60.csv")){
      if(std::find(PDFweightsLHEsets_.begin(), PDFweightsLHEsets_.end(), 260400) != PDFweightsLHEsets_.end()){
        // 260400	NNPDF30_nlo_as_0118_nf_4	101	
        pdfWeightLHAnumber_ = 260400;
      }else{
        // 292000	NNPDF30_nlo_nf_4_pdfas	103	
        pdfWeightLHAnumber_ = 292000;
      }
      nPdfWeights_        = 100;
      nPdfEigWeights_     = 60;
      
    }else if(mc2hessianCSV_str.Contains("NNPDF30_lo_as_0130_nf_4_hessian_60.csv")){
      // 263400	NNPDF30_lo_as_0130_nf_4	  101	
      pdfWeightLHAnumber_ = 263400;
      nPdfWeights_        = 100;
      nPdfEigWeights_     = 60;
      
    }else if(mc2hessianCSV_str.Contains("NNPDF30_lo_as_0130_hessian_60.csv")){
      // 263000	NNPDF30_lo_as_0130	      101	
      pdfWeightLHAnumber_ = 263000;
      nPdfWeights_        = 100;
      nPdfEigWeights_     = 60;
    }
    cout 
    << "pdfWeightLHAnumber= " << pdfWeightLHAnumber_
    << " nPdfWeights= " << nPdfWeights_
    << " nPdfEigWeights= " << nPdfEigWeights_
    << endl;
  }
  
  pdfweightshelper_.Init(nPdfWeights_,nPdfEigWeights_,mc2hessianCSV);
  
}

void PDFWeightsProducer::produce( edm::Event & evt, const edm::EventSetup & ){

  auto_ptr<std::vector<float>> pdfeigweights_ptr( new std::vector<float>() );
  auto_ptr<unsigned int > pdfWeightLHAnumber_ptr( new unsigned int() );
  Handle<LHEEventProduct> lheInfo;
  evt.getByLabel( edm::InputTag("externalLHEProducer"), lheInfo );
  
  cout << "lheInfo->weights().size()= " << lheInfo->weights().size() << endl;
  cout << "PDFweightsLHEorder_.size()= " << PDFweightsLHEorder_.size() << endl;

  double nomlheweight = lheInfo->weights()[0].wgt;
  Handle<GenEventInfoProduct> genInfo;
  evt.getByToken( srcTokenGen_, genInfo );
  
  weight_ = genInfo->weight();
  
  //get the original mc replica weights
  std::vector<double> inpdfweights;
  for (unsigned int ipdf=0; ipdf<lheInfo->weights().size();ipdf++){
    
    if( PDFweightsLHEorder_[ipdf] > pdfWeightLHAnumber_
        && PDFweightsLHEorder_[ipdf] <= pdfWeightLHAnumber_+nPdfWeights_ ){

      //this is the weight to be used for evaluating uncertainties with mc replica weights
      // pdfweights_.push_back(lheInfo->weights()[ipdf].wgt*weight_/nomlheweight);
      
      //this is the raw weight to be fed to the mc2hessian convertor
      inpdfweights.push_back(lheInfo->weights()[ipdf].wgt);
      // foo = std::vector(myVec.begin () + 100000, myVec.begin () + 150000);
      
    }
    
    cout << "weight " << inpdfweights.size() 
    // << " --- PDFmap " << PDFweightsLHEorder_[ipdf] 
    << " --- PDFmap atoi " << PDFweightsLHEorder_[ipdf]
    << " --- lheInfo " << lheInfo->weights()[ipdf].wgt 
    << endl;
    
    if(inpdfweights.size() == nPdfWeights_) break;
  }
  
  std::vector<double> outpdfweights(nPdfEigWeights_);
  //do the actual conversion, where the nominal lhe weight is needed as the reference point for the linearization
  pdfweightshelper_.DoMC2Hessian(nomlheweight,inpdfweights.data(),outpdfweights.data());
  
  for (unsigned int iwgt=0; iwgt<nPdfEigWeights_; ++iwgt){
    double wgtval = outpdfweights[iwgt];
    
    //the is the weight to be used for evaluating uncertainties with hessian weights
    pdfeigweights_ptr->push_back(wgtval*weight_/nomlheweight);
    cout << "weight " << pdfeigweights_ptr->size() << " value " << (wgtval*weight_/nomlheweight) << endl;
  }

  *pdfWeightLHAnumber_ptr = pdfWeightLHAnumber_;
  
  evt.put( pdfeigweights_ptr, "outputHessianWeights" );
  evt.put( pdfWeightLHAnumber_ptr, "inputPDFset" );

}

void PDFWeightsProducer::endRun(edm::Run const& iRun, edm::EventSetup const& es){
}

void PDFWeightsProducer::endJob(){
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( PDFWeightsProducer );
