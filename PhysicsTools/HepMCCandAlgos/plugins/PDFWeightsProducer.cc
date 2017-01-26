/* \class PDFWeightsProducer
*
* \author Luca Perrozzi, ETHZ
*
* Convert PDF replicas in Hessian representation
* Originally developed by Josh BEndavid, adapted by Luca Perrozzi
*
*
*/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
// #include "FWCore/Framework/interface/EDProducer.h"
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

class PDFWeightsProducer : public edm::one::EDProducer<edm::BeginRunProducer,
edm::EndRunProducer> {
public:
  
  explicit PDFWeightsProducer( const ParameterSet & cfg ) :
  lheRunInfo_(consumes<LHERunInfoProduct,edm::InRun>(edm::InputTag("externalLHEProducer"))),
  srcTokenlheEvent_( consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"))),
  srcTokenGen_( consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  pdfWeightLHAnumber_(cfg.getParameter<unsigned int>("pdfWeightLHAnumber")),
  nPdfWeights_(cfg.getParameter<unsigned int>("nPdfWeights")),
  nPdfEigWeights_(cfg.getParameter<unsigned int>("nPdfEigWeights"))
  {
    produces<std::vector<float>>("outputHessianWeights");
    produces<unsigned int, edm::InRun>("pdfWeightLHAnumber");
    mc2hessianCSV = cfg.getParameter<edm::FileInPath>("mc2hessianCSV");
  }

  ~PDFWeightsProducer();
  
private:
  
  virtual void beginRunProduce(edm::Run & run, edm::EventSetup const& es) override;
  virtual void endRunProduce(edm::Run & run, edm::EventSetup const& es) override;
  virtual void produce( edm::Event&, const edm::EventSetup&) override;
  virtual void beginJob();
  virtual void endJob();
  
  edm::FileInPath mc2hessianCSV;
  PDFWeightsHelper pdfweightshelper_;
  EDGetTokenT<LHERunInfoProduct> lheRunInfo_;
  EDGetTokenT<LHEEventProduct> srcTokenlheEvent_;
  EDGetTokenT<GenEventInfoProduct> srcTokenGen_;
  
  unsigned int pdfWeightLHAnumber_;
  unsigned int nPdfWeights_;
  unsigned int nPdfEigWeights_;
  unsigned int PDFweightsLHEoffset_;
  
  Handle<LHERunInfoProduct> run;

  float weight_;
  std::vector<unsigned int> PDFweightsLHEorder_;
  std::vector<unsigned int> PDFweightsLHEoffsetBuffer_;
  std::vector<unsigned int> PDFweightsLHEsets_;
  
};

PDFWeightsProducer::~PDFWeightsProducer(){
}

void PDFWeightsProducer::beginJob(){
}

void PDFWeightsProducer::beginRunProduce(edm::Run & iRun, edm::EventSetup const& es) {

  // Access LHERunInfoProduct (i.e. LHE header) to retrieve the list of weights in the sample
  bool product_exists = iRun.getByLabel( edm::InputTag("externalLHEProducer"), run );
  if( product_exists ){
    TString mc2hessianCSV_str = mc2hessianCSV.fullPath();

    typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
    LHERunInfoProduct myLHERunInfoProduct = *(run.product());
    for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
      std::vector<std::string> lines = iter->lines();
      for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
        TString line_tstr =lines.at(iLine).c_str();
        // cout<<"before weight id= " << iLine<< " "<<line_tstr.ReplaceAll("\n","").Data()<<endl;
        if (line_tstr.Contains("<weight id=")) {
          // cout<<iLine<< " "<<line_tstr.ReplaceAll("\n","").Data()<<endl;
          line_tstr = line_tstr.ReplaceAll("\n","").ReplaceAll("<weight id=\"","").ReplaceAll("</weight>","");
          TObjArray *tx = line_tstr.Tokenize("> ");
          line_tstr = ((TObjString *)(tx->At(tx->GetEntries()-1)))->String();
          unsigned int pdfnumber = line_tstr.ReplaceAll("PDF set = ","").ReplaceAll("pdfset=","").Atoi();
          PDFweightsLHEorder_.push_back(pdfnumber);
          if(pdfnumber%100==1){
            unsigned int LHEset = pdfnumber-1;
            unsigned int LHEposition = PDFweightsLHEorder_.size()-1;
            cout << LHEposition << " WEIGHT - PDF SET CANDIDATE FOR REWEIGHTING IN THIS SAMPLE= " << (LHEset) << endl;
            PDFweightsLHEsets_.push_back(LHEset);
            PDFweightsLHEoffsetBuffer_.push_back(LHEposition);
            
            // Unless details are specified in the cfg, check whether the sample actually contains 
            // the weights of the appropriate PDF set to be converted to Hessian.
            // A check is performed on the NLO NNPDF samples to catch the _pdfas PDF set in case is used
            if(pdfWeightLHAnumber_<1 || nPdfWeights_<1 || nPdfEigWeights_<1){
              
              if(
                  // 260000	NNPDF30_nlo_as_0118	      101
                  // 292200	NNPDF30_nlo_nf_5_pdfas	  103	
                  (mc2hessianCSV_str.Contains("NNPDF30_nlo_as_0118_hessian_60.csv") && (LHEset == 260000 || LHEset == 292200))
                  // 260400	NNPDF30_nlo_as_0118_nf_4	101	
                  // 292000	NNPDF30_nlo_nf_4_pdfas	  103	
                  || (mc2hessianCSV_str.Contains("NNPDF30_nlo_as_0118_nf_4_hessian_60.csv") && (LHEset == 260400 || LHEset == 292000))
                  // 263400	NNPDF30_lo_as_0130_nf_4	  101	
                  || (mc2hessianCSV_str.Contains("NNPDF30_lo_as_0130_nf_4_hessian_60.csv") && (LHEset == 263400))
                  // 263000	NNPDF30_lo_as_0130	      101	
                  || (mc2hessianCSV_str.Contains("NNPDF30_lo_as_0130_hessian_60.csv") && (LHEset == 263000))
                  ){
                cout << "Either pdfWeightLHAnumber_<1 or nPdfWeights_<1 or nPdfEigWeights_<1: using hardcoded settings for" << endl;
                cout << "mc2hessianCSV " << mc2hessianCSV_str << endl;
                pdfWeightLHAnumber_  = LHEset;
                PDFweightsLHEoffset_ = LHEposition;
                nPdfWeights_         = 100;
                nPdfEigWeights_      = 60;
                cout 
                << "LHEset= " << LHEset
                << " nPdfWeights= " << nPdfWeights_
                << " nPdfEigWeights= " << nPdfEigWeights_
                << endl;
                break;
              }
            }else if (LHEset == pdfWeightLHAnumber_){
              PDFweightsLHEoffset_ = LHEposition;
              break;
            }
            
          }
        }
      }
    }
    if(pdfWeightLHAnumber_<1 || nPdfWeights_<1 || nPdfEigWeights_<1){
      cout << "NO SUITABLE SET FOUND FOR MC2HESSIAN PDF CONVERSION!" << endl;
    }

    pdfweightshelper_.Init(nPdfWeights_,nPdfEigWeights_,mc2hessianCSV);

    auto_ptr<unsigned int > p(new unsigned int() );
    *p = pdfWeightLHAnumber_;
    iRun.put(p, "pdfWeightLHAnumber");
  }
}

void PDFWeightsProducer::produce( edm::Event & evt, const edm::EventSetup & ){

  auto_ptr<std::vector<float>> pdfeigweights_ptr( new std::vector<float>() );
  // auto_ptr<unsigned int > pdfWeightLHAnumber_ptr( new unsigned int() );
  Handle<LHEEventProduct> lheInfo;
  bool product_exists = evt.getByLabel(edm::InputTag("externalLHEProducer"), lheInfo );
  if( product_exists && !(pdfWeightLHAnumber_<1 || nPdfWeights_<1 || nPdfEigWeights_<1)){
    
    // cout << "lheInfo->weights().size()= " << lheInfo->weights().size() << endl;
    // cout << "PDFweightsLHEorder_.size()= " << PDFweightsLHEorder_.size() << endl;

    double nomlheweight = lheInfo->weights()[0].wgt;
    Handle<GenEventInfoProduct> genInfo;
    evt.getByToken( srcTokenGen_, genInfo );
    
    weight_ = genInfo->weight();
    
    //get the original mc replica weights
    std::vector<double> inpdfweights;
    // for (unsigned int ipdf=0; ipdf<lheInfo->weights().size();ipdf++){
    for (unsigned int ipdf=PDFweightsLHEoffset_; ipdf<PDFweightsLHEoffset_+nPdfWeights_;ipdf++){
      
      if( PDFweightsLHEorder_[ipdf] > pdfWeightLHAnumber_
          && PDFweightsLHEorder_[ipdf] <= pdfWeightLHAnumber_+nPdfWeights_ ){

        //this is the raw weight to be fed to the mc2hessian convertor
        inpdfweights.push_back(lheInfo->weights()[ipdf].wgt);
        
      }
      
      // cout << "weight " << ipdf << " --> " << inpdfweights.size() 
      // << " --- PDFmap atoi " << PDFweightsLHEorder_[ipdf]
      // << " --- lheInfo " << lheInfo->weights()[ipdf].wgt 
      // << endl;
      
      if(inpdfweights.size() == nPdfWeights_) break;
    }
    
    std::vector<double> outpdfweights(nPdfEigWeights_);
    //do the actual conversion, where the nominal lhe weight is needed as the reference point for the linearization
    pdfweightshelper_.DoMC2Hessian(nomlheweight,inpdfweights.data(),outpdfweights.data());
    
    for (unsigned int iwgt=0; iwgt<nPdfEigWeights_; ++iwgt){
      double wgtval = outpdfweights[iwgt];
      
      //the is the weight to be used for evaluating uncertainties with hessian weights
      pdfeigweights_ptr->push_back(wgtval*weight_/nomlheweight);
      // cout << "weight " << pdfeigweights_ptr->size() << " value " << (wgtval*weight_/nomlheweight) << endl;
    }

    // *pdfWeightLHAnumber_ptr = pdfWeightLHAnumber_;
    
    evt.put( pdfeigweights_ptr, "outputHessianWeights" );
    // evt.put( pdfWeightLHAnumber_ptr, "inputPDFset" );

  }
}

void PDFWeightsProducer::endRunProduce(edm::Run& run, edm::EventSetup const& es){
}

void PDFWeightsProducer::endJob(){
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( PDFWeightsProducer );
