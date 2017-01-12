/* \class HTXSRivetProducer
 *
 * \author David Sperka, University of Florida
 *
 * $Id: HTXSRivetProducer.cc,v 1.1 2016/09/27 13:07:29 dsperka Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Run.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

#include "Rivet/AnalysisHandler.hh"
#include "GeneratorInterface/RivetInterface/src/HiggsTemplateCrossSections.cc"
#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"

#include <vector>
#include <stdio.h>
#include <string.h>


using namespace Rivet;
using namespace edm;
using namespace std;

class HTXSRivetProducer : public edm::EDProducer {
public:
    
    explicit HTXSRivetProducer(const edm::ParameterSet& cfg) : 
        _hepmcCollection(consumes<HepMCProduct>(cfg.getParameter<edm::InputTag>("HepMCCollection"))),
        _lheRunInfo(consumes<LHERunInfoProduct,edm::InRun>(cfg.getParameter<edm::InputTag>("LHERunInfo")))
    {
        _analysisHandler = new Rivet::AnalysisHandler();
        _HTXS = new Rivet::HiggsTemplateCrossSections();
        
        _isFirstEvent = true;
        _prodMode = cfg.getParameter<string>("ProductionMode");
        m_HiggsProdMode = HTXS::UNKNOWN;
        
        produces<HTXS::HiggsClassification>("HiggsClassification").setBranchAlias("HiggsClassification");

    }
    ~HTXSRivetProducer();
    
private:
    
    virtual void beginJob();
    virtual void produce( edm::Event&, const edm::EventSetup&);
    virtual void endJob();
    
    virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& es);
    virtual void endRun(edm::Run const& iRun, edm::EventSetup const& es);
    
    edm::EDGetTokenT<edm::HepMCProduct> _hepmcCollection;
    edm::EDGetTokenT<LHERunInfoProduct> _lheRunInfo;
    
    Rivet::AnalysisHandler* _analysisHandler;
    Rivet::HiggsTemplateCrossSections* _HTXS;
    
    bool _isFirstEvent;
    std::string _prodMode;
    HTXS::HiggsProdMode m_HiggsProdMode;
    
    HTXS::HiggsClassification cat_;
    
};

HTXSRivetProducer::~HTXSRivetProducer(){
}

void HTXSRivetProducer::beginJob(){
    _analysisHandler->addAnalysis(_HTXS);
}

void HTXSRivetProducer::produce( edm::Event & iEvent, const edm::EventSetup & ) {
    
    //get the hepmc product from the event
    edm::Handle<HepMCProduct> evt;
    iEvent.getByToken(_hepmcCollection, evt);
    
    // get HepMC GenEvent
    const HepMC::GenEvent *myGenEvent = evt->GetEvent();
    
    if (_prodMode == "AUTO") {
        
        // for these prod modes, don't change what is set in BeginRun
        if (m_HiggsProdMode != HTXS::GGF && m_HiggsProdMode != HTXS::VBF && m_HiggsProdMode != HTXS::GG2ZH) {
            
            unsigned nWs = 0;
            unsigned nZs = 0;
            unsigned nTs = 0;
            unsigned nBs = 0;
            unsigned nHs = 0;
            
            HepMC::GenVertex *HSvtx = myGenEvent->signal_process_vertex();
            
            for (auto ptcl:particles(HSvtx,HepMC::children)) {
                if (ptcl->pdg_id() == 24 || ptcl->pdg_id() == -24) ++nWs;
                if (ptcl->pdg_id() == 23) ++nZs;
                if (abs(ptcl->pdg_id()) == 6) ++nTs;
                if (abs(ptcl->pdg_id()) == 5) ++nBs;
                if (ptcl->pdg_id() == 25) ++nHs;
            }
                        
            if (nZs==1 && nHs==1 && (nWs+nTs)==0) {
                m_HiggsProdMode = HTXS::QQ2ZH;
            } else if (nWs==1 && nHs==1 && (nZs+nTs)==0) {
                m_HiggsProdMode = HTXS::WH; 
            } else if (nTs==2 && nHs==1 && nZs==0) {
                m_HiggsProdMode = HTXS::TTH; 
            } else if (nTs==1 && nHs==1 && nZs==0) {
                m_HiggsProdMode = HTXS::TH; 
            } else if (nBs==2 && nHs==1 && nZs==0) {
                m_HiggsProdMode = HTXS::BBH; 
            }
            
            _HTXS->setHiggsProdMode(m_HiggsProdMode);
            
        }
    }
    

    if (_isFirstEvent){

        // set the production mode if not done already
        if      ( _prodMode == "GGF"   ) m_HiggsProdMode = HTXS::GGF;
        else if ( _prodMode == "VBF"   ) m_HiggsProdMode = HTXS::VBF;
        else if ( _prodMode == "WH"    ) m_HiggsProdMode = HTXS::WH;
        else if ( _prodMode == "ZH"    ) m_HiggsProdMode = HTXS::QQ2ZH;
        else if ( _prodMode == "QQ2ZH" ) m_HiggsProdMode = HTXS::QQ2ZH;
        else if ( _prodMode == "GG2ZH" ) m_HiggsProdMode = HTXS::GG2ZH;
        else if ( _prodMode == "TTH"   ) m_HiggsProdMode = HTXS::TTH;
        else if ( _prodMode == "BBH"   ) m_HiggsProdMode = HTXS::BBH;
        else if ( _prodMode == "TH"    ) m_HiggsProdMode = HTXS::TH;
        else if ( _prodMode == "AUTO"  ) {
            cout<<"Using AUTO for HiggsProdMode, found it to be: "<<m_HiggsProdMode<< "\n";
            cout<<"(UNKNOWN=0, GGF=1, VBF=2, WH=3, QQ2ZH=4, GG2ZH=5, TTH=6, BBH=7, TH=8)"<<endl;
        } else {
            throw cms::Exception("HTXSRivetProducer") << "ProductionMode must be one of: GGF,VBF,WH,ZH,QQ2ZH,GG2ZH,TTH,BBH,TH,AUTO " ; 
        }
        _HTXS->setHiggsProdMode(m_HiggsProdMode);

        // at this point the production mode must be known
        if (m_HiggsProdMode == HTXS::UNKNOWN) {
            throw cms::Exception("HTXSRivetProducer") << "HiggsProduction mode is UNKNOWN";
        }            

        // initialize rivet analysis
        _analysisHandler->init(*myGenEvent);
        _isFirstEvent = false;

    }

    // classify the event
    Rivet::HiggsClassification rivet_cat = _HTXS->classifyEvent(*myGenEvent,m_HiggsProdMode);
    cat_ = HTXS::Rivet2Root(rivet_cat);

    unique_ptr<HTXS::HiggsClassification> cat( new HTXS::HiggsClassification( cat_ ) ); 

    iEvent.put(std::move(cat),"HiggsClassification");

}

void HTXSRivetProducer::endJob(){
    _HTXS->printClassificationSummary();
}

void HTXSRivetProducer::endRun(edm::Run const& iRun, edm::EventSetup const& es) 
{
}

void HTXSRivetProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& es) {

    if (_prodMode == "AUTO") {    

        edm::Handle<LHERunInfoProduct> run;
        iRun.getByLabel( edm::InputTag("externalLHEProducer"), run );
        
        typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
        LHERunInfoProduct myLHERunInfoProduct = *(run.product());
        for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
            std::vector<std::string> lines = iter->lines();
            for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
                std::string line=lines.at(iLine);
                // POWHEG
                if (strstr(line.c_str(),"gg_H_quark-mass-effects")) {
                    std::cout<<iLine<< " "<<line<<std::endl;
                    m_HiggsProdMode = HTXS::GGF;
                    break;
                }
                if (strstr(line.c_str(),"VBF_H")) {
                    std::cout<<iLine<< " "<<line<<std::endl;
                    m_HiggsProdMode = HTXS::VBF;
                    break;
                }
                if (strstr(line.c_str(),"HZJ")) {
                    std::cout<<iLine<< " "<<line<<std::endl;
                    m_HiggsProdMode = HTXS::QQ2ZH;
                    break;
                }
                if (strstr(line.c_str(),"ggHZ")) {
                    std::cout<<iLine<< " "<<line<<std::endl;
                    m_HiggsProdMode = HTXS::GG2ZH;
                    break;
                }
                // MC@NLO
                if (strstr(line.c_str(),"ggh012j")) {
                    std::cout<<iLine<< " "<<line<<std::endl;
                    m_HiggsProdMode = HTXS::GGF;
                    break;
                }
                if (strstr(line.c_str(),"vbfh")) {
                    std::cout<<iLine<< " "<<line<<std::endl;
                    m_HiggsProdMode = HTXS::VBF;
                    break;
                }
                if (strstr(line.c_str(),"zh012j")) {
                    std::cout<<iLine<< " "<<line<<std::endl;
                    m_HiggsProdMode = HTXS::QQ2ZH;
                    break;
                }
                if (strstr(line.c_str(),"ggzh01j")) {
                    std::cout<<iLine<< " "<<line<<std::endl;
                    m_HiggsProdMode = HTXS::GG2ZH;
                    break;
                }
            }
              
            if ( m_HiggsProdMode != HTXS::UNKNOWN) break;
        }
    }
    
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( HTXSRivetProducer );
