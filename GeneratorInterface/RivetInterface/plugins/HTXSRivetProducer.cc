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
//#include "GeneratorInterface/RivetInterface/src/HiggsTemplateCrossSections.h"
#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"

#include <vector>
#include <stdio.h>
#include <string.h>


using namespace Rivet;
using namespace edm;
using namespace std;
using namespace lhef;


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
        
        produces<int>("stage0cat").setBranchAlias("stage0cat");
        produces<int>("stage1cat").setBranchAlias("stage1cat");
        produces<int>("njets").setBranchAlias("njets");
        produces<float>("pTH").setBranchAlias("pTH");
        produces<float>("pTV").setBranchAlias("pTV");
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
    
    int stage0cat_;
    int stage1cat_;
    int njets_;
    float pTH_;
    float pTV_;
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
    
    if ( _prodMode == "AUTO" || _prodMode == "VH") {
        
        unsigned nWs = 0;
        unsigned nZs = 0;
        unsigned nTs = 0;
        unsigned nHs = 0;
        
        HepMC::GenVertex *HSvtx = myGenEvent->signal_process_vertex();

        for (auto ptcl:particles(HSvtx,HepMC::children)) {
            if (ptcl->pdg_id() == 24 || ptcl->pdg_id() == -24) ++nWs;
            if (ptcl->pdg_id() == 23) ++nZs;
            if (abs(ptcl->pdg_id()) == 6) ++nTs;
            if (ptcl->pdg_id() == 25) ++nHs;
        }

        if (_prodMode == "VH" && nZs == 0 && nWs == 0) {
            throw cms::Exception("HTXSRivetProducer") << "ProductionMode set to ZH but no W or Z is found ";
        }
        if (_prodMode == "VH" && nZs > 0 && nWs > 0) {
            throw cms::Exception("HTXSRivetProducer") << "ProductionMode set to ZH and both Ws and Zs are found ";
        }

        if (nZs==1 && nHs==1 && (nWs+nTs)==0) {
            m_HiggsProdMode = HTXS::QQ2ZH;
        } else if (nWs==1 && nHs==1 && (nZs+nTs)==0) {
            m_HiggsProdMode = HTXS::WH; 
        } else if (nTs==2 && nHs==1) {
            m_HiggsProdMode = HTXS::TTH; 
        }

        _HTXS->setHiggsProdMode(m_HiggsProdMode);

    }

    if (_isFirstEvent){

        // set the production mode
        if      ( _prodMode == "GGF"   ) m_HiggsProdMode = HTXS::GGF;
        else if ( _prodMode == "VBF"   ) m_HiggsProdMode = HTXS::VBF;
        else if ( _prodMode == "WH"    ) m_HiggsProdMode = HTXS::WH;
        else if ( _prodMode == "ZH"    ) m_HiggsProdMode = HTXS::QQ2ZH;
        else if ( _prodMode == "VH"    ) m_HiggsProdMode = HTXS::WH;  // But we will check every event and switch on the fly
        else if ( _prodMode == "QQ2ZH" ) m_HiggsProdMode = HTXS::QQ2ZH;
        else if ( _prodMode == "GG2ZH" ) m_HiggsProdMode = HTXS::GG2ZH;
        else if ( _prodMode == "TTH"   ) m_HiggsProdMode = HTXS::TTH;
        else if ( _prodMode == "BBH"   ) m_HiggsProdMode = HTXS::BBH;
        else if ( _prodMode == "TH"    ) m_HiggsProdMode = HTXS::TH;
        else if ( _prodMode == "AUTO"  ) cout<<"Using AUTO for HiggsProdMode, found it to be: "<<m_HiggsProdMode<<endl;
        else {
            throw cms::Exception("HTXSRivetProducer") << "ProductionMode must be one of: GGF,VBF,WH,ZH,VH,QQ2ZH,GG2ZH,TTH,BBH,TH,AUTO " ; 
        }
        _HTXS->setHiggsProdMode(m_HiggsProdMode);

        // initialize rivet analysis
        _analysisHandler->init(*myGenEvent);
        _isFirstEvent = false;

    }

    // classify the event
    Rivet::HiggsClassification rivet_cat = _HTXS->classifyEvent(*myGenEvent,m_HiggsProdMode);
    cat_ = HTXS::Rivet2Root(rivet_cat);

    // print out some results
    //cout<<"HTXSRivetProducer rivet_cat.prodMode "<<rivet_cat.prodMode <<endl;
    //cout<<"HTXSRivetProducer rivet_cat.errorCode "<< rivet_cat.errorCode <<endl;
    //cout<<"HTXSRivetProducer rivet_cat.stage0_cat "<< rivet_cat.stage0_cat <<endl;
    //cout<<"HTXSRivetProducer rivet_cat.stage1_rivet_catpTjet30GeV "<< rivet_cat.stage1_cat_pTjet30GeV <<endl;
    //cout<<"HTXSRivetProducer njets30 " << rivet_cat.jets30.size() <<endl;
    //for (unsigned int i=0; i<rivet_cat.jets30.size(); i++) {
    //    cout<<"pTj"<<i<<" "<< rivet_cat.jets30[i].momentum().pt()<<" eta "<<rivet_cat.jets30[i].momentum().eta()<<endl;
    //    Rivet::Particles &jetparticles = rivet_cat.jets30[i].particles();
    //    for ( auto p:jetparticles ) {
    //        cout<<p.pdgId()<<" "<<p.pt()<<" "<<p.eta()<<endl; 
    //    }
    //}
    //if (rivet_cat.jets30.size()>=2) {
    //    const Rivet::FourMomentum &j1=rivet_cat.jets30[0].momentum(), &j2=rivet_cat.jets30[1].momentum();
    //    cout<<"mass "<<(j1+j2).mass()<<" delta_rap: "<<std::abs(j1.rapidity()-j2.rapidity())<<endl;
    //    cout<<"pthjj: "<<(j1+j2+rivet_cat.higgs.momentum()).pt()<<endl;
    //}

    stage0cat_ = cat_.stage0_cat;
    stage1cat_ = cat_.stage1_cat_pTjet30GeV;
    njets_ = cat_.jets30.size();
    pTH_ = cat_.higgs.Pt();
    pTV_ = cat_.V.Pt();

    unique_ptr<int> stage0cat( new int( stage0cat_ ) );
    unique_ptr<int> stage1cat( new int( stage1cat_ ) );
    unique_ptr<int> njets( new int( njets_ ) );
    unique_ptr<float> pTH( new float( pTH_ ) );
    unique_ptr<float> pTV( new float( pTV_ ) );
    unique_ptr<HTXS::HiggsClassification> cat( new HTXS::HiggsClassification( cat_ ) ); 

    iEvent.put(std::move(stage0cat),"stage0cat");
    iEvent.put(std::move(stage1cat),"stage1cat");
    iEvent.put(std::move(njets),"njets");
    iEvent.put(std::move(pTH),"pTH");
    iEvent.put(std::move(pTV),"pTV");
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
        //iRun.getByToken( _lheRunInfo, run );
        iRun.getByLabel( edm::InputTag("externalLHEProducer"), run );
        
        typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
        LHERunInfoProduct myLHERunInfoProduct = *(run.product());
        for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
            std::cout << iter->tag() << std::endl;
            std::vector<std::string> lines = iter->lines();
            for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
                std::string line=lines.at(iLine);
                //std::cout<<iLine<< " "<<line<<std::endl;
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

            }
              
            if ( m_HiggsProdMode != HTXS::UNKNOWN) break;
        }
    }
    
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( HTXSRivetProducer );
