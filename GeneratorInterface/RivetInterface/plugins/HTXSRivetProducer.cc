/* \class HTXSRivetProducer
 *
 * \author David Sperka, University of Florida
 *
 * $Id: HTXSRivetProducer.cc,v 1.1 2016/09/27 13:07:29 dsperka Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "Rivet/AnalysisHandler.hh"
#include "GeneratorInterface/RivetInterface/src/HiggsTemplateCrossSections.cc"
#include "GeneratorInterface/RivetInterface/src/HiggsTemplateCrossSections.h"

#include <vector>
#include <string>

class HTXSRivetProducer : public edm::EDProducer {
public:

    explicit HTXSRivetProducer(const edm::ParameterSet&);
    ~HTXSRivetProducer();

private:

    virtual void beginJob();
    virtual void produce( edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    edm::EDGetTokenT<edm::HepMCProduct> _hepmcCollection;
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

};

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

using namespace Rivet;
using namespace edm;
using namespace std;

HTXSRivetProducer::HTXSRivetProducer( const ParameterSet & cfg ) {

    _hepmcCollection = consumes<HepMCProduct>(cfg.getParameter<edm::InputTag>("HepMCCollection"));

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

}

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

    if (_isFirstEvent){

        // set the production mode
        if      ( _prodMode == "GGF"   ) m_HiggsProdMode = HTXS::GGF;
        else if ( _prodMode == "VBF"   ) m_HiggsProdMode = HTXS::VBF;
        else if ( _prodMode == "WH"    ) m_HiggsProdMode = HTXS::WH;
        else if ( _prodMode == "ZH"    ) m_HiggsProdMode = HTXS::QQ2ZH;
        else if ( _prodMode == "QQ2ZH" ) m_HiggsProdMode = HTXS::QQ2ZH;
        else if ( _prodMode == "GG2ZH" ) m_HiggsProdMode = HTXS::GG2ZH;
        else if ( _prodMode == "TTH"   ) m_HiggsProdMode = HTXS::TTH;
        else if ( _prodMode == "BBH"   ) m_HiggsProdMode = HTXS::BBH;
        else if ( _prodMode == "TH"    ) m_HiggsProdMode = HTXS::TH;
        else {
            throw cms::Exception("HTXSRivetProducer") << "ProductionMode must be one of: GGF,VBF,WH,ZH,QQ2ZH,GG2ZH,TTH,BBH,TH " ; 
        }
        _HTXS->setHiggsProdMode(m_HiggsProdMode);

        // initialize rivet analysis
        _analysisHandler->init(*myGenEvent);
        _isFirstEvent = false;

    }

    // classify the event
    //cout<<"HTXSRivetProducer categorize... "<<endl;
    HiggsClassification cat = _HTXS->classifyEvent(*myGenEvent,m_HiggsProdMode);

    // print out some results
    //cout<<"HTXSRivetProducer cat.prodMode "<<cat.prodMode <<endl;
    //cout<<"HTXSRivetProducer cat.errorCode "<< cat.errorCode <<endl;
    //cout<<"HTXSRivetProducer cat.stage0_cat "<< cat.stage0_cat <<endl;
    //cout<<"HTXSRivetProducer cat.stage1_cat_pTjet30GeV "<< cat.stage1_cat_pTjet30GeV <<endl;
    //cout<<"HTXSRivetProducer njets30 " << cat.jets30.size() <<endl;
    //for (unsigned int i=0; i<cat.jets30.size(); i++) {
    //    cout<<"pTj"<<i<<" "<< cat.jets30[i].momentum().pt()<<" eta "<<cat.jets30[i].momentum().eta()<<endl;
    //    Rivet::Particles &jetparticles = cat.jets30[i].particles();
    //    for ( auto p:jetparticles ) {
    //        cout<<p.pdgId()<<" "<<p.pt()<<" "<<p.eta()<<endl; 
    //    }
    //}
    //if (cat.jets30.size()>=2) {
    //    const Rivet::FourMomentum &j1=cat.jets30[0].momentum(), &j2=cat.jets30[1].momentum();
    //    cout<<"mass "<<(j1+j2).mass()<<" delta_rap: "<<std::abs(j1.rapidity()-j2.rapidity())<<endl;
    //    cout<<"pthjj: "<<(j1+j2+cat.higgs.momentum()).pt()<<endl;
    //}

    stage0cat_ = cat.stage0_cat;
    stage1cat_ = cat.stage1_cat_pTjet30GeV;
    njets_ = cat.jets30.size();
    pTH_ = cat.higgs.momentum().pt();
    pTV_ = cat.V.momentum().pt();

    unique_ptr<int> stage0cat( new int( stage0cat_ ) );
    unique_ptr<int> stage1cat( new int( stage1cat_ ) );
    unique_ptr<int> njets( new int( njets_ ) );
    unique_ptr<float> pTH( new float( pTH_ ) );
    unique_ptr<float> pTV( new float( pTV_ ) );

    iEvent.put(std::move(stage0cat),"stage0cat");
    iEvent.put(std::move(stage1cat),"stage1cat");
    iEvent.put(std::move(njets),"njets");
    iEvent.put(std::move(pTH),"pTH");
    iEvent.put(std::move(pTV),"pTV");

}

void HTXSRivetProducer::endJob(){
    _HTXS->printClassificationSummary();
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( HTXSRivetProducer );
