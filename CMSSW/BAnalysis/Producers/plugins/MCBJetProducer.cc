// -*- C++ -*-
//
// Package:    MCBJetProducer
// Class:      MCBJetProducer
// 
/**\class MCBJetProducer MCBJetProducer.cc BAnalsysis/Producers/src/MCBJetProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Lukas WEHRLI
//         Created:  Wed Apr  1 15:21:27 CET 2009
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/Vector3D.h"
// Math
#include "Math/GenVector/VectorUtil.h"

#include "AnalysisDataFormats/BAnalysis/interface/SimBHadron.h"
#include "AnalysisDataFormats/BAnalysis/interface/SimBJet.h"

 #include "DataFormats/JetReco/interface/JetCollection.h"


//
// class decleration
//


class MCBJetProducer : public edm::EDProducer {
   public:
      explicit MCBJetProducer(const edm::ParameterSet&);
  ~MCBJetProducer();
  
   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob();


      // ----------member data ---------------------------

      //edm::InputTag jetsrc_;
      //int noEvent;


};

//constants, enums and typedefs

//typedef std::set<int>::const_iterator intsetit;
//typedef SimBJetCollection::iterator SBHit;
//typedef std::set<const reco::Candidate *>::iterator candsetit; 
//typedef std::vector<reco::Jet> JetCollection;

//static data member definitions



//constructors and destructor

MCBJetProducer::MCBJetProducer(const edm::ParameterSet& iConfig)// :
  //jetsrc_(iConfig.getParameter<edm::InputTag>("jetsrc"))

{
 produces<JetCollection>(); 


   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
*/
   //now do what ever other initialization is needed

}


MCBJetProducer::~MCBJetProducer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MCBJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;



   //////////////////////////////////////////
   //--------------SIMBHADRONS----------------------"
   //SimBHadron
   Handle<SimBHadronCollection> SBHC;
   iEvent.getByLabel("BHadrons", SBHC); 
   const SimBHadronCollection sbhcoll = *(SBHC.product());
   //////////////////////////////////////////

   //output collection
   JetCollection * output = new JetCollection();

/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
*/ 

   for(unsigned int b=0; b<sbhcoll.size(); b++){
     
     //math::XYZPoint p(sbhcoll[b].vx(), sbhcoll[b].vy(), sbhcoll[b].vz());
     //temp: 
     Jet jet(sbhcoll[b].p4(), sbhcoll[b].vertex()); 
     
     output->push_back(jet);
   }



   std::auto_ptr<JetCollection> pOut(output);
   iEvent.put(pOut);

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/

}

// ------------ method called once each job just before starting event loop  ------------
void 
MCBJetProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MCBJetProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MCBJetProducer);
