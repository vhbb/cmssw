// -*- C++ -*-
//
// Package:    SubProducer
// Class:      SubProducer
// 
/**\class SubProducer SubProducer.cc Prod/SubProducer/src/SubProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  David Lopes Pegna,Address unknown,NONE,
//         Created:  Wed Jul 22 19:17:38 EDT 2009
// $Id: SubProducer.cc,v 1.1 2011/02/22 21:52:54 dlopes Exp $
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

#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"

//
// class decleration
//

class SubProducer : public edm::EDProducer {
   public:
      explicit SubProducer(const edm::ParameterSet&);
      ~SubProducer();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------

  edm::InputTag jetLabel_;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
SubProducer::SubProducer(const edm::ParameterSet& iConfig):
jetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("jetTag"))
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
*/
   //now do what ever other initialization is needed
produces<reco::BasicJetCollection>();
  
}


SubProducer::~SubProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
SubProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;

   std::auto_ptr<BasicJetCollection> SubJetList( new BasicJetCollection() );

   edm::Handle<BasicJetCollection> jetHandle;
   iEvent.getByLabel(jetLabel_,jetHandle);

     for(BasicJetCollection::const_iterator it(jetHandle->begin()), itEnd(jetHandle->end());  it!=itEnd;++it)
       {
      reco::Jet::Constituents constituents = it->getJetConstituents();
 
        for (unsigned int iJC(0); iJC<constituents.size(); ++iJC ){
        Jet::Constituent icandJet = constituents[iJC];

  //      std::cout << "and here : " << icandJet->pt() << "\n";
        const BasicJet * subjet = static_cast<const BasicJet * >(&*icandJet);

  //      if(subjet!=0) std::cout << "here I am: " << subjet->pt() << "\n"; 
        if(subjet!=0) SubJetList->push_back(*subjet);

        }  

      }

     iEvent.put(SubJetList);

}

// ------------ method called once each job just before starting event loop  ------------
void 
SubProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SubProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SubProducer);
