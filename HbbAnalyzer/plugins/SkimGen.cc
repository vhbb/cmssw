// -*- C++ -*-
//
// Package:    SkimGen
// Class:      SkimGen
// 
/**\class SkimGen SkimGen.cc myGEN/SkimGen/src/SkimGen.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Leonard Apanasevich
//         Created:  Thu Sep 15 13:04:24 CDT 2011
// $Id: SkimGen.cc,v 1.1 2012/12/01 20:41:12 apana Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

//#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
//#include "TH1F.h"

//#include "HepMC/WeightContainer.h"
//#include "HepMC/GenEvent.h"
//#include "HepMC/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

typedef std::vector<int> vint;

//
// class declaration
//

class SkimGen : public edm::EDFilter {
   public:
      explicit SkimGen(const edm::ParameterSet&);
      ~SkimGen();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

  edm::InputTag mcTruthCollection_;
  vint genparts_;
  double genpt_;
  bool Debug_;
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
SkimGen::SkimGen(const edm::ParameterSet& iConfig)
  : mcTruthCollection_(iConfig.getParameter<edm::InputTag>("mcTruthCollection"))
  , genparts_(iConfig.getParameter<vint>("genparts"))
  , genpt_(iConfig.getParameter<double>("genpt"))
  , Debug_(iConfig.getParameter<bool>("Debug"))
{
   //now do what ever initialization is needed

}


SkimGen::~SkimGen()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SkimGen::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   bool accept(false);

   Handle<reco::GenParticleCollection> genParticles ;
   iEvent.getByLabel(mcTruthCollection_, genParticles) ;

   // if (Debug_) std::cout << "GenPart: " << genpart_ << " pt cut:" << genpt_ << std::endl;
   reco::GenParticleCollection::const_iterator mcIter ;
   int iline=0;
   for ( mcIter=genParticles->begin() ; mcIter!=genParticles->end() ; mcIter++ ){
     iline++;
     if(// mcIter->status() == 1 &&  
	find(genparts_.begin(), genparts_.end(), abs((mcIter->pdgId()))) != genparts_.end()){
       // found generated particle
       double pt=mcIter->pt();
       if (Debug_) std::cout <<"Event: "    << iEvent.id().event()
		 <<"  Found generated particle " << mcIter->pdgId()
		 <<"  with pt: "<< pt 
		 <<"  mass: "   << mcIter->mass() 
		 <<"  Status: " << mcIter->status()
		 <<"  Line: " << iline
		 << std::endl; 
       if (pt>genpt_ ){
	   accept=true;
	   break;
       }
     }
   } // end loop over "gen"
   if (accept && Debug_) std::cout << "Accepted event" << std::endl;
   return accept;
}

// ------------ method called once each job just before starting event loop  ------------
void 
SkimGen::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SkimGen::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
SkimGen::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
SkimGen::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
SkimGen::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
SkimGen::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SkimGen::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(SkimGen);
