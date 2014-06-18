// -*- C++ -*-
//
// Package:    MCBHadronProducer
// Class:      MCBHadronProducer
// 
/**\class MCBHadronProducer MCBHadronProducer.cc BAnalsysis/MCBHadronProducer/src/MCBHadronProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Andrea RIZZI
//         Created:  Wed Dec  3 12:14:27 CET 2008
// $Id: MCBHadronProducer.cc,v 1.2 2011/10/31 11:17:12 arizzi Exp $
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

//#include "AnalysisDataFormats/BAnalysis/interface/SimBHadron.h"
//#include "AnalysisDataFormats/BAnalysis/interface/SimSecondaryVertex.h"

#include "ZSV/BAnalysis/interface/SimBHadron.h"
#include "ZSV/BAnalysis/interface/SimSecondaryVertex.h"

//
// class decleration
//

class MCBHadronProducer : public edm::EDProducer {
   public:
      explicit MCBHadronProducer(const edm::ParameterSet&);
  ~MCBHadronProducer();
  
   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob();
  
      void removeMothers(const reco::Candidate &gp, SimBHadronCollection &out);

      void setMoDauQuantities(reco::GenParticle gp, SimBHadron &bh);
      void findbquarks(const reco::Candidate &gp, std::set<const reco::Candidate *> &set); 


      // ----------member data ---------------------------

      int noEvent;


};

//constants, enums and typedefs

typedef std::set<int>::const_iterator intsetit;
typedef SimBHadronCollection::iterator SBHit;
typedef std::set<const reco::Candidate *>::iterator candsetit; 

//static data member definitions



//constructors and destructor

MCBHadronProducer::MCBHadronProducer(const edm::ParameterSet& iConfig) 

{
 produces<SimBHadronCollection>(); 


   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
*/
   //now do what ever other initialization is needed

}


MCBHadronProducer::~MCBHadronProducer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MCBHadronProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  noEvent++;
  ////std::cout << "Event NO " << noEvent << std::endl;
   using namespace edm;
   using namespace reco;

   //output collection
   SimBHadronCollection * output = new SimBHadronCollection();


  /////////////////////////////////////////////
   //--------------GenParticles
   Handle<GenParticleCollection> GPCollection;
   iEvent.getByLabel("genParticles", GPCollection);
   const GenParticleCollection & gpc = *(GPCollection.product());

   for(unsigned int i=0; i<gpc.size(); i++){
     int pdgid = gpc[i].pdgId();
     int pdid = abs(pdgid);
     if(((pdid/100)==5 || (pdid/1000)==5) && ((pdid%10)==1 || (pdid%10)==2 )) {
       //std::cout << "B found at " << i << " pdgid " << pdgid << std::endl;
///////////////////TEST FOR NON-PYTHIA
//        bool btake = true;
//        for(unsigned int dn=0; dn<gpc[i].numberOfDaughters(); dn++){
//       int pdidDaught = (*gpc[i].daughter(dn)).pdgId(); 
//       int apdidDaught = abs(pdidDaught);
//       std::cout << "mother " << pdgid << " daughter " << pdidDaught << std::endl;
//       if((apdidDaught/100)==5 || (apdidDaught/1000)==5) btake = false;
//        }
//        if(btake){
/////////////////////////////////////
//        Particle::LorentzVector lv = gpc[i].p4();

       //check for B hadron mothers of pgc[i] already in output 
       //collection (->remove them)
       if((*output).size()>0) removeMothers(gpc[i], *output);

       //std::cout << "from Bproducer " << (*gpc[i].daughter(0)).vx() << std::endl;



       SimBHadron sbh(gpc[i]); 
       sbh.decPosition.SetXYZ((*gpc[i].daughter(0)).vx(),(*gpc[i].daughter(0)).vy(),(*gpc[i].daughter(0)).vz());
       setMoDauQuantities(gpc[i], sbh);
       output->push_back(sbh);
     }
   } 
   //std::cout << "number of b " << output->size() << std::endl;
   std::auto_ptr<SimBHadronCollection> pOut(output);
   iEvent.put(pOut);

}

// ------------ method called once each job just before starting event loop  ------------
void 
MCBHadronProducer::beginJob()
{
  noEvent = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MCBHadronProducer::endJob() {
}


void 
MCBHadronProducer::removeMothers(const reco::Candidate &gp, SimBHadronCollection &out){
  for(SBHit it = out.begin(); it != out.end(); it++){
    if(gp.p4()==(*it).p4()) {
      //int t3 = (*it).pdgId();
      SBHit temp = it;
      out.erase(temp);
      it--;
      //std::cout << "removed " << t3 << " B hadron from output\n";
    }
  }
  if(gp.status() != 3){
    for(unsigned int m=0; m<gp.numberOfMothers(); m++) {
      removeMothers((*gp.mother(m)), out);
    }
  }
}
      
void 
MCBHadronProducer::setMoDauQuantities(reco::GenParticle gp, SimBHadron &bh){
   //This function is buggy, it make a mess on single top

  bh.quarkstatus = -1;
  bh.brotherstatus = -1;
  
  std::set<const reco::Candidate *> bquarks; 
  findbquarks(gp, bquarks);
  int id2 = gp.pdgId(); 
  bool bmes = false; 
  if(id2/1000==0) bmes = true;
  if(bquarks.size()>1){
    for(candsetit it = bquarks.begin(); it !=bquarks.end(); it++){
      int id1 = (*it)->pdgId();
      //same sign 
      if(id1/abs(id1)==id2/abs(id2)){
	if(bmes){
	  //remove
	  candsetit temp = it;
	  it--;
	  bquarks.erase(temp);
	}
      }
      else{
	if(!bmes){
	  //remove
	  candsetit temp = it;
	  it--;
	  bquarks.erase(temp);
	}
      }
    }
  }
  if(bquarks.size()!=1) {
    ////std::cout << "MORE THAN ONE B/BBAR!!!\n";
  }
  else{
    candsetit it = bquarks.begin();
    bh.quarkstatus = (*it)->status(); 
    for(unsigned int m=0; m<(*it)->numberOfMothers(); m++){
      bh.motherids.push_back((*it)->mother(m)->pdgId());
    }
    if((*it)->numberOfMothers()>0){
      for(unsigned int d=0; d<((*it)->mother())->numberOfDaughters(); d++){
	int id1 = ((*it)->mother())->daughter(d)->pdgId();
	if(abs(id1)==5 && ((*it)->mother())->daughter(d)->p4() != (*it)->p4()){
	  if(id1!=(*it)->pdgId()) bh.brotherstatus = ((*it)->mother())->daughter(d)->status();
	}
      }
    }
  }
  

}

void 
MCBHadronProducer::findbquarks(const reco::Candidate &gp, std::set<const reco::Candidate *> &set){
  for(unsigned int m=0; m<gp.numberOfMothers(); m++){
    if(abs((*gp.mother(m)).pdgId())==5) {
      //const reco::Candidate *temp = gp.mother(m);
      bool hasbmother = false;
      for(unsigned int m1=0; m1<(*gp.mother(m)).numberOfMothers(); m1++){
	if(abs((*(*gp.mother(m)).mother(m1)).pdgId())==5) hasbmother = true;
      }
      if(!hasbmother) set.insert(gp.mother(m));
    }
    if((*gp.mother(m)).status()!=3 || (*gp.mother(m)).pdgId()!=2212) findbquarks((*gp.mother(m)), set);
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(MCBHadronProducer);
