// -*- C++ -*-
//
// Package:    MCSecVertexProducer
// Class:      MCSecVertexProducer
// 
/**\class MCSecVertexProducer MCSecVertexProducer.cc BAnalsysis/MCSecVertexProducer/src/MCSecVertexProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Andrea RIZZI
//         Created:  Wed Dec  3 12:14:27 CET 2008
// $Id: MCSecVertexProducer.cc,v 1.8 2010/04/16 08:58:51 wehrlilu Exp $
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

#include "AnalysisDataFormats/BAnalysis/interface/SimSecondaryVertex.h"
//
// class decleration
//


class MCSecVertexProducer : public edm::EDProducer {
   public:
      explicit MCSecVertexProducer(const edm::ParameterSet&);
  ~MCSecVertexProducer();
  
   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob();
  
      int vecInJet(reco::CaloJetCollection jc, TVector3 v, double &dR);

      // ----------member data ---------------------------

      edm::InputTag jetsrc_;


};


//constants, enums and typedefs

typedef std::set<int>::const_iterator intsetit;


//static data member definitions



//constructors and destructor

MCSecVertexProducer::MCSecVertexProducer(const edm::ParameterSet& iConfig) :
  jetsrc_(iConfig.getParameter<edm::InputTag>("jetsrc"))

{
 produces<SimSecondaryVertexCollection>(); 
 produces<SimPrimaryVertex>();

   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
*/
   //now do what ever other initialization is needed

}


MCSecVertexProducer::~MCSecVertexProducer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MCSecVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
*/ 

   /////////////////////////////////////////////
   //--------------GenParticles
   //Find simulated primary vertex
   Handle<GenParticleCollection> GPCollection; 
   iEvent.getByLabel("genParticles", GPCollection);
   const GenParticleCollection & gpc = *(GPCollection.product());
   Vertex::Point simpvpos;
   bool secvfound = false; 
   Vertex::Point p0; p0.SetXYZ(0,0,0);
   for(unsigned int i=0; i<gpc.size(); i++){
     //simulated primary vertex: 
     if(!secvfound && (gpc[i].vertex()!=p0)){
       simpvpos = gpc[i].vertex(); 
       secvfound = true;
     }
   }
   SimPrimaryVertex * spv = new SimPrimaryVertex();
   (*spv).position.SetXYZ(simpvpos.X(), simpvpos.Y(), simpvpos.Z());
   /////////////////////////////////////////////


   /////////////////////////////////////////////
   //--------------CALOJETS-----------------------
   //Collection of CaloJets
   Handle<CaloJetCollection> CJCollection;
   iEvent.getByLabel(jetsrc_, CJCollection);
   const CaloJetCollection cjc = *(CJCollection.product());
   /////////////////////////////////////////////

   /////////////////////////////////////////////
   //------------------Secondary Vertices (MC)-----------------------
   Handle<SimVertexContainer> SVC; 
   iEvent.getByLabel("g4SimHits", SVC); 
   const SimVertexContainer simverts = *(SVC.product());
   std::set<int> simVertAfterCuts;
   //std::vector<TVector3> simSVdir, simSVpos;
   /////////////////////////////////////////////

   /////////////////////////////////////////////
   //------------------SimTracks (MC)-----------------------
   //std::set<int> simVwithctracks; 
   Handle<SimTrackContainer> STC; 
   iEvent.getByLabel("g4SimHits", STC); 
   const SimTrackContainer simtracks = *(STC.product());
   //std::cout << "start simtracks\n";
   for(unsigned int i=0; i<simtracks.size(); i++){
     //insert all simsecvertices with at least one charged track
     if(fabs(simtracks[i].charge())>0) {
       simVertAfterCuts.insert(simtracks[i].vertIndex());
     }
   }
   //std::cout << "end simtracks\n";
   /////////////////////////////////////////////


   ////
   //remove all simvertices with distance (in transv. plane!) to simprimavy vertex > 10cm //was 5cm
   //and all with zcoordinate of distance > 30cm
   for(intsetit it=simVertAfterCuts.begin(); it!=simVertAfterCuts.end(); it++){
     int i = (*it);
     math::XYZTLorentzVectorD pos = simverts[i].position();
     TVector3 displ;
     displ.SetXYZ(pos.x()-simpvpos.x(),pos.y()-simpvpos.y(),pos.z()-simpvpos.z());
     //if(displ.Mag()>10){
     if((displ.x()*displ.x()+displ.y()*displ.y())>100 || displ.z()>30){
       simVertAfterCuts.erase(*it);
     }

   }
   //reduce number of simsecvert by removing all sv too close to pv 
   //and grouping sv together
   for(intsetit it=simVertAfterCuts.begin(); it!=simVertAfterCuts.end(); it++){
     math::XYZTLorentzVectorD pos = simverts[*it].position();
     //too near to primary vertex
     TVector3 v2pv; 
     v2pv.SetXYZ(pos.x()-simpvpos.x(),pos.y()-simpvpos.y(),pos.z()-simpvpos.z());
     if(v2pv.Mag()<0.001){
       simVertAfterCuts.erase(*it);
       continue; 
     }

     //erase other vectors too near to that one
     intsetit it2=it;
     for(it2++; it2!=simVertAfterCuts.end(); it2++){
       math::XYZTLorentzVectorD pos2 = simverts[*it2].position();
       TVector3 distance;
       distance.SetXYZ(pos2.x()-pos.x(),pos2.y()-pos.y(),pos2.z()-pos.z());

       if(distance.Mag()<0.001) {
	 simVertAfterCuts.erase(*it2);
	 it2--;
       }
     }
   }

   //output
   SimSecondaryVertexCollection * output = new SimSecondaryVertexCollection();
   for(intsetit it=simVertAfterCuts.begin(); it!=simVertAfterCuts.end(); it++){
     SimSecondaryVertex ssv;
     math::XYZTLorentzVectorD p = simverts[(*it)].position();
     TVector3 temp; 
     temp.SetXYZ(p.x(), p.y(), p.z());
     ssv.position = temp;
     temp.SetXYZ(p.x()-simpvpos.x(), p.y()-simpvpos.y(), p.z()-simpvpos.z());
     ssv.direction = temp;  
     //double dR;
     ssv.jet = edm::Ptr<reco::Jet>(CJCollection,vecInJet(cjc,temp,ssv.dR_jet));
//      for(int i=0; i<5; i++){
//      edm::Ptr<reco::Jet> j = edm::Ptr<reco::Jet>(CJCollection,i);
     
//      std::cout << "JETpt " << cjc[i].pt() << " is " << (*(j.get())).pt() << "\n";
//    }

     output->push_back(ssv);
   }


   
//    //output
//    SimSecondaryVertexCollection * output = new SimSecondaryVertexCollection();
//    SimSecondaryVertex v;
//    output->push_back(v); 
 
   std::auto_ptr<SimPrimaryVertex> pOut1(spv);
   iEvent.put(pOut1);

   std::auto_ptr<SimSecondaryVertexCollection> pOut(output);
   iEvent.put(pOut);

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/

}

// ------------ method called once each job just before starting event loop  ------------
void 
MCSecVertexProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MCSecVertexProducer::endJob() {
}

int
MCSecVertexProducer::vecInJet(reco::CaloJetCollection jc, TVector3 v, double &dR){
  dR = 100;
  double temp;
  int returnIndex = -1; 
  for(unsigned int i=0; i<jc.size(); i++){
    reco::Particle::LorentzVector lv = jc[i].p4();
    //momentum of jet
    TVector3 vjet;
    vjet.SetPtEtaPhi(lv.pt(), lv.eta(), lv.phi());
    temp = vjet.DeltaR(v);
    if(temp<dR) {
      returnIndex = i;
      dR = temp;
    }
  }
  return returnIndex;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MCSecVertexProducer);
