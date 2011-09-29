// -*- C++ -*-
//
// Package:    BCandidateProducer
// Class:      BCandidateProducer
// 
/**\class BCandidateProducer BCandidateProducer.cc BAnalysis/BCandidateProducer/src/BCandidateProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lukas Wehrli (IPP/ETHZ) [wehrlilu]
//         Created:  Wed May 12 09:34:25 CEST 2010
// $Id: BCandidateProducer.cc,v 1.2 2010/11/25 09:25:27 wehrlilu Exp $
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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

#include <algorithm>

struct ownVertex{
  unsigned int intNo; 
  reco::Vertex vert; 
  double invm;
  unsigned int ntracks;
};

inline bool operator<(ownVertex v1, ownVertex v2){
  if(v1.ntracks>2 && v2.ntracks<3) return true;
  if(v1.ntracks<3 && v2.ntracks>2) return false;
  return (v1.invm>v2.invm);
}


class BCandidateProducer : public edm::EDProducer {
   public:
      explicit BCandidateProducer(const edm::ParameterSet&);
      ~BCandidateProducer();

   private:
      virtual void beginJob();
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob(); 


      bool testOnBCchain(unsigned int &bCand, unsigned int &i, std::vector<ownVertex> &ownVertices, reco::Vertex &pv,  math::XYZTLorentzVector &bCandmom, std::set<reco::TrackRef> &bCandtracks); 
      void addTracks(std::set<reco::TrackRef> &trackrefs, reco::Vertex &v, double weight=0.5);
      GlobalVector flightDirection(reco::Vertex &pv, reco::Vertex &sv);
      math::XYZTLorentzVector sum4momentaTracks(std::set<reco::TrackRef> &trackrefs);

      inline double deltaPhi(double phi1, double phi2) {
        double result = phi1 - phi2;
        while (result > M_PI) result -= 2.*M_PI;
        while (result <= -M_PI) result += 2.*M_PI;
        return result;
      }

      inline double deltaR2(GlobalVector &v1, GlobalVector &v2) {
        double deta = v1.eta() - v2.eta();
        double dphi = deltaPhi(v1.phi(), v2.phi());
        return deta*deta + dphi*dphi;
      }

      inline double deltaR(GlobalVector &v1, GlobalVector &v2) {
        return sqrt(deltaR2 (v1,v2));
      }
      inline double deltaR(math::XYZVector v1, GlobalVector v2) {
        GlobalVector v(v1.x(),v1.y(),v1.z());
        return sqrt(deltaR2 (v,v2));
      }
      inline double deltaR(math::XYZVector v1, math::XYZVector v2) {
        GlobalVector v(v1.x(),v1.y(),v1.z());
        GlobalVector w(v2.x(),v2.y(),v2.z());
        return sqrt(deltaR2 (v,w));
      }


      
      edm::InputTag vertexsrc_, offlinePrimaryVertices_;
      double minDRForUnique, vecSumIMCUTForUnique, minCosPAtomerge, maxPtreltomerge; 
      bool debug; 
};

BCandidateProducer::BCandidateProducer(const edm::ParameterSet& iConfig):
  vertexsrc_(iConfig.getParameter<edm::InputTag>("src")),
  offlinePrimaryVertices_(iConfig.getParameter<edm::InputTag>("primaryVertices")),
  minDRForUnique(iConfig.getUntrackedParameter<double>("minDRUnique",0.4)),
  vecSumIMCUTForUnique(iConfig.getUntrackedParameter<double>("minvecSumIMifsmallDRUnique",5.5)), 
  minCosPAtomerge(iConfig.getUntrackedParameter<double>("minCosPAtomerge",0.99)),
  maxPtreltomerge(iConfig.getUntrackedParameter<double>("maxPtreltomerge",7777.0))

{

  //produces<reco::VertexCollection>("Selected");
  produces<reco::VertexCollection>();
  //produces<reco::VertexCollection>("NotSelected");
   produces<std::vector<reco::LeafCandidate> >();
   produces<std::vector<int> >();
  
}


BCandidateProducer::~BCandidateProducer()
{
 
}


// ------------ method called to produce the data  ------------
void
BCandidateProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;

   VertexCollection *bvertices = new VertexCollection();
   VertexCollection *nonbvertices = new VertexCollection();
   std::vector<reco::LeafCandidate> *blc = new std::vector<reco::LeafCandidate>;
   

   ///////////////////////////////////////////
   //------------------Primary Vertices-----------------------
   edm::Handle<reco::VertexCollection> primaryVertexCollection;
   iEvent.getByLabel(offlinePrimaryVertices_, primaryVertexCollection);
   const reco::VertexCollection pvc = *(primaryVertexCollection.product());
   Vertex pv = pvc[0];
   //////////////////////////////////////////


   ///////////////////////////////////////////
   //------------------Secondary Vertices-----------------------
   Handle<std::vector<reco::Vertex> > SVC;
   iEvent.getByLabel(vertexsrc_, SVC);
   const std::vector<reco::Vertex> svc = *(SVC.product());
   ///////////////////////////////////////////
   
   std::vector<ownVertex> ownVertices; 
   for(unsigned int i=0; i<svc.size(); i++){

        reco::Vertex sv = svc[i];
	ownVertex ov = {i,sv,sv.p4(0.13957,0.5).M(),sv.tracksSize()};
	ownVertices.push_back(ov);
   }
   sort(ownVertices.begin(), ownVertices.end()); 
   
   if(debug){
     for(unsigned int i=0; i<ownVertices.size(); i++){
       std::cout << i << " ntracks " << ownVertices[i].ntracks << " invm " << ownVertices[i].invm << "\n";
     }
     std::cout << "\n";
   }

    
   std::vector<math::XYZTLorentzVectorD> bCandmoms; 
   math::XYZTLorentzVectorD bCandmom; 
   std::vector<std::set<TrackRef> > bCandtrackSets; 
   std::set<TrackRef> bCandtracks;
   
   unsigned int doLater = 0; 
   for(unsigned int bCand=0; bCand<ownVertices.size(); bCand++){
     bCandmom = ownVertices[bCand].vert.p4(0.13957,0.5);
     bCandtracks.clear(); 
     addTracks(bCandtracks, ownVertices[bCand].vert, 0.5); 
     for(unsigned int i=bCand+1; i<ownVertices.size(); i++){
       if(ownVertices[i].invm>1.5){
	 doLater++; 
       }
       else{
	 if(testOnBCchain(bCand, i, ownVertices, pv, bCandmom, bCandtracks)){
	   i--;
	 }
       }
     }
     for(unsigned int j=bCand+1; j<bCand+1+doLater; j++){
       if(testOnBCchain(bCand, j, ownVertices, pv, bCandmom, bCandtracks))
	 {
	   j--; doLater--;
	 }
     }
     doLater=0;
     bCandmoms.push_back(bCandmom); 
     bCandtrackSets.push_back(bCandtracks); 
   }


   if(ownVertices.size()!=bCandmoms.size()) std::cout << "SIZE MISMATCH\n\n\n";
   else if(debug) std::cout << "NUMBER BCANDID: " << ownVertices.size() << "\n";

   //store number of tracks of bcandidates as ints in vector
   std::auto_ptr<std::vector<int> > bcandNtracks( new std::vector<int>() );

   for(unsigned int v=0; v<ownVertices.size(); v++){
     if(bCandmoms[v].M() < ownVertices[v].invm) std::cout << "MASS ERROR";
    
     bvertices->push_back(ownVertices[v].vert); 
     Vertex::Point p(ownVertices[v].vert.x(),ownVertices[v].vert.y(),ownVertices[v].vert.z());
 
     double charge = 0.0;
     int ntracksBC = 0; 
     for(std::set<reco::TrackRef>::const_iterator it=bCandtrackSets[v].begin(); it!=bCandtrackSets[v].end(); it++){
       charge += (*it)->charge();
       ntracksBC++;  
     }

     LeafCandidate lc(charge, bCandmoms[v], p); 
     blc->push_back(lc); 
     bcandNtracks->push_back(ntracksBC); 
   }

 
   std::auto_ptr<VertexCollection> bvertColl(bvertices);
   std::auto_ptr<VertexCollection> nonbvertColl(nonbvertices);
   std::auto_ptr<std::vector<reco::LeafCandidate> > bcandColl(blc); 
   //iEvent.put(bvertColl,"Selected");
   iEvent.put(bvertColl);
   iEvent.put(bcandColl);
   iEvent.put(bcandNtracks); 


 
}

// ------------ method called once each job just before starting event loop  ------------
void 
BCandidateProducer::beginJob()
{
  debug = false;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BCandidateProducer::endJob() {
}

bool
BCandidateProducer::testOnBCchain(unsigned int &bCand, unsigned int &i, std::vector<ownVertex> &ownVertices, reco::Vertex &pv, math::XYZTLorentzVector &bCandmom, std::set<reco::TrackRef> &bCandtracks){
  bool debug = false; 
  std::set<reco::TrackRef> trackrefs;

  //add all tracks
  trackrefs.clear();
  for(std::set<reco::TrackRef>::const_iterator it=bCandtracks.begin(); it!=bCandtracks.end(); it++) 
    trackrefs.insert(*it);
  addTracks(trackrefs, ownVertices[i].vert, 0.5); 
  
  
  //compute pair inv mass and DR
  //GlobalVector gv1(bCandmom.X(),bCandmom.Y(),bCandmom.Z()); 
  GlobalVector gv1 = flightDirection(pv, ownVertices[bCand].vert);
  GlobalVector gv2 = flightDirection(pv, ownVertices[i].vert);
 
  double dR = deltaR(gv1,gv2); 
  if(debug) std::cout << "DR " << dR << std::endl;
  math::XYZTLorentzVector pairTrkSum = sum4momentaTracks(trackrefs);
  if(debug) std::cout << "InvM cand" << bCandmom.M() << " + add? " << ownVertices[i].invm << " pair " << pairTrkSum.M() << std::endl;
  
  //pointing angle and ptrel
  GlobalVector axis1 = GlobalVector(ownVertices[bCand].vert.p4().X(),ownVertices[bCand].vert.p4().Y(),ownVertices[bCand].vert.p4().Z());
  GlobalVector axis2 = GlobalVector(ownVertices[i].vert.p4().X(),ownVertices[i].vert.p4().Y(),ownVertices[i].vert.p4().Z());
  reco::SecondaryVertex sv1(pv,ownVertices[bCand].vert, axis1, true);
  reco::SecondaryVertex sv2(pv,ownVertices[i].vert, axis2, true);
  double cosPA, ptrel; 
  GlobalVector pvToNear; 
  if(sv1.dist3d().value()<sv2.dist3d().value()){
    GlobalVector nearToFar = flightDirection(ownVertices[bCand].vert,ownVertices[i].vert);
    cosPA = nearToFar.dot(axis2)/axis2.mag()/nearToFar.mag(); 
    pvToNear = flightDirection(pv, ownVertices[bCand].vert); 
    double cosa = pvToNear.dot(axis2)/pvToNear.mag()/axis2.mag(); 
    ptrel = sqrt(1.0 - cosa*cosa)*axis2.mag();  
  }
  else {
    GlobalVector nearToFar = flightDirection(ownVertices[i].vert,ownVertices[bCand].vert);
    cosPA = nearToFar.dot(axis1)/axis1.mag()/nearToFar.mag(); 
    pvToNear = flightDirection(pv, ownVertices[i].vert); 
    //ptrel = pvToNear.dot(axis1)/pvToNear.mag();
    double cosa = pvToNear.dot(axis1)/pvToNear.mag()/axis1.mag(); 
    ptrel = sqrt(1.0 - cosa*cosa)*axis1.mag();  
  }

  if (debug){
    std::cout << "pv pos " << pv.position().x() << " " << pv.position().y() << " " << pv.position().z() << std::endl;
    std::cout << "V1 pos " << ownVertices[bCand].vert.position().x() << " " 
	      << ownVertices[bCand].vert.position().y() << " " 
	      << ownVertices[bCand].vert.position().z() << std::endl;
    std::cout << "V1 mom " << axis1.x() << " " << axis1.y() << " " << axis1.z() << std::endl;
    
    std::cout << "V2 pos " << ownVertices[i].vert.position().x() << " " 
	      << ownVertices[i].vert.position().y() << " " 
	      << ownVertices[i].vert.position().z() << std::endl;
    std::cout << "V2 mom " << axis2.x() << " " << axis2.y() << " " << axis2.z() << std::endl;
  }
    

  if(debug) std::cout << "dR " << dR << " pair sum " << pairTrkSum.M() << " cosPA " << cosPA << " ptrel " << ptrel << std::endl;

  //decide on b-c-chain  
  if(dR<minDRForUnique && pairTrkSum.M()<vecSumIMCUTForUnique && cosPA>minCosPAtomerge && fabs(ptrel)<maxPtreltomerge) {
    if(debug) std::cout << "REMOVE\n";
    addTracks(bCandtracks, ownVertices[i].vert, 0.5); 
    bCandmom = sum4momentaTracks(bCandtracks);
    ownVertices.erase(ownVertices.begin()+i);
    //i--;
    return true;
  }
  return false;
}

void 
BCandidateProducer::addTracks(std::set<reco::TrackRef> &trackrefs, reco::Vertex &v, double weight){

  for(reco::Vertex::trackRef_iterator ti = v.tracks_begin(); ti!=v.tracks_end(); ti++){
    if(v.trackWeight(*ti)>weight){
      reco::TrackRef t = ti->castTo<reco::TrackRef>();
      trackrefs.insert(t);
    }

  }

}

GlobalVector
BCandidateProducer::flightDirection(reco::Vertex &pv, reco::Vertex &sv){
  GlobalVector res(sv.position().X() - pv.position().X(),
                    sv.position().Y() - pv.position().Y(),
                    sv.position().Z() - pv.position().Z());
  return res;
}

math::XYZTLorentzVector
BCandidateProducer::sum4momentaTracks(std::set<reco::TrackRef> &trackrefs){
  math::XYZTLorentzVector res(0,0,0,0);
  for(std::set<reco::TrackRef>::const_iterator it = trackrefs.begin(); it!= trackrefs.end(); it++){
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> >  temp;
    temp.SetPx((*it)->px()); 
    temp.SetPy((*it)->py()); 
    temp.SetPz((*it)->pz()); 
    temp.SetM(0.13957); 
    res += temp;
  }
  return res;
}


//define this as a plug-in
DEFINE_FWK_MODULE(BCandidateProducer);
