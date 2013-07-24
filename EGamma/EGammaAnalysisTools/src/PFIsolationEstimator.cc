#include <TFile.h>
#include "EGamma/EGammaAnalysisTools/interface/PFIsolationEstimator.h"
#include <cmath>
#include "DataFormats/Math/interface/deltaR.h"
using namespace std;

#ifndef STANDALONE
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#endif

using namespace reco;



//--------------------------------------------------------------------------------------------------
PFIsolationEstimator::PFIsolationEstimator() :
fisInitialized(kFALSE)
{
  // Constructor.
}



//--------------------------------------------------------------------------------------------------
PFIsolationEstimator::~PFIsolationEstimator()
{

}

//--------------------------------------------------------------------------------------------------
void PFIsolationEstimator::initialize( Bool_t  bApplyVeto, int iParticleType ) {

  setParticleType(iParticleType);

  //By default check for an option vertex association
  checkClosestZVertex = kTRUE;
  
  //Apply vetoes
  setApplyVeto(bApplyVeto);
  
  setDeltaRVetoBarrelPhotons();
  setDeltaRVetoBarrelNeutrals();
  setDeltaRVetoBarrelCharged();
  setDeltaRVetoEndcapPhotons();
  setDeltaRVetoEndcapNeutrals();
  setDeltaRVetoEndcapCharged();

  
  setRectangleDeltaPhiVetoBarrelPhotons();
  setRectangleDeltaPhiVetoBarrelNeutrals();
  setRectangleDeltaPhiVetoBarrelCharged();
  setRectangleDeltaPhiVetoEndcapPhotons();
  setRectangleDeltaPhiVetoEndcapNeutrals();
  setRectangleDeltaPhiVetoEndcapCharged();
  

  setRectangleDeltaEtaVetoBarrelPhotons();
  setRectangleDeltaEtaVetoBarrelNeutrals();
  setRectangleDeltaEtaVetoBarrelCharged();
  setRectangleDeltaEtaVetoEndcapPhotons();
  setRectangleDeltaEtaVetoEndcapNeutrals();
  setRectangleDeltaEtaVetoEndcapCharged();


  if(bApplyVeto && iParticleType==kElectron){
    //Setup veto conditions for electrons
    setDeltaRVetoBarrel(kFALSE);
    setDeltaRVetoEndcap(kTRUE);
    setRectangleVetoBarrel(kFALSE);
    setRectangleVetoEndcap(kFALSE);
    
    //Current recommended default value for the electrons
    setDeltaRVetoEndcapPhotons(0.08);
    setDeltaRVetoEndcapCharged(0.015);


  }else{
    //Setup veto conditions for photons
    setDeltaRVetoBarrel(kFALSE);
    setDeltaRVetoEndcap(kFALSE);
    setRectangleVetoBarrel(kFALSE);
    setRectangleVetoEndcap(kFALSE);

  }


}

//--------------------------------------------------------------------------------------------------
void PFIsolationEstimator::initializeElectronIsolation( Bool_t  bApplyVeto,float  fConeSize){
  initialize(bApplyVeto,kElectron);
  initializeRings(1, fConeSize);

}

//--------------------------------------------------------------------------------------------------
void PFIsolationEstimator::initializePhotonIsolation( Bool_t  bApplyVeto, float fConeSize ){
  initialize(bApplyVeto,kPhoton);
  initializeRings(1, fConeSize);
}


//--------------------------------------------------------------------------------------------------
void PFIsolationEstimator::initializeElectronIsolationInRings( Bool_t  bApplyVeto, int iNumberOfRings, float fRingSize ){
  initialize(bApplyVeto,kElectron);
  initializeRings(iNumberOfRings, fRingSize);
}

//--------------------------------------------------------------------------------------------------
void PFIsolationEstimator::initializePhotonIsolationInRings( Bool_t  bApplyVeto, int iNumberOfRings, float fRingSize  ){
  initialize(bApplyVeto,kPhoton);
  initializeRings(iNumberOfRings, fRingSize);
}


//--------------------------------------------------------------------------------------------------
void PFIsolationEstimator::initializeRings(int iNumberOfRings, float fRingSize){
 
  setRingSize(fRingSize);
  setNumbersOfRings(iNumberOfRings);
 
  fIsolationInRings.clear();
  for(int isoBin =0;isoBin<iNumberOfRings;isoBin++){
    float fTemp = 0.0;
    fIsolationInRings.push_back(fTemp);
    
    float fTempPhoton = 0.0;
    fIsolationInRingsPhoton.push_back(fTempPhoton);

    float fTempNeutral = 0.0;
    fIsolationInRingsNeutral.push_back(fTempNeutral);

    float fTempCharged = 0.0;
    fIsolationInRingsCharged.push_back(fTempCharged);

    float fTempChargedAll = 0.0;
    fIsolationInRingsChargedAll.push_back(fTempChargedAll);

  }

  fConeSize = fRingSize * (float)iNumberOfRings;

}
  
 
//--------------------------------------------------------------------------------------------------
float PFIsolationEstimator::fGetIsolation(const reco::PFCandidate * pfCandidate, const reco::PFCandidateCollection* pfParticlesColl, reco::Vertex& vtx, edm::Handle< reco::VertexCollection > vertices) {
 
  fGetIsolationInRings( pfCandidate, pfParticlesColl, vtx, vertices);
  fIsolation = fIsolationInRings[0];
  
  return fIsolation;
}


//--------------------------------------------------------------------------------------------------
vector<float >  PFIsolationEstimator::fGetIsolationInRings(const reco::PFCandidate * pfCandidate, const reco::PFCandidateCollection* pfParticlesColl, reco::Vertex& vtx, edm::Handle< reco::VertexCollection > vertices) {

  int isoBin;

  
  for(isoBin =0;isoBin<iNumberOfRings;isoBin++){
    fIsolationInRings[isoBin]=0.;
    fIsolationInRingsPhoton[isoBin]=0.;
    fIsolationInRingsNeutral[isoBin]=0.;
    fIsolationInRingsCharged[isoBin]=0.;
    fIsolationInRingsChargedAll[isoBin]=0.;
  }
  

  fEta =  pfCandidate->eta();
  fPhi =  pfCandidate->phi();
  fPt =  pfCandidate->pt();
  fVx =  pfCandidate->vx();
  fVy =  pfCandidate->vy();
  fVz =  pfCandidate->vz();

  for(unsigned iPF=0; iPF<pfParticlesColl->size(); iPF++) {

    const reco::PFCandidate& pfParticle= (*pfParticlesColl)[iPF]; 

    if(&pfParticle==(pfCandidate))
      continue;

    if(pfParticle.pdgId()==22){
      
      if(isPhotonParticleVetoed( &pfParticle)>=0.){
	isoBin = (int)(fDeltaR/fRingSize);
	fIsolationInRingsPhoton[isoBin]  = fIsolationInRingsPhoton[isoBin] + pfParticle.pt();
      }
      
    }else if(abs(pfParticle.pdgId())==130){
        
      if(isNeutralParticleVetoed(  &pfParticle)>=0.){
       	isoBin = (int)(fDeltaR/fRingSize);
	fIsolationInRingsNeutral[isoBin]  = fIsolationInRingsNeutral[isoBin] + pfParticle.pt();
      }
    

      //}else if(abs(pfParticle.pdgId()) == 11 ||abs(pfParticle.pdgId()) == 13 || abs(pfParticle.pdgId()) == 211){
    }else if(abs(pfParticle.pdgId()) == 211){
      if(isChargedParticleVetoed( &pfParticle, vtx, vertices)>=0.){
	isoBin = (int)(fDeltaR/fRingSize);
	fIsolationInRingsCharged[isoBin]  = fIsolationInRingsCharged[isoBin] + pfParticle.pt();
      }

    }
  }

 
  for(int isoBin =0;isoBin<iNumberOfRings;isoBin++){
    fIsolationInRings[isoBin]= fIsolationInRingsPhoton[isoBin]+ fIsolationInRingsNeutral[isoBin] +  fIsolationInRingsCharged[isoBin];
  }

  return fIsolationInRings;
}


//--------------------------------------------------------------------------------------------------
float PFIsolationEstimator::fGetIsolation(const reco::Photon * photon, const reco::PFCandidateCollection* pfParticlesColl, reco::Vertex& vtx, edm::Handle< reco::VertexCollection > vertices) {
 
  fGetIsolationInRings( photon, pfParticlesColl, vtx, vertices);
  fIsolation = fIsolationInRings[0];
  
  return fIsolation;
}


//--------------------------------------------------------------------------------------------------
vector<float >  PFIsolationEstimator::fGetIsolationInRings(const reco::Photon * photon, const reco::PFCandidateCollection* pfParticlesColl, reco::Vertex& vtx, edm::Handle< reco::VertexCollection > vertices) {

  int isoBin;
  
  for(isoBin =0;isoBin<iNumberOfRings;isoBin++){
    fIsolationInRings[isoBin]=0.;
    fIsolationInRingsPhoton[isoBin]=0.;
    fIsolationInRingsNeutral[isoBin]=0.;
    fIsolationInRingsCharged[isoBin]=0.;
    fIsolationInRingsChargedAll[isoBin]=0.;
  }
  
  int iMatch =  matchPFObject(photon,pfParticlesColl);


  fEta =  photon->eta();
  fPhi =  photon->phi();
  fPt =  photon->pt();
  fVx =  photon->vx();
  fVy =  photon->vy();
  fVz =  photon->vz();

  
  for(unsigned iPF=0; iPF<pfParticlesColl->size(); iPF++) {

    const reco::PFCandidate& pfParticle= (*pfParticlesColl)[iPF]; 
    

    if(iMatch == (int)iPF)
      continue;

    if(pfParticle.pdgId()==22){
    
      if(isPhotonParticleVetoed(&pfParticle)>=0.){
	isoBin = (int)(fDeltaR/fRingSize);
	fIsolationInRingsPhoton[isoBin]  = fIsolationInRingsPhoton[isoBin] + pfParticle.pt();

      }
      
    }else if(abs(pfParticle.pdgId())==130){
        
      if(isNeutralParticleVetoed( &pfParticle)>=0.){
       	isoBin = (int)(fDeltaR/fRingSize);
	fIsolationInRingsNeutral[isoBin]  = fIsolationInRingsNeutral[isoBin] + pfParticle.pt();
      }

      //}else if(abs(pfParticle.pdgId()) == 11 ||abs(pfParticle.pdgId()) == 13 || abs(pfParticle.pdgId()) == 211){
    }else if(abs(pfParticle.pdgId()) == 211){
      if(isChargedParticleVetoed(  &pfParticle, vtx, vertices)>=0.){
	isoBin = (int)(fDeltaR/fRingSize);
	fIsolationInRingsCharged[isoBin]  = fIsolationInRingsCharged[isoBin] + pfParticle.pt();
      }

    }
  }

 
  for(int isoBin =0;isoBin<iNumberOfRings;isoBin++){
    fIsolationInRings[isoBin]= fIsolationInRingsPhoton[isoBin]+ fIsolationInRingsNeutral[isoBin] +  fIsolationInRingsCharged[isoBin];
    }
  
  return fIsolationInRings;
}


//--------------------------------------------------------------------------------------------------
float  PFIsolationEstimator::isPhotonParticleVetoed( const reco::PFCandidate* pfIsoCand ){
  
  
  fDeltaR = deltaR(fEta,fPhi,pfIsoCand->eta(),pfIsoCand->phi()); 

  if(fDeltaR > fConeSize)
    return -999.;
  
  fDeltaPhi = deltaPhi(fPhi,pfIsoCand->phi()); 
  fDeltaEta = fEta-pfIsoCand->eta(); 

  if(!bApplyVeto)
    return fDeltaR;
 

  if(abs(fEta)<1.44442){
    if(bDeltaRVetoBarrel){
      if(fDeltaR < fDeltaRVetoBarrelPhotons)
        return -999.;
    }
    
    if(bRectangleVetoBarrel){
      if(abs(fDeltaEta) < fRectangleDeltaEtaVetoBarrelPhotons && abs(fDeltaPhi) < fRectangleDeltaPhiVetoBarrelPhotons){
	return -999.;
      }
    }
  }else{
    if(bDeltaRVetoEndcap){
      if(fDeltaR < fDeltaRVetoEndcapPhotons)
	return -999.;
    }
    if(bRectangleVetoEndcap){
      if(abs(fDeltaEta) < fRectangleDeltaEtaVetoEndcapPhotons && abs(fDeltaPhi) < fRectangleDeltaPhiVetoEndcapPhotons){
	 return -999.;
      }
    }
  }

  return fDeltaR;
}

//--------------------------------------------------------------------------------------------------
float  PFIsolationEstimator::isNeutralParticleVetoed( const reco::PFCandidate* pfIsoCand ){

  fDeltaR = deltaR(fEta,fPhi,pfIsoCand->eta(),pfIsoCand->phi()); 
  
  if(fDeltaR > fConeSize)
    return -999;
  
  fDeltaPhi = deltaPhi(fPhi,pfIsoCand->phi()); 
  fDeltaEta = fEta-pfIsoCand->eta(); 

  if(!bApplyVeto)
    return fDeltaR;

  if(abs(fEta)<1.44442){
    if(!bDeltaRVetoBarrel&&!bRectangleVetoBarrel){
      return fDeltaR;
    }
    
    if(bDeltaRVetoBarrel){
	if(fDeltaR < fDeltaRVetoBarrelNeutrals)
	  return -999.;
      }
      if(bRectangleVetoBarrel){
	if(abs(fDeltaEta) < fRectangleDeltaEtaVetoBarrelNeutrals && abs(fDeltaPhi) < fRectangleDeltaPhiVetoBarrelNeutrals){
	    return -999.;
	}
      }
      
    }else{
     if(!bDeltaRVetoEndcap&&!bRectangleVetoEndcap){
       return fDeltaR;
     }
      if(bDeltaRVetoEndcap){
	if(fDeltaR < fDeltaRVetoEndcapNeutrals)
	  return -999.;
      }
      if(bRectangleVetoEndcap){
	if(abs(fDeltaEta) < fRectangleDeltaEtaVetoEndcapNeutrals && abs(fDeltaPhi) < fRectangleDeltaPhiVetoEndcapNeutrals){
	  return -999.;
	}
      }
  }

  return fDeltaR;
}


//----------------------------------------------------------------------------------------------------
float  PFIsolationEstimator::isChargedParticleVetoed(const reco::PFCandidate* pfIsoCand, edm::Handle< reco::VertexCollection >  vertices  ){
  //need code to handle special conditions
  
  return -999;
}

//-----------------------------------------------------------------------------------------------------
float  PFIsolationEstimator::isChargedParticleVetoed(const reco::PFCandidate* pfIsoCand, reco::Vertex& vtxMain, edm::Handle< reco::VertexCollection >  vertices  ){
  

  
  
  VertexRef vtx = chargedHadronVertex(vertices,  *pfIsoCand );
  if(vtx.isNull())
    return -999.;
  
  float fVtxMainX = vtxMain.x();
  float fVtxMainY = vtxMain.y();
  float fVtxMainZ = vtxMain.z();

  if(iParticleType==kPhoton){

    
    //this piece of code does not use the chargedhadronvertex function. 
    /*  float dz = fabs(pfIsoCand->vz() - fVtxMainZ);
    if (dz > 1.)
      return -999.;
    
    double dxy = ( -(pfIsoCand->vx() - fVtxMainX)*pfIsoCand->py() + (pfIsoCand->vy() - fVtxMainY)*pfIsoCand->px()) / pfIsoCand->pt();

    if(fabs(dxy) > 0.1)
      return -999.;
    */
    
    float dz = fabs(vtx->z() - fVtxMainZ);
    if (dz > 1.)
      return -999.;
    
   
    double dxy = ( -(vtx->x() - fVtxMainX)*pfIsoCand->py() + (vtx->y() - fVtxMainY)*pfIsoCand->px()) / pfIsoCand->pt();

    if(fabs(dxy) > 0.2)
      return -999.;
    
  }else{
  

    float dz = fabs(vtx->z() - fVtxMainZ);
    if (dz > 1.)
      return -999.;
    
    double dxy = ( -(vtx->x() - fVx)*pfIsoCand->py() + (vtx->y() - fVy)*pfIsoCand->px()) / pfIsoCand->pt();
    if(fabs(dxy) > 0.1)
      return -999.;
  }
    
  fDeltaR = deltaR(pfIsoCand->eta(),pfIsoCand->phi(),fEta,fPhi); 

  if(fDeltaR > fConeSize)
    return -999.;

  fDeltaPhi = deltaPhi(fPhi,pfIsoCand->phi()); 
  fDeltaEta = fEta-pfIsoCand->eta(); 
  
  if(!bApplyVeto)
    return fDeltaR;  
  
  if(abs(fEta)<1.44442){
    if(!bDeltaRVetoBarrel&&!bRectangleVetoBarrel){
      return fDeltaR;
    }
    
    if(bDeltaRVetoBarrel){
	if(fDeltaR < fDeltaRVetoBarrelCharged)
	  return -999.;
      }
      if(bRectangleVetoBarrel){
	if(abs(fDeltaEta) < fRectangleDeltaEtaVetoBarrelCharged && abs(fDeltaPhi) < fRectangleDeltaPhiVetoBarrelCharged){
	    return -999.;
	}
      }
      
    }else{
     if(!bDeltaRVetoEndcap&&!bRectangleVetoEndcap){
       return fDeltaR;
     }
      if(bDeltaRVetoEndcap){
	if(fDeltaR < fDeltaRVetoEndcapCharged)
	  return -999.;
      }
      if(bRectangleVetoEndcap){
	if(abs(fDeltaEta) < fRectangleDeltaEtaVetoEndcapCharged && abs(fDeltaPhi) < fRectangleDeltaPhiVetoEndcapCharged){
	  return -999.;
	}
      }
  }
		   
  

  return fDeltaR;
}


//--------------------------------------------------------------------------------------------------
 VertexRef  PFIsolationEstimator::chargedHadronVertex(  edm::Handle< reco::VertexCollection > verticesColl, const reco::PFCandidate& pfcand ){

  //code copied from Florian's PFNoPU class
    
  reco::TrackBaseRef trackBaseRef( pfcand.trackRef() );

  size_t  iVertex = 0;
  unsigned index=0;
  unsigned nFoundVertex = 0;

  float bestweight=0;
  
  const reco::VertexCollection& vertices = *(verticesColl.product());

  for( reco::VertexCollection::const_iterator iv=vertices.begin(); iv!=vertices.end(); ++iv, ++index) {
    
    const reco::Vertex& vtx = *iv;
    
    // loop on tracks in vertices
    for(reco::Vertex::trackRef_iterator iTrack=vtx.tracks_begin();iTrack!=vtx.tracks_end(); ++iTrack) {

      const reco::TrackBaseRef& baseRef = *iTrack;

      // one of the tracks in the vertex is the same as 
      // the track considered in the function
      if(baseRef == trackBaseRef ) {
        float w = vtx.trackWeight(baseRef);
        //select the vertex for which the track has the highest weight
        if (w > bestweight){
          bestweight=w;
          iVertex=index;
          nFoundVertex++;
        }
      }
    }
    
  }
 
 
  
  if (nFoundVertex>0){
    if (nFoundVertex!=1)
      edm::LogWarning("TrackOnTwoVertex")<<"a track is shared by at least two verteces. Used to be an assert";
    return  VertexRef( verticesColl, iVertex);
  }
  // no vertex found with this track. 

  // optional: as a secondary solution, associate the closest vertex in z
  if ( checkClosestZVertex ) {

    double dzmin = 10000.;
    double ztrack = pfcand.vertex().z();
    bool foundVertex = false;
    index = 0;
    for( reco::VertexCollection::const_iterator  iv=vertices.begin(); iv!=vertices.end(); ++iv, ++index) {

      double dz = fabs(ztrack - iv->z());
      if(dz<dzmin) {
        dzmin = dz;
        iVertex = index;
        foundVertex = true;
      }
    }

    if( foundVertex ) 
      return  VertexRef( verticesColl, iVertex);  
  
  }
   
  return  VertexRef( );
}



int PFIsolationEstimator::matchPFObject(const reco::Photon* photon, const reco::PFCandidateCollection* Candidates ){
  
  Int_t iMatch = -1;

  int i=0;
  for(reco::PFCandidateCollection::const_iterator iPF=Candidates->begin();iPF !=Candidates->end();iPF++){
    const reco::PFCandidate& pfParticle = (*iPF);
    //    if((((pfParticle.pdgId()==22 && pfParticle.mva_nothing_gamma()>0.01) || TMath::Abs(pfParticle.pdgId())==11) )){
    if((((pfParticle.pdgId()==22 ) || TMath::Abs(pfParticle.pdgId())==11) )){
     
      if(pfParticle.superClusterRef()==photon->superCluster())
	iMatch= i;
     
    }
    
    i++;
  }
  
  if(iMatch == -1){
    i=0;
    float fPt = -1;
    for(reco::PFCandidateCollection::const_iterator iPF=Candidates->begin();iPF !=Candidates->end();iPF++){
      const reco::PFCandidate& pfParticle = (*iPF);
      if((((pfParticle.pdgId()==22 ) || TMath::Abs(pfParticle.pdgId())==11) )){
	if(pfParticle.pt()>fPt){
	  fDeltaR = deltaR(pfParticle.eta(),pfParticle.phi(),photon->eta(),photon->phi());
	  if(fDeltaR<0.1){
	    iMatch = i;
	    fPt = pfParticle.pt();
	  }
	}
      }
      i++;
    }
  }
  
  return iMatch;

}
