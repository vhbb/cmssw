#ifndef SimSecVertex_H
#define SimSecVertex_H
//#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/Jet.h" 
#include <vector>
#include "DataFormats/Math/interface/Vector3D.h"

#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/LorentzVector.h"
//#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include <TVector3.h>

class SimPrimaryVertex {
public:
  SimPrimaryVertex(){}
  TVector3 position;
};

class SimSecondaryVertex {
public:

  SimSecondaryVertex(){}
  TVector3 direction;
  TVector3 position;
  //math::XYZVector direction; 
  //math::XYZVector position; 
  edm::Ptr<reco::Jet> jet; 
  double dR_jet;
  bool fromB;
 
};

typedef  std::vector<SimSecondaryVertex> SimSecondaryVertexCollection;
typedef  edm::Ref<SimSecondaryVertexCollection> SimSecondaryVertexRef;
typedef  edm::RefProd<SimSecondaryVertexCollection> SimSecondaryVertexRefProd;
typedef  edm::RefVector<SimSecondaryVertexCollection> SimSecondaryVertexRefVector;

#endif
