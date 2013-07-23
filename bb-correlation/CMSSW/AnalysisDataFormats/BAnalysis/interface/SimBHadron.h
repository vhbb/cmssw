#ifndef SimBHadron_H
#define SimBHadron_H
//#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include <vector>
#include "DataFormats/Math/interface/Vector3D.h"

//#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
//#include "SimDataFormats/Track/interface/SimTrackContainer.h"
//include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/LorentzVector.h"
//#include "DataFormats/Math/interface/deltaR.h"
//#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

//#include <TVector3.h>

#include "AnalysisDataFormats/BAnalysis/interface/SimSecondaryVertex.h"

class SimBHadron : public reco::GenParticle{
public:

  SimBHadron(){ boolInfoSet_=false;}
  SimBHadron(const reco::GenParticle& gen) :
             reco::GenParticle(gen){ boolInfoSet_=false;}

  //TVector3 decPosition; 
  math::XYZVector decPosition; 

//  edm::Ptr<SimSecondaryVertex> ssv; 
//  double dist_ssv;
  
//  edm::Ptr<reco::Jet> jet;
//  double dR_jet;
  
  int quarkstatus; 
  int brotherstatus; 
  std::vector<int> motherids; 

  //set info
  void infoSet(bool val){boolInfoSet_=val;}
  void dist3D(double val){dist3D_ = val;}
  void dist2D(double val){dist2D_ = val;}
  void ncST(int val){ncST_ = val;}
  void addcST(){ncST_++;}
  void nRT(int val){nRT_ = val;}
  void nGoodRT(int val){nGoodRT_ = val;}
  void addRT(){nRT_++;}
  void addGoodRT(){nGoodRT_++;}
  void nmcVert(int val){nmcVert_=val;}
  void nmatmcVert(int val){nmatmcVert_=val;}
  void ncutsOkmcVert(int val){ncutsOkmcVert_=val;}
  void etaCutTracks(double val){etaCutTracks_ = val;}
  void ptCutTracks(double val){ptCutTracks_ = val;}
  void etaOk(bool val){etaOk_ = val;}
  void simOk(bool val){simOk_ = val;}
  void trackOk(bool val){trackOk_ = val;}
  void trackGoodOk(bool val){trackGoodOk_ = val;}
  void mcvertOk(bool val){mcvertOk_ = val;}
  void vertOk(bool val){vertOk_ = val;}
  void mccutsOk(bool val){mccutsOk_ = val;}
  void cutsOk(bool val){cutsOk_ = val;}
  void cutsOkNEW(bool val){cutsOkNEW_ = val;}
  void mcmatchOk(bool val){mcmatchOk_ = val;}
  void matchOk(bool val){matchOk_ = val;}
  void mccutsmatchOk(bool val){mccutsmatchOk_ = val;}
  void cutsmatchOk(bool val){cutsmatchOk_ = val;}
  void cutsmatchOkNEW(bool val){cutsmatchOkNEW_ = val;}

  //get info 
  double dist3D() const {return dist3D_;} 
  double dist2D() const {return dist2D_;} 
  int ncST() const {return ncST_;} 
  int nRT() const {return nRT_;} 
  int nGoodRT() const {return nGoodRT_;} 
  int nmcVert() const {return nmcVert_;} 
  int nmatmcVert() const {return nmatmcVert_;} 
  int ncutsOkmcVert() const {return ncutsOkmcVert_;} 
  double etaCutTracks() const {return etaCutTracks_;} 
  double ptCutTracks() const {return ptCutTracks_;} 
  bool infoSet() const {return boolInfoSet_;} 
  bool etaOk() const {return etaOk_;} 
  bool simOk() const {return simOk_;} 
  bool trackOk() const {return trackOk_;} 
  bool trackGoodOk() const {return trackGoodOk_;} 
  bool mcvertOk() const {return mcvertOk_;} 
  bool vertOk() const {return vertOk_;} 
  bool mccutsOk() const {return mccutsOk_;} 
  bool cutsOk() const {return cutsOk_;} 
  bool cutsOkNEW() const {return cutsOkNEW_;} 
  bool mcmatchOk() const {return mcmatchOk_;} 
  bool matchOk() const {return matchOk_;} 
  bool mccutsmatchOk() const {return mccutsmatchOk_;} 
  bool cutsmatchOk() const {return cutsmatchOk_;} 
  bool cutsmatchOkNEW() const {return cutsmatchOkNEW_;} 

private:
  double dist3D_; 
  double dist2D_; 
  int ncST_;
  int nRT_;
  int nGoodRT_;
  int nmcVert_;
  int nmatmcVert_;
  int ncutsOkmcVert_;
  double etaCutTracks_;
  double ptCutTracks_; 
  bool boolInfoSet_; 
  bool etaOk_;
  bool simOk_; 
  bool trackOk_; 
  bool trackGoodOk_; 
  bool mcvertOk_;
  bool vertOk_;
  bool mccutsOk_; //can not be filled because in Step2 (EDM) no cuts defined *
  bool cutsOk_; //can not be filled because in Step2 (EDM) no cuts defined *
  bool cutsOkNEW_; //can not be filled because in Step2 (EDM) no cuts defined *
  bool mcmatchOk_; //not filled at the moment (only important on Event level -> filled directly) *
  bool matchOk_; //not filled at the moment (only important on Event level -> filled directly) *
  bool mccutsmatchOk_; //not filled at the moment (only important on Event level -> filled directly) *
  bool cutsmatchOk_; //not filled at the moment (only important on Event level -> filled directly) *
  bool cutsmatchOkNEW_; //not filled at the moment (only important on Event level -> filled directly) *
  //* filled in Step 3

};

typedef  std::vector<SimBHadron> SimBHadronCollection;
typedef  edm::Ref<SimBHadronCollection> SimBHadronRef;
typedef  edm::RefProd<SimBHadronCollection> SimBHadronRefProd;
typedef  edm::RefVector<SimBHadronCollection> SimBHadronRefVector;

#endif

