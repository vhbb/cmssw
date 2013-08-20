#ifndef RecoBVertexCand_H
#define RecoBVertexCand_H
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include <vector>
#include "DataFormats/Math/interface/Vector3D.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/LorentzVector.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

/* #include "AnalysisDataFormats/BAnalysis/interface/SimSecondaryVertex.h" */

class RecoBVertexCand : public reco::Vertex{
public:

  RecoBVertexCand(){}
  RecoBVertexCand(const reco::Vertex& vert) :
             reco::Vertex(vert){}

  //TVector3 decPosition; 
  //math::XYZVector decPosition; 

  //set info
    //void infoSet(bool val){boolInfoSet_=val;} 
    void invM(float val){invM_=val;}
    void gamma(float val){gamma_=val;}
    void eta(float val){eta_=val;}
    void pt(float val){pt_=val;}
    void noftr(int val){noftr_=val;}
    void noftrWok(int val){noftrWok_=val;}
    void dist2D(float val){dist2D_=val;}
    void distSig2D(float val){distSig2D_=val;}
    void dist3D(float val){dist3D_=val;}
    void dist3D_norm(float val){dist3D_norm_=val;}
    void distSig3D(float val){distSig3D_=val;}
    void cosFdTdsum(float val){cosFdTdsum_=val;}
    void normChi2(float val){normChi2_=val;}
    void coneSize(float val){coneSize_=val;}
    void coneSize_norm(float val){coneSize_norm_=val;}
    void dRtoClosestB(float val){dRtoClosestB_=val;}
    void dRtoMatchedB(float val){dRtoMatchedB_=val;}
    void massOfMatchedB(float val){massOfMatchedB_=val;}
    void ptOfMatchedB(float val){ptOfMatchedB_=val;}
    void nofMatchedB(int val){nofMatchedB_=val;}
    void processCode(int val){processCode_=val;}
    void eventNo(int val){eventNo_=val;}
    void runNo(int val){runNo_=val;}
    void matched(bool val){matched_=val;}
    void selected(bool val){selected_=val;}
    void unique(bool val){unique_=val;}
    void primV(bool val){primV_=val;}
    void secV(bool val){secV_=val;}
    void tertV(bool val){tertV_=val;}
    void BweakDec(bool val){BweakDec_=val;}
    void DweakDec(bool val){DweakDec_=val;}
    void LightDec(bool val){LightDec_=val;}
    void KsDec(bool val){KsDec_=val;}
    void LongLivedDec(bool val){LongLivedDec_=val;}
    void fakeV(bool val){fakeV_=val;}
    void flightDirection(GlobalVector val){flightDirection_=val;}
    void trackSumDirection(GlobalVector val){trackSumDirection_=val;}
    void trackSumDirectionWok(GlobalVector val){trackSumDirectionWok_=val;}

  //get info 
    float invM() const {return invM_;}
    float gamma() const {return gamma_;}
    float eta() const {return eta_;}
    float pt() const {return pt_;}
    int noftr() const {return noftr_;}
    int noftrWok() const {return noftrWok_;}
    float dist2D() const {return dist2D_;}
    float dist3D() const {return dist3D_;}
    float dist3D_norm() const {return dist3D_norm_;}
    float distSig2D() const {return distSig2D_;}
    float distSig3D() const {return distSig3D_;}
    float cosFdTdsum() const {return cosFdTdsum_;}
    float normChi2() const {return normChi2_;}
    float coneSize() const {return coneSize_;}
    float coneSize_norm() const {return coneSize_norm_;}
    float dRtoClosestB() const {return dRtoClosestB_;}
    float dRtoMatchedB() const {return dRtoMatchedB_;}
    float massOfMatchedB() const {return massOfMatchedB_;}
    float ptOfMatchedB() const {return ptOfMatchedB_;}
    int nofMatchedB() const {return nofMatchedB_;}
    int processCode() const {return processCode_;}
    int eventNo() const {return eventNo_;}
    int runNo() const {return runNo_;}
    bool matched() const {return matched_;}
    bool selected() const {return selected_;}
    bool unique() const {return unique_;}
    bool primV() const {return primV_;}
    bool secV() const {return secV_;}
    bool tertV() const {return tertV_;} 
    bool BweakDec() const {return BweakDec_;}
    bool DweakDec() const {return DweakDec_;}
    bool LightDec() const {return LightDec_;}
    bool KsDec() const {return KsDec_;}
    bool LongLivedDec() const {return LongLivedDec_;}
    bool fakeV() const {return fakeV_;}
    GlobalVector flightDirection() const {return flightDirection_;}
    GlobalVector trackSumDirection() const {return trackSumDirection_;}
    GlobalVector trackSumDirectionWok() const {return trackSumDirectionWok_;}

private:
    GlobalVector flightDirection_;
    GlobalVector trackSumDirection_;
    GlobalVector trackSumDirectionWok_;
    //GlobalPoint position; 

    float invM_; 
    float gamma_; 
    float eta_; 
/*     float ptOfMatchedB;  */
/*     float massOfMatchedB;  */
    int noftr_;
    int noftrWok_;
    float dist2D_; 
    float dist3D_; 
    float dist3D_norm_; 
    float distSig2D_; 
    float distSig3D_;
    float cosFdTdsum_; 
/*     float sigmaPhi; */
/*     float sigmaEta;  */
/*     float chi2; */
/*     float ndof;  */
    float normChi2_; // --> normalizedChi2() makes errors when plotting with fwl
    float coneSize_; 
    float coneSize_norm_; 
    float dRtoClosestB_; //MC
    float dRtoMatchedB_; //MC
    float massOfMatchedB_; //MC
    float ptOfMatchedB_; //MC
    int nofMatchedB_; //MC
    int processCode_; //MC
    int eventNo_; 
    int runNo_; 
/*     int nofSimB;  */
/*     int nofSelRSV;  */
    float pt_; 
    bool matched_; 
    bool selected_; 
    bool unique_; 

    //Victor
    bool primV_;
    bool secV_;
    bool tertV_; 
    bool BweakDec_;
    bool DweakDec_;
    bool LightDec_;
    bool KsDec_;
    bool LongLivedDec_;
    bool fakeV_;    

};

typedef  std::vector<RecoBVertexCand> RecoBVertexCandCollection;
typedef  edm::Ref<RecoBVertexCandCollection> RecoBVertexCandRef;
typedef  edm::RefProd<RecoBVertexCandCollection> RecoBVertexCandRefProd;
typedef  edm::RefVector<RecoBVertexCandCollection> RecoBVertexCandRefVector;

#endif

