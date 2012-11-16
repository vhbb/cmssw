#ifndef VHbbCandidate__H 
#define VHbbCandidate__H 

#include <TLorentzVector.h>
#include <TVector2.h>
#include <vector>

#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"

class VHbbCandidate {
 public:
   //Zmumu = 0 
   //Zee = 1 
   //Wmun = 2
   //Wen = 3
   //Znn = 4

  enum CandidateType{Zmumu, Zee, Wmun, Wen, Znn, UNKNOWN};

    VHbbCandidate(){candidateType=UNKNOWN;}

  class VectorCandidate {
  public:
  VectorCandidate() : firstLepton(99999),secondLepton(99999) {}
    double Mt(CandidateType candidateType) const {
    if(candidateType==Wen)
      {
       float ptl=electrons[0].p4.Pt();
       float met=mets[0].p4.Pt();
       float et=ptl+met;
       return sqrt(et*et - p4.Pt()*p4.Pt()  );
      }
    if(candidateType==Wmun)
      {
       float ptl=muons[0].p4.Pt();
       float met=mets[0].p4.Pt();
       float et=ptl+met;
       return sqrt(et*et - p4.Pt()*p4.Pt()  );
      }
    return 0;
   }
    TLorentzVector p4;
    std::vector<VHbbEvent::MuonInfo> muons;
    std::vector<VHbbEvent::ElectronInfo> electrons;
    std::vector<VHbbEvent::TauInfo> taus;
    std::vector<VHbbEvent::METInfo> mets;
    
    unsigned int firstLepton,secondLepton;
    
  };
  
  class HiggsCandidate {
  public:
   TLorentzVector p4;
    std::vector<VHbbEvent::SimpleJet> jets;
    bool HiggsFlag;
    float deltaTheta;
    std::vector <float> helicities;
 public:
    VHbbEvent::SimpleJet& firstJet(){return jets[0];}
    VHbbEvent::SimpleJet& secondJet(){return jets[1];}
 };
  
  class FatHiggsCandidate {
  public:
   TLorentzVector p4;
    std::vector<VHbbEvent::SimpleJet> jets;
    bool FatHiggsFlag;
    int subjetsSize;
    float deltaTheta;
    std::vector <float> helicities;
 public:
    VHbbEvent::SimpleJet& firstJet(){return jets[0];}
    VHbbEvent::SimpleJet& secondJet(){return jets[1];}
 };
    
    // *****************************
    // added by Nhan
class RawFatHiggsCandidate {
public:
    
    TLorentzVector p4;
    bool RawFatHiggsFlag;
    
    float prunedJetMass;
    float trimmedJetMass;
    float filteredJetMass;
    float ungroomedJetMass;
    
    float pT_pr;
    float eta_pr;    
    float px_pr;        
    float py_pr;        
    float pz_pr;        
    float e_pr;        
    float tau1;
    float tau2;
    float tau3;
    float tau4;
    
    float sj1_px_pr;
    float sj1_py_pr;
    float sj1_pz_pr;
    float sj1_e_pr;    
    float sj2_px_pr;
    float sj2_py_pr;
    float sj2_pz_pr;
    float sj2_e_pr;
    
public:
    TLorentzVector getp4(){ return p4; }
};    
    // *****************************    
    
  void setCandidateType (CandidateType c){candidateType = c;}


  double deltaPhi() const {
   return V.p4.DeltaPhi(H.p4);
  }

  double Mt() const {
   return V.Mt(candidateType);
  }
  
 int additionalLeptons() const {
   int expectedLeptons = 0;
   if(  candidateType == Wmun ||  candidateType == Wen) expectedLeptons =1;
   if(  candidateType == Zmumu ||  candidateType == Zee) expectedLeptons =2;

   return ( V.muons.size() + V.electrons.size() - expectedLeptons);

 } 

 public:
  TLorentzVector p4(){return V.p4+H.p4;}
  CandidateType candidateType;
  HiggsCandidate H;
  FatHiggsCandidate FatH;
    RawFatHiggsCandidate RawFatH;
  VectorCandidate V;
  std::vector<VHbbEvent::SimpleJet> additionalJets;
  std::vector<VHbbEvent::SimpleJet> additionalJetsFat;
};


#endif
