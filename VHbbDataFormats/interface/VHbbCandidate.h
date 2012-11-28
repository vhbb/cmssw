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
class FatHiggsCandidateFJ3 {
public:
    
    TLorentzVector p4;
    bool FatHiggsFJ3Flag;
    
    TLorentzVector p4_pr;
    TLorentzVector p4_ft;
    TLorentzVector p4_tr;
    TLorentzVector area_pr;
    TLorentzVector area_ft;
    TLorentzVector area_tr;
    
    float tau1;
    float tau2;
    float tau3;
    float tau4;
    float qjetVol;
        
        // mdft matching
    bool matchedMDFTCandidate;
    TLorentzVector matchedMDFTCandidate_p4;
    std::vector< TLorentzVector > subjets_mdft_p4;
    std::vector< TLorentzVector > subjets_mdft_csv;

    std::vector< TLorentzVector > subjets_pr_p4;
    std::vector< TLorentzVector > subjets_pr_csv;

    std::vector< TLorentzVector > subjets_ft_p4;
    std::vector< TLorentzVector > subjets_ft_csv;
    
    
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
    FatHiggsCandidateFJ3 FatHFJ3;
  VectorCandidate V;
  std::vector<VHbbEvent::SimpleJet> additionalJets;
  std::vector<VHbbEvent::SimpleJet> additionalJetsFat;
};


#endif
