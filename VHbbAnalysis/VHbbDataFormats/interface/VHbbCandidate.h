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

  enum CandidateType{Zmumu, Zee, Wmun, Wen, Znn,  Zemu, Ztaumu, Ztaue, Wtaun, Ztautau, Zbb, UNKNOWN};

    VHbbCandidate(){candidateType=UNKNOWN;}

  class VectorCandidate {
  public:
  VectorCandidate() : firstLepton(0),secondLepton(1),firstLeptonOrig(99),secondLeptonOrig(99) {}
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
    
    double MtTau(CandidateType candidateTypeWithTau) const {
      if(candidateTypeWithTau == Wtaun)
	{
	  float ptl=taus[0].p4.Pt();
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
    
    unsigned int firstLepton,secondLepton; // position in the above lepton collections
    unsigned int firstLeptonOrig,secondLeptonOrig; //position in the original VHEvent collections
    
  };
  
  class HiggsCandidate {
  public:
   TLorentzVector p4;
    std::vector<VHbbEvent::SimpleJet> jets;
    size_t indices[2];
    bool HiggsFlag;
    float deltaTheta;
    std::vector <float> helicities;
 public:
    VHbbEvent::SimpleJet& firstJet(){return jets[0];}
    VHbbEvent::SimpleJet& secondJet(){return jets[1];}

    size_t firstJetIndex(){return indices[0];}
    size_t secondJetIndex(){return indices[1];}
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

  void setCandidateType (CandidateType c){candidateType = c;}


  double deltaPhi() const {
   return V.p4.DeltaPhi(H.p4);
  }

  double Mt() const {
   return V.Mt(candidateType);
  }

  double MtTau() const {
   return V.MtTau(candidateTypeWithTau);
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
  CandidateType candidateTypeWithTau;
  HiggsCandidate H;
  FatHiggsCandidate FatH;
  VectorCandidate V;
  VectorCandidate VTau;
  std::vector<VHbbEvent::SimpleJet> additionalJets;
  std::vector<VHbbEvent::SimpleJet> additionalJetsFat;
};


#endif
