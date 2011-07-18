#ifndef VHbbEventAuxInfo__H 
#define VHbbEventAuxInfo__H 

#include <TLorentzVector.h>
#include <TVector2.h>
#include <vector>

class VHbbEventAuxInfo{
 public:

  class TriggerInfo {
  public:
    TriggerInfo() :  triggerMu9(-99),
      triggerIsoMu9(-99),
      triggerIsoMu13_3(-99),
      triggerMu11(-99),
      triggerDoubleMu3(-99),
      triggerDoubleMu3_2(-99),
      triggerMu15(-99),
      triggerMu15_1(-99),
      triggerDoubleElec10(-99),
      triggerDoubleElec15_1(-99),
      triggerDoubleElec17_1(-99),
      triggerMet100_1(-99),
      triggerSingleEle1(-99),
      triggerSingleEle2(-99),
      triggerSingleEle3(-99),
      triggerSingleEle4(-99),
      triggerBtagMu1(-99),
      triggerBtagMu2(-99),
      triggerBtagMu0(-99),
      triggerBtagMu11(-99),
      triggerBtagMuJet1(-99),
      triggerBtagMuJet2(-99),
      triggerBtagMuJet3(-99),
      triggerBtagMuJet4(-99),
      triggerIsoMu15(-99),
      triggerIsoMu17v5(-99),
      triggerIsoMu17v6(-99) {
      for (unsigned int i=0; i< 500; ++i){
	flag[i]= -99;
      }
    }
  public:
    int flag[500];
    int triggerMu9,
      triggerIsoMu9,
      triggerIsoMu13_3,
      triggerMu11,
      triggerDoubleMu3,
      triggerDoubleMu3_2,
      triggerMu15,
      triggerMu15_1,
      triggerDoubleElec10,
      triggerDoubleElec15_1,
      triggerDoubleElec17_1,
      triggerMet100_1,
      triggerSingleEle1,
      triggerSingleEle2,
      triggerSingleEle3,
      triggerSingleEle4,
      triggerBtagMu1,
      triggerBtagMu2,
      triggerBtagMu0,
      triggerBtagMu11,
      triggerBtagMuJet1,
      triggerBtagMuJet2,
      triggerBtagMuJet3,
      triggerBtagMuJet4,
      triggerIsoMu15,
      triggerIsoMu17v5,
      triggerIsoMu17v6;
  };

  class PrimaryVertexInfo {
  public:
    TVector3 firstPVInPT2,firstPVInProb;
  };

  class PUInfo{
  public:
    PUInfo(): rho(-99) {}
  public:
    float rho;
  };


 public:
  PUInfo puInfo;
  
  TriggerInfo triggerInfo;
  
  PrimaryVertexInfo pvInfo;
  

};
#endif
