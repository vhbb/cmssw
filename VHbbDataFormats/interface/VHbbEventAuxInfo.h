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
    PrimaryVertexInfo() : nVertices (-99){}
    TVector3 firstPVInPT2,firstPVInProb;
    int nVertices;
  };

  class PUInfo{
  public:
    PUInfo(): rho(-99) {}
  public:
    float rho;
  };


  class ParticleMCInfo {
  public:
    ParticleMCInfo(): status(-99), momid(-99), gmomid(-99), charge(-99){}
  public:
    int status;
    int momid;
    int gmomid;
    float charge;
    TLorentzVector p4;
    //    int ndau;
    std::vector<int> dauid;
    std::vector<TLorentzVector> dauFourMomentum;    
  };

 double genBBDeltaR() const
 {
   if(mcB.size() > 0 && mcBbar.size() > 0)
   return mcB[0].p4.DeltaR(mcBbar[0].p4);
   else return -99;
 }
 double genCCDeltaR() const
 {
   if(mcC.size() >=2)
   return mcC[0].p4.DeltaR(mcC[1].p4);
   else return -99;
 }
 
 public:
  PUInfo puInfo;
  
  TriggerInfo triggerInfo;
  
  PrimaryVertexInfo pvInfo;
  std::vector<ParticleMCInfo> mcH;
  std::vector<ParticleMCInfo> mcW;
  std::vector<ParticleMCInfo> mcZ;
  
  std::vector<ParticleMCInfo> mcBbar;
  std::vector<ParticleMCInfo> mcB;
  std::vector<ParticleMCInfo> mcC;       

};
#endif
