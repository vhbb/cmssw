#ifndef VHbbEvent__H 
#define VHbbEvent__H 

#include <TLorentzVector.h>
#include <TVector2.h>
#include <vector>

class VHbbEvent{
 public:

  class ParticleMCInfo {
  public:
    ParticleMCInfo(): status(-99), momid(-99), gmomid(-99), charge(-1){}
  public:
    int status;
    int momid;
    int gmomid;
    int charge;
    TLorentzVector fourMomentum;
    //    int ndau;
    std::vector<int> dauid;
    std::vector<TLorentzVector> dauFourMomentum;    
  };

  class SimpleJet {
    public:
    SimpleJet(): flavour(-99), tche(-99), tchp(-99), jpb(-99), jp(-99), 
      ssvhe(-99), csv(-99), csvmva(-99), ntracks(-99), charge(-99),
      bestMCid(-99), bestMCmomid(-99){}
  public:
    int flavour;
    float tche,tchp, jpb,jp , ssvhe, csv, csvmva;
    int ntracks;
    float charge;
    TLorentzVector fourMomentum;
    TLorentzVector chargedTracksFourMomentum;

    int bestMCid, bestMCmomid;
    // new
    TVector2 tVector;
  };


  class HardJet {
  public:
    HardJet(): constituents(-99){}
  public:
    int constituents;
    TLorentzVector fourMomentum;
    std::vector<TLorentzVector> subFourMomentum;
    std::vector<float> etaSub, phiSub;
  };


  class METInfo {
  public:
    METInfo(): sumEt(-99), metSig(-99), eLong(-99){}
  public:
    float sumEt, metSig, eLong;
    TLorentzVector fourMomentum;
  };

  class MuonInfo {
  public:
    MuonInfo(): charge(-99),tIso(-99), eIso(-99), hIso(-99), 
      acop(-99), ipDb(-99), ipErrDb(-99), zPVPt(-99),zPVProb(-99), chi2(-99), globChi2(-99),
      cat(-99), nHits(-99), nPixelHits(-99), globNHits(-99),validMuStations(-99),
      mcId(-99), mcMomId(-99), mcgMomId(-99){}
  public:
    TLorentzVector fourMomentum;
    int charge;
    float tIso, eIso, hIso, acop, ipDb, ipErrDb, zPVPt,zPVProb, chi2, globChi2;
    int cat, nHits, nPixelHits, globNHits, validMuStations;
    TLorentzVector mcFourMomentum;
    int mcId, mcMomId, mcgMomId;
  };

  class ElectronInfo {
  public:
    ElectronInfo() : scEta(-99), scPhi(-99), charge(-99), 
    tIso(-99), eIso(-99), hIso(-99), 
      acop(-99), id95(-99),id85(-99),id70(-99),id95r(-99), 
      id70r(-99), id85r(-99),mcId(-99), mcMomId(-99), mcgMomId (-99){}
  public:
    TLorentzVector fourMomentum;
    float scEta, scPhi;
    int charge;
    float tIso, eIso, hIso, acop;
    float  id95,id85,id70,id95r, id70r, id85r;
    TLorentzVector mcFourMomentum;
    int mcId, mcMomId, mcgMomId;
  };

  class TauInfo{
  public:
    TauInfo()  : charge(-99), tIso(-99), eIso(-99), hIso(-99), acop(-99), 
    idbyIso(-99),idbyTrackIso(-99),idbyTaNCfrOnePercent(-99),
    idbyTaNCfrHalfPercent(-99), idbyTaNCfrQuarterPercent(-99), 
      idbyTaNCfrTenthPercent(-99), idbyTaNC(-99), mcId(-99), mcMomId(-99), mcgMomId(-99) {}
  public:
    TLorentzVector fourMomentum;
    int charge;
    float tIso, eIso, hIso, acop;
    float  idbyIso,idbyTrackIso,idbyTaNCfrOnePercent,idbyTaNCfrHalfPercent, idbyTaNCfrQuarterPercent, idbyTaNCfrTenthPercent, idbyTaNC;
    TLorentzVector mcFourMomentum;
    int mcId, mcMomId, mcgMomId;
  };


  class DiMuonInfo  {
  public:
    TLorentzVector fourMomentum;
    MuonInfo daughter1, daughter2;
  };


  class DiElectronInfo {
  public:
    TLorentzVector fourMomentum;
    ElectronInfo daughter1, daughter2;
  };


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
  ParticleMCInfo mcH;
  ParticleMCInfo mcW;
  ParticleMCInfo mcZ;
  std::vector<SimpleJet> simpleJets;
  std::vector<SimpleJet> simpleJets2; //???
  std::vector<SimpleJet> subJets; //???
  std::vector<HardJet> hardJets;
  
  ParticleMCInfo mcBbar;
  ParticleMCInfo mcB;
  ParticleMCInfo mcC;    
  
  PUInfo puInfo;
  
  
  
  METInfo calomet;
  METInfo tcmet;
  METInfo pfmet;
  
  std::vector<MuonInfo> muInfo;
  std::vector<ElectronInfo> eleInfo;
  std::vector<TauInfo> tauInfo;
  
  TriggerInfo triggerInfo;
  
  PrimaryVertexInfo pvInfo;
  
  std::vector<DiMuonInfo> diMuonInfo;
  std::vector<DiElectronInfo> diElectronInfo;
  
  //    HiggsCandidate dijetHiggs;
  //    HiggsCandidate fatjetHiggs;
  
  // to be decided.. one single or many V? 
  //  VCandidate vcand;
  //  WCandidate wele;
  //  WCandidate wmu;
  //WCandidate wtau;
  //  ZCandidate zmu;
  //  ZCandidate zele;
  //ZInvCandidate znunu;
  
  //
  // really needed????
  //

};
#endif
