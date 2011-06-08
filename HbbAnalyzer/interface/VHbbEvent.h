#ifndef VHbbEvent__H 
#define VHbbEvent__H 

#include <TLorentzVector.h>
#include <TMatrix.h>
#include <vector>

class VHbbEvent{
 public:

  class ParticleMCInfo {
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
    int flavour;
    float tche,tchp, jpb,jp , ssvhe, csv, csvmva;
    int ntracks;
    float charge;
    TLorentzVector fourMomentum;
    int b1BestMCid, b1BestMCmomid;
  };


  class HardJet{
  public:
    int constituents;
    TLorentzVector fourMomentum;
    std::vector<TLorentzVector> subFourMomentum;
    std::vector<float> etaSub, phiSub;
  };


  class METInfo{
  public:
    float sumEt, metSig, eLong;
    TLorentzVector fourMomentum;
  };

  class MuonInfo {
  public:
    TLorentzVector fourMomentum;
    int charge;
    float tIso, eIso, hIso, acop, ipDb, ipErrDb, zPVPt,zPVProb, chi2, globChi2;
    int cat, nHits, nPixelHits, globNHits;
    TLorentzVector mcFourMomentum;
    int mcId, mcMomId, mcgMomId;
  };


  class ElectronInfo {
  public:
    TLorentzVector fourMomentum;
    float scEta, scPhi;
    int charge;
    float tIso, eIso, hIso, acop;
    float  id95,id85,id70,id95r, id70r;
    TLorentzVector mcFourMomentum;
    int mcId, mcMomId, mcgMomId;
  };

  class TauInfo {
  public:
    TLorentzVector fourMomentum;
    int charge;
    float tIso, eIso, hIso, acop;
    float  idbyIso,idbyTrackIso,idbyTaNCfrOnePercent,idbyTaNCfrHalfPercent, idbyTaNCfrQuarterPercent, idbyTaNCfrTenthPercent, idbyTaNC;
    TLorentzVector mcFourMomentum;
    int mcId, mcMomId, mcgMomId;
  };


  class DiMuonInfo {
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

  class PUInfo {
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
