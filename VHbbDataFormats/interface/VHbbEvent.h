#ifndef VHbbEvent__H 
#define VHbbEvent__H 

#include <TLorentzVector.h>
#include <TVector2.h>
#include <vector>

class VHbbEvent{
 public:


  class SimpleJet {
    public:
    SimpleJet(): flavour(-99), tche(-99), tchp(-99), jpb(-99), jp(-99), 
      ssvhe(-99), csv(-99), csvmva(-99), ntracks(-99), charge(-99),
      bestMCid(-99), bestMCmomid(-99){}
  public:
    double Pt() {return p4.Pt();}
    int flavour;
    float tche,tchp, jpb,jp , ssvhe, csv, csvmva;
    int ntracks;
    float charge;
    TLorentzVector p4;
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
    TLorentzVector p4;
    std::vector<TLorentzVector> subFourMomentum;
    std::vector<float> etaSub, phiSub;
  };


  class METInfo {
  public:
    METInfo(): sumEt(-99), metSig(-99), eLong(-99){}
  public:
    float sumEt, metSig, eLong;
    TLorentzVector p4;
  };

  class MuonInfo {
  public:
    MuonInfo(): charge(-99),tIso(-99), eIso(-99), hIso(-99),pfChaIso(-99), pfPhoIso(-99), pfNeuIso(-99),
      acop(-99), ipDb(-99), ipErrDb(-99), zPVPt(-99),zPVProb(-99), chi2(-99), globChi2(-99),
      cat(-99), nHits(-99), nPixelHits(-99), globNHits(-99),validMuStations(-99),
      mcId(-99), mcMomId(-99), mcgMomId(-99){}
  public:
    TLorentzVector p4;
    int charge;
    float tIso, eIso, hIso, pfChaIso,pfPhoIso,pfNeuIso,acop, ipDb, ipErrDb, zPVPt,zPVProb, chi2, globChi2;
    int cat, nHits, nPixelHits, globNHits, validMuStations;
    TLorentzVector mcFourMomentum;
    int mcId, mcMomId, mcgMomId;
  };

  class ElectronInfo {
  public:
    ElectronInfo() : scEta(-99), scPhi(-99), charge(-99), 
    tIso(-99), eIso(-99), hIso(-99),pfChaIso(-99), pfPhoIso(-99), pfNeuIso(-99), acop(-99),
      id95(-99),id85(-99),id80(-99),id70(-99),
      id95r(-99),id85r(-99),id80r(-99),id70r(-99),
      mcId(-99), mcMomId(-99), mcgMomId (-99){}
  public:
    TLorentzVector p4;
    float scEta, scPhi;
    int charge;
    float tIso, eIso, hIso, pfChaIso,pfPhoIso,pfNeuIso, acop;
    float  id95,id85,id80,id70,id95r, id85r,id80r, id70r;
    TLorentzVector mcFourMomentum;
    int mcId, mcMomId, mcgMomId;
  };

  class TauInfo{
  public:
    TauInfo()  : charge(-99), tIso(-99), eIso(-99), hIso(-99),pfChaIso(-99), pfPhoIso(-99), pfNeuIso(-99), acop(-99), 
    idbyIso(-99),idbyTrackIso(-99),idbyTaNCfrOnePercent(-99),
    idbyTaNCfrHalfPercent(-99), idbyTaNCfrQuarterPercent(-99), 
      idbyTaNCfrTenthPercent(-99), idbyTaNC(-99), mcId(-99), mcMomId(-99), mcgMomId(-99) {}
  public:
    TLorentzVector p4;
    int charge;
    float tIso, eIso, hIso,pfChaIso,pfPhoIso,pfNeuIso, acop;
    float  idbyIso,idbyTrackIso,idbyTaNCfrOnePercent,idbyTaNCfrHalfPercent, idbyTaNCfrQuarterPercent, idbyTaNCfrTenthPercent, idbyTaNC;
    TLorentzVector mcFourMomentum;
    int mcId, mcMomId, mcgMomId;
  };


  class DiMuonInfo  {
  public:
    TLorentzVector p4;
    MuonInfo daughter1, daughter2;
  };


  class DiElectronInfo {
  public:
    TLorentzVector p4;
    ElectronInfo daughter1, daughter2;
  };



 public:
  std::vector<SimpleJet> simpleJets;
  std::vector<SimpleJet> simpleJets2; //???
  std::vector<SimpleJet> subJets; //???
  std::vector<HardJet> hardJets;
  
 METInfo calomet;
  METInfo tcmet;
  METInfo pfmet;
  
  std::vector<MuonInfo> muInfo;
  std::vector<ElectronInfo> eleInfo;
  std::vector<TauInfo> tauInfo;
  
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
