#ifndef VHbbEvent__H 
#define VHbbEvent__H 

#include <TLorentzVector.h>
#include <TVector2.h>
#include <vector>

class VHbbEvent{
 public:


  class SimpleJet {
    public:
    SimpleJet(): flavour(-99), isSemiLept(-99), isSemiLeptMCtruth(-99), SoftLeptpdgId(-99), SoftLeptIdlooseMu(-99), SoftLeptId95(-99), SoftLeptPt(-99), SoftLeptdR(-99), SoftLeptptRel(-99),  SoftLeptRelCombIso(-99),   tche(-99), tchp(-99), jpb(-99), jp(-99), 
      ssvhe(-99), csv(-99), csvmva(-99), csvivf(-99), cmva(-99),
      vtxMass(-99), vtx3dL(-99), vtx3deL(-99), vtxNTracks(-99), vtxProb(-99),
      ntracks(-99), charge(-99),jecunc(-99),
      chargedHadronEFraction(-99), neutralHadronEFraction(-99), chargedEmEFraction(-99), neutralEmEFraction(-99),nConstituents(-99), ptRaw(-99), ptLeadTrack(-99), jetArea(-99),
      bestMCid(-99), bestMCmomid(-99){}
  public:
    double Pt() {return p4.Pt();}
    int flavour, isSemiLept, isSemiLeptMCtruth, SoftLeptpdgId , SoftLeptIdlooseMu, SoftLeptId95;
    float SoftLeptPt, SoftLeptdR, SoftLeptptRel, SoftLeptRelCombIso; 
    float tche,tchp, jpb,jp , ssvhe, csv, csvmva,csvivf,cmva;
    // CSV vertex info
    float vtxMass, vtx3dL, vtx3deL, vtxNTracks, vtxProb;
    TVector3 vtxPosition;
    TLorentzVector  vtxP4;
    std::vector<unsigned int> vtxTrackIds, btagTrackIds, csvTrackIds;
    int csvNTracks, btagNTracks;
    int ntracks;
    float charge;
    TLorentzVector p4;
    TLorentzVector chargedTracksFourMomentum;

    // jec uncertainty

    float jecunc;
    // energy fractions as in 
    /*    //chargedHadronEnergyFraction (relative to uncorrected jet energy)
	  float chargedHadronEnergyFraction() const {return chargedHadronEnergy()/((jecSetsAvailable() ? jecFactor(0) : 1.)*energy());}
	  /// neutralHadronEnergyFraction (relative to uncorrected jet energy)
	  float neutralHadronEnergyFraction() const {return neutralHadronEnergy()/((jecSetsAvailable() ? jecFactor(0) : 1.)*energy());}
	  /// chargedEmEnergyFraction (relative to uncorrected jet energy)
	  float chargedEmEnergyFraction()     const {return chargedEmEnergy()/((jecSetsAvailable() ? jecFactor(0) : 1.)*energy());}
	  /// neutralEmEnergyFraction (relative to uncorrected jet energy)
	  float neutralEmEnergyFraction()     const {return neutralEmEnergy()/((jecSetsAvailable() ? jecFactor(0) : 1.)*energy());}
    */

    float chargedHadronEFraction, neutralHadronEFraction, chargedEmEFraction, neutralEmEFraction;
    int nConstituents;
    float ptRaw,ptLeadTrack;  
    float jetArea;  

    //MC parton matching
    TLorentzVector bestMCp4,bestMCp4mom;
    int bestMCid, bestMCmomid;
    // new
    TVector2 tVector;
    //Scale factor
    float SF_CSVL;
    float SF_CSVM;
    float SF_CSVT; 
    float SF_CSVLerr;
    float SF_CSVMerr;
    float SF_CSVTerr; 
    
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
    MuonInfo(): charge(-99),tIso(-99), eIso(-99), hIso(-99),pfChaIso(-99),pfChaPUIso(99), pfPhoIso(-99), pfNeuIso(-99),
      acop(-99), ipDb(-99), ipErrDb(-99), zPVPt(-99),zPVProb(-99), chi2(-99), globChi2(-99),
      cat(-99), nValidTracker(-99), nValidPixel(-99), nMatches(-99),nHits(-99), nPixelHits(-99), globNHits(-99),validMuStations(-99),emEnergy(-99), hadEnergy(-99),nValidLayers(-99),isPF(-99),
      mcId(-99), mcMomId(-99), mcgMomId(-99){}
  public:
    TLorentzVector p4;
    int charge;
    float tIso, eIso, hIso, pfChaIso,pfChaPUIso,  pfPhoIso,pfNeuIso,acop, ipDb, ipErrDb, zPVPt,zPVProb, chi2, globChi2;
    int cat, nValidTracker, nValidPixel, nMatches,nHits, nPixelHits, globNHits, validMuStations;
    float emEnergy, hadEnergy;
    TLorentzVector mcFourMomentum;
    int mcId, mcMomId, mcgMomId;
    int isPF,nValidLayers;
    std::vector< int > hltMatchedBits;
  };

  class ElectronInfo {
  public:
    ElectronInfo() : scEta(-99), scPhi(-99), charge(-99), 
      tIso(-99), eIso(-99), hIso(-99),pfChaIso(-99),pfChaPUIso(-99), pfPhoIso(-99), pfNeuIso(-99), acop(-99),
      sihih(-99), Dphi(-99), Deta(-99), HoE(-99), convDist(-99), convDcot(-99), innerHits(-99),
      isEB(false),isEE(false),
      ipDb(-99), ipErrDb(-99),
      id95(-99),id85(-99),id80(-99),id70(-99),
      id95r(-99),id85r(-99),id80r(-99),id70r(-99),
      mcId(-99), mcMomId(-99), mcgMomId (-99),
      dxy(-99), dz(-99),
     //MVA ELEID 2012 input vars (keep same names to make life easier in integrating with MVA evaluation code
      fMVAVar_EoP(-99), fMVAVar_HoE(-99), fMVAVar_IoEmIoP(-99), fMVAVar_PreShowerOverRaw(-99), fMVAVar_R9(-99), fMVAVar_d0(-99), fMVAVar_deta(-99), fMVAVar_detacalo(-99), fMVAVar_dphi(-99), fMVAVar_e1x5e5x5(-99), fMVAVar_eleEoPout(-99), fMVAVar_eta(-99), fMVAVar_etawidth(-99), fMVAVar_fbrem(-99), fMVAVar_gsfchi2(-99), fMVAVar_ip3d(-99), fMVAVar_kfchi2(-99), fMVAVar_kfhits(-99), fMVAVar_phiwidth(-99), fMVAVar_pt(-99), fMVAVar_see(-99), fMVAVar_spp(-99), mvaOut(-99),mvaOutTrig(-99)
{}
  public:
    TLorentzVector p4;
    float scEta, scPhi;
    int charge;
    float tIso, eIso, hIso, pfChaIso,pfChaPUIso, pfPhoIso,pfNeuIso, acop;
    float sihih, Dphi, Deta, HoE, convDist, convDcot;
    int innerHits;
    bool isEB,isEE;
    float       ipDb, ipErrDb, dxy,dz;
    float  id95,id85,id80,id70,id95r, id85r,id80r, id70r;
    float fMVAVar_EoP, fMVAVar_HoE, fMVAVar_IoEmIoP, fMVAVar_PreShowerOverRaw, fMVAVar_R9, fMVAVar_d0, fMVAVar_deta, fMVAVar_detacalo, fMVAVar_dphi, fMVAVar_e1x5e5x5, fMVAVar_eleEoPout, fMVAVar_eta, fMVAVar_etawidth, fMVAVar_fbrem, fMVAVar_gsfchi2, fMVAVar_ip3d, fMVAVar_kfchi2, fMVAVar_kfhits, fMVAVar_phiwidth, fMVAVar_pt, fMVAVar_see, fMVAVar_spp, mvaOut,mvaOutTrig;
    TLorentzVector mcFourMomentum;
    int mcId, mcMomId, mcgMomId;
    std::vector< int > hltMatchedBits;
  };

  class TauInfo{
  public:
    TauInfo()  : charge(-99), tIso(-99), eIso(-99), hIso(-99),pfChaIso(-99), pfChaPUIso(-99),pfPhoIso(-99), pfNeuIso(-99), acop(-99), 
      byIsolation(-99),trackIsolation(-99),byTaNCfrOnePercent(-99),
      byTaNCfrHalfPercent(-99), byTaNCfrQuarterPercent(-99), 
      byTaNCfrTenthPercent(-99), byTaNC(-99), mcId(-99), mcMomId(-99), mcgMomId(-99),
      againstElectronLoose(-99),  againstElectronMedium(-99),  againstElectronTight(-99),  againstMuonLoose(-99),  againstMuonTight(-99),  
      byLooseIsolation(-99),  byMediumIsolation(-99),  byTightIsolation(-99),  byVLooseIsolation(-99),  decayModeFinding(-99),
      isolationPFChargedHadrCandsPtSum(-999999.), isolationPFGammaCandsEtSum(-999999.),leadPFChargedHadrCandPt(-999999.),
      NsignalPFChargedHadrCands (-99),NsignalPFGammaCands(-99), byLooseCombinedIsolationDeltaBetaCorr(-99), againstElectronMVA(-99) {}
  public:
    TLorentzVector p4;
    int charge;
    float tIso, eIso, hIso,pfChaIso,pfChaPUIso,pfPhoIso,pfNeuIso, acop;
    float  byIsolation,trackIsolation,byTaNCfrOnePercent,byTaNCfrHalfPercent, byTaNCfrQuarterPercent, byTaNCfrTenthPercent, byTaNC;
    int mcId, mcMomId, mcgMomId;
    float againstElectronLoose, againstElectronMedium, againstElectronTight, againstMuonLoose, againstMuonTight;
    float byLooseIsolation, byMediumIsolation, byTightIsolation, byVLooseIsolation, decayModeFinding;

    float isolationPFChargedHadrCandsPtSum,isolationPFGammaCandsEtSum,leadPFChargedHadrCandPt;
    int NsignalPFChargedHadrCands, NsignalPFGammaCands;
    float byLooseCombinedIsolationDeltaBetaCorr, againstElectronMVA;


    TLorentzVector mcFourMomentum;
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
  std::vector<SimpleJet> simpleJets3;
  std::vector<SimpleJet> simpleJets4; //???
  std::vector<SimpleJet> subJets; //???
  std::vector<SimpleJet> filterJets; //???
  std::vector<HardJet> hardJets;
  
  METInfo calomet;
  METInfo tcmet;
  METInfo pfmet;
  METInfo pfmetType1corr;
  METInfo pfmetType1p2corr;
  METInfo pfmetNoPUType1corr;
  METInfo pfmetNoPUType1p2corr;
  METInfo mht;
  METInfo metNoPU;
  METInfo metCh;
  std::vector<METInfo> metUncInfo;

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
