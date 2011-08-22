#ifndef CUTS200X_H
#define CUTS200X_H
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbProxy.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/CutsAndHistos.h"
#include <TH1F.h>
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include <sstream>
#include "TKey.h"


#define CSVM 0.679
#define CSVL 0.244
#define CSVT 0.898

struct CompareJetPt {
  bool operator()( const VHbbEvent::SimpleJet& j1, const  VHbbEvent::SimpleJet& j2 ) const {
    return j1.p4.Pt() > j2.p4.Pt();
  }
};

struct CompareBTag {
  bool operator()(const  VHbbEvent::SimpleJet& j1, const  VHbbEvent::SimpleJet& j2 ) const {
    return j1.csv > j2.csv;
  }
};


// New implementations of the control region
// The signal regions must be implemented incrementally since cutflow is needed

class VlightRegionHWmun: public Cut {
  std::string name() {return "VlightRegionHWmun";};
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = p.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = p.getVHbbCandidate()->at(0).H;
  
  return (  iCand->at(0).candidateType == VHbbCandidate::Wmun
        && V.muons[0].p4.Pt() > 20
        && H.jets.size() >= 2 
        && H.jets.at(0).Pt() > 30
        && H.jets.at(1).Pt() > 30
        && H.p4.Pt() > 150
        && V.p4.Pt() > 150
        && V.Mt(VHbbCandidate::Wmun) < 160
        && ( H.jets.at(0).csv < CSVM)
        && ( H.jets.at(1).csv < CSVM)
        && iCand->at(0).additionalJets.size() < 2
        && V.mets[0].metSig > 2.5
//	&& TMath::Abs( Geom::deltaPhi(H.p4.Phi(), V.p4.Phi()) ) > 2.5 
        && p.trigger()->accept("HLT_IsoMu17_v.*") 
        );
  }
};


class VlightRegionHWen: public Cut {
  std::string name() {return "VlightRegionHWen";};
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = p.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = p.getVHbbCandidate()->at(0).H;
  
  return (  iCand->at(0).candidateType == VHbbCandidate::Wen
        && H.jets.size() >= 2 
        && H.jets.at(0).Pt() > 30
        && H.jets.at(1).Pt() > 30
        && H.p4.Pt() > 150
        && V.p4.Pt() > 150
        && ( H.jets.at(0).csv < CSVM)
        && ( H.jets.at(1).csv < CSVM)
        && iCand->at(0).additionalJets.size() < 5
        && V.mets[0].metSig > 2
	&& TMath::Abs( Geom::deltaPhi(H.p4.Phi(), V.p4.Phi()) ) > 2.5  
        && p.trigger()->accept("HLT_Ele\\(\\(27\\)\\|\\(32\\)\\)_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v.")
        );
  }
};

class VlightRegionHZmumu: public Cut {
  std::string name() {return "VlightRegionHZmumu";};
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = p.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = p.getVHbbCandidate()->at(0).H;

  return (  iCand->at(0).candidateType == VHbbCandidate::Zmumu
        && V.muons[0].p4.Pt() > 20
        && H.jets.size() >= 2 
        && H.jets.at(0).Pt() > 20
        && H.jets.at(1).Pt() > 20
//        && H.p4.Pt() > 100
        && V.p4.Pt() > 100
        && ( H.jets.at(0).csv < CSVL)
        && ( H.jets.at(1).csv < CSVL)
        && iCand->at(0).additionalJets.size() < 2
        && p.trigger()->accept("HLT_IsoMu17_v.*")
        && V.p4.M() > 75
	&& V.p4.M() < 105
        );
  }
};

class VlightRegionHZee: public Cut {
  std::string name() {return "VlightRegionHZee";};
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = p.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = p.getVHbbCandidate()->at(0).H;

  return (  iCand->at(0).candidateType == VHbbCandidate::Zee
        && H.jets.size() >= 2 
        && H.jets.at(0).Pt() > 20
        && H.jets.at(1).Pt() > 20
  //      && H.p4.Pt() > 100
        && V.p4.Pt() > 100
        && ( H.jets.at(0).csv < CSVL)
        && ( H.jets.at(1).csv < CSVL)
        && iCand->at(0).additionalJets.size() < 2
        && p.trigger()->accept("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v.*")       
        && V.p4.M() > 75
	&& V.p4.M() < 105
        );
  }
};



class TTbarRegionHWmun: public Cut {
  std::string name() {return "TTbarRegionHWmun";};
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = p.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = p.getVHbbCandidate()->at(0).H;
  
  return (  iCand->at(0).candidateType == VHbbCandidate::Wmun
        && V.muons[0].p4.Pt() > 20
        && H.jets.size() >= 2
        && H.jets.at(0).Pt() > 30
        && H.jets.at(1).Pt() > 30
        && H.p4.Pt() > 100
        && V.p4.Pt() > 100
      //  && V.Mt(VHbbCandidate::Wmun) > 40
      //  && V.Mt(VHbbCandidate::Wmun) < 120
        && ( H.jets.at(0).csv > CSVT ||  H.jets.at(1).csv > CSVT)
        && iCand->at(0).additionalJets.size() > 1
        && p.trigger()->accept("HLT_IsoMu17_v.*")
        );
  }
};

class TTbarRegionHWen: public Cut {
  std::string name() {return "TTbarRegionHWen";};
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = p.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = p.getVHbbCandidate()->at(0).H;
  
  return (  iCand->at(0).candidateType == VHbbCandidate::Wen
        && H.jets.size() >= 2
        && H.jets.at(0).Pt() > 30
        && H.jets.at(1).Pt() > 30
        && H.p4.Pt() > 75
        && V.p4.Pt() > 75
        && ( H.jets.at(0).csv > CSVT ||  H.jets.at(1).csv > CSVT)
        && TMath::Abs( Geom::deltaPhi(H.p4.Phi(), V.p4.Phi()) ) > 2.5
	&& iCand->at(0).additionalJets.size() > 0
        && iCand->at(0).additionalJets.at(0).p4.Pt() > 35
        && V.mets[0].metSig > 2
        &&  ( TMath::Abs( Geom::deltaPhi( V.mets[0].p4.Phi(), H.jets.at(0).p4.Phi())) > 1.5
            || TMath::Abs( Geom::deltaPhi(V.mets[0].p4.Phi(), H.jets.at(1).p4.Phi())) > 1.5  )
        && p.trigger()->accept("HLT_Ele\\(\\(27\\)\\|\\(32\\)\\)_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v.")


        );
  }
};


class TTbarRegionHZmumu: public Cut {
  std::string name() {return "TTbarRegionHZmumu";};
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = p.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = p.getVHbbCandidate()->at(0).H;

  return (  iCand->at(0).candidateType == VHbbCandidate::Zmumu
        && V.muons[0].p4.Pt() > 20
        && H.jets.size() >= 2
        && H.jets.at(0).Pt() > 20
        && H.jets.at(1).Pt() > 20
        && ( H.jets.at(0).csv > CSVT ||  H.jets.at(1).csv > CSVT)
        && ( V.mets.size() >0 && V.mets.at(0).p4.Pt() > 50)
        && iCand->at(0).additionalJets.size() > 1
        && V.p4.M() > 120
        && p.trigger()->accept("HLT_IsoMu17_v.*")
	);
  }
};


class TTbarRegionHZee: public Cut {
  std::string name() {return "TTbarRegionHZee";};
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = p.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = p.getVHbbCandidate()->at(0).H;

  return (  iCand->at(0).candidateType == VHbbCandidate::Zee
        && H.jets.size() >= 2
        && H.jets.at(0).Pt() > 20
        && H.jets.at(1).Pt() > 20
      //  && ( H.jets.at(0).csv > CSVT ||  H.jets.at(1).csv > CSVT)
        && ( V.mets.size() >0 && V.mets.at(0).p4.Pt() > 50)
    //    && iCand->at(0).additionalJets.size() > 1
        && p.trigger()->accept("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v.*")

	);
  }
};


class VbbRegionHWmun: public Cut {
  std::string name() {return "VbbRegionHWmun";};
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = p.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = p.getVHbbCandidate()->at(0).H;

  return (  iCand->at(0).candidateType == VHbbCandidate::Wmun
        && V.muons[0].p4.Pt() > 20
        && H.jets.size() >= 2
        && H.jets.at(0).Pt() > 30
        && H.jets.at(1).Pt() > 30
        && H.p4.Pt() < 150
        && V.p4.Pt() < 150 
        && V.Mt(VHbbCandidate::Wmun) < 120
        && V.Mt(VHbbCandidate::Wmun) > 40
        && ( H.jets.at(0).csv > CSVT ||  H.jets.at(1).csv > CSVT)
        && TMath::Abs( Geom::deltaPhi(H.p4.Phi(), V.p4.Phi()) ) > 2.5
        && iCand->at(0).additionalJets.size() ==0
        && V.mets[0].metSig > 2.5
        && p.trigger()->accept("HLT_IsoMu17_v.*")
        );
  }
};

class VbbRegionHWen: public Cut {
  std::string name() {return "VbbRegionHWen";};
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = p.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = p.getVHbbCandidate()->at(0).H;

  return (  iCand->at(0).candidateType == VHbbCandidate::Wen
        && H.jets.size() >= 2
        && H.jets.at(0).Pt() > 30
        && H.jets.at(1).Pt() > 30
        && H.p4.Pt() < 150
        && V.p4.Pt() < 150 
        && V.p4.M() > 50
        && ( H.jets.at(0).csv > CSVT ||  H.jets.at(1).csv > CSVT)
        && TMath::Abs( Geom::deltaPhi(H.p4.Phi(), V.p4.Phi()) ) > 2.4
        && iCand->at(0).additionalJets.size() < 2
        && V.mets[0].metSig > 2
        && p.trigger()->accept("HLT_Ele\\(\\(27\\)\\|\\(32\\)\\)_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v.*")

        );
  }
};

class VbbRegionHZmumu: public Cut {
  std::string name() {return "VbbRegionHZmumu";};
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = p.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = p.getVHbbCandidate()->at(0).H;

  return (  iCand->at(0).candidateType == VHbbCandidate::Zmumu
        && V.muons[0].p4.Pt() > 20
        && H.jets.size() >= 2
        && H.jets.at(0).Pt() > 20
        && H.jets.at(1).Pt() > 20
        &&  ( H.p4.M() < 100 ||  H.p4.M() > 140)
//	        && ( H.jets.at(0).csv > CSVT ||  H.jets.at(1).csv > CSVT)
        && ( H.jets.at(0).csv > CSVT && H.jets.at(1).csv > CSVT)
        && TMath::Abs( Geom::deltaPhi(H.p4.Phi(), V.p4.Phi()) ) > 2.9
        && iCand->at(0).additionalJets.size() < 2
        && ( V.mets.size() ==0 || V.mets.at(0).p4.Pt() < 30)
        && p.trigger()->accept("HLT_IsoMu17_v.*")
        && V.p4.M() > 75 
	&& V.p4.M() < 105
  
      );
  }
};

class VbbRegionHZee: public Cut {
  std::string name() {return "VbbRegionHZee";};
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = p.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = p.getVHbbCandidate()->at(0).H;

  return (  iCand->at(0).candidateType == VHbbCandidate::Zee
        && H.jets.size() >= 2
        && H.jets.at(0).Pt() > 20
        && H.jets.at(1).Pt() > 20
        && ( H.p4.M() < 95 ||  H.p4.M() > 145)
        && ( H.jets.at(0).csv > CSVT ||  H.jets.at(1).csv > CSVT)
        && ( H.jets.at(0).csv > 0.5 && H.jets.at(1).csv > 0.5)
//        && TMath::Abs( Geom::deltaPhi(H.p4.Phi(), V.p4.Phi()) ) > 2.9
        && iCand->at(0).additionalJets.size() < 2
  //      && V.mets[0].p4.Pt() < 30
        && p.trigger()->accept("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v.*")
	&& V.p4.M() > 75
	&& V.p4.M() < 105
        );
  }
};


class SignalPreSelectionWen : public Cut {
  std::string name() {return "SignalPreSelWen";};

  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  return (  iCand->at(0).candidateType == VHbbCandidate::Wen);
  }
};

class SignalPreSelectionWmun : public Cut {
  std::string name() {return "SignalPreSelWmun";};
  
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  return (  iCand->at(0).candidateType == VHbbCandidate::Wmun);
  }
};

class SignalPreSelectionZee : public Cut {
  std::string name() {return "SignalPreSelZee";};
  
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  return (  iCand->at(0).candidateType == VHbbCandidate::Zee);
  }
};

class SignalPreSelectionZmumu : public Cut {
  std::string name() {return "SignalPreSelZmumu";};
  
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  return (  iCand->at(0).candidateType == VHbbCandidate::Zmumu);
  }
};

class SignalPreSelectionZnn : public Cut {
  std::string name() {return "SignalPreSelZnn";};
  
  Bool_t pass(VHbbProxy &p){
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  return (  iCand->at(0).candidateType == VHbbCandidate::Znn);
  }
};

class HPtCut : public PCut
{
 public:
  HPtCut(double ptMin):PCut(ptMin){}
  bool pass(VHbbProxy &p) {
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  return (  iCand->at(0).H.p4.Pt() > m_cut);
  }
  virtual std::string name()  {return "Higgs_Pt_Gt_"+cutValueString(); } 
};

class VPtCut : public PCut
{
 public:
  VPtCut(double ptMin):PCut(ptMin){}
  bool pass(VHbbProxy &p) {
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  return (  iCand->at(0).V.p4.Pt() > m_cut);
  }
  virtual std::string name()  {return "Vector_Pt_Gt_"+cutValueString(); }
};


class DoubleBTagCut : public PCut
{
public:
  DoubleBTagCut(double csvMin):PCut(csvMin){}
  bool pass(VHbbProxy &p) {
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::HiggsCandidate  H = p.getVHbbCandidate()->at(0).H;
  return (  
        H.jets.size() >= 2
        && ( H.jets.at(0).csv > m_cut && H.jets.at(1).csv > m_cut)
      );
  } 
  virtual std::string name()  {return "DoubleCSV_"+cutValueString(); }
};

class SingleBTagCut : public PCut
{
  SingleBTagCut(double csvMin):PCut(csvMin){}
  bool pass(VHbbProxy &p) {
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  const VHbbCandidate::HiggsCandidate & H = p.getVHbbCandidate()->at(0).H;
  return (
        H.jets.size() >= 2
        && ( H.jets.at(0).csv > m_cut || H.jets.at(1).csv > m_cut)
      );
  }
  virtual std::string name()  {return "SingleCSV_"+cutValueString(); }

};

class VHDeltaPhiCut : public PCut
{
  VHDeltaPhiCut(double deltaPhiMin):PCut(deltaPhiMin){}
  bool pass(VHbbProxy &p) {
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  const VHbbCandidate::HiggsCandidate & H = p.getVHbbCandidate()->at(0).H;
  const VHbbCandidate::VectorCandidate & V = p.getVHbbCandidate()->at(0).V;
 
 return (TMath::Abs( Geom::deltaPhi(H.p4.Phi(), V.p4.Phi()) ) > m_cut
        );
  }
  virtual std::string name()  {return "VHDeltaPhi_Gt_"+cutValueString(); }
};


class AdditionalJetsCut : public PCut
{
 /// Filtering "< n"
  AdditionalJetsCut(int n):PCut(n) {}
  bool pass(VHbbProxy &p) {
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
 return (iCand->at(0).additionalJets.size() < m_cut);
  }
  virtual std::string name()  {return "NumAddJets_Lt_"+cutValueString(); }
};

class AdditionalLeptonsCut : public PCut
{
 /// Filtering "< n"
  AdditionalLeptonsCut(int n):PCut(n) {}
  bool pass(VHbbProxy &p) {
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  const VHbbCandidate::VectorCandidate & V = p.getVHbbCandidate()->at(0).V;
  int expectedLeptons = 0;
  if(  iCand->at(0).candidateType == VHbbCandidate::Wmun ||  iCand->at(0).candidateType == VHbbCandidate::Wen) expectedLeptons =1;
  if(  iCand->at(0).candidateType == VHbbCandidate::Zmumu ||  iCand->at(0).candidateType == VHbbCandidate::Zee) expectedLeptons =2;
  
   return ( V.muons.size() + V.electrons.size() - expectedLeptons  < m_cut);
  }
  virtual std::string name()  {return "NumAddLeptons_Lt_"+cutValueString(); }
};

class METCut : public PCut
{
 public:
  METCut(double metMin):PCut(metMin){}
  bool pass(VHbbProxy &p) {
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  return (  iCand->at(0).V.mets.at(0).p4.Pt() > m_cut);
  }
  virtual std::string name()  {return "MET_Gt_"+cutValueString(); }
};

class METSignificanceCut : public PCut
{
 public:
  METSignificanceCut(double metMin):PCut(metMin){}
  bool pass(VHbbProxy &p) {
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  return (   iCand->at(0).V.mets[0].metSig  > m_cut);
  }
  virtual std::string name()  {return "METSig_Gt_"+cutValueString(); }
};

class JetMETDeltaPhiCut : public PCut
{
 public:
  JetMETDeltaPhiCut(double deltaPhiMin):PCut(deltaPhiMin){}
  bool pass(VHbbProxy &p) {
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  const VHbbCandidate::VectorCandidate V = p.getVHbbCandidate()->at(0).V;
  const VHbbCandidate::HiggsCandidate H = p.getVHbbCandidate()->at(0).H;

  return ( TMath::Abs( Geom::deltaPhi( V.mets[0].p4.Phi(), H.jets.at(0).p4.Phi())) > m_cut
            && TMath::Abs( Geom::deltaPhi(V.mets[0].p4.Phi(), H.jets.at(1).p4.Phi())) > m_cut  );
}
  virtual std::string name()  {return "JetMETDeltaPhi_Gt_"+cutValueString(); }
};
class DiJetMassMinCut : public PCut
{
 public:
  DiJetMassMinCut(double mMin):PCut(mMin){}
  bool pass(VHbbProxy &p) {
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  return (  iCand->at(0).H.p4.M() > m_cut);
  }
  virtual std::string name()  {return "Higgs_M_Gt_"+cutValueString(); }
};


class DiJetMassMaxCut : public PCut
{
 public:
  DiJetMassMaxCut(double mMax):PCut(mMax){}
  bool pass(VHbbProxy &p) {
  const std::vector<VHbbCandidate> *iCand = p.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  return (  iCand->at(0).H.p4.M() <= m_cut);
  }
  virtual std::string name()  {return "Higgs_M_Lt_"+cutValueString(); }
};


#endif
