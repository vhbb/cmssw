#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbProxy.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/CutsAndHistos.h"
#include <TH1F.h>
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include <sstream>
#include "TKey.h"


struct CompareJetPt {
  bool operator()( const VHbbEvent::SimpleJet& j1, const  VHbbEvent::SimpleJet& j2 ) const {
    return j1.fourMomentum.Pt() > j2.fourMomentum.Pt();
  }
};

struct CompareBTag {
  bool operator()(const  VHbbEvent::SimpleJet& j1, const  VHbbEvent::SimpleJet& j2 ) const {
    return j1.csv > j2.csv;
  }
};



class SignalRegion: public Cut {
  std::string name() {return "SignalRegion";};
  
  Bool_t pass(VHbbProxy &iProxy){
    
    const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();

    if(iCand->size() < 1)
      return false;

    VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
    VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;

    if(iCand->at(0).candidateType == VHbbCandidate::Zmumu || iCand->at(0).candidateType == VHbbCandidate::Zee){ 
      btag_csv_min = 0.5;
      btag_csv_max = 0.9;
      Higgs_pt = 150;
      V_pt = 150;
      VH_deltaPhi = 2.70;
      nOfAdditionalJets = 2;
      pullAngle = 1.57;
      helicityAngle = 0.8;
    }
    else if(iCand->at(0).candidateType == VHbbCandidate::Wmun || iCand->at(0).candidateType == VHbbCandidate::Wen){
      btag_csv_min = 0.5;
      btag_csv_max = 0.9;
      Higgs_pt = 150;
      V_pt = 150;
      VH_deltaPhi = 2.95;
      nOfAdditionalJets = 2;
      pullAngle = 1.57;
      helicityAngle = 0.8;
    }
    else if(iCand->at(0).candidateType == VHbbCandidate::Znn){
      btag_csv_min = 0.5;
      btag_csv_max = 0.9;
      Higgs_pt = 150;
      V_pt = 150;
      VH_deltaPhi = 2.95;
      nOfAdditionalJets = 2;
      pullAngle = 1.57;
      helicityAngle = 0.8;
    }
    else
      std::cerr << "No vector boson reconstructed. No histos will be filled." << std::endl;
    
    Bool_t go = false;
    if( H.fourMomentum.Pt() > Higgs_pt 
	&& V.fourMomentum.Pt() > V_pt 
	&& TMath::Abs( Geom::deltaPhi(H.fourMomentum.Phi(), V.fourMomentum.Phi()) ) > VH_deltaPhi  
	&& ( H.jets.at(0).csv > btag_csv_min && H.jets.at(1).csv > btag_csv_min )
	&& ( H.jets.at(0).csv > btag_csv_max || H.jets.at(1).csv > btag_csv_max )
	&& iCand->at(0).additionalJets.size() < nOfAdditionalJets 
	&& TMath::Abs(H.deltaTheta) < pullAngle 
	//	&& TMath::Abs(H.helicityAngle) < helicityAngle 
	)
      go = true;
     return go;

  }

private:

  Double_t btag_csv_min;
  Double_t btag_csv_max;
  Double_t Higgs_pt;
  Double_t V_pt;
  Double_t VH_deltaPhi;
  unsigned int nOfAdditionalJets;
  unsigned int nOfAdditionalLeptons;
  Double_t pullAngle;
  Double_t helicityAngle;

};


class ControlRegion_Vusdg: public Cut {
  std::string name() {return "ControlRegion_Vusdg";};
  
  Bool_t pass(VHbbProxy &iProxy){
    
    const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();

    if(iCand->size() < 1)
      return false;

    VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
    VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;

    CSVL = 0.244;
    CSVM = 0.679;
    CSVT = 0.898;

    if(iCand->at(0).candidateType == VHbbCandidate::Zmumu){
      pt_b1 = 20;
      pt_b2 = 20;
      btag_csv_min = CSVL;
      btag_csv_max = CSVL;
      Higgs_pt = 150;
      V_pt = 150;
      VH_deltaPhi = 2.70;
      nOfAdditionalJets = 2;
      pullAngle = 1.57;
      helicityAngle = 0.8;
    }
    else if( iCand->at(0).candidateType == VHbbCandidate::Zee){ 
      pt_b1 = 20;
      pt_b2 = 20;
      btag_csv_min = CSVL;
      btag_csv_max = CSVL;
      Higgs_pt = 150;
      V_pt = 150;
      VH_deltaPhi = 2.70;
      nOfAdditionalJets = 2;
      pullAngle = 1.57;
      helicityAngle = 0.8;
    }
    else if(iCand->at(0).candidateType == VHbbCandidate::Wmun){
      pt_b1 = 30;
      pt_b2 = 30;
      pt_b3 = 30;
      Higgs_pt = 150;
      btag_csv_min = CSVM;
      btag_csv_max = CSVT;
      btag_csv_additional = CSVL;
      V_pt = 150;
      VH_deltaPhi = 2.95;
      nOfAdditionalJets = 2;
      pullAngle = 1.57;
      helicityAngle = 0.8;
    }
    else if(iCand->at(0).candidateType == VHbbCandidate::Wen){
      pt_b1 = 30;
      pt_b2 = 30;
      pt_b3 = 30;
      btag_csv_min = CSVM;
      btag_csv_max = CSVT;
      btag_csv_additional = CSVL;
      Higgs_pt = 150;
      V_pt = 150;
      VH_deltaPhi = 2.95;
      nOfAdditionalJets = 2;
      pullAngle = 1.57;
      helicityAngle = 0.8;
    }
    else
      std::cerr << "No vector boson reconstructed for this control region." << std::endl;

    CompareBTag bTagComparator;
    additionalJetsBtag = iProxy.getVHbbEvent().additionalJets;
    std::sort( additionalJetsBtag.begin(), additionalJetsBtag.end(), bTagComparator );
    
    Bool_t go = false;
    if( H.fourMomentum.Pt() > Higgs_pt 
	&& H.jets.at(0).Pt() > pt_b1 
	&& H.jets.at(1).Pt() > pt_b2 
	&& H.additionalJetsBtag.at(0).Pt() > pt_b3
	&& V.fourMomentum.Pt() > V_pt 
	&& TMath::Abs( Geom::deltaPhi(H.fourMomentum.Phi(), V.fourMomentum.Phi()) ) > VH_deltaPhi  
	&& ( H.jets.at(0).csv > btag_csv_min && H.jets.at(1).csv > btag_csv_min )
	&& ( H.jets.at(0).csv > btag_csv_max || H.jets.at(1).csv > btag_csv_max )
	&& iCand->at(0).additionalJets.size() < nOfAdditionalJets 
	&& TMath::Abs(H.deltaTheta) < pullAngle 
	//	&& TMath::Abs(H.helicityAngle) < helicityAngle 
	)
      go = true;
     return go;

  }

private:

  Double_t CSVL;
  Double_t CSVM;  
  Double_t CSVT;
  Double_t btag_csv_additional;
  Double_t btag_csv_min;
  Double_t btag_csv_max;
  Double_t Higgs_pt;
  Double_t V_pt;
  Double_t VH_deltaPhi;
  unsigned int nOfAdditionalJets;
  unsigned int nOfAdditionalLeptons;
  Double_t pullAngle;
  Double_t helicityAngle;
  std::vector<VHbbEvent::SimpleJet> additionalJetsBtag;

};


class ControlRegion_Top: public Cut {
  std::string name() {return "ControlRegion_Top";};
  
  Bool_t pass(VHbbProxy &iProxy){
    
    const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();

    if(iCand->size() < 1)
      return false;

    VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
    VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;

    if(iCand->at(0).candidateType == VHbbCandidate::Zmumu){
      pt_b1 = 20;
      pt_b2 = 20;
      btag_csv_min = 0.5;
      btag_csv_max = 0.9;
      Higgs_pt = 150;
      V_pt = 150;
      VH_deltaPhi = 2.70;
      nOfAdditionalJets = 2;
      pullAngle = 1.57;
      helicityAngle = 0.8;
    }
    else if( iCand->at(0).candidateType == VHbbCandidate::Zee){ 
      pt_b1 = 20;
      pt_b2 = 20;
      btag_csv_min = 0.5;
      btag_csv_max = 0.9;
      Higgs_pt = 150;
      V_pt = 150;
      VH_deltaPhi = 2.70;
      nOfAdditionalJets = 2;
      pullAngle = 1.57;
      helicityAngle = 0.8;
    }
    else if(iCand->at(0).candidateType == VHbbCandidate::Wmun){
      pt_b1 = 30;
      pt_b2 = 30;
      pt_b3 = 30;
      btag_csv_min = 0.5;
      btag_csv_max = 0.9;
      Higgs_pt = 150;
      V_pt = 150;
      VH_deltaPhi = 2.95;
      nOfAdditionalJets = 2;
      pullAngle = 1.57;
      helicityAngle = 0.8;
    }
    else if(iCand->at(0).candidateType == VHbbCandidate::Wen){
      pt_b1 = 30;
      pt_b2 = 30;
      pt_b3 = 30;
      btag_csv_min = 0.5;
      btag_csv_max = 0.9;
      Higgs_pt = 150;
      V_pt = 150;
      VH_deltaPhi = 2.95;
      nOfAdditionalJets = 2;
      pullAngle = 1.57;
      helicityAngle = 0.8;
    }
    else
      std::cerr << "No vector boson reconstructed for this control region." << std::endl;
    
    Bool_t go = false;
    if( H.fourMomentum.Pt() > Higgs_pt 
	&& H.jets.at(0).Pt() > pt_b1 
	&& H.jets.at(1).Pt() > pt_b2 
	//	&& H.additionalJets.at(0).Pt() > pt_b3
	&& V.fourMomentum.Pt() > V_pt 
	&& TMath::Abs( Geom::deltaPhi(H.fourMomentum.Phi(), V.fourMomentum.Phi()) ) > VH_deltaPhi  
	&& ( H.jets.at(0).csv > btag_csv_min && H.jets.at(1).csv > btag_csv_min )
	&& ( H.jets.at(0).csv > btag_csv_max || H.jets.at(1).csv > btag_csv_max )
	&& iCand->at(0).additionalJets.size() < nOfAdditionalJets 
	&& TMath::Abs(H.deltaTheta) < pullAngle 
	//	&& TMath::Abs(H.helicityAngle) < helicityAngle 
	)
      go = true;
     return go;

  }

private:

  Double_t btag_csv_min;
  Double_t btag_csv_max;
  Double_t Higgs_pt;
  Double_t V_pt;
  Double_t VH_deltaPhi;
  unsigned int nOfAdditionalJets;
  unsigned int nOfAdditionalLeptons;
  Double_t pullAngle;
  Double_t helicityAngle;

};



class ControlRegion_Vbb: public Cut {
  std::string name() {return "ControlRegion_Vbb";};
  
  Bool_t pass(VHbbProxy &iProxy){
    
    const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();

    if(iCand->size() < 1)
      return false;

    VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
    VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;

    if(iCand->at(0).candidateType == VHbbCandidate::Zmumu){
      pt_b1 = 20;
      pt_b2 = 20;
      btag_csv_min = 0.5;
      btag_csv_max = 0.9;
      Higgs_pt = 150;
      V_pt = 150;
      VH_deltaPhi = 2.70;
      nOfAdditionalJets = 2;
      pullAngle = 1.57;
      helicityAngle = 0.8;
    }
    else if( iCand->at(0).candidateType == VHbbCandidate::Zee){ 
      pt_b1 = 20;
      pt_b2 = 20;
      btag_csv_min = 0.5;
      btag_csv_max = 0.9;
      Higgs_pt = 150;
      V_pt = 150;
      VH_deltaPhi = 2.70;
      nOfAdditionalJets = 2;
      pullAngle = 1.57;
      helicityAngle = 0.8;
    }
    else if(iCand->at(0).candidateType == VHbbCandidate::Wmun){
      pt_b1 = 30;
      pt_b2 = 30;
      pt_b3 = 30;
      btag_csv_min = 0.5;
      btag_csv_max = 0.9;
      Higgs_pt = 150;
      V_pt = 150;
      VH_deltaPhi = 2.95;
      nOfAdditionalJets = 2;
      pullAngle = 1.57;
      helicityAngle = 0.8;
    }
    else if(iCand->at(0).candidateType == VHbbCandidate::Wen){
      pt_b1 = 30;
      pt_b2 = 30;
      pt_b3 = 30;
      btag_csv_min = 0.5;
      btag_csv_max = 0.9;
      Higgs_pt = 150;
      V_pt = 150;
      VH_deltaPhi = 2.95;
      nOfAdditionalJets = 2;
      pullAngle = 1.57;
      helicityAngle = 0.8;
    }
    else
      std::cerr << "No vector boson reconstructed for this control region." << std::endl;
    
    Bool_t go = false;
    if( H.fourMomentum.Pt() > Higgs_pt 
	&& H.jets.at(0).Pt() > pt_b1 
	&& H.jets.at(1).Pt() > pt_b2 
	//	&& H.additionalJets.at(0).Pt() > pt_b3
	&& V.fourMomentum.Pt() > V_pt 
	&& TMath::Abs( Geom::deltaPhi(H.fourMomentum.Phi(), V.fourMomentum.Phi()) ) > VH_deltaPhi  
	&& ( H.jets.at(0).csv > btag_csv_min && H.jets.at(1).csv > btag_csv_min )
	&& ( H.jets.at(0).csv > btag_csv_max || H.jets.at(1).csv > btag_csv_max )
	&& iCand->at(0).additionalJets.size() < nOfAdditionalJets 
	&& TMath::Abs(H.deltaTheta) < pullAngle 
	//	&& TMath::Abs(H.helicityAngle) < helicityAngle 
	)
      go = true;
     return go;

  }

private:

  Double_t btag_csv_min;
  Double_t btag_csv_max;
  Double_t Higgs_pt;
  Double_t V_pt;
  Double_t VH_deltaPhi;
  unsigned int nOfAdditionalJets;
  unsigned int nOfAdditionalLeptons;
  Double_t pullAngle;
  Double_t helicityAngle;

};


