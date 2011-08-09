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
  Bool_t pass(VHbbProxy &iProxy){
  const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;
  
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
        && iProxy.trigger()->accept("HLT_IsoMu17_v.*") 
        );
  }
};


class VlightRegionHWen: public Cut {
  std::string name() {return "VlightRegionHWen";};
  Bool_t pass(VHbbProxy &iProxy){
  const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;
  
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
        && iProxy.trigger()->accept("HLT_Ele\\(\\(27\\)\\|\\(32\\)\\)_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v.")
        );
  }
};

class VlightRegionHZmumu: public Cut {
  std::string name() {return "VlightRegionHZmumu";};
  Bool_t pass(VHbbProxy &iProxy){
  const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;

  return (  iCand->at(0).candidateType == VHbbCandidate::Zmumu
        && H.jets.size() >= 2 
        && H.jets.at(0).Pt() > 20
        && H.jets.at(1).Pt() > 20
//        && H.p4.Pt() > 100
        && V.p4.Pt() > 100
        && ( H.jets.at(0).csv < CSVL)
        && ( H.jets.at(1).csv < CSVL)
        && iCand->at(0).additionalJets.size() < 2
        && iProxy.trigger()->accept("HLT_IsoMu17_v.*")
        && V.p4.M() > 75
	&& V.p4.M() < 105
        );
  }
};

class VlightRegionHZee: public Cut {
  std::string name() {return "VlightRegionHZee";};
  Bool_t pass(VHbbProxy &iProxy){
  const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;

  return (  iCand->at(0).candidateType == VHbbCandidate::Zee
        && H.jets.size() >= 2 
        && H.jets.at(0).Pt() > 20
        && H.jets.at(1).Pt() > 20
  //      && H.p4.Pt() > 100
        && V.p4.Pt() > 100
        && ( H.jets.at(0).csv < CSVL)
        && ( H.jets.at(1).csv < CSVL)
        && iCand->at(0).additionalJets.size() < 2
        && iProxy.trigger()->accept("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v.*")       
        && V.p4.M() > 75
	&& V.p4.M() < 105
        );
  }
};



class TTbarRegionHWmun: public Cut {
  std::string name() {return "TTbarRegionHWmun";};
  Bool_t pass(VHbbProxy &iProxy){
  const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;
  
  return (  iCand->at(0).candidateType == VHbbCandidate::Wmun
        && V.muons[0].p4.Pt() > 20
        && H.jets.size() >= 2
        && H.jets.at(0).Pt() > 30
        && H.jets.at(1).Pt() > 30
        && H.p4.Pt() > 100
        && V.p4.Pt() > 100
        && V.Mt(VHbbCandidate::Wmun) > 40
        && V.Mt(VHbbCandidate::Wmun) < 120
        && ( H.jets.at(0).csv > CSVT ||  H.jets.at(1).csv > CSVT)
        && iCand->at(0).additionalJets.size() > 1
        && iProxy.trigger()->accept("HLT_IsoMu17_v.*")
        );
  }
};

class TTbarRegionHWen: public Cut {
  std::string name() {return "TTbarRegionHWen";};
  Bool_t pass(VHbbProxy &iProxy){
  const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;
  
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
        && iProxy.trigger()->accept("HLT_Ele\\(\\(27\\)\\|\\(32\\)\\)_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v.")


        );
  }
};


class TTbarRegionHZmumu: public Cut {
  std::string name() {return "TTbarRegionHZmumu";};
  Bool_t pass(VHbbProxy &iProxy){
  const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;

  return (  iCand->at(0).candidateType == VHbbCandidate::Zmumu
        && H.jets.size() >= 2
        && H.jets.at(0).Pt() > 20
        && H.jets.at(1).Pt() > 20
//        && ( H.jets.at(0).csv > CSVT ||  H.jets.at(1).csv > CSVT)
        && ( V.mets.size() >0 && V.mets.at(0).p4.Pt() > 50)
  //      && iCand->at(0).additionalJets.size() > 1
        && V.p4.M() > 120
        && iProxy.trigger()->accept("HLT_IsoMu17_v.*")
	);
  }
};


class TTbarRegionHZee: public Cut {
  std::string name() {return "TTbarRegionHZee";};
  Bool_t pass(VHbbProxy &iProxy){
  const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;

  return (  iCand->at(0).candidateType == VHbbCandidate::Zee
        && H.jets.size() >= 2
        && H.jets.at(0).Pt() > 20
        && H.jets.at(1).Pt() > 20
      //  && ( H.jets.at(0).csv > CSVT ||  H.jets.at(1).csv > CSVT)
        && ( V.mets.size() >0 && V.mets.at(0).p4.Pt() > 50)
    //    && iCand->at(0).additionalJets.size() > 1
        && iProxy.trigger()->accept("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v.*")

	);
  }
};


class VbbRegionHWmun: public Cut {
  std::string name() {return "VbbRegionHWmun";};
  Bool_t pass(VHbbProxy &iProxy){
  const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;

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
        && V.mets[0].metSig > 1
        && iProxy.trigger()->accept("HLT_IsoMu17_v.*")
        );
  }
};

class VbbRegionHWen: public Cut {
  std::string name() {return "VbbRegionHWen";};
  Bool_t pass(VHbbProxy &iProxy){
  const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;

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
        && iProxy.trigger()->accept("HLT_Ele\\(\\(27\\)\\|\\(32\\)\\)_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v.*")

        );
  }
};

class VbbRegionHZmumu: public Cut {
  std::string name() {return "VbbRegionHZmumu";};
  Bool_t pass(VHbbProxy &iProxy){
  const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;

  return (  iCand->at(0).candidateType == VHbbCandidate::Zmumu
        && H.jets.size() >= 2
        && H.jets.at(0).Pt() > 20
        && H.jets.at(1).Pt() > 20
        &&  ( H.p4.M() < 100 ||  H.p4.M() > 140)
        && ( H.jets.at(0).csv > CSVT ||  H.jets.at(1).csv > CSVT)
        && ( H.jets.at(0).csv > 0.5 && H.jets.at(1).csv > 0.5)
        && TMath::Abs( Geom::deltaPhi(H.p4.Phi(), V.p4.Phi()) ) > 2.9
        && iCand->at(0).additionalJets.size() < 2
        && ( V.mets.size() ==0 || V.mets.at(0).p4.Pt() < 30)
        && iProxy.trigger()->accept("HLT_IsoMu17_v.*")
        && V.p4.M() > 75 
	&& V.p4.M() < 105
  
      );
  }
};

class VbbRegionHZee: public Cut {
  std::string name() {return "VbbRegionHZee";};
  Bool_t pass(VHbbProxy &iProxy){
  const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();
  if(iCand->size() < 1) return false;
  VHbbCandidate::VectorCandidate V = iProxy.getVHbbCandidate()->at(0).V;
  VHbbCandidate::HiggsCandidate H = iProxy.getVHbbCandidate()->at(0).H;

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
        && iProxy.trigger()->accept("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v.*")
	&& V.p4.M() > 75
	&& V.p4.M() < 105
        );
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
    if( H.jets.size() >= 2 && H.p4.Pt() > Higgs_pt 
	&& V.p4.Pt() > V_pt 
	&& TMath::Abs( Geom::deltaPhi(H.p4.Phi(), V.p4.Phi()) ) > VH_deltaPhi  
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

  Double_t pt_b1;
  Double_t pt_b2;
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




