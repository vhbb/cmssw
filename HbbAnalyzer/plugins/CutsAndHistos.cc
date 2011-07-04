#include "VHbbAnalysis/HbbAnalyzer/interface/VHbbEvent.h"
#include "VHbbAnalysis/HbbAnalyzer/interface/VHbbCandidate.h"
#include "VHbbAnalysis/HbbAnalyzer/interface/VHbbProxy.h"
#include "VHbbAnalysis/HbbAnalyzer/interface/CutsAndHistos.h"
#include <TH1F.h>
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include <sstream>
#include "TKey.h"

//class LeptonSelection : public Cut 

class mcHPt70Cut :  public Cut {
  std::string name() {return "mcHpT70";};
  bool pass(VHbbProxy &iProxy) {
    return iProxy.getVHbbEvent()->mcH.fourMomentum.Pt() > 70;
  } 
};

class mcHPt100Cut :  public Cut {
  std::string name() {return "mcHpT100";};
  bool pass(VHbbProxy &iProxy) {
    return iProxy.getVHbbEvent()->mcH.fourMomentum.Pt() > 100;
  } 
};

class mcHPt150Cut :  public Cut {
  std::string name() {return "mcHpT150";};
  bool pass(VHbbProxy &iProxy) {
    return iProxy.getVHbbEvent()->mcH.fourMomentum.Pt() > 150;
  } 
};

class mcWPt50To150Cut:  public Cut {
  std::string name() {return "mcWpT50To150";};
  bool pass(VHbbProxy &iProxy) {
    const VHbbEvent* iEvent = iProxy.getVHbbEvent();
    return ( iEvent->mcW.fourMomentum.Pt() > 50 
	     && iEvent->mcW.fourMomentum.Pt() < 150 );
  } 
};

class mcZPt50To150Cut:  public Cut {
  std::string name() {return "mcZpT50To150";};
  bool pass(VHbbProxy &iProxy) {
    const VHbbEvent* iEvent = iProxy.getVHbbEvent();
    return ( iEvent->mcZ.fourMomentum.Pt() > 50 
	     && iEvent->mcZ.fourMomentum.Pt() < 150 );
  } 
};

class VPt150Cut: public Cut {
  std::ostringstream  V_pt;
  std::string name() {
    V_pt << 150;
    return "VpT"+V_pt.str();
  }
  bool pass(VHbbProxy &iProxy) {
    if(iProxy.getVHbbCandidate()->size() < 1)
      return false;
    else 
      return ((iProxy.getVHbbCandidate()->at(0)).V.fourMomentum.Pt() > 150);
  } 
};

class HPt150Cut: public Cut {
  std::ostringstream H_pt;
  std::string name() {
    H_pt << 150;
    return "HpT"+H_pt.str();
  }
  bool pass(VHbbProxy &iProxy) {
    if(iProxy.getVHbbCandidate()->size() < 1)
      return false;
    else 
      return ((iProxy.getVHbbCandidate()->at(0)).H.fourMomentum.Pt() > 150);
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
	//	&& V.additionalLeptons.size() < nOfAddionalLeptons
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




class MCHistos : public Histos {
  
public:

  virtual void book(TFile &f, std::string suffix) {

    TDirectory *subDir;

    if( ! f.GetDirectory(suffix.c_str()) )
      subDir = f.mkdir(suffix.c_str());
    else 
      subDir = f.GetDirectory(suffix.c_str());

    subDir->cd();

    bin_mass = 500;
    min_mass = 0;
    max_mass = 300;
    
    bin_pt = 500;
    min_pt = 0;
    max_pt = 500;
    
    bin_hel = 50;
    min_hel = 0;
    max_hel = 1;

    //from MC
    h_simHMass = new TH1F(("simHiggsMass"+suffix).c_str(),("Sim Higgs Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    h_simHPt = new TH1F(("simHiggsPt"+suffix).c_str(),("Sim Higgs Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    h_simZMass = new TH1F(("simZMass"+suffix).c_str(),("Sim Z Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    h_simZPt = new TH1F(("simZPt"+suffix).c_str(),("Sim Z Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    h_simWMass = new TH1F(("simWMass"+suffix).c_str(),("Sim W Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    h_simWPt = new TH1F(("simWPt"+suffix).c_str(),("Sim W Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );

  }
  
  virtual void fill(VHbbProxy &iProxy,float w) {

    const VHbbEvent *iEvent = iProxy.getVHbbEvent();

    //from MC    
    h_simHMass->Fill(iEvent->mcH.fourMomentum.M(), w); 
    h_simHPt->Fill(iEvent->mcH.fourMomentum.Pt(), w); 
    h_simZMass->Fill(iEvent->mcZ.fourMomentum.M(), w); 
    h_simZPt->Fill(iEvent->mcZ.fourMomentum.Pt(), w); 
    h_simWMass->Fill(iEvent->mcW.fourMomentum.M(), w); 
    h_simWPt->Fill(iEvent->mcW.fourMomentum.Pt(), w); 

  }


  TH1F * h_simHMass;
  TH1F * h_simHPt;
  TH1F * h_simWMass;
  TH1F * h_simWPt;
  TH1F * h_simZMass;
  TH1F * h_simZPt;
 
private:

  Int_t bin_mass;
  Double_t min_mass;
  Double_t max_mass;

  Int_t bin_pt;
  Double_t min_pt;
  Double_t max_pt;

  Int_t bin_hel;
  Double_t min_hel;
  Double_t max_hel;


};

class StandardHistos : public Histos {
  
public:

  virtual void book(TFile &f, std::string suffix) {


    TDirectory *subDir = f.mkdir(suffix.c_str());
    subDir->cd();

    bin_mass = 500;
    min_mass = 0;
    max_mass = 300;
    
    bin_pt = 500;
    min_pt = 0;
    max_pt = 500;
    
    bin_hel = 50;
    min_hel = 0;
    max_hel = 1;

    //Candidates
    h_HMass = new TH1F(("HiggsMass"+suffix).c_str(),(" Higgs Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    h_HPt = new TH1F(("HiggsPt"+suffix).c_str(),(" Higgs Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    h_ZMass = new TH1F(("ZMass"+suffix).c_str(),(" Z Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    h_ZPt = new TH1F(("ZPt"+suffix).c_str(),(" Z Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    h_WMass = new TH1F(("WMass"+suffix).c_str(),(" W Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    h_WPt = new TH1F(("WPt"+suffix).c_str(),(" W Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    h_HHel = new TH1F(("HiggsHel"+suffix).c_str(),("Higgs helicity angle ("+suffix+")").c_str(), bin_hel, min_hel, max_hel );
  }
  
  virtual void fill(VHbbProxy &iProxy,float w) {

    const VHbbEvent *iEvent = iProxy.getVHbbEvent();
    const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();

    //Candidates
    if(iCand->size() > 0){
      VHbbCandidate::CandidateType iCandType = iCand->at(0).candidateType;
      VHbbCandidate::HiggsCandidate H = iCand->at(0).H;
      VHbbCandidate::VectorCandidate V = iCand->at(0).V;
      
      h_HMass->Fill(H.fourMomentum.M(), w); 
      h_HPt->Fill(H.fourMomentum.Pt(), w); 
      if( iCandType == VHbbCandidate::Zmumu || iCandType == VHbbCandidate::Zee || iCandType == VHbbCandidate::Znn ){
	h_ZMass->Fill(V.fourMomentum.M(), w); 
	h_ZPt->Fill(V.fourMomentum.Pt(), w);
      } 
      else if(iCandType == VHbbCandidate::Wen || iCandType == VHbbCandidate::Wmun){
	h_WMass->Fill(V.fourMomentum.M(), w); 
	h_WPt->Fill(V.fourMomentum.Pt(), w); 
      }
      //    h_HHel->Fill(H.hel(), w); 
    }
  }
  
  TH1F * h_HMass;
  TH1F * h_HPt;
  TH1F * h_HHel;
  TH1F * h_WMass;
  TH1F * h_WPt;
  TH1F * h_ZMass;
  TH1F * h_ZPt;

  TH1F * h_pull;
  TH1F * h_helicity;
  TH1F * h_dPhiVH;
  TH1F * h_dRVH;
 
private:

  Int_t bin_mass;
  Double_t min_mass;
  Double_t max_mass;

  Int_t bin_pt;
  Double_t min_pt;
  Double_t max_pt;

  Int_t bin_hel;
  Double_t min_hel;
  Double_t max_hel;


};
