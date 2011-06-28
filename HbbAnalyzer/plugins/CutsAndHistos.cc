#include "VHbbAnalysis/HbbAnalyzer/interface/VHbbEvent.h"
#include "VHbbAnalysis/HbbAnalyzer/interface/CutsAndHistos.hpp"
#include <TH1F.h>
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

//class LeptonSelection : public Cut 

class HPt70Cut :  public Cut {
  std::string name() {return "HpT70";};
  bool pass(const VHbbEvent &e) {return e.mcH.fourMomentum.Pt() > 70;} 
};

class WPt50To150Cut:  public Cut {
  std::string name() {return "WpT50To150";};
  bool pass(const VHbbEvent &e) {return e.mcW.fourMomentum.Pt() > 50 && e.mcW.fourMomentum.Pt() < 150;} 
};

class ZPt50To150Cut:  public Cut {
  std::string name() {return "ZpT50To150";};
  bool pass(const VHbbEvent &e) {return e.mcZ.fourMomentum.Pt() > 50 && e.mcZ.fourMomentum.Pt() < 150;} 
};

// class VBoost: public Cut {
//   std::stringstream  V_pt = 150;
//   std::tsringstream V_type;
//   std::string name() {return V_type+V_pt.str()}
//   bool pass(VHbbEvent &e) {return e.mcZ.fourMomentum.Pt() > 50 && e.mcZ.fourMomentum.Pt() < 150;} 

// };

// class SignalRegion: public Cut {
//   std::string name() {return "SignalRegion";};

//   VectorCandidate V = e.VHbbCandidate.VectorCandidate;
//   HiggsCandidate H = e.VHbbCandidate.HiggsCandidate;

//   if(V.isZLL){ 
//     btag_min = 0.5;
//     btag_max = 0.9;
//     Higgs_pt = 150;
//     V_pt = 150;
//     VH_deltaPhi = 2.70;
//     nOfAdditionalJets = 2;
//     pullAngle = 1.57;
//     helicityAngle = 0.8;
//   }
//   else if(V.isWLNu){
//     btag_min = 0.5;
//     btag_max = 0.9;
//     Higgs_pt = 150;
//     V_pt = 150;
//     VH_deltaPhi = 2.95;
//     nOfAdditionalJets = 2;
//     pullAngle = 1.57;
//     helicityAngle = 0.8;
//   }
//   else if(V.isZNuNu){
//     btag_min = 0.5;
//     btag_max = 0.9;
//     Higgs_pt = 150;
//     V_pt = 150;
//     VH_deltaPhi = 2.95;
//     nOfAdditionalJets = 2;
//     pullAngle = 1.57;
//     helicityAngle = 0.8;
//   }
//   else
//     std::cerr << "No vector boson reconstructed. 
// No histos will be filled." << std::endl;
  
//   Bool_t pass(VHbbEvent &e){
//     Bool_t go = false;
//     if( H.fourMomentum.Pt() > Higgs_pt 
// 	&& e.VCandidate.fourMomentum.Pt() > V_pt 
// 	&& TMath::Abs( Geom::deltaPhi(H.fourMomentum.Phi(), e.ZCandidate.fourMomentum.Phi()) ) > VH_deltaPhi  
// 	&& ( H.jets.at(0).csv > btag_csv_min && H.jets.at(1).csv > btag_csv_min )
// 	&& ( H.jets.at(0).csv > btag_csv_max || H.jets.at(1).csv > btag_csv_max )
// 	&& H.additionalJets.size() < nOfAddionalJets 
// 	&& e.VCandidate.additionalLeptons.size() < nOfAddionalLeptons
// 	&& TMath::Abs(H.pullAngle) < pullAngle 
// 	&& TMath::Abs(H.helicityAngle) < helicityAngle )
//       go = true;
//      return go;
//   }

// private:

//   Double_t btag_min;
//   Double_t btag_max;
//   Double_t Higgs_pt;
//   Double_t V_pt;
//   Double_t VH_deltaPhi;
//   Int_t nOfAdditionalJets;
//   Int_t nOfAdditionalLeptons;
//   Double_t pullAngle;
//   Double_t helicityAngle;

// };


// class ControlRegion_noBoost_bTag: public Cut {
//   std::string name() {return "ControlRegion_noBoost_bTag";};

//   VectorCandidate V = e.VHbbCandidate.VectorCandidate;
//   HiggsCandidate H = e.VHbbCandidate.HiggsCandidate;

//   if(V.isZLL){ 
//     btag_min = 0.5;
//     btag_max = 0.9;
//     Higgs_pt = 150;
//     V_pt = 150;
//     VH_deltaPhi = 2.70;
//     nOfAdditionalJets = 2;
//     pullAngle = 1.57;
//     helicityAngle = 0.8;
//   }
//   else if(V.isWLNu){
//     btag_min = 0.5;
//     btag_max = 0.9;
//     Higgs_pt = 150;
//     V_pt = 150;
//     VH_deltaPhi = 2.95;
//     nOfAdditionalJets = 2;
//     pullAngle = 1.57;
//     helicityAngle = 0.8;
//   }
//   else if(V.isZNuNu){
//     btag_min = 0.5;
//     btag_max = 0.9;
//     Higgs_pt = 150;
//     V_pt = 150;
//     VH_deltaPhi = 2.95;
//     nOfAdditionalJets = 2;
//     pullAngle = 1.57;
//     helicityAngle = 0.8;
//   }
//   else
//     std::cerr << "No vector boson reconstructed. 
// No histos will be filled." << std::endl;

//   Bool_t pass(VHbbEvent &e){
//     Bool_t go = false;
//     if( H.fourMomentum.Pt() < Higgs_pt 
// 	&& e.VCandidate.fourMomentum.Pt() < V_pt 
// 	&& TMath::Abs( Geom::deltaPhi(H.fourMomentum.Phi(), e.ZCandidate.fourMomentum.Phi()) ) > VH_deltaPhi  
// 	&& ( H.jets.at(0).csv > btag_csv_min && H.jets.at(1).csv > btag_csv_min )
// 	&& ( H.jets.at(0).csv > btag_csv_max || H.jets.at(1).csv > btag_csv_max )
// 	&& H.additionalJets.size() < nOfAddionalJets 
// 	&& e.VCandidate.additionalLeptons.size() < nOfAddionalLeptons
// 	&& TMath::Abs(H.pullAngle) < pullAngle 
// 	&& TMath::Abs(H.helicityAngle) < helicityAngle )
//       go = true;
//      return go;
//   }

// private:
  
//   Double_t btag_min;
//   Double_t btag_max;
//   Double_t Higgs_pt;
//   Double_t V_pt;
//   Double_t VH_deltaPhi;
//   Int_t nOfAdditionalJets;
//   Int_t nOfAdditionalLeptons;
//   Double_t pullAngle;
//   Double_t helicityAngle;

// };


// class ControlRegion_Boost_nobTag: public Cut {
//   std::string name() {return "ControlRegion_Boost_nobTag";};

//   VectorCandidate V = e.VHbbCandidate.VectorCandidate;
//   HiggsCandidate H = e.VHbbCandidate.HiggsCandidate;

//   if(V.isZLL){ 
//     btag_min = 0.;
//     btag_max = 0.;
//     Higgs_pt = 150;
//     V_pt = 150;
//     VH_deltaPhi = 2.70;
//     nOfAdditionalJets = 2;
//     pullAngle = 1.57;
//     helicityAngle = 0.8;
//   }
//   else if(V.isWLNu){
//     btag_min = 0.;
//     btag_max = 0.;
//     Higgs_pt = 150;
//     V_pt = 150;
//     VH_deltaPhi = 2.95;
//     nOfAdditionalJets = 2;
//     pullAngle = 1.57;
//     helicityAngle = 0.8;
//   }
//   else if(V.isZNuNu){
//     btag_min = 0.;
//     btag_max = 0.;
//     Higgs_pt = 150;
//     V_pt = 150;
//     VH_deltaPhi = 2.95;
//     nOfAdditionalJets = 2;
//     pullAngle = 1.57;
//     helicityAngle = 0.8;
//   }
//   else
//     std::cerr << "No vector boson reconstructed. 
// No histos will be filled." << std::endl;

//   Bool_t pass(VHbbEvent &e){
//     Bool_t go = false;
//     if( H.fourMomentum.Pt() > Higgs_pt 
// 	&& e.VCandidate.fourMomentum.Pt() > V_pt 
// 	&& TMath::Abs( Geom::deltaPhi(H.fourMomentum.Phi(), e.ZCandidate.fourMomentum.Phi()) ) > VH_deltaPhi  
// 	&& ( H.jets.at(0).csv > btag_csv_min && H.jets.at(1).csv > btag_csv_min )
// 	&& ( H.jets.at(0).csv > btag_csv_max || H.jets.at(1).csv > btag_csv_max )
// 	&& H.additionalJets.size() < nOfAddionalJets 
// 	&& e.VCandidate.additionalLeptons.size() < nOfAddionalLeptons
// 	&& TMath::Abs(H.pullAngle) < pullAngle 
// 	&& TMath::Abs(H.helicityAngle) < helicityAngle )
//       go = true;
//      return go;
//   }

// private:
  
//   Double_t btag_min;
//   Double_t btag_max;
//   Double_t Higgs_pt;
//   Double_t V_pt;
//   Double_t VH_deltaPhi;
//   Int_t nOfAdditionalJets;
//   Int_t nOfAdditionalLeptons;
//   Double_t pullAngle;
//   Double_t helicityAngle;

// };


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

    //from MC
    h_simHMass = new TH1F(("simHiggsMass"+suffix).c_str(),("Sim Higgs Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    h_simHPt = new TH1F(("simHiggsPt"+suffix).c_str(),("Sim Higgs Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    h_simZMass = new TH1F(("simZMass"+suffix).c_str(),("Sim Z Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    h_simZPt = new TH1F(("simZPt"+suffix).c_str(),("Sim Z Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    h_simWMass = new TH1F(("simWMass"+suffix).c_str(),("Sim W Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    h_simWPt = new TH1F(("simWPt"+suffix).c_str(),("Sim W Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );

    //Candidates
//     h_HMass = new TH1F(("HiggsMass"+suffix).c_str(),(" Higgs Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
//     h_HPt = new TH1F(("HiggsPt"+suffix).c_str(),(" Higgs Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
//     h_ZMass = new TH1F(("ZMass"+suffix).c_str(),(" Z Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
//     h_ZPt = new TH1F(("ZPt"+suffix).c_str(),(" Z Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
//     h_WMass = new TH1F(("WMass"+suffix).c_str(),(" W Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
//     h_WPt = new TH1F(("WPt"+suffix).c_str(),(" W Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );

    //    h_HHel = new TH1F(("HiggsHel"+suffix).c_str(),("Higgs helicity angle ("+suffix+")").c_str(), bin_hel, min_hel, max_hel );
  }
  
  virtual void fill(const VHbbEvent &e,float w) {
    //from MC

//     HiggsCandidate H = e.HiggsCandidate;
//     VectorCandidate V = e.VectorCandidate;
    
    h_simHMass->Fill(e.mcH.fourMomentum.M(), w); 
    h_simHPt->Fill(e.mcH.fourMomentum.Pt(), w); 
    h_simZMass->Fill(e.mcZ.fourMomentum.M(), w); 
    h_simZPt->Fill(e.mcZ.fourMomentum.Pt(), w); 
    h_simWMass->Fill(e.mcW.fourMomentum.M(), w); 
    h_simWPt->Fill(e.mcW.fourMomentum.Pt(), w); 

    //Candidates
//     h_HMass->Fill(H.fourMomentum.M(), w); 
//     h_HPt->Fill(H.fourMomentum.Pt(), w); 
//     if(V.isZLL || V.isZNuNu){
//       h_ZMass->Fill(V.fourMomentum.M(), w); 
//       h_ZPt->Fill(V.fourMomentum.Pt(), w);
//     } 
//     else if(V.isWLNu){
//       h_WMass->Fill(V.fourMomentum.M(), w); 
//       h_WPt->Fill(V.fourMomentum.Pt(), w); 
//     }

    //    h_HHel->Fill(e.Candidate.Higgs.hel(), w); 
  }
  

  TH1F * h_simHMass;
  TH1F * h_simHPt;
  TH1F * h_simWMass;
  TH1F * h_simWPt;
  TH1F * h_simZMass;
  TH1F * h_simZPt;
 
  TH1F * h_HMass;
  TH1F * h_HPt;
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

