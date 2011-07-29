#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbProxy.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/CutsAndHistos.h"
#include <TH1F.h>
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include <sstream>
#include "TKey.h"

//class LeptonSelection : public Cut 

class mcHPtCut :  public Cut {
public:
  mcHPtCut(Double_t mcHPtCutMin_):
    mcHPtCutMin(mcHPtCutMin_) {}
  mcHPtCut(Double_t mcHPtCutMin_, Double_t mcHPtCutMax_):
    mcHPtCutMin(mcHPtCutMin_),
    mcHPtCutMax(mcHPtCutMax_) {}
  std::string name() {
    oss_mcHPtCutMin << mcHPtCutMin;
    if(!mcHPtCutMax)
      return "mcHpT"+oss_mcHPtCutMin.str();
    else{
      oss_mcHPtCutMax << mcHPtCutMax;
      return "mcHpT"+oss_mcHPtCutMin.str()+"To"+oss_mcHPtCutMax.str();
    }
  }
  bool pass(VHbbProxy &iProxy) {
    if(!mcHPtCutMax)
      return ( iProxy.getVHbbEventAuxInfo()->mcH.size() !=0 && ((iProxy.getVHbbEventAuxInfo()->mcH)[0]).p4.Pt() > mcHPtCutMin );
    else
      return ( iProxy.getVHbbEventAuxInfo()->mcH.size() !=0 && ((iProxy.getVHbbEventAuxInfo()->mcH)[0]).p4.Pt() > mcHPtCutMin 
	       && iProxy.getVHbbEventAuxInfo()->mcH.size() !=0 && ((iProxy.getVHbbEventAuxInfo()->mcH)[0]).p4.Pt() < mcHPtCutMax );
  } 
private:
  Double_t mcHPtCutMin;
  Double_t mcHPtCutMax;
  std::ostringstream oss_mcHPtCutMin;
  std::ostringstream oss_mcHPtCutMax;
};


class mcWPtCut :  public Cut {
public:
  mcWPtCut(Double_t mcWPtCutMin_):
    mcWPtCutMin(mcWPtCutMin_) {}
  mcWPtCut(Double_t mcWPtCutMin_, Double_t mcWPtCutMax_):
    mcWPtCutMin(mcWPtCutMin_),
    mcWPtCutMax(mcWPtCutMax_) {}
  std::string name() {
    oss_mcWPtCutMin << mcWPtCutMin;
    if(!mcWPtCutMax)
      return "mcWpT"+oss_mcWPtCutMin.str();
    else{
      oss_mcWPtCutMax << mcWPtCutMax;
      return "mcWpT"+oss_mcWPtCutMin.str()+"To"+oss_mcWPtCutMax.str();
    }
  }
  bool pass(VHbbProxy &iProxy) {
    if(!mcWPtCutMax)
      return ( iProxy.getVHbbEventAuxInfo()->mcW.size() !=0 && ((iProxy.getVHbbEventAuxInfo()->mcW)[0]).p4.Pt() > mcWPtCutMin );
    else
      return ( iProxy.getVHbbEventAuxInfo()->mcW.size() !=0 && ((iProxy.getVHbbEventAuxInfo()->mcW)[0]).p4.Pt() > mcWPtCutMin 
	       && iProxy.getVHbbEventAuxInfo()->mcW.size() !=0 && ((iProxy.getVHbbEventAuxInfo()->mcW)[0]).p4.Pt() < mcWPtCutMax );
  } 
private:
  Double_t mcWPtCutMin;
  Double_t mcWPtCutMax;
  std::ostringstream oss_mcWPtCutMin;
  std::ostringstream oss_mcWPtCutMax;
};



class mcZPtCut :  public Cut {
public:
  mcZPtCut(Double_t mcZPtCutMin_):
    mcZPtCutMin(mcZPtCutMin_) {}
  mcZPtCut(Double_t mcZPtCutMin_, Double_t mcZPtCutMax_):
    mcZPtCutMin(mcZPtCutMin_),
    mcZPtCutMax(mcZPtCutMax_) {}
  std::string name() {
    oss_mcZPtCutMin << mcZPtCutMin;
    if(!mcZPtCutMax)
      return "mcZpT"+oss_mcZPtCutMin.str();
    else{
      oss_mcZPtCutMax << mcZPtCutMax;
      return "mcZpT"+oss_mcZPtCutMin.str()+"To"+oss_mcZPtCutMax.str();
    }
  }
  bool pass(VHbbProxy &iProxy) {
    if(!mcZPtCutMax)
      return ( iProxy.getVHbbEventAuxInfo()->mcZ.size() !=0 && ((iProxy.getVHbbEventAuxInfo()->mcZ)[0]).p4.Pt() > mcZPtCutMin );
    else
      return ( iProxy.getVHbbEventAuxInfo()->mcZ.size() !=0 && ((iProxy.getVHbbEventAuxInfo()->mcZ)[0]).p4.Pt() > mcZPtCutMin 
	       && iProxy.getVHbbEventAuxInfo()->mcZ.size() !=0 && ((iProxy.getVHbbEventAuxInfo()->mcZ)[0]).p4.Pt() < mcZPtCutMax );
  } 
private:
  Double_t mcZPtCutMin;
  Double_t mcZPtCutMax;
  std::ostringstream oss_mcZPtCutMin;
  std::ostringstream oss_mcZPtCutMax;
};


class SimpleJet1PtCut: public Cut {
public:
  SimpleJet1PtCut(Double_t simpleJet1PtCutMin_):
    simpleJet1PtCutMin(simpleJet1PtCutMin_) {}
  SimpleJet1PtCut(Double_t simpleJet1PtCutMin_, Double_t simpleJet1PtCutMax_):
    simpleJet1PtCutMin(simpleJet1PtCutMin_),
    simpleJet1PtCutMax(simpleJet1PtCutMax_) {}
  std::string name() {
    oss_simpleJet1PtCutMin << simpleJet1PtCutMin;
    if(!simpleJet1PtCutMax){
      return "SimpleJet1PtCut"+oss_simpleJet1PtCutMin.str();
    }
    else{
      oss_simpleJet1PtCutMax << simpleJet1PtCutMax;
      return "SimpleJet1PtCut"+oss_simpleJet1PtCutMin.str()+"To"+oss_simpleJet1PtCutMax.str();
    }
  }
  bool pass(VHbbProxy &iProxy) {
    if(iProxy.getVHbbCandidate()->size() < 1)
      return false;
    else 
      if(!simpleJet1PtCutMin)
	return ((iProxy.getVHbbCandidate()->at(0)).H.jets.at(0).p4.Pt() > simpleJet1PtCutMin);
      else
	return ((iProxy.getVHbbCandidate()->at(0)).H.jets.at(0).p4.Pt() > simpleJet1PtCutMin
		&& (iProxy.getVHbbCandidate()->at(0)).H.jets.at(0).p4.Pt() < simpleJet1PtCutMax);
  } 
private:
  Double_t simpleJet1PtCutMin;
  Double_t simpleJet1PtCutMax;
  std::ostringstream  oss_simpleJet1PtCutMin;
  std::ostringstream  oss_simpleJet1PtCutMax;
};

class VPtCut: public Cut {
public:
  VPtCut(Double_t VptCutMin_):
    VptCutMin(VptCutMin_) {}
  VPtCut(Double_t VptCutMin_, Double_t VptCutMax_):
    VptCutMin(VptCutMin_),
    VptCutMax(VptCutMax_) {}
  std::string name() {
    oss_VptMin << VptCutMin;
    if(! VptCutMax ){
      return "VpT"+oss_VptMin.str();
    }
    else{
      oss_VptMax << VptCutMax;
      return "VpT"+oss_VptMin.str()+"To"+oss_VptMax.str();
    }
  }
  bool pass(VHbbProxy &iProxy) {
    if(iProxy.getVHbbCandidate()->size() < 1)
      return false;
    else 
      if(VptCutMax != 1e10)
	return ((iProxy.getVHbbCandidate()->at(0)).V.p4.Pt() > VptCutMin 
		&& (iProxy.getVHbbCandidate()->at(0)).V.p4.Pt() < VptCutMax);
      else
	return ((iProxy.getVHbbCandidate()->at(0)).V.p4.Pt() > VptCutMin );
		
  } 
private:
  Double_t VptCutMin;
  Double_t VptCutMax;
  std::ostringstream oss_VptMin;
  std::ostringstream oss_VptMax;
};


class HPtCut: public Cut {
public:
  HPtCut(Double_t hPtCutMin_):
    hPtCutMin(hPtCutMin_) {}
  HPtCut(Double_t hPtCutMin_, Double_t hPtCutMax_):
    hPtCutMin(hPtCutMin_),
    hPtCutMax(hPtCutMax_) {}

  std::string name() {
    oss_HPtCutMin << hPtCutMin;
    if(!hPtCutMax){
      return "HpT"+oss_HPtCutMin.str();
    }
    else{
      oss_HPtCutMax << hPtCutMax;
      return "HpT"+oss_HPtCutMin.str()+"To"+oss_HPtCutMax.str();
    }
  }
  bool pass(VHbbProxy &iProxy) {
    if(iProxy.getVHbbCandidate()->size() < 1)
      return false;
    else 
      if(!hPtCutMax)
	return ((iProxy.getVHbbCandidate()->at(0)).H.p4.Pt() > hPtCutMin);
      else 
	return ((iProxy.getVHbbCandidate()->at(0)).H.p4.Pt() > hPtCutMin
		&& (iProxy.getVHbbCandidate()->at(0)).H.p4.Pt() < hPtCutMax);
  }
private:
  Double_t hPtCutMin;
  Double_t hPtCutMax;
  std::ostringstream oss_HPtCutMin;
  std::ostringstream oss_HPtCutMax;
};


class SignalRegionOld: public Cut {
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
    if( H.p4.Pt() > Higgs_pt 
	&& V.p4.Pt() > V_pt 
	&& TMath::Abs( Geom::deltaPhi(H.p4.Phi(), V.p4.Phi()) ) > VH_deltaPhi  
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
    McH_simHMass = new TH1F(("simHiggsMass"+suffix).c_str(),("Sim Higgs Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    McH_simHPt = new TH1F(("simHiggsPt"+suffix).c_str(),("Sim Higgs Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    McH_simZMass = new TH1F(("simZMass"+suffix).c_str(),("Sim Z Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    McH_simZPt = new TH1F(("simZPt"+suffix).c_str(),("Sim Z Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    McH_simWMass = new TH1F(("simWMass"+suffix).c_str(),("Sim W Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    McH_simWPt = new TH1F(("simWPt"+suffix).c_str(),("Sim W Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );

  }
  
  virtual void fill(VHbbProxy &iProxy,float w) {

    //    const VHbbEvent *iEvent = iProxy.getVHbbEvent();  
    const VHbbEventAuxInfo *iAuxInfo = iProxy.getVHbbEventAuxInfo();
   if(iAuxInfo) {
    //from MC    
     if (iAuxInfo->mcH.size() )
       McH_simHMass->Fill(iAuxInfo->mcH[0].p4.M(), w); 
     if (iAuxInfo->mcH.size() )
    McH_simHPt->Fill(iAuxInfo->mcH[0].p4.Pt(), w); 
     if (iAuxInfo->mcZ.size() )
    McH_simZMass->Fill(iAuxInfo->mcZ[0].p4.M(), w); 
     if (iAuxInfo->mcZ.size() )
    McH_simZPt->Fill(iAuxInfo->mcZ[0].p4.Pt(), w); 
     if (iAuxInfo->mcW.size() )
    McH_simWMass->Fill(iAuxInfo->mcW[0].p4.M(), w); 
     if (iAuxInfo->mcW.size() )
    McH_simWPt->Fill(iAuxInfo->mcW[0].p4.Pt(), w); 
   }
  }


  TH1F * McH_simHMass;
  TH1F * McH_simHPt;
  TH1F * McH_simWMass;
  TH1F * McH_simWPt;
  TH1F * McH_simZMass;
  TH1F * McH_simZPt;
 
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


class BTagHistos : public Histos {

public:
  
  virtual void book(TFile &f, std::string suffix) {
    
    TDirectory *subDir;
    if( ! f.GetDirectory(suffix.c_str()) )
      subDir = f.mkdir(suffix.c_str());
    else 
      subDir = f.GetDirectory(suffix.c_str());
    subDir->cd();
    
    bin_btag_prob = 20;
    min_btag_prob = 0;
    max_btag_prob = 1;

    bin_btag_count = 10;
    min_btag_count = 0;
    max_btag_count = 10;

    //Candidates
    BTagH_bJet1_csv = new TH1F(("BJet1_CSV"+suffix).c_str(),("BJet1 CSV ("+suffix+")").c_str(), bin_btag_prob, min_btag_prob, max_btag_prob );
    BTagH_bJet2_csv = new TH1F(("BJet2_CSV"+suffix).c_str(),("BJet2 CSV ("+suffix+")").c_str(), bin_btag_prob, min_btag_prob, max_btag_prob );
    BTagH_bJet1_csvmva = new TH1F(("BJet1_CSVMVA"+suffix).c_str(),("BJet1 CSVMVA ("+suffix+")").c_str(), bin_btag_prob, min_btag_prob, max_btag_prob );
    BTagH_bJet2_csvmva = new TH1F(("BJet2_CSVMVA"+suffix).c_str(),("BJet2 CSVMVA ("+suffix+")").c_str(), bin_btag_prob, min_btag_prob, max_btag_prob );
    BTagH_bJet1_ssvhe = new TH1F(("BJet1_SSVHE"+suffix).c_str(),("BJet1 SSVHE ("+suffix+")").c_str(), bin_btag_prob, min_btag_prob, max_btag_prob );
    BTagH_bJet2_ssvhe = new TH1F(("BJet2_SSVHE"+suffix).c_str(),("BJet2 SSVHE ("+suffix+")").c_str(), bin_btag_prob, min_btag_prob, max_btag_prob );
    BTagH_bJet1_jpb = new TH1F(("BJet1_JPB"+suffix).c_str(),("BJet1 JPB ("+suffix+")").c_str(), bin_btag_prob, min_btag_prob, max_btag_prob );
    BTagH_bJet2_jpb = new TH1F(("BJet2_JPB"+suffix).c_str(),("BJet2 JPB ("+suffix+")").c_str(), bin_btag_prob, min_btag_prob, max_btag_prob );
    BTagH_bJet1_tche = new TH1F(("BJet1_TCHE"+suffix).c_str(),("BJet1 TCHE ("+suffix+")").c_str(), bin_btag_count, min_btag_count, max_btag_count );
    BTagH_bJet2_tche = new TH1F(("BJet2_TCHE"+suffix).c_str(),("BJet2 TCHE ("+suffix+")").c_str(), bin_btag_count, min_btag_count, max_btag_count );
    BTagH_bJet1_jp = new TH1F(("BJet1_JP"+suffix).c_str(),("BJet1 JP ("+suffix+")").c_str(), bin_btag_prob, min_btag_prob, max_btag_prob );
    BTagH_bJet2_jp = new TH1F(("BJet2_JP"+suffix).c_str(),("BJet2 JP ("+suffix+")").c_str(), bin_btag_prob, min_btag_prob, max_btag_prob );
    BTagH_bJet1_tchp = new TH1F(("BJet1_TCHP"+suffix).c_str(),("BJet1 TCHP ("+suffix+")").c_str(), bin_btag_count, min_btag_count, max_btag_count );
    BTagH_bJet2_tchp = new TH1F(("BJet2_TCHP"+suffix).c_str(),("BJet2 TCHP ("+suffix+")").c_str(), bin_btag_count, min_btag_count, max_btag_count );

  }
  
  virtual void fill(VHbbProxy &iProxy,float w) {

    const VHbbEvent *iEvent = iProxy.getVHbbEvent();
    const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();

    //Candidates
    if(iCand->size() > 0){
      VHbbCandidate::CandidateType iCandType = iCand->at(0).candidateType;
      VHbbCandidate::HiggsCandidate H = iCand->at(0).H;
      
      BTagH_bJet1_csv->Fill(H.jets.at(0).csv, w);
      BTagH_bJet2_csv->Fill(H.jets.at(1).csv, w);
      BTagH_bJet1_csvmva->Fill(H.jets.at(0).csvmva, w);
      BTagH_bJet2_csvmva->Fill(H.jets.at(1).csvmva, w);
      BTagH_bJet1_ssvhe->Fill(H.jets.at(0).ssvhe, w);
      BTagH_bJet2_ssvhe->Fill(H.jets.at(1).ssvhe, w);
      BTagH_bJet1_tche->Fill(H.jets.at(0).tche, w);
      BTagH_bJet2_tche->Fill(H.jets.at(1).tche, w);
      BTagH_bJet1_tchp->Fill(H.jets.at(0).tchp, w);
      BTagH_bJet2_tchp->Fill(H.jets.at(1).tchp, w);
      BTagH_bJet1_jpb->Fill(H.jets.at(0).jpb, w);
      BTagH_bJet2_jpb->Fill(H.jets.at(1).jpb, w);
      BTagH_bJet1_jp->Fill(H.jets.at(0).jp, w);
      BTagH_bJet2_jp->Fill(H.jets.at(1).jp, w);
      
    }
  }
  
  TH1F * BTagH_bJet1_csv;
  TH1F * BTagH_bJet2_csv;
  TH1F * BTagH_bJet1_csvmva;
  TH1F * BTagH_bJet2_csvmva;
  TH1F * BTagH_bJet1_ssvhe;
  TH1F * BTagH_bJet2_ssvhe;
  TH1F * BTagH_bJet1_jpb;
  TH1F * BTagH_bJet2_jpb;
  TH1F * BTagH_bJet1_tche;
  TH1F * BTagH_bJet2_tche;
  TH1F * BTagH_bJet1_jp;
  TH1F * BTagH_bJet2_jp;
  TH1F * BTagH_bJet1_tchp;
  TH1F * BTagH_bJet2_tchp;
 
private:

  Int_t bin_btag_prob;
  Double_t min_btag_prob;
  Double_t max_btag_prob;
  
  Int_t bin_btag_count;
  Double_t min_btag_count;
  Double_t max_btag_count;
  
};
 
class CountHisto : public Histos {
   virtual void book(TFile &f, std::string suffix) { 
      
    TDirectory *subDir;
    if( ! f.GetDirectory(suffix.c_str()) )
      subDir = f.mkdir(suffix.c_str());
    else 
      subDir = f.GetDirectory(suffix.c_str());
    subDir->cd();
     count = new TH1F(("Count"+suffix).c_str(),("Count ("+suffix+")").c_str(), 1,0,2 );

   }
   virtual void fill(VHbbProxy &iProxy,float w) {
    count->Fill(1,w);
   }

   TH1F * count;
};
 
class StandardHistos : public Histos {
    
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

    bin_btag = 20;
    min_btag = 0;
    max_btag = 1;

    bin_deltaR = bin_deltaPhi = bin_deltaEta = 20;
    min_deltaR = min_deltaPhi = min_deltaEta = 0;
    max_deltaR = max_deltaPhi = max_deltaEta = 5;

    //Candidates
    StH_simpleJet1_pt = new TH1F(("SimpleJet1_pt"+suffix).c_str(),("Simple Jet1 pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    StH_simpleJet2_pt = new TH1F(("SimpleJet2_pt"+suffix).c_str(),("Simple Jet2 pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    StH_simpleJet1_bTag = new TH1F(("SimpleJet1_bTag"+suffix).c_str(),("Simple Jet1 bTag ("+suffix+")").c_str(), bin_btag, min_btag, max_btag );
    StH_simpleJet2_bTag = new TH1F(("SimpleJet2_bTag"+suffix).c_str(),("Simple Jet2 bTag ("+suffix+")").c_str(), bin_btag, min_btag, max_btag );
    StH_simpleJets_dR = new TH1F(("SimpleJets_dR"+suffix).c_str(),("Simple Jets deltaR ("+suffix+")").c_str(), bin_deltaR, min_deltaR, max_deltaR );
    StH_simpleJets_dPhi = new TH1F(("SimpleJets_dPhi"+suffix).c_str(),("Simple Jets deltaPhi ("+suffix+")").c_str(), bin_deltaPhi, min_deltaPhi, max_deltaPhi );
    StH_simpleJets_dEta = new TH1F(("SimpleJets_dEta"+suffix).c_str(),("Simple Jets deltaEta ("+suffix+")").c_str(), bin_deltaEta, min_deltaEta, max_deltaEta );

    StH_HMass = new TH1F(("HiggsMass"+suffix).c_str(),(" Higgs Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    StH_HPt = new TH1F(("HiggsPt"+suffix).c_str(),(" Higgs Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    StH_HHel = new TH1F(("HiggsHel"+suffix).c_str(),("Higgs helicity angle ("+suffix+")").c_str(), bin_hel, min_hel, max_hel );
    StH_HPullAngle = new TH1F(("HiggsPullAngle"+suffix).c_str(),("Higgs pull angle ("+suffix+")").c_str(), bin_deltaPhi, min_deltaPhi, max_deltaPhi );

    StH_ZMass = new TH1F(("ZMass"+suffix).c_str(),(" Z Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    StH_ZPt = new TH1F(("ZPt"+suffix).c_str(),(" Z Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    StH_ZH_dPhi = new TH1F(("ZH_dPhi"+suffix).c_str(),(" ZH delta Phi ("+suffix+")").c_str(), bin_deltaPhi, min_deltaPhi, max_deltaPhi );

    StH_WMass = new TH1F(("WMass"+suffix).c_str(),(" W Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    StH_WPt = new TH1F(("WPt"+suffix).c_str(),(" W Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    StH_WH_dPhi = new TH1F(("WH_dPhi"+suffix).c_str(),(" WH delta Phi ("+suffix+")").c_str(), bin_deltaPhi, min_deltaPhi, max_deltaPhi );

  }
  
  virtual void fill(VHbbProxy &iProxy,float w) {

    const VHbbEvent *iEvent = iProxy.getVHbbEvent();
    const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();

    //Candidates
    if(iCand->size() > 0){
      VHbbCandidate::CandidateType iCandType = iCand->at(0).candidateType;
      VHbbCandidate::HiggsCandidate H = iCand->at(0).H;
      VHbbCandidate::VectorCandidate V = iCand->at(0).V;
      
      StH_simpleJet1_pt->Fill(H.jets.at(0).p4.Pt(), w);
      StH_simpleJet2_pt->Fill(H.jets.at(1).p4.Pt(), w);
      StH_simpleJet1_bTag->Fill(H.jets.at(0).csv, w);
      StH_simpleJet2_bTag->Fill(H.jets.at(1).csv, w);
      StH_simpleJets_dR->Fill(H.jets.at(0).p4.DeltaR(H.jets.at(1).p4), w);
      StH_simpleJets_dPhi->Fill(H.jets.at(0).p4.DeltaPhi(H.jets.at(1).p4), w);
      StH_simpleJets_dEta->Fill(TMath::Abs(H.jets.at(0).p4.Eta()-H.jets.at(1).p4.Eta()), w);

      StH_HMass->Fill(H.p4.M(), w); 
      StH_HPt->Fill(H.p4.Pt(), w); 
     //    StH_HHel->Fill(H.hel(), w); 
      StH_HPullAngle->Fill(H.deltaTheta, w); 
      if( iCandType == VHbbCandidate::Zmumu || iCandType == VHbbCandidate::Zee || iCandType == VHbbCandidate::Znn ){
	StH_ZMass->Fill(V.p4.M(), w); 
	StH_ZPt->Fill(V.p4.Pt(), w);
	StH_ZH_dPhi->Fill(V.p4.DeltaPhi(H.p4.Phi()), w); 
      } 
      else if(iCandType == VHbbCandidate::Wen || iCandType == VHbbCandidate::Wmun){
	StH_WMass->Fill(V.p4.M(), w); 
	StH_WPt->Fill(V.p4.Pt(), w); 
	StH_WH_dPhi->Fill(V.p4.DeltaPhi(H.p4.Phi()), w); 
      }
 
    }
  }
  
  TH1F * StH_simpleJet1_pt;
  TH1F * StH_simpleJet2_pt;
  TH1F * StH_simpleJet1_bTag;
  TH1F * StH_simpleJet2_bTag;
  TH1F * StH_simpleJets_dR;
  TH1F * StH_simpleJets_dPhi;
  TH1F * StH_simpleJets_dEta;

  TH1F * StH_HMass;
  TH1F * StH_HPt;
  TH1F * StH_HHel;
  TH1F * StH_HPullAngle;
  TH1F * StH_WMass;
  TH1F * StH_WPt;
  TH1F * StH_WH_dPhi;
  TH1F * StH_ZMass;
  TH1F * StH_ZPt;
  TH1F * StH_ZH_dPhi;
 
private:

  Int_t bin_btag;
  Double_t min_btag;
  Double_t max_btag;

  Int_t bin_deltaEta;
  Double_t min_deltaEta;
  Double_t max_deltaEta;

  Int_t bin_deltaPhi;
  Double_t min_deltaPhi;
  Double_t max_deltaPhi;

  Int_t bin_deltaR;
  Double_t min_deltaR;
  Double_t max_deltaR;

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

class HardJetHistos : public Histos {
    
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

    bin_btag = 20;
    min_btag = 0;
    max_btag = 1;

    bin_deltaR = bin_deltaPhi = bin_deltaEta = 20;
    min_deltaR = min_deltaPhi = min_deltaEta = 0;
    max_deltaR = max_deltaPhi = max_deltaEta = 5;

    //Candidates

    HardJetH_subJet1_pt = new TH1F(("SubJet1_pt"+suffix).c_str(),("Sub Jet1 pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    HardJetH_subJet2_pt = new TH1F(("SubJet2_pt"+suffix).c_str(),("Sub Jet2 pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    HardJetH_subJet1_bTag = new TH1F(("SubJet1_bTag"+suffix).c_str(),("Sub Jet1 bTag ("+suffix+")").c_str(), bin_btag, min_btag, max_btag );
    HardJetH_subJet2_bTag = new TH1F(("SubJet2_bTag"+suffix).c_str(),("Sub Jet2 bTag ("+suffix+")").c_str(), bin_btag, min_btag, max_btag );
    HardJetH_subJets_dR = new TH1F(("SubJets_dR"+suffix).c_str(),("Sub Jets deltaR ("+suffix+")").c_str(), bin_deltaR, min_deltaR, max_deltaR );
    HardJetH_subJets_dPhi = new TH1F(("SubJets_dPhi"+suffix).c_str(),("Sub Jets deltaPhi ("+suffix+")").c_str(), bin_deltaPhi, min_deltaPhi, max_deltaPhi );
    HardJetH_subJets_dEta = new TH1F(("SubJets_dEta"+suffix).c_str(),("Sub Jets deltaEta ("+suffix+")").c_str(), bin_deltaEta, min_deltaEta, max_deltaEta );

    HardJetH_HMass = new TH1F(("HiggsMass"+suffix).c_str(),(" Higgs Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    HardJetH_HPt = new TH1F(("HiggsPt"+suffix).c_str(),(" Higgs Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    HardJetH_HHel = new TH1F(("HiggsHel"+suffix).c_str(),("Higgs helicity angle ("+suffix+")").c_str(), bin_hel, min_hel, max_hel );

    HardJetH_ZMass = new TH1F(("ZMass"+suffix).c_str(),(" Z Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    HardJetH_ZPt = new TH1F(("ZPt"+suffix).c_str(),(" Z Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    HardJetH_ZH_dPhi = new TH1F(("ZH_dPhi"+suffix).c_str(),(" ZH delta Phi ("+suffix+")").c_str(), bin_deltaPhi, min_deltaPhi, max_deltaPhi );

    HardJetH_WMass = new TH1F(("WMass"+suffix).c_str(),(" W Mass ("+suffix+")").c_str(), bin_mass, min_mass, max_mass );
    HardJetH_WPt = new TH1F(("WPt"+suffix).c_str(),(" W Pt ("+suffix+")").c_str(), bin_pt, min_pt, max_pt );
    HardJetH_WH_dPhi = new TH1F(("WH_dPhi"+suffix).c_str(),(" WH delta Phi ("+suffix+")").c_str(), bin_deltaPhi, min_deltaPhi, max_deltaPhi );

  }
  
  virtual void fill(VHbbProxy &iProxy,float w) {

    const VHbbEvent *iEvent = iProxy.getVHbbEvent();
if(iEvent)
{    const std::vector<VHbbCandidate> *iCand = iProxy.getVHbbCandidate();

    //Candidates
    if(iCand->size() > 0){
      VHbbCandidate::CandidateType iCandType = iCand->at(0).candidateType;
      VHbbCandidate::HiggsCandidate H = iCand->at(0).H;
      VHbbCandidate::VectorCandidate V = iCand->at(0).V;
      std::vector<VHbbEvent::HardJet> iHardJets = iEvent->hardJets;
      VHbbEvent::HardJet iHardJet = iHardJets.at(0);

      HardJetH_subJet1_pt->Fill(iHardJet.subFourMomentum.at(0).Pt(), w);
      HardJetH_subJet2_pt->Fill(iHardJet.subFourMomentum.at(1).Pt(), w);
      //SubJet information on btag missing
//       HardJetH_subJet1_bTag->Fill(iHardJet.at(0).csv, w);
//       HardJetH_subJet2_bTag->Fill(iHardJet.at(1).csv, w);
      HardJetH_subJets_dR->Fill(iHardJet.subFourMomentum.at(0).DeltaR(iHardJet.subFourMomentum.at(1)), w);
      HardJetH_subJets_dPhi->Fill(iHardJet.subFourMomentum.at(0).DeltaPhi(iHardJet.subFourMomentum.at(1)), w);
      HardJetH_subJets_dEta->Fill(TMath::Abs(iHardJet.subFourMomentum.at(0).Eta()-iHardJet.subFourMomentum.at(1).Eta()), w);

      //Here there should be the higgs candidate from HardJet
//       HardJetH_HMass->Fill(H.p4.M(), w); 
//       HardJetH_HPt->Fill(H.p4.Pt(), w); 
//      //    HardJetH_HHel->Fill(H.hel(), w); 
//       if( iCandType == VHbbCandidate::Zmumu || iCandType == VHbbCandidate::Zee || iCandType == VHbbCandidate::Znn ){
// 	HardJetH_ZMass->Fill(V.p4.M(), w); 
// 	HardJetH_ZPt->Fill(V.p4.Pt(), w);
// 	HardJetH_ZH_dPhi->Fill(V.p4.DeltaPhi(H.p4.Phi()), w); 
//       } 
//       else if(iCandType == VHbbCandidate::Wen || iCandType == VHbbCandidate::Wmun){
// 	HardJetH_WMass->Fill(V.p4.M(), w); 
// 	HardJetH_WPt->Fill(V.p4.Pt(), w); 
// 	HardJetH_WH_dPhi->Fill(V.p4.DeltaPhi(H.p4.Phi()), w); 
//       }
 }
    }
  }
  
  TH1F * HardJetH_subJet1_pt;
  TH1F * HardJetH_subJet2_pt;
  TH1F * HardJetH_subJet1_bTag;
  TH1F * HardJetH_subJet2_bTag;
  TH1F * HardJetH_subJets_dR;
  TH1F * HardJetH_subJets_dPhi;
  TH1F * HardJetH_subJets_dEta;

  TH1F * HardJetH_HMass;
  TH1F * HardJetH_HPt;
  TH1F * HardJetH_HHel;
  TH1F * HardJetH_WMass;
  TH1F * HardJetH_WPt;
  TH1F * HardJetH_WH_dPhi;
  TH1F * HardJetH_ZMass;
  TH1F * HardJetH_ZPt;
  TH1F * HardJetH_ZH_dPhi;

private:

  Int_t bin_btag;
  Double_t min_btag;
  Double_t max_btag;

  Int_t bin_deltaEta;
  Double_t min_deltaEta;
  Double_t max_deltaEta;

  Int_t bin_deltaPhi;
  Double_t min_deltaPhi;
  Double_t max_deltaPhi;

  Int_t bin_deltaR;
  Double_t min_deltaR;
  Double_t max_deltaR;

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
