#ifndef VHBBCANDIDATETOOLS_H
#define VHBBCANDIDATETOOLS_H

#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"

#include <iostream>

struct CompareJetPtMuons {
  bool operator()( const VHbbEvent::MuonInfo& j1, const  VHbbEvent::MuonInfo& j2 ) const {
    return j1.p4.Pt() > j2.p4.Pt();
  }
};
struct CompareJetPtElectrons {
  bool operator()( const VHbbEvent::ElectronInfo& j1, const  VHbbEvent::ElectronInfo& j2 ) const {
    return j1.p4.Pt() > j2.p4.Pt();
  }
};



class VHbbCandidateTools {
 public:

 VHbbCandidateTools(bool verbose = false): verbose_(verbose){}

  float deltaPhi(float in2, float in1){
    float dphi = in2-in1;
    if ( dphi > M_PI ) {
      dphi -= 2.0*M_PI;
    } else if ( dphi <= -M_PI ) {
      dphi += 2.0*M_PI;
    }
    return dphi;
  }
  
  
  VHbbCandidate getHZtaumuCandidate(const VHbbCandidate & in, bool & ok, std::vector<unsigned int>& muPos, std::vector<unsigned int>& tauPos){
    if (verbose_) std::cout <<" getHZtaumuCandidate input mu "<<in.V.muons.size()<<" e "<<in.V.electrons.size()<< " tau " << in.V.taus.size() << std::endl;
    ok = false;
    VHbbCandidate temp=in;
    // Require exactly one tau and muon, and no electrons
    if (temp.V.taus.size()!=1) return in;
    if (temp.V.muons.size()!=1) return in ;
    if (temp.V.electrons.size()!=0) return in ;
    temp.V.p4 = temp.V.taus[0].p4 + temp.V.muons[0].p4;
    temp.V.firstLepton = muPos[0];
    temp.V.secondLepton = tauPos[0];
    
    ok = true;
    return temp;
    
  }
  VHbbCandidate getHWtaunCandidate(const VHbbCandidate & in, bool & ok , std::vector<unsigned int>& pos){
    
    if (verbose_){
      std::cout <<" getHWtaunCandidate input mu "<<in.V.muons.size()<<" e "<<in.V.electrons.size()<< " tau " << in.V.taus.size() << std::endl;
      std::cout << " pos.size()=" << pos.size() << std::endl;
    }
    
    ok = false;
    VHbbCandidate temp=in;
    // require a tau and no electrons or muons
    if (temp.V.taus.size()!=1) return in;
    if (temp.V.muons.size()!=0) return in ;
    if (temp.V.electrons.size()!=0) return in ;
    if (temp.V.mets.size()<1) return in ;
    temp.V.p4 = temp.V.taus[0].p4+temp.V.mets[0].p4;
    temp.V.firstLepton=pos[0];
    
    if (verbose_){
      std::cout << "Done with getHWtaunCandidate" << std::endl;
    }
    
    ok=true;
    return temp;
  }
 

  
  VHbbCandidate getHZmumuCandidate(const VHbbCandidate & in, bool & ok, std::vector<unsigned int>& pos){
    if (verbose_){
      std::cout <<" getHZmumuCandidate input mu "<<in.V.muons.size()<<" e "<<in.V.electrons.size()<<std::endl;
    }
    ok = false;
    VHbbCandidate temp=in;
    
    //
    // change: allow for additional leptons; by definition 
    //
    if (temp.V.muons.size()<2) return in ;
    //    if (temp.V.electrons.size()!=0) return in ;
    std::vector<VHbbEvent::MuonInfo> muons_ = temp.V.muons;

    // beware: assumes already sorted!!!!

    //    CompareJetPtMuons ptComparator;
    //    std::sort(muons_.begin(), muons_.end(), ptComparator);
    if (muons_[0].p4.Pt()<20 || muons_[1].p4.Pt()<20 ) return in;
    
    //
    // now I need to ask also for the charge
    //
    int selectMu2=1;
    
    if (muons_[0].charge* muons_[selectMu2].charge > 0){
      if (muons_.size() ==2) return in;
      //
      // i need to find a proper pair
      //
      
      for (unsigned int it=2; it< muons_.size(); ++it){
	if (  muons_[it].charge * muons_[0].charge < 0) {
	  selectMu2 = it;
	  break;
	}
	if (selectMu2 == 1) return in;
      }  
    }
    temp.V.p4 = muons_[0].p4+muons_[selectMu2].p4;
    std::vector<VHbbEvent::MuonInfo> muons2_;
    for (std::vector<VHbbEvent::MuonInfo>::const_iterator it = muons_.begin(); it!= muons_.end(); ++it){
      if (it->p4.Pt()>20) muons2_.push_back(*it);
    }
    temp.V.muons = muons2_;
    
    // the same for electrons
    std::vector<VHbbEvent::ElectronInfo> electrons_ = temp.V.electrons;

    // beware; assumes already sorted

    //    CompareJetPtElectrons ptComparator2;
    //    std::sort(electrons_.begin(), electrons_.end(), ptComparator2);
    std::vector<VHbbEvent::ElectronInfo> electrons2_;
    for (std::vector<VHbbEvent::ElectronInfo>::const_iterator it = electrons_.begin(); it!= electrons_.end(); ++it){
      if (it->p4.Pt()>20) electrons2_.push_back(*it);
    }
    temp.V.electrons = electrons2_;
    
    //
    // consider all 
    //
    
    
    //    if (temp.V.Pt()<150 ) return in;
    //    if (temp.H.Pt()<150) return in;
    //    if (temp.H.firstJet().csv< 0.9) return in;
    //    if (temp.H.secondJet().csv<0.5) return in;
    //    if (deltaPhi(temp.V.Phi(),temp.H.Phi())<2.7) return in;
    //    if (temp.V.FourMomentum.Mass()<75 || temp.V.FourMomentum.Mass()>105) return in; 
    //    if (temp.additionalJets.size()>0) return in;
    //    if (std::Abs(deltaTheta) ????

    temp.V.firstLepton = 0;
    temp.V.secondLepton = selectMu2;
    temp.V.firstLeptonOrig = pos[0];
    temp.V.secondLeptonOrig = pos[selectMu2];

   ok = true;
    return temp;
  }  
  VHbbCandidate getHZeeCandidate(const VHbbCandidate & in, bool & ok, std::vector<unsigned int>& pos){
    if (verbose_){
      std::cout <<" getHZeeCandidate input mu "<<in.V.electrons.size()<<" e "<<in.V.muons.size()<<std::endl;
    }
    ok = false;
    VHbbCandidate temp=in;
    
    //
    // change: allow for additional leptons; by definition 
    //
    if (temp.V.electrons.size()<2) return in ;
    //    if (temp.V.electrons.size()!=0) return in ;
    std::vector<VHbbEvent::ElectronInfo> electrons_ = temp.V.electrons;

    // beware assumes already sorted

    //    CompareJetPtElectrons ptComparator;
    //    std::sort(electrons_.begin(), electrons_.end(), ptComparator);
    if (electrons_[0].p4.Pt()<20 || electrons_[1].p4.Pt()<20 ) return in;
    //
    // now I need to ask also for the charge
    //
    int selectE2=1;
    
    if (electrons_[0].charge* electrons_[selectE2].charge > 0){
      if (electrons_.size() ==2) return in;
      //
      // i need to find a proper pair
      //
      
      for (unsigned int it=2; it< electrons_.size(); ++it){
	if (  electrons_[it].charge * electrons_[0].charge < 0) {
	  selectE2 = it;
	  break;
	}
	if (selectE2 == 1) return in;
      }  
    }
   temp.V.p4 = electrons_[0].p4+electrons_[selectE2].p4;
    std::vector<VHbbEvent::ElectronInfo> electrons2_;
    for (std::vector<VHbbEvent::ElectronInfo>::const_iterator it = electrons_.begin(); it!= electrons_.end(); ++it){
      if (it->p4.Pt()>20) electrons2_.push_back(*it);
    }
    temp.V.electrons = electrons2_;
  
    // the same for muonss
    std::vector<VHbbEvent::MuonInfo> muons_ = temp.V.muons;

    // beware assumes already sorted

    //    CompareJetPtMuons ptComparator2;
    //    std::sort(muons_.begin(), muons_.end(), ptComparator2);
    std::vector<VHbbEvent::MuonInfo> muons2_;
    for (std::vector<VHbbEvent::MuonInfo>::const_iterator it = muons_.begin(); it!= muons_.end(); ++it){
      if (it->p4.Pt()>20) muons2_.push_back(*it);
    }
    temp.V.muons = muons2_;

    temp.V.firstLepton = 0;
    temp.V.secondLepton = selectE2;
    temp.V.firstLeptonOrig = pos[0];
    temp.V.secondLeptonOrig = pos[selectE2];

    
   ok = true;
    return temp;
  }  
  VHbbCandidate getHZnnCandidate(const VHbbCandidate & in, bool & ok){
    if (verbose_){
      std::cout <<" getHZnnCandidate input mu "<<in.V.muons.size()<<" e "<<in.V.electrons.size()<<" met "<<in.V.mets.size()<<std::endl;
    }
    
    ok = false;
    VHbbCandidate temp=in;
    if (temp.V.mets.size()!=1) return in;
//always build a NuNu candidate, exclusion come from if/else series on the caller side
// this allow to still build candidates for mu+e  and same sign dilept
//    if (temp.V.muons.size()!=0) return in ;
//    if (temp.V.electrons.size()!=0) return in ;
    
    temp.V.p4 = temp.V.mets[0].p4;
    if (verbose_) {
      std::cout <<" debug met "<< temp.V.mets[0].metSig << " " <<  temp.V.mets[0].sumEt<< std::endl;
    }   
//    if (temp.V.mets[0].metSig<5) return in;
    if (temp.V.mets[0].p4.Pt()<80) return in;
    //    if (temp.H.p4.Pt()<150)return in;
    //    if (temp.H.firstJet().csv< 0.9) return in;
    //    if (temp.H.secondJet().csv<0.5) return in;
    //    if (deltaPhi(temp.V.p4.Phi(),temp.H.p4.Phi())<2.95) return in;
    //    if (temp.V.electrons.size()>0 || temp.V.muons.size()>0 ) return in;
    //    if (temp.additionalJets.size()>0) return in;
    //    if (std::Abs(deltaTheta) ????
    
    ok = true;
    return temp;
  }

  VHbbCandidate getHWmunCandidate(const VHbbCandidate & in, bool & ok , std::vector<unsigned int>& pos){
    ok = false;
    VHbbCandidate temp=in;
    // require a muon and no electrons
    if (temp.V.muons.size()!=1) return in ;
    if (temp.V.electrons.size()!=0) return in ;
    if (temp.V.mets.size()<1) return in ;
    //
    /*pT(W) > 150 GeV (pt(W) computed using lepton px,y and PF MET x and y components)
      pT(H) > 150 GeV
      best btag, CSV > 0.90
      second-best btag, CSV > 0.50
      Dphi(W,H) > 2.95
      no additional isolated leptons (pT > 15 GeV)
      
      same lepton definition as in the physics objects section 
      
      No additional ak5PFjets (pT > 30 GeV; |eta| < 2.4)
      MET>35. for the electron BDT analysis
      |cos(theta * )| (TBO)
      color flow pull angle (TBO)
      We don't cut on the transverse mass (for boosted objects cutting on the transverse mass introduces an inefficiency due to the angle between the MET and the lepton being close to 0.) 
    */
    
    temp.V.p4 = temp.V.muons[0].p4+temp.V.mets[0].p4;
    temp.V.firstLeptonOrig=pos[0];

    ok=true;
    return temp;
  }

  VHbbCandidate getHWenCandidate(const VHbbCandidate & in, bool & ok, std::vector<unsigned int>& pos){
    ok = false;
    VHbbCandidate temp=in;
    if (temp.V.electrons.size()!=1) return in ;
    if (temp.V.muons.size()!=0) return in ;
    if (temp.V.mets.size()<1) return in ;
    //
    /*pT(W) > 150 GeV (pt(W) computed using lepton px,y and PF MET x and y components)
      pT(H) > 150 GeV
      best btag, CSV > 0.90
      second-best btag, CSV > 0.50
      Dphi(W,H) > 2.95
      no additional isolated leptons (pT > 15 GeV)
      
      same lepton definition as in the physics objects section 
      
      No additional ak5PFjets (pT > 30 GeV; |eta| < 2.4)
      MET>35. for the electron BDT analysis
      |cos(theta * )| (TBO)
      color flow pull angle (TBO)
      We don't cut on the transverse mass (for boosted objects cutting on the transverse mass introduces an inefficiency due to the angle between the MET and the lepton being close to 0.) 
    */
    
    ok=true;
    temp.V.firstLeptonOrig=pos[0];
    temp.V.p4 = temp.V.electrons[0].p4+temp.V.mets[0].p4;
    return temp;
  }

 public: 
  bool verbose_;
  
 public:
  static float getDeltaTheta( const VHbbEvent::SimpleJet & j1, const VHbbEvent::SimpleJet & j2 ) {
    double deltaTheta = 1e10;
    TLorentzVector pi(0,0,0,0);
    TLorentzVector v_j1 = j1.chargedTracksFourMomentum;
    TLorentzVector v_j2 = j2.chargedTracksFourMomentum;
    
    if( v_j2.Mag() == 0 
	|| v_j1.Mag() == 0 )
      return deltaTheta = 1e10;
    
    //use j1 to calculate the pull vector
    TVector2 t = j1.tVector;
    
    if( t.Mod() == 0 )
      return deltaTheta = 1e10;
    
    Double_t dphi =  v_j2.Phi()- v_j1.Phi();
    if ( dphi > M_PI ) {
      dphi -= 2.0*M_PI;
    } else if ( dphi <= -M_PI ) {
      dphi += 2.0*M_PI;
    }
    Double_t deltaeta = v_j2.Rapidity() - v_j1.Rapidity();
    TVector2 BBdir( deltaeta, dphi );
    
    deltaTheta = t.DeltaPhi(BBdir);
    
    return deltaTheta;
    
  }
  
  float getHelicity(const VHbbEvent::SimpleJet& j, TVector3 boost) const {
    double hel = 1e10;
    TLorentzVector jet = j.p4;
    jet.Boost( -boost );
    hel = TMath::Cos( jet.Vect().Angle( boost ) );
    return hel;
  }


};



#endif






