#ifndef VHBBCANDIDATETOOLS_H
#define VHBBCANDIDATETOOLS_H

#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"

#include <iostream>

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

  VHbbCandidate getHZmumuCandidate(const VHbbCandidate & in, bool & ok){
    if (verbose_){
      std::cout <<" getHZmumuCandidate input mu "<<in.V.muons.size()<<" e "<<in.V.electrons.size()<<std::endl;
    }
    ok = false;
    VHbbCandidate temp=in;
    if (temp.V.muons.size()!=2) return in ;
    if (temp.V.electrons.size()!=0) return in ;
    temp.V.p4 = temp.V.muons[0].p4+temp.V.muons[1].p4;
   
    if (temp.V.muons[0].p4.Pt()<20 || temp.V.muons[1].p4.Pt()<20 ) return in;
    
    //    if (temp.V.Pt()<150 ) return in;
    //    if (temp.H.Pt()<150) return in;
    //    if (temp.H.firstJet().csv< 0.9) return in;
    //    if (temp.H.secondJet().csv<0.5) return in;
    //    if (deltaPhi(temp.V.Phi(),temp.H.Phi())<2.7) return in;
    //    if (temp.V.FourMomentum.Mass()<75 || temp.V.FourMomentum.Mass()>105) return in; 
    //    if (temp.additionalJets.size()>0) return in;
    //    if (std::Abs(deltaTheta) ????
    ok = true;
    return temp;
  }  
  VHbbCandidate getHZeeCandidate(const VHbbCandidate & in, bool & ok){
    if (verbose_){
      std::cout <<" getHZeeCandidate input mu "<<in.V.muons.size()<<" e "<<in.V.electrons.size()<<std::endl;
    }

    ok = false;
    VHbbCandidate temp=in;
    if (temp.V.electrons.size()!=2) return in ;
    if (temp.V.muons.size()!=0) return in ;

    temp.V.p4 = temp.V.electrons[0].p4+temp.V.electrons[1].p4;
   
    //
    // i need to ask VBTF and pt NEEDS ADJUSTING!!!!!
    //
    if (temp.V.electrons[0].p4.Pt()<20 ||temp.V.electrons[1].p4.Pt()<20  ) return in;
    if (temp.V.electrons[0].id95r < -100000 ||temp.V.electrons[1].id95r < -100000) return in;

    //    if (temp.V.p4.Pt()<150 ) return in;
    //    if (temp.H.p4.Pt()<150) return in;
    //    if (temp.H.firstJet().csv< 0.9) return in;
    //    if (temp.H.secondJet().csv<0.5) return in;
    //    if (deltaPhi(temp.V.p4.Phi(),temp.H.p4.Phi())<2.7) return in;
    //    if (temp.V.p4.M()<75 || temp.V.p4.M()>105) return in; 
    //    if (temp.additionalJets.size()>0) return in;
    //    if (std::Abs(deltaTheta) ????
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
    if (temp.V.muons.size()!=0) return in ;
    if (temp.V.electrons.size()!=0) return in ;
    
    temp.V.p4 = temp.V.mets[0].p4;
    if (verbose_) {
      std::cout <<" debug met "<< temp.V.mets[0].metSig << " " <<  temp.V.mets[0].sumEt<< std::endl;
    }   
    if (temp.V.mets[0].metSig<5) return in;
    if (temp.V.mets[0].sumEt<50) return in;
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

  VHbbCandidate getHWmunCandidate(const VHbbCandidate & in, bool & ok){
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

    ok=true;
    return temp;
  }

  VHbbCandidate getHWenCandidate(const VHbbCandidate & in, bool & ok){
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

    temp.V.p4 = temp.V.electrons[0].p4+temp.V.mets[0].p4;
    return temp;
  }

 public: 
  bool verbose_;
  
 public:
  float getDeltaTheta( const VHbbEvent::SimpleJet & j1, const VHbbEvent::SimpleJet & j2 ) const {
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






