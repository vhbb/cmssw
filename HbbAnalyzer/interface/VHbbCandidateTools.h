#ifndef VHBBCANDIDATETOOLS_H
#define VHBBCANDIDATETOOLS_H

#include "VHbbAnalysis/HbbAnalyzer/interface/VHbbCandidate.h"

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
    temp.V.fourMomentum = temp.V.muons[0].fourMomentum+temp.V.muons[1].fourMomentum;
   
    if (temp.V.muons[0].fourMomentum.Pt()<20 || temp.V.muons[1].fourMomentum.Pt()<20 ) return in;
    
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

    temp.V.fourMomentum = temp.V.electrons[0].fourMomentum+temp.V.electrons[1].fourMomentum;
   
    //
    // i need to ask VBTF and pt NEEDS ADJUSTING!!!!!
    //
    if (temp.V.electrons[0].fourMomentum.Pt()<20 ||temp.V.electrons[1].fourMomentum.Pt()<20  ) return in;
    if (temp.V.electrons[0].id95r < -100000 ||temp.V.electrons[1].id95r < -100000) return in;

    //    if (temp.V.fourMomentum.Pt()<150 ) return in;
    //    if (temp.H.fourMomentum.Pt()<150) return in;
    //    if (temp.H.firstJet().csv< 0.9) return in;
    //    if (temp.H.secondJet().csv<0.5) return in;
    //    if (deltaPhi(temp.V.fourMomentum.Phi(),temp.H.fourMomentum.Phi())<2.7) return in;
    //    if (temp.V.fourMomentum.M()<75 || temp.V.fourMomentum.M()>105) return in; 
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
    
    temp.V.fourMomentum = temp.V.mets[0].fourMomentum;
    if (verbose_) {
      std::cout <<" debug met "<< temp.V.mets[0].metSig << " " <<  temp.V.mets[0].sumEt<< std::endl;
    }   
    if (temp.V.mets[0].metSig<5) return in;
    if (temp.V.mets[0].sumEt<50) return in;
    //    if (temp.H.fourMomentum.Pt()<150)return in;
    //    if (temp.H.firstJet().csv< 0.9) return in;
    //    if (temp.H.secondJet().csv<0.5) return in;
    //    if (deltaPhi(temp.V.fourMomentum.Phi(),temp.H.fourMomentum.Phi())<2.95) return in;
    //    if (temp.V.electrons.size()>0 || temp.V.muons.size()>0 ) return in;
    //    if (temp.additionalJets.size()>0) return in;
    //    if (std::Abs(deltaTheta) ????
    
    ok = true;
    return temp;
  }

  VHbbCandidate getHWmunCandidate(const VHbbCandidate & in, bool & ok){
    ok = false;
    VHbbCandidate temp=in;
    if (temp.V.muons.size()!=1) return in ;
    if (temp.V.electrons.size()!=0) return in ;
    //
    // not clear to me, where is MET? it is not used????
    //
    return in;
  }

  VHbbCandidate getHWenCandidate(const VHbbCandidate & in, bool & ok){
    ok = false;
    VHbbCandidate temp=in;
    if (temp.V.electrons.size()!=1) return in ;
    if (temp.V.muons.size()!=0) return in ;
    //
    // not clear to me, where is MET? it is not used????
    //
    return in;
  }

 public: 
  bool verbose_;

};



#endif






