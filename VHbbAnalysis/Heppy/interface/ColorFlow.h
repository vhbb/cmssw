// author : Pierluigi Bortignon
// email : pierluigi.bortignon@gmail.com
// date : 16.10.2015 - Fermilab
// version : 0

#ifndef PhysicsTools_Heppy_ColorFlow_h
#define PhysicsTools_Heppy_ColorFlow_h

#include <vector>
#include <iostream>
#include <cmath>
#include <TVector2.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <boost/python.hpp>

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"


namespace heppy{
class ColorFlow {
    
 public:
 // This class cluster a set of variables that could have sensitivity to color flow connection between two jets
 // The variable is inspired from this paper: http://arxiv.org/abs/1001.5027
 //
 // The class defines a 2D vector (a pull vector) per each jet. The pull vector represents the direction of the radiation in the rapiity-phi plane.
 // The interface suggests to use the phi of the pull vector (in the rapidity-phi plane) and its magnitude as sensitive variables to color connection.


// ColorFlow(const std::vector<pat::Jet> &jets);
 ColorFlow(boost::python::object *jets);
 
 float getPullVectorPhi() { return r_mag_; };
 float getPullVectorMag() { return r_phi_; };

private:

  reco::PFCandidateCollection * pfCands;

  TVector2 calculate_pull_vector_(const std::vector<pat::Jet> &jets) ;

  TVector2 pull_;
  TVector2 null_;
  TVector2 ci_;
  TVector2 r_;
  TVector2 t_vector_;
  TLorentzVector pi_;
  TLorentzVector J_;
  float r_mag_;
  float r_phi_;
  float patJetpfcPt_;
  float max_pfCand_pt;
  float pfCand_pt_;
  unsigned int nOfConst_;

//  const std::vector<pat::Jet> jets_;
  
};
}
#endif   
 
