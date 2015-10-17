// author : Pierluigi Bortignon
// email : pierluigi.bortignon@gmail.com
// date : 16.10.2015 - Fermilab
// version : 0

//#include "PhysicsTools/Heppy/interface/ReclusterJets.h"
#include "VHbbAnalysis/Heppy/interface/ColorFlow.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

using namespace std;

namespace heppy{

  // ColorFlow::ColorFlow(const std::vector<pat::Jet> &jets):
  // jets_(jets)

  ColorFlow::ColorFlow(boost::python::object *jets)
    {
      // initialisation - for the moment here, then I'll move it out.
      max_pfCand_pt = 1;

      // loop over the jets
      for( unsigned int i = 0; i<jets->attr("size"); i++ ){
        // get the number of pfCandidates of a jet
        unsigned int nOfCands = boost::python::extract<unsigned int>(jets->attr("at")(i).attr("numberOfSourceCandidatePtrs"));
        // loop over the PF candidate of a jet
        for ( unsigned int cand_j=0; cand_j < nOfCands; cand_j++){
          // get a pf candidate
          reco::CandidatePtr pfCand = boost::python::extract<reco::CandidatePtr> (jets->attr("at")(i).attr("sourceCandidatePtr")(cand_j));
          if (!pfCand) continue; // check if it is available (?)
          if (pfCand->charge()==0 or pfCand->pt() < max_pfCand_pt) continue; // calculating a tvector using only charged tracks. This is how it was in the Run1 analysis.

          pi_.SetPtEtaPhiE( pfCand->pt(), pfCand->eta(), pfCand->phi(), pfCand->energy() );
          J_+=pi_; // building the jet using only yhe chared pf candidates four momentum
          nOfConst_+=1; // count the number of charge constituents

        } // candidate loop
        
        //if there are no more than 2 consituents there is no way to calculate the pull Angle - maybe in the future this can be revised.
        if (nOfConst_ < 2) t_vector_.Set(0.,0.);
        else{ // else calculate the t_vector

          for ( unsigned int cand_j=0; cand_j < nOfCands; cand_j++){
          // get a pf candidate
          reco::CandidatePtr pfCand = boost::python::extract<reco::CandidatePtr> (jets->attr("at")(i).attr("sourceCandidatePtr")(cand_j));
          //reco::CandidatePtr pfCand = jets->at(i).sourceCandidatePtr(cand_j);
          if (!pfCand) continue; // check if it is available (?)
          if (pfCand->charge()==0 or pfCand->pt() < max_pfCand_pt) continue; // calculating a tvector using only charged tracks. This is how it was in the Run1 analysis.

          pfCand_pt_ = pfCand->pt();
          pi_.SetPtEtaPhiE(pfCand->pt(),pfCand->eta(),pfCand->phi(),pfCand->energy());
          r_.Set( pi_.Rapidity() - J_.Rapidity(), deltaPhi( pfCand->phi(), J_.Phi() ) );
          r_mag_ = r_.Mod();
          t_vector_ += ( pfCand_pt_ / J_.Pt() ) * r_mag_ * r_;

          } // candidate loop

        } // else

        // now fill the jet information
        // here I have to ind a way to add the property to the jets from the c++. which kind of objects are these?

      } // jet loop

    } // constructor


} // namespace heppy

