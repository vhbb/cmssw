//----------------------------------------------------------------------
//  This file is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  This file is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  The GNU General Public License is available at
//  http://www.gnu.org/licenses/gpl.html or you can write to the Free Software
//  Foundation, Inc.:
//      59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//----------------------------------------------------------------------

#include "RecoJets/JetAlgorithms/interface/MultiRHEPTopTaggerWrapper.h"

#include <fastjet/Error.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"

#include <math.h>
#include <limits>
#include <cassert>
using namespace std;

#include "RecoJets/JetAlgorithms/interface/HEPTopTagger.h"
#include "RecoJets/JetAlgorithms/interface/MultiRHEPTopTagger.h"

FASTJET_BEGIN_NAMESPACE


// Expected R_min for tops (as function of filtered initial fatjet pT in GeV (using CA, R=0.2, n=10)
// From ttbar sample, matched to hadronically decaying top with delta R < 0.8 and true top pT > 200
// Cuts are: fW < 0.175 and  m_top = 120..170
// Input objects are packed pfCandidates (wo/ filtering)
// IMPORTANT: this might need to be changed when using CHS or very different other cuts
double R_min_expected_function(double x){
  if (x<300)
    return 1.17 + 1.91e-03*x - 6.45e-06*x*x;
  else if (x<500)
    return 1.89 - 2.89e-03*x + 1.55e-06*x*x;
  else
    return 1.86 - 2.78e-03*x + 1.44e-06*x*x;
}



//------------------------------------------------------------------------
// returns the tagged PseudoJet if successful, 0 otherwise
//  - jet   the PseudoJet to tag
PseudoJet MultiRHEPTopTagger::result(const PseudoJet & jet) const{

  // make sure that there is a "regular" cluster sequence associated
  // with the jet. Note that we also check it is valid (to avoid a
  // more criptic error later on)
  if (!jet.has_valid_cluster_sequence()){
    throw Error("HEPTopTagger can only be applied on jets having an associated (and valid) ClusterSequence");
  }

  double m_W = 80.4;
  double m_top = 172.3;

  external::MultiR_TopTagger tagger(R_max_,   // R_max
				    R_min_,   // R_min
				    0.1,   // R_step (using a stepsize smaller than 0.1 would currently not work)
				    0.2,   // Mass-drop threshold
				    false, // use_dR_max_triplet
				    *(jet.associated_cluster_sequence()),
				    jet,
				    m_top,
				    m_W);
   
  // Unclustering, Filtering & Subjet Settings
  tagger.set_max_subjet_mass(subjetMass_);
  tagger.set_mass_drop_threshold(muCut_); 
  tagger.set_Rfilt(filtR_);
  tagger.set_nfilt(filtN_);
  tagger.set_minpt_subjet(minSubjetPt_);

  // How to select among candidates
  tagger.set_mode(mode_);
  
  // Requirements to accept a candidate
  tagger.set_minpt_tag(minCandPt_); 
  tagger.set_top_range(minCandMass_, maxCandMass_); 
  tagger.set_mass_ratio_cut(minM23Cut_, minM13Cut_, maxM13Cut_);
  tagger.set_f_W(massRatioWidth_/100.);

  // Set function to calculate R_min_expected
  tagger.set_r_min_exp_function(R_min_expected_function);

  tagger.run_tagger();
  
  // MultiR tagging requirements:
  //   - found a valid R_min
  //   - HTT top mass window
  //   - HTT mass ratio cuts
  //   - HTT minimal candidate pT
  // If this is not intended: use loose top mass and ratio windows
  if ( (tagger.Rmin_raw()==0) || (! tagger.cand_Rmin().is_tagged()))
      return PseudoJet();
  
  // create the result and its structure
  const JetDefinition::Recombiner *rec
    = jet.associated_cluster_sequence()->jet_def().recombiner();

  const vector<PseudoJet>& subjets = tagger.cand_Rmin().top_subjets();
  assert(subjets.size() == 3);

  PseudoJet non_W = subjets[0];
  PseudoJet W1 = subjets[1];
  PseudoJet W2 = subjets[2];
  PseudoJet W = join(subjets[1], subjets[2], *rec);

  PseudoJet result = join<MultiRHEPTopTaggerStructure>( W1, W2, non_W, *rec);
  MultiRHEPTopTaggerStructure *s = (MultiRHEPTopTaggerStructure*) result.structure_non_const_ptr();

  s->_fj_mass  = jet.m();
  s->_fj_pt    = jet.perp();
  s->_fj_eta   = jet.eta();
  s->_fj_phi   = jet.phi();

  s->_top_mass = tagger.cand_Rmin().t().m();
  s->_pruned_mass = tagger.cand_Rmin().pruned_mass();
  s->_unfiltered_mass = tagger.cand_Rmin().unfiltered_mass();
  s->_fW = tagger.cand_Rmin().fW();
  s->_mass_ratio_passed = tagger.cand_Rmin().is_masscut_passed();
  s->_Rmin = tagger.Rmin();
  s->_ptFiltForRminExp = tagger.pt_for_exp(); // CA, R=0.2, n=10 is the current default in the tagger
  s->_RminExpected = tagger.R_min_exp();
  
  return result;
}

FASTJET_END_NAMESPACE
