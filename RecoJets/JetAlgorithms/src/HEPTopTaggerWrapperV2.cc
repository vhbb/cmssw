// 2011 Christopher Vermilion
//
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

#include "RecoJets/JetAlgorithms/interface/HEPTopTaggerWrapperV2.h"

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

#include "RecoJets/JetAlgorithms/interface/HEPTopTaggerV2.h"


FASTJET_BEGIN_NAMESPACE

//------------------------------------------------------------------------
// returns the tagged PseudoJet if successful, 0 otherwise
//  - jet   the PseudoJet to tag
PseudoJet HEPTopTaggerV2::result(const PseudoJet & jet) const{

  // make sure that there is a "regular" cluster sequence associated
  // with the jet. Note that we also check it is valid (to avoid a
  // more criptic error later on)
  if (!jet.has_valid_cluster_sequence()){
    throw Error("HEPTopTagger can only be applied on jets having an associated (and valid) ClusterSequence");
  }

  external::HEPTopTaggerV2 tagger(jet);

  // translate the massRatioWidth (which should be the half-width given in %) 
  // to values useful for the A-shape cuts
  double mw_over_mt = 80.4/172.3;
  double ratio_min = mw_over_mt * (100.-massRatioWidth_)/100.;
  double ratio_max = mw_over_mt * (100.+massRatioWidth_)/100.;
 
  
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
  tagger.set_mass_ratio_range(ratio_min, ratio_max);

  tagger.run_tagger();

  // Requires:
  //   - top mass window
  //   - mass ratio cuts
  //   - minimal candidate pT
  // If this is not intended: use loose top mass and ratio windows
  if (!tagger.is_tagged())
    return PseudoJet();
  
  // create the result and its structure
  const JetDefinition::Recombiner *rec
    = jet.associated_cluster_sequence()->jet_def().recombiner();

  const vector<PseudoJet>& subjets = tagger.top_subjets();
  assert(subjets.size() == 3);

  PseudoJet non_W = subjets[0];
  PseudoJet W1 = subjets[1];
  PseudoJet W2 = subjets[2];
  PseudoJet W = join(subjets[1], subjets[2], *rec);

  PseudoJet result = join<HEPTopTaggerV2Structure>( W1, W2, non_W, *rec);
  HEPTopTaggerV2Structure *s = (HEPTopTaggerV2Structure*) result.structure_non_const_ptr();

  s->_fj_mass  = jet.m();
  s->_fj_pt    = jet.perp();
  s->_fj_eta   = jet.eta();
  s->_fj_phi   = jet.phi();

  s->_top_mass = tagger.t().m();
  s->_pruned_mass = tagger.pruned_mass();
  s->_unfiltered_mass = tagger.unfiltered_mass();
  s->_fW = tagger.fW();
  s->_mass_ratio_passed = tagger.is_masscut_passed();

  // Removed selectors as all cuts are applied ion HTT
  return result;
}

FASTJET_END_NAMESPACE
