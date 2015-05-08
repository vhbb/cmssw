#include "../interface/MultiRHEPTopTagger.h"

// Do not change next line, it's needed by the sed-code that makes the tagger CMSSW-compatible.
namespace external {


//uncluster a fat jet to subjets of given cone size
void MultiR_TopTagger::UnclusterFatjets(const vector<fastjet::PseudoJet> & big_fatjets, 
					vector<fastjet::PseudoJet> & small_fatjets, 
					const ClusterSequence & cs, 
					const double small_radius) {
  for (unsigned i=0; i < big_fatjets.size(); i++) {
    PseudoJet this_jet = big_fatjets[i];
    PseudoJet parent1(0, 0, 0, 0), parent2(0, 0, 0, 0);
    bool test = cs.has_parents(this_jet, parent1, parent2);
    double dR = sqrt(parent1.squared_distance(parent2));

    if (!test || dR<small_radius) {
      small_fatjets.push_back(this_jet);
    } else {
      vector<fastjet::PseudoJet> parents;
      parents.push_back(parent1);
      parents.push_back(parent2);
      UnclusterFatjets(parents, small_fatjets, cs, small_radius);
    }
  }
}

//////MultiR_TopTagger/////////////////////////
//Start with a fatjet clustered with radius max_fatjet_R
//Then run the HEPTopTagger for fatjets with smaller radii
//Quantities stored as map<double,...> where the key is the fatjet radius

MultiR_TopTagger::MultiR_TopTagger(double max_fatjet_R,
				   double min_fatjet_R,
				   double step_R,
				   double multiR_threshold,
				   bool use_dR_max_triplet,
				   const fastjet::ClusterSequence & cs, 
				   const fastjet::PseudoJet & jet, 
				   double mtmass, double mwmass
				   ) : _cs(&cs),  _jet(&jet),
				       _mtmass(mtmass),	_mwmass(mwmass), _mass_drop_threshold(0.8), _minpt_tag(200.), _minpt_subjet(0.), _mode(1), _n_filt_exp(10), _R_filt_exp(0.2),
				       _max_fatjet_R(max_fatjet_R), _min_fatjet_R(min_fatjet_R), _step_R(step_R), _multiR_threshold(multiR_threshold), 
				       _use_dR_max_triplet(use_dR_max_triplet), _debug(false)
{}

void MultiR_TopTagger::run_tagger() {
  if (_debug) {
    cout << "============================="  << endl 
	 << "new MultiR" << endl;
  }

  // Before doing the MultiR procedure:
  // determine the filtered fatjet pT that can be used for 
  // fitting R_min_expected 
  fastjet::Filter filter(JetDefinition(fastjet::cambridge_algorithm, _R_filt_exp), SelectorNHardest(_n_filt_exp));
  _pt_for_exp = filter.result(*_jet).perp();

   
  // Do MultiR procedure  
  vector<fastjet::PseudoJet> big_fatjets;
  vector<fastjet::PseudoJet> small_fatjets;
  
  big_fatjets.push_back(* _jet);
  _Rmin = 0;
  _mass_Rmin = 0.;
  _pt_Rmin = 0.;

  int maxR = int(_max_fatjet_R * 10);
  int minR = int(_min_fatjet_R * 10);
  int stepR = int(_step_R * 10);
    
  for (int R = maxR; R >= minR; R -= stepR) {
    
    // TODO: check!
    // Deactivated this piece for now as dR_max_triplet is not used.
    // Problem in lines below?
    //  
    //float dR_max_triplet = 9999;
    //if (_use_dR_max_triplet) {
    //  small_fatjets = big_fatjets;
    //  dR_max_triplet = R / 10.;   
    //} else {

    UnclusterFatjets(big_fatjets, small_fatjets, *_cs, R / 10.);
    //}
    
    if (_debug) {cout << "R = " << R << " -> n_small_fatjets = " << small_fatjets.size();}
    
    _n_small_fatjets[R] = small_fatjets.size();

    // We are sorting by pt - so start with a negative dummy
    double dummy = -99999;

    for (unsigned i = 0; i < small_fatjets.size(); i++) {
      external::HEPTopTaggerV2 htt(small_fatjets[i], _mtmass, _mwmass);
      htt.set_top_range(_top_range[0], _top_range[1]);
      htt.set_mass_ratio_cut(_mass_ratios[0], _mass_ratios[1], _mass_ratios[2]);
      htt.set_max_subjet_mass(_subjet_mass);
      htt.set_minpt_subjet(_minpt_subjet);
      htt.set_minpt_tag(_minpt_tag);
      htt.set_mass_drop_threshold(_mass_drop_threshold);
      htt.set_nfilt(_n_filt);
      htt.set_Rfilt(_R_filt);
      htt.set_mass_ratio_range((1.-_f_W)*_mwmass/_mtmass, (1.+_f_W)*_mwmass/_mtmass); 
      htt.set_mode(_mode); 
      
      htt.run_tagger();
     
      if (htt.top_candidate().perp() > dummy) {
	dummy = htt.top_candidate().perp();
	_HEPTopTaggerV2[R] = htt;
      }
    } //End of loop over small_fatjets
    
    // Only check if we have not found Rmin yet
    if (_Rmin == 0 && R < maxR) {                 
      // If the new mass is OUTSIDE the window ..
      if (_HEPTopTaggerV2[R].top_candidate().m() < (1-_multiR_threshold)*_HEPTopTaggerV2[maxR].top_candidate().m())
	// .. set _Rmin to the previous mass 
	_Rmin = R + stepR;
    }
    
    big_fatjets = small_fatjets;
    small_fatjets.clear();
  }//End of loop over R

  // if we did not find Rmin in the loop, pick the last value
  if (_Rmin == 0 && _HEPTopTaggerV2[maxR].top_candidate().m() > 0)
    _Rmin = minR;

  //for the case that there is no tag at all (< 3 hard substructures)
  if (_Rmin == 0 && _HEPTopTaggerV2[maxR].top_candidate().m() == 0)
    _Rmin = maxR;

  _mass_Rmin = _HEPTopTaggerV2[_Rmin].top_candidate().m();
  _pt_Rmin = _HEPTopTaggerV2[_Rmin].top_candidate().perp();
 
  if (_debug) {cout << "MultiR done" << endl;}
}

MultiR_TopTagger::~MultiR_TopTagger(){}

// Do not change next line, it's needed by the sed-code that makes the tagger CMSSW-compatible.
};
