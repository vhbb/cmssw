//Adapted Version based on the work of Gregor.

#ifndef __MULTIR_TOPTAGGERV2_HH__
#define __MULTIR_TOPTAGGERV2_HH__

#include <vector>
#include <algorithm>  // for swap
#include <math.h>
#include "../interface/HEPTopTaggerV2.h"

using namespace std;
using namespace fastjet;

// Do not change next line, it's needed by the sed-code that makes the tagger CMSSW-compatible.
namespace external {


class MultiR_TopTagger {
public:
  MultiR_TopTagger(double max_fatjet_R,
		   double min_fatjet_R,
		   double step_R,
		   double multiR_threshold,
		   bool use_dR_max_triplet,
		   const fastjet::ClusterSequence & cs, 
		   const fastjet::PseudoJet & jet, 
		   double mtmass, double mwmass
		   );

  ~MultiR_TopTagger();

  //run tagger
  void run_tagger();

  // Return the candidate (and some properties) at R=R_min
  HEPTopTaggerV2 cand_Rmin(){return _HEPTopTaggerV2[_Rmin];}
  const int & Rmin_raw() const {return _Rmin;}
  const double Rmin() const {return _Rmin/10.;}
  const double & mass_Rmin() const {return _mass_Rmin;}
  const double & pt_Rmin() const {return _pt_Rmin;}
  
  // Return pT for usage in determining R_min expected filtered pT of
  // the initial fatjet, but with different filtering settings then
  // used for the HTT
  const double & pt_for_exp() const {return _pt_for_exp;}
  
  double R_min_exp() {return _r_min_exp_function(_pt_for_exp);}

  // Access to all candidates and number-of-small-fatjets
  HEPTopTaggerV2 HTTagger(int i)  {return _HEPTopTaggerV2[i];}
  const double n_small_fatjets(int i) {return _n_small_fatjets[i];}

  void set_mode(int mode) {_mode = mode;}

  void set_max_subjet_mass(double x) {_subjet_mass = x;}
  void set_mass_drop_threshold(double x) {_mass_drop_threshold = x;}
  void set_minpt_subjet(double x) {_minpt_subjet = x;}
  void set_minpt_tag(double x) {_minpt_tag = x;}
  void set_top_range(double top_range_min, double top_range_max) {_top_range[0] = top_range_min; _top_range[1] = top_range_max;}
  void set_f_W(double f_W) {_f_W = f_W;}
  void set_mass_ratio_cut(double mass_ratios_0, double mass_ratios_1, double mass_ratios_2) {_mass_ratios[0] = mass_ratios_0; _mass_ratios[1] = mass_ratios_1; _mass_ratios[2] = mass_ratios_2;}
  void set_nfilt(unsigned nfilt) {_n_filt = nfilt;}
  void set_Rfilt(double Rfilt) {_R_filt = Rfilt;}
  void set_nfilt_for_exp_Rmin(unsigned nfilt_exp) {_n_filt_exp = nfilt_exp;}
  void set_Rfilt_for_exp_Rmin(double Rfilt_exp) {_R_filt_exp = Rfilt_exp;}
  void set_debug(bool debug) {_debug = debug;}
  void set_r_min_exp_function(double (*f)(double)) {_r_min_exp_function = f;}
 

private:
  const ClusterSequence * _cs;
  const PseudoJet *       _jet;
  double _mtmass, _mwmass;
  double _mass_drop_threshold;
  double _subjet_mass;
  double _minpt_tag;
  double _minpt_subjet;
  double _pt_for_exp;
  map<int,external::HEPTopTaggerV2> _HEPTopTaggerV2;
  map<int,int> _n_small_fatjets;
  int _Rmin;
  int _mode;
  double _mass_Rmin, _pt_Rmin;
  double _mass_mean, _mass_width;
  double _top_range[2];
  unsigned _n_filt;
  double _R_filt;
  unsigned _n_filt_exp;
  double _R_filt_exp;
  double _f_W;
  double _mass_ratios[3];
  double _max_fatjet_R, _min_fatjet_R, _step_R, _multiR_threshold;
  bool _use_dR_max_triplet;
  bool _debug;
  double (*_r_min_exp_function)(double);

  void UnclusterFatjets(const vector<fastjet::PseudoJet> & big_fatjets, vector<fastjet::PseudoJet> & small_fatjets, const ClusterSequence & cs, const double small_radius);

};
//--------------------------------------------------------------------
// Do not change next line, it's needed by the sed-code that makes the tagger CMSSW-compatible.
};

#endif // __MULTIR_TOPTAGGER_HH__

