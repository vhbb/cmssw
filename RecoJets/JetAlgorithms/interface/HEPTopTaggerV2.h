//HEPTopTaggerV2: 
// Modes:
// 0 = EARLY_MASSRATIO_SORT_MASS (apply massratio, then sort by distance to true top mass)
// 1 = LATE_MASSRATIO_SORT_MASS (sort by distance to true top mass, then apply mass-ratio. Old HTT)
// 2 = EARLY_MASSRATIO_SORT_MODDJADE (apply massratio, then sort by modified d-jade)
// 3 = LATE_MASSRATIO_SORT_MODDJADE (sort by modified d-jade, then apply mass-ratio)
// 4 = TWO_STEP_FILTER (take three highest pT objects after unclustering and apply mass ratio)


#ifndef __HEPTOPTAGGERV2_HH__
#define __HEPTOPTAGGERV2_HH__

#include <math.h>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"

// Allow putting evertything into a separate namepsace
// Do not change next line, it's needed by the sed-code that makes the tagger CMSSW-compatible.
namespace external {


class HEPTopTaggerV2 {

public:

  enum Mode {EARLY_MASSRATIO_SORT_MASS, 
	     LATE_MASSRATIO_SORT_MASS, 
	     EARLY_MASSRATIO_SORT_MODDJADE,
	     LATE_MASSRATIO_SORT_MODDJADE,
	     TWO_STEP_FILTER};

  typedef fastjet::ClusterSequence ClusterSequence;
  typedef fastjet::JetAlgorithm JetAlgorithm;
  typedef fastjet::JetDefinition JetDefinition;
  typedef fastjet::PseudoJet PseudoJet;

  HEPTopTaggerV2();
  
  HEPTopTaggerV2(fastjet::PseudoJet jet);
  
  HEPTopTaggerV2(fastjet::PseudoJet jet,
	       double mtmass, double mwmass);

  //run tagger
  void run_tagger();
  
  //get information
  bool is_maybe_top() const {return _is_maybe_top;}
  bool is_masscut_passed() const {return _is_masscut_passed;}
  bool is_minptcut_passed() const {return _is_ptmincut_passed;}
  bool is_tagged() const {return (_is_masscut_passed && _is_ptmincut_passed);}
  const PseudoJet & top_candidate() const {return _top_candidate;}
  const std::vector<PseudoJet> & top_subjets() const {return _top_subjets;}
  const std::vector<PseudoJet> & top_hadrons() const {return _top_hadrons;}
  const std::vector<PseudoJet> & hardparts() const {return _top_parts;}
  unsigned parts_size() const {return _parts_size;}
  double delta_top() const {return _delta_top;}
  double djsum() const {return _djsum;}
  double pruned_mass() const {return _pruned_mass;}
  double unfiltered_mass() const {return _unfiltered_mass;}
  double fW();
  void get_setting() const;
  void get_info() const;
  const PseudoJet & t() const {return _top_candidate;}
  const PseudoJet & b() const {return _top_subjets[0];}
  const PseudoJet & W() const {return _W;}
  const PseudoJet & W1() const {return _top_subjets[1];}
  const PseudoJet & W2() const {return _top_subjets[2];}
  const PseudoJet & j1() const {return _top_subs[0];}
  const PseudoJet & j2() const {return _top_subs[1];}
  const PseudoJet & j3() const {return _top_subs[2];}
  const PseudoJet & fat() {return _fat;}
 
  //set parameters
  void set_max_subjet_mass(double x) {_max_subjet_mass = x;}
  void set_mass_drop_threshold(double x) {_mass_drop_threshold = x;}
  void set_top_range(double xmin, double xmax) {_mtmin = xmin; _mtmax = xmax;}
  void set_mass_ratio_range(double rmin, double rmax) {_rmin = rmin; _rmax = rmax;}
  void set_mass_ratio_cut(double m23cut, double m13cutmin,double m13cutmax) {_m23cut = m23cut; _m13cutmin = m13cutmin; _m13cutmax = m13cutmax;}
  void set_nfilt(unsigned nfilt) {_nfilt = nfilt;}
  void set_Rfilt(double Rfilt) {_Rfilt = Rfilt;}
  void set_filtering_jetalgorithm(JetAlgorithm jet_algorithm) {_jet_algorithm_filter = jet_algorithm;}
  void set_reclustering_jetalgorithm(JetAlgorithm jet_algorithm) {_jet_algorithm_recluster = jet_algorithm;}
  void set_pruner_cuts(double zcut, double rcut_factor) {_zcut = zcut; _rcut_factor = rcut_factor;}
  void set_mode(int mode) {_mode = Mode(mode);}
  void set_debug(bool debug) {_debug = debug;}
  void set_minpt_tag(double x) {_minpt_tag = x;}
  void set_minpt_subjet(double x) {_minpt_subjet = x;}
  
private:
  const PseudoJet* _jet;
  double _mtmass, _mwmass;
  double _mass_drop_threshold;
  double _max_subjet_mass;
  double _mtmin, _mtmax;
  double _rmin, _rmax;
  double _m23cut, _m13cutmin, _m13cutmax;
  size_t _nfilt;
  double _Rfilt;
  double _Rprun;
  JetAlgorithm _jet_algorithm_filter;
  JetAlgorithm _jet_algorithm_recluster;
  double _zcut;
  double _rcut_factor;
  Mode _mode;
  double _minpt_tag;
  double _minpt_subjet;
  bool _debug;
  PseudoJet _fat;
  
  bool _is_masscut_passed;
  bool _is_ptmincut_passed;
  bool _is_maybe_top;
  double _djsum;
  double _delta_top;
  double _pruned_mass;
  double _unfiltered_mass;
  double _fw;
  unsigned _parts_size;
  PseudoJet _top_candidate;
  PseudoJet _W;
  std::vector<PseudoJet> _top_subs;
  std::vector<PseudoJet> _top_subjets;
  std::vector<PseudoJet> _top_hadrons;
  std::vector<PseudoJet> _top_parts;
  static bool _first_time;

  //internal functions
  void FindHardSubst(const PseudoJet& jet, std::vector<fastjet::PseudoJet>& t_parts);
  std::vector<PseudoJet> Filtering(const std::vector <PseudoJet> & top_constits, const JetDefinition & filtering_def);
  void store_topsubjets(const std::vector<PseudoJet>& top_subs);
  bool check_mass_criteria(const std::vector<fastjet::PseudoJet> & top_subs) const;
  void print_banner();
  double perp(const PseudoJet & vec, const fastjet::PseudoJet & ref);
  double djademod (const fastjet::PseudoJet & subjet_i, const fastjet::PseudoJet & subjet_j, const fastjet::PseudoJet & ref);
};
//--------------------------------------------------------------------
// Do not change next line, it's needed by the sed-code that makes the tagger CMSSW-compatible.
};
#endif // __HEPTOPTAGGER_HH__
