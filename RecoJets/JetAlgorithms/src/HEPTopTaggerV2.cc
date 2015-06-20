
#include "../interface/HEPTopTaggerV2.h"

// Do not change next line, it's needed by the sed-code that makes the tagger CMSSW-compatible.
namespace external {

bool HEPTopTaggerV2::_first_time = true;

void HEPTopTaggerV2::print_banner() {
  if (!_first_time) {return;}
  _first_time = false;

  std::cout << "#--------------------------------------------------------------------------\n";
  std::cout << "#                   HEPTopTaggerV2 - under construction                      \n";
  std::cout << "#                                                                          \n";
  std::cout << "# Please cite JHEP 1010 (2010) 078 [arXiv:1006.2833 [hep-ph]]              \n";
  std::cout << "# and Phys.Rev. D89 (2014) 074047 [arXiv:1312.1504 [hep-ph]]               \n";
  std::cout << "#--------------------------------------------------------------------------\n";
  get_setting();
}

double HEPTopTaggerV2::perp(const PseudoJet & vec, const fastjet::PseudoJet & ref) {
  double ref_ref = ref.px() * ref.px() + ref.py() * ref.py() + ref.pz() * ref.pz();
  double vec_ref = vec.px() * ref.px() + vec.py() * ref.py() + vec.pz() * ref.pz();
  double per_per = vec.px() * vec.px() + vec.py() * vec.py() + vec.pz() * vec.pz();
  if (ref_ref > 0.) 
    per_per -= vec_ref * vec_ref / ref_ref;
  if (per_per < 0.) 
    per_per = 0.;
  return sqrt(per_per);
}

double HEPTopTaggerV2::djademod (const fastjet::PseudoJet& subjet_i, const fastjet::PseudoJet& subjet_j, const fastjet::PseudoJet& ref) {
  double dj = -1.0;
  double delta_phi = subjet_i.delta_phi_to(subjet_j);
  double delta_eta = subjet_i.eta() - subjet_j.eta();
  double delta_R = sqrt(delta_eta * delta_eta + delta_phi * delta_phi);	
  dj = perp(subjet_i, ref) * perp(subjet_j, ref) * pow(delta_R, 4.);
  return dj;
}

double HEPTopTaggerV2::fW() {
  // Minimal:
  // |(m_ij / m_123) / (m_w/ m_t) - 1|

  double m12 = (_top_subs[0] + _top_subs[1]).m();
  double m13 = (_top_subs[0] + _top_subs[2]).m();
  double m23 = (_top_subs[1] + _top_subs[2]).m();
  double m123 = (_top_subs[0] + _top_subs[1] + _top_subs[2]).m();

  double fw12 = fabs( (m12/m123) / (_mwmass/_mtmass) - 1);
  double fw13 = fabs( (m13/m123) / (_mwmass/_mtmass) - 1);
  double fw23 = fabs( (m23/m123) / (_mwmass/_mtmass) - 1);
  
  return std::min(fw12, std::min(fw13, fw23));  
}

//Find hard substructures
void HEPTopTaggerV2::FindHardSubst(const PseudoJet & this_jet, std::vector<fastjet::PseudoJet> & t_parts) {
  PseudoJet parent1(0, 0, 0, 0), parent2(0, 0, 0, 0);
  if (this_jet.m() < _max_subjet_mass || !this_jet.validated_cs()->has_parents(this_jet, parent1, parent2)) {
    t_parts.push_back(this_jet);
  } else {
    if (parent1.m() < parent2.m()) 
      std::swap(parent1, parent2);   
    FindHardSubst(parent1, t_parts);
    if (parent1.m() < _mass_drop_threshold * this_jet.m())
      FindHardSubst(parent2, t_parts);   
  }
}

//store subjets as vector<PseudoJet> with [0]->b [1]->W-jet 1 [2]->W-jet 2
void HEPTopTaggerV2::store_topsubjets(const std::vector<PseudoJet>& top_subs) {
  _top_subjets.resize(0);
  double m12 = (top_subs[0] + top_subs[1]).m();
  double m13 = (top_subs[0] + top_subs[2]).m();
  double m23 = (top_subs[1] + top_subs[2]).m();
  double dm12 = fabs(m12 - _mwmass);
  double dm13 = fabs(m13 - _mwmass);
  double dm23 = fabs(m23 - _mwmass);
  
  if (dm23 <= dm12 && dm23 <= dm13) {
    _top_subjets.push_back(top_subs[0]); 
    _top_subjets.push_back(top_subs[1]); 
    _top_subjets.push_back(top_subs[2]);	
  } else if (dm13 <= dm12 && dm13 < dm23) {
    _top_subjets.push_back(top_subs[1]);
    _top_subjets.push_back(top_subs[0]);
    _top_subjets.push_back(top_subs[2]);
  } else if (dm12 < dm23 && dm12 < dm13) {
    _top_subjets.push_back(top_subs[2]);
    _top_subjets.push_back(top_subs[0]);
    _top_subjets.push_back(top_subs[1]);
  }
  _W = _top_subjets[1] + _top_subjets[2];
  return;
}

//check mass plane cuts
bool HEPTopTaggerV2::check_mass_criteria(const std::vector<PseudoJet> & top_subs) const {
  bool is_passed = false;
  double m12 = (top_subs[0] + top_subs[1]).m();
  double m13 = (top_subs[0] + top_subs[2]).m();
  double m23 = (top_subs[1] + top_subs[2]).m();
  double m123 = (top_subs[0] + top_subs[1] + top_subs[2]).m();
  if (
      (atan(m13/m12) > _m13cutmin && _m13cutmax > atan(m13/m12)
       && (m23/m123 > _rmin && _rmax > m23/m123))
      ||
      (((m23/m123) * (m23/m123) < 1 - _rmin * _rmin* (1 + (m13/m12) * (m13/m12)))
       &&
       ((m23/m123) * (m23/m123) > 1 - _rmax * _rmax * (1 + (m13/m12) * (m13/m12)))
       && 
       (m23/m123 > _m23cut))
      ||
      (((m23/m123) * (m23/m123) < 1 - _rmin * _rmin * (1 + (m12/m13) * (m12/m13)))
       &&
       ((m23/m123) * (m23/m123) > 1 - _rmax * _rmax * (1 + (m12/m13) * (m12/m13)))
       && 
       (m23/m123 > _m23cut))
      ) { 
    is_passed = true;
  }
  return is_passed;
}

HEPTopTaggerV2::HEPTopTaggerV2() {}

HEPTopTaggerV2::HEPTopTaggerV2(fastjet::PseudoJet jet) : 
  _jet(&jet), _mtmass(172.3), _mwmass(80.4), 
  _mass_drop_threshold(0.8), _max_subjet_mass(30.), 
  _mtmin(150.), _mtmax(200.), _rmin(0.85*80.4/172.3), _rmax(1.15*80.4/172.3), 
  _m23cut(0.35), _m13cutmin(0.2), _m13cutmax(1.3), 
  _nfilt(5), _Rfilt(0.3), _Rprun(1.5), _jet_algorithm_filter(fastjet::cambridge_algorithm), _jet_algorithm_recluster(fastjet::cambridge_algorithm), _zcut(0.1),
  _rcut_factor(0.5), _mode(Mode(0)), _minpt_tag(200.), _minpt_subjet(0.), _debug(false), _fat(jet)
{}

HEPTopTaggerV2::HEPTopTaggerV2(fastjet::PseudoJet jet, 
			   double mtmass, double mwmass
			   ) : 
  _jet(&jet), _mtmass(mtmass), _mwmass(mwmass), 
  _mass_drop_threshold(0.8), _max_subjet_mass(30.), 
  _mtmin(150.), _mtmax(200.), _rmin(0.85*mwmass/mtmass), _rmax(1.15*mwmass/mtmass), 
  _m23cut(0.35), _m13cutmin(0.2), _m13cutmax(1.3), 
  _nfilt(5), _Rfilt(0.3), _Rprun(1.5), _jet_algorithm_filter(fastjet::cambridge_algorithm), _jet_algorithm_recluster(fastjet::cambridge_algorithm), _zcut(0.1),
  _rcut_factor(0.5), _mode(Mode(0)), _minpt_tag(200.), _minpt_subjet(0.), _debug(false), _fat(jet)
{}

void HEPTopTaggerV2::run_tagger() {
  print_banner();

  if ((_mode != EARLY_MASSRATIO_SORT_MASS) 
      && (_mode != LATE_MASSRATIO_SORT_MASS) 
      && (_mode != EARLY_MASSRATIO_SORT_MODDJADE)
      && (_mode != LATE_MASSRATIO_SORT_MODDJADE)
      && (_mode != TWO_STEP_FILTER) ) {
    std::cout << "ERROR: UNKNOWN MODE" << std::endl;
    return;
  }
  
  //initialization
  _djsum = 0.;
  _delta_top = 1000000000000.0;
  _pruned_mass = 0.;
  _unfiltered_mass = 0.;
  _top_candidate.reset(0., 0., 0., 0.);
  _parts_size = 0;
  _is_maybe_top = _is_masscut_passed = _is_ptmincut_passed = false;
  _top_subs.clear();
  _top_subjets.clear();
  _top_hadrons.clear();
  _top_parts.clear();
  
  //find hard substructures
  FindHardSubst(*_jet, _top_parts);
  _parts_size = _top_parts.size();
  
  if (_top_parts.size() < 3) { 
    if (_debug) {std::cout << "< 3 hard substructures " << std::endl;}
    return; //such events are not interesting   
  }
  
  // Sort subjets-after-unclustering by pT.
  // Necessary so that two-step-filtering can use the leading-three.
  _top_parts=sorted_by_pt(_top_parts);

  // loop over triples
  _top_parts = sorted_by_pt(_top_parts);
  for (unsigned rr = 0; rr < _top_parts.size(); rr++) {
    for (unsigned ll = rr + 1; ll < _top_parts.size(); ll++) {
      for (unsigned kk = ll + 1; kk < _top_parts.size(); kk++) {
	
	// two-step filtering 
	// This means that we only look at the triplet formed by the
	// three leading-in-pT subjets-after-unclustering.
	if((_mode==TWO_STEP_FILTER) && rr>0)
	  continue;
	if((_mode==TWO_STEP_FILTER) && ll>1)
	  continue;
	if((_mode==TWO_STEP_FILTER) && kk>2)
	  continue;

      	//pick triple
	PseudoJet triple = join(_top_parts[rr], _top_parts[ll], _top_parts[kk]);
	
	//filtering 
	double filt_top_R 
	  = std::min(_Rfilt, 0.5*sqrt(std::min(_top_parts[kk].squared_distance(_top_parts[ll]), 
				     std::min(_top_parts[rr].squared_distance(_top_parts[ll]), 
					 _top_parts[kk].squared_distance(_top_parts[rr])))));
	JetDefinition filtering_def(_jet_algorithm_filter, filt_top_R);
	fastjet::Filter filter(filtering_def, fastjet::SelectorNHardest(_nfilt) * fastjet::SelectorPtMin(_minpt_subjet));
	PseudoJet topcandidate = filter(triple);


	//mass window cut
  	if (topcandidate.m() < _mtmin || _mtmax < topcandidate.m()) continue;

	// Sanity cut: can't recluster less than 3 objects into three subjets
	if (topcandidate.pieces().size() < 3)
	  continue;
       
	// Recluster to 3 subjets and apply mass plane cuts
	// Use a self-deleting CS-pointer. Taken from CMSSW version of HTT
	// Initial CMSSW edit by CKV, suggested by G. P. Salam)	
	JetDefinition reclustering(_jet_algorithm_recluster, 3.14);
	ClusterSequence *  cs_top_sub = new ClusterSequence(topcandidate.pieces(), reclustering);
        std::vector <PseudoJet> top_subs = sorted_by_pt(cs_top_sub->exclusive_jets(3));         
	cs_top_sub->delete_self_when_unused();

	// Require the third subjet to be above the pT threshold
	if (top_subs[2].perp() < _minpt_subjet)
	  continue;

	// Modes with early 2d-massplane cuts
	if (_mode == EARLY_MASSRATIO_SORT_MASS      && !check_mass_criteria(top_subs)) {continue;}
	if (_mode == EARLY_MASSRATIO_SORT_MODDJADE  && !check_mass_criteria(top_subs)) {continue;}

	//is this candidate better than the other? -> update
	double deltatop = fabs(topcandidate.m() - _mtmass);
	double djsum = djademod(top_subs[0], top_subs[1], topcandidate) 
			       + djademod(top_subs[0], top_subs[2], topcandidate)
			       + djademod(top_subs[1], top_subs[2], topcandidate);
	bool better = false;

	// Modes 0 and 1 sort by top mass
	if ( (_mode == EARLY_MASSRATIO_SORT_MASS) 
	     || (_mode == LATE_MASSRATIO_SORT_MASS)) {
	  if (deltatop < _delta_top) 
	    better = true;
	}
	// Modes 2 and 3 sort by modified jade distance
	else if ( (_mode == EARLY_MASSRATIO_SORT_MODDJADE) 
		  || (_mode == LATE_MASSRATIO_SORT_MODDJADE)) {
	  if (djsum > _djsum) 
	    better = true;
	}
	// Mode 4 is the two-step filtering. No sorting necessary as
	// we just look at the triplet of highest pT objects after
	// unclustering
	else if (_mode == TWO_STEP_FILTER) {
	  better = true;
	} 
	else {
	  std::cout << "ERROR: UNKNOWN MODE (IN DISTANCE MEASURE SELECTION)" << std::endl;
	  return;
	}

	if (better) {
	  _djsum = djsum;
	  _delta_top = deltatop; 
	  _is_maybe_top = true;
	  _top_candidate = topcandidate;
	  _top_subs = top_subs;
	  store_topsubjets(top_subs);
	  _top_hadrons = topcandidate.constituents();
	  //pruning
	  JetDefinition jet_def_prune(fastjet::cambridge_algorithm, _Rprun);
	  fastjet::Pruner pruner(jet_def_prune, _zcut, _rcut_factor);
	  PseudoJet prunedjet = pruner(triple);
	  _pruned_mass = prunedjet.m();
	  _unfiltered_mass = triple.m();
	  
	  //are all criteria fulfilled?
	  _is_masscut_passed = false;
	  if (check_mass_criteria(top_subs)) {
	    _is_masscut_passed = true;
	  }
	  _is_ptmincut_passed = false;
	  if (_top_candidate.pt() > _minpt_tag) {
	    _is_ptmincut_passed = true;
	  }
	}//end better
      }//end kk
    }//end ll
  }//end rr
  return;
}

void HEPTopTaggerV2::get_info() const {  
  std::cout << "#--------------------------------------------------------------------------\n";
  std::cout << "#                          HEPTopTaggerV2 Result" << std::endl;
  std::cout << "#" << std::endl;
  std::cout << "# is top candidate: " << _is_maybe_top << std::endl;
  std::cout << "# mass plane cuts passed: " << _is_masscut_passed << std::endl;
  std::cout << "# top candidate mass: " << _top_candidate.m() << std::endl;
  std::cout << "# top candidate (pt, eta, phi): (" 
       << _top_candidate.perp() << ", "
       << _top_candidate.eta() << ", "
       << _top_candidate.phi_std() << ")" << std::endl;
  std::cout << "# top hadrons: " << _top_hadrons.size() << std::endl;
  std::cout << "# hard substructures: " << _parts_size << std::endl;
  std::cout << "# |m - mtop| : " << _delta_top << std::endl; 
  std::cout << "# djsum : " << _djsum << std::endl;
  std::cout << "# is consistency cut passed: " << _is_ptmincut_passed << std::endl; 
  std::cout << "#--------------------------------------------------------------------------\n";
  return;
}

void HEPTopTaggerV2::get_setting() const {
  std::cout << "#--------------------------------------------------------------------------\n";
  std::cout << "#                         HEPTopTaggerV2 Settings" << std::endl;
  std::cout << "#" << std::endl;
  std::cout << "# mode: " << _mode << " (0 = EARLY_MASSRATIO_SORT_MASS) " << std::endl;
  std::cout << "#        "         << " (1 = LATE_MASSRATIO_SORT_MASS)  " << std::endl;
  std::cout << "#        "         << " (2 = EARLY_MASSRATIO_SORT_MODDJADE)  " << std::endl;
  std::cout << "#        "         << " (3 = LATE_MASSRATIO_SORT_MODDJADE)  " << std::endl;
  std::cout << "#        "         << " (4 = TWO_STEP_FILTER)  " << std::endl;
  std::cout << "# top mass: " << _mtmass << "    ";
  std::cout << "W mass: " << _mwmass << std::endl;
  std::cout << "# top mass window: [" << _mtmin << ", " << _mtmax << "]" << std::endl;
  std::cout << "# W mass ratio: [" << _rmin << ", " << _rmax << "] (["
       <<_rmin*_mtmass/_mwmass<< "%, "<< _rmax*_mtmass/_mwmass << "%])"<< std::endl;
  std::cout << "# mass plane cuts: (m23cut, m13min, m13max) = (" 
       << _m23cut << ", " << _m13cutmin << ", " << _m13cutmax << ")" << std::endl;
  std::cout << "# mass_drop_threshold: " << _mass_drop_threshold << "    ";
  std::cout << "max_subjet_mass: " << _max_subjet_mass << std::endl;
  std::cout << "# R_filt: " << _Rfilt << "    ";
  std::cout << "n_filt: " << _nfilt << std::endl;
  std::cout << "# minimal subjet pt: " << _minpt_subjet << std::endl;
  std::cout << "# minimal reconstructed pt: " << _minpt_tag << std::endl;
  std::cout << "# internal jet algorithms (0 = kt, 1 = C/A, 2 = anti-kt): " << std::endl; 
  std::cout << "#   filtering: "<< _jet_algorithm_filter << std::endl;
  std::cout << "#   reclustering: "<< _jet_algorithm_recluster << std::endl;
  std::cout << "#--------------------------------------------------------------------------\n";
  
  return;
}
// Do not change next line, it's needed by the sed-code that makes the tagger CMSSW-compatible.
};
