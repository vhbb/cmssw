#ifndef AnalysisDataFormats_TopObjects_interface_HTTTopJetTagInfo_h
#define AnalysisDataFormats_TopObjects_interface_HTTTopJetTagInfo_h


// \class HTTTopJetTagInfo
// 
// \short specific tag info for HEPTopTagger tagging algorithm
// HTTTopJetTagInfo is a class to hold the discriminator variables for the
// HEPTopTagger algorithm.
// 
//
// \author Gregor Kasieczka (based on  CATopJetTagInfo by Salvatore Rappoccio)
// \version first version on 25 Sep 2014

#include "DataFormats/BTauReco/interface/RefMacros.h"
#include "DataFormats/BTauReco/interface/JetTagInfo.h"
#include "DataFormats/BTauReco/interface/TaggingVariable.h"
#include <vector>

namespace reco {
 
class HTTTopJetProperties {
public:
  HTTTopJetProperties() {
    fjPt             = 0.;
    fjMass           = 0.;
    fjEta            = 0.;
    fjPhi            = 0.;
    topMass          = 0.;
    unfilteredMass   = 0.;
    prunedMass	     = 0.;
    fW		     = 0.;
    massRatioPassed  = 0.;
    isMultiR	     = 0;
    Rmin	     = 0.;
    RminExpected     = 0.;
    ptFiltForRminExp = 0.;
  }
  double              fjPt;             //<! Mass of the inital Fatjet passed to the TT
  double              fjMass;           //<! Mass of the inital Fatjet passed to the TT
  double              fjEta;            //<! Mass of the inital Fatjet passed to the TT
  double              fjPhi;            //<! Mass of the inital Fatjet passed to the TT
  double              topMass;          //<! Mass of the HTT top quark candidate [GeV] (at R=Rmin for MultiR)
  double              unfilteredMass;   //<! Unfiltered mass of the triplet [GeV] (at R=Rmin for MultiR)
  double              prunedMass;       //<! Mass of the pruned fat jet [GeV] (at R=Rmin for MultiR)
  double              fW;               //<! Minimum distance of m_ij/m_123 from m_W/m_top (at R=Rmin for MultiR)
  double              massRatioPassed;  //<! Did the candidate pass the default mass ratio? Can be used instead of fW (at R=Rmin for MultiR)
  bool                isMultiR;         //<! Tagger operated in MultiR mode
  double              Rmin;             //<! R_min found in MultiR procedure. Set to -1 for non-MultiR mode.
  double              RminExpected;     //<! R_min expected for a top quark based on filtered fat-jet pT. Set to -1 for non-MultiR mode.
  double              ptFiltForRminExp; //<! Filtered initial fatjet pT for re-doing Rmin(expected) fit  Set to -1 for non-MultiR mode.
};

 class HTTTopJetTagInfo : public JetTagInfo {
public:
  typedef edm::RefToBase<Jet> jet_type;
  typedef HTTTopJetProperties  properties_type;
    
    HTTTopJetTagInfo(void) {}

    virtual ~HTTTopJetTagInfo(void) {}
  
    virtual HTTTopJetTagInfo* clone(void) const { return new HTTTopJetTagInfo(*this); }
    
    const properties_type & properties() const {
      return properties_;
    }

    void insert(const edm::RefToBase<Jet> & jet, const HTTTopJetProperties & properties) {
      setJetRef(jet);
      properties_ = properties;
    }

protected:
    properties_type properties_;

};

DECLARE_EDM_REFS( HTTTopJetTagInfo )

}

#endif // AnalysisDataFormats_TopObjects_interface_HTTTopJetTagInfo_h
