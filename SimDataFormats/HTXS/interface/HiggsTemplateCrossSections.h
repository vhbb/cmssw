#ifndef SimDataFormats_GeneratorProducts_HiggsTemplateCrossSections_h
#define SimDataFormats_GeneratorProducts_HiggsTemplateCrossSections_h

#include "DataFormats/Math/interface/LorentzVector.h"
#include <vector>

/// Higgs Template Cross Section namespace
namespace HTXS {

  /// Error code: whether the classification was successful or failed
  enum ErrorCode {
    UNDEFINED=-99,
    SUCCESS = 0,               ///< successful classification
    PRODMODE_DEFINED = 1,      ///< production mode not defined
    MOMENTUM_CONSERVATION = 2, ///< failed momentum conservation
    HIGGS_IDENTIFICATION = 3,  ///< failed to identify Higgs boson
    HS_VTX_IDENTIFICATION = 4, ///< failed to identify hard scatter vertex
    VH_IDENTIFICATION = 5,     ///< failed to identify associated vector boson
    TOP_W_IDENTIFICATION = 6   ///< failed to identify top decay
  };

  /// Higgs production modes, corresponding to input sample
  enum HiggsProdMode {
    UNKNOWN = 0,
    GGF = 1, VBF = 2, WH = 3, QQ2ZH = 4, GG2ZH = 5,
    TTH = 6, BBH = 7, TH = 8 
  };
  
  ///   Two digit number of format PF
  ///   P is digit for the physics process
  ///   and F is 0 for |yH|>2.5 and 11 for |yH|<2.5 ("in fiducial")

  /// Namespace for Stage0 categorization
  namespace Stage0 {
    /// @enum Stage-0 ategorization: Two-digit number of format PF, with P for process and F being 0 for |yH|>2.5 and 1 for |yH|<2.5
    enum Category {
      UNKNOWN  = 0, GG2H_FWDH = 10, GG2H = 11, VBF_FWDH = 20, VBF = 21, VH2HQQ_FWDH = 22, VH2HQQ = 23,
      QQ2HLNU_FWDH = 30, QQ2HLNU = 31, QQ2HLL_FWDH = 40, QQ2HLL = 41, GG2HLL_FWDH = 50, GG2HLL = 51,
      TTH_FWDH = 60, TTH = 61, BBH_FWDH = 70, BBH = 71, TH_FWDH = 80, TH = 81 };
  }
 
  /// Categorization Stage 1:
  /// Three digit integer of format PF
  /// Where P is a digit representing the process
  /// F is a unique integer ( F < 99 ) corresponding to each Stage1 phase-space region (bin)
  namespace Stage1 {
    enum Category {
      UNKNOWN  = 0,
      // Gluon fusion
      GG2H_FWDH = 100,
      GG2H_VBFTOPO_JET3VETO = 101, GG2H_VBFTOPO_JET3 = 102,
      GG2H_0J   = 103,
      GG2H_1J_PTH_0_60 = 104,      GG2H_1J_PTH_60_120 = 105,
      GG2H_1J_PTH_120_200 = 106,   GG2H_1J_PTH_GT200 = 107,
      GG2H_GE2J_PTH_0_60 = 108,      GG2H_GE2J_PTH_60_120 = 109,
      GG2H_GE2J_PTH_120_200 = 110,   GG2H_GE2J_PTH_GT200 = 111,
      // "VBF"
      QQ2HQQ_FWDH = 200,
      QQ2HQQ_VBFTOPO_JET3VETO = 201, QQ2HQQ_VBFTOPO_JET3 = 202,
      QQ2HQQ_VH2JET = 203, QQ2HQQ_REST = 204, QQ2HQQ_PTJET1_GT200 = 205,
      // qq -> WH
      QQ2HLNU_FWDH = 300,
      QQ2HLNU_PTV_0_150 = 301,
      QQ2HLNU_PTV_150_250_0J = 302,
      QQ2HLNU_PTV_150_250_GE1J = 303,
      QQ2HLNU_PTV_GT250 = 304,
      // qq -> ZH
      QQ2HLL_FWDH = 400,
      QQ2HLL_PTV_0_150 = 401,
      QQ2HLL_PTV_150_250_0J = 402,
      QQ2HLL_PTV_150_250_GE1J = 403,
      QQ2HLL_PTV_GT250 = 404,
      // gg -> ZH
      GG2HLL_FWDH = 500,
      GG2HLL_PTV_0_150 = 501,
      GG2HLL_PTV_GT150_0J = 502,
      GG2HLL_PTV_GT150_GE1J = 503,
      // ttH
      TTH_FWDH = 600, TTH = 601,
      // bbH
      BBH_FWDH = 700, BBH = 701,
      // tH
      TH_FWDH = 800, TH = 801
    };
  } // namespace Stage1

  
//#ifdef ROOT_TLorentzVector
    //typedef TLorentzVector TLV;
    typedef math::XYZTLorentzVectorD TLV;
    typedef std::vector<TLV> TLVs;
    
    template <class vec4>
      TLV MakeTLV(vec4 const p) { return TLV(p.px(),p.py(),p.pz(),p.E()); }
    
    template <class Vvec4>
      inline TLVs MakeTLVs(Vvec4 const &rivet_jets){ 
      TLVs jets; for ( auto jet:rivet_jets ) jets.push_back(MakeTLV(jet)); 
      return jets; 
    }
    
    // Structure holding information about the current event:
    // Four-momenta and event classification according to the
    // Higgs Template Cross Section
    struct HiggsClassification {
      // Higgs production mode
      HTXS::HiggsProdMode prodMode;
      // The Higgs boson
      TLV higgs;
      // The Higgs boson decay products
      TLV p4decay_higgs;
      // Associated vector bosons
      TLV V;
      // The V-boson decay products
      TLV p4decay_V;
      // Jets are built ignoring Higgs decay products and leptons from V decays
      // jets with pT > 25 GeV and 30 GeV
      TLVs jets25, jets30;
      // Event categorization according to YR4 wrtietup
      // https://cds.cern.ch/record/2138079
      HTXS::Stage0::Category stage0_cat;
      HTXS::Stage1::Category stage1_cat_pTjet25GeV;
      HTXS::Stage1::Category stage1_cat_pTjet30GeV;
      // Error code :: classification was succesful or some error occured
      HTXS::ErrorCode errorCode;
    };
    
    template <class category>
      inline HTXS::HiggsClassification Rivet2Root(category const &htxs_cat_rivet){
      HTXS::HiggsClassification cat;
      cat.prodMode = htxs_cat_rivet.prodMode;
      cat.errorCode = htxs_cat_rivet.errorCode;
      cat.higgs = MakeTLV(htxs_cat_rivet.higgs);
      cat.V = MakeTLV(htxs_cat_rivet.V);
      cat.p4decay_higgs = MakeTLV(htxs_cat_rivet.p4decay_higgs);
      cat.p4decay_V = MakeTLV(htxs_cat_rivet.p4decay_V);
      cat.jets25 = MakeTLVs(htxs_cat_rivet.jets25);
      cat.jets30 = MakeTLVs(htxs_cat_rivet.jets30);
      cat.stage0_cat = htxs_cat_rivet.stage0_cat;
      cat.stage1_cat_pTjet25GeV = htxs_cat_rivet.stage1_cat_pTjet25GeV;
      cat.stage1_cat_pTjet30GeV = htxs_cat_rivet.stage1_cat_pTjet30GeV;
      return cat;    
    }
    
//#endif 

} // namespace HTXS


#ifdef RIVET_Particle_HH
//#ifdef HIGGSTRUTHCLASSIFIER_HIGGSTRUTHCLASSIFIER_CC
//#include "Rivet/Particle.hh"
namespace Rivet {

  /// @struct HiggsClassification
  /// @brief Structure holding information about the current event:
  ///        Four-momenta and event classification according to the
  ///        Higgs Template Cross Section
  struct HiggsClassification {
    /// Higgs production mode
    HTXS::HiggsProdMode prodMode;
    /// The Higgs boson
    Rivet::Particle higgs;
    /// Vector boson produced in asscoiation with the Higgs
    Rivet::Particle V;
    /// The four momentum sum of all stable decay products orignating from the Higgs boson
    Rivet::FourMomentum p4decay_higgs;
    /// The four momentum sum of all stable decay products orignating from the vector boson in associated production
    Rivet::FourMomentum p4decay_V;
    /// Jets built ignoring Higgs decay products and leptons from V decays, pT thresholds at 25 GeV and 30 GeV
    Rivet::Jets jets25, jets30;
    /// Stage-0 HTXS event classifcation, see: https://cds.cern.ch/record/2138079
    HTXS::Stage0::Category stage0_cat;
    /// Stage-1 HTXS event classifcation, see: https://cds.cern.ch/record/2138079
    HTXS::Stage1::Category stage1_cat_pTjet25GeV;
    /// Stage-1 HTXS event classifcation, see: https://cds.cern.ch/record/2138079
    HTXS::Stage1::Category stage1_cat_pTjet30GeV;
    /// Error code: Whether classification was succesful or some error occured
    HTXS::ErrorCode errorCode;
  };
} // namespace Rivet
#endif



#endif
