// -*- C++ -*-
//
// Package:    BCorrAnalyzer
// Class:      BCorrAnalyzer
// 
/**\class BCorrAnalyzer BCorrAnalyzer.cc BAnalysis/BCorrAnalyzer/src/BCorrAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lukas Wehrli (IPP/ETHZ) [wehrlilu]
//         Created:  Wed May 26 10:41:33 CEST 2010
// $Id: BCorrAnalyzer.cc,v 1.10 2010/07/19 12:15:34 wehrlilu Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//#include "DataFormats/PatCandidates/interface/Jet.h"

//Jet headers
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "RecoJets/JetAlgorithms/interface/JetIDHelper.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

//Trigger headers
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"


#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include "TH1F.h"

//sim headers
#include "AnalysisDataFormats/BAnalysis/interface/SimBHadron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

#include <vector>

#define MAXJETS 100
#define MAXPDF 100

#define runPrescale 136393

#define DEBUG 1

using namespace std;


struct bcandid{
  float massBcand;
  float massVert;
  float gamma;
  float pt; 
  float eta; 
  float phi;
  float dist3D; 
  float distSig3D;
  float dist2D; 
  float distSig2D;
  float dist3D_norm;
  float eventNo; 
  float runNo;
  float pthat; 
  float ptHardestPJ;
  float etaHardestPJ;
  int nvert; 
  int nSelected; 
  int ndRmatched;
  int flavor; //dR matching to SimBHadrons, dRmin<0.1 -> flavor 5
              //else dR matching to D (c)?, dRmin<0.1 -> flavor 4
              //else flavor 1 (light)
  int jet6; 
  int jet10; 
  int jet10nobptx; 
  int jet15; 
  int jet15hcalnf; 
  int jet30; 
  int jet50; 
  int process; //see event
  int eventProcId;
  bool selected; 
  bool dRmatched; 
}; 
/*
struct bcandid{
  float massBcand;
  float massVert;
  float gamma;
  float pt; 
  float eta; 
  float phi;
  float dist3D; 
  float distSig3D;
  float dist2D; 
  float distSig2D;
  float dist3D_norm;
  float eventNo; 
  float runNo;
  float ptHardestPJ; 
  float etaHardestPJ; 
  int nvert; 
  int nSelected; 
  int jet6; 
  int jet10; 
  int jet10nobptx; 
  int jet15; 
  int jet15hcalnf; 
  int jet30; 
  int jet50; 
  bool selected; 
}; 
*/

struct event{
  int nB; 
  int nV; 
  int nMat; 
  int process; //0 if not defined, 1 FCR, 2 FEX, 3 GSP
  int eventProcId;
  float dRvv;
  float dEtavv; 
  float dPhivv;
  float dRbb;
  float dEtabb; 
  float dPhibb;
  float massV1;
  float massV2;
  float massV3;
  float ptV1;
  float ptV2;
  float ptV3;
  float etaV1;
  float etaV2;
  float etaV3;
  float phiV1;
  float phiV2;
  float phiV3;
  float dist3D1; 
  float dist3D2; 
  float dist3D3; 
  float distSig3D1; 
  float distSig3D2;
  float distSig3D3;
  float dist2D1; 
  float dist2D2; 
  float dist2D3; 
  float distSig2D1; 
  float distSig2D2; 
  float distSig2D3;  
  float massB1;
  float massB2;
  float ptB1;
  float ptB2;
  float etaB1;
  float etaB2;
  float phiB1;
  float phiB2;
  float ptHardestGJ; 
  float ptHardestPJ; 
  float ptHardestPGJ; 
  float eventNo;
  float runNo;
  float pthat; 
  int flavors; 
  int jet6; 
  int jet10; 
  int jet10nobptx; 
  int jet15; 
  int jet15hcalnf; 
  int jet30; 
  int jet50; 

};

struct brsvpair{
  unsigned int index1;
  unsigned int index2;
  double discr;
};

inline void removePairs(std::set<brsvpair> &res){
  std::set<brsvpair>::const_iterator it1, it2, temp;
  for(it1=res.begin(); it1!=res.end(); it1++){
    it2=it1;
    for(it2++; it2!=res.end(); it2++){
      if((*it2).index1==(*it1).index1 || (*it2).index2==(*it1).index2){
        temp = it2; it2--;
        res.erase(*temp);
      }
    }
  }
}

inline bool operator<(brsvpair b1, brsvpair b2){return (b1.discr<b2.discr);}

class BCorrAnalyzer : public edm::EDAnalyzer {
   public:
      explicit BCorrAnalyzer(const edm::ParameterSet&);
      ~BCorrAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      GlobalVector flightDirection(reco::Vertex &pv, reco::Vertex &sv);

      // ----------member data ---------------------------
      edm::Service<TFileService> fs_;
      TTree *tEvents, *tBCands;
      event sevents;
      bcandid sbcands; 
      edm::InputTag vertexsrc_, simbsrc_, offlinePrimaryVertices_, trigressrc_;

      double minInvMass, distSig3Dmin, maxEtaVertex;
      TH1F *hdR, *hdEta, *hdPhi; 

      float count[8][8][8][8][8][8]; 
      float norm[8][8][8]; 

      float countJet[8][8][8][8];
      float normJet[8][8];

  
  ////Bjet analysis var
  //input var
  int isData_;
  edm::InputTag JetCollection_;
  edm::InputTag GenJetCollection_;
  edm::InputTag PFJetCollection_;
  edm::InputTag trigTag_;
  const edm::ParameterSet   jetID_;

  reco::helper::JetIDHelper jetIDHelper;

  std::string JetAlgoLabel_;

  float CaloJetSelection_minPt_, CaloJetSelection_maxEta_,  CaloJetSelection_EMF_, CaloJetSelection_fHPD_, CaloJetSelection_n90Hits_;
  float PFJetSelection_minPt_, PFJetSelection_maxEta_,  PFJetSelection_NHF_, PFJetSelection_NEF_,
    PFJetSelection_CEF_, PFJetSelection_CHF_, PFJetSelection_CM_,PFJetSelection_NCONST_ ;
  float BTag_tchel_, BTag_tchem_;
  float BTag_tchpm_, BTag_tchpt_;
  float BTag_ssvm_, BTag_ssvt_;
  float BTag_ssvhp_;

  //an var
  int PtHatEff;
  int TotalEvts;

  int eventRunNumber, eventNumber, eventLumiSection;
  int eventProcId;

  double pthat;

  double pdf_x1, pdf_x2, pdf_pdf1, pdf_pdf2, pdf_Q;
  int pdf_id1, pdf_id2;

  double weight;

  int CaloJet_nJets;
  float CaloJet_pt[MAXJETS], CaloJet_eta[MAXJETS], CaloJet_phi[MAXJETS], CaloJet_E[MAXJETS];
  float CaloJet_px[MAXJETS], CaloJet_py[MAXJETS],CaloJet_pz[MAXJETS];
  int CaloJet_JetIDLoose[MAXJETS], CaloJet_JetIDTight[MAXJETS];
  int CaloJet_JetFlavour[MAXJETS];
  float CaloJet_emEF[MAXJETS], CaloJet_fHPD[MAXJETS], CaloJet_n90Hits[MAXJETS], CaloJet_fRBX[MAXJETS];

  int PFJet_nJets;
  float PFJet_pt[MAXJETS], PFJet_eta[MAXJETS], PFJet_phi[MAXJETS], PFJet_E[MAXJETS];
  float PFJet_px[MAXJETS], PFJet_py[MAXJETS],PFJet_pz[MAXJETS];
  int PFJet_JetIDLoose[MAXJETS], PFJet_JetFlavour[MAXJETS];
  float PFJet_nHEF[MAXJETS], PFJet_nEmEF[MAXJETS], PFJet_cHEF[MAXJETS], PFJet_cEmEF[MAXJETS];
  int PFJet_cMult[MAXJETS], PFJet_nConst[MAXJETS];

  vector<string> WP;
  vector<string> BALGO;
  map<string, float> WPVALUES;
  map<string, float> varFloat;
  map<string, float*> varFloatArr;
  map<string, int*> varIntArr;
  map<string, int> varInt;
  map<string, int> varTrigger;



};

BCorrAnalyzer::BCorrAnalyzer(const edm::ParameterSet& iConfig):
  vertexsrc_(iConfig.getParameter<edm::InputTag>("src")),
  simbsrc_(iConfig.getParameter<edm::InputTag>("simbsrc")),
  offlinePrimaryVertices_(iConfig.getParameter<edm::InputTag>("primaryVertices")),
  trigressrc_(iConfig.getParameter<edm::InputTag>("trigressrc")),
  minInvMass(iConfig.getUntrackedParameter<double>("minInvM",1.4)),
  distSig3Dmin(iConfig.getUntrackedParameter<double>("mindistSig3D",5.0)),
  maxEtaVertex(iConfig.getUntrackedParameter<double>("maxEtaVertex",2.0))

{

  varTrigger["HLT_L1_BscMinBiasOR_BptxPlusORMinus"] = 0;
  varTrigger["HLT_L1_BscMinBiasOR_BptxPlusORMinus_NoBPTX"] = 0;
  varTrigger["HLT_MinBias"] = 0;
  varTrigger["HLT_MinBiasBSC"] = 0;
  varTrigger["HLT_MinBiasBSC_NoBPTX"] = 0;

  varTrigger["HLT_L1_BscMinBiasOR_BptxPlusORMinus_prescale"] = -1;
  varTrigger["HLT_L1_BscMinBiasOR_BptxPlusORMinus_NoBPTX_prescale"] = -1;
  varTrigger["HLT_MinBias_prescale"] = -1;
  varTrigger["HLT_MinBiasBSC_prescale"] = -1;
  varTrigger["HLT_MinBiasBSC_NoBPTX_prescale"] = -1;

  varTrigger["HLT_L1Jet6U"] = 0;
  varTrigger["HLT_L1Jet10U"] = 0;
  varTrigger["HLT_HT100U"] = 0;
  varTrigger["HLT_Jet15U"] = 0;
  varTrigger["HLT_Jet30U"] = 0;
  varTrigger["HLT_Jet50U"] = 0;
  varTrigger["HLT_HT100U"] = 0;

  varTrigger["HLT_L1Jet6U_prescale"] = -1;
  varTrigger["HLT_L1Jet10U_prescale"] = -1;
  varTrigger["HLT_HT100U_prescale"] = -1;
  varTrigger["HLT_Jet15U_prescale"] = -1;
  varTrigger["HLT_Jet30U_prescale"] = -1;
  varTrigger["HLT_Jet50U_prescale"] = -1;
  varTrigger["HLT_HT100U_prescale"] = -1;

  BALGO.push_back("tche");
  BALGO.push_back("tchp");
  BALGO.push_back("ssvhe");
  BALGO.push_back("ssvhp");
  WP.push_back("tchel");
  WP.push_back("tchem");
  WP.push_back("tchpm");
  WP.push_back("tchpt");
  WP.push_back("ssvhem");
  WP.push_back("ssvhet");
  WP.push_back("ssvhp");

  CaloJetSelection_minPt_ = iConfig.getUntrackedParameter<double>("CaloJetSelection_minPt" );
  CaloJetSelection_maxEta_ = iConfig.getUntrackedParameter<double>("CaloJetSelection_maxEta" );
  //not used (CaloJetid selector used instead)
  //CaloJetSelection_EMF_ = iConfig.getUntrackedParameter<double>("CaloJetSelection_EMF" );
  //CaloJetSelection_fHPD_ = iConfig.getUntrackedParameter<double>("CaloJetSelection_fHPD" );
  //CaloJetSelection_n90Hits_ = iConfig.getUntrackedParameter<double>("CaloJetSelection_n90Hits" );

  PFJetSelection_minPt_ = iConfig.getUntrackedParameter<double>("PFJetSelection_minPt" );
  PFJetSelection_maxEta_= iConfig.getUntrackedParameter<double>("PFJetSelection_maxEta" );
  PFJetSelection_NHF_ = iConfig.getUntrackedParameter<double>("PFJetSelection_NHF" );
  PFJetSelection_NEF_ = iConfig.getUntrackedParameter<double>("PFJetSelection_NEF" );
  PFJetSelection_CEF_ = iConfig.getUntrackedParameter<double>("PFJetSelection_CEF" );
  PFJetSelection_CHF_= iConfig.getUntrackedParameter<double>("PFJetSelection_CHF" );
  PFJetSelection_CM_ = iConfig.getUntrackedParameter<double>("PFJetSelection_CM" );
  PFJetSelection_NCONST_ = iConfig.getUntrackedParameter<double>("PFJetSelection_NCONST" );

  for(vector<string>::iterator wp = WP.begin(); wp!=WP.end(); wp++){
    string wpn = (*wp);
    WPVALUES[wpn] = iConfig.getUntrackedParameter<double>("BTag_"+wpn );
  }

  isData_ = iConfig.getUntrackedParameter<int>("isData" );

  JetCollection_ = iConfig.getUntrackedParameter<edm::InputTag>("JetCollection" );
  PFJetCollection_ = iConfig.getUntrackedParameter<edm::InputTag>("PFJetCollection" );

  for(vector<string>::iterator wp = WP.begin(); wp!=WP.end(); wp++){
    string wpn = (*wp);
    varFloatArr["CaloJet_btag_"+wpn+"_pt"] = new float[MAXJETS];
    varFloatArr["CaloJet_btag_"+wpn+"_eta"] = new float[MAXJETS];
    varFloatArr["CaloJet_btag_"+wpn+"_phi"] = new float[MAXJETS];
    varFloatArr["CaloJet_btag_"+wpn+"_discr"] = new float[MAXJETS];
    varIntArr["CaloJet_btag_"+wpn+"_JetFlavour"] = new int[MAXJETS];

    varInt["CaloJet_btag_"+wpn+"_nsel"] = 0;
    varFloat["CaloJet_deltaPhi_2bjet_"+wpn] = -99;
    varFloat["CaloJet_deltaR_2bjet_"+wpn] = -99;

    varFloatArr["PFJet_btag_"+wpn+"_pt"] = new float[MAXJETS];
    varFloatArr["PFJet_btag_"+wpn+"_eta"] = new float[MAXJETS];
    varFloatArr["PFJet_btag_"+wpn+"_phi"] = new float[MAXJETS];
    varFloatArr["PFJet_btag_"+wpn+"_discr"] = new float[MAXJETS];
    varIntArr["PFJet_btag_"+wpn+"_JetFlavour"] = new int[MAXJETS];

    varInt["PFJet_btag_"+wpn+"_nsel"] = 0;
    varFloat["PFJet_deltaPhi_2bjet_"+wpn] = -99;
    varFloat["PFJet_deltaR_2bjet_"+wpn] = -99;

  }

  for(vector<string>::iterator balg = BALGO.begin(); balg!=BALGO.end(); balg++){
    varFloatArr["CaloJet_btag_"+(*balg)] = new float[MAXJETS];
    varFloatArr["PFJet_btag_"+(*balg)] = new float[MAXJETS];
  }

}


BCorrAnalyzer::~BCorrAnalyzer()
{
 
}


// ------------ method called to for each event  ------------
void
BCorrAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;

   //Run info
   eventLumiSection = iEvent.luminosityBlock();
   EventID eventId = iEvent.id();
   eventRunNumber = eventId.run ();
   eventNumber = eventId.event ();

   //Gen info
   Handle<SimBHadronCollection> SBHC;
   SimBHadronCollection sbhc ;
   Handle<GenParticleCollection> GPC;
   GenParticleCollection gpc;
   Handle<GenJetCollection> GJC;
   GenJetCollection gjc;

   double ptHardestGJ = -1;
   int proc1 = -1, proc2 = -1;  //1FCR, 2FEX, 3GSP, 4FourB

   if(isData_==0){
     Handle< GenEventInfoProduct > GenInfoHandle;
     iEvent.getByLabel( "generator", GenInfoHandle );
     eventProcId   = GenInfoHandle->signalProcessID();
     pthat = ( GenInfoHandle->hasBinningValues() ?
               (GenInfoHandle->binningValues())[0] : 0.0);

     //weight = 1;
     //double event_weight = GenInfoHandle->weight();
     //if(event_weight!=1) weight = event_weight;

     pdf_Q = GenInfoHandle->pdf()->scalePDF;
     pdf_id1 = GenInfoHandle->pdf()->id.first;
     pdf_x1 = GenInfoHandle->pdf()->x.first;
     pdf_pdf1 = GenInfoHandle->pdf()->xPDF.first;
     pdf_id2 = GenInfoHandle->pdf()->id.second;
     pdf_x2 = GenInfoHandle->pdf()->x.second;
     pdf_pdf2 = GenInfoHandle->pdf()->xPDF.second;

     //------------------SimBHadrons------------------------------
     if(DEBUG==1) cout << "[DEBUG] SimH init" << endl;
     iEvent.getByLabel(simbsrc_, SBHC);
     sbhc = *(SBHC.product());
     if(sbhc.size()==2){
       if(sbhc[0].quarkstatus + sbhc[0].brotherstatus==6) proc1 = 1;
       if(sbhc[0].quarkstatus + sbhc[0].brotherstatus==5) proc1 = 2;
       if(sbhc[0].quarkstatus + sbhc[0].brotherstatus==4) proc1 = 3;
       if(sbhc[1].quarkstatus + sbhc[1].brotherstatus==6) proc2 = 1;
       if(sbhc[1].quarkstatus + sbhc[1].brotherstatus==5) proc2 = 2;
       if(sbhc[1].quarkstatus + sbhc[1].brotherstatus==4) proc2 = 3;
       if(proc1 != proc2) std::cout << "PROBLEM WITH PROCESS\n\n";
     }
     if(DEBUG==1) cout << "[DEBUG] SimH end" << endl;


     //------------------GenParticles------------------------------
     //Handle<GenParticleCollection> GPC;
     iEvent.getByLabel("genParticles", GPC);
     gpc = *(GPC.product());
     
     //------------------GenJets------------------------------
     iEvent.getByLabel("ak5GenJets", GJC);
     gjc = *(GJC.product());
     ///////////////////////////////////////////
     for(unsigned int i=0; i<gjc.size(); i++)
       if(gjc[i].pt()>ptHardestGJ) ptHardestGJ = gjc[i].pt();
     
   }




   //------------------Trigger Results -----------------------
   edm::Handle<edm::TriggerResults> hltTR;
   iEvent.getByLabel(trigressrc_, hltTR);

   bool wantTable = false; 

   edm::TriggerNames const&  triggerNames = iEvent.triggerNames(*hltTR);
   int number6 = 777, number10 = 777, number10nobptx = 777, number15hcalnf = 777, number15 = 777, number30 = 777, number50 = 777; 
   unsigned int b6 = 0, b10 = 0, b10nobptx = 0, b15hcalnf = 0, b15 = 0, b30 = 0, b50 = 0; 

   for (unsigned i = 0; i < triggerNames.size(); ++i) {
     std::string trigName = triggerNames.triggerName(i);
     varTrigger[trigName] = 0;
     if ( hltTR->accept(i)) varTrigger[trigName] = 1;

     if(trigName=="HLT_L1Jet6U") number6 = i;
     if(trigName=="HLT_L1Jet10U") number10 = i;
     if(trigName=="HLT_L1Jet10U_NoBPTX") number10nobptx = i;
     if(trigName=="HLT_Jet15U_HcalNoiseFiltered") number15hcalnf = i;
     if(trigName=="HLT_Jet15U") number15 = i;
     if(trigName=="HLT_Jet30U") number30 = i;
     if(trigName=="HLT_Jet50U") number50 = i;

     if(wantTable){
       std::cout << i << " " << triggerNames.triggerName(i) << " status: " << hltTR->accept(i) << std::endl;
     }
   }

   //prescales
   if(isData_==1){
     Handle<pat::TriggerEvent> triggerEvent;
     iEvent.getByLabel("patTriggerEvent", triggerEvent);
     pat::TriggerEvent const * trig = &*triggerEvent;

     if ( trig->wasRun() && trig->wasAccept() ) {
       for(map<string,int>::iterator hlt=varTrigger.begin(); hlt!=varTrigger.end(); hlt++){
	 pat::TriggerPath const * Path = trig->path( hlt->first );
	 if(!Path) continue;
	 varTrigger[hlt->first+"_prescale"] = -1;
	 if(!Path->wasRun()) continue;
	 if(eventRunNumber < runPrescale || Path->wasError() ) continue;
	 varTrigger[ hlt->first+"_prescale"] = Path->prescale();
       }
     }
   }
     
   if(number6==777) b6 = 2; 
   else if(hltTR->accept(number6)) b6 = 1; 
   if(number10==777) b10 = 2; 
   else if(hltTR->accept(number10)) b10 = 1; 
   if(number10nobptx==777) b10nobptx = 2; 
   else if(hltTR->accept(number10nobptx)) b10nobptx = 1; 
   if(number15==777) b15 = 2; 
   else if(hltTR->accept(number15)) b15 = 1; 
   if(number15hcalnf==777) b15hcalnf = 2; 
   else if(hltTR->accept(number15hcalnf)) b15hcalnf = 1; 
   if(number30==777) b30 = 2; 
   else if(hltTR->accept(number30)) b30 = 1; 
   if(number50==777) b50 = 2; 
   else if(hltTR->accept(number50)) b50 = 1; 


   //------------------Primary Vertices-----------------------
   edm::Handle<reco::VertexCollection> primaryVertexCollection;
   iEvent.getByLabel(offlinePrimaryVertices_, primaryVertexCollection);
   const reco::VertexCollection pvc = *(primaryVertexCollection.product());
   Vertex pv = pvc[0];
   

   //------------------PFJets-----------------------
   double ptHardestPJ = -1, ptHardestPGJ = -1;
   double etaHardestPJ = -1, etaHardestPGJ = -1;

   Handle< View< pat::Jet > > pfjets;
   iEvent.getByLabel(PFJetCollection_,pfjets);

   PFJet_nJets=0;
   map<string, vector<pat::Jet> > bpfjets;

   for(edm::View<pat::Jet>::const_iterator jet=pfjets->begin(); jet!=pfjets->end(); ++jet){

     double jetPt=jet->pt(); double jetEta=jet->eta();

     PFJet_JetIDLoose[PFJet_nJets] = 0;
     if( jetPt <= PFJetSelection_minPt_) continue;
     else if( fabs(jetEta) >= PFJetSelection_maxEta_ ) continue ;

     //JetID
     bool selectedJet = true;
     if( fabs(jetEta)<2.4){
       if(jet->neutralHadronEnergyFraction() >= PFJetSelection_NHF_ ) selectedJet = false;
       else if(jet->neutralEmEnergyFraction() >= PFJetSelection_NEF_ ) selectedJet = false;
       else if(jet->chargedHadronEnergyFraction() <= PFJetSelection_CHF_   ) selectedJet = false;
       else if(jet->chargedEmEnergyFraction() >=  PFJetSelection_CEF_   ) selectedJet = false;
       else if(jet->chargedMultiplicity() <= PFJetSelection_CM_  ) selectedJet = false;
       else if(jet->nConstituents() <= PFJetSelection_NCONST_  ) selectedJet = false;
     }
     else{
       if(jet->neutralHadronEnergyFraction() >= PFJetSelection_NHF_ ) selectedJet = false;
       else if(jet->neutralEmEnergyFraction() >= PFJetSelection_NEF_ ) selectedJet = false;
       else if(jet->nConstituents() <= PFJetSelection_NCONST_  ) selectedJet = false;
     }
     
     if(selectedJet)PFJet_JetIDLoose[PFJet_nJets] = 1;
     
     //GenJet maxPt check
     if(jet->genJet())
       if(jet->genJet()->pt()>ptHardestPGJ){ ptHardestPGJ = jet->genJet()->pt(); etaHardestPGJ = jet->genJet()->eta();} 
     if(jet->pt()>ptHardestPJ) {ptHardestPJ = jet->pt(); etaHardestPJ = jet->eta();}
       
     
     //Var filling
     PFJet_nHEF[PFJet_nJets] = jet->neutralHadronEnergyFraction();
     PFJet_nEmEF[PFJet_nJets] =  jet->neutralEmEnergyFraction();
     PFJet_cHEF[PFJet_nJets] = jet->chargedHadronEnergyFraction();
     PFJet_cEmEF[PFJet_nJets] =  jet->chargedEmEnergyFraction();
     PFJet_cMult[PFJet_nJets] =  jet->chargedMultiplicity();
     PFJet_nConst[PFJet_nJets] =  jet->nConstituents();
     
     PFJet_px[PFJet_nJets] = jet->px();  PFJet_py[PFJet_nJets] = jet->py(); PFJet_pz[PFJet_nJets] = jet->pz();
     PFJet_pt[PFJet_nJets] = jetPt; PFJet_eta[PFJet_nJets] = jetEta; PFJet_phi[PFJet_nJets] = jet->phi();
     PFJet_E[PFJet_nJets] = jet->energy();
     
     if(isData_==0) PFJet_JetFlavour[PFJet_nJets] =jet->partonFlavour();
     
     //BTag filling
     for(vector<string>::iterator balg = BALGO.begin(); balg!=BALGO.end(); balg++){
       string balgo = (*balg);
       if(balgo=="tche") varFloatArr["PFJet_btag_"+balgo][PFJet_nJets] = jet->bDiscriminator("trackCountingHighEffBJetTags");
       else if(balgo=="tchp") varFloatArr["PFJet_btag_"+balgo][PFJet_nJets] = jet->bDiscriminator("trackCountingHighPurBJetTags");
	 else if(balgo=="ssvhe") varFloatArr["PFJet_btag_"+balgo][PFJet_nJets] = jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
	 else if(balgo=="ssvhp") varFloatArr["PFJet_btag_"+balgo][PFJet_nJets] = jet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
	 else cout << "WTF???" << endl;
     }
     
     string balgo = "";
     for(vector<string>::iterator wp = WP.begin(); wp!=WP.end(); wp++){
       if((*wp)=="tchem" || (*wp)=="tchel") balgo = "tche";
       else if ((*wp)=="tchpm" || (*wp)=="tchpt") balgo = "tchp";
       else if ((*wp)=="ssvhem" || (*wp)=="sshhet") balgo = "ssvhe";
       else if ((*wp)=="ssvhp" ) balgo = "ssvhp";
       
       if( varFloatArr["PFJet_btag_"+balgo][PFJet_nJets] > WPVALUES[(*wp)]) bpfjets[(*wp)].push_back((*jet));
       }
   
     PFJet_nJets++;
   }//end jets


   //bjets quantity
   for(vector<string>::iterator wp = WP.begin(); wp!=WP.end(); wp++){
     string wpn = (*wp);
     
     for(int b=0;b<varInt["PFJet_btag_"+wpn+"_nsel"];b++){
       varFloatArr["PFJet_btag_"+wpn+"_pt"][b] = bpfjets[wpn][b].pt();
       varFloatArr["PFJet_btag_"+wpn+"_eta"][b] = bpfjets[wpn][b].eta();
       varFloatArr["PFJet_btag_"+wpn+"_phi"][b] = bpfjets[wpn][b].phi();
       varIntArr["PFJet_btag_"+wpn+"_JetFlavour"][b] = bpfjets[wpn][b].partonFlavour();
       string discr = "";
       if(wpn=="tchel" || wpn=="tchem") discr ="trackCountingHighEffBJetTags";
       else if(wpn=="tchpm" || wpn=="tchpt") discr = "trackCountingHighPurBJetTags";
       else if(wpn=="ssvpm" || wpn=="ssvpt") discr = "simpleSecondaryVertexHighEffBJetTags";
       else if(wpn=="ssvhp") discr = "simpleSecondaryVertexHighPurBJetTags";
       varFloatArr["PFJet_btag_"+wpn+"_discr"][b] = bpfjets[wpn][b].bDiscriminator(discr.c_str());
     }
   }


   //------------------CaloJets-----------------------
   JetIDSelectionFunctor jetIDLoose( JetIDSelectionFunctor::PURE09,
				     JetIDSelectionFunctor::LOOSE );
   pat::strbitset ret = jetIDLoose.getBitTemplate();

   JetIDSelectionFunctor jetIDTight( JetIDSelectionFunctor::PURE09,
				     JetIDSelectionFunctor::TIGHT );
   pat::strbitset retTight = jetIDTight.getBitTemplate();

   Handle< View< pat::Jet > > jets;
   iEvent.getByLabel(JetCollection_,jets);
   CaloJet_nJets=0;
   map<string, vector<pat::Jet> > bjets;
   for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
     double jetPt=jet->pt(); double jetEta=jet->eta();

     //JetID
     CaloJet_JetIDLoose[CaloJet_nJets] = 0;
     CaloJet_JetIDTight[CaloJet_nJets] = 0;

     if( jetPt <= CaloJetSelection_minPt_) continue;
     else if( fabs(jetEta) >= CaloJetSelection_maxEta_ ) continue ;

     ret.set(false);
     bool selectedJet = jetIDLoose(*jet, ret);
     if(selectedJet)CaloJet_JetIDLoose[CaloJet_nJets] = 1;
     
     retTight.set(false);
     if( jetIDTight(*jet, retTight)) CaloJet_JetIDTight[CaloJet_nJets] = 1;;

     //Var Filling
     CaloJet_emEF[CaloJet_nJets] = jet->emEnergyFraction();
     CaloJet_fHPD[CaloJet_nJets] = jet->jetID().fHPD;
     CaloJet_n90Hits[CaloJet_nJets] = jet->jetID().n90Hits;
     CaloJet_fRBX[ CaloJet_nJets] = jet->jetID().fRBX;
     CaloJet_px[CaloJet_nJets] = jet->px();  CaloJet_py[CaloJet_nJets] = jet->py(); CaloJet_pz[CaloJet_nJets] = jet->pz();
     CaloJet_pt[CaloJet_nJets] = jetPt; CaloJet_eta[CaloJet_nJets] = jetEta; CaloJet_phi[CaloJet_nJets] = jet->phi();
     CaloJet_E[CaloJet_nJets] = jet->energy();
     CaloJet_emEF[CaloJet_nJets] = jet->emEnergyFraction(); CaloJet_fHPD[CaloJet_nJets] = jet->jetID().fHPD; CaloJet_n90Hits[CaloJet_nJets] = jet->jetID().n90Hits;
     
     if(isData_==0) CaloJet_JetFlavour[CaloJet_nJets] =jet->partonFlavour();
     else CaloJet_JetFlavour[CaloJet_nJets] = 0;

     //Btag
     for(vector<string>::iterator balg = BALGO.begin(); balg!=BALGO.end(); balg++){
       string balgo = (*balg);
       if(balgo=="tche") {varFloatArr["CaloJet_btag_"+balgo][CaloJet_nJets] = jet->bDiscriminator("trackCountingHighEffBJetTags");}
       else if(balgo=="tchp"){ varFloatArr["CaloJet_btag_"+balgo][CaloJet_nJets] = jet->bDiscriminator("trackCountingHighPurBJetTags");}
       else if(balgo=="ssvhe"){ varFloatArr["CaloJet_btag_"+balgo][CaloJet_nJets] = jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");}
       else if(balgo=="ssvhp"){ varFloatArr["CaloJet_btag_"+balgo][CaloJet_nJets] = jet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");}
       else
	 cout << "WTF???" << endl;
     }
     string balgo = "";
     for(vector<string>::iterator wp = WP.begin(); wp!=WP.end(); wp++){
       if((*wp)=="tchem" || (*wp)=="tchel") balgo = "tche";
       else if ((*wp)=="tchpm" || (*wp)=="tchpt") balgo = "tchp";
       else if ((*wp)=="ssvhem" || (*wp)=="sshhet") balgo = "ssvhe";
       else if ((*wp)=="ssvhp" ) balgo = "ssvhp";
       
       if( varFloatArr["CaloJet_btag_"+balgo][CaloJet_nJets] > WPVALUES[(*wp)]) bjets[(*wp)].push_back((*jet));
     }
     CaloJet_nJets++;
   }//end jets
   

   //------------------Secondary Vertices-----------------------
   Handle<std::vector<reco::Vertex> > SVC;
   iEvent.getByLabel(vertexsrc_, SVC);
   const std::vector<reco::Vertex> svc = *(SVC.product());

   //------------------Leaf Candidates--------------------------
   Handle<std::vector<reco::LeafCandidate> > LCC;
   iEvent.getByLabel(vertexsrc_, LCC);
   const std::vector<reco::LeafCandidate> lcc = *(LCC.product());


   if(lcc.size()!=svc.size()){
	std::cout << "COLLECTIONS SIZE MISMATCH\n\n";
	return; 
   }

   int goodVert=0, index1 = -1, index2 = -1, index3 = -1; 
   //index3: if three selected it is third selected, otherwise if two selected it is the first remaining one
   std::set<brsvpair> pairs;

   //   std::vector<bcandidmc> bcandids;
   std::vector<bcandid> bcandids;
   
   for(unsigned int i=0; i<lcc.size(); i++){
     //if(isData_==0)bcandidmc curbc;
     //    else 
     bcandid curbc;
     
     GlobalVector axis = GlobalVector(svc[i].p4().X(),svc[i].p4().Y(),svc[i].p4().Z());
     Vertex rv = svc[i]; 
     SecondaryVertex sv(pv,rv,axis,true);
     curbc.selected=false; 
     curbc.dRmatched = false;
     if(lcc[i].p4().M()>minInvMass && sv.dist3d().significance()>distSig3Dmin && fabs(flightDirection(pv,rv).eta())<maxEtaVertex) {
       goodVert++; 
       if(index1==-1) index1=i; 
       else {
	 if(index2==-1) index2=i;
	 else index3=i; 
       }
       curbc.selected = true;
     }
     GlobalVector fldir = flightDirection(pv,rv); 
     curbc.massBcand = lcc[i].p4().M(); 
     curbc.massVert = svc[i].p4().M();
     curbc.gamma = lcc[i].p4().Gamma(); 
     curbc.pt = lcc[i].p4().Pt();  
     curbc.eta = fldir.eta();  
     curbc.phi = fldir.phi(); 
     curbc.dist3D = sv.dist3d().value();  
     curbc.distSig3D = sv.dist3d().significance(); 
     curbc.dist2D = sv.dist2d().value();  
     curbc.distSig2D = sv.dist2d().significance(); 
     curbc.dist3D_norm = sv.dist3d().value()/curbc.gamma; 
     curbc.eventNo = iEvent.id().event();  
     curbc.runNo = iEvent.id().run(); 
     curbc.pthat = pthat; 
     curbc.nvert = (int)lcc.size();
     curbc.flavor = 1;
     //flavor:
     
     for(unsigned int b=0; b<sbhc.size(); b++){
       double dRbr = deltaR(sbhc[b].p4(),fldir); 
       //add an entry in the vecindpair
       brsvpair pair = {b,i,dRbr};
       pairs.insert(pair);
     } 
     bcandids.push_back(curbc);
     
   }
   if(index3==-1 && lcc.size()>2 && goodVert==2){
     for(int i=0; i<(int)lcc.size(); i++){
       if(index3==-1 && index2!=i && index1!=i){
	 index3 = i; 
       }
     }
   }
   //flavor: 
   removePairs(pairs);
   for(std::set<brsvpair>::const_iterator it = pairs.begin(); it!=pairs.end(); it++){
     if(it->discr<0.1) bcandids[it->index2].flavor = 5; 
     //std::cout << "b " << it->index1 << " v " << it->index2 << " dr " << it->discr << "\n";
   }
   
   //charm? 
   if(DEBUG==1) cout << "[DEBUG] charm init " << endl;

   for(unsigned int i=0; i<bcandids.size(); i++){
     if(bcandids[i].flavor!=5){
       for(unsigned int gp=0; gp<gpc.size(); gp++){
	 unsigned int id = abs(gpc[gp].pdgId());
	 if(id==4 || (id/100)==4 || (id/1000)==4) {
	   double dR = deltaR(gpc[gp].p4(),lcc[i].p4());
	     if(dR<0.1) bcandids[i].flavor=4;
	 }
       }
     }
   }
   if(DEBUG==1) cout << "[DEBUG] charm end" << endl;





   sevents.nB =sbhc.size();
   sevents.nV = goodVert; 
   sevents.nMat = -1;
   sevents.process = proc1; 
   sevents.eventProcId = eventProcId; 
   sevents.dRvv = -777.0;
   sevents.dEtavv = -777.0; 
   sevents.dPhivv = -777.0;
   sevents.dRbb = -777.0;
   sevents.dEtabb = -777.0; 
   sevents.dPhibb = -777.0;
   sevents.massV1 = -777.0;
   sevents.massV2 = -777.0;
   sevents.massV3 = -777.0;
   sevents.ptV1 = -777.0;
   sevents.ptV2 = -777.0;
   sevents.ptV3 = -777.0;
   sevents.etaV1 = -777.0;
   sevents.etaV2 = -777.0;
   sevents.etaV3 = -777.0;
   sevents.phiV1 = -777.0;
   sevents.phiV2 = -777.0;
   sevents.phiV3 = -777.0;
   sevents.dist3D1 = -777.0;
   sevents.dist3D2 = -777.0;
   sevents.dist3D3 = -777.0;
   sevents.distSig3D1 = -777.0;
   sevents.distSig3D2 = -777.0;
   sevents.distSig3D3 = -777.0;
   sevents.dist2D1 = -777.0;
   sevents.dist2D2 = -777.0;
   sevents.dist2D3 = -777.0;
   sevents.distSig2D1 = -777.0;
   sevents.distSig2D2 = -777.0;
   sevents.distSig2D3 = -777.0;
   sevents.massB1 = -777.0;
   sevents.massB2 = -777.0;
   sevents.ptB1 = -777.0;
   sevents.ptB2 = -777.0;
   sevents.etaB1 = -777.0;
   sevents.etaB2 = -777.0;
   sevents.phiB1 = -777.0;
   sevents.phiB2 = -777.0;
   sevents.ptHardestGJ = ptHardestGJ; 
   sevents.ptHardestPJ = ptHardestPJ; 
   sevents.ptHardestPGJ = ptHardestPGJ; 
   sevents.eventNo = iEvent.id().event();
   sevents.runNo = iEvent.id().run();
   sevents.pthat = pthat; 
   sevents.flavors = 0; 
   sevents.jet6 = b6; 
   sevents.jet10 = b10; 
   sevents.jet10nobptx = b10nobptx; 
   sevents.jet15 = b15; 
   sevents.jet15hcalnf = b15hcalnf; 
   sevents.jet30 = b30; 
   sevents.jet50 = b50; 


   if(sbhc.size()==2){
     sevents.dRbb = deltaR(sbhc[0].p4(), sbhc[1].p4());
     sevents.dEtabb = sbhc[0].p4().Eta()-sbhc[1].p4().Eta();
     sevents.dPhibb = deltaPhi(sbhc[0].p4().Phi(), sbhc[1].p4().Phi());
     sevents.massB1 = sbhc[0].p4().M();
     sevents.massB2 = sbhc[1].p4().M();
     sevents.ptB1 = sbhc[0].p4().Pt();
     sevents.ptB2 = sbhc[1].p4().Pt();
     sevents.etaB1 = sbhc[0].p4().Eta();
     sevents.etaB2 = sbhc[1].p4().Eta();
     sevents.phiB1 = sbhc[0].p4().Phi();
     sevents.phiB2 = sbhc[1].p4().Phi();

   }
   int nmat = -777; 
   
   if(goodVert==2){
     Vertex sv1 = svc[index1]; 
     Vertex sv2 = svc[index2];
     GlobalVector flightDir1 = flightDirection(pv,sv1); 
     GlobalVector flightDir2 = flightDirection(pv,sv2); 

     sevents.dRvv = deltaR(flightDir1,flightDir2);
     sevents.dEtavv = flightDir1.eta()-flightDir2.eta(); 
     sevents.dPhivv = deltaPhi(flightDir1.phi(),flightDir2.phi());
     sevents.massV1 = lcc[index1].p4().M();
     sevents.massV2 = lcc[index2].p4().M();
     sevents.ptV1 = lcc[index1].p4().Pt();
     sevents.ptV2 = lcc[index2].p4().Pt();
     sevents.etaV1 = flightDir1.eta();
     sevents.etaV2 = flightDir2.eta();
     sevents.phiV1 = flightDir1.phi();
     sevents.phiV2 = flightDir2.phi();
     sevents.dist3D1 = bcandids[index1].dist3D; 
     sevents.dist3D2 = bcandids[index2].dist3D; 
     sevents.distSig3D1 = bcandids[index1].distSig3D; 
     sevents.distSig3D2 = bcandids[index2].distSig3D; 
     sevents.dist2D1 = bcandids[index1].dist2D; 
     sevents.dist2D2 = bcandids[index2].dist2D; 
     sevents.distSig2D1 = bcandids[index1].distSig2D; 
     sevents.distSig2D2 = bcandids[index2].distSig2D; 

     //flavor: 
     int flav1 = bcandids[index1].flavor; 
     int flav2 = bcandids[index2].flavor; 
     if(flav1<flav2){
       int temp = flav2; flav2 = flav1; flav1 = temp; 
     }
     sevents.flavors = 10*flav1 + flav2; 

     hdR->Fill(sevents.dRvv); 
     hdEta->Fill(sevents.dEtavv); 
     hdPhi->Fill(sevents.dPhivv); 

     if(sbhc.size()==2){
       double min = -1; 
       int pos = 0;
       nmat = 0; 
       double dR11 = deltaR(sbhc[0].p4(),flightDir1); 
       min=dR11;
       double dR12 = deltaR(sbhc[0].p4(),flightDir2); 
       if(dR12<min){min=dR12; pos=1;}
       double dR21 = deltaR(sbhc[1].p4(),flightDir1);
       if(dR21<min){min=dR21; pos=2;}
       double dR22 = deltaR(sbhc[1].p4(),flightDir2); 
       if(dR22<min){min=dR22; pos=3;}
       if(min<0.1){
	 if(pos==0 || pos==2) bcandids[index1].dRmatched = true;
	 else bcandids[index2].dRmatched = true;
	 nmat++;
	 switch(pos){
	 case 0: {if(dR22<0.1) {nmat++; bcandids[index2].dRmatched = true;} break;}
	 case 1: {if(dR21<0.1) {nmat++; bcandids[index1].dRmatched = true;} break;}
	 case 2: {if(dR12<0.1) {nmat++; bcandids[index2].dRmatched = true;} break;}
	 case 3: {if(dR11<0.1) {nmat++; bcandids[index1].dRmatched = true;} break;}
	 }
       }
       sevents.nMat = nmat; 
     }
  
   }

   if(index3!=-1){
     Vertex sv3 = svc[index3]; 
     GlobalVector flightDir3 = flightDirection(pv,sv3); 

     sevents.massV3 = lcc[index3].p4().M();
     sevents.ptV3 = lcc[index3].p4().Pt();
     sevents.etaV3 = flightDir3.eta();
     sevents.phiV3 = flightDir3.phi();
     sevents.dist3D3 = bcandids[index3].dist3D; 
     sevents.distSig3D3 = bcandids[index3].distSig3D; 
     sevents.dist2D3 = bcandids[index3].dist2D; 
     sevents.distSig2D3 = bcandids[index3].distSig2D; 

   }

   tEvents->Fill();

   for(unsigned int i=0; i<bcandids.size(); i++){
     sbcands = bcandids[i]; 
     sbcands.nSelected = goodVert; 
     sbcands.ndRmatched = nmat; 
     sbcands.jet6 = b6; 
     sbcands.jet10 = b10; 
     sbcands.jet10nobptx = b10nobptx; 
     sbcands.jet15 = b15; 
     sbcands.jet15hcalnf = b15hcalnf; 
     sbcands.jet30 = b30; 
     sbcands.jet50 = b50;
     sbcands.process = proc1;
     sbcands.eventProcId = eventProcId;
     tBCands->Fill();
   }

   //EFFIC CONSTRUCT
   if(sbhc.size()==2 && fabs(sbhc[0].eta())<2.4 && fabs(sbhc[1].eta())<2.4){
     unsigned int binSDR = 7;
     if(sevents.dRbb<3.0){
       binSDR = 6; 
       if(sevents.dRbb<2.5){
	 binSDR = 5;
	 if(sevents.dRbb<1.2){
	   binSDR = 4;
	   if(sevents.dRbb<0.8){
	     binSDR = 3;
	     if(sevents.dRbb<0.4){
	       binSDR = 2;
	       if(sevents.dRbb<0.2){
		 binSDR = 1;
		 if(sevents.dRbb<0.1){
		   binSDR = 0;
		 }
	       }
	     }
	   }
	 }
       }
     }
     unsigned int binSpts = 7; 
     if(sevents.ptB1<100 || sevents.ptB2<100){
       binSpts = 6; 
       if(sevents.ptB1<70 || sevents.ptB2<70){
	 binSpts = 5; 
	 if(sevents.ptB1<50 || sevents.ptB2<50){
	   binSpts = 4; 
	   if(sevents.ptB1<40 || sevents.ptB2<40){
	     binSpts = 3; 
	     if(sevents.ptB1<30 || sevents.ptB2<30){
	       binSpts = 2; 
	       if(sevents.ptB1<20 || sevents.ptB2<20){
		 binSpts = 1; 
		 if(sevents.ptB1<10 || sevents.ptB2<10){
		   binSpts = 0; 
		 }
	       }
	     }
	   }
	 }
       }
     }
     unsigned int binSpth = 7; 
     if(sevents.ptB1<100 && sevents.ptB2<100){
       binSpth = 6; 
       if(sevents.ptB1<70 && sevents.ptB2<70){
	 binSpth = 5; 
	 if(sevents.ptB1<50 && sevents.ptB2<50){
	   binSpth = 4; 
	   if(sevents.ptB1<40 && sevents.ptB2<40){
	     binSpth = 3; 
	     if(sevents.ptB1<30 && sevents.ptB2<30){
	       binSpth = 2; 
	       if(sevents.ptB1<20 && sevents.ptB2<20){
		 binSpth = 1; 
		 if(sevents.ptB1<10 && sevents.ptB2<10){
		   binSpth = 0; 
		 }
	       }
	     }
	   }
	 }
       }
     }
     unsigned int binSptJET = 7; 
     if(ptHardestGJ<100){
       binSptJET = 6; 
       if(ptHardestGJ<70){
	 binSptJET = 5; 
	 if(ptHardestGJ<50){
	   binSptJET = 4; 
	   if(ptHardestGJ<40){
	     binSptJET = 3; 
	     if(ptHardestGJ<30){
	       binSptJET = 2; 
	       if(ptHardestGJ<20){
		 binSptJET = 1; 
		 if(ptHardestGJ<10){
		   binSptJET = 0; 
		 }
	       }
	     }
	   }
	 }
       }
     }

//      if(sevents.ptB1<20 && sevents.ptB2<20){
//        binSpth = 1; 
//        if(sevents.ptB1<10 && sevents.ptB2<10){
// 	 binSpth = 0; 
//        }
//      }  
     norm[binSDR][binSpts][binSpth]++;
     normJet[binSDR][binSptJET]++;

     //VERTICES
     if(goodVert==2){
       unsigned int binRDR = 7;
       if(sevents.dRvv<3.0){
	 binRDR = 6; 
	 if(sevents.dRvv<2.5){
	   binRDR = 5;
	   if(sevents.dRvv<1.2){
	     binRDR = 4;
	     if(sevents.dRvv<0.8){
	       binRDR = 3;
	       if(sevents.dRvv<0.4){
		 binRDR = 2;
		 if(sevents.dRvv<0.2){
		   binRDR = 1;
		   if(sevents.dRvv<0.1){
		     binRDR = 0;
		   }
		 }
	       }
	     }
	   }
	 }
       }
       
       unsigned int binRpts = 7; 
       if(sevents.ptV1<100 || sevents.ptV2<100){
	 binRpts = 6; 
	 if(sevents.ptV1<70 || sevents.ptV2<70){
	   binRpts = 5; 
	   if(sevents.ptV1<50 || sevents.ptV2<50){
	     binRpts = 4; 
	     if(sevents.ptV1<40 || sevents.ptV2<40){
	       binRpts = 3; 
	       if(sevents.ptV1<30 || sevents.ptV2<30){
		 binRpts = 2; 
		 if(sevents.ptV1<20 || sevents.ptV2<20){
		   binRpts = 1; 
		   if(sevents.ptV1<10 || sevents.ptV2<10){
		     binRpts = 0; 
		   }
		 }
	       }
	     }
	   }
	 }
       }
       unsigned int binRpth = 7; 
       if(sevents.ptV1<100 && sevents.ptV2<100){
	 binRpth = 6; 
	 if(sevents.ptV1<70 && sevents.ptV2<70){
	   binRpth = 5; 
	   if(sevents.ptV1<50 && sevents.ptV2<50){
	     binRpth = 4; 
	     if(sevents.ptV1<40 && sevents.ptV2<40){
	       binRpth = 3; 
	       if(sevents.ptV1<30 && sevents.ptV2<30){
		 binRpth = 2; 
		 if(sevents.ptV1<20 && sevents.ptV2<20){
		   binRpth = 1; 
		   if(sevents.ptV1<10 && sevents.ptV2<10){
		     binRpth = 0; 
		   }
		 }
	       }
	     }
	   }
	 }
       }
       unsigned int binRptJET = 7; 
       if(ptHardestPJ<100){
	 binRptJET = 6; 
	 if(ptHardestPJ<70){
	   binRptJET = 5; 
	   if(ptHardestPJ<50){
	     binRptJET = 4; 
	     if(ptHardestPJ<40){
	       binRptJET = 3; 
	       if(ptHardestPJ<30){
		 binRptJET = 2; 
		 if(ptHardestPJ<20){
		   binRptJET = 1; 
		   if(ptHardestPJ<10){
		     binRptJET = 0; 
		   }
		 }
	       }
	     }
	   }
	 }
       }
       
//        unsigned int binRpts = 2; 
//        if(sevents.ptV1<20 || sevents.ptV2<20){
// 	 binRpts = 1; 
// 	 if(sevents.ptV1<10 || sevents.ptV2<10){
// 	   binRpts = 0; 
// 	 }
//        }
//        unsigned int binRpth = 2; 
//        if(sevents.ptV1<20 && sevents.ptV2<20){
// 	 binRpth = 1; 
// 	 if(sevents.ptV1<10 && sevents.ptV2<10){
// 	   binRpth = 0; 
// 	 }
//        }  
//        std::cout << "fill " << binSDR << binSpts << binSpth << binRDR << binRpts << binRpth << std::endl;
       count[binSDR][binSpts][binSpth][binRDR][binRpts][binRpth]++;
       countJet[binSDR][binSptJET][binRDR][binRptJET]++;
     }//end VERTICES
   }

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
BCorrAnalyzer::beginJob()
{
  //  TFileDirectory sdtEvents = fs_->mkdir("tevents","tree for events");
  //TFileDirectory sdpEvents = fs_->mkdir("pevents","plots for events");
  //TFileDirectory sdtBCands = fs_->mkdir("tbcandidates","tree for bcandidates");

  tEvents = /*(&sdtEvents)*/ fs_->make<TTree>("tEvents","event properties");
  tEvents->Branch("BEvents", &sevents.nB, "nB/I:nV/I:nMat/I:process/I:eventProcId/I:dRvv/F:dEtavv/F:dPhivv/F:dRbb/F:dEtabb/F:dPhibb/F:massV1/F:massV2/F:massV3/F:ptV1/F:ptV2/F:ptV3/F:etaV1/F:etaV2/F:etaV3/F:phiV1/F:phiV2/F:phiV3/F:dist3D1/F:dist3D2/F:dist3D3/F:distSig3D1/F:distSig3D2/F:distSig3D3/F:dist2D1/F:dist2D2/F:dist2D3/F:distSig2D1/F:distSig2D2/F:distSig2D3/F:massB1/F:massB2/F:ptB1/F:ptB2/F:etaB1/F:etaB2/F:phiB1/F:phiB2/F:ptHardestGJ/F:ptHardestPJ/F:ptHardestPGJ/F:ptHardestGJ/F:etaHardestPJ/F:etaHardestPGJ/F:eventNo/F:runNo/F:pthat/F:flavors/I:jet6/I:jet10/I:jet10nobptx/I:jet15/I:jet15hcalnf/I:jet30/I:jet50/I", 128000);

  //run info
  tEvents->Branch("pthat",&pthat,"pthat/D");
  tEvents->Branch("eventNumber",&eventNumber,"eventNumber/I");
  tEvents->Branch("eventRunNumber",&eventRunNumber,"eventRunNumber/I");
  tEvents->Branch("eventLumiSection",&eventLumiSection,"eventLumiSection/I");

  //pdf info
  tEvents->Branch("pdf_Q",&pdf_Q,"pdf_Q/D");
  tEvents->Branch("pdf_x1",&pdf_x1,"pdf_x1/D");
  tEvents->Branch("pdf_x2",&pdf_x2,"pdf_x2/D");
  tEvents->Branch("pdf_pdf1",&pdf_pdf1,"pdf_pdf1/D");
  tEvents->Branch("pdf_pdf2",&pdf_pdf2,"pdf_pdf2/D");
  tEvents->Branch("pdf_id1",&pdf_id1,"pdf_id1/I");
  tEvents->Branch("pdf_id2",&pdf_id2,"pdf_id2/I");


  //trigger info
  for( map<string, int>::iterator t = varTrigger.begin(); t!=varTrigger.end(); t++){
    tEvents->Branch(( t->first ).c_str(), &(varTrigger[ t->first ]),(t->first+"/I").c_str());
  }

  //PFJets
  tEvents->Branch("PFJet_nJets",&PFJet_nJets,"PFJet_nJets/I");
  tEvents->Branch("PFJet_px",PFJet_px,"PFJet_px[PFJet_nJets]/F");
  tEvents->Branch("PFJet_py",PFJet_py,"PFJet_py[PFJet_nJets]/F");
  tEvents->Branch("PFJet_pz",PFJet_pz,"PFJet_pz[PFJet_nJets]/F");
  tEvents->Branch("PFJet_pt",PFJet_pt,"PFJet_pt[PFJet_nJets]/F");
  tEvents->Branch("PFJet_eta",PFJet_eta,"PFJet_eta[PFJet_nJets]/F");
  tEvents->Branch("PFJet_phi",PFJet_phi,"PFJet_phi[PFJet_nJets]/F");
  tEvents->Branch("PFJet_E",PFJet_E,"PFJet_E[PFJet_nJets]/F");
  tEvents->Branch("PFJet_JetIDLoose",PFJet_JetIDLoose,"PFJet_JetIDLoose[PFJet_nJets]/I");
  tEvents->Branch("PFJet_JetFlavour",PFJet_JetFlavour,"PFJet_JetFlavour[PFJet_nJets]/I");

  tEvents->Branch("PFJet_nHEF",PFJet_nHEF,"PFJet_nHEF[PFJet_nJets]/F");
  tEvents->Branch("PFJet_nEmEF",PFJet_nEmEF,"PFJet_nEmEF[PFJet_nJets]/F");
  tEvents->Branch("PFJet_cHEF",PFJet_cHEF,"PFJet_cHEF[PFJet_nJets]/F");
  tEvents->Branch("PFJet_cEmEF",PFJet_cEmEF,"PFJet_cEmEF[PFJet_nJets]/F");
  tEvents->Branch("PFJet_cMult",PFJet_cMult,"PFJet_cMult[PFJet_nJets]/I");
  tEvents->Branch("PFJet_nConst",PFJet_nConst,"PFJet_nConst[PFJet_nJets]/I");


  for(vector<string>::iterator balg = BALGO.begin(); balg!=BALGO.end(); balg++){
    string balgo= (*balg);
    tEvents->Branch(("PFJet_btag_"+balgo).c_str(), varFloatArr["PFJet_btag_"+balgo],("PFJet_btag_"+balgo+"[PFJet_nJets]/F").c_str());
  }

  for(vector<string>::iterator wp = WP.begin(); wp!=WP.end(); wp++){
    string wpn = (*wp);
    tEvents->Branch(("PFJet_btag_"+wpn+"_nsel").c_str(),&(varInt["PFJet_btag_"+wpn+"_nsel"]),("PFJet_btag_"+wpn+"_nsel/I").c_str() );
    tEvents->Branch(("PFJet_btag_"+wpn+"_pt").c_str(),(varFloatArr["PFJet_btag_"+wpn+"_pt"]),("PFJet_btag_"+wpn+"_pt[PFJet_btag_"+wpn+"_nsel]/F").c_str() );
    tEvents->Branch(("PFJet_btag_"+wpn+"_eta").c_str(),(varFloatArr["PFJet_btag_"+wpn+"_eta"]),("PFJet_btag_"+wpn+"_eta[PFJet_btag_"+wpn+"_nsel]/F").c_str() );
    tEvents->Branch(("PFJet_btag_"+wpn+"_phi").c_str(),(varFloatArr["PFJet_btag_"+wpn+"_phi"]),("PFJet_btag_"+wpn+"_phi[PFJet_btag_"+wpn+"_nsel]/F").c_str() );
    tEvents->Branch(("PFJet_btag_"+wpn+"_discr").c_str(),(varFloatArr["PFJet_btag_"+wpn+"_discr"]),("PFJet_btag_"+wpn+"_discr[PFJet_btag_"+wpn+"_nsel]/F").c_str() );
    tEvents->Branch(("PFJet_btag_"+wpn+"_JetFlavour").c_str(),(varIntArr["PFJet_btag_"+wpn+"_JetFlavour"]),("PFJet_btag_"+wpn+"_JetFlavour[PFJet_btag_"+wpn+"_nsel]/I").c_str() );
  }

  tEvents->Branch("CaloJet_nJets",&CaloJet_nJets,"CaloJet_nJets/I");
  tEvents->Branch("CaloJet_px",CaloJet_px,"CaloJet_px[CaloJet_nJets]/F");
  tEvents->Branch("CaloJet_py",CaloJet_py,"CaloJet_py[CaloJet_nJets]/F");
  tEvents->Branch("CaloJet_pz",CaloJet_pz,"CaloJet_pz[CaloJet_nJets]/F");
  tEvents->Branch("CaloJet_pt",CaloJet_pt,"CaloJet_pt[CaloJet_nJets]/F");
  tEvents->Branch("CaloJet_eta",CaloJet_eta,"CaloJet_eta[CaloJet_nJets]/F");
  tEvents->Branch("CaloJet_phi",CaloJet_phi,"CaloJet_phi[CaloJet_nJets]/F");
  tEvents->Branch("CaloJet_E",CaloJet_E,"CaloJet_E[CaloJet_nJets]/F");
  tEvents->Branch("CaloJet_JetIDLoose",CaloJet_JetIDLoose,"CaloJet_JetIDLoose[CaloJet_nJets]/I");
  tEvents->Branch("CaloJet_JetIDTight",CaloJet_JetIDTight,"CaloJet_JetIDTight[CaloJet_nJets]/I");
  tEvents->Branch("CaloJet_JetFlavour",CaloJet_JetFlavour,"CaloJet_JetFlavour[CaloJet_nJets]/I");
  tEvents->Branch("CaloJet_emEF",CaloJet_emEF,"CaloJet_emEF[CaloJet_nJets]/F");
  tEvents->Branch("CaloJet_fHPD",CaloJet_fHPD,"CaloJet_fHPD[CaloJet_nJets]/F");
  tEvents->Branch("CaloJet_n90Hits",CaloJet_n90Hits,"CaloJet_n90Hits[CaloJet_nJets]/F");
  tEvents->Branch("CaloJet_fRBX",CaloJet_fRBX,"CaloJet_fRBX[CaloJet_nJets]/F");
  
  for(vector<string>::iterator balg = BALGO.begin(); balg!=BALGO.end(); balg++){
    string balgo= (*balg);
    tEvents->Branch(("CaloJet_btag_"+balgo).c_str(), varFloatArr["CaloJet_btag_"+balgo],("CaloJet_btag_"+balgo+"[CaloJet_nJets]/F").c_str());
  }
  
  for(vector<string>::iterator wp = WP.begin(); wp!=WP.end(); wp++){
    string wpn = (*wp);
    tEvents->Branch(("CaloJet_btag_"+wpn+"_nsel").c_str(),&(varInt["CaloJet_btag_"+wpn+"_nsel"]),("CaloJet_btag_"+wpn+"_nsel/I").c_str() );
    tEvents->Branch(("CaloJet_btag_"+wpn+"_pt").c_str(),(varFloatArr["CaloJet_btag_"+wpn+"_pt"]),("CaloJet_btag_"+wpn+"_pt[CaloJet_btag_"+wpn+"_nsel]/F").c_str() );
    tEvents->Branch(("CaloJet_btag_"+wpn+"_eta").c_str(),(varFloatArr["CaloJet_btag_"+wpn+"_eta"]),("CaloJet_btag_"+wpn+"_eta[CaloJet_btag_"+wpn+"_nsel]/F").c_str() );
    tEvents->Branch(("CaloJet_btag_"+wpn+"_phi").c_str(),(varFloatArr["CaloJet_btag_"+wpn+"_phi"]),("CaloJet_btag_"+wpn+"_phi[CaloJet_btag_"+wpn+"_nsel]/F").c_str() );
    tEvents->Branch(("CaloJet_btag_"+wpn+"_discr").c_str(),(varFloatArr["CaloJet_btag_"+wpn+"_discr"]),("CaloJet_btag_"+wpn+"_discr[CaloJet_btag_"+wpn+"_nsel]/F").c_str() );
    tEvents->Branch(("CaloJet_btag_"+wpn+"_JetFlavour").c_str(),(varIntArr["CaloJet_btag_"+wpn+"_JetFlavour"]),("CaloJet_btag_"+wpn+"_JetFlavour[CaloJet_btag_"+wpn+"_nsel]/I").c_str() );
  }




  tBCands = /*(&sdtBCands)*/ fs_->make<TTree>("tBCands","bcandidate properties");
  tBCands->Branch("BCands", &sbcands.massBcand, "massBcand/F:massVert/F:gamma/F:pt/F:eta/F:phi/F:dist3D/F:distSig3D/F:dist2D/F:distSig2D/F:dist3D_norm/F:eventNo/F:runNo/F:pthat/F:nvert/I:nSelected/I:ndRmatched/I:flavor/I:jet6/I:jet10/I:jet10nobptx/I:jet15/I:jet15hcalnf/I:jet30/I:jet50/I:process/I:eventProcId/I:selected/O:dRmatched/O", 128000);


  hdR = /*(&sdpEvents)*/fs_->make<TH1F>("ivf_hdR","hdR",20,0,10);
  hdEta = /*(&sdpEvents)*/fs_->make<TH1F>("ivf_hdEta","hdEta",20,-5,5);
  hdPhi = /*(&sdpEvents)*/fs_->make<TH1F>("ivf_hdPhi","hdPhi",20,-4,4);
  
  for(unsigned int i1=0; i1<8; i1++){
    for(unsigned int i2=0; i2<8; i2++){
      normJet[i1][i2] = 0; 
      for(unsigned int i3=0; i3<8; i3++){
	norm[i1][i2][i3] = 0;
	for(unsigned int i4=0; i4<8; i4++){
	  countJet[i1][i2][i3][i4] = 0;
	  for(unsigned int i5=0; i5<8; i5++){
	    for(unsigned int i6=0; i6<8; i6++){
	      count[i1][i2][i3][i4][i5][i6] = 0; 
	    }
	  }
	}
      }
    }
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
BCorrAnalyzer::endJob() {
  /*
  std::cout << "efficMatBV1 "; 
  for(unsigned int i1=0; i1<8; i1++){
    for(unsigned int i2=0; i2<8; i2++){
      for(unsigned int i3=0; i3<8; i3++){
	for(unsigned int i4=0; i4<8; i4++){
	  for(unsigned int i5=0; i5<8; i5++){
	    for(unsigned int i6=0; i6<8; i6++){
	      std::cout <<  norm[i1][i2][i3] << " ";

	    }
	  }
	}
      }
    }
  }
  
  std::cout << std::endl;
  std::cout << "efficMatBV2 "; 

  for(unsigned int i1=0; i1<8; i1++){
    for(unsigned int i2=0; i2<8; i2++){
      for(unsigned int i3=0; i3<8; i3++){
// 	std::cout << "NORM " << norm[i1][i2][i3] << "\n";
	for(unsigned int i4=0; i4<8; i4++){
	  for(unsigned int i5=0; i5<8; i5++){
	    for(unsigned int i6=0; i6<8; i6++){
// 	      if(count[i1][i2][i3][i4][i5][i6]!=0) std::cout << i1 << i2 << i3 << i4 << i5 << i6 << "    " << count[i1][i2][i3][i4][i5][i6] << std::endl;
	      std::cout << count[i1][i2][i3][i4][i5][i6] << " "; 
	    }
	  }
	}
      }
    }
  }
  std::cout << std::endl;

  std::cout << "efficMatJET1 "; 
  for(unsigned int i1=0; i1<8; i1++){
    for(unsigned int i2=0; i2<8; i2++){
      for(unsigned int i3=0; i3<8; i3++){
	for(unsigned int i4=0; i4<8; i4++){
	  std::cout << normJet[i1][i2] << " ";
	}
      }
    }
  }
  std::cout << "\nefficMatJET2 "; 
  for(unsigned int i1=0; i1<8; i1++){
    for(unsigned int i2=0; i2<8; i2++){
      for(unsigned int i3=0; i3<8; i3++){
	for(unsigned int i4=0; i4<8; i4++){
	  std::cout << countJet[i1][i2][i3][i4] << " ";
	}
      }
    }
  }
  */
  std::cout << std::endl;
  std::cout << std::endl;

}

GlobalVector
BCorrAnalyzer::flightDirection(reco::Vertex &pv, reco::Vertex &sv){
  GlobalVector res(sv.position().X() - pv.position().X(),
                    sv.position().Y() - pv.position().Y(),
                    sv.position().Z() - pv.position().Z());
  return res;
}

//define this as a plug-in
DEFINE_FWK_MODULE(BCorrAnalyzer);
