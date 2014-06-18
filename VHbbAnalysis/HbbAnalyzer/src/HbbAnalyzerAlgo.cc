// -*- C++ -*-
//
// Package:    HbbAnalyzerAlgo
// Class:      HbbAnalyzerAlgo
// 
/**\class HbbAnalyzerAlgo HbbAnalyzerAlgo.cc Analysis/HbbAnalyzer/src/HbbAnalyzerAlgo.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  David Lopes Pegna,Address unknown,NONE,
//         Created:  Thu Mar  5 13:51:28 EST 2009
// $Id: HbbAnalyzerAlgo.cc,v 1.79 2013/01/25 01:56:57 jiafulow Exp $
//
//


//uncomment to save also jet collections 1 and 4
//#define ENABLE_SIMPLEJETS1
//#define ENABLE_SIMPLEJETS4

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h" 
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "VHbbAnalysis/HbbAnalyzer/interface/HbbAnalyzerAlgo.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "Math/GenVector/PxPyPzM4D.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

//TOM #include "CMGTools/External/interface/PileupJetIdentifier.h"
//TOM #include "CMGTools/External/interface/PileupJetIdAlgo.h"
//#include "CMGTools/External/interface/PileupJetIdProducer.h"

#include <cmath>




#define GENPTOLOR(a) TLorentzVector((a).px(), (a).py(), (a).pz(), (a).energy())
#define GENPTOLORP(a) TLorentzVector((a)->px(), (a)->py(), (a)->pz(), (a)->energy())



struct CompareJetPtMuons {
  bool operator()( const VHbbEvent::MuonInfo& j1, const  VHbbEvent::MuonInfo& j2 ) const {
    return j1.p4.Pt() > j2.p4.Pt();
  }
};
struct CompareJetPtElectrons {
  bool operator()( const VHbbEvent::ElectronInfo& j1, const  VHbbEvent::ElectronInfo& j2 ) const {
    return j1.p4.Pt() > j2.p4.Pt();
  }
};
struct CompareJetPtTaus {
  bool operator()( const VHbbEvent::TauInfo& j1, const  VHbbEvent::TauInfo& j2 ) const {
    return j1.p4.Pt() > j2.p4.Pt();
  }
};



HbbAnalyzerAlgo::HbbAnalyzerAlgo(){
  /*process.HbbAnalyzerNew = cms.EDProducer("HbbAnalyzerNew",
    runOnMC = cms.bool(isMC),
    hltResultsTag = cms.InputTag("TriggerResults::HLT"),
    lep_ptCutForBjets = cms.double(5),
    electronNoCutsTag = cms.InputTag("gsfElectrons"),
#    electronTag = cms.InputTag("selectedElectronsMatched"),
    electronTag = cms.InputTag("selectedPatElectrons"),
    tauTag = cms.InputTag("patTaus"),

    muonNoCutsTag = cms.InputTag("muons"),
   muonTag = cms.InputTag("selectedPatMuons"),
#    muonTag = cms.InputTag("selectedMuonsMatched"),
    jetTag = cms.InputTag("selectedPatJetsCAVHFatPF"),
    subjetTag = cms.InputTag("selectedPatJetsCAVHSubPF"),
    filterjetTag = cms.InputTag("selectedPatJetsCAVHFilterPF"),
    simplejet2Tag = cms.InputTag("cleanPatJets"),
    simplejet3Tag = cms.InputTag("selectedPatJetsAK7PF"),
    photonTag = cms.InputTag("selectedPatPhotons"),
    metTag = cms.InputTag("met"), #this input tag is used to fill calo MET 
    verbose = cms.untracked.bool(False),

    #TODO: clean up the analyzer
    simplejet1Tag = cms.InputTag("UNUSED_WAS_selectedPatJets"),
    simplejet4Tag = cms.InputTag("UNUSED_WAS_selectedPatJetsAK7Calo"),

)
  */

  verbose_ = false;
  eleLabel_ = edm::InputTag("slimmedElectrons","","");
  muoLabel_ = edm::InputTag("slimmedMuons","","");
  lep_ptCutForBjets_ = 5;
  // TOM remove
  //  elenoCutsLabel_(iConfig.getParameter<edm::InputTag>("electro")),
  //  muonoCutsLabel_(iConfig.getParameter<edm::InputTag>("muonNoCutsTag")),
  jetLabel_ = edm::InputTag("slimmedJets","","");
  //  subjetLabel_(iConfig.getParameter<edm::InputTag>("subjetTag")),
  //  filterjetLabel_(iConfig.getParameter<edm::InputTag>("filterjetTag")),
  // simplejet1Label_(iConfig.getParameter<edm::InputTag>("simplejet1Tag")),
  simplejet2Label_=  edm::InputTag("slimmedJets","","");
  //simplejet3Label_(iConfig.getParameter<edm::InputTag>("simplejet3Tag")),
  //simplejet4Label_(iConfig.getParameter<edm::InputTag>("simplejet4Tag")),
  tauLabel_ = edm::InputTag("slimmedTaus","","");
  metLabel_ = edm::InputTag("slimmedMETs","","");
  phoLabel_ = edm::InputTag("slimmedPhotons");
  hltResults_ = edm::InputTag("TriggerResults","","HLT");
  runOnMC_ = true;


  /*
  eleLabel_(iConfig.getParameter<edm::InputTag>("electronTag")),

  muoLabel_(iConfig.getParameter<edm::InputTag>("muonTag")),
  lep_ptCutForBjets_(iConfig.getParameter<double>("lep_ptCutForBjets")),
  elenoCutsLabel_(iConfig.getParameter<edm::InputTag>("electronNoCutsTag")),
  muonoCutsLabel_(iConfig.getParameter<edm::InputTag>("muonNoCutsTag")),
  jetLabel_(iConfig.getParameter<edm::InputTag>("jetTag")),
  subjetLabel_(iConfig.getParameter<edm::InputTag>("subjetTag")),
  filterjetLabel_(iConfig.getParameter<edm::InputTag>("filterjetTag")),
  simplejet1Label_(iConfig.getParameter<edm::InputTag>("simplejet1Tag")),
  simplejet2Label_(iConfig.getParameter<edm::InputTag>("simplejet2Tag")),
  simplejet3Label_(iConfig.getParameter<edm::InputTag>("simplejet3Tag")),
  simplejet4Label_(iConfig.getParameter<edm::InputTag>("simplejet4Tag")),
  tauLabel_(iConfig.getParameter<edm::InputTag>("tauTag")),
  metLabel_(iConfig.getParameter<edm::InputTag>("metTag")),
  phoLabel_(iConfig.getParameter<edm::InputTag>("photonTag")),
  hltResults_(iConfig.getParameter<edm::InputTag>("hltResultsTag")),
  runOnMC_(iConfig.getParameter<bool>("runOnMC")), verbose_(iConfig.getUntrackedParameter<bool>("verbose")) {
  */
  //
  // put the setwhatproduced etc etc
 
  //  produces<VHbbEvent>();
  //  produces<VHbbEventAuxInfo>();

 
}


HbbAnalyzerAlgo::~HbbAnalyzerAlgo(){
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HbbAnalyzerAlgo::produce( fwlite::Event & ev,VHbbEvent& vhbbevent_, VHbbEventAuxInfo& vhbbeventauxinfo_ ){
  using namespace edm;
  using namespace reco;
  

  
  VHbbEvent* hbbInfo= &vhbbevent_;
  VHbbEventAuxInfo* auxInfo = &vhbbeventauxinfo_;
  fwlite::Event& iEvent = ev;

  if (runOnMC_){
    fwlite::Handle<GenEventInfoProduct> evt_info;
    evt_info.getByLabel(ev,"generator");
    auxInfo->weightMCProd = evt_info->weight(); 
  }
  else
    { auxInfo->weightMCProd =1.;}
  //
  // ??
  
  // trigger

  // trigger
  fwlite::Handle<edm::TriggerResults>  hltresults;
  //iEvent.getByLabel("TriggerResults", hltresults);
   
  //edm::InputTag tag("TriggerResults::HLT");
  //  edm::InputTag tag("TriggerResults::HLT0");
  hltresults.getByLabel(iEvent, hltResults_.label().c_str());
   
  const edm::TriggerNames & triggerNames_ = iEvent.triggerNames(*hltresults);

  int ntrigs = hltresults->size();
  if (ntrigs==0){std::cerr << "%HLTInfo -- No trigger name given in TriggerResults of the input " << std::endl;}

  BeamSpot vertexBeamSpot;
  fwlite::Handle<reco::BeamSpot> recoBeamSpotHandle;
  recoBeamSpotHandle.getByLabel(iEvent,"offlineBeamSpot","");
  vertexBeamSpot = *recoBeamSpotHandle;
  /*
    double BSx=vertexBeamSpot.x0();
    double BSy=vertexBeamSpot.y0();
    double BSz=vertexBeamSpot.z0();
  */

  double MinVtxProb=-999.;
  int VtxIn=-99;
  
  fwlite::Handle<reco::VertexCollection> recVtxs;
  recVtxs.getByLabel(iEvent,"offlineSlimmedPrimaryVertices");
  
  auxInfo->pvInfo.nVertices = recVtxs->size();

  for(size_t i = 0; i < recVtxs->size(); ++ i) {
    const Vertex &vtx = (*recVtxs)[i];
    double RecVtxProb=TMath::Prob(vtx.chi2(),vtx.ndof());
    if(RecVtxProb>MinVtxProb){
      VtxIn=i;
      MinVtxProb=RecVtxProb;      
    }
  }

  const Vertex &RecVtx = (*recVtxs)[VtxIn];
  const Vertex &RecVtxFirst = (*recVtxs)[0];
  const Vertex &vertex = RecVtxFirst; //used in ele id 2012
  
  auxInfo->pvInfo.firstPVInPT2 = TVector3(RecVtxFirst.x(), RecVtxFirst.y(), RecVtxFirst.z());
  auxInfo->pvInfo.firstPVInProb = TVector3(RecVtx.x(), RecVtx.y(), RecVtx.z());
  
  (auxInfo->pvInfo).efirstPVInPT2 = (RecVtxFirst.error());
  (auxInfo->pvInfo).efirstPVInProb = RecVtx.error();
   
  fwlite::Handle<double> rhoHandle;
  rhoHandle.getByLabel(iEvent,"fixedGridRhoFastjetAll");
  auxInfo->puInfo.rho = *rhoHandle;
  //  std::cout <<" RHO "<< auxInfo->puInfo.rho<<std::endl;

  //AR change wrt 53X 
  fwlite::Handle<double> rho25Handle;
  rho25Handle.getByLabel(iEvent,"fixedGridRhoFastjetCentralChargedPileUp");  // was kt6PFJets25 
  auxInfo->puInfo.rho25 = *rho25Handle;
  fwlite::Handle<double> rho25HandleIso;
  rho25HandleIso.getByLabel(iEvent,"fixedGridRhoFastjetAll"); //was kt6PFJetsForIsolation", "rho");   
  auxInfo->puInfo.rho25Iso = *rho25HandleIso;
 
  fwlite::Handle<double> rhoNeutralHandle;
  rhoNeutralHandle.getByLabel(iEvent,"fixedGridRhoFastjetCentralNeutral"); // was kt6PFJetsCentralNeutral", "rho");   
  auxInfo->puInfo.rhoNeutral = *rhoNeutralHandle;

  fwlite::Handle<std::vector< PileupSummaryInfo> > puHandle;

  if (runOnMC_){
    puHandle.getByLabel(iEvent,"");
    if (puHandle.isValid()){
      
      std::vector< PileupSummaryInfo> pu = (*puHandle); 
      for (std::vector<PileupSummaryInfo>::const_iterator it= pu.begin(); it!=pu.end(); ++it){
	 int bx = (*it).getBunchCrossing();
         if(bx == 0) { auxInfo->puInfo.truePU = (*it).getTrueNumInteractions();}
	unsigned int num = (*it).getPU_NumInteractions();
	//	std::cout <<" PU PUSHING "<<bx<<" " <<num<<std::endl;
	auxInfo->puInfo.pus[bx]  =num;
      }
    }
  }

  //// real start
  
  
  fwlite::Handle<GenParticleCollection> genParticles;
  
  //  bool printJet=0;
  
  
  if(runOnMC_){
    
    genParticles.getByLabel(iEvent,"prunedGenParticles");
    
    for(size_t i = 0; i < genParticles->size(); ++ i) {
     
      const GenParticle & p = (*genParticles)[i];
      int id = p.pdgId();
      int st = p.status();  
      
      if(id==25){
	
	VHbbEventAuxInfo::ParticleMCInfo htemp;
	htemp.status=st;
	htemp.charge=p.charge();
	if(p.mother(0)!=0) htemp.momid=p.mother(0)->pdgId();
	if(p.mother(0)!=0 && p.mother(0)->mother(0)!=0) htemp.gmomid=p.mother(0)->mother(0)->pdgId(); 
	htemp.p4 = GENPTOLOR(p);
	
	int ndau = p.numberOfDaughters();
	for(int j = 0; j < ndau; ++ j) {
	  const Candidate * Hdau = p.daughter( j );
	  htemp.dauid.push_back(Hdau->pdgId());
	  htemp.dauFourMomentum.push_back(GENPTOLORP(Hdau));
	}
	(auxInfo->mcH).push_back(htemp);
      }
      
      
      if(abs(id)==24){
	
	VHbbEventAuxInfo::ParticleMCInfo wtemp;
	wtemp.status=st;
	wtemp.charge=p.charge();
	if(p.mother(0)!=0) wtemp.momid=p.mother(0)->pdgId();
	if(p.mother(0)!=0 && p.mother(0)->mother(0)!=0) wtemp.gmomid=p.mother(0)->mother(0)->pdgId();
	wtemp.p4=GENPTOLOR(p);
	
	int ndau = p.numberOfDaughters();
	for(int j = 0; j < ndau; ++ j) {
	  const Candidate * Wdau = p.daughter( j );
	  wtemp.dauid.push_back(Wdau->pdgId());
	  wtemp.dauFourMomentum.push_back(GENPTOLORP(Wdau));
	}
	auxInfo->mcW.push_back(wtemp);
      }

      if(abs(id)==15) {
	VHbbEventAuxInfo::ParticleMCInfo tautemp;
        tautemp.status=st;
        tautemp.charge=p.charge();
        if(p.mother(0)!=0) tautemp.momid=p.mother(0)->pdgId();
        if(p.mother(0)!=0 && p.mother(0)->mother(0)!=0) tautemp.gmomid=p.mother(0)->mother(0)->pdgId();
        tautemp.p4=GENPTOLOR(p);

        int ndau = p.numberOfDaughters();
        for(int j = 0; j < ndau; ++ j) {
          const Candidate * Taudau = p.daughter( j );
          tautemp.dauid.push_back(Taudau->pdgId());
          tautemp.dauFourMomentum.push_back(GENPTOLORP(Taudau));
        }
        auxInfo->mcTau.push_back(tautemp);
      }
	
      if(abs(id)==23){


	VHbbEventAuxInfo::ParticleMCInfo ztemp;
	ztemp.status=st;
	ztemp.charge=p.charge();
	if(p.mother(0)!=0) ztemp.momid=p.mother(0)->pdgId();
	if(p.mother(0)!=0 && p.mother(0)->mother(0)!=0) ztemp.gmomid=p.mother(0)->mother(0)->pdgId();
	ztemp.p4=GENPTOLOR(p);

	int ndau = p.numberOfDaughters();
	for(int j = 0; j < ndau; ++ j) {
	  const Candidate * Zdau = p.daughter( j );
	  ztemp.dauid.push_back(Zdau->pdgId());
	  ztemp.dauFourMomentum.push_back(GENPTOLORP(Zdau));
	}
	auxInfo->mcZ.push_back(ztemp);
      }
      //
      // binfo
      //
      
      
      if(id==5){

	VHbbEventAuxInfo::ParticleMCInfo btemp;
	btemp.status=st;
	btemp.charge=p.charge();
	if(p.mother(0)!=0) btemp.momid=p.mother(0)->pdgId();
	if(p.mother(0)!=0 && p.mother(0)->mother(0)!=0) btemp.gmomid=p.mother(0)->mother(0)->pdgId(); 

	btemp.p4=GENPTOLOR(p);
	
	int nHDaubdau = p.numberOfDaughters();
	for(int j = 0; j < nHDaubdau; ++ j) {
	  const Candidate * Bdau = p.daughter( j );
	  btemp.dauid.push_back(Bdau->pdgId());
	}
	auxInfo->mcB.push_back(btemp);
      }
      
      if(id==-5){

	VHbbEventAuxInfo::ParticleMCInfo bbtemp;
	
	bbtemp.status=st;
	bbtemp.charge=p.charge();
	if(p.mother(0)!=0) bbtemp.momid=p.mother(0)->pdgId();
	if(p.mother(0)!=0 && p.mother(0)->mother(0)!=0) bbtemp.gmomid=p.mother(0)->mother(0)->pdgId(); 

	bbtemp.p4=GENPTOLOR(p);
	
	int nHDaubdau = p.numberOfDaughters();
	for(int j = 0; j < nHDaubdau; ++ j) {
	  const Candidate * Bdau = p.daughter( j );
	  bbtemp.dauid.push_back(Bdau->pdgId());
	}


	auxInfo->mcBbar.push_back(bbtemp);
     }
      
      if(abs(id)==4){
	VHbbEventAuxInfo::ParticleMCInfo ctemp;
	ctemp.status=st;
	ctemp.charge=p.charge();
	if(p.mother(0)!=0) ctemp.momid=p.mother(0)->pdgId();
	if(p.mother(0)!=0 && p.mother(0)->mother(0)!=0) ctemp.gmomid=p.mother(0)->mother(0)->pdgId();

	ctemp.p4=GENPTOLOR(p);
	
	int nHDaubdau = p.numberOfDaughters();
	for(int j = 0; j < nHDaubdau; ++ j) {
	  const Candidate * Bdau = p.daughter( j );
	  ctemp.dauid.push_back(Bdau->pdgId());
	}

	auxInfo->mcC.push_back(ctemp);	

      }

    }

  }   // isMC

  /////// end generator block    

/// photon used in isolation
  /* TOM
 * AR: use packedPfCands 	
 *
  fwlite::Handle<edm::View<reco::PFCandidate> > photonIsoH;
  photonIsoH.getByLabel(iEvent,"pfAllPhotons");
  edm::View<reco::PFCandidate> photonsForIso = *photonIsoH;
 */
  fwlite::Handle<std::vector<pat::Muon> > muonHandle;
  muonHandle.getByLabel(iEvent,muoLabel_.label().c_str());
  std::vector<pat::Muon> muons = *muonHandle;

  // hard jet   
  fwlite::Handle<std::vector<pat::Jet> > fatjetHandle;
  fatjetHandle.getByLabel(iEvent,jetLabel_.label().c_str());
  std::vector<pat::Jet> fatjets = *fatjetHandle;
  /* TOM

  // sub jet   
  fwlite::Handle<std::vector<pat::Jet> > subjetHandle;
  subjetHandle.getByLabel(iEvent,subjetLabel_.label().c_str());
  std::vector<pat::Jet> subjets = *subjetHandle;

  // filter jet   
  fwlite::Handle<std::vector<pat::Jet> > filterjetHandle;
  filterjetHandle.getByLabel(iEvent,filterjetLabel_.label().c_str());
  std::vector<pat::Jet> filterjets = *filterjetHandle;

  // standard jets

  */
  fwlite::Handle<std::vector<pat::Jet> > simplejet2Handle;
  simplejet2Handle.getByLabel(iEvent,simplejet2Label_.label().c_str());
  std::vector<pat::Jet> simplejets2 = *simplejet2Handle;
  /* TOM
  fwlite::Handle<std::vector<pat::Jet> > simplejet3Handle;
  simplejet3Handle.getByLabel(iEvent,simplejet3Label_.label().c_str());
  std::vector<pat::Jet> simplejets3 = *simplejet3Handle;

  */

  fwlite::Handle<std::vector<pat::Electron> > electronHandle;
  electronHandle.getByLabel(iEvent,eleLabel_.label().c_str());
  std::vector<pat::Electron> electrons = *electronHandle;


  //   edm::Handle<std::vector<pat::Photon> > phoHandle;
  //   iEvent.getByLabel(phoLabel_,phoHandle);
  //   std::vector<pat::Photon> photons = *phoHandle;

  fwlite::Handle<std::vector<pat::Tau> > tauHandle;
  tauHandle.getByLabel(iEvent,tauLabel_.label().c_str());
  std::vector<pat::Tau> taus = *tauHandle;
  
  /* TOM - CANNOT WORK HERE

  //Get the computer for the CSV
  ESHandle<JetTagComputer> handle;
  iSetup.get<JetTagComputerRecord>().get("combinedSecondaryVertex", handle);
  computer = dynamic_cast<const GenericMVAJetTagComputer*>(handle.product());

  //BTAGGING SCALE FACTOR FROM DATABASE
  //Combined Secondary Vertex Loose
  edm::ESHandle<BtagPerformance> bTagSF_CSVL_;
  iSetup.get<BTagPerformanceRecord>().get("BTAGCSVL",bTagSF_CSVL_);
  //Combined Secondary Vertex Medium
  edm::ESHandle<BtagPerformance> bTagSF_CSVM_;
  iSetup.get<BTagPerformanceRecord>().get("BTAGCSVM",bTagSF_CSVM_);
  //Combined Secondary Vertex Tight
  edm::ESHandle<BtagPerformance> bTagSF_CSVT_;
  iSetup.get<BTagPerformanceRecord>().get("BTAGCSVT",bTagSF_CSVT_);

  edm::ESHandle<BtagPerformance> mistagSF_CSVL_;
  iSetup.get<BTagPerformanceRecord>().get("MISTAGCSVL",mistagSF_CSVL_);
  //Combined Secondary Vertex Medium
  edm::ESHandle<BtagPerformance> mistagSF_CSVM_;
  iSetup.get<BTagPerformanceRecord>().get("MISTAGCSVM",mistagSF_CSVM_);
  //Combined Secondary Vertex Tight
  edm::ESHandle<BtagPerformance> mistagSF_CSVT_;
  iSetup.get<BTagPerformanceRecord>().get("MISTAGCSVT",mistagSF_CSVT_);

  */

  //ccla moved out to SimpleJets2 loop
  // TOM
  /* fwlite::Handle<std::vector<reco::Candidate> >  muonNoCutsHandle;
  muonNoCutsHandle.getByLabel(iEvent,muonoCutsLabel_.label().c_str());
  std::vector<reco::Candidate> muonsNoCuts = *muonNoCutsHandle; 

  fwlite::Handle<std::vector<reco::Candidate> > eleNoCutsHandle;
  eleNoCutsHandle.getByLabel(iEvent,elenoCutsLabel_.label().c_str());
  std::vector<reco::Candidate> elesNoCuts = *eleNoCutsHandle; 
*/
  /* TOM
BTagSFContainer btagSFs;
  btagSFs.BTAGSF_CSVL = (bTagSF_CSVL_.product());
  btagSFs.BTAGSF_CSVM = (bTagSF_CSVM_.product());
  btagSFs.BTAGSF_CSVT = (bTagSF_CSVT_.product());
  btagSFs.MISTAGSF_CSVL = (mistagSF_CSVL_.product());
  btagSFs.MISTAGSF_CSVM = (mistagSF_CSVM_.product());
  btagSFs.MISTAGSF_CSVT = (mistagSF_CSVT_.product());

#ifdef ENABLE_SIMPLEJETS1
  fwlite::Handle<std::vector<pat::Jet> > simplejet1Handle;
  iEvent.getByLabel(simplejet1Label_,simplejet1Handle);
  std::vector<pat::Jet> simplejets1 = *simplejet1Handle;
  for(std::vector<pat::Jet>::const_iterator jet_iter = simplejets1.begin(); jet_iter!=simplejets1.end(); ++jet_iter){
    //     if(jet_iter->pt()>50)
    //       njetscounter++;
    VHbbEvent::SimpleJet sj;
    //    std::cout <<" sj1"<<std::endl;
    fillSimpleJet(sj,jet_iter);
    //    if(!runOnMC_) 

   setJecUnc(sj,jecUnc);

    Particle::LorentzVector p4Jet = jet_iter->p4();

    if(runOnMC_){

      fillScaleFactors(sj, btagSFs);

      //PAT genJet matching
      //genJet
      const reco::GenJet *gJ = jet_iter->genJet();
      //physical parton for mother info ONLY
      if( (jet_iter->genParton()) ){
	sj.bestMCid = jet_iter->genParton()->pdgId();
	if( (jet_iter->genParton()->mother()) )
	  sj.bestMCmomid=jet_iter->genParton()->mother()->pdgId();
      }
      TLorentzVector gJp4;
      if(gJ){
	gJp4.SetPtEtaPhiE(gJ->pt(),gJ->eta(),gJ->phi(),gJ->energy());
	sj.bestMCp4 = gJp4;
	if(verbose_){
	  std::clog << "genJet matched Pt = " << gJp4.Pt() << std::endl;
	  std::clog << "genJet matched eta = " << gJp4.Eta() << std::endl;
	  std::clog << "genJet matched deltaR = " <<gJp4.DeltaR(sj.p4) << std::endl;
	  std::clog << "genJet matched mother id = " << sj.bestMCmomid << std::endl;
	}
      }
      
    } //isMC
    hbbInfo->simpleJets.push_back(sj);
    
  }
#endif //ENABLE_SIMPLEJETS1


  */

  /* TOM
     
     for(std::vector<pat::Jet>::const_iterator jet_iter = simplejets3.begin(); jet_iter!=simplejets3.end(); ++jet_iter){
     //     if(jet_iter->pt()>50)
    //       njetscounter++;
    VHbbEvent::SimpleJet sj;
    //    std::cout <<" sj3"<<std::endl;
    fillSimpleJet(sj,jet_iter);
    //    if(!runOnMC_)  
    //  setJecUnc(sj,jecUnc);

   Particle::LorentzVector p4Jet = jet_iter->p4();

    if(runOnMC_){

      //  TOM 
      // TOM     fillScaleFactors(sj, btagSFs);

      //PAT genJet matching
      //genJet
      const reco::GenJet *gJ = jet_iter->genJet();
      //physical parton for mother info ONLY
      if( (jet_iter->genParton()) ){
	sj.bestMCid = jet_iter->genParton()->pdgId();
	if( (jet_iter->genParton()->mother()) )
	  sj.bestMCmomid=jet_iter->genParton()->mother()->pdgId();
      }
      TLorentzVector gJp4;
      if(gJ){
	gJp4.SetPtEtaPhiE(gJ->pt(),gJ->eta(),gJ->phi(),gJ->energy());
	sj.bestMCp4 = gJp4;
	if(verbose_){
	  std::clog << "genJet matched Pt = " << gJp4.Pt() << std::endl;
	  std::clog << "genJet matched eta = " << gJp4.Eta() << std::endl;
	  std::clog << "genJet matched deltaR = " <<gJp4.DeltaR(sj.p4) << std::endl;
	  std::clog << "genJet matched mother id = " << sj.bestMCmomid << std::endl;
	}
      }
      
    } //isMC
    // 


    hbbInfo->simpleJets3.push_back(sj);
    
 }
  */


  for(std::vector<pat::Jet>::const_iterator jet_iter = simplejets2.begin(); jet_iter!=simplejets2.end(); ++jet_iter){
    


    VHbbEvent::SimpleJet sj;
    

    fillSimpleJet(sj,jet_iter);  


    ///###########          PU JET ID #################
   // add puId...
/*    fwlite::Handle<edm::ValueMap<float> > puJetIdMVA;
    puJetIdMVA.getByLabel(iEvent,"puJetMva","fullDiscriminant");

    edm::Handle<edm::ValueMap<int> > puJetIdFlag;
    puJetIdFlag.getByLabel(iEvent,"puJetMva", "fullId");

    //    cout  << " pt " << jet_iter->pt() << " eta " << jet_iter->eta() << std::endl;
    unsigned int idx = jet_iter - simplejets2.begin();



    sj.puJetIdMva   = (*puJetIdMVA)[simplejets2.refAt(idx)];
    int    idflag = (*puJetIdFlag)[simplejets2.refAt(idx)];
  */
     sj.puJetIdMva = jet_iter->userFloat("pileupJetId:fullDiscriminant");

/* TOM 
 *
 * //     cout << " PU JetID MVA " << mva; 
    if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose )) {
      //cout << " pass loose wp";
      sj.puJetIdL =1;
	}
    if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kMedium )) {
      //    cout << " pass medium wp";
      sj.puJetIdM =1;
    }
    if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kTight )) {
      //    cout << " pass tight wp";
      sj.puJetIdT =1;
    }
    //    cout << endl;
    //  #############  END OF PU JET ID ######################

*/
    

//TOM  
//    if(!runOnMC_)  
//	sj.jecunc=
//          setJecUnc(sj,jecUnc);
   if(runOnMC_)
    sj.flavour = jet_iter->jetFlavourInfo().getPartonFlavour();
    
    
    sj.tche=jet_iter->bDiscriminator("trackCountingHighEffBJetTags");
    sj.tchp=jet_iter->bDiscriminator("trackCountingHighPurBJetTags");
    sj.jp=jet_iter->bDiscriminator("jetProbabilityBJetTags");
    sj.jpb=jet_iter->bDiscriminator("jetBProbabilityBJetTags");
    sj.ssvhe=jet_iter->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    sj.csv=jet_iter->bDiscriminator("combinedSecondaryVertexBJetTags");
    sj.csvmva=jet_iter->bDiscriminator("combinedSecondaryVertexMVABJetTags");
    sj.charge=jet_iter->jetCharge();
    sj.ntracks=jet_iter->associatedTracks().size();
    sj.p4=GENPTOLORP(jet_iter);
    sj.chargedTracksFourMomentum=(getChargedTracksMomentum(&*(jet_iter)));
    sj.SF_CSVL=1;
    sj.SF_CSVM=1;
    sj.SF_CSVT=1;
    sj.SF_CSVLerr=0;
    sj.SF_CSVMerr=0;
    sj.SF_CSVTerr=0;
   
 
    //
    // addtaginfo for csv
    //

    if (jet_iter->hasTagInfo("SimpleSecondaryVertex")) {

      const reco::SecondaryVertexTagInfo * tf = jet_iter->tagInfoSecondaryVertex();
      sj.vtxMass = tf->secondaryVertex(0).p4().mass();
      sj.vtxNTracks = tf->secondaryVertex(0).nTracks();
      Measurement1D m = tf->flightDistance(0);
      sj.vtx3dL = m.value();
      sj.vtx3deL = m.error();
    } else {
		sj.vtxMass=jet_iter->userFloat("vtxMass");
		sj.vtxNTracks=jet_iter->userFloat("vtxNtracks");
		sj.vtx3dL=jet_iter->userFloat("vtx3DVal");
		sj.vtx3deL= sj.vtx3dL/jet_iter->userFloat("vtx3DSig");
	}


    //
    // add tVector
    //
    sj.tVector = getTvect(&(*jet_iter));
   
    Particle::LorentzVector p4Jet = jet_iter->p4();

    if(runOnMC_){

      //BTV scale factors
      // TOM       fillScaleFactors(sj, btagSFs);

      //PAT genJet matching
      //genJet
      const reco::GenJet *gJ = jet_iter->genJet();
      //physical parton for mother info ONLY
      if( (jet_iter->genParton()) ){
	sj.bestMCid = jet_iter->genParton()->pdgId();
	if( (jet_iter->genParton()->mother()) )
	  sj.bestMCmomid=jet_iter->genParton()->mother()->pdgId();
      }
      TLorentzVector gJp4;
      if(gJ){
	gJp4.SetPtEtaPhiE(gJ->pt(),gJ->eta(),gJ->phi(),gJ->energy());
	sj.bestMCp4 = gJp4;
	if(verbose_){
	  std::clog << "genJet matched Pt = " << gJp4.Pt() << std::endl;
	  std::clog << "genJet matched eta = " << gJp4.Eta() << std::endl;
	  std::clog << "genJet matched deltaR = " << gJp4.DeltaR(sj.p4) << std::endl;
	  std::clog << "genJet matched mother id = " << sj.bestMCmomid << std::endl;
	}
      }

        // add flag if a mc lepton is find inside a cone around the jets... 
      genParticles.getByLabel(iEvent,"prunedGenParticles" );
      
      for(size_t i = 0; i < genParticles->size(); ++ i) {
     
      const GenParticle & p = (*genParticles)[i];
      int id = 0; 
      p.pt()> lep_ptCutForBjets_ ? id= p.pdgId(): 0;
   
      //      std::cout<< "found a muon with pt " << mu->pt()   << std::endl;
      if   ((abs(id)==13 || abs(id)==11) && deltaR(p.eta(), p.phi(), sj.p4.Eta(), sj.p4.Phi() ) <0.5)  sj.isSemiLeptMCtruth=1;
      }

    }  //isMC

        // add flag if a reco lepton is find inside a cone around the jets... 
    
    

    

    for(std::vector<pat::Muon>::const_iterator mu = muons.begin(); mu!=muons.end() && sj.isSemiLept!=1; ++mu){
      //      std::cout<< "found a muon with pt " << mu->pt()   << std::endl;
      const pat::Muon& m = *mu; 
      float Smpt = m.pt(); 
      float Smeta = m.eta();
      float Smphi = m.phi();
      
      float SmJdR = deltaR(Smeta, Smphi, sj.p4.Eta(), sj.p4.Phi());
      
      if   ( Smpt> lep_ptCutForBjets_ && SmJdR <0.5)  {
	sj.isSemiLept=1;
	//isSemiLept(-99), isSemiLeptMCtruth(-99), SoftLeptPt(-99), SoftLeptdR(-99), SoftLeptptRel(-99), SoftLeptpdgId(-99), SoftLeptIdlooseMu(-99), SoftLeptId95(-99), SoftLeptRelCombIso(-99),  
	sj.SoftLeptpdgId =13;
	sj.SoftLeptdR= SmJdR;
	sj.SoftLeptPt=Smpt;
	TVector3 mvec ( m.p4().Vect().X(), m.p4().Vect().Y(), m.p4().Vect().Z()  ); 
	sj.SoftLeptptRel=  sj.p4.Perp(  mvec );
	sj.SoftLeptRelCombIso = (m.trackIso() + m.ecalIso() + m.hcalIso() ) / Smpt ;
	sj.SoftLeptIdlooseMu=m.muonID("TMLastStationLoose");
      }
    }
    
     

      for(std::vector<pat::Electron>::const_iterator ele = electrons.begin(); ele!=electrons.end() && sj.isSemiLept!=1; ++ele){
    
	const pat::Electron& e = *ele; 
	float Smpt = e.pt(); 
	float Smeta = e.eta();
	float Smphi = e.phi();
       
	float SmJdR = deltaR(Smeta, Smphi, sj.p4.Eta(), sj.p4.Phi());
	if   ( Smpt> lep_ptCutForBjets_ && SmJdR <0.5)  {
	 sj.isSemiLept=1;
	 sj.SoftLeptpdgId =11;
	 sj.SoftLeptdR= SmJdR;
	 sj.SoftLeptPt=Smpt;
	 TVector3 mvec ( e.p4().Vect().X(), e.p4().Vect().Y(), e.p4().Vect().Z()  ); 
	 sj.SoftLeptptRel=  sj.p4.Perp(  mvec );
	 sj.SoftLeptRelCombIso = (e.trackIso() + e.ecalIso() + e.hcalIso() ) / Smpt ;
	 //	 sj.SoftLeptId95=e.electronID("eidVBTFCom95");
	 //std::cout << "before ele id " << std::endl;      
	 // std::cout << " e.e.sigmaIetaIeta " << e.sigmaIetaIeta() <<  std::endl;
	 //std::cout << " e.isEB() " << e.isEB() << std::endl;
	 if (  
	     ( fabs(Smeta)<2.5 && !( abs(Smeta)>1.4442 && abs(Smeta)<1.566))  && 
              
	     (( abs(Smeta)>1.566  && (e.sigmaIetaIeta()<0.01) && ( e.deltaPhiSuperClusterTrackAtVtx()<0.8  && e.deltaPhiSuperClusterTrackAtVtx()>-0.8) && ( e.deltaEtaSuperClusterTrackAtVtx()<0.007 && e.deltaEtaSuperClusterTrackAtVtx()>-0.007 )  )
	      || ( abs(Smeta)<1.4442  && (e.sigmaIetaIeta()<0.03) && ( e.deltaPhiSuperClusterTrackAtVtx()<0.7 && e.deltaPhiSuperClusterTrackAtVtx()>-0.7 ) && ( e.deltaEtaSuperClusterTrackAtVtx()<0.01 && e.deltaEtaSuperClusterTrackAtVtx()>-0.01 ) ))
	     )
	   sj.SoftLeptId95=1;
       }
     }
  
	
 
    
    
    hbbInfo->simpleJets2.push_back(sj);
    
  }


  /*   const GenJet* jet1Mc = bjet1.genJet();
       const GenJet* jet2Mc = bjet2.genJet();
       if(jet1Mc!=0){
       MCbJet1MomId=jet1Mc->mother()->pdgId();
       MCbJet1GMomId=jet1Mc->mother()->mother()->pdgId();
       }

       if(jet2Mc!=0){
       MCbJet2MomId=jet2Mc->mother()->pdgId();
       MCbJet2GMomId=jet2Mc->mother()->mother()->pdgId();
       }
  */
   
  /* TOM

    // S U B S T R U C T U R E
    // (FastJet3)
    fwlite::Handle<std::vector<pat::Jet> > ca12jetHandle;
    ca12jetHandle.getByLabel(iEvent,"selectedPatJetsCA12PF");
    std::vector<pat::Jet> ca12jets = *ca12jetHandle;
    fwlite::Handle<std::vector<pat::Jet> > camdft12jetHandle;
    camdft12jetHandle.getByLabel(iEvent,"selectedPatJetsCA12MassDropFilteredPF");
    std::vector<pat::Jet> camdft12jets = *camdft12jetHandle;

    fwlite::Handle<std::vector<pat::Jet> > camdft12subjetHandle;
    camdft12subjetHandle.getByLabel(iEvent,"selectedPatJetsCA12MassDropFilteredSubjetsPF");
    std::vector<pat::Jet> camdft12subjets = *camdft12subjetHandle;
    fwlite::Handle<std::vector<pat::Jet> > caft12subjetHandle;
    caft12subjetHandle.getByLabel(iEvent,"selectedPatJetsCA12FilteredSubjetsPF");
    std::vector<pat::Jet> caft12subjets = *caft12subjetHandle;
    fwlite::Handle<std::vector<pat::Jet> > capr12subjetHandle;
    capr12subjetHandle.getByLabel(iEvent,"selectedPatJetsCA12PrunedSubjetsPF");
    std::vector<pat::Jet> capr12subjets = *capr12subjetHandle;

    // C A 1 2   R A W   J E T S
    //std::cout << "Fill CA12 Raw jet!" << std::endl;
    for(std::vector<pat::Jet>::const_iterator jet_iter = ca12jets.begin(); jet_iter!=ca12jets.end(); ++jet_iter){
        const std::vector<reco::PFCandidatePtr>& constituents = jet_iter->getPFConstituents();
        //std::cout << "size: " << constituents.size() << std::endl;
        VHbbEvent::RawJet rj;
        rj.p4 = GENPTOLORP(jet_iter);
        rj.Nconstituents = constituents.size();
        for (unsigned int iJC=0; iJC<constituents.size(); ++iJC) {
            rj.constituents_px.push_back( constituents.at(iJC)->px() );
            rj.constituents_py.push_back( constituents.at(iJC)->py() );
            rj.constituents_pz.push_back( constituents.at(iJC)->pz() );
            rj.constituents_e.push_back( constituents.at(iJC)->energy() );
            rj.constituents_pdgId.push_back( constituents.at(iJC)->pdgId() );
        }
        
        if (jet_iter->genJet()){
            const std::vector <const reco::GenParticle*>& genconstituents = jet_iter->genJet()->getGenConstituents();
            rj.genP4 = GENPTOLORP(jet_iter->genJet());
            rj.Ngenconstituents = genconstituents.size();
            for (unsigned int iJC=0; iJC<genconstituents.size(); ++iJC) {
                rj.genconstituents_px.push_back( genconstituents.at(iJC)->px() );
                rj.genconstituents_py.push_back( genconstituents.at(iJC)->py() );
                rj.genconstituents_pz.push_back( genconstituents.at(iJC)->pz() );
                rj.genconstituents_e.push_back( genconstituents.at(iJC)->energy() );
                rj.genconstituents_pdgId.push_back( genconstituents.at(iJC)->pdgId() );
            }
        }
        hbbInfo->CA12_rawJets.push_back(rj);
    }
    
    // C A 1 2   M A S S D R O P / F I L T E R E D   J E T S
    //std::cout << "Fill CA12 MDFT jet!" << std::endl;
    for(std::vector<pat::Jet>::const_iterator jet_iter = camdft12jets.begin(); jet_iter!=camdft12jets.end(); ++jet_iter){
        if(printJet) {std::cout << "Jet Pt: " << jet_iter->pt() << " E,M: " << jet_iter->p4().E() << " " << jet_iter->p4().M() << "\n";}
        const reco::Jet::Constituents& constituents = jet_iter->getJetConstituents();
        //std::cout << "size: " << constituents.size() << std::endl;
        VHbbEvent::HardJet hj;
        hj.p4 = GENPTOLORP(jet_iter);
        hj.constituents = constituents.size();
        
        for (unsigned int iJC=0; iJC<constituents.size(); ++iJC) {
            const reco::Jet::Constituent& icandJet = constituents.at(iJC);
            if(printJet) {std::cout << "subJet Pt: " << icandJet->pt() << " subJet E,M,eta,phi: " <<  icandJet->p4().E() << "," << icandJet->p4().M() << "," << icandJet->eta() << "," << icandJet->phi() << "\n"; } 
            hj.subFourMomentum.push_back(GENPTOLORP(icandJet));
            hj.etaSub.push_back(icandJet->eta());
            hj.phiSub.push_back(icandJet->phi());
        }
        hbbInfo->CA12mdft_hardJets.push_back(hj);
    }
    
    // C A 1 2   M A S S D R O P / F I L T E R E D   S U B J E T S
    //std::cout << "Fill CA12 MDFT subjet!" << std::endl;    
    for(std::vector<pat::Jet>::const_iterator subjet_iter = camdft12subjets.begin(); subjet_iter!=camdft12subjets.end(); ++subjet_iter) {
        VHbbEvent::SimpleJet sj;
        fillSimpleJet(sj,subjet_iter);
        //setJecUnc(sj,jecUnc);
        hbbInfo->CA12mdft_subJets.push_back(sj);
    }
    
    // C A 1 2   F I L T E R E D   S U B J E T S
    //std::cout << "Fill CA12 FT subjet!" << std::endl;    
    for(std::vector<pat::Jet>::const_iterator subjet_iter = caft12subjets.begin(); subjet_iter!=caft12subjets.end(); ++subjet_iter) {
        VHbbEvent::SimpleJet sj;
        fillSimpleJet(sj,subjet_iter);
        //setJecUnc(sj,jecUnc);
        hbbInfo->CA12ft_subJets.push_back(sj);
    }
    
    // C A 1 2   P R U N E D   S U B J E T S
    //std::cout << "Fill CA12 PR subjet!" << std::endl;    
    for(std::vector<pat::Jet>::const_iterator subjet_iter = capr12subjets.begin(); subjet_iter!=capr12subjets.end(); ++subjet_iter) {
        VHbbEvent::SimpleJet sj;
        fillSimpleJet(sj,subjet_iter);
        //setJecUnc(sj,jecUnc);
        hbbInfo->CA12pr_subJets.push_back(sj);
    }

    // S U B S T R U C T U R E
    // (FastJet2)
    for(std::vector<pat::Jet>::const_iterator jet_iter = fatjets.begin(); jet_iter!=fatjets.end(); ++jet_iter){
        if(printJet) {std::cout << "Jet Pt: " << jet_iter->pt() << " E,M: " << jet_iter->p4().E() << " " << jet_iter->p4().M() << "\n";}
        const reco::Jet::Constituents& constituents = jet_iter->getJetConstituents();
        //std::cout << "size: " << constituents.size() << std::endl;
        VHbbEvent::HardJet hj;
        hj.p4 = GENPTOLORP(jet_iter);
        hj.constituents=constituents.size();
        
        for (unsigned int iJC=0; iJC<constituents.size(); ++iJC ){
            const Jet::Constituent& icandJet = constituents.at(iJC);
            if(printJet) {std::cout << "subJet Pt: " << icandJet->pt() << " subJet E,M,eta,phi: " <<  icandJet->p4().E() << "," << icandJet->p4().M() << "," << icandJet->eta() << "," << icandJet->phi() << "\n"; } 
            hj.subFourMomentum.push_back(GENPTOLORP(icandJet));
            hj.etaSub.push_back(icandJet->eta());
            hj.phiSub.push_back(icandJet->phi());
        }
        hbbInfo->hardJets.push_back(hj);
    }

    //double matEta[1000*30],matPhi[1000*30];
    //for(int i=0;i<1000;i++){for(int j=0;j<30;j++){matEta[i*j]=-99.;matPhi[i*j]=-99.;}}
    //HardJetSubEta2.SetMatrixArray(matEta);
    //HardJetSubPhi2.SetMatrixArray(matPhi);
    //TMatrixDRow a1(HardJetSubEta2,0);
    //for(int i=0;i<30;i++){
    //  std::cout << "test: " << a1[i] << "\n";  
    //}

    */

    /* TOM

    for(std::vector<pat::Jet>::const_iterator subjet_iter = subjets.begin(); subjet_iter!=subjets.end(); ++subjet_iter) {
        if(printJet) {std::cout << "SubJetTagged Pt: " << subjet_iter->pt() << " E,M,eta,phi,Btag: " << subjet_iter->p4().E() 
            << "," << subjet_iter->p4().M() << "," << subjet_iter->eta() << "," << subjet_iter->phi() 
            << "," << subjet_iter->bDiscriminator("combinedSecondaryVertexBJetTags") << "\n";}
        VHbbEvent::SimpleJet sj;
        fillSimpleJet(sj,subjet_iter);
	// TOM         setJecUnc(sj,jecUnc);
        hbbInfo->subJets.push_back(sj);
    }

   TOM  */

    /* TOM
    for(std::vector<pat::Jet>::const_iterator subjet_iter = filterjets.begin(); subjet_iter!=filterjets.end(); ++subjet_iter) {
        if(printJet) {std::cout << "FilterjetTagged Pt: " << subjet_iter->pt() << " E,M,eta,phi,Btag: " << subjet_iter->p4().E() 
            << "," << subjet_iter->p4().M() << "," << subjet_iter->eta() << "," << subjet_iter->phi() 
            << "," << subjet_iter->bDiscriminator("combinedSecondaryVertexBJetTags") << "\n";}
        VHbbEvent::SimpleJet sj;
        fillSimpleJet(sj,subjet_iter);
	// TOM         setJecUnc(sj,jecUnc);
        

        if(runOnMC_) {
            // BTV scale factors
            //fillScaleFactors(sj, btagSFs);

            // physical parton for mother info ONLY
            if( (subjet_iter->genParton()) ) {
                sj.bestMCid = subjet_iter->genParton()->pdgId();
                if( (subjet_iter->genParton()->mother()) )
                    sj.bestMCmomid=subjet_iter->genParton()->mother()->pdgId();
            }
            // genJet
            const reco::GenJet * gJ = subjet_iter->genJet();
            TLorentzVector gJp4;
            if(gJ) {
                gJp4.SetPtEtaPhiE(gJ->pt(),gJ->eta(),gJ->phi(),gJ->energy());
                sj.bestMCp4 = gJp4;
            }
        }

	//ccla    
    

	for(std::vector<reco::Candidate>::const_iterator mu = muonsNoCuts.begin(); mu!=muonsNoCuts.end() && sj.isSemiLept!=1; ++mu){
	  //      std::cout<< "found a muon with pt " << mu->pt()   << std::endl;
	  const pat::Muon& m = static_cast <const pat::Muon&> (*mu); 
	  float Smpt = m.pt(); 
	  float Smeta = m.eta();
	  float Smphi = m.phi();
      
	  float SmJdR = deltaR(Smeta, Smphi, sj.p4.Eta(), sj.p4.Phi());
      
	  if   ( Smpt> lep_ptCutForBjets_ && SmJdR <0.5)  {
	    sj.isSemiLept=1;
	    //isSemiLept(-99), isSemiLeptMCtruth(-99), SoftLeptPt(-99), SoftLeptdR(-99), SoftLeptptRel(-99), SoftLeptpdgId(-99), SoftLeptIdlooseMu(-99), SoftLeptId95(-99), SoftLeptRelCombIso(-99),  
	    sj.SoftLeptpdgId =13;
	    sj.SoftLeptdR= SmJdR;
	    sj.SoftLeptPt=Smpt;
	    TVector3 mvec ( m.p4().Vect().X(), m.p4().Vect().Y(), m.p4().Vect().Z()  ); 
	    sj.SoftLeptptRel=  sj.p4.Perp(  mvec );
	    sj.SoftLeptRelCombIso = (m.trackIso() + m.ecalIso() + m.hcalIso() ) / Smpt ;
	    sj.SoftLeptIdlooseMu=m.muonID("TMLastStationLoose");
	  }
	}
    
    
	//edm::Handle<std::vector<reco::Candidate> > eleNoCutsHandle;
	//iEvent.getByLabel(elenoCutsLabel_,eleNoCutsHandle);
	//std::vector<reco::Candidate> elesNoCuts = *eleNoCutsHandle; 
 
    

	for(std::vector<reco::Candidate>::const_iterator ele = elesNoCuts.begin(); ele!=elesNoCuts.end() && sj.isSemiLept!=1; ++ele){
    
	  const pat::Electron& e = static_cast <const pat::Electron&> (*ele); 
	  float Smpt = e.pt(); 
	  float Smeta = e.eta();
	  float Smphi = e.phi();
       
	  float SmJdR = deltaR(Smeta, Smphi, sj.p4.Eta(), sj.p4.Phi());
	  if   ( Smpt> lep_ptCutForBjets_ && SmJdR <0.5)  {
	    sj.isSemiLept=1;
	    sj.SoftLeptpdgId =11;
	    sj.SoftLeptdR= SmJdR;
	    sj.SoftLeptPt=Smpt;
	    TVector3 mvec ( e.p4().Vect().X(), e.p4().Vect().Y(), e.p4().Vect().Z()  ); 
	    sj.SoftLeptptRel=  sj.p4.Perp(  mvec );
	    sj.SoftLeptRelCombIso = (e.trackIso() + e.ecalIso() + e.hcalIso() ) / Smpt ;
	    //	 sj.SoftLeptId95=e.electronID("eidVBTFCom95");
	    //std::cout << "before ele id " << std::endl;      
	    // std::cout << " e.e.sigmaIetaIeta " << e.sigmaIetaIeta() <<  std::endl;
	    //std::cout << " e.isEB() " << e.isEB() << std::endl;
	    if (  
		( fabs(Smeta)<2.5 && !( abs(Smeta)>1.4442 && abs(Smeta)<1.566))  && 
              
		(( abs(Smeta)>1.566  && (e.sigmaIetaIeta()<0.01) && ( e.deltaPhiSuperClusterTrackAtVtx()<0.8  && e.deltaPhiSuperClusterTrackAtVtx()>-0.8) && ( e.deltaEtaSuperClusterTrackAtVtx()<0.007 && e.deltaEtaSuperClusterTrackAtVtx()>-0.007 )  )
		 || ( abs(Smeta)<1.4442  && (e.sigmaIetaIeta()<0.03) && ( e.deltaPhiSuperClusterTrackAtVtx()<0.7 && e.deltaPhiSuperClusterTrackAtVtx()>-0.7 ) && ( e.deltaEtaSuperClusterTrackAtVtx()<0.01 && e.deltaEtaSuperClusterTrackAtVtx()>-0.01 ) ))
		  )
	      sj.SoftLeptId95=1;
	  }
	}


        hbbInfo->filterJets.push_back(sj);
    }

    TOM */

  //
  // add charged met
  //
  
  /* TOM

  fwlite::Handle<std::vector<reco::MET> > metChargedHandle;
  metChargedHandle.getByLabel(iEvent,"pfMETNoPUCharge");
  std::vector<reco::MET> metsCh = *metChargedHandle;
  if(metsCh.size()){
    hbbInfo->metCh.sumEt=(metsCh[0]).sumEt();
    hbbInfo->metCh.metSig=metSignificance(& (metsCh[0]));
    hbbInfo->metCh.eLong=(metsCh[0]).e_longitudinal();
    hbbInfo->metCh.p4=GENPTOLOR((metsCh[0]));
    if (verbose_)     std::cout <<" METCharged "<<     hbbInfo->metCh.metSig <<" " <<     hbbInfo->metCh.sumEt<<std::endl;
  }

  */


//AR: MET and its uncerts
//  
  //default met is type1 corrected
  
  edm::Handle<std::vector<pat::MET> > metHandle;
  iEvent.getByLabel(metLabel_,metHandle);
  const pat::MET & met = (*metHandle)[0];
  hbbInfo->pfmetType1corr.sumEt=met.sumEt();
  hbbInfo->pfmetType1corr.metSig=metSignificance(&met);
  hbbInfo->pfmetType1corr.eLong=met.e_longitudinal();
  hbbInfo->pfmetType1corr.p4=GENPTOLOR(met);
  if (verbose_)     std::cout <<" type 1 corrected pfMET "<<     hbbInfo->pfmetType1corr.metSig <<" " <<     hbbInfo->pfmetType1corr.sumEt<<std::endl;

/*TOM AR not available
  // type 1 + 2 corr met
  edm::Handle<std::vector<reco::MET> > pfmetType1p2corrHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMet",pfmetType1p2corrHandle);
  std::vector<reco::MET> pfmetsType1p2corr = *pfmetType1p2corrHandle;
  if(pfmetsType1p2corr.size()){
    hbbInfo->pfmetType1p2corr.sumEt=(pfmetsType1p2corr[0]).sumEt();
    hbbInfo->pfmetType1p2corr.metSig=metSignificance(& (pfmetsType1p2corr[0]));
    hbbInfo->pfmetType1p2corr.eLong=(pfmetsType1p2corr[0]).e_longitudinal();
    hbbInfo->pfmetType1p2corr.p4=GENPTOLOR((pfmetsType1p2corr[0]));
    if (verbose_)     std::cout <<" type 1 +2 corrected pfMET "<<     hbbInfo->pfmetType1p2corr.metSig <<" " <<     hbbInfo->pfmetType1p2corr.sumEt<<std::endl;
  }

  // type 1 corr met NoPU
  /--*  edm::Handle<std::vector<reco::MET> > pfmetNoPUType1corrHandle;
  iEvent.getByLabel("patType1CorrectedPFMetNoPU",pfmetNoPUType1corrHandle);
  std::vector<reco::MET> pfmetsNoPUType1corr = *pfmetNoPUType1corrHandle;
  if(pfmetsNoPUType1corr.size()){
    hbbInfo->pfmetNoPUType1corr.sumEt=(pfmetsNoPUType1corr[0]).sumEt();
    hbbInfo->pfmetNoPUType1corr.metSig=metSignificance(& (pfmetsNoPUType1corr[0]));
    hbbInfo->pfmetNoPUType1corr.eLong=(pfmetsNoPUType1corr[0]).e_longitudinal();
    hbbInfo->pfmetNoPUType1corr.p4=GENPTOLOR((pfmetsNoPUType1corr[0]));
    if (verbose_)     std::cout <<" type 1 corrected pfMET NoPU"<<     hbbInfo->pfmetNoPUType1corr.metSig <<" " <<     hbbInfo->pfmetNoPUType1corr.sumEt<<std::endl;
  }


  // type 1 + 2 corr met
  edm::Handle<std::vector<reco::MET> > pfmetNoPUType1p2corrHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMetNoPU",pfmetNoPUType1p2corrHandle);
  std::vector<reco::MET> pfmetsNoPUType1p2corr = *pfmetNoPUType1p2corrHandle;
  if(pfmetsNoPUType1p2corr.size()){
    hbbInfo->pfmetNoPUType1p2corr.sumEt=(pfmetsNoPUType1p2corr[0]).sumEt();
    hbbInfo->pfmetNoPUType1p2corr.metSig=metSignificance(& (pfmetsNoPUType1p2corr[0]));
    hbbInfo->pfmetNoPUType1p2corr.eLong=(pfmetsNoPUType1p2corr[0]).e_longitudinal();
    hbbInfo->pfmetNoPUType1p2corr.p4=GENPTOLOR((pfmetsNoPUType1p2corr[0]));
    if (verbose_)     std::cout <<" type 1 +2 corrected pfMET "<<     hbbInfo->pfmetNoPUType1p2corr.metSig <<" " <<     hbbInfo->pfmetNoPUType1p2corr.sumEt<<std::endl;
  }

*--/

/--*
  // MET uncertainty vector
    vector<pat::MET>                 "patType1CorrectedPFMet"    ""                "VH"      
    vector<pat::MET>                 "patType1CorrectedPFMetElectronEnDown"   ""                "VH"      
    vector<pat::MET>                 "patType1CorrectedPFMetElectronEnUp"   ""                "VH"      
    vector<pat::MET>                 "patType1CorrectedPFMetJetEnDown"   ""                "VH"      
    vector<pat::MET>                 "patType1CorrectedPFMetJetEnUp"   ""                "VH"      
    vector<pat::MET>                 "patType1CorrectedPFMetJetResDown"   ""                "VH"      
    vector<pat::MET>                 "patType1CorrectedPFMetJetResUp"   ""                "VH"      
    vector<pat::MET>                 "patType1CorrectedPFMetMuonEnDown"   ""                "VH"      
    vector<pat::MET>                 "patType1CorrectedPFMetMuonEnUp"   ""                "VH"      
    vector<pat::MET>                 "patType1CorrectedPFMetNoPU"   ""                "VH"      
    vector<pat::MET>                 "patType1CorrectedPFMetTauEnDown"   ""                "VH"      
    vector<pat::MET>                 "patType1CorrectedPFMetTauEnUp"   ""                "VH"      
    vector<pat::MET>                 "patType1CorrectedPFMetUnclusteredEnDown"   ""                "VH"      
    vector<pat::MET>                 "patType1CorrectedPFMetUnclusteredEnUp"   ""                "VH"      
    vector<pat::MET>                 "patType1p2CorrectedPFMet"   ""                "VH"      
    vector<pat::MET>                 "patType1p2CorrectedPFMetElectronEnDown"   ""                "VH"      
    vector<pat::MET>                 "patType1p2CorrectedPFMetElectronEnUp"   ""                "VH"      
    vector<pat::MET>                 "patType1p2CorrectedPFMetJetEnDown"   ""                "VH"      
    vector<pat::MET>                 "patType1p2CorrectedPFMetJetEnUp"   ""                "VH"      
    vector<pat::MET>                 "patType1p2CorrectedPFMetJetResDown"   ""                "VH"      
    vector<pat::MET>                 "patType1p2CorrectedPFMetJetResUp"   ""                "VH"      
    vector<pat::MET>                 "patType1p2CorrectedPFMetMuonEnDown"   ""                "VH"      
    vector<pat::MET>                 "patType1p2CorrectedPFMetMuonEnUp"   ""                "VH"      
    vector<pat::MET>                 "patType1p2CorrectedPFMetNoPU"   ""                "VH"      
    vector<pat::MET>                 "patType1p2CorrectedPFMetTauEnDown"   ""                "VH"      
    vector<pat::MET>                 "patType1p2CorrectedPFMetTauEnUp"   ""                "VH"      
    vector<pat::MET>                 "patType1p2CorrectedPFMetUnclusteredEnDown"   ""                "VH"      
    vector<pat::MET>                 "patType1p2CorrectedPFMetUnclusteredEnUp"   ""                "VH"      
      *--/


*/
//      enum METUncertainty {
  //      JetEnUp=0, JetEnDown=1, JetResUp=2, JetResDown=3,
    //    MuonEnUp=4, MuonEnDown=5, ElectronEnUp=6, ElectronEnDown=7, TauEnUp=8,TauEnDown=9,
      //  UnclusteredEnUp=10,UnclusteredEnDown=11, METUncertaintySize=12
     // };

    for( int i=0; i < 12;i++)   {
	    VHbbEvent::METInfo metunc;
	    metunc.sumEt =met.shiftedSumEt(pat::MET::METUncertainty(i)); 
	    metunc.metSig=0; //AR not available: metSignificance(& (patType1CorrectedPFMetsElectronEnDown[0]));
	    metunc.eLong=0; //AR not available: (patType1CorrectedPFMetsElectronEnDown[0]).e_longitudinal();
	    metunc.p4=GENPTOLOR(met.shiftedP4(pat::MET::METUncertainty(i)));
	    hbbInfo->metUncInfo.push_back(metunc); 
    }
    for( int i=0; i < 12;i++)   {
            VHbbEvent::METInfo metunc;
            metunc.sumEt =met.shiftedSumEt(pat::MET::METUncertainty(i),pat::MET::Type1p2);
            metunc.metSig=0; //AR not available: metSignificance(& (patType1CorrectedPFMetsElectronEnDown[0]));
            metunc.eLong=0; //AR not available: (patType1CorrectedPFMetsElectronEnDown[0]).e_longitudinal();
            metunc.p4=GENPTOLOR(met.shiftedP4(pat::MET::METUncertainty(i),pat::MET::Type1p2));
            hbbInfo->metUncInfo.push_back(metunc);
    }

  


  /* TOM missign the other METs... compute on the fly?

  edm::Handle<std::vector<pat::MET> > metTCHandle;
  iEvent.getByLabel("patMETsTC",metTCHandle);
  std::vector<pat::MET> metsTC = *metTCHandle;
  if(metsTC.size()){
    hbbInfo->tcmet.sumEt=(metsTC[0]).sumEt();
    hbbInfo->tcmet.metSig=metSignificance(&(metsTC[0]));
    hbbInfo->tcmet.eLong=(metsTC[0]).e_longitudinal();
    hbbInfo->tcmet.p4=GENPTOLOR((metsTC[0]));
    if (verbose_)     std::cout <<" METTC "<<     hbbInfo->tcmet.metSig <<" " <<     hbbInfo->tcmet.sumEt<<std::endl;
  }
  
  edm::Handle<std::vector<reco::MET> > pfMETNoPUHandle;
  iEvent.getByLabel("pfMETNoPU",pfMETNoPUHandle);
  std::vector<reco::MET> metspfMETNoPU = *pfMETNoPUHandle;
  if(metspfMETNoPU.size()){
    hbbInfo->metNoPU.sumEt=(metspfMETNoPU[0]).sumEt();
    hbbInfo->metNoPU.metSig=metSignificance(&(metspfMETNoPU[0]));                
    hbbInfo->metNoPU.eLong=(metspfMETNoPU[0]).e_longitudinal();
    hbbInfo->metNoPU.p4=GENPTOLOR((metspfMETNoPU[0]));
    if (verbose_)     std::cout <<" pfMETNoPU "<<     hbbInfo->metNoPU.metSig <<" " <<     hbbInfo->metNoPU.sumEt<<std::endl;
  }
 
  edm::Handle<std::vector<reco::MET> > mHTHandle;
  iEvent.getByLabel("patMETsHT",mHTHandle);
  std::vector<reco::MET> metsHT = *mHTHandle;
  if(metsHT.size()){
    hbbInfo->mht.sumEt=(metsHT[0]).sumEt();
    hbbInfo->mht.metSig=metSignificance(&(metsHT[0]));
    hbbInfo->mht.eLong=(metsHT[0]).e_longitudinal();
    hbbInfo->mht.p4=GENPTOLOR((metsHT[0]));
    if (verbose_)     std::cout <<" METHT "<<     hbbInfo->mht.metSig <<" " <<     hbbInfo->mht.sumEt<<std::endl;
  }
  
  edm::Handle<std::vector<reco::MET> > metHandle;
  iEvent.getByLabel(metLabel_,metHandle);
  std::vector<reco::MET> mets = *metHandle;
  
  if(mets.size()){
    hbbInfo->calomet.sumEt=(mets[0]).sumEt();
    hbbInfo->calomet.metSig=metSignificance(&(mets[0]));
    hbbInfo->calomet.eLong=(mets[0]).e_longitudinal();
    hbbInfo->calomet.p4=GENPTOLOR((mets[0]));
    if (verbose_)     std::cout <<" METCALO "<<     hbbInfo->calomet.metSig <<" " <<     hbbInfo->calomet.sumEt<<std::endl;
  }
  
  edm::Handle<std::vector<pat::MET> > metPFHandle;
  iEvent.getByLabel("patMETsPFlow",metPFHandle);
  std::vector<pat::MET> metsPF = *metPFHandle;
 */ 

  hbbInfo->pfmet.sumEt=met.corSumEt();//uncorrect
  hbbInfo->pfmet.metSig=0; //AR n.a.
  hbbInfo->pfmet.eLong=0; //AR n.a
  hbbInfo->pfmet.p4=TLorentzVector(met.px()*met.uncorrectedPt()/met.pt(), met.py()*met.uncorrectedPt()/met.pt(), met.pz(), met.energy());
  
  

//  if(verbose_)
//    std::cout << " INPUT MUONS "<<muons.size()<<std::endl;



  for(std::vector<pat::Muon>::const_iterator mu = muons.begin(); mu!=muons.end(); ++mu){
    VHbbEvent::MuonInfo mf;
    mf.p4 =GENPTOLORP( mu);
    mf.charge=mu->charge();
    mf.tIso=mu->trackIso();
    mf.eIso=mu->ecalIso();
    mf.hIso=mu->hcalIso();
    mf.pfChaIso=mu->chargedHadronIso();
    mf.pfChaPUIso=mu->puChargedHadronIso(); //userIso(5);
    mf.pfPhoIso=mu->photonIso();
    mf.pfNeuIso=mu->neutralHadronIso(); 
    Geom::Phi<double> deltaphi(mu->phi()-atan2(mf.p4.Px(), mf.p4.Py()));
    double acop = deltaphi.value();
    mf.acop=acop;

    mf.emEnergy = mu->calEnergy().em;
    mf.hadEnergy = mu->calEnergy().had;

    mf.nMatches = mu->numberOfMatches();

    mf.ipDb=mu->dB();
    mf.ipErrDb=mu->edB();
    mf.cat=0;

    if(mu->isGlobalMuon()) mf.cat|=1;
    if(mu->isTrackerMuon()) mf.cat|=2;
    if(mu->isStandAloneMuon()) mf.cat|=4;
    TrackRef trkMu1Ref = mu->get<TrackRef>();
    if(trkMu1Ref.isNonnull()){
      const Track* MuTrk1 = mu->get<TrackRef>().get();
      mf.zPVPt=MuTrk1->dz(RecVtxFirst.position());
      mf.zPVProb=MuTrk1->dz(RecVtx.position());
      mf.nHits=MuTrk1->numberOfValidHits();
      mf.chi2=MuTrk1->normalizedChi2();
      TrackRef iTrack1 = mu->innerTrack();
      const reco::HitPattern& p1 = iTrack1->hitPattern();
      mf.nPixelHits=p1.pixelLayersWithMeasurement();

      mf.nValidTracker = p1.numberOfValidTrackerHits(); 
      mf.nValidPixel = p1.numberOfValidPixelHits(); 
      mf.nValidLayers = p1.trackerLayersWithMeasurement();
      mf.isPF = mu->isPFMuon();



   } 
    if(mu->isGlobalMuon()){
      TrackRef gTrack = mu->globalTrack();
      const reco::HitPattern& q = gTrack->hitPattern();
      mf.globChi2=gTrack.get()->normalizedChi2();
      mf.globNHits=q.numberOfValidMuonHits();
      mf.validMuStations = q.muonStationsWithValidHits();
    }else{
      mf.globChi2=-99;
      mf.globNHits=-99;
    }

    //Muon trigger matching
    for (int itrig = 0; itrig != ntrigs; ++itrig){
      std::string trigName=triggerNames_.triggerName(itrig);
      if( (mu->triggerObjectMatchesByPath(trigName,false,false).size() != 0) ){
	mf.hltMatchedBits.push_back(itrig);
	if(verbose_){
	  std::clog << "Trigger Matching box" << std::endl;
	  std::clog << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	  std::clog << "Matching parameters are defined in the cfg" << std::endl;
	  std::clog << "Trigger bit = " << itrig << std::endl;
	  std::clog << "Trigger name = " << trigName << std::endl;
	  std::clog << "Trigger object matched collection size = " << mu->triggerObjectMatchesByPath(trigName,false,false).size() << std::endl;
	  std::clog << "Pat Muon pt = " << mf.p4.Pt() << " HLT object matched = " << mu->triggerObjectMatch(0)->p4().Pt() << std::endl;
	  std::clog << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	}
      }
    }
    //

    // add stamuon 

    //    if (mu->isStandAloneMuon()) {
    //      reco::TrackRef sta = mu->standAloneMuon();
    //      
    //    }


    //     int muInfo[12];
    //     fillMuBlock(mu,  muInfo);
    if(runOnMC_){
      const GenParticle* muMc = mu->genLepton();
      if(muMc!=0){
	mf.mcId=muMc->pdgId();
	mf.mcFourMomentum=GENPTOLORP(muMc);
	if(muMc->mother()!=0) mf.mcMomId=muMc->mother()->pdgId();
	if(muMc->mother()!=0 && muMc->mother()->mother()!=0) mf.mcgMomId=muMc->mother()->mother()->pdgId();
      } }
    hbbInfo->muInfo.push_back(mf);
  }

 if(verbose_)
    std::cout << " INPUT electrons "<<electrons.size()<<std::endl;

 // TOM WHAT TO DO HERE?

 /*  InputTag  reducedEBRecHitCollection(string("reducedEcalRecHitsEB"));
  InputTag  reducedEERecHitCollection(string("reducedEcalRecHitsEE"));
  EcalClusterLazyTools lazyTools(iEvent, iSetup, reducedEBRecHitCollection, reducedEERecHitCollection);
  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  const TransientTrackBuilder & transientTrackBuilder= *(builder.product());
 */


 for(std::vector<pat::Electron>::const_iterator elec = electrons.begin(); elec!=electrons.end(); ++elec){
    VHbbEvent::ElectronInfo ef;
    ef.p4=GENPTOLORP(elec);
    ef.scEta =elec->superCluster()->eta();
    ef.scPhi =elec->superCluster()->phi();
    //    if(ElecEta[eleccont]!=0) ElecEt[eleccont]=elec->superCluster()->energy()/cosh(elec->superCluster()->eta());
    ef.charge=elec->charge();
    ef.tIso=elec->trackIso();
    ef.eIso=elec->ecalIso();
    ef.hIso=elec->hcalIso();
    ef.pfChaIso=elec->chargedHadronIso();
    ef.pfChaPUIso=elec->puChargedHadronIso();//userIso(5);
    ef.pfPhoIso=elec->photonIso();
    ef.pfPhoIsoDoubleCounted=0;

/* Check if there are photons sharing the super cluster*/
/* TOM
    for(size_t k=0;k<photonsForIso.size();k++) {
       if(deltaR(elec->eta(),elec->phi(),photonsForIso[k].eta(),photonsForIso[k].phi()) < 0.05 && abs(photonsForIso[k].pt()-elec->pt())/(photonsForIso[k].pt()+elec->pt()) < 0.05 ) {
          std::cout << "Double counting of supercluster!" << std::endl;
	  ef.pfPhoIsoDoubleCounted+=photonsForIso[k].pt(); 
     }
    }
 TOM */
    ef.pfNeuIso=elec->neutralHadronIso();

    //
    // ip info
    //

    ef.ipDb=elec->dB();
    ef.ipErrDb=elec->edB();
    


    Geom::Phi<double> deltaphi(elec->superCluster()->phi()-atan2(hbbInfo->pfmet.p4.Py(),hbbInfo->pfmet.p4.Px()));
    ef.acop = deltaphi.value();
    //
    ef.sihih = elec->sigmaIetaIeta();
    ef.Dphi = elec->deltaPhiSuperClusterTrackAtVtx();
    ef.Deta = elec->deltaEtaSuperClusterTrackAtVtx();
    ef.HoE = elec->hadronicOverEm();
    ef.convDist = elec->convDist();
    ef.convDcot = elec->convDcot();
    if(elec->gsfTrack().isNonnull()) 
    {
     ef.innerHits = elec->gsfTrack()->trackerExpectedHitsInner().numberOfHits();   
    }
    ef.isEB = elec->isEB();
    ef.isEE = elec->isEE();
/* 2012 ELEID*/

  const pat::Electron & ele = *elec;
  bool validKF= false; 
  reco::TrackRef myTrackRef = ele.closestCtfTrackRef();
  validKF = (myTrackRef.isAvailable());
  validKF = (myTrackRef.isNonnull());  

  // Pure tracking variables
  ef.fMVAVar_fbrem           =  ele.fbrem();
  ef.fMVAVar_kfchi2          =  (validKF) ? myTrackRef->normalizedChi2() : 0 ;
  ef.fMVAVar_kfhits          =  (validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. ; 
  //  fMVAVar_kfhitsall          =  (validKF) ? myTrackRef->numberOfValidHits() : -1. ;   //  save also this in your ntuple as possible alternative
  ef.fMVAVar_gsfchi2         =  ele.gsfTrack()->normalizedChi2();  

  
  // Geometrical matchings
  ef.fMVAVar_deta            =  ele.deltaEtaSuperClusterTrackAtVtx();
  ef.fMVAVar_dphi            =  ele.deltaPhiSuperClusterTrackAtVtx();
  ef.fMVAVar_detacalo        =  ele.deltaEtaSeedClusterTrackAtCalo();
  // fMVAVar_dphicalo        =  ele.deltaPhiSeedClusterTrackAtCalo();   //  save also this in your ntuple 


  // Pure ECAL -> shower shapes
  ef.fMVAVar_see             =  ele.userFloat("sigmaIetaIeta_NoZS"); //AR: ele.sigmaIetaIeta();    //EleSigmaIEtaIEta
  ef.fMVAVar_spp 	     =  ele.userFloat("sigmaIetaIphi_NoZS");

  // fMVAVar_sigmaIEtaIPhi = vCov[1];  //  save also this in your ntuple 

  ef.fMVAVar_etawidth        =  ele.superCluster()->etaWidth();
  ef.fMVAVar_phiwidth        =  ele.superCluster()->phiWidth();
  ef.fMVAVar_e1x5e5x5        =  ele.userFloat("e1x5_over_e5x5_NoZS"); //AR:(ele.e5x5()) !=0. ? 1.-(ele.e1x5()/ele.e5x5()) : -1. ;
  ef.fMVAVar_R9		     = ele.userFloat("r9_NoZS");
  //fMVAVar_nbrems          =  fabs(ele.numberOfBrems());    //  save also this in your ntuple 

  // Energy matching
  ef.fMVAVar_HoE             =  ele.hadronicOverEm();
  ef.fMVAVar_EoP             =  ele.eSuperClusterOverP();
  // fMVAVar_IoEmIoP         =  (1.0/(ele.superCluster()->energy())) - (1.0 / ele.p());  // in the future to be changed with ele.gsfTrack()->p()
  ef.fMVAVar_IoEmIoP         =  (1.0/ele.ecalEnergy()) - (1.0 / ele.p());  // in the future to be changed with ele.gsfTrack()->p()   // 24/04/2012 changed to correctly access the   corrected supercluster energy from CMSSW_52X

  ef.fMVAVar_eleEoPout       =  ele.eEleClusterOverPout();
  ef.fMVAVar_PreShowerOverRaw=  ele.superCluster()->preshowerEnergy() / ele.superCluster()->rawEnergy();
  // fMVAVar_EoPout          =  ele.eSeedClusterOverPout();     //  save also this in your ntuple 


  // Spectators
  ef.fMVAVar_eta             =  ele.superCluster()->eta();         
  ef.fMVAVar_pt              =  ele.pt();                          

  //additional for cut based
  ef.dxy = elec->gsfTrack()->dxy(vertex.position());
  ef.dz  = elec->gsfTrack()->dz(vertex.position());


    //d0
    if (ele.gsfTrack().isNonnull()) {
      ef.fMVAVar_d0 = (-1.0)*ele.gsfTrack()->dxy(vertex.position()); 
    } else if (ele.closestCtfTrackRef().isNonnull()) {
      ef.fMVAVar_d0 = (-1.0)*ele.closestCtfTrackRef()->dxy(vertex.position()); 
    } else {
      ef.fMVAVar_d0 = -9999.0;
    
    //default values for IP3D
    ef.fMVAVar_ip3d = ele.ip3d();
 
    }
  
  

/* end of 2012 ELEID*/

    //
    // fill eleids
    //    
/*    ef.id95 = elec->electronID("simpleEleId95cIso");
    ef.id85 = elec->electronID("simpleEleId85cIso");
    ef.id70 = elec->electronID("simpleEleId70cIso");
    ef.id95r = elec->electronID("simpleEleId95relIso");
    ef.id70r = elec->electronID("simpleEleId70relIso");
    ef.id85r = elec->electronID("simpleEleId85relIso");
*/
    ef.loose =elec->electronID("eidLoose");
    ef.robustLoose=elec->electronID("eidRobustLoose");
    ef.tight =elec->electronID("eidTight");
    ef.robustTight=elec->electronID("eidRobustTight");
    ef.robustHighEnergy=elec->electronID("eidRobustHighEnergy");

    //Electron trigger matching
    for (int itrig = 0; itrig != ntrigs; ++itrig){
      std::string trigName=triggerNames_.triggerName(itrig);
      if( (elec->triggerObjectMatchesByPath(trigName).size() != 0) ){
	ef.hltMatchedBits.push_back(itrig);
	if(verbose_){
	  std::clog << "Trigger Matching box" << std::endl;
	  std::clog << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	  std::clog << "Matching parameters are defined in the cfg" << std::endl;
	  std::clog << "Trigger bit = " << itrig << std::endl;
	  std::clog << "Trigger name = " << trigName << std::endl;
	  std::clog << "Trigger object matched collection size = " << elec->triggerObjectMatchesByPath(trigName).size() << std::endl;
	  std::clog << "Pat Electron pt = " << ef.p4.Pt() << " HLT object matched = " << elec->triggerObjectMatch(0)->p4().Pt() << std::endl;
	  std::clog << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	}
      }
    }

    if(runOnMC_){
      const GenParticle* elecMc = elec->genLepton();
      if(elecMc!=0){
	ef.mcId=elecMc->pdgId();
	ef.mcFourMomentum=GENPTOLORP(elecMc);
	if(elecMc->mother()!=0) ef.mcMomId=elecMc->mother()->pdgId();
	if(elecMc->mother()!=0 && elecMc->mother()->mother()!=0) ef.mcgMomId=elecMc->mother()->mother()->pdgId();
      }}
    hbbInfo->eleInfo.push_back(ef);
  }

  if(verbose_)
    std::cout << " INPUT taus "<<taus.size()<<std::endl;
  for(std::vector<pat::Tau>::const_iterator tau = taus.begin(); tau!=taus.end(); ++tau){
    VHbbEvent::TauInfo tf;
    tf.p4=GENPTOLORP(tau);
    tf.charge=tau->charge();
    tf.tIso=tau->trackIso();
    tf.eIso=tau->ecalIso();
    tf.hIso=tau->hcalIso();
    Geom::Phi<double> deltaphi(tau->phi()-atan2(hbbInfo->pfmet.p4.Py(),hbbInfo->pfmet.p4.Px()));
    double acop = deltaphi.value();
    tf.acop=acop;
    if (tau->isTauIDAvailable("againstElectronLoose")) tf.againstElectronLoose=tau->tauID("againstElectronLoose");
    if (tau->isTauIDAvailable("againstElectronMedium")) tf.againstElectronMedium=tau->tauID("againstElectronMedium");
    if (tau->isTauIDAvailable("againstElectronTight")) tf.againstElectronTight=tau->tauID("againstElectronTight");
    if (tau->isTauIDAvailable("againstMuonLoose")) tf.againstMuonLoose=tau->tauID("againstMuonLoose");
    if (tau->isTauIDAvailable("againstMuonTight")) tf.againstMuonTight=tau->tauID("againstMuonTight");
    if (tau->isTauIDAvailable("byLooseIsolation")) tf.byLooseIsolation=tau->tauID("byLooseIsolation");
    if (tau->isTauIDAvailable("byMediumIsolation")) tf.byMediumIsolation=tau->tauID("byMediumIsolation");
    if (tau->isTauIDAvailable("byTightIsolation")) tf.byTightIsolation=tau->tauID("byTightIsolation");
    if (tau->isTauIDAvailable("byVLooseIsolation")) tf.byVLooseIsolation=tau->tauID("byVLooseIsolation");
    if (tau->isTauIDAvailable("decayModeFinding")) tf.decayModeFinding=tau->tauID("decayModeFinding");
    if (tau->isTauIDAvailable("byIsolation")) tf.byIsolation=tau->tauID("byIsolation");
    if (tau->isTauIDAvailable("trackIsolation")) tf.trackIsolation=tau->tauID("trackIsolation");
    if (tau->isTauIDAvailable("byTaNCfrOnePercent")) tf.byTaNCfrOnePercent=tau->tauID("byTaNCfrOnePercent");
    if (tau->isTauIDAvailable("byTaNCfrHalfPercent")) tf.byTaNCfrHalfPercent=tau->tauID("byTaNCfrHalfPercent");
    if (tau->isTauIDAvailable("byTaNCfrQuarterPercent")) tf.byTaNCfrQuarterPercent=tau->tauID("byTaNCfrQuarterPercent");
    if (tau->isTauIDAvailable("byTaNCfrTenthPercent")) tf.byTaNCfrTenthPercent=tau->tauID("byTaNCfrTenthPercent");
    if (tau->isTauIDAvailable("byTaNC")) tf.byTaNC=tau->tauID("byTaNC");
    if (tau->isTauIDAvailable("byLooseCombinedIsolationDeltaBetaCorr")) tf.byLooseCombinedIsolationDeltaBetaCorr=tau->tauID("byLooseCombinedIsolationDeltaBetaCorr");
    if (tau->isTauIDAvailable("againstElectronMVA")) tf.againstElectronMVA=tau->tauID("againstElectronMVA");
    if (tau->isPFTau()) {
      tf.isolationPFChargedHadrCandsPtSum = tau->isolationPFChargedHadrCandsPtSum();
      tf.isolationPFGammaCandsEtSum = tau->isolationPFGammaCandsEtSum();
      if (tau->leadPFChargedHadrCand().isAvailable()) tf.leadPFChargedHadrCandPt = tau->leadPFChargedHadrCand()->pt();
      tf.NsignalPFChargedHadrCands = tau->signalPFChargedHadrCands().size();
      tf.NsignalPFGammaCands = tau->signalPFGammaCands().size();
    }
    hbbInfo->tauInfo.push_back(tf);
    if (verbose_) {
      std::cout << "SCZ DEBUG: againstElectronLoose is " << tf.againstElectronLoose << std::endl;
      std::cout << "SCZ DEBUG: againstElectronMedium is " << tf.againstElectronMedium << std::endl;
      std::cout << "SCZ DEBUG: againstElectronTight is " << tf.againstElectronTight << std::endl;
      std::cout << "SCZ DEBUG: againstMuonLoose is " << tf.againstMuonLoose << std::endl;
      std::cout << "SCZ DEBUG: againstMuonTight is " << tf.againstMuonTight << std::endl;
      std::cout << "SCZ DEBUG: byLooseIsolation is " << tf.byLooseIsolation << std::endl;
      std::cout << "SCZ DEBUG: byMediumIsolation is " << tf.byMediumIsolation << std::endl;
      std::cout << "SCZ DEBUG: byTightIsolation is " << tf.byTightIsolation << std::endl;
      std::cout << "SCZ DEBUG: byVLooseIsolation is " << tf.byVLooseIsolation << std::endl;
      std::cout << "SCZ DEBUG: decayModeFinding is " << tf.decayModeFinding << std::endl;
      std::cout << "SCZ DEBUG: byIsolation is " << tf.byIsolation<< std::endl;
      std::cout << "SCZ DEBUG: trackIsolation is " << tf.trackIsolation << std::endl;
      std::cout << "SCZ DEBUG: byTaNCfrOnePercent is " << tf.byTaNCfrOnePercent << std::endl;
      std::cout << "SCZ DEBUG: byTaNCfrHalfPercent is " << tf.byTaNCfrHalfPercent << std::endl;
      std::cout << "SCZ DEBUG: byTaNCfrQuarterPercent is " << tf.byTaNCfrQuarterPercent << std::endl;
      std::cout << "SCZ DEBUG: byTaNCfrTenthPercent is " << tf.byTaNCfrTenthPercent << std::endl;
      std::cout << "SCZ DEBUG: byTaNC is " << tf.byTaNC << std::endl;
      std::cout << "SCZ DEBUG: isolationPFChargedHadrCandsPtSum is " << tf.isolationPFChargedHadrCandsPtSum << std::endl;
      std::cout << "SCZ DEBUG: isolationPFGammaCandsEtSum is " << tf.isolationPFGammaCandsEtSum << std::endl;
      std::cout << "SCZ DEBUG: isolationPFGammaCandsEtSum is " << tf.leadPFChargedHadrCandPt << std::endl;
      std::cout << "SCZ DEBUG: NsignalPFChargedHadrCands is " << tf.NsignalPFChargedHadrCands << std::endl;
      std::cout << "SCZ DEBUG: NsignalPFGammaCands is " << tf.NsignalPFGammaCands << std::endl;
      std::cout << "SCZ DEBUG: byLooseCombinedIsolationDeltaBetaCorr is " << tf.byLooseCombinedIsolationDeltaBetaCorr << std::endl;
      std::cout << "SCZ DEBUG: againstElectronMVA is " << tf.againstElectronMVA << std::endl;
    }
  }

  CompareJetPtMuons ptComparatorMu;
  CompareJetPtElectrons ptComparatorE;
  CompareJetPtTaus ptComparatorTau;

  std::sort(hbbInfo->muInfo.begin(), hbbInfo->muInfo.end(), ptComparatorMu);
  std::sort(hbbInfo->eleInfo.begin(), hbbInfo->eleInfo.end(), ptComparatorE);
  std::sort(hbbInfo->tauInfo.begin(), hbbInfo->tauInfo.end(), ptComparatorTau);


    
  
   if (verbose_){
     std::cout <<" Pushing hbbInfo "<<std::endl;
     std::cout <<" SimpleJets1 = "<<hbbInfo->simpleJets.size()<<std::endl<<
       " SimpleJets2 = "<<hbbInfo->simpleJets2.size()<<std::endl<<
       " SubJets = "<<hbbInfo->subJets.size()<<std::endl<<
       " HardJets = "<<hbbInfo->hardJets.size()<<std::endl<<
       " FilterJets = "<<hbbInfo->filterJets.size()<<std::endl<<
       " Muons = "<<hbbInfo->muInfo.size()<<std::endl<<
       " Electrons = "<<hbbInfo->eleInfo.size()<<std::endl<<
       " Taus = "<<hbbInfo->tauInfo.size()<<std::endl<<
       " Electrons = "<<hbbInfo->eleInfo.size()<<std::endl<<
       "--------------------- "<<std::endl;
  }


   //  iEvent.put(hbbInfo);
   //  iEvent.put(auxInfo);


   //  delete jecUnc;

}
  
void
HbbAnalyzerAlgo::fillMuBlock(std::vector<pat::Muon>::const_iterator mu, int *muInfo)
{
  if(muon::isGoodMuon(*mu,muon::TMLastStationLoose)) muInfo[0]=1;
  if(muon::isGoodMuon(*mu,muon::TMLastStationTight)) muInfo[1]=1;
  if(muon::isGoodMuon(*mu,muon::TM2DCompatibilityLoose)) muInfo[2]=1;
  if(muon::isGoodMuon(*mu,muon::TM2DCompatibilityTight)) muInfo[3]=1;
  if(muon::isGoodMuon(*mu,muon::TMOneStationLoose)) muInfo[4]=1;
  if(muon::isGoodMuon(*mu,muon::TMOneStationTight)) muInfo[5]=1;
  if(muon::isGoodMuon(*mu,muon::TMLastStationOptimizedLowPtLoose)) muInfo[6]=1;
  if(muon::isGoodMuon(*mu,muon::TMLastStationOptimizedLowPtTight))muInfo[7]=1;
  if(muon::isGoodMuon(*mu,muon::TMOneStationAngLoose)) muInfo[8]=1;
  if(muon::isGoodMuon(*mu,muon::TMOneStationAngTight)) muInfo[9]=1;
  if(muon::isGoodMuon(*mu,muon::TMLastStationAngLoose)) muInfo[10]=1;
  if(muon::isGoodMuon(*mu,muon::TMLastStationAngTight)) muInfo[11]=1;
  if(muon::isGoodMuon(*mu,muon::GMTkChiCompatibility)) muInfo[12]=1;
  if(muon::isGoodMuon(*mu,muon::GMStaChiCompatibility)) muInfo[13]=1;
  if(muon::isGoodMuon(*mu,muon::GMTkKinkTight)) muInfo[14]=1;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HbbAnalyzerAlgo::beginJob(){
}


// ------------ method called once each job just after ending the event loop  ------------
void 
HbbAnalyzerAlgo::endJob() {
}

TVector2 HbbAnalyzerAlgo::getTvect( const pat::Jet* patJet ){

  TVector2 t_Vect(0,0);
  TVector2 null(0,0);
  TVector2 ci(0,0);
  TLorentzVector pi(0,0,0,0);
  TLorentzVector J(0,0,0,0);
  TVector2 r(0,0);
  double patJetpfcPt = 1e10;
  double r_mag = 1e10;
  unsigned int nOfconst = 0;


  //  std::cout <<" OOOOOOOOOOO"<<std::endl;


  //re-reconstruct the jet direction with the charged tracks
  std::vector<reco::CandidatePtr>
    patJetpfc = patJet->getJetConstituents();
  for(size_t idx = 0; idx < patJetpfc.size(); idx++){
    if( patJetpfc.at(idx)->charge() != 0 ){
      pi.SetPtEtaPhiE( patJetpfc.at(idx)->pt(), patJetpfc.at(idx)->eta(), patJetpfc.at(idx)->phi(), patJetpfc.at(idx)->energy() );
      J += pi;
      nOfconst++;
    }
  }
// if there are less than two charged tracks do not calculate the pull (there is not enough info). It returns a null vector

  if( nOfconst < 2 )
    return null;
  


  TVector2 v_J( J.Rapidity(), J.Phi() );
//calculate TVector using only charged tracks
  for(size_t idx = 0; idx < patJetpfc.size(); idx++){
    if( patJetpfc.at(idx)->charge() != 0  ){
      patJetpfcPt = patJetpfc.at(idx)->pt();
      pi.SetPtEtaPhiE( patJetpfc.at(idx)->pt(), patJetpfc.at(idx)->eta(), patJetpfc.at(idx)->phi(), patJetpfc.at(idx)->energy() );
      //TOM WHAT TO DO 
      r.Set( pi.Rapidity() - J.Rapidity(), Geom::deltaBarePhi( patJetpfc.at(idx)->phi(), J.Phi() ) );
      r_mag = r.Mod();
      t_Vect += ( patJetpfcPt / J.Pt() ) * r_mag * r;
    }
  }

  
  return t_Vect;
  
}

TLorentzVector HbbAnalyzerAlgo::getChargedTracksMomentum(const pat::Jet* patJet ){

  // std::cout <<"fff ECCCCCCOOOOO "<<patJet->isPFJet()<<std::endl;
  

  //  return TLorentzVector();
  TLorentzVector pi(0,0,0,0);
  TLorentzVector v_j1(0,0,0,0);



  if (patJet->isPFJet() == false ){
      v_j1 = GENPTOLORP(patJet);
      return v_j1;
  }

  //  std::cout <<" DDDDD" <<std::endl;

  std::vector<reco::CandidatePtr>
    j1pfc = patJet->daughterPtrVector();
  //  std::cout<<"SSSS"<<std::endl;
  //  std::cout<<"SSSS "<<j1pfc.size()<<std::endl;
  for(size_t idx = 0; idx < j1pfc.size(); idx++){
    //    std::cout <<" AT "<<idx<<std::endl;
    if( j1pfc.at(idx)->charge() != 0 ){
      pi.SetPtEtaPhiE( j1pfc.at(idx)->pt(), j1pfc.at(idx)->eta(), j1pfc.at(idx)->phi(), j1pfc.at(idx)->energy() );
      v_j1 += pi;
    }
  }
  return v_j1;
  //re-
}


//Btagging scale factors
void HbbAnalyzerAlgo::fillScaleFactors(VHbbEvent::SimpleJet& sj, BTagSFContainer iSF){


  BinningPointByMap measurePoint;
  //for a USDG
  //for CB jets
  //scale factor 1 for CB jets over 240GeV/c
  if( TMath::Abs(sj.flavour) == 4 or TMath::Abs(sj.flavour) == 5 ){
    measurePoint.insert( BinningVariables::JetEt, sj.p4.Et() );
    measurePoint.insert( BinningVariables::JetAbsEta, fabs(sj.p4.Eta()) );
    if( iSF.BTAGSF_CSVL->isResultOk(PerformanceResult::BTAGBEFFCORR , measurePoint) ){
      sj.SF_CSVL = iSF.BTAGSF_CSVL->getResult(PerformanceResult::BTAGBEFFCORR , measurePoint);
      sj.SF_CSVLerr = iSF.BTAGSF_CSVL->getResult(PerformanceResult::BTAGBERRCORR , measurePoint);	  
      if(verbose_){
	std::clog << "C/B Jet flavour = " << sj.flavour << std::endl;
	std::clog << "C/B Jet Et = " << sj.p4.Et() << std::endl;
	std::clog << "C/B Jet eta = " << sj.p4.Eta() << std::endl;	    
	std::clog << "C/B CSVL Scale Factor = " << sj.SF_CSVL << std::endl; 
	std::clog << "C/B CSVL Scale Factor error = " << sj.SF_CSVLerr << std::endl; 
      }
    }
    if( iSF.BTAGSF_CSVM->isResultOk(PerformanceResult::BTAGBEFFCORR , measurePoint) ){
      sj.SF_CSVM = iSF.BTAGSF_CSVM->getResult(PerformanceResult::BTAGBEFFCORR , measurePoint);
      sj.SF_CSVMerr = iSF.BTAGSF_CSVM->getResult(PerformanceResult::BTAGBERRCORR , measurePoint);	  
    }
    if( iSF.BTAGSF_CSVT->isResultOk(PerformanceResult::BTAGBEFFCORR , measurePoint) ){
      sj.SF_CSVT = iSF.BTAGSF_CSVT->getResult(PerformanceResult::BTAGBEFFCORR , measurePoint);
      sj.SF_CSVTerr = iSF.BTAGSF_CSVT->getResult(PerformanceResult::BTAGBERRCORR , measurePoint);	  
    }
    else{
      if(verbose_){
	std::cerr << "No SF found in the database for this jet" << std::endl;
	std::clog << "No SF found: Jet flavour = " << sj.flavour << std::endl;
	std::clog << "No SF found: Jet Et = " << sj.p4.Et() << std::endl;
	std::clog << "No SF found: Jet eta = " << sj.p4.Eta() << std::endl;
      }
    }
  }
  else {
    measurePoint.insert( BinningVariables::JetEt, sj.p4.Et() );
    measurePoint.insert( BinningVariables::JetAbsEta, fabs(sj.p4.Eta()) );
    if( iSF.MISTAGSF_CSVL->isResultOk(PerformanceResult::BTAGLEFFCORR , measurePoint) ){
      sj.SF_CSVL = iSF.MISTAGSF_CSVL->getResult(PerformanceResult::BTAGLEFFCORR , measurePoint);
      sj.SF_CSVLerr = iSF.MISTAGSF_CSVL->getResult(PerformanceResult::BTAGLERRCORR , measurePoint);
      if(verbose_){
	std::clog << "Light Jet flavour = " << sj.flavour << std::endl;
	std::clog << "Light Jet Et = " << sj.p4.Et() << std::endl;
	std::clog << "Light Jet eta = " << sj.p4.Eta() << std::endl;	    
	std::clog << "Light CSVL Scale Factor = " << sj.SF_CSVL << std::endl; 
	std::clog << "Light CSVL Scale Factor error = " << sj.SF_CSVLerr << std::endl; 
      }
    }
    if( iSF.MISTAGSF_CSVM->isResultOk(PerformanceResult::BTAGLEFFCORR , measurePoint) ){
      sj.SF_CSVM = iSF.MISTAGSF_CSVM->getResult(PerformanceResult::BTAGLEFFCORR , measurePoint);
      sj.SF_CSVMerr = iSF.MISTAGSF_CSVM->getResult(PerformanceResult::BTAGLERRCORR , measurePoint);
    }
    if( iSF.MISTAGSF_CSVT->isResultOk(PerformanceResult::BTAGLEFFCORR , measurePoint) ){
      sj.SF_CSVT = iSF.MISTAGSF_CSVT->getResult(PerformanceResult::BTAGLEFFCORR , measurePoint);
      sj.SF_CSVTerr = iSF.MISTAGSF_CSVT->getResult(PerformanceResult::BTAGLERRCORR , measurePoint);
    }
    else{
      if(verbose_){
	std::cerr << "No SF found in the database for this jet" << std::endl;
	std::clog << "No SF found: Jet flavour = " << sj.flavour << std::endl;
	std::clog << "No SF found: Jet Et = " << sj.p4.Et() << std::endl;
	std::clog << "No SF found: Jet eta = " << sj.p4.Eta() << std::endl;
      }
    }
  }
    
}

void HbbAnalyzerAlgo::setJecUnc(VHbbEvent::SimpleJet& sj,JetCorrectionUncertainty* jecunc){
  //
  // test 
  //

  //  return;
  double eta = sj.p4.Eta();
  double pt = sj.p4.Pt();
  
  jecunc->setJetEta(eta);
  jecunc->setJetPt(pt); // here you must use the CORRECTED jet pt
  double unc = jecunc->getUncertainty(true);
  sj.jecunc= unc;
}


void HbbAnalyzerAlgo ::fillSimpleJet (VHbbEvent::SimpleJet& sj, std::vector<pat::Jet>::const_iterator jet_iter){
      sj.flavour = jet_iter->partonFlavour();
    sj.tche=jet_iter->bDiscriminator("trackCountingHighEffBJetTags");
    sj.tchp=jet_iter->bDiscriminator("trackCountingHighPurBJetTags");
    sj.jp=jet_iter->bDiscriminator("jetProbabilityBJetTags");
    sj.jpb=jet_iter->bDiscriminator("jetBProbabilityBJetTags");
    sj.ssvhe=jet_iter->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    sj.csv=jet_iter->bDiscriminator("combinedSecondaryVertexBJetTags");
    sj.csvmva=jet_iter->bDiscriminator("combinedSecondaryVertexMVABJetTags");
    sj.csvivf=jet_iter->bDiscriminator("combinedInclusiveSecondaryVertexBJetTags");
    sj.cmva=jet_iter->bDiscriminator("combinedMVABJetTags");


    sj.charge=jet_iter->jetCharge();
    sj.ntracks=jet_iter->associatedTracks().size();
    sj.p4=GENPTOLORP(jet_iter);
    //    std::cout << " ECCO "<<sj.csv<< " "<< sj.p4.Pt()<<std::endl;
         sj.chargedTracksFourMomentum=(getChargedTracksMomentum(&*(jet_iter)));
	 //    std::cout <<" OUT"<<std::endl;
    sj.SF_CSVL=1;
    sj.SF_CSVM=1;
    sj.SF_CSVT=1;
    sj.SF_CSVLerr=0;
    sj.SF_CSVMerr=0;
    sj.SF_CSVTerr=0;

    //    std::cout <<" PIPPO"<<std::endl;

    //    std::cout <<" TOMMMMMMM "<<jet_iter->partonFlavour()<<std::endl;

    /* TOM

    //for quark-gluon tagger
    sj.constituentPtDistribution = jet_iter->constituentPtDistribution();
    std::cout <<" PIPPOasdasd"<<std::endl;
    sj.constituentEtaPhiSpread = jet_iter->constituentEtaPhiSpread(); 
    std::cout <<" PIPPO"<<std::endl;
    
    TOM */

    if (jet_iter->isPFJet() == true) {
            sj.chargedHadronEFraction = jet_iter-> chargedHadronEnergyFraction();
      sj.neutralHadronEFraction = jet_iter-> neutralHadronEnergyFraction ();

      sj.chargedEmEFraction = jet_iter-> chargedEmEnergyFraction ();
      sj.neutralEmEFraction = jet_iter-> neutralEmEnergyFraction ();


      sj.nConstituents = jet_iter->daughterPtrVector().size();
      
    }
    //    std::cout <<" PIPPO"<<std::endl;
    sj.jetArea = jet_iter->jetArea();
    //
    // addtaginfo for csv
    //

    //    if (jet_iter->hasTagInfo("SimpleSecondaryVertex")) {
    //    std::cout <<" PIPPOff"<<std::endl;
    const reco::SecondaryVertexTagInfo * tf = jet_iter->tagInfoSecondaryVertex();
   if (tf){
      math::XYZTLorentzVectorD vertexSum;
      for(size_t vi=0;vi< tf->nVertices();vi++)
      {
        vertexSum+=tf->secondaryVertex(vi).p4();
      }
      sj.vtxP4 = GENPTOLOR(vertexSum);

     if (tf->nVertices() >0){
	sj.vtxPosition = TVector3(tf->secondaryVertex(0).position().x(),tf->secondaryVertex(0).position().y(),tf->secondaryVertex(0).position().z());
	sj.vtxMass =  tf->secondaryVertex(0).p4().mass();
	sj.vtxNTracks = tf->secondaryVertex(0).nTracks();
	std::vector<reco::TrackBaseRef >::const_iterator tit =  tf->secondaryVertex(0).tracks_begin();
	for (; tit<  tf->secondaryVertex(0).tracks_end(); ++tit){
	  sj.vtxTrackIds.push_back(tit->key());
	}
	Measurement1D m = tf->flightDistance(0);
	sj.vtx3dL = m.value();
	sj.vtx3deL = m.error();
     }
    }
    
   //   std::cout <<" SONO QUI"<<std::endl;
   /* TOM

   // CSV track info
   const reco::SecondaryVertexTagInfo * svTagInfos = jet_iter->tagInfoSecondaryVertex();
    const reco::TrackIPTagInfo * ipTagInfos = jet_iter->tagInfoTrackIP();
    for (edm::RefVector<reco::TrackCollection>::const_iterator t = ipTagInfos->selectedTracks().begin(); t != ipTagInfos->selectedTracks().end(); t++){
      sj.btagTrackIds.push_back(t->key());
    }// all btag IP selected tracks    
    std::vector<const reco::BaseTagInfo*> tagInfos;
    tagInfos.push_back(dynamic_cast<const reco::BaseTagInfo*>(ipTagInfos));
    tagInfos.push_back(dynamic_cast<const reco::BaseTagInfo*>(svTagInfos));
    JetTagComputer::TagInfoHelper helper(tagInfos);
    reco::TaggingVariableList varList = computer->taggingVariables(helper); // computer for getting CSV variables



      
    for(reco::TaggingVariableList::const_iterator iter = varList.begin(); iter != varList.end(); ++iter)
    {
      //std::cout << reco::TaggingVariableTokens[iter->first] << " = " << iter->second << std::endl;
      for (edm::RefVector<reco::TrackCollection>::const_iterator t = ipTagInfos->selectedTracks().begin(); t != ipTagInfos->selectedTracks().end(); t++){
        
        if (strcmp(reco::TaggingVariableTokens[iter->first], "trackMomentum") == 0 && (fabs((float)iter->second - (float)(*t)->p()) < 0.0001) ){
          sj.csvTrackIds.push_back(t->key());
        }// if tagged track
      }// loop on IPtracks        
    }// loop on CSV variables

    
    sj.btagNTracks= ipTagInfos->selectedTracks().size();
    sj.csvNTracks = sj.csvTrackIds.size();
 TOM */
   //
    // add tVector
    //
    sj.tVector = getTvect(&(*jet_iter));

    //    std::cout <<"DFGHJKL"<<std::endl;

    sj.ptRaw = jet_iter->correctedJet(0).pt();
    //    std::cout <<" OOOOO"<<std::endl;
    sj.ptLeadTrack =-9999.;
    if (jet_iter->isPFJet() == true) {
       std::vector <reco::CandidatePtr> constituents = jet_iter->getJetConstituents ();
       for (unsigned ic = 0; ic < constituents.size (); ++ic) {
	 //	 std::cout <<" SJSJSJS "<<std::endl;
	 //        if ( constituents[ic]->particleId() > 3 ) continue;
	 //         reco::TrackRef trackRef = ->trackRef();
	 if (constituents[ic]->charge()!=0  ) { if(constituents[ic]->pt() > sj.ptLeadTrack) sj.ptLeadTrack=constituents[ic]->pt(); }
      }
     }


}

float HbbAnalyzerAlgo::metSignificance(const reco::MET * met)
{
double sigmaX2= met->getSignificanceMatrix()(0,0);
double sigmaY2= met->getSignificanceMatrix()(1,1);
double significance = 0;
try {
if(sigmaX2<1.e10 && sigmaY2<1.e10) significance = met->significance();
}
catch(...)
{
 std::cout << "PROBLEM WITH MET SIGNIFICANCE sigma X2 and Y2 are: " << sigmaX2 << " " << sigmaY2 << std::endl;
}
return significance;
}


