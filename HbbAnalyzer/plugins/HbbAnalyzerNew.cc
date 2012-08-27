// -*- C++ -*-
//
// Package:    HbbAnalyzerNew
// Class:      HbbAnalyzerNew
// 
/**\class HbbAnalyzerNew HbbAnalyzerNew.cc Analysis/HbbAnalyzer/src/HbbAnalyzerNew.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  David Lopes Pegna,Address unknown,NONE,
//         Created:  Thu Mar  5 13:51:28 EST 2009
// $Id: HbbAnalyzerNew.cc,v 1.76 2012/06/12 19:40:16 arizzi Exp $
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

#include "VHbbAnalysis/HbbAnalyzer/interface/HbbAnalyzerNew.h"
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

#include "CMGTools/External/interface/PileupJetIdentifier.h"
#include "CMGTools/External/interface/PileupJetIdAlgo.h"
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



HbbAnalyzerNew::HbbAnalyzerNew(const edm::ParameterSet& iConfig):
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

  //
  // put the setwhatproduced etc etc
 
  produces<VHbbEvent>();
  produces<VHbbEventAuxInfo>();

 
}


HbbAnalyzerNew::~HbbAnalyzerNew(){
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HbbAnalyzerNew::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  using namespace reco;
 
 
  // JEC Uncertainty

  //  JetCorrectionUncertainty *jecUnc=0;
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get("AK5PFchs",JetCorParColl); 
  JetCorrectionUncertainty *jecUnc=0;
  //  if (!runOnMC_){
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  jecUnc = new JetCorrectionUncertainty(JetCorPar);
  //  }
  
  std::auto_ptr<VHbbEvent> hbbInfo( new VHbbEvent() );  
  std::auto_ptr<VHbbEventAuxInfo> auxInfo( new VHbbEventAuxInfo() );
 

  if (runOnMC_){
    Handle<GenEventInfoProduct> evt_info;
    iEvent.getByType(evt_info);
    auxInfo->weightMCProd = evt_info->weight(); 
  }
  else
    { auxInfo->weightMCProd =1.;}
  //
  // ??
  
  // trigger

  // trigger
  edm::Handle<edm::TriggerResults>  hltresults;
  //iEvent.getByLabel("TriggerResults", hltresults);
   
  //edm::InputTag tag("TriggerResults::HLT");
  //  edm::InputTag tag("TriggerResults::HLT0");
  iEvent.getByLabel(hltResults_, hltresults);
   
  const edm::TriggerNames & triggerNames_ = iEvent.triggerNames(*hltresults);

  int ntrigs = hltresults->size();
  if (ntrigs==0){std::cerr << "%HLTInfo -- No trigger name given in TriggerResults of the input " << std::endl;}

  BeamSpot vertexBeamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot",recoBeamSpotHandle);
  vertexBeamSpot = *recoBeamSpotHandle;
  /*
    double BSx=vertexBeamSpot.x0();
    double BSy=vertexBeamSpot.y0();
    double BSz=vertexBeamSpot.z0();
  */

  double MinVtxProb=-999.;
  int VtxIn=-99;
  
  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("offlinePrimaryVertices", recVtxs);
  
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
    
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(edm::InputTag("kt6PFJets", "rho"),rhoHandle);   
  auxInfo->puInfo.rho = *rhoHandle;
   
  edm::Handle<double> rho25Handle;
  iEvent.getByLabel(edm::InputTag("kt6PFJets25", "rho"),rho25Handle);   
  auxInfo->puInfo.rho25 = *rho25Handle;
  edm::Handle<double> rho25HandleIso;
  iEvent.getByLabel(edm::InputTag("kt6PFJetsForIsolation", "rho"),rho25HandleIso);   
  auxInfo->puInfo.rho25Iso = *rho25HandleIso;
 
  edm::Handle<double> rhoNeutralHandle;
  iEvent.getByLabel(edm::InputTag("kt6PFJetsCentralNeutral", "rho"),rhoNeutralHandle);   
  auxInfo->puInfo.rhoNeutral = *rhoNeutralHandle;


  edm::Handle<std::vector< PileupSummaryInfo> > puHandle;

  if (runOnMC_){
    iEvent.getByType(puHandle);
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
  
  
  Handle<GenParticleCollection> genParticles;
  
  bool printJet=0;
  
  
  if(runOnMC_){
    
   iEvent.getByLabel("genParticles", genParticles);
    
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
  edm::Handle<edm::View<reco::PFCandidate> > photonIsoH;
  iEvent.getByLabel("pfAllPhotons",photonIsoH);
  edm::View<reco::PFCandidate> photonsForIso = *photonIsoH;

  edm::Handle<edm::View<pat::Muon> > muonHandle;
  iEvent.getByLabel(muoLabel_,muonHandle);
  edm::View<pat::Muon> muons = *muonHandle;

  // hard jet   
  edm::Handle<edm::View<pat::Jet> > jetHandle;
  iEvent.getByLabel(jetLabel_,jetHandle);
  edm::View<pat::Jet> jets = *jetHandle;

  // sub jet   
  edm::Handle<edm::View<pat::Jet> > subjetHandle;
  iEvent.getByLabel(subjetLabel_,subjetHandle);
  edm::View<pat::Jet> subjets = *subjetHandle;

  // filter jet   
  edm::Handle<edm::View<pat::Jet> > filterjetHandle;
  iEvent.getByLabel(filterjetLabel_,filterjetHandle);
  edm::View<pat::Jet> filterjets = *filterjetHandle;

  // standard jets


  edm::Handle<edm::View<pat::Jet> > simplejet2Handle;
  iEvent.getByLabel(simplejet2Label_,simplejet2Handle);
  edm::View<pat::Jet> simplejets2 = *simplejet2Handle;

  edm::Handle<edm::View<pat::Jet> > simplejet3Handle;
  iEvent.getByLabel(simplejet3Label_,simplejet3Handle);
  edm::View<pat::Jet> simplejets3 = *simplejet3Handle;



  edm::Handle<edm::View<pat::Electron> > electronHandle;
  iEvent.getByLabel(eleLabel_,electronHandle);
  edm::View<pat::Electron> electrons = *electronHandle;


  //   edm::Handle<edm::View<pat::Photon> > phoHandle;
  //   iEvent.getByLabel(phoLabel_,phoHandle);
  //   edm::View<pat::Photon> photons = *phoHandle;

  edm::Handle<edm::View<pat::Tau> > tauHandle;
  iEvent.getByLabel(tauLabel_,tauHandle);
  edm::View<pat::Tau> taus = *tauHandle;
  
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

BTagSFContainer btagSFs;
  btagSFs.BTAGSF_CSVL = (bTagSF_CSVL_.product());
  btagSFs.BTAGSF_CSVM = (bTagSF_CSVM_.product());
  btagSFs.BTAGSF_CSVT = (bTagSF_CSVT_.product());
  btagSFs.MISTAGSF_CSVL = (mistagSF_CSVL_.product());
  btagSFs.MISTAGSF_CSVM = (mistagSF_CSVM_.product());
  btagSFs.MISTAGSF_CSVT = (mistagSF_CSVT_.product());

#ifdef ENABLE_SIMPLEJETS1
  edm::Handle<edm::View<pat::Jet> > simplejet1Handle;
  iEvent.getByLabel(simplejet1Label_,simplejet1Handle);
  edm::View<pat::Jet> simplejets1 = *simplejet1Handle;
  for(edm::View<pat::Jet>::const_iterator jet_iter = simplejets1.begin(); jet_iter!=simplejets1.end(); ++jet_iter){
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

  for(edm::View<pat::Jet>::const_iterator jet_iter = simplejets3.begin(); jet_iter!=simplejets3.end(); ++jet_iter){
    //     if(jet_iter->pt()>50)
    //       njetscounter++;
    VHbbEvent::SimpleJet sj;
    //    std::cout <<" sj3"<<std::endl;
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
    // 


    hbbInfo->simpleJets3.push_back(sj);
    
  }

#ifdef ENABLE_SIMPLEJETS4
  edm::Handle<edm::View<pat::Jet> > simplejet4Handle;
  iEvent.getByLabel(simplejet4Label_,simplejet4Handle);
  edm::View<pat::Jet> simplejets4 = *simplejet4Handle;
  for(edm::View<pat::Jet>::const_iterator jet_iter = simplejets4.begin(); jet_iter!=simplejets4.end(); ++jet_iter){
    //     if(jet_iter->pt()>50)
    //       njetscounter++;
  
   VHbbEvent::SimpleJet sj;
    //    std::cout <<" sj4"<<std::endl;
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
    hbbInfo->simpleJets4.push_back(sj);
    
  }
#endif //ENABLE SIMPLEJETS4


  for(edm::View<pat::Jet>::const_iterator jet_iter = simplejets2.begin(); jet_iter!=simplejets2.end(); ++jet_iter){
    


    VHbbEvent::SimpleJet sj;
    //    std::cout <<" sj2"<<std::endl;
    fillSimpleJet(sj,jet_iter);  

    ///###########          PU JET ID #################
   // add puId...
    edm::Handle<edm::ValueMap<float> > puJetIdMVA;
    iEvent.getByLabel("puJetMva","fullDiscriminant", puJetIdMVA);

    edm::Handle<edm::ValueMap<int> > puJetIdFlag;
    iEvent.getByLabel("puJetMva", "fullId", puJetIdFlag);

    //    cout  << " pt " << jet_iter->pt() << " eta " << jet_iter->eta() << std::endl;
    unsigned int idx = jet_iter - simplejets2.begin();



    sj.puJetIdMva   = (*puJetIdMVA)[simplejets2.refAt(idx)];
    int    idflag = (*puJetIdFlag)[simplejets2.refAt(idx)];
  
    
    //     cout << " PU JetID MVA " << mva; 
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

  
    //  if(!runOnMC_)  
 setJecUnc(sj,jecUnc);
    /*    sj.flavour = jet_iter->partonFlavour();
    
    
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
    }


    //
    // add tVector
    //
    sj.tVector = getTvect(&(*jet_iter));
    */
    Particle::LorentzVector p4Jet = jet_iter->p4();

    if(runOnMC_){

      //BTV scale factors
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
	  std::clog << "genJet matched deltaR = " << gJp4.DeltaR(sj.p4) << std::endl;
	  std::clog << "genJet matched mother id = " << sj.bestMCmomid << std::endl;
	}
      }

        // add flag if a mc lepton is find inside a cone around the jets... 
      iEvent.getByLabel("genParticles", genParticles);
      
      for(size_t i = 0; i < genParticles->size(); ++ i) {
     
      const GenParticle & p = (*genParticles)[i];
      int id = 0; 
      p.pt()> lep_ptCutForBjets_ ? id= p.pdgId(): 0;
   
      //      std::cout<< "found a muon with pt " << mu->pt()   << std::endl;
      if   ((abs(id)==13 || abs(id)==11) && deltaR(p.eta(), p.phi(), sj.p4.Eta(), sj.p4.Phi() ) <0.5)  sj.isSemiLeptMCtruth=1;
      }

    }  //isMC

        // add flag if a reco lepton is find inside a cone around the jets... 
    edm::Handle<edm::View<reco::Candidate> > muonNoCutsHandle;
    iEvent.getByLabel(muonoCutsLabel_,muonNoCutsHandle);
    edm::View<reco::Candidate> muonsNoCuts = *muonNoCutsHandle; 
    
    

    for(edm::View<reco::Candidate>::const_iterator mu = muonsNoCuts.begin(); mu!=muonsNoCuts.end() && sj.isSemiLept!=1; ++mu){
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
    
    
    edm::Handle<edm::View<reco::Candidate> > eleNoCutsHandle;
      iEvent.getByLabel(elenoCutsLabel_,eleNoCutsHandle);
      edm::View<reco::Candidate> elesNoCuts = *eleNoCutsHandle; 
 
    

      for(edm::View<reco::Candidate>::const_iterator ele = elesNoCuts.begin(); ele!=elesNoCuts.end() && sj.isSemiLept!=1; ++ele){
    
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
   


  /////// hard jet

  
  double matEta[1000*30],matPhi[1000*30];
  for(int i=0;i<1000;i++){for(int j=0;j<30;j++){matEta[i*j]=-99.;matPhi[i*j]=-99.;}}

  for(edm::View<pat::Jet>::const_iterator jet_iter = jets.begin(); jet_iter!=jets.end(); ++jet_iter){
    
    if(printJet) {std::cout << "Jet Pt: " << jet_iter->pt() << " E,M: " << jet_iter->p4().E() << " " << jet_iter->p4().M() << "\n";} 
    
    reco::Jet::Constituents constituents = jet_iter->getJetConstituents();
    
    //    if(printJet) {std::cout << "NsubJets: " << constituents.size() << "\n";} 
    VHbbEvent::HardJet hj;
    hj.constituents=constituents.size();
    hj.p4 =GENPTOLORP(jet_iter);
    
    for (unsigned int iJC(0); iJC<constituents.size(); ++iJC ){
      Jet::Constituent icandJet = constituents[iJC];

      if(printJet) {std::cout << "subJet Pt: " << icandJet->pt() << " subJet E,M,eta,phi: " <<  icandJet->p4().E() << "," 
			      << icandJet->p4().M() << "," << icandJet->eta() << "," << icandJet->phi() << "\n"; } 


      hj.subFourMomentum.push_back(GENPTOLORP(icandJet));
      hj.etaSub.push_back(icandJet->eta());
      hj.phiSub.push_back(icandJet->phi());
     
    }
    hbbInfo->hardJets.push_back(hj);
  }
    
  //  HardJetSubEta2.SetMatrixArray(matEta);
  //  HardJetSubPhi2.SetMatrixArray(matPhi);
  //  TMatrixDRow a1(HardJetSubEta2,0);
  //  for(int i=0;i<30;i++){
  //   std::cout << "test: " << a1[i] << "\n";  
  //  }

  // pat subJets with Btag


  for(edm::View<pat::Jet>::const_iterator subjet_iter = subjets.begin(); subjet_iter!=subjets.end(); ++subjet_iter){

    if(printJet) {std::cout << "SubJetTagged Pt: " << subjet_iter->pt() << " E,M,eta,phi,Btag: " << subjet_iter->p4().E() 
			    << "," << subjet_iter->p4().M() << "," << subjet_iter->eta() << "," << subjet_iter->phi()  
			    << "," << subjet_iter->bDiscriminator("combinedSecondaryVertexBJetTags") << "\n";}

    VHbbEvent::SimpleJet sj;
    //    std::cout <<" sub jet "<<std::endl;
    fillSimpleJet(sj,subjet_iter);
    //  if(!runOnMC_)  
    setJecUnc(sj,jecUnc);
    /*    sj.flavour = subjet_iter->partonFlavour();
    sj.tVector = getTvect(&(*subjet_iter));
    sj.tche=subjet_iter->bDiscriminator("trackCountingHighEffBJetTags");
    sj.tchp=subjet_iter->bDiscriminator("trackCountingHighPurBJetTags");
    sj.jp=subjet_iter->bDiscriminator("jetProbabilityBJetTags");
    sj.jpb=subjet_iter->bDiscriminator("jetBProbabilityBJetTags");
    sj.ssvhe=subjet_iter->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    sj.csv=subjet_iter->bDiscriminator("combinedSecondaryVertexBJetTags");
    sj.csvmva=subjet_iter->bDiscriminator("combinedSecondaryVertexMVABJetTags");
    sj.charge=subjet_iter->jetCharge();
    sj.ntracks=subjet_iter->associatedTracks().size();
    sj.p4=GENPTOLORP(subjet_iter);
    sj.p4=(getChargedTracksMomentum(&*(subjet_iter)));

    //
    // addtaginfo for csv
    //

    if (subjet_iter->hasTagInfo("SimpleSecondaryVertex")) {

      const reco::SecondaryVertexTagInfo * tf = subjet_iter->tagInfoSecondaryVertex();
      sj.vtxMass = tf->secondaryVertex(0).p4().mass();
      sj.vtxNTracks = tf->secondaryVertex(0).nTracks();
      Measurement1D m = tf->flightDistance(0);
      sj.vtx3dL = m.value();
      sj.vtx3deL = m.error();
    }
    */
    hbbInfo->subJets.push_back(sj);

  }

  for(edm::View<pat::Jet>::const_iterator filterjet_iter = filterjets.begin(); filterjet_iter!=filterjets.end(); ++filterjet_iter){
 
    if(printJet) {std::cout << "FilterjetTagged Pt: " << filterjet_iter->pt() << " E,M,eta,phi,Btag: " << filterjet_iter->p4().E()  << "," << filterjet_iter->p4().M() << "," << filterjet_iter->eta() << "," << filterjet_iter->phi()  << "," << filterjet_iter->bDiscriminator("combinedSecondaryVertexBJetTags") << "\n";}
 
    VHbbEvent::SimpleJet fj;
    //    std::cout <<" sub jet "<<std::endl;
    fillSimpleJet(fj,filterjet_iter);
    //  if(!runOnMC_)  
    setJecUnc(fj,jecUnc);

    if(runOnMC_){

      //BTV scale factors
     // fillScaleFactors(sj, btagSFs);

      //PAT genJet matching
      //genJet
      const reco::GenJet *gJ = filterjet_iter->genJet();
      //physical parton for mother info ONLY
      if( (filterjet_iter->genParton()) ){
        fj.bestMCid = filterjet_iter->genParton()->pdgId();
        if( (filterjet_iter->genParton()->mother()) )
          fj.bestMCmomid=filterjet_iter->genParton()->mother()->pdgId();
      }
      TLorentzVector gJp4;
      if(gJ){
        gJp4.SetPtEtaPhiE(gJ->pt(),gJ->eta(),gJ->phi(),gJ->energy());
        fj.bestMCp4 = gJp4;
        if(verbose_){
          std::clog << "filter genJet matched Pt = " << gJp4.Pt() << std::endl;
          std::clog << "filter genJet matched eta = " << gJp4.Eta() << std::endl;
          std::clog << "filter genJet matched deltaR = " << gJp4.DeltaR(fj.p4) << std::endl;
          std::clog << "filter genJet matched mother id = " << fj.bestMCmomid << std::endl;
        }
      }
    }

    hbbInfo->filterJets.push_back(fj);
    

  }

  //
  // add charged met
  //
  
  edm::Handle<edm::View<reco::MET> > metChargedHandle;
  iEvent.getByLabel("pfMETNoPUCharge",metChargedHandle);
  edm::View<reco::MET> metsCh = *metChargedHandle;
  if(metsCh.size()){
    hbbInfo->metCh.sumEt=(metsCh[0]).sumEt();
    hbbInfo->metCh.metSig=metSignificance(& (metsCh[0]));
    hbbInfo->metCh.eLong=(metsCh[0]).e_longitudinal();
    hbbInfo->metCh.p4=GENPTOLOR((metsCh[0]));
    if (verbose_)     std::cout <<" METCharged "<<     hbbInfo->metCh.metSig <<" " <<     hbbInfo->metCh.sumEt<<std::endl;
  }

  // type 1 corr met
  edm::Handle<edm::View<reco::MET> > pfmetType1corrHandle;
  iEvent.getByLabel("patType1CorrectedPFMet",pfmetType1corrHandle);
  edm::View<reco::MET> pfmetsType1corr = *pfmetType1corrHandle;
  if(pfmetsType1corr.size()){
    hbbInfo->pfmetType1corr.sumEt=(pfmetsType1corr[0]).sumEt();
    hbbInfo->pfmetType1corr.metSig=metSignificance(& (pfmetsType1corr[0]));
    hbbInfo->pfmetType1corr.eLong=(pfmetsType1corr[0]).e_longitudinal();
    hbbInfo->pfmetType1corr.p4=GENPTOLOR((pfmetsType1corr[0]));
    if (verbose_)     std::cout <<" type 1 corrected pfMET "<<     hbbInfo->pfmetType1corr.metSig <<" " <<     hbbInfo->pfmetType1corr.sumEt<<std::endl;
  }


  // type 1 + 2 corr met
  edm::Handle<edm::View<reco::MET> > pfmetType1p2corrHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMet",pfmetType1p2corrHandle);
  edm::View<reco::MET> pfmetsType1p2corr = *pfmetType1p2corrHandle;
  if(pfmetsType1p2corr.size()){
    hbbInfo->pfmetType1p2corr.sumEt=(pfmetsType1p2corr[0]).sumEt();
    hbbInfo->pfmetType1p2corr.metSig=metSignificance(& (pfmetsType1p2corr[0]));
    hbbInfo->pfmetType1p2corr.eLong=(pfmetsType1p2corr[0]).e_longitudinal();
    hbbInfo->pfmetType1p2corr.p4=GENPTOLOR((pfmetsType1p2corr[0]));
    if (verbose_)     std::cout <<" type 1 +2 corrected pfMET "<<     hbbInfo->pfmetType1p2corr.metSig <<" " <<     hbbInfo->pfmetType1p2corr.sumEt<<std::endl;
  }

  // type 1 corr met NoPU
/*  edm::Handle<edm::View<reco::MET> > pfmetNoPUType1corrHandle;
  iEvent.getByLabel("patType1CorrectedPFMetNoPU",pfmetNoPUType1corrHandle);
  edm::View<reco::MET> pfmetsNoPUType1corr = *pfmetNoPUType1corrHandle;
  if(pfmetsNoPUType1corr.size()){
    hbbInfo->pfmetNoPUType1corr.sumEt=(pfmetsNoPUType1corr[0]).sumEt();
    hbbInfo->pfmetNoPUType1corr.metSig=metSignificance(& (pfmetsNoPUType1corr[0]));
    hbbInfo->pfmetNoPUType1corr.eLong=(pfmetsNoPUType1corr[0]).e_longitudinal();
    hbbInfo->pfmetNoPUType1corr.p4=GENPTOLOR((pfmetsNoPUType1corr[0]));
    if (verbose_)     std::cout <<" type 1 corrected pfMET NoPU"<<     hbbInfo->pfmetNoPUType1corr.metSig <<" " <<     hbbInfo->pfmetNoPUType1corr.sumEt<<std::endl;
  }


  // type 1 + 2 corr met
  edm::Handle<edm::View<reco::MET> > pfmetNoPUType1p2corrHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMetNoPU",pfmetNoPUType1p2corrHandle);
  edm::View<reco::MET> pfmetsNoPUType1p2corr = *pfmetNoPUType1p2corrHandle;
  if(pfmetsNoPUType1p2corr.size()){
    hbbInfo->pfmetNoPUType1p2corr.sumEt=(pfmetsNoPUType1p2corr[0]).sumEt();
    hbbInfo->pfmetNoPUType1p2corr.metSig=metSignificance(& (pfmetsNoPUType1p2corr[0]));
    hbbInfo->pfmetNoPUType1p2corr.eLong=(pfmetsNoPUType1p2corr[0]).e_longitudinal();
    hbbInfo->pfmetNoPUType1p2corr.p4=GENPTOLOR((pfmetsNoPUType1p2corr[0]));
    if (verbose_)     std::cout <<" type 1 +2 corrected pfMET "<<     hbbInfo->pfmetNoPUType1p2corr.metSig <<" " <<     hbbInfo->pfmetNoPUType1p2corr.sumEt<<std::endl;
  }

*/

  /*
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
      */

    VHbbEvent::METInfo metunc;
  edm::Handle<edm::View<reco::MET> > patType1CorrectedPFMetElectronEnDownHandle;
  iEvent.getByLabel("patType1CorrectedPFMetElectronEnDown",patType1CorrectedPFMetElectronEnDownHandle);
  edm::View<reco::MET> patType1CorrectedPFMetsElectronEnDown = *patType1CorrectedPFMetElectronEnDownHandle;
  if(patType1CorrectedPFMetsElectronEnDown.size()){
    metunc.sumEt =(patType1CorrectedPFMetsElectronEnDown[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1CorrectedPFMetsElectronEnDown[0]));
    metunc.eLong=(patType1CorrectedPFMetsElectronEnDown[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1CorrectedPFMetsElectronEnDown[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }

  edm::Handle<edm::View<reco::MET> > patType1CorrectedPFMetElectronEnUpHandle;
  iEvent.getByLabel("patType1CorrectedPFMetElectronEnUp",patType1CorrectedPFMetElectronEnUpHandle);
  edm::View<reco::MET> patType1CorrectedPFMetsElectronEnUp = *patType1CorrectedPFMetElectronEnUpHandle;
  if(patType1CorrectedPFMetsElectronEnUp.size()){
    metunc.sumEt =(patType1CorrectedPFMetsElectronEnUp[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1CorrectedPFMetsElectronEnUp[0]));
    metunc.eLong=(patType1CorrectedPFMetsElectronEnUp[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1CorrectedPFMetsElectronEnUp[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }



  edm::Handle<edm::View<reco::MET> > patType1CorrectedPFMetMuonEnDownHandle;
  iEvent.getByLabel("patType1CorrectedPFMetMuonEnDown",patType1CorrectedPFMetMuonEnDownHandle);
  edm::View<reco::MET> patType1CorrectedPFMetsMuonEnDown = *patType1CorrectedPFMetMuonEnDownHandle;
  if(patType1CorrectedPFMetsMuonEnDown.size()){
    metunc.sumEt =(patType1CorrectedPFMetsMuonEnDown[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1CorrectedPFMetsMuonEnDown[0]));
    metunc.eLong=(patType1CorrectedPFMetsMuonEnDown[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1CorrectedPFMetsMuonEnDown[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }

  edm::Handle<edm::View<reco::MET> > patType1CorrectedPFMetMuonEnUpHandle;
  iEvent.getByLabel("patType1CorrectedPFMetMuonEnUp",patType1CorrectedPFMetMuonEnUpHandle);
  edm::View<reco::MET> patType1CorrectedPFMetsMuonEnUp = *patType1CorrectedPFMetMuonEnUpHandle;
  if(patType1CorrectedPFMetsMuonEnUp.size()){
    metunc.sumEt =(patType1CorrectedPFMetsMuonEnUp[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1CorrectedPFMetsMuonEnUp[0]));
    metunc.eLong=(patType1CorrectedPFMetsMuonEnUp[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1CorrectedPFMetsMuonEnUp[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }



  edm::Handle<edm::View<reco::MET> > patType1CorrectedPFMetTauEnDownHandle;
  iEvent.getByLabel("patType1CorrectedPFMetTauEnDown",patType1CorrectedPFMetTauEnDownHandle);
  edm::View<reco::MET> patType1CorrectedPFMetsTauEnDown = *patType1CorrectedPFMetTauEnDownHandle;
  if(patType1CorrectedPFMetsTauEnDown.size()){
    metunc.sumEt =(patType1CorrectedPFMetsTauEnDown[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1CorrectedPFMetsTauEnDown[0]));
    metunc.eLong=(patType1CorrectedPFMetsTauEnDown[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1CorrectedPFMetsTauEnDown[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }

  edm::Handle<edm::View<reco::MET> > patType1CorrectedPFMetTauEnUpHandle;
  iEvent.getByLabel("patType1CorrectedPFMetTauEnUp",patType1CorrectedPFMetTauEnUpHandle);
  edm::View<reco::MET> patType1CorrectedPFMetsTauEnUp = *patType1CorrectedPFMetTauEnUpHandle;
  if(patType1CorrectedPFMetsTauEnUp.size()){
    metunc.sumEt =(patType1CorrectedPFMetsTauEnUp[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1CorrectedPFMetsTauEnUp[0]));
    metunc.eLong=(patType1CorrectedPFMetsTauEnUp[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1CorrectedPFMetsTauEnUp[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }


  edm::Handle<edm::View<reco::MET> > patType1CorrectedPFMetJetEnDownHandle;
  iEvent.getByLabel("patType1CorrectedPFMetJetEnDown",patType1CorrectedPFMetJetEnDownHandle);
  edm::View<reco::MET> patType1CorrectedPFMetsJetEnDown = *patType1CorrectedPFMetJetEnDownHandle;
  if(patType1CorrectedPFMetsJetEnDown.size()){
    metunc.sumEt =(patType1CorrectedPFMetsJetEnDown[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1CorrectedPFMetsJetEnDown[0]));
    metunc.eLong=(patType1CorrectedPFMetsJetEnDown[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1CorrectedPFMetsJetEnDown[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }

  edm::Handle<edm::View<reco::MET> > patType1CorrectedPFMetJetEnUpHandle;
  iEvent.getByLabel("patType1CorrectedPFMetJetEnUp",patType1CorrectedPFMetJetEnUpHandle);
  edm::View<reco::MET> patType1CorrectedPFMetsJetEnUp = *patType1CorrectedPFMetJetEnUpHandle;
  if(patType1CorrectedPFMetsJetEnUp.size()){
    metunc.sumEt =(patType1CorrectedPFMetsJetEnUp[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1CorrectedPFMetsJetEnUp[0]));
    metunc.eLong=(patType1CorrectedPFMetsJetEnUp[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1CorrectedPFMetsJetEnUp[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }


  edm::Handle<edm::View<reco::MET> > patType1CorrectedPFMetJetResDownHandle;
  iEvent.getByLabel("patType1CorrectedPFMetJetResDown",patType1CorrectedPFMetJetResDownHandle);
  edm::View<reco::MET> patType1CorrectedPFMetsJetResDown = *patType1CorrectedPFMetJetResDownHandle;
  if(patType1CorrectedPFMetsJetResDown.size()){
    metunc.sumEt =(patType1CorrectedPFMetsJetResDown[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1CorrectedPFMetsJetResDown[0]));
    metunc.eLong=(patType1CorrectedPFMetsJetResDown[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1CorrectedPFMetsJetResDown[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }

  edm::Handle<edm::View<reco::MET> > patType1CorrectedPFMetJetResUpHandle;
  iEvent.getByLabel("patType1CorrectedPFMetJetResUp",patType1CorrectedPFMetJetResUpHandle);
  edm::View<reco::MET> patType1CorrectedPFMetsJetResUp = *patType1CorrectedPFMetJetResUpHandle;
  if(patType1CorrectedPFMetsJetResUp.size()){
    metunc.sumEt =(patType1CorrectedPFMetsJetResUp[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1CorrectedPFMetsJetResUp[0]));
    metunc.eLong=(patType1CorrectedPFMetsJetResUp[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1CorrectedPFMetsJetResUp[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }


  edm::Handle<edm::View<reco::MET> > patType1CorrectedPFMetUnclusteredEnDownHandle;
  iEvent.getByLabel("patType1CorrectedPFMetUnclusteredEnDown",patType1CorrectedPFMetUnclusteredEnDownHandle);
  edm::View<reco::MET> patType1CorrectedPFMetsUnclusteredEnDown = *patType1CorrectedPFMetUnclusteredEnDownHandle;
  if(patType1CorrectedPFMetsUnclusteredEnDown.size()){
    metunc.sumEt =(patType1CorrectedPFMetsUnclusteredEnDown[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1CorrectedPFMetsUnclusteredEnDown[0]));
    metunc.eLong=(patType1CorrectedPFMetsUnclusteredEnDown[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1CorrectedPFMetsUnclusteredEnDown[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }

  edm::Handle<edm::View<reco::MET> > patType1CorrectedPFMetUnclusteredEnUpHandle;
  iEvent.getByLabel("patType1CorrectedPFMetUnclusteredEnUp",patType1CorrectedPFMetUnclusteredEnUpHandle);
  edm::View<reco::MET> patType1CorrectedPFMetsUnclusteredEnUp = *patType1CorrectedPFMetUnclusteredEnUpHandle;
  if(patType1CorrectedPFMetsUnclusteredEnUp.size()){
    metunc.sumEt =(patType1CorrectedPFMetsUnclusteredEnUp[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1CorrectedPFMetsUnclusteredEnUp[0]));
    metunc.eLong=(patType1CorrectedPFMetsUnclusteredEnUp[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1CorrectedPFMetsUnclusteredEnUp[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }


  edm::Handle<edm::View<reco::MET> > patType1p2CorrectedPFMetElectronEnDownHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMetElectronEnDown",patType1p2CorrectedPFMetElectronEnDownHandle);
  edm::View<reco::MET> patType1p2CorrectedPFMetsElectronEnDown = *patType1p2CorrectedPFMetElectronEnDownHandle;
  if(patType1p2CorrectedPFMetsElectronEnDown.size()){
    metunc.sumEt =(patType1p2CorrectedPFMetsElectronEnDown[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1p2CorrectedPFMetsElectronEnDown[0]));
    metunc.eLong=(patType1p2CorrectedPFMetsElectronEnDown[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1p2CorrectedPFMetsElectronEnDown[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }

  edm::Handle<edm::View<reco::MET> > patType1p2CorrectedPFMetElectronEnUpHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMetElectronEnUp",patType1p2CorrectedPFMetElectronEnUpHandle);
  edm::View<reco::MET> patType1p2CorrectedPFMetsElectronEnUp = *patType1p2CorrectedPFMetElectronEnUpHandle;
  if(patType1p2CorrectedPFMetsElectronEnUp.size()){
    metunc.sumEt =(patType1p2CorrectedPFMetsElectronEnUp[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1p2CorrectedPFMetsElectronEnUp[0]));
    metunc.eLong=(patType1p2CorrectedPFMetsElectronEnUp[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1p2CorrectedPFMetsElectronEnUp[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }



  edm::Handle<edm::View<reco::MET> > patType1p2CorrectedPFMetMuonEnDownHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMetMuonEnDown",patType1p2CorrectedPFMetMuonEnDownHandle);
  edm::View<reco::MET> patType1p2CorrectedPFMetsMuonEnDown = *patType1p2CorrectedPFMetMuonEnDownHandle;
  if(patType1p2CorrectedPFMetsMuonEnDown.size()){
    metunc.sumEt =(patType1p2CorrectedPFMetsMuonEnDown[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1p2CorrectedPFMetsMuonEnDown[0]));
    metunc.eLong=(patType1p2CorrectedPFMetsMuonEnDown[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1p2CorrectedPFMetsMuonEnDown[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }

  edm::Handle<edm::View<reco::MET> > patType1p2CorrectedPFMetMuonEnUpHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMetMuonEnUp",patType1p2CorrectedPFMetMuonEnUpHandle);
  edm::View<reco::MET> patType1p2CorrectedPFMetsMuonEnUp = *patType1p2CorrectedPFMetMuonEnUpHandle;
  if(patType1p2CorrectedPFMetsMuonEnUp.size()){
    metunc.sumEt =(patType1p2CorrectedPFMetsMuonEnUp[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1p2CorrectedPFMetsMuonEnUp[0]));
    metunc.eLong=(patType1p2CorrectedPFMetsMuonEnUp[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1p2CorrectedPFMetsMuonEnUp[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }



  edm::Handle<edm::View<reco::MET> > patType1p2CorrectedPFMetTauEnDownHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMetTauEnDown",patType1p2CorrectedPFMetTauEnDownHandle);
  edm::View<reco::MET> patType1p2CorrectedPFMetsTauEnDown = *patType1p2CorrectedPFMetTauEnDownHandle;
  if(patType1p2CorrectedPFMetsTauEnDown.size()){
    metunc.sumEt =(patType1p2CorrectedPFMetsTauEnDown[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1p2CorrectedPFMetsTauEnDown[0]));
    metunc.eLong=(patType1p2CorrectedPFMetsTauEnDown[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1p2CorrectedPFMetsTauEnDown[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }

  edm::Handle<edm::View<reco::MET> > patType1p2CorrectedPFMetTauEnUpHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMetTauEnUp",patType1p2CorrectedPFMetTauEnUpHandle);
  edm::View<reco::MET> patType1p2CorrectedPFMetsTauEnUp = *patType1p2CorrectedPFMetTauEnUpHandle;
  if(patType1p2CorrectedPFMetsTauEnUp.size()){
    metunc.sumEt =(patType1p2CorrectedPFMetsTauEnUp[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1p2CorrectedPFMetsTauEnUp[0]));
    metunc.eLong=(patType1p2CorrectedPFMetsTauEnUp[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1p2CorrectedPFMetsTauEnUp[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }


  edm::Handle<edm::View<reco::MET> > patType1p2CorrectedPFMetJetEnDownHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMetJetEnDown",patType1p2CorrectedPFMetJetEnDownHandle);
  edm::View<reco::MET> patType1p2CorrectedPFMetsJetEnDown = *patType1p2CorrectedPFMetJetEnDownHandle;
  if(patType1p2CorrectedPFMetsJetEnDown.size()){
    metunc.sumEt =(patType1p2CorrectedPFMetsJetEnDown[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1p2CorrectedPFMetsJetEnDown[0]));
    metunc.eLong=(patType1p2CorrectedPFMetsJetEnDown[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1p2CorrectedPFMetsJetEnDown[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }

  edm::Handle<edm::View<reco::MET> > patType1p2CorrectedPFMetJetEnUpHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMetJetEnUp",patType1p2CorrectedPFMetJetEnUpHandle);
  edm::View<reco::MET> patType1p2CorrectedPFMetsJetEnUp = *patType1p2CorrectedPFMetJetEnUpHandle;
  if(patType1p2CorrectedPFMetsJetEnUp.size()){
    metunc.sumEt =(patType1p2CorrectedPFMetsJetEnUp[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1p2CorrectedPFMetsJetEnUp[0]));
    metunc.eLong=(patType1p2CorrectedPFMetsJetEnUp[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1p2CorrectedPFMetsJetEnUp[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }


  edm::Handle<edm::View<reco::MET> > patType1p2CorrectedPFMetJetResDownHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMetJetResDown",patType1p2CorrectedPFMetJetResDownHandle);
  edm::View<reco::MET> patType1p2CorrectedPFMetsJetResDown = *patType1p2CorrectedPFMetJetResDownHandle;
  if(patType1p2CorrectedPFMetsJetResDown.size()){
    metunc.sumEt =(patType1p2CorrectedPFMetsJetResDown[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1p2CorrectedPFMetsJetResDown[0]));
    metunc.eLong=(patType1p2CorrectedPFMetsJetResDown[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1p2CorrectedPFMetsJetResDown[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }

  edm::Handle<edm::View<reco::MET> > patType1p2CorrectedPFMetJetResUpHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMetJetResUp",patType1p2CorrectedPFMetJetResUpHandle);
  edm::View<reco::MET> patType1p2CorrectedPFMetsJetResUp = *patType1p2CorrectedPFMetJetResUpHandle;
  if(patType1p2CorrectedPFMetsJetResUp.size()){
    metunc.sumEt =(patType1p2CorrectedPFMetsJetResUp[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1p2CorrectedPFMetsJetResUp[0]));
    metunc.eLong=(patType1p2CorrectedPFMetsJetResUp[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1p2CorrectedPFMetsJetResUp[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }


  edm::Handle<edm::View<reco::MET> > patType1p2CorrectedPFMetUnclusteredEnDownHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMetUnclusteredEnDown",patType1p2CorrectedPFMetUnclusteredEnDownHandle);
  edm::View<reco::MET> patType1p2CorrectedPFMetsUnclusteredEnDown = *patType1p2CorrectedPFMetUnclusteredEnDownHandle;
  if(patType1p2CorrectedPFMetsUnclusteredEnDown.size()){
    metunc.sumEt =(patType1p2CorrectedPFMetsUnclusteredEnDown[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1p2CorrectedPFMetsUnclusteredEnDown[0]));
    metunc.eLong=(patType1p2CorrectedPFMetsUnclusteredEnDown[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1p2CorrectedPFMetsUnclusteredEnDown[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }

  edm::Handle<edm::View<reco::MET> > patType1p2CorrectedPFMetUnclusteredEnUpHandle;
  iEvent.getByLabel("patType1p2CorrectedPFMetUnclusteredEnUp",patType1p2CorrectedPFMetUnclusteredEnUpHandle);
  edm::View<reco::MET> patType1p2CorrectedPFMetsUnclusteredEnUp = *patType1p2CorrectedPFMetUnclusteredEnUpHandle;
  if(patType1p2CorrectedPFMetsUnclusteredEnUp.size()){
    metunc.sumEt =(patType1p2CorrectedPFMetsUnclusteredEnUp[0]).sumEt();
    metunc.metSig=metSignificance(& (patType1p2CorrectedPFMetsUnclusteredEnUp[0]));
    metunc.eLong=(patType1p2CorrectedPFMetsUnclusteredEnUp[0]).e_longitudinal();
    metunc.p4=GENPTOLOR((patType1p2CorrectedPFMetsUnclusteredEnUp[0]));
    hbbInfo->metUncInfo.push_back(metunc); 
  }



  //
  // met is calomet
  //

  edm::Handle<edm::View<pat::MET> > metTCHandle;
  iEvent.getByLabel("patMETsTC",metTCHandle);
  edm::View<pat::MET> metsTC = *metTCHandle;
  if(metsTC.size()){
    hbbInfo->tcmet.sumEt=(metsTC[0]).sumEt();
    hbbInfo->tcmet.metSig=metSignificance(&(metsTC[0]));
    hbbInfo->tcmet.eLong=(metsTC[0]).e_longitudinal();
    hbbInfo->tcmet.p4=GENPTOLOR((metsTC[0]));
    if (verbose_)     std::cout <<" METTC "<<     hbbInfo->tcmet.metSig <<" " <<     hbbInfo->tcmet.sumEt<<std::endl;
  }
  
  edm::Handle<edm::View<reco::MET> > pfMETNoPUHandle;
  iEvent.getByLabel("pfMETNoPU",pfMETNoPUHandle);
  edm::View<reco::MET> metspfMETNoPU = *pfMETNoPUHandle;
  if(metspfMETNoPU.size()){
    hbbInfo->metNoPU.sumEt=(metspfMETNoPU[0]).sumEt();
    hbbInfo->metNoPU.metSig=metSignificance(&(metspfMETNoPU[0]));                
    hbbInfo->metNoPU.eLong=(metspfMETNoPU[0]).e_longitudinal();
    hbbInfo->metNoPU.p4=GENPTOLOR((metspfMETNoPU[0]));
    if (verbose_)     std::cout <<" pfMETNoPU "<<     hbbInfo->metNoPU.metSig <<" " <<     hbbInfo->metNoPU.sumEt<<std::endl;
  }
 
  edm::Handle<edm::View<reco::MET> > mHTHandle;
  iEvent.getByLabel("patMETsHT",mHTHandle);
  edm::View<reco::MET> metsHT = *mHTHandle;
  if(metsHT.size()){
    hbbInfo->mht.sumEt=(metsHT[0]).sumEt();
    hbbInfo->mht.metSig=metSignificance(&(metsHT[0]));
    hbbInfo->mht.eLong=(metsHT[0]).e_longitudinal();
    hbbInfo->mht.p4=GENPTOLOR((metsHT[0]));
    if (verbose_)     std::cout <<" METHT "<<     hbbInfo->mht.metSig <<" " <<     hbbInfo->mht.sumEt<<std::endl;
  }
  
  edm::Handle<edm::View<reco::MET> > metHandle;
  iEvent.getByLabel(metLabel_,metHandle);
  edm::View<reco::MET> mets = *metHandle;
  
  if(mets.size()){
    hbbInfo->calomet.sumEt=(mets[0]).sumEt();
    hbbInfo->calomet.metSig=metSignificance(&(mets[0]));
    hbbInfo->calomet.eLong=(mets[0]).e_longitudinal();
    hbbInfo->calomet.p4=GENPTOLOR((mets[0]));
    if (verbose_)     std::cout <<" METCALO "<<     hbbInfo->calomet.metSig <<" " <<     hbbInfo->calomet.sumEt<<std::endl;
  }
  
  edm::Handle<edm::View<pat::MET> > metPFHandle;
  iEvent.getByLabel("patMETsPFlow",metPFHandle);
  edm::View<pat::MET> metsPF = *metPFHandle;
  
  if(metsPF.size()){
    hbbInfo->pfmet.sumEt=(metsPF[0]).sumEt();
    hbbInfo->pfmet.metSig=metSignificance(&(metsPF[0]));
    hbbInfo->pfmet.eLong=(metsPF[0]).e_longitudinal();
    hbbInfo->pfmet.p4=GENPTOLOR((metsPF[0]));
    if (verbose_)     std::cout <<" METPF "<<     hbbInfo->pfmet.metSig <<" " <<     hbbInfo->pfmet.sumEt<<std::endl;
  }
  
  
  if(verbose_){
    std::cout << "METs: calomet "<<mets.size()<<" tcmet"<<metsTC.size()<<" pfmet "<<metsPF.size()<<" MHT" <<metsHT.size()<<std::endl;  
  }

  if(verbose_)
    std::cout << " INPUT MUONS "<<muons.size()<<std::endl;

  for(edm::View<pat::Muon>::const_iterator mu = muons.begin(); mu!=muons.end(); ++mu){
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
  InputTag  reducedEBRecHitCollection(string("reducedEcalRecHitsEB"));
  InputTag  reducedEERecHitCollection(string("reducedEcalRecHitsEE"));
  EcalClusterLazyTools lazyTools(iEvent, iSetup, reducedEBRecHitCollection, reducedEERecHitCollection);
  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  const TransientTrackBuilder & transientTrackBuilder= *(builder.product());

  for(edm::View<pat::Electron>::const_iterator elec = electrons.begin(); elec!=electrons.end(); ++elec){
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
    for(size_t k=0;k<photonsForIso.size();k++) {
       if(deltaR(elec->eta(),elec->phi(),photonsForIso[k].eta(),photonsForIso[k].phi()) < 0.05 && abs(photonsForIso[k].pt()-elec->pt())/(photonsForIso[k].pt()+elec->pt()) < 0.05 ) {
          std::cout << "Double counting of supercluster!" << std::endl;
	  ef.pfPhoIsoDoubleCounted+=photonsForIso[k].pt(); 
     }
    }
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
  ef.fMVAVar_see             =  ele.sigmaIetaIeta();    //EleSigmaIEtaIEta
  std::vector<float> vCov = lazyTools.localCovariances(*(ele.superCluster()->seed())) ;
  if (!isnan(vCov[2])) ef.fMVAVar_spp = sqrt (vCov[2]);   //EleSigmaIPhiIPhi
  else ef.fMVAVar_spp = 0.;    
  // fMVAVar_sigmaIEtaIPhi = vCov[1];  //  save also this in your ntuple 

  ef.fMVAVar_etawidth        =  ele.superCluster()->etaWidth();
  ef.fMVAVar_phiwidth        =  ele.superCluster()->phiWidth();
  ef.fMVAVar_e1x5e5x5        =  (ele.e5x5()) !=0. ? 1.-(ele.e1x5()/ele.e5x5()) : -1. ;
  ef.fMVAVar_R9              =  lazyTools.e3x3(*(ele.superCluster()->seed())) / ele.superCluster()->rawEnergy();
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
    ef.fMVAVar_ip3d = -999.0; 
    // fMVAVar_ip3dSig = 0.0;
    if (ele.gsfTrack().isNonnull()) {
      const double gsfsign   = ( (-ele.gsfTrack()->dxy(vertex.position()))   >=0 ) ? 1. : -1.;
      
      const reco::TransientTrack &tt = transientTrackBuilder.build(ele.gsfTrack()); 
      const std::pair<bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt,vertex);
      if (ip3dpv.first) {
	double ip3d = gsfsign*ip3dpv.second.value();
	//double ip3derr = ip3dpv.second.error();  
	ef.fMVAVar_ip3d = ip3d; 
	// fMVAVar_ip3dSig = ip3d/ip3derr;
      }
    }
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
    ef.id95 =elec->electronID("eidVBTFCom95");
    ef.id95r=elec->electronID("eidVBTFRel95");
    ef.id85 =elec->electronID("eidVBTFCom85");
    ef.id85r=elec->electronID("eidVBTFRel85");
    ef.id80 =elec->electronID("eidVBTFCom80");
    ef.id80r=elec->electronID("eidVBTFRel80");
    ef.id70 =elec->electronID("eidVBTFCom70");
    ef.id70r=elec->electronID("eidVBTFRel70");
    ef.mvaOut=elec->electronID("mvaNonTrigV0");
    ef.mvaOutTrig=elec->electronID("mvaTrigV0");

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
  for(edm::View<pat::Tau>::const_iterator tau = taus.begin(); tau!=taus.end(); ++tau){
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


  iEvent.put(hbbInfo);
  iEvent.put(auxInfo);


  delete jecUnc;

}
  
void
HbbAnalyzerNew::fillMuBlock(edm::View<pat::Muon>::const_iterator mu, int muInfo[15])
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
HbbAnalyzerNew::beginJob(){
}


// ------------ method called once each job just after ending the event loop  ------------
void 
HbbAnalyzerNew::endJob() {
}

TVector2 HbbAnalyzerNew::getTvect( const pat::Jet* patJet ){

  TVector2 t_Vect(0,0);
  TVector2 null(0,0);
  TVector2 ci(0,0);
  TLorentzVector pi(0,0,0,0);
  TLorentzVector J(0,0,0,0);
  TVector2 r(0,0);
  double patJetpfcPt = 1e10;
  double r_mag = 1e10;
  unsigned int nOfconst = 0;


  if (patJet->isPFJet() == false) {
    return t_Vect;
  }
  

  //re-reconstruct the jet direction with the charged tracks
  std::vector<reco::PFCandidatePtr>
    patJetpfc = patJet->getPFConstituents();
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
      r.Set( pi.Rapidity() - J.Rapidity(), Geom::deltaPhi( patJetpfc.at(idx)->phi(), J.Phi() ) );
      r_mag = r.Mod();
      t_Vect += ( patJetpfcPt / J.Pt() ) * r_mag * r;
    }
  }

  
  return t_Vect;
  
}

TLorentzVector HbbAnalyzerNew::getChargedTracksMomentum(const pat::Jet* patJet ){
  //  return TLorentzVector();
  TLorentzVector pi(0,0,0,0);
  TLorentzVector v_j1(0,0,0,0);


  //  std::cout <<"fff ECCCCCCOOOOO "<<patJet->isPFJet()<<std::endl;

  if (patJet->isPFJet() == false ){
      v_j1 = GENPTOLORP(patJet);
      return v_j1;
  }
  std::vector<reco::PFCandidatePtr>
    j1pfc = patJet->getPFConstituents();
  for(size_t idx = 0; idx < j1pfc.size(); idx++){
    if( j1pfc.at(idx)->charge() != 0 ){
      pi.SetPtEtaPhiE( j1pfc.at(idx)->pt(), j1pfc.at(idx)->eta(), j1pfc.at(idx)->phi(), j1pfc.at(idx)->energy() );
      v_j1 += pi;
    }
  }
  return v_j1;
  //re-
}


//Btagging scale factors
void HbbAnalyzerNew::fillScaleFactors(VHbbEvent::SimpleJet& sj, BTagSFContainer iSF){


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

void HbbAnalyzerNew::setJecUnc(VHbbEvent::SimpleJet& sj,JetCorrectionUncertainty* jecunc){
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


void HbbAnalyzerNew ::fillSimpleJet (VHbbEvent::SimpleJet& sj, edm::View<pat::Jet>::const_iterator jet_iter){
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
    sj.SF_CSVL=1;
    sj.SF_CSVM=1;
    sj.SF_CSVT=1;
    sj.SF_CSVLerr=0;
    sj.SF_CSVMerr=0;
    sj.SF_CSVTerr=0;

    


    
    if (jet_iter->isPFJet() == true) {

      sj.chargedHadronEFraction = jet_iter-> chargedHadronEnergyFraction();
      sj.neutralHadronEFraction = jet_iter-> neutralHadronEnergyFraction ();
      sj.chargedEmEFraction = jet_iter-> chargedEmEnergyFraction ();
      sj.neutralEmEFraction = jet_iter-> neutralEmEnergyFraction ();
      sj.nConstituents = jet_iter->getPFConstituents().size();
      
    }
    sj.jetArea = jet_iter->jetArea();
    //
    // addtaginfo for csv
    //

    //    if (jet_iter->hasTagInfo("SimpleSecondaryVertex")) {

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

   //
    // add tVector
    //
    sj.tVector = getTvect(&(*jet_iter));

    sj.ptRaw = jet_iter->correctedJet(0).pt();

    sj.ptLeadTrack =-9999.;
    if (jet_iter->isPFJet() == true) {
       std::vector <reco::PFCandidatePtr> constituents = jet_iter->getPFConstituents ();
       for (unsigned ic = 0; ic < constituents.size (); ++ic) {
         if ( constituents[ic]->particleId() > 3 ) continue;
         reco::TrackRef trackRef = constituents[ic]->trackRef();
       if ( trackRef.isNonnull() ) { if(trackRef->pt() > sj.ptLeadTrack) sj.ptLeadTrack=trackRef->pt(); }
      }
     }


}

float HbbAnalyzerNew::metSignificance(const reco::MET * met)
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


//define this as a plug-in
DEFINE_FWK_MODULE(HbbAnalyzerNew);
