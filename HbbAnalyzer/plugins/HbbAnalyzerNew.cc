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
// $Id: HbbAnalyzerNew.cc,v 1.35 2011/09/14 14:15:04 tboccali Exp $
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
  jetLabel_(iConfig.getParameter<edm::InputTag>("jetTag")),
  subjetLabel_(iConfig.getParameter<edm::InputTag>("subjetTag")),
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
  iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl); 
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
  
  auxInfo->pvInfo.firstPVInPT2 = TVector3(RecVtxFirst.x(), RecVtxFirst.y(), RecVtxFirst.z());
  auxInfo->pvInfo.firstPVInProb = TVector3(RecVtx.x(), RecVtx.y(), RecVtx.z());

    
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(edm::InputTag("kt6PFJets", "rho"),rhoHandle);   
  auxInfo->puInfo.rho = *rhoHandle;
  
  edm::Handle<std::vector< PileupSummaryInfo> > puHandle;

  if (runOnMC_){
    iEvent.getByType(puHandle);
    if (puHandle.isValid()){
      
      std::vector< PileupSummaryInfo> pu = (*puHandle); 
      for (std::vector<PileupSummaryInfo>::const_iterator it= pu.begin(); it!=pu.end(); ++it){
	 int bx = (*it).getBunchCrossing();
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
	  ztemp.dauFourMomentum.push_back(GENPTOLOR(p));
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

  // standard jets

  edm::Handle<edm::View<pat::Jet> > simplejet1Handle;
  iEvent.getByLabel(simplejet1Label_,simplejet1Handle);
  edm::View<pat::Jet> simplejets1 = *simplejet1Handle;

  edm::Handle<edm::View<pat::Jet> > simplejet2Handle;
  iEvent.getByLabel(simplejet2Label_,simplejet2Handle);
  edm::View<pat::Jet> simplejets2 = *simplejet2Handle;

  edm::Handle<edm::View<pat::Jet> > simplejet3Handle;
  iEvent.getByLabel(simplejet3Label_,simplejet3Handle);
  edm::View<pat::Jet> simplejets3 = *simplejet3Handle;

  edm::Handle<edm::View<pat::Jet> > simplejet4Handle;
  iEvent.getByLabel(simplejet4Label_,simplejet4Handle);
  edm::View<pat::Jet> simplejets4 = *simplejet4Handle;


  edm::Handle<edm::View<pat::Electron> > electronHandle;
  iEvent.getByLabel(eleLabel_,electronHandle);
  edm::View<pat::Electron> electrons = *electronHandle;


  //   edm::Handle<edm::View<pat::Photon> > phoHandle;
  //   iEvent.getByLabel(phoLabel_,phoHandle);
  //   edm::View<pat::Photon> photons = *phoHandle;

  edm::Handle<edm::View<pat::Tau> > tauHandle;
  iEvent.getByLabel(tauLabel_,tauHandle);
  edm::View<pat::Tau> taus = *tauHandle;


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
  for(edm::View<pat::Jet>::const_iterator jet_iter = simplejets1.begin(); jet_iter!=simplejets1.end(); ++jet_iter){
    //     if(jet_iter->pt()>50)
    //       njetscounter++;
    VHbbEvent::SimpleJet sj;
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
      if( (jet_iter->genParton())
	  and (jet_iter->genParton()->mother()) )
	sj.bestMCmomid=jet_iter->genParton()->mother()->pdgId();
      TLorentzVector gJp4;
      if(gJ){
	gJp4.SetPtEtaPhiE(gJ->pt(),gJ->eta(),gJ->phi(),gJ->energy());
	sj. bestMCp4mom = gJp4;
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
      if( (jet_iter->genParton())
	  and (jet_iter->genParton()->mother()) )
	sj.bestMCmomid=jet_iter->genParton()->mother()->pdgId();
      TLorentzVector gJp4;
      if(gJ){
	gJp4.SetPtEtaPhiE(gJ->pt(),gJ->eta(),gJ->phi(),gJ->energy());
	sj. bestMCp4mom = gJp4;
	if(verbose_){
	  std::clog << "genJet matched Pt = " << gJp4.Pt() << std::endl;
	  std::clog << "genJet matched eta = " << gJp4.Eta() << std::endl;
	  std::clog << "genJet matched deltaR = " <<gJp4.DeltaR(sj.p4) << std::endl;
	  std::clog << "genJet matched mother id = " << sj.bestMCmomid << std::endl;
	}
      }
      
    } //isMC
    hbbInfo->simpleJets3.push_back(sj);
    
  }

#ifdef ENABLE_SIMPLEJETS4
  for(edm::View<pat::Jet>::const_iterator jet_iter = simplejets4.begin(); jet_iter!=simplejets4.end(); ++jet_iter){
    //     if(jet_iter->pt()>50)
    //       njetscounter++;
    VHbbEvent::SimpleJet sj;
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
      if( (jet_iter->genParton())
	  and (jet_iter->genParton()->mother()) )
	sj.bestMCmomid=jet_iter->genParton()->mother()->pdgId();
      TLorentzVector gJp4;
      if(gJ){
	gJp4.SetPtEtaPhiE(gJ->pt(),gJ->eta(),gJ->phi(),gJ->energy());
	sj. bestMCp4mom = gJp4;
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
    fillSimpleJet(sj,jet_iter);    
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
      if( (jet_iter->genParton())
	  and (jet_iter->genParton()->mother()) )
	sj.bestMCmomid=jet_iter->genParton()->mother()->pdgId();
      TLorentzVector gJp4;
      if(gJ){
	gJp4.SetPtEtaPhiE(gJ->pt(),gJ->eta(),gJ->phi(),gJ->energy());
	sj. bestMCp4mom = gJp4;
	if(verbose_){
	  std::clog << "genJet matched Pt = " << gJp4.Pt() << std::endl;
	  std::clog << "genJet matched eta = " << gJp4.Eta() << std::endl;
	  std::clog << "genJet matched deltaR = " << gJp4.DeltaR(sj.p4) << std::endl;
	  std::clog << "genJet matched mother id = " << sj.bestMCmomid << std::endl;
	}
      }

    }   //isMC
    
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


  //
  // met is calomet
  //

  edm::Handle<edm::View<pat::MET> > metTCHandle;
  iEvent.getByLabel("patMETsTC",metTCHandle);
  edm::View<pat::MET> metsTC = *metTCHandle;
  if(metsTC.size()){
    hbbInfo->tcmet.sumEt=(metsTC[0]).sumEt();
    hbbInfo->tcmet.metSig=(metsTC[0]).significance();
    hbbInfo->tcmet.eLong=(metsTC[0]).e_longitudinal();
    hbbInfo->tcmet.p4=GENPTOLOR((metsTC[0]));
    if (verbose_)     std::cout <<" METTC "<<     hbbInfo->tcmet.metSig <<" " <<     hbbInfo->tcmet.sumEt<<std::endl;
  }
  
 
  edm::Handle<edm::View<reco::MET> > mHTHandle;
  iEvent.getByLabel("patMETsHT",mHTHandle);
  edm::View<reco::MET> metsHT = *mHTHandle;
  if(metsHT.size()){
    hbbInfo->mht.sumEt=(metsHT[0]).sumEt();
    hbbInfo->mht.metSig=(metsHT[0]).significance();
    hbbInfo->mht.eLong=(metsHT[0]).e_longitudinal();
    hbbInfo->mht.p4=GENPTOLOR((metsHT[0]));
    if (verbose_)     std::cout <<" METHT "<<     hbbInfo->mht.metSig <<" " <<     hbbInfo->mht.sumEt<<std::endl;
  }
  
  edm::Handle<edm::View<pat::MET> > metHandle;
  iEvent.getByLabel(metLabel_,metHandle);
  edm::View<pat::MET> mets = *metHandle;
  
  if(mets.size()){
    hbbInfo->calomet.sumEt=(mets[0]).sumEt();
    hbbInfo->calomet.metSig=(mets[0]).significance();
    hbbInfo->calomet.eLong=(mets[0]).e_longitudinal();
    hbbInfo->calomet.p4=GENPTOLOR((mets[0]));
    if (verbose_)     std::cout <<" METCALO "<<     hbbInfo->calomet.metSig <<" " <<     hbbInfo->calomet.sumEt<<std::endl;
  }
  
  edm::Handle<edm::View<pat::MET> > metPFHandle;
  iEvent.getByLabel("patMETsPF",metPFHandle);
  edm::View<pat::MET> metsPF = *metPFHandle;
  
  if(metsPF.size()){
    hbbInfo->pfmet.sumEt=(metsPF[0]).sumEt();
    hbbInfo->pfmet.metSig=(metsPF[0]).significance();
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
    mf.pfPhoIso=mu->photonIso();
    mf.pfNeuIso=mu->neutralHadronIso(); 
    Geom::Phi<double> deltaphi(mu->phi()-atan2(mf.p4.Px(), mf.p4.Py()));
    double acop = deltaphi.value();
    mf.acop=acop;

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
    } 
    if(mu->isGlobalMuon()){
      TrackRef gTrack = mu->globalTrack();
      const reco::HitPattern& q = gTrack->hitPattern();
      mf.globChi2=gTrack.get()->normalizedChi2();
      mf.globNHits=q.numberOfValidMuonHits();
      mf. validMuStations = q. muonStationsWithValidHits();
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
    ef.pfPhoIso=elec->photonIso();
    ef.pfNeuIso=elec->neutralHadronIso();

    Geom::Phi<double> deltaphi(elec->superCluster()->phi()-atan2(hbbInfo->calomet.p4.Py(),hbbInfo->calomet.p4.Px()));
    ef.acop = deltaphi.value();
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
    Geom::Phi<double> deltaphi(tau->phi()-atan2(hbbInfo->calomet.p4.Py(),hbbInfo->calomet.p4.Px()));
    double acop = deltaphi.value();
    tf.acop=acop;
    tf.idbyIso=tau->tauID("byIsolation");
    tf.idbyTrackIso=tau->tauID("trackIsolation");
    tf.idbyTaNCfrOnePercent=tau->tauID("byTaNCfrOnePercent");
    tf.idbyTaNCfrHalfPercent=tau->tauID("byTaNCfrHalfPercent");
    tf.idbyTaNCfrQuarterPercent=tau->tauID("byTaNCfrQuarterPercent");
    tf.idbyTaNCfrTenthPercent=tau->tauID("byTaNCfrTenthPercent");  
    tf.idbyTaNC=tau->tauID("byTaNC");
    hbbInfo->tauInfo.push_back(tf);
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
       " Muons = "<<hbbInfo->muInfo.size()<<std::endl<<
       " Electrons = "<<hbbInfo->eleInfo.size()<<std::endl<<
       " Taus = "<<hbbInfo->tauInfo.size()<<std::endl<<
       " Electrons = "<<hbbInfo->eleInfo.size()<<std::endl<<
       "--------------------- "<<std::endl;
  }


  iEvent.put(hbbInfo);
  iEvent.put(auxInfo);


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

  return;
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

    
    if (jet_iter->isPFJet() == true) {

      sj.chargedHadronEFraction = jet_iter-> chargedHadronEnergyFraction();
      sj.neutralHadronEFraction = jet_iter-> neutralHadronEnergyFraction ();
      sj.chargedEmEFraction = jet_iter-> chargedEmEnergyFraction ();
      sj.neutralEmEFraction = jet_iter-> neutralEmEnergyFraction ();
      sj. nConstituents = jet_iter->getPFConstituents().size();
      
    }
    //
    // addtaginfo for csv
    //

    //    if (jet_iter->hasTagInfo("SimpleSecondaryVertex")) {

    const reco::SecondaryVertexTagInfo * tf = jet_iter->tagInfoSecondaryVertex();
   if (tf){
     if (tf->nVertices() >0){
	sj.vtxMass = tf->secondaryVertex(0).p4().mass();
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
   //
    // add tVector
    //
    sj.tVector = getTvect(&(*jet_iter));
}


//define this as a plug-in
DEFINE_FWK_MODULE(HbbAnalyzerNew);
