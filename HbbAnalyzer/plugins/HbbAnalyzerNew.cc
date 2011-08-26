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
// $Id: HbbAnalyzerNew.cc,v 1.25 2011/08/23 14:43:03 bortigno Exp $
//
//

#include "VHbbAnalysis/HbbAnalyzer/interface/HbbAnalyzerNew.h"

#define GENPTOLOR(a) TLorentzVector((a).px(), (a).py(), (a).pz(), (a).energy())
#define GENPTOLORP(a) TLorentzVector((a)->px(), (a)->py(), (a)->pz(), (a)->energy())

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
  dimuLabel_(iConfig.getParameter<edm::InputTag>("dimuTag")),
  dielecLabel_(iConfig.getParameter<edm::InputTag>("dielecTag")),
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
  
  
  std::auto_ptr<VHbbEvent> hbbInfo( new VHbbEvent() );  
  std::auto_ptr<VHbbEventAuxInfo> auxInfo( new VHbbEventAuxInfo() );
  
  //
  // ??
   
  // trigger
  edm::Handle<edm::TriggerResults>  hltresults;
  //iEvent.getByLabel("TriggerResults", hltresults);
   
  //edm::InputTag tag("TriggerResults::HLT");
  //  edm::InputTag tag("TriggerResults::HLT0");
  iEvent.getByLabel(hltResults_, hltresults);
   
  const edm::TriggerNames & triggerNames_ = iEvent.triggerNames(*hltresults);
   
  int ntrigs = hltresults->size();
  if (ntrigs==0){std::cerr << "%HLTInfo -- No trigger name given in TriggerResults of the input " << std::endl;}

  for (int itrig = 0; itrig != ntrigs; ++itrig){

    TString trigName=triggerNames_.triggerName(itrig);
    bool accept = hltresults->accept(itrig);

    if (accept){(auxInfo->triggerInfo.flag)[itrig] = 1;}
    else { (auxInfo->triggerInfo.flag)[itrig] = 0;}

    //    std::cout << "%HLTInfo --  Number of HLT Triggers: " << ntrigs << std::endl;
    //    std::cout << "%HLTInfo --  HLTTrigger(" << itrig << "): " << trigName << " = " << accept << std::endl;
  }

//   //
//   // big bloat
//   //


  int goodDoubleMu3=0,goodDoubleMu3_2=0,  
    goodMu9=0, goodIsoMu9=0, goodMu11=0, goodIsoMu13_3=0, goodMu15=0, goodMu15_1=0; 
  int goodDoubleElec10=0,goodDoubleElec15_1=0,goodDoubleElec17_1=0;
  int goodMet100=0;
  int goodSingleEle1=0,goodSingleEle2=0,goodSingleEle3=0,goodSingleEle4=0;
  int goodBtagMu1=0,goodBtagMu2=0,goodBtagMu0=0,goodBtagMu11=0;  
  int goodBtagMuJet1=0, goodBtagMuJet2=0, goodBtagMuJet3=0, goodBtagMuJet4=0; 
  int goodIsoMu15=0,goodIsoMu17v5=0,goodIsoMu17v6=0;

  for (int itrig = 0; itrig != ntrigs; ++itrig){
    TString trigName=triggerNames_.triggerName(itrig);
    if(strcmp(trigName,"HLT_Mu9")==0) goodMu9++;
    if(strcmp(trigName,"HLT_IsoMu9")==0) goodIsoMu9++;
    if(strcmp(trigName,"HLT_IsoMu13_v3")==0) goodIsoMu13_3++;
    if(strcmp(trigName,"HLT_Mu11")==0) goodMu11++;
    if(strcmp(trigName,"HLT_DoubleMu3")==0) goodDoubleMu3++;
    if(strcmp(trigName,"HLT_DoubleMu3_v2")==0) goodDoubleMu3_2++;
    if(strcmp(trigName,"HLT_Mu15")==0) goodMu15++;
    if(strcmp(trigName,"HLT_Mu15_v1")==0) goodMu15_1++;

    if(strcmp(trigName,"HLT_DoubleEle10_SW_L1R")==0) goodDoubleElec10++;
    if(strcmp(trigName,"HLT_DoubleEle15_SW_L1R_v1")==0) goodDoubleElec15_1++;
    if(strcmp(trigName,"HLT_DoubleEle17_SW_L1R_v1")==0) goodDoubleElec17_1++;
    if(strcmp(trigName,"HLT_MET100_v1")==0) goodMet100++;
    if(strcmp(trigName,"HLT_Ele15_SW_L1R")==0) goodSingleEle1++;
    if(strcmp(trigName,"HLT_Ele17_SW_TightEleId_L1R")==0) goodSingleEle2++;
    if(strcmp(trigName,"HLT_Ele17_SW_TighterEleIdIsol_L1R_v2")==0) goodSingleEle3++;
    if(strcmp(trigName,"HLT_Ele17_SW_TighterEleIdIsol_L1R_v3")==0) goodSingleEle4++;
    if(strcmp(trigName,"HLT_BTagMu_DiJet20U_v1")==0) goodBtagMu1++;
    if(strcmp(trigName,"HLT_BTagMu_DiJet30U_Mu5_v3")==0) goodBtagMu2++;
    if(strcmp(trigName,"HLT_BTagMu_Jet20U")==0) goodBtagMu0++;
    if(strcmp(trigName,"HLT_BTagMu_DiJet20U_Mu5_v1")==0) goodBtagMu11++;
    if(strcmp(trigName,"HLT_Mu17_CentralJet30_BTagIP_v2")==0) goodBtagMuJet1++;
    if(strcmp(trigName,"HLT_Mu17_CentralJet30_v2")==0) goodBtagMuJet2++;
    if(strcmp(trigName,"HLT_HT200_Mu5_PFMHT35_v2")==0) goodBtagMuJet3++;
    if(strcmp(trigName,"HLT_Mu12_CentralJet30_BTagIP_v2")==0) goodBtagMuJet4++;  

    if(strcmp(trigName,"HLT_IsoMu15_v5")==0) goodIsoMu15++; 
    if(strcmp(trigName,"HLT_IsoMu17_v5")==0) goodIsoMu17v5++;  
    if(strcmp(trigName,"HLT_IsoMu17_v6")==0) goodIsoMu17v6++;  
  }
  int itrig1=-99;
  if(goodMu9!=0) itrig1 = triggerNames_.triggerIndex("HLT_Mu9");
  int itrig2=-99;
  if(goodIsoMu9!=0) itrig2 = triggerNames_.triggerIndex("HLT_IsoMu9");
  int itrig3=-99;
  if(goodIsoMu13_3!=0) itrig3 = triggerNames_.triggerIndex("HLT_IsoMu13_v3"); 
  int itrig4=-99;
  if(goodMu11!=0) itrig4 = triggerNames_.triggerIndex("HLT_Mu11"); 
  int itrig5=-99;  
  if(goodDoubleMu3!=0) itrig5 = triggerNames_.triggerIndex("HLT_DoubleMu3"); 
  int itrig6=-99;
  if(goodDoubleMu3_2!=0) itrig6 = triggerNames_.triggerIndex("HLT_DoubleMu3_v2");
  int itrig7=-99;
  if(goodMu15!=0) itrig7 = triggerNames_.triggerIndex("HLT_Mu15");
  int itrig8=-99;
  if(goodMu15_1!=0) itrig8 = triggerNames_.triggerIndex("HLT_Mu15_v1"); 

  int itrig9=-99;
  if(goodDoubleElec10!=0) itrig9 = triggerNames_.triggerIndex("HLT_DoubleEle10_SW_L1R"); 
  int itrig10=-99;
  if(goodDoubleElec15_1!=0) itrig10 = triggerNames_.triggerIndex("HLT_DoubleEle15_SW_L1R_v1"); 
  int itrig11=-99;
  if(goodDoubleElec17_1!=0) itrig11 = triggerNames_.triggerIndex("HLT_DoubleEle17_SW_L1R_v1"); 
  int itrig12=-99;
  if(goodMet100!=0) itrig12 = triggerNames_.triggerIndex("HLT_MET100_v1"); 

  int itrig13=-99;
  if(goodSingleEle1!=0) itrig13 = triggerNames_.triggerIndex("HLT_Ele15_SW_L1R");
  int itrig14=-99;
  if(goodSingleEle2!=0) itrig14 = triggerNames_.triggerIndex("HLT_Ele17_SW_TightEleId_L1R");
  int itrig15=-99;
  if(goodSingleEle3!=0) itrig15 = triggerNames_.triggerIndex("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2");
  int itrig16=-99;
  if(goodSingleEle4!=0) itrig16 = triggerNames_.triggerIndex("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3");

  int itrig17=-99; 
  if(goodBtagMu1!=0) itrig17 = triggerNames_.triggerIndex("HLT_BTagMu_DiJet20U_v1");      
  int itrig18=-99;
  if(goodBtagMu2!=0) itrig18 = triggerNames_.triggerIndex("HLT_BTagMu_DiJet30U_Mu5_v3");    
  int itrig19=-99;
  if(goodBtagMu0!=0) itrig19 = triggerNames_.triggerIndex("HLT_BTagMu_Jet20U");
  int itrig20=-99;
  if(goodBtagMu11!=0) itrig20 = triggerNames_.triggerIndex("HLT_BTagMu_DiJet20U_Mu5_v1");
  int itrig21=-99;
  if(goodBtagMuJet1!=0) itrig21 = triggerNames_.triggerIndex("HLT_Mu17_CentralJet30_BTagIP_v2");
  int itrig22=-99;
  if(goodBtagMuJet2!=0) itrig22 = triggerNames_.triggerIndex("HLT_Mu17_CentralJet30_v2"); 
  int itrig23=-99;
  if(goodBtagMuJet3!=0) itrig23 = triggerNames_.triggerIndex("HLT_HT200_Mu5_PFMHT35_v2");
  int itrig231=-99;
  if(goodBtagMuJet4!=0) itrig231 = triggerNames_.triggerIndex("HLT_Mu12_CentralJet30_BTagIP_v2");  

  int itrig24=-99;
  if(goodIsoMu15!=0) itrig24 = triggerNames_.triggerIndex("HLT_IsoMu15_v5"); 
  int itrig25=-99;
  if(goodIsoMu17v5!=0) itrig25 = triggerNames_.triggerIndex("HLT_IsoMu17_v5");
  int itrig26=-99;
  if(goodIsoMu17v6!=0) itrig26 = triggerNames_.triggerIndex("HLT_IsoMu17_v6");
   
  if (itrig1!=-99 && hltresults->accept(itrig1))  auxInfo->triggerInfo.triggerMu9=1; else auxInfo->triggerInfo.triggerMu9=0;
  if (itrig2!=-99 && hltresults->accept(itrig2))  auxInfo->triggerInfo.triggerIsoMu9=1; else auxInfo->triggerInfo.triggerIsoMu9=0;
  if (itrig3!=-99 && hltresults->accept(itrig3))  auxInfo->triggerInfo.triggerIsoMu13_3=1; else auxInfo->triggerInfo.triggerIsoMu13_3=0;
  if (itrig4!=-99 && hltresults->accept(itrig4))  auxInfo->triggerInfo.triggerMu11=1; else auxInfo->triggerInfo.triggerMu11=0;
  if (itrig5!=-99 && hltresults->accept(itrig5))  auxInfo->triggerInfo.triggerDoubleMu3=1; else auxInfo->triggerInfo.triggerDoubleMu3=0; 
  if (itrig6!=-99 && hltresults->accept(itrig6))  auxInfo->triggerInfo.triggerDoubleMu3_2=1; else auxInfo->triggerInfo.triggerDoubleMu3_2=0;
  if (itrig7!=-99 && hltresults->accept(itrig7))  auxInfo->triggerInfo.triggerMu15=1; else auxInfo->triggerInfo.triggerMu15=0;
  if (itrig8!=-99 && hltresults->accept(itrig8))  auxInfo->triggerInfo.triggerMu15_1=1; else auxInfo->triggerInfo.triggerMu15_1=0;  
   
  if (itrig9!=-99 && hltresults->accept(itrig9))  auxInfo->triggerInfo.triggerDoubleElec10=1; else auxInfo->triggerInfo.triggerDoubleElec10=0;  
  if (itrig10!=-99 && hltresults->accept(itrig10))  auxInfo->triggerInfo.triggerDoubleElec15_1=1; else auxInfo->triggerInfo.triggerDoubleElec15_1=0;  
  if (itrig11!=-99 && hltresults->accept(itrig11))  auxInfo->triggerInfo.triggerDoubleElec17_1=1; else auxInfo->triggerInfo.triggerDoubleElec17_1=0;  
  if (itrig12!=-99 && hltresults->accept(itrig12))  auxInfo->triggerInfo.triggerMet100_1=1; else auxInfo->triggerInfo.triggerMet100_1=0;
   
  if (itrig13!=-99 && hltresults->accept(itrig13))  auxInfo->triggerInfo.triggerSingleEle1=1; else auxInfo->triggerInfo.triggerSingleEle1=0;
  if (itrig14!=-99 && hltresults->accept(itrig14))  auxInfo->triggerInfo.triggerSingleEle2=1; else auxInfo->triggerInfo.triggerSingleEle2=0;
  if (itrig15!=-99 && hltresults->accept(itrig15))  auxInfo->triggerInfo.triggerSingleEle3=1; else auxInfo->triggerInfo.triggerSingleEle3=0;
  if (itrig16!=-99 && hltresults->accept(itrig16))  auxInfo->triggerInfo.triggerSingleEle4=1; else auxInfo->triggerInfo.triggerSingleEle4=0;
   
  if (itrig17!=-99 && hltresults->accept(itrig17))  auxInfo->triggerInfo.triggerBtagMu1=1; else auxInfo->triggerInfo.triggerBtagMu1=0;   
  if (itrig18!=-99 && hltresults->accept(itrig18))  auxInfo->triggerInfo.triggerBtagMu2=1; else auxInfo->triggerInfo.triggerBtagMu2=0;
  if (itrig19!=-99 && hltresults->accept(itrig19))  auxInfo->triggerInfo.triggerBtagMu0=1; else auxInfo->triggerInfo.triggerBtagMu0=0;
  if (itrig20!=-99 && hltresults->accept(itrig20))  auxInfo->triggerInfo.triggerBtagMu11=1; else auxInfo->triggerInfo.triggerBtagMu11=0;
  if (itrig21!=-99 && hltresults->accept(itrig21))  auxInfo->triggerInfo.triggerBtagMuJet1=1; else auxInfo->triggerInfo.triggerBtagMuJet1=0;
  if (itrig22!=-99 && hltresults->accept(itrig22))  auxInfo->triggerInfo.triggerBtagMuJet2=1; else auxInfo->triggerInfo.triggerBtagMuJet2=0;
  if (itrig23!=-99 && hltresults->accept(itrig23))  auxInfo->triggerInfo.triggerBtagMuJet3=1; else auxInfo->triggerInfo.triggerBtagMuJet3=0;
  if (itrig231!=-99 && hltresults->accept(itrig231))  auxInfo->triggerInfo.triggerBtagMuJet4=1; else auxInfo->triggerInfo.triggerBtagMuJet4=0;
   
  if (itrig24!=-99 && hltresults->accept(itrig24))  auxInfo->triggerInfo.triggerIsoMu15=1; else auxInfo->triggerInfo.triggerIsoMu15=0;
  if (itrig25!=-99 && hltresults->accept(itrig25))  auxInfo->triggerInfo.triggerIsoMu17v5=1; else auxInfo->triggerInfo.triggerIsoMu17v5=0;
  if (itrig26!=-99 && hltresults->accept(itrig26))  auxInfo->triggerInfo.triggerIsoMu17v6=1; else auxInfo->triggerInfo.triggerIsoMu17v6=0;
   
  if (itrig1==-99)  auxInfo->triggerInfo.triggerMu9=-99; 
  if (itrig2==-99)  auxInfo->triggerInfo.triggerIsoMu9=-99; 
  if (itrig3==-99)  auxInfo->triggerInfo.triggerIsoMu13_3=-99; 
  if (itrig4==-99)  auxInfo->triggerInfo.triggerMu11=-99; 
  if (itrig5==-99)  auxInfo->triggerInfo.triggerDoubleMu3=-99;
  if (itrig6==-99)  auxInfo->triggerInfo.triggerDoubleMu3_2=-99; 
  if (itrig7==-99)  auxInfo->triggerInfo.triggerMu15=-99;
  if (itrig8==-99)  auxInfo->triggerInfo.triggerMu15_1=-99; 

  if (itrig9==-99)  auxInfo->triggerInfo.triggerDoubleElec10=-99; 
  if (itrig10==-99)  auxInfo->triggerInfo.triggerDoubleElec15_1=-99; 
  if (itrig11==-99)  auxInfo->triggerInfo.triggerDoubleElec17_1=-99; 
  if (itrig12==-99) auxInfo->triggerInfo.triggerMet100_1=-99;

  if (itrig13==-99) auxInfo->triggerInfo.triggerSingleEle1=-99;
  if (itrig14==-99) auxInfo->triggerInfo.triggerSingleEle2=-99;
  if (itrig15==-99) auxInfo->triggerInfo.triggerSingleEle3=-99;
  if (itrig16==-99) auxInfo->triggerInfo.triggerSingleEle4=-99;

  if(itrig17==-99)  auxInfo->triggerInfo.triggerBtagMu1=-99;   
  if(itrig18==-99)  auxInfo->triggerInfo.triggerBtagMu2=-99;
  if(itrig19==-99)  auxInfo->triggerInfo.triggerBtagMu0=-99;
  if(itrig20==-99)  auxInfo->triggerInfo.triggerBtagMu11=-99;
  if(itrig21==-99)  auxInfo->triggerInfo.triggerBtagMuJet1=-99;
  if(itrig22==-99)  auxInfo->triggerInfo.triggerBtagMuJet2=-99;
  if(itrig23==-99)  auxInfo->triggerInfo.triggerBtagMuJet3=-99;
  if(itrig231==-99)  auxInfo->triggerInfo.triggerBtagMuJet4=-99;

  if(itrig24==-99)  auxInfo->triggerInfo.triggerIsoMu15=-99;
  if(itrig25==-99)  auxInfo->triggerInfo.triggerIsoMu17v5=-99;
  if(itrig26==-99)  auxInfo->triggerInfo.triggerIsoMu17v6=-99;

  //  MinDRMu=-99.,MCBestMuId=-99,MCBestMuMomId=-99,MCBestMuGMomId=-99;
  //  for(int i=0;i<50;i++) {DeltaRMu[i]=-99;}



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
  iEvent.getByLabel(edm::InputTag("kt6PFJets", "rho"),rhoHandle);   auxInfo->puInfo.rho = *rhoHandle;
  
  //// real start
  
  
  Handle<GenParticleCollection> genParticles;
  
  bool printJet=0;
  
  
  if(runOnMC_){
    
    int Hin=-99,Win=-99,Zin=-99,bin=-99,bbarin=-99;
    iEvent.getByLabel("genParticles", genParticles);
    
    for(size_t i = 0; i < genParticles->size(); ++ i) {
      int Hin=-99,Win=-99,Zin=-99,bin=-99,bbarin=-99;
      
      const GenParticle & p = (*genParticles)[i];
      int id = p.pdgId();
      int st = p.status();  
      
      if(id==25){
	
	/*          int wh=0;
		    int nMoth = p.numberOfMothers();
		    for(size_t jj = 0; jj < nMoth; ++ jj) {
		    const Candidate * mom = p.mother( jj );
		    int nmomdau= mom->numberOfDaughters();
		    for(size_t j = 0; j < nmomdau; ++ j) {
		    const Candidate * Wmomdau = mom->daughter( j );
		    if(abs(Wmomdau->pdgId())==24) wh++;
		    }
		    }
		    
		    if(wh==0) continue;
	*/
	VHbbEventAuxInfo::ParticleMCInfo htemp;
	Hin=i; 
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
	  
	Win=i;
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

	Zin=i;
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
	bin=i; 
	VHbbEventAuxInfo::ParticleMCInfo btemp;
	btemp.status=st;
	btemp.charge=p.charge();
	if(p.mother(0)!=0) btemp.momid=p.mother(0)->pdgId();
	if(p.mother(0)!=0 && p.mother(0)->mother(0)!=0) btemp.gmomid=p.mother(0)->mother(0)->pdgId(); 
	auxInfo->mcB.push_back(btemp);
      }
      
      if(id==-5){
	bbarin=i; 
	VHbbEventAuxInfo::ParticleMCInfo bbtemp;
	
	bbtemp.status=st;
	bbtemp.charge=p.charge();
	if(p.mother(0)!=0) bbtemp.momid=p.mother(0)->pdgId();
	if(p.mother(0)!=0 && p.mother(0)->mother(0)!=0) bbtemp.gmomid=p.mother(0)->mother(0)->pdgId(); 
	auxInfo->mcBbar.push_back(bbtemp);
     }
      
      if(abs(id)==4){
	VHbbEventAuxInfo::ParticleMCInfo ctemp;
	ctemp.status=st;
	ctemp.charge=p.charge();
	if(p.mother(0)!=0) ctemp.momid=p.mother(0)->pdgId();
	if(p.mother(0)!=0 && p.mother(0)->mother(0)!=0) ctemp.gmomid=p.mother(0)->mother(0)->pdgId();
	auxInfo->mcC.push_back(ctemp);	
      }
    
      if(bin!=-99 && bbarin!=-99){
      const Candidate & bGen = (*genParticles)[bin];
      const Candidate & bbarGen = (*genParticles)[bbarin]; 
      ((auxInfo->mcB).back()).p4=GENPTOLOR(bGen);
     ((auxInfo->mcBbar).back()).p4=GENPTOLOR(bbarGen);
      
      int nHDaubdau = bGen.numberOfDaughters();
      for(int j = 0; j < nHDaubdau; ++ j) {
	const Candidate * Bdau = bGen.daughter( j );
	((auxInfo->mcB).back()).dauid.push_back(Bdau->pdgId());
	((auxInfo->mcB).back()).dauFourMomentum.push_back(GENPTOLORP(Bdau));
      }
      int nHDaubbardau = bbarGen.numberOfDaughters();
      for(int j = 0; j < nHDaubbardau; ++ j) {
	const Candidate * Bbardau = bbarGen.daughter( j );
	((auxInfo->mcBbar).back()).dauid.push_back(Bbardau->pdgId());
	((auxInfo->mcBbar).back()).dauFourMomentum.push_back(GENPTOLORP(Bbardau));
      }
      
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

  edm::Handle<edm::View<pat::MET> > metHandle;
  iEvent.getByLabel(metLabel_,metHandle);
  edm::View<pat::MET> mets = *metHandle;

  //   edm::Handle<edm::View<pat::Photon> > phoHandle;
  //   iEvent.getByLabel(phoLabel_,phoHandle);
  //   edm::View<pat::Photon> photons = *phoHandle;

  edm::Handle<edm::View<pat::Tau> > tauHandle;
  iEvent.getByLabel(tauLabel_,tauHandle);
  edm::View<pat::Tau> taus = *tauHandle;

  edm::Handle<CandidateView> dimuons;
  iEvent.getByLabel(dimuLabel_,dimuons);

  edm::Handle<CandidateView> dielectrons;
  iEvent.getByLabel(dielecLabel_,dielectrons);

  //BTAGGING SCALE FACTOR FROM DATABASE
  //Combined Secondary Vertex Loose
  edm::ESHandle<BtagPerformance> bTagSF_CSVL_;
  iSetup.get<BTagPerformanceRecord>().get("BTAGCSVL",bTagSF_CSVL_);
  const BtagPerformance & bTagSF_CSVL = *(bTagSF_CSVL_.product());
  //Combined Secondary Vertex Medium
  edm::ESHandle<BtagPerformance> bTagSF_CSVM_;
  iSetup.get<BTagPerformanceRecord>().get("BTAGCSVM",bTagSF_CSVM_);
  const BtagPerformance & bTagSF_CSVM = *(bTagSF_CSVM_.product());
  //Combined Secondary Vertex Tight
  edm::ESHandle<BtagPerformance> bTagSF_CSVT_;
  iSetup.get<BTagPerformanceRecord>().get("BTAGCSVT",bTagSF_CSVT_);
  const BtagPerformance & bTagSF_CSVT = *(bTagSF_CSVT_.product());

  edm::ESHandle<BtagPerformance> mistagSF_CSVL_;
  iSetup.get<BTagPerformanceRecord>().get("MISTAGCSVL",mistagSF_CSVL_);
  const BtagPerformance & mistagSF_CSVL = *(mistagSF_CSVL_.product());
  //Combined Secondary Vertex Medium
  edm::ESHandle<BtagPerformance> mistagSF_CSVM_;
  iSetup.get<BTagPerformanceRecord>().get("MISTAGCSVM",mistagSF_CSVM_);
  const BtagPerformance & mistagSF_CSVM = *(mistagSF_CSVM_.product());
  //Combined Secondary Vertex Tight
  edm::ESHandle<BtagPerformance> mistagSF_CSVT_;
  iSetup.get<BTagPerformanceRecord>().get("MISTAGCSVT",mistagSF_CSVT_);
  const BtagPerformance & mistagSF_CSVT = *(mistagSF_CSVT_.product());

  iBTV.BTAGSF_CSVL = (bTagSF_CSVL_.product());
  iBTV.BTAGSF_CSVM = (bTagSF_CSVM_.product());
  iBTV.BTAGSF_CSVT = (bTagSF_CSVT_.product());
  iBTV.MISTAGSF_CSVL = (mistagSF_CSVL_.product());
  iBTV.MISTAGSF_CSVM = (mistagSF_CSVM_.product());
  iBTV.MISTAGSF_CSVT = (mistagSF_CSVT_.product());

  for(edm::View<pat::Jet>::const_iterator jet_iter = simplejets1.begin(); jet_iter!=simplejets1.end(); ++jet_iter){
    //     if(jet_iter->pt()>50)
    //       njetscounter++;
    VHbbEvent::SimpleJet sj;
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
    //
    // add tVector
    //
    sj.tVector = getTvect(&(*jet_iter));

    Particle::LorentzVector p4Jet = jet_iter->p4();

    if(runOnMC_){

      fillScaleFactors(sj, iBTV);

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
	if(verbose_){
	  std::clog << "genJet matched Pt = " << gJp4.Pt() << std::endl;
	  std::clog << "genJet matched eta = " << gJp4.Eta() << std::endl;
	  std::clog << "genJet matched deltaR = " << gJp4.DeltaR(sj.p4) << std::endl;
	  std::clog << "genJet matched mother id = " << sj.bestMCmomid << std::endl;
	}
      }
      
    } //isMC
    hbbInfo->simpleJets.push_back(sj);
    
  }
  
  
  for(edm::View<pat::Jet>::const_iterator jet_iter = simplejets2.begin(); jet_iter!=simplejets2.end(); ++jet_iter){
    
    VHbbEvent::SimpleJet sj;
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
    //
    // add tVector
    //
    sj.tVector = getTvect(&(*jet_iter));

    Particle::LorentzVector p4Jet = jet_iter->p4();


    if(runOnMC_){

      //BTV scale factors
      fillScaleFactors(sj, iBTV);

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

    sj.flavour = subjet_iter->partonFlavour();
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
    hbbInfo->tcmet.metSig=(metsTC[0]).mEtSig();
    hbbInfo->tcmet.eLong=(metsTC[0]).e_longitudinal();
    hbbInfo->tcmet.p4=GENPTOLOR((metsTC[0]));
    if (verbose_)     std::cout <<" METTC "<<     hbbInfo->tcmet.metSig <<" " <<     hbbInfo->tcmet.sumEt<<std::endl;
  }
  
  if(mets.size()){
    hbbInfo->calomet.sumEt=(mets[0]).sumEt();
    hbbInfo->calomet.metSig=(mets[0]).mEtSig();
    hbbInfo->calomet.eLong=(mets[0]).e_longitudinal();
    hbbInfo->calomet.p4=GENPTOLOR((mets[0]));
    if (verbose_)     std::cout <<" METTC "<<     hbbInfo->calomet.metSig <<" " <<     hbbInfo->calomet.sumEt<<std::endl;
  }

  
  edm::Handle<edm::View<pat::MET> > metPFHandle;
  iEvent.getByLabel("patMETsPF",metPFHandle);
  edm::View<pat::MET> metsPF = *metPFHandle;

  if(metsPF.size()){
    hbbInfo->pfmet.sumEt=(metsPF[0]).sumEt();
    hbbInfo->pfmet.metSig=(metsPF[0]).mEtSig();
    hbbInfo->pfmet.eLong=(metsPF[0]).e_longitudinal();
    hbbInfo->pfmet.p4=GENPTOLOR((metsPF[0]));
    if (verbose_)     std::cout <<" METTC "<<     hbbInfo->pfmet.metSig <<" " <<     hbbInfo->pfmet.sumEt<<std::endl;
  }


  if(verbose_){
    std::cout << "METs: calomet "<<mets.size()<<" tcmet "<<metsTC.size()<<" pfmet "<<metsPF.size()<<std::endl;  
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


  // dimuons and dielectrons

  for( size_t i = 0; i < dimuons->size(); i++ ) {
    VHbbEvent::DiMuonInfo df;
    const Candidate & dimuonCand = (*dimuons)[ i ];
    df.p4= GENPTOLOR(dimuonCand);
    const Candidate * lep0 = dimuonCand.daughter( 0 );
    const Candidate * lep1 = dimuonCand.daughter( 1 );
    // needed to access specific methods of pat::Muon
    const pat::Muon & muonDau0 = dynamic_cast<const pat::Muon &>(*lep0->masterClone());
    const pat::Muon & muonDau1 = dynamic_cast<const pat::Muon &>(*lep1->masterClone());
    
    df.daughter1.p4=GENPTOLOR(muonDau0);
    df.daughter2.p4=GENPTOLOR(muonDau1);
    
    df.daughter1.tIso= muonDau0.trackIso();
    df.daughter2.tIso= muonDau1.trackIso();

    df.daughter1.eIso= muonDau0.ecalIso();
    df.daughter2.eIso= muonDau1.ecalIso();

    df.daughter1.hIso= muonDau0.hcalIso();
    df.daughter2.hIso= muonDau1.hcalIso();

    df.daughter1.ipDb=muonDau0.dB();
    df.daughter2.ipDb=muonDau1.dB();

    df.daughter1.ipErrDb=muonDau0.edB();
    df.daughter2.ipErrDb=muonDau1.edB();

    df.daughter1.cat=0;
    if(muonDau0.isGlobalMuon()) df.daughter1.cat|=1;
    if(muonDau0.isTrackerMuon()) df.daughter1.cat|=2;
    if(muonDau0.isStandAloneMuon()) df.daughter1.cat|=4;
    df.daughter2.cat=0;
    if(muonDau1.isGlobalMuon()) df.daughter2.cat|=1;
    if(muonDau1.isTrackerMuon()) df.daughter2.cat|=2;
    if(muonDau1.isStandAloneMuon()) df.daughter2.cat|=4;

    TrackRef trkMu1Ref = muonDau0.get<TrackRef>();
    TrackRef trkMu2Ref = muonDau1.get<TrackRef>();

    if(trkMu1Ref.isNonnull() && trkMu2Ref.isNonnull()){
      const Track* MuTrk1 = muonDau0.get<TrackRef>().get();
      const Track* MuTrk2 = muonDau1.get<TrackRef>().get();
      df.daughter1.zPVPt=MuTrk1->dz(RecVtxFirst.position());
      df.daughter1.zPVProb=MuTrk1->dz(RecVtx.position());
      df.daughter2.zPVPt=MuTrk2->dz(RecVtxFirst.position());
      df.daughter2.zPVProb=MuTrk2->dz(RecVtx.position());

      df.daughter1.nHits = MuTrk1->numberOfValidHits();
      df.daughter2.nHits = MuTrk2->numberOfValidHits();

      df.daughter1.chi2 = MuTrk1->normalizedChi2();
      df.daughter2.chi2 = MuTrk2->normalizedChi2();

      TrackRef iTrack1 = muonDau0.innerTrack();
      const reco::HitPattern& p1 = iTrack1->hitPattern();
      TrackRef iTrack2 = muonDau1.innerTrack();
      const reco::HitPattern& p2 = iTrack2->hitPattern();

      df.daughter1.nPixelHits = p1.pixelLayersWithMeasurement();
      df.daughter2.nPixelHits = p2.pixelLayersWithMeasurement();

      if(muonDau0.isGlobalMuon()){
	TrackRef gTrack = muonDau0.globalTrack();
	const reco::HitPattern& q = gTrack->hitPattern();
	df.daughter1.globNHits=q.numberOfValidMuonHits();
	df.daughter1.globChi2=gTrack.get()->normalizedChi2();
        df.daughter1.validMuStations = q. muonStationsWithValidHits();
      }
      if(muonDau1.isGlobalMuon()){
	TrackRef gTrack = muonDau1.globalTrack();
	const reco::HitPattern& q = gTrack->hitPattern();
	df.daughter2.globNHits=q.numberOfValidMuonHits();
	df.daughter2.globChi2=gTrack.get()->normalizedChi2();
        df.daughter2.validMuStations = q. muonStationsWithValidHits();
      }
  
    }
      
    hbbInfo->diMuonInfo.push_back(df);
  }

  for( size_t i = 0; i < dielectrons->size(); i++ ) {
    VHbbEvent::DiElectronInfo df;
    const Candidate & dielecCand = (*dielectrons)[ i ];

    df.p4=GENPTOLOR(dielecCand);

    // accessing the daughters of the dimuon candidate
    const Candidate * lep0 = dielecCand.daughter( 0 );
    const Candidate * lep1 = dielecCand.daughter( 1 );
    // needed to access specific methods of pat::Muon
    const pat::Electron & elecDau0 = dynamic_cast<const pat::Electron &>(*lep0->masterClone());
    const pat::Electron & elecDau1 = dynamic_cast<const pat::Electron &>(*lep1->masterClone());

    df.daughter1.p4 = GENPTOLOR(elecDau0);
    df.daughter2.p4= GENPTOLOR(elecDau1);

    df.daughter1.tIso = elecDau0.trackIso();
    df.daughter2.tIso = elecDau1.trackIso();

    df.daughter1.eIso = elecDau0.ecalIso();
    df.daughter2.eIso = elecDau1.ecalIso();

    df.daughter1.hIso = elecDau0.hcalIso();
    df.daughter2.hIso = elecDau1.hcalIso();
    
    // ids
    /*df.daughter1.id95 = elecDau0.electronID("simpleEleId95cIso");
    df.daughter1.id85 = elecDau0.electronID  ("simpleEleId85cIso");
    df.daughter1.id70 = elecDau0.electronID  ("simpleEleId70cIso");
    df.daughter1.id95r = elecDau0.electronID ("simpleEleId95relIso");
    df.daughter1.id85r = elecDau0.electronID ("simpleEleId85relIso");
    df.daughter1.id70r = elecDau0.electronID ("simpleEleId70relIso");


    df.daughter2.id95 = elecDau1.electronID("simpleEleId95cIso");
    df.daughter2.id85 = elecDau1.electronID  ("simpleEleId85cIso");
    df.daughter2.id70 = elecDau1.electronID  ("simpleEleId70cIso");
    df.daughter2.id95r = elecDau1.electronID ("simpleEleId95relIso");
    df.daughter2.id85r = elecDau1.electronID ("simpleEleId85relIso");
    df.daughter2.id70r = elecDau1.electronID ("simpleEleId70relIso");
*/
    hbbInfo->diElectronInfo.push_back(df);
    
  }
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
void HbbAnalyzerNew::fillScaleFactors(VHbbEvent::SimpleJet sj, BTV_SF iSF){


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

//define this as a plug-in
DEFINE_FWK_MODULE(HbbAnalyzerNew);
