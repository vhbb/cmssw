// -*- C++ -*-
//
// Package:    MCatNLOAnalyzer
// Class:      MCatNLOAnalyzer
// 
/**\class MCatNLOAnalyzer MCatNLOAnalyzer.cc BAnalysis/MCatNLOAnalyzer/src/MCatNLOAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lukas Wehrli (IPP/ETHZ) [wehrlilu]
//         Created:  Thu Aug 19 08:50:46 CEST 2010
// $Id: MCatNLOAnalyzer.cc,v 1.1 2010/09/08 11:24:40 wehrlilu Exp $
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

#include "AnalysisDataFormats/BAnalysis/interface/SimBHadron.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TH1F.h"
#include "TH2F.h"

class MCatNLOAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MCatNLOAnalyzer(const edm::ParameterSet&);
      ~MCatNLOAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      edm::Service<TFileService> fs_;
      edm::InputTag simbsrc_;
      double maxEta, minPt;
      double maxPThat, minPThat;  
      double ptcutLeadingJet; 
      TH1F *hdR, *hdEta, *hdPhi, *hpthat;
      

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MCatNLOAnalyzer::MCatNLOAnalyzer(const edm::ParameterSet& iConfig):
  simbsrc_(iConfig.getParameter<edm::InputTag>("simbsrc")),
  maxEta(iConfig.getUntrackedParameter<double>("maxEta",2.4)),
  minPt(iConfig.getUntrackedParameter<double>("minPt",15)),
  maxPThat(iConfig.getUntrackedParameter<double>("maxPThat",0.0)),
  minPThat(iConfig.getUntrackedParameter<double>("minPThat",10000.0)),
  ptcutLeadingJet(iConfig.getUntrackedParameter<double>("ptcutLeadingJet",84.0))

{
   //now do what ever initialization is needed

}


MCatNLOAnalyzer::~MCatNLOAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
MCatNLOAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;

   bool takeEvent = true; 

   /////////////////////////////////////////////
   //--------------GenParticles
   Handle<GenParticleCollection> GPCollection;
   iEvent.getByLabel("genParticles", GPCollection);
   const GenParticleCollection & gpc = *(GPCollection.product());
   

   //pthat defined here as pt of hardest b quark
   double thispt = -1; 
   for(unsigned int i=0; i<gpc.size(); i++){
     int pdgid = gpc[i].pdgId();
     int pdid = abs(pdgid);
     if(pdid==5){
	double pt = gpc[i].pt(); 
	if(pt>thispt) thispt = pt;
     }
   }

   //pthat cuts (for mcatnlo samples)
   if(thispt<minPThat) {
     std::cout << "WARNING ABOUT LOWER PT BOUNDARY\n";
     std::cout << "minpthat " << minPThat << " thispt " << thispt << "\n";
     takeEvent = false;
   }
   if(thispt>maxPThat) takeEvent = false;

   if(!takeEvent) {
// 	std::cout << "cont\n";
	return; 
   } 

   hpthat->Fill(thispt);

   //------------------SimBHadrons------------------------------
   Handle<SimBHadronCollection> SBHC;
   iEvent.getByLabel(simbsrc_, SBHC);
   SimBHadronCollection sbhc = *(SBHC.product());


   //------------------GenJets------------------------------
   Handle<GenJetCollection> GJC;
   iEvent.getByLabel("ak5GenJets", GJC);
   GenJetCollection gjc = *(GJC.product());
   double ptHardestGJ = -1, etaHardestGJ=-1;
   for(unsigned int i=0; i<gjc.size(); i++){
       if(gjc[i].pt()>ptHardestGJ && fabs(gjc[i].eta())<3.0) {
         ptHardestGJ = gjc[i].pt();
         etaHardestGJ = gjc[i].eta();
       }
   }

   //leading jet cut at ptcutLeadingJet ((84 corresponding to cut of hlt_jet30u trigger))
   if(ptHardestGJ>ptcutLeadingJet){
// 	std::cout << "LJ ok";
	if(sbhc.size()==2){
	  //kinematic cuts on b: as done for data and other mc event generators
	   if(fabs(sbhc[0].eta())<maxEta && fabs(sbhc[1].eta())<maxEta && sbhc[0].pt()>minPt && sbhc[1].pt()>minPt){
 	        double dRbb = deltaR(sbhc[0].p4(),sbhc[1].p4());
		double dEtabb = sbhc[0].p4().Eta()-sbhc[1].p4().Eta();
		double dPhibb = deltaPhi(sbhc[0].p4().Phi(), sbhc[1].p4().Phi());
		hdR->Fill(dRbb); 
	        hdPhi->Fill(fabs(dPhibb)); 
		hdEta->Fill(fabs(dEtabb));
// 		std::cout << " 2B ok";
	   }
	}
   }
//    std::cout << "\n";

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
MCatNLOAnalyzer::beginJob()
{
  hdR = fs_->make<TH1F>("hdR_mcatnlo","#Delta R between two B hadrons",60,0,6);
  hdEta = fs_->make<TH1F>("hdEta_mcatnlo","#Delta #eta between two B hadrons",100,0,10);
  hdPhi = fs_->make<TH1F>("hdPhi_mcatnlo","#Delta #phi between two B hadrons",32,0,3.2);
  hpthat = fs_->make<TH1F>("hpthat","pthat after cut",100,0,400);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MCatNLOAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MCatNLOAnalyzer);
