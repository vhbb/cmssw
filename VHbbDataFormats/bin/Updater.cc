#include <TH1F.h>
#include <TH3F.h>
#include <TH2F.h>
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//btagging
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
//#include "VHbbAnalysis/VHbbDataFormats/interface/BTagWeight.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/HbbCandidateFinderAlgo.h"
#include "VHbbAnalysis/VHbbDataFormats/src/HbbCandidateFinderAlgo.cc"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidateTools.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/TriggerWeight.h"

#include <sstream>
#include <string>
#define nhJets_MAX 2
#define naJets_MAX 20

#define MAXJ 30
#define MAXL 10
#define MAXB 10
typedef struct
{
  float et;
  float sumet;
  float sig;
  float phi;
} METInfo;


int main(int argc, char* argv[])
{
  gROOT->Reset();

  // parse arguments
  if ( argc < 2 ) {
    return 0;
  }
  // get the python configuration
  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in  = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput" );
  const edm::ParameterSet& out = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteOutput");
  const edm::ParameterSet& ana = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("Analyzer");
  // now get each parameter
//  int maxEvents_( in.getParameter<int>("maxEvents") );
//  int skipEvents_( in.getParameter<int>("skipEvents") );
//  unsigned int outputEvery_( in.getParameter<unsigned int>("outputEvery") );
  std::string outputFile_( out.getParameter<std::string>("fileName" ) );
  std::string inputFile_( in.getParameter<std::string>("fileName" ) );
  bool replaceWeights( ana.getParameter<bool>("replaceWeights" ) );
  bool redoPU( ana.getParameter<bool>("redoPU" ) );
  bool redoTrigger( ana.getParameter<bool>("redoTrigger" ) );
  bool redoHiggs( ana.getParameter<bool>("redoHiggs" ) );

  std::string Weight3DfileName_ = in.getParameter<std::string> ("Weight3DfileName") ;

  TriggerWeight triggerWeight(ana);

  edm::LumiReWeighting   lumiWeights;
  edm::Lumi3DReWeighting   lumiWeights2011B;

  HbbCandidateFinderAlgo finder(false, 20, true);
  VHbbCandidateTools selector(false);


  if(redoPU)
  {
  std::string PUmcfileName_ = in.getParameter<std::string> ("PUmcfileName") ;
  std::string PUmcfileName2011B_ = in.getParameter<std::string> ("PUmcfileName2011B") ;
  std::string PUdatafileName_ = in.getParameter<std::string> ("PUdatafileName") ;
  std::string PUdatafileName2011B_ = in.getParameter<std::string> ("PUdatafileName2011B") ;
  lumiWeights = edm::LumiReWeighting(PUmcfileName_,PUdatafileName_ , "pileup", "pileup");
//  lumiWeights2011B = edm::Lumi3DReWeighting(PUmcfileName2011B_,PUdatafileName2011B_ , "pileup", "pileup");


  //lumiWeights2011B.weight3D_init(); // generate the weights the fisrt time;
  /*                 if(Weight3DfileName_!="")
                      { lumiWeights2011B.weight3D_init(Weight3DfileName_.c_str()); }
                   else
                      {
                        lumiWeights2011B.weight3D_init(73.5/68.); // generate the weights the fisrt time;
                      }

*/
  }



   //Get old file, old tree and set top branch address
   TFile *oldfile = new TFile(inputFile_.c_str());
   TTree *oldtree = (TTree*)oldfile->Get("tree");
   TH1F *  count = (TH1F*)oldfile->Get("Count");
   TH1F *  countWithPU = (TH1F*)oldfile->Get("CountWithPU");
   TH1F *  countWithPU2011B = (TH1F*)oldfile->Get("CountWithPU2011B");
   TH3F *  input3DPU = (TH3F*)oldfile->Get("Input3DPU");


   Int_t nentries = (Int_t)oldtree->GetEntries();
   Int_t           nvlep;
   Float_t         vLepton_pt[50];
   Float_t         vLepton_eta[50];
/*   Int_t           nhJets;
   Float_t         hJet_pt[50];
   Float_t         hJet_eta[50];
   Int_t           naJets;
   Float_t         aJet_pt[50];
   Float_t         aJet_eta[50];
*/   Int_t           Vtype;
   METInfo         MET;
   Int_t           nhJets;
   Int_t           naJets;
   Float_t         hJet_pt[nhJets_MAX];
   Float_t         hJet_eta[nhJets_MAX];
   Float_t         hJet_phi[nhJets_MAX];
   Float_t         hJet_e[nhJets_MAX];
   Float_t         hJet_csv[nhJets_MAX];
   Float_t         hJet_cosTheta[nhJets_MAX];
   Int_t           hJet_numTracksSV[nhJets_MAX];
   Float_t         hJet_chf[nhJets_MAX];
   Float_t         hJet_nhf[nhJets_MAX];
   Float_t         hJet_cef[nhJets_MAX];
   Float_t         hJet_nef[nhJets_MAX];
   Float_t         hJet_nch[nhJets_MAX];
   Float_t         hJet_nconstituents[nhJets_MAX];
   Float_t         hJet_flavour[nhJets_MAX];
   Float_t         hJet_genPt[nhJets_MAX];
   Float_t         hJet_genEta[nhJets_MAX];
   Float_t         hJet_genPhi[nhJets_MAX];
   Float_t         hJet_JECUnc[nhJets_MAX];
   Float_t         hJet_vtxMass[nhJets_MAX];
   Float_t         hJet_vtx3dL[nhJets_MAX];
   Float_t         hJet_vtx3deL[nhJets_MAX];
   UChar_t         hJet_id[nhJets_MAX];
   Float_t         aJet_pt[naJets_MAX];
   Float_t         aJet_eta[naJets_MAX];
   Float_t         aJet_phi[naJets_MAX];
   Float_t         aJet_e[naJets_MAX];
   Float_t         aJet_csv[naJets_MAX];
   Float_t         aJet_cosTheta[naJets_MAX];
   Int_t           aJet_numTracksSV[naJets_MAX];
   Float_t         aJet_chf[naJets_MAX];
   Float_t         aJet_nhf[naJets_MAX];
   Float_t         aJet_cef[naJets_MAX];
   Float_t         aJet_nef[naJets_MAX];
   Float_t         aJet_nch[naJets_MAX];
   Float_t         aJet_nconstituents[naJets_MAX];
   Float_t         aJet_flavour[naJets_MAX];
   Float_t         aJet_genPt[naJets_MAX];
   Float_t         aJet_genEta[naJets_MAX];
   Float_t         aJet_genPhi[naJets_MAX];
   Float_t         aJet_JECUnc[naJets_MAX];
   Float_t         aJet_vtxMass[naJets_MAX];
   Float_t         aJet_vtx3dL[naJets_MAX];
   Float_t         aJet_vtx3deL[naJets_MAX];
   UChar_t         aJet_id[naJets_MAX];

   oldtree->SetBranchAddress("nvlep", &nvlep);
   oldtree->SetBranchAddress("nhJets", &nhJets);
   oldtree->SetBranchAddress("naJets", &naJets);

   oldtree->SetBranchAddress("vLepton_pt", vLepton_pt);
   oldtree->SetBranchAddress("vLepton_eta", vLepton_eta);

   oldtree->SetBranchAddress("hJet_pt", hJet_pt);
   oldtree->SetBranchAddress("hJet_eta", hJet_eta);

   oldtree->SetBranchAddress("aJet_pt", aJet_pt);
   oldtree->SetBranchAddress("aJet_eta", aJet_eta);

   oldtree->SetBranchAddress("Vtype", &Vtype);
   oldtree->SetBranchAddress("MET", &MET);

   // Trigger weights
   Float_t         weightTrig;
   Float_t         weightTrigMay;
   Float_t         weightTrigV4;
   Float_t         weightTrigMET;
   Float_t         weightTrigOrMu30;
   Float_t         weightEleRecoAndId;
   Float_t         weightEleTrigJetMETPart;
   Float_t         weightEleTrigElePart;
   if(replaceWeights)
   {
   std::cout << "Replacing the weights in the same branch names" << std::endl;
   oldtree->SetBranchAddress("weightTrig", &weightTrig);
   oldtree->SetBranchAddress("weightTrigMay", &weightTrigMay);
   oldtree->SetBranchAddress("weightTrigV4", &weightTrigV4);
   oldtree->SetBranchAddress("weightTrigMET", &weightTrigMET);
   oldtree->SetBranchAddress("weightTrigOrMu30", &weightTrigOrMu30);
   oldtree->SetBranchAddress("weightEleRecoAndId", &weightEleRecoAndId);
   oldtree->SetBranchAddress("weightEleTrigJetMETPart", &weightEleTrigJetMETPart);
   oldtree->SetBranchAddress("weightEleTrigElePart", &weightEleTrigElePart);
   } 
 //Jet variables
 //  bool redoHiggs = true;
  if(redoHiggs)
  {
//   oldtree->SetBranchAddress("nhJets", &nhJets);
//   oldtree->SetBranchAddress("naJets", &naJets);
//   oldtree->SetBranchAddress("hJet_pt", hJet_pt);
//   oldtree->SetBranchAddress("hJet_eta", hJet_eta);
   oldtree->SetBranchAddress("hJet_phi", hJet_phi);
   oldtree->SetBranchAddress("hJet_e", hJet_e);
   oldtree->SetBranchAddress("hJet_csv", hJet_csv);
   oldtree->SetBranchAddress("hJet_cosTheta", hJet_cosTheta);
   oldtree->SetBranchAddress("hJet_numTracksSV", hJet_numTracksSV);
   oldtree->SetBranchAddress("hJet_chf", hJet_chf);
   oldtree->SetBranchAddress("hJet_nhf", hJet_nhf);
   oldtree->SetBranchAddress("hJet_cef", hJet_cef);
   oldtree->SetBranchAddress("hJet_nef", hJet_nef);
   oldtree->SetBranchAddress("hJet_nch", hJet_nch);
   oldtree->SetBranchAddress("hJet_nconstituents", hJet_nconstituents);
   oldtree->SetBranchAddress("hJet_flavour", hJet_flavour);
   oldtree->SetBranchAddress("hJet_genPt", hJet_genPt);
   oldtree->SetBranchAddress("hJet_genEta", hJet_genEta);
   oldtree->SetBranchAddress("hJet_genPhi", hJet_genPhi);
   oldtree->SetBranchAddress("hJet_JECUnc", hJet_JECUnc);
   oldtree->SetBranchAddress("hJet_vtxMass", hJet_vtxMass);
   oldtree->SetBranchAddress("hJet_vtx3dL", hJet_vtx3dL);
   oldtree->SetBranchAddress("hJet_vtx3deL", hJet_vtx3deL);
   oldtree->SetBranchAddress("hJet_id", hJet_id);

//   oldtree->SetBranchAddress("aJet_pt", aJet_pt);
//   oldtree->SetBranchAddress("aJet_eta", aJet_eta);
   oldtree->SetBranchAddress("aJet_phi", aJet_phi);
   oldtree->SetBranchAddress("aJet_e", aJet_e);
   oldtree->SetBranchAddress("aJet_csv", aJet_csv);
   oldtree->SetBranchAddress("aJet_cosTheta", aJet_cosTheta);
   oldtree->SetBranchAddress("aJet_numTracksSV", aJet_numTracksSV);
   oldtree->SetBranchAddress("aJet_chf", aJet_chf);
   oldtree->SetBranchAddress("aJet_nhf", aJet_nhf);
   oldtree->SetBranchAddress("aJet_cef", aJet_cef);
   oldtree->SetBranchAddress("aJet_nef", aJet_nef);
   oldtree->SetBranchAddress("aJet_nch", aJet_nch);
   oldtree->SetBranchAddress("aJet_nconstituents", aJet_nconstituents);
   oldtree->SetBranchAddress("aJet_flavour", aJet_flavour);
   oldtree->SetBranchAddress("aJet_genPt", aJet_genPt);
   oldtree->SetBranchAddress("aJet_genEta", aJet_genEta);
   oldtree->SetBranchAddress("aJet_genPhi", aJet_genPhi);
   oldtree->SetBranchAddress("aJet_JECUnc", aJet_JECUnc);
   oldtree->SetBranchAddress("aJet_vtxMass", aJet_vtxMass);
   oldtree->SetBranchAddress("aJet_vtx3dL", aJet_vtx3dL);
   oldtree->SetBranchAddress("aJet_vtx3deL", aJet_vtx3deL);
   oldtree->SetBranchAddress("aJet_id", aJet_id);

  } 
 
   //Pileup Info
   Float_t         PUweight;
   Float_t         PUweight2011B;
   Float_t  PU0,PUp1,PUm1;
   if(redoPU)
   {
    oldtree->SetBranchAddress("PU0", &PU0);
    oldtree->SetBranchAddress("PUp1", &PUp1);
    oldtree->SetBranchAddress("PUm1", &PUm1);
    oldtree->SetBranchAddress("PUweight", &PUweight);
    oldtree->SetBranchAddress("PUweight2011B", &PUweight2011B);
   }


   //Create a new file + a clone of old tree in new file + clone of the histos + additional/updated weights

   TFile *newfile = new TFile(outputFile_.c_str(),"RECREATE");
   TTree *newtree = oldtree->CloneTree(0);
   if(count) count->Clone()->Write();
   if(input3DPU) input3DPU->Clone()->Write();

   if(redoPU)
   { 
     if(countWithPU) countWithPU->Clone("CountWithPU_OLD")->Write();
     if(countWithPU2011B) countWithPU2011B->Clone("CountWithPU2011B_OLD")->Write();

//recompute the normalization
     countWithPU = new TH1F("CountWithPU","CountWithPU", 1,0,2 );
     countWithPU2011B = new TH1F("CountWithPU2011B","CountWithPU2011B", 1,0,2 );
     for(int ix=1;ix<=input3DPU->GetNbinsX();ix++)
      for(int iy=1;iy<=input3DPU->GetNbinsY();iy++)
       for(int iz=1;iz<=input3DPU->GetNbinsZ();iz++)
        {
          Float_t nev=input3DPU->GetBinContent(ix,iy,iz);
          PUweight =  lumiWeights.weight( iy-1 );  // bin 1 is [-0.5,0.5]
          PUweight2011B = lumiWeights2011B.weight3D( ix-1, iy-1, iz-1);
          countWithPU->Fill(1,PUweight*nev);
          countWithPU2011B->Fill(1,PUweight2011B*nev);
        }

     countWithPU->Write();
     countWithPU2011B->Write();
 
   }
   else
   { //Just clone the old ones
     if(countWithPU) countWithPU->Clone()->Write();
     if(countWithPU2011B) countWithPU2011B->Clone()->Write();
   }
    
  if(!replaceWeights)
   {
   std::cout << "Creating new branch names with _up postfix" << std::endl;
   newtree->Branch("weightTrig_up", &weightTrig, "weightTrig_up/F");
   newtree->Branch("weightTrigMay_up",  &weightTrigMay,"weightTrigMay/F");
   newtree->Branch("weightTrigV4_up",  &weightTrigV4,"weightTrigV4/F");
   newtree->Branch("weightTrigMET_up",  &weightTrigMET,"weightTrigMET/F");
   newtree->Branch("weightTrigOrMu30_up",  &weightTrigOrMu30,"weightTrigOrMu30/F");
   newtree->Branch("weightEleRecoAndId_up",  &weightEleRecoAndId,"weightEleRecoAndId/F");
   newtree->Branch("weightEleTrigJetMETPart_up",  &weightEleTrigJetMETPart,"weightEleTrigJetMETPart/F");
   newtree->Branch("weightEleTrigElePart_up",  &weightEleTrigElePart,"weightEleTrigElePart/F");

/* TODO:
  _outTree->Branch("weightEleTrigEleAugPart"        , &weightEleTrigEleAugPart          ,  "weightEleTrigEleAugPart/F");
  
  _outTree->Branch("weightTrigMET80"        , &weightTrigMET80          , "weightTrigMET80/F");
  _outTree->Branch("weightTrigMET100"        , &weightTrigMET100          , "weightTrigMET100/F");
  _outTree->Branch("weightTrig2CJet20"        , &weightTrig2CJet20          , "weightTrig2CJet20/F");
  _outTree->Branch("weightTrigMET150"        , &weightTrigMET150          , "weightTrigMET150/F");
  _outTree->Branch("weightTrigMET802CJet"        , &weightTrigMET802CJet          , "weightTrigMET802CJet/F");
  _outTree->Branch("weightTrigMET1002CJet"        , &weightTrigMET1002CJet          , "weightTrigMET1002CJet/F");
  _outTree->Branch("weightTrigMETLP"        , &weightTrigMETLP          , "weightTrigMETLP/F");
*/


   }


   
//   Float_t weightTrigUpdate;
//   newtree->Branch("weightTrigUpdate"        , &weightTrigUpdate          ,  "weightTrigUpdate/F");

   for (Int_t i=0;i<nentries; i++) {
         oldtree->GetEntry(i);

         if(redoPU) {
          PUweight =  lumiWeights.weight( PU0 );
          PUweight2011B = lumiWeights2011B.weight3D( PUm1, PU0, PUp1);
        }
  
        if(redoHiggs)
        {

std::vector<VHbbEvent::SimpleJet> jets;
for(int i =0; i < nhJets; i++)
{
  VHbbEvent::SimpleJet j;
  float scale=1.0;
  j.p4.SetPtEtaPhiE(hJet_pt[i]*scale,hJet_eta[i],hJet_phi[i],hJet_e[i]*scale);
  j.csv=hJet_csv[i];
  j.vtxMass=hJet_numTracksSV[i] ;
  j.vtxMass=hJet_vtxMass[i];
  j.vtx3dL=hJet_vtx3dL[i] ;
  j.vtx3deL=hJet_vtx3deL[i] ;
  j.chargedHadronEFraction=hJet_chf[i];
  j.neutralHadronEFraction=hJet_nhf[i]  ;
  j.chargedEmEFraction=hJet_cef[i]  ;
  j.neutralEmEFraction=hJet_nef[i]  ;
  j.nConstituents=hJet_nconstituents[i]  ;
  j.ntracks=hJet_nch[i];
/*  j.SF_CSVL=hJet_SF_CSVL[i];
  j.SF_CSVM=hJet_SF_CSVM[i];
  j.SF_CSVT=hJet_SF_CSVT[i];
  j.SF_CSVLerr=hJet_SF_CSVLerr[i];
  j.SF_CSVMerr=hJet_SF_CSVMerr[i];
  j.SF_CSVTerr=hJet_SF_CSVTerr[i];
*/
  j.flavour=hJet_flavour[i];
  j.bestMCp4.SetPtEtaPhiE(hJet_genPt[i], hJet_genEta[i], hJet_genPhi[i],  hJet_genPt[i]); //E not saved!
//j.bestMCp4.Phi()=hJet_genPhi[i];
  j.jecunc=hJet_JECUnc[i];
  jets.push_back(j);
}


for(int i =0; i < naJets; i++)
{
  VHbbEvent::SimpleJet j;
  float scale=1.0;
  j.p4.SetPtEtaPhiE(aJet_pt[i]*scale,aJet_eta[i],aJet_phi[i],aJet_e[i]*scale);
  j.csv=aJet_csv[i];
  j.vtxMass=aJet_numTracksSV[i] ;
  j.vtxMass=aJet_vtxMass[i];
  j.vtx3dL=aJet_vtx3dL[i] ;
  j.vtx3deL=aJet_vtx3deL[i] ;
  j.chargedHadronEFraction=aJet_chf[i];
  j.neutralHadronEFraction=aJet_nhf[i]  ;
  j.chargedEmEFraction=aJet_cef[i]  ;
  j.neutralEmEFraction=aJet_nef[i]  ;
  j.nConstituents=aJet_nconstituents[i]  ;
  j.ntracks=aJet_nch[i];
/*  j.SF_CSVL=aJet_SF_CSVL[i];
  j.SF_CSVM=aJet_SF_CSVM[i];
  j.SF_CSVT=aJet_SF_CSVT[i];
  j.SF_CSVLerr=aJet_SF_CSVLerr[i];
  j.SF_CSVMerr=aJet_SF_CSVMerr[i];
  j.SF_CSVTerr=aJet_SF_CSVTerr[i];
*/
  j.flavour=aJet_flavour[i];
  j.bestMCp4.SetPtEtaPhiE(aJet_genPt[i],aJet_genEta[i],aJet_genPhi[i],aJet_genPt[i]);
  j.jecunc=aJet_JECUnc[i];
  jets.push_back(j);
}

/*if (useHighestPtHiggs == false){
    foundJets = findDiJets(jets,j1,j2,addJets) ;
  }else{*/
  std::vector<VHbbEvent::SimpleJet> addJets;
  VHbbEvent::SimpleJet j1;
  VHbbEvent::SimpleJet j2;
  
   bool foundJets = finder.findDiJetsHighestPt(jets,j1,j2,addJets) ;
/*  TVector3 higgsBoost;
  higgsBoost = (temp.H.p4).BoostVector();
  temp.H.helicities.clear();
  temp.H.helicities.push_back(selector.getHelicity(j1,higgsBoost));
  temp.H.helicities.push_back(selector.getHelicity(j2,higgsBoost));
  temp.H.deltaTheta = selector.getDeltaTheta(j1,j2);
  temp.additionalJets = addJets;
          hJets.cosTheta[0]=  vhCand.H.helicities[0];
        hJets.cosTheta[1]=  vhCand.H.helicities[1];*/

  }
 
 if(redoTrigger) 
 {    
         std::vector<float> jet30eta;
         std::vector<float> jet30pt;
         for( int j = 0 ; j < nhJets; j++) if(hJet_pt[j]>30 ) { jet30eta.push_back(hJet_eta[j]); jet30pt.push_back(hJet_pt[j]); } 
         for( int j = 0 ; j < naJets; j++) if(aJet_pt[j]>30 ) { jet30eta.push_back(aJet_eta[j]); jet30pt.push_back(aJet_pt[j]); } 

	 if(Vtype == 0 ){
          float cweightID = triggerWeight.scaleMuID(vLepton_pt[0],vLepton_eta[0]) * triggerWeight.scaleMuID(vLepton_pt[1],vLepton_eta[1]) ;
          float weightTrig1 = triggerWeight.scaleMuIsoHLT(vLepton_pt[0],vLepton_eta[0]);
          float weightTrig2 = triggerWeight.scaleMuIsoHLT(vLepton_pt[1],vLepton_eta[1]);
          float cweightTrig = weightTrig1 + weightTrig2 - weightTrig1*weightTrig2;
          weightTrig = cweightID * cweightTrig;

        }
        if( Vtype == 1 ){
          std::vector<float> pt,eta;
          pt.push_back(vLepton_pt[0]); eta.push_back(vLepton_eta[0]);
          pt.push_back(vLepton_pt[1]); eta.push_back(vLepton_eta[1]);
          weightEleRecoAndId=triggerWeight.scaleID95Ele(vLepton_pt[0],vLepton_eta[0]) * triggerWeight.scaleRecoEle(vLepton_pt[0],vLepton_eta[0]) *
            triggerWeight.scaleID95Ele(vLepton_pt[1],vLepton_eta[1]) * triggerWeight.scaleRecoEle(vLepton_pt[1],vLepton_eta[1]);
          weightEleTrigElePart = triggerWeight.scaleDoubleEle17Ele8(pt,eta);
          //REMOVE FOR "float" for newer ntupler and add branch
          float weightEleTrigEleAugPart = triggerWeight.scaleDoubleEle17Ele8Aug(pt,eta);
          weightTrig = (weightEleTrigElePart*1.14+weightEleTrigEleAugPart*0.98 )/2.12 * weightEleRecoAndId;



        }
        if(Vtype == 2 ){
          float cweightID = triggerWeight.scaleMuID(vLepton_pt[0],vLepton_eta[0]);
          float weightTrig1 = triggerWeight.scaleMuIsoHLT(vLepton_pt[0],vLepton_eta[0]);
          float cweightTrig = weightTrig1;
          weightTrig = cweightID * cweightTrig;
          float weightTrig1OrMu30 = triggerWeight.scaleMuOr30IsoHLT(vLepton_pt[0],vLepton_eta[0]);
          weightTrigOrMu30 = cweightID*weightTrig1OrMu30;

        }
        if( Vtype == 3 ){
          weightTrigMay = triggerWeight.scaleSingleEleMay(vLepton_pt[0],vLepton_eta[0]);
          weightTrigV4 = triggerWeight.scaleSingleEleV4(vLepton_pt[0],vLepton_eta[0]);
          weightEleRecoAndId=triggerWeight.scaleID80Ele(vLepton_pt[0],vLepton_eta[0]) * triggerWeight.scaleRecoEle(vLepton_pt[0],vLepton_eta[0]);
          weightEleTrigJetMETPart=triggerWeight.scaleJet30Jet25(jet30pt,jet30eta)*triggerWeight.scalePFMHTEle(MET.et);
          weightEleTrigElePart= weightTrigV4; //this is for debugging only, checking only the V4 part

          weightTrigMay*=weightEleRecoAndId;
          weightTrigV4*=weightEleRecoAndId;
          weightTrigV4*=weightEleTrigJetMETPart;
//        weightTrig = weightTrigMay * 0.187 + weightTrigV4 * (1.-0.187); //FIXME: use proper lumi if we reload 2.fb
          weightTrig = (weightTrigMay * 0.215 + weightTrigV4 * 1.915)/ 2.13; //FIXME: use proper lumi if we reload 2.fb


        }
/* TODO
 if(isMC_)
{
        weightTrigMET80 =  triggerWeight.scaleMET80(MET.et);
        weightTrigMET100 =  triggerWeight.scaleMET80(MET.et);
        weightTrig2CJet20 = triggerWeight.scale2CentralJet( jet10pt, jet10eta);
        weightTrigMET150 = triggerWeight.scaleMET150(MET.et);
        weightTrigMET802CJet= weightTrigMET80 * weightTrig2CJet20;
        weightTrigMET1002CJet= weightTrigMET100 * weightTrig2CJet20;
}
        if( Vtype == VHbbCandidate::Znn ){
          nvlep=0;
          float weightTrig1 = triggerWeight.scaleMetHLT(vhCand.V.mets.at(0).p4.Pt());
          weightTrigMETLP = weightTrig1;
          weightTrig = weightTrigMET150 + weightTrigMET802CJet  - weightTrigMET802CJet*weightTrigMET150;
        }
*/

  }


      newtree->Fill();
   }
   newtree->Print();
   newtree->AutoSave();
   delete oldfile;
   delete newfile;
}
