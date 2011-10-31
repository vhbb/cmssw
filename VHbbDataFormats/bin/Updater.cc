#include <TH1F.h>
#include <TH3F.h>
#include <TH2F.h>
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//btagging
//#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
//#include "VHbbAnalysis/VHbbDataFormats/interface/BTagWeight.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/TriggerWeight.h"

#include <sstream>
#include <string>

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

  TriggerWeight triggerWeight(ana);

  edm::LumiReWeighting   lumiWeights;
  edm::LumiReWeighting   lumiWeights2011B;

  if(redoPU)
  {
  std::string PUmcfileName_ = in.getParameter<std::string> ("PUmcfileName") ;
  std::string PUmcfileName2011B_ = in.getParameter<std::string> ("PUmcfileName2011B") ;
  std::string PUdatafileName_ = in.getParameter<std::string> ("PUdatafileName") ;
  std::string PUdatafileName2011B_ = in.getParameter<std::string> ("PUdatafileName2011B") ;
  lumiWeights = edm::LumiReWeighting(PUmcfileName_,PUdatafileName_ , "pileup", "pileup");
  lumiWeights2011B = edm::LumiReWeighting(PUmcfileName2011B_,PUdatafileName2011B_ , "pileup", "pileup");


  //lumiWeights2011B.weight3D_init(); // generate the weights the fisrt time;
  lumiWeights2011B.weight3D_init("Weight3D.root");

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
   Int_t           nhJets;
   Float_t         hJet_pt[50];
   Float_t         hJet_eta[50];
   Int_t           naJets;
   Float_t         aJet_pt[50];
   Float_t         aJet_eta[50];
   Int_t           Vtype;
   METInfo         MET;

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

   }


   
//   Float_t weightTrigUpdate;
//   newtree->Branch("weightTrigUpdate"        , &weightTrigUpdate          ,  "weightTrigUpdate/F");

   for (Int_t i=0;i<nentries; i++) {
         oldtree->GetEntry(i);

         if(redoPU) {
          PUweight =  lumiWeights.weight( PU0 );
          PUweight2011B = lumiWeights2011B.weight3D( PUm1, PU0, PUp1);
        }
  


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
        if( Vtype == 4 ){
          float weightTrig1 = triggerWeight.scaleMetHLT(MET.et);
          weightTrig = weightTrig1;
          weightTrigMET = weightTrig1;

        }




      newtree->Fill();
   }
   newtree->Print();
   newtree->AutoSave();
   delete oldfile;
   delete newfile;
}
