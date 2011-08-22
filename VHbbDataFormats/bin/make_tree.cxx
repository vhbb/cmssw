#include <TSystem.h>
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/TriggerReader.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/Histos.h"
#include <iostream>
#include <fstream>
#include "VHbbAnalysis/VHbbDataFormats/interface/Cuts200X.h"
#include "DataFormats/Math/interface/deltaR.h"


int main( int argc, char ** argv ){

  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  gSystem->Load("libDataFormatsFWLite");
  AutoLibraryLoader::enable();

  std::string name(argv[2]);

  TFile * TreeFile = new TFile((name+"TreeFile.root").c_str(),"RECREATE");
  TTree * tree = new TTree("treeMVA","treeMVA");
  std::clog << "Booking the tree ... "<< std::endl;

  //book the tree for MVA
  float bbMass,bbPt,btag1,btag2,NaddJet,DeltaRbb,helicity,DeltaPhiVH,bPt1,bPt2,VMass,VPt,pullAngle,DeltaEtabb,deltaPhipfMETjet1,deltaPhipfMETjet2,pfMET,pfMETsig;
  
  tree->Branch("bbMass",&bbMass,"bbMass/F");
  tree->Branch("bbPt",&bbPt,"bbPt/F");
  tree->Branch("btag1",&btag1,"btag1/F");
  tree->Branch("btag2",&btag2,"btag2/F");
  tree->Branch("NaddJet",&NaddJet,"NaddJet/F");
  tree->Branch("DeltaRbb", &DeltaRbb, "DeltaRbb/F" );
  tree->Branch("helicity", &helicity, "helicity/F" );
  tree->Branch("DeltaPhiVH", &DeltaPhiVH ,"DeltaPhiVH/F");
  tree->Branch("bPt1", &bPt1, "bPt1/F" );
  tree->Branch("bPt2", &bPt2, "bPt2/F" );
  tree->Branch("VMass", &VMass, "VMass/F" );
  tree->Branch("VPt", &VPt, "VPt/F" );
  tree->Branch("pullAngle", &pullAngle, "pullAngle/F" );
  tree->Branch("DeltaEtabb",&DeltaEtabb, "DeltaEtabb/F");
  //for the precuts
  tree->Branch("deltaPhipfMETjet1",&deltaPhipfMETjet1,"deltaPhipfMETjet1/F");
  tree->Branch("deltaPhipfMETjet2",&deltaPhipfMETjet2,"deltaPhipfMETjet2/F");
  tree->Branch("pfMET",&pfMET,"pfMET/F");
  tree->Branch("pfMETsig",&pfMETsig,"pfMETsig/F");

  std::string fl(argv[1]);
  std::vector<std::string> inputFiles_;
  std::ifstream in(fl.c_str());
  while(!in.eof())
  {
   std::string line;
   in>> line;
   inputFiles_.push_back(line);
  }
  TriggerReader trigger;
 
  //Loop on all files
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
    std::cout << iFile << std::endl;
    // open input file (can be located on castor)
    TFile* f = TFile::Open(inputFiles_[iFile].c_str());
    if(f==0) continue;
    
    fwlite::Event ev(f);
    unsigned int event = 0; 
    //Loop on all events of this file
    for( ev.toBegin(); !ev.atEnd(); ++ev) {
      
      event++;
      //      std::clog << "Entering event "<< event << std::endl;
      fwlite::Handle< std::vector<VHbbCandidate> > vhbbCandHandle; 
      vhbbCandHandle.getByLabel(ev,"hbbCandidates");
      const std::vector<VHbbCandidate> & iCand = *vhbbCandHandle.product();
      
      fwlite::Handle< VHbbEventAuxInfo > vhbbAuxHandle; 
      vhbbAuxHandle.getByLabel(ev,"HbbAnalyzerNew");
      const VHbbEventAuxInfo & aux = *vhbbAuxHandle.product();
      
      /*  fwlite::Handle< VHbbEvent > vhbbHandle; 
	  vhbbHandle.getByLabel(ev,"HbbAnalyzerNew");
	  const VHbbEvent iEvent = *vhbbHandle.product();
      */
      trigger.setEvent(&ev);
      VHbbProxy iProxy(0,0, &iCand,&trigger);
      
      //      std::clog << "Filling tree "<< std::endl;
      

      const std::vector<VHbbCandidate> *pcand = iProxy.getVHbbCandidate();
      //////////////
      // TMVA ADD //
      //////////////
      if(pcand->size() > 0 and pcand->at(0).H.jets.size() > 1){
	//variables for tmva
	bbMass = pcand->at(0).H.p4.M();
	bbPt = pcand->at(0).H.p4.Pt();
	btag1 = pcand->at(0).H.jets[0].csv;
	btag2 = pcand->at(0).H.jets[1].csv;
	NaddJet = pcand->at(0).additionalJets.size();
	DeltaRbb = deltaR(pcand->at(0).H.jets[0].p4.Eta(),pcand->at(0).H.jets[0].p4.Phi(),pcand->at(0).H.jets[1].p4.Eta(),pcand->at(0).H.jets[1].p4.Phi());
	helicity = pcand->at(0).H.helicities[0];
	DeltaPhiVH = Geom::deltaPhi(pcand->at(0).H.p4.Phi(),pcand->at(0).V.p4.Phi()) ;
	bPt1 = pcand->at(0).H.jets[0].p4.Pt();
	bPt2 = pcand->at(0).H.jets[1].p4.Pt();
	VMass = pcand->at(0).V.p4.M();
	VPt = pcand->at(0).V.p4.Pt();
	pullAngle = pcand->at(0).H.deltaTheta;
	DeltaEtabb = TMath::Abs( pcand->at(0).H.jets[0].p4.Eta() - pcand->at(0).H.jets[1].p4.Eta() );
	deltaPhipfMETjet1 = Geom::deltaPhi( pcand->at(0).V.mets.at(0).p4.Phi(), pcand->at(0).H.jets[0].p4.Phi() );
	deltaPhipfMETjet2 = Geom::deltaPhi( pcand->at(0).V.mets.at(0).p4.Phi(), pcand->at(0).H.jets[1].p4.Phi() );
	pfMET = pcand->at(0).V.mets.at(0).p4.Pt();
	pfMETsig = pcand->at(0).V.mets.at(0).metSig;
	
	//fill the tree for tmva    
	tree->Fill();
      }

    }
    
    f->Close();
    
  }//loop on files
  
  TreeFile->Write();
  TreeFile->Close();

  return 0;
  
}

