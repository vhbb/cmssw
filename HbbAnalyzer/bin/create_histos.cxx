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
#include "VHbbAnalysis/HbbAnalyzer/interface/HbbAnalyzerNew.h"
#include "VHbbAnalysis/HbbAnalyzer/interface/HbbCandidateFinder.h"
#include "VHbbAnalysis/HbbAnalyzer/src/classes.h"
#include "VHbbAnalysis/HbbAnalyzer/plugins/CutsAndHistos.cc"

int main( void ){

  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  gSystem->Load("libDataFormatsFWLite");
  AutoLibraryLoader::enable();

  //  mkdir("./Histograms",755);

  //TFile f("../test/PAT.edm.root");
  TFile f("../test/VHbbPAT.edm.root");
  fwlite::Event ev(&f);

  TFile *fout = new TFile("fout.root","RECREATE");

  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetNumberContours(5);

  CutSet noBoostControlRegion;
  noBoostControlRegion.add(new HPt70Cut);
  CutSet boostRegion;
  boostRegion.add(new ZPt50To150Cut);

  std::vector<Histos *> histosForControlRegions;
  histosForControlRegions.push_back(new StandardHistos);

  std::vector<Histos *> histosForBoostedRegions;
  histosForBoostedRegions.push_back(new StandardHistos);

  CutsAndHistos  noBoost(noBoostControlRegion,histosForControlRegions);
  CutsAndHistos  boost(boostRegion,histosForBoostedRegions);

  noBoost.book(*fout);
  boost.book(*fout);

  unsigned int event = 0; 
  for( ev.toBegin(); !ev.atEnd(); ++ev) {
    event++;
    std::cout << " --- Processing event " << event << " --- " << std::endl;

    fwlite::Handle< VHbbEvent > vhbbHandle; 
    vhbbHandle.getByLabel(ev,"HbbAnalyzerNew");
    const VHbbEvent iEvent = *vhbbHandle.product();

    fwlite::Handle< std::vector<VHbbCandidate> > vhbbCandHandle; 
    vhbbCandHandle.getByLabel(ev,"hbbCandidates");
    const std::vector<VHbbCandidate> iCand = *vhbbCandHandle.product();

    VHbbProxy iProxy(&iEvent, &iCand);

    const std::vector<VHbbCandidate> *pcand = iProxy.getVHbbCandidate();

    std::cout << "pt MC higgs from proxy = " << iProxy.getVHbbEvent()->mcH.fourMomentum.Pt() << std::endl;
    std::cout << "pt higgs = " << iEvent.mcH.fourMomentum.Pt() << std::endl;
    if(iCand.size()>0){
      std::cout << "pt Candidate higgs from proxy = " << pcand->at(0).H.fourMomentum.Pt() << std::endl;
      std::cout << "pt Candidate higgs from proxy = " << iProxy.getVHbbCandidate()->at(0).H.fourMomentum.Pt() << std::endl;
      std::cout << "pt Candidate higgs = " << iCand.at(0).H.fourMomentum.Pt() << std::endl;
    }

    noBoost.process(iProxy, 1);
    boost.process(iProxy, 1);

  }

  fout->Write();
  fout->Close();

  return 0;

}

