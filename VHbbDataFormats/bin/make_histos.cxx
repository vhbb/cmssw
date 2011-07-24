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
//#include "VHbbAnalysis/HbbAnalyzer/interface/HbbAnalyzerNew.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
//#include "VHbbAnalysis/HbbAnalyzer/interface/HbbCandidateFinder.h"
#include "VHbbAnalysis/VHbbDataFormats/src/classes.h"
#include "VHbbAnalysis/VHbbDataFormats/src/CutsAndHistos.cc"
#include <iostream>
#include <fstream>
#include "VHbbAnalysis/VHbbDataFormats/src/Cuts200X.cc"
int main( int argc, char ** argv ){
int event_all=0;
float btag_csv_min = 0.69;
float btag_csv_max = 0.25;
int s=0;
int b=0;
int c=0;
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  gSystem->Load("libDataFormatsFWLite");
  AutoLibraryLoader::enable();

  //  mkdir("./Histograms",755);

  //TFile f("../test/PAT.edm.root");
  //  TFile f("../test/VHbbPAT.edm_6_0_HAB.root");
  //ddTFile &f =* TFile::Open("/scratch/arizzi/VHbbPAT.edm_9_0_tPZ.root");
 TFile *fout = new TFile("VHbbHistos.root","RECREATE");

  std::vector<Histos *> histosForSignalRegions;
  histosForSignalRegions.push_back(new StandardHistos);
  histosForSignalRegions.push_back(new MCHistos);


  std::vector<CutsAndHistos *> allHistos;

  allHistos.push_back(new CutsAndHistos(new SignalRegion,histosForSignalRegions));
  allHistos.push_back(new CutsAndHistos(new VlightRegionHZee,new StandardHistos));
  allHistos.push_back(new CutsAndHistos(new VlightRegionHZmumu,new StandardHistos));
  allHistos.push_back(new CutsAndHistos(new VlightRegionHWen,new StandardHistos));
  allHistos.push_back(new CutsAndHistos(new VlightRegionHWmun,new StandardHistos));
  allHistos.push_back(new CutsAndHistos(new VbbRegionHZee,new StandardHistos));
  allHistos.push_back(new CutsAndHistos(new VbbRegionHZmumu,new StandardHistos));
  allHistos.push_back(new CutsAndHistos(new VbbRegionHWen,new StandardHistos));
  allHistos.push_back(new CutsAndHistos(new VbbRegionHWmun,new StandardHistos));
  allHistos.push_back(new CutsAndHistos(new TTbarRegionHZee,new StandardHistos));
  allHistos.push_back(new CutsAndHistos(new TTbarRegionHZmumu,new StandardHistos));
  allHistos.push_back(new CutsAndHistos(new TTbarRegionHWen,new StandardHistos));
  allHistos.push_back(new CutsAndHistos(new TTbarRegionHWmun,new StandardHistos));

for(size_t a=0;a < allHistos.size(); a++)
{
 allHistos[a]->book(*fout);
}



  std::string fl(argv[1]);
  std::vector<std::string> inputFiles_;
  std::ifstream in(fl.c_str());
  while(!in.eof())
  {
   std::string line;
   in>> line;
   inputFiles_.push_back(line);
  }
   
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
    event_all++;


    fwlite::Handle< std::vector<VHbbCandidate> > vhbbCandHandle; 
    vhbbCandHandle.getByLabel(ev,"hbbCandidates");
    const std::vector<VHbbCandidate> iCand = *vhbbCandHandle.product();

  /*  fwlite::Handle< VHbbEvent > vhbbHandle; 
    vhbbHandle.getByLabel(ev,"HbbAnalyzerNew");
    const VHbbEvent iEvent = *vhbbHandle.product();
*/

    VHbbProxy iProxy(0, &iCand);

    const std::vector<VHbbCandidate> *pcand = iProxy.getVHbbCandidate();

    if(iCand.size()>0)
    {
	for(size_t a=0;a < allHistos.size(); a++)
	{
	 allHistos[a]->process(iProxy, 1);
	}
    } 
  }
  f->Close();
 }//loop on files

  fout->Write();
  fout->Close();
 std::cout << event_all << std::endl;
  return 0;

}

