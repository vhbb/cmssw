#include <TH1F.h>
#include <TH3F.h>
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "VHbbAnalysis/VHbbDataFormats/interface/HbbCandidateFinderAlgo.h" 
#include "VHbbAnalysis/VHbbDataFormats/src/HbbCandidateFinderAlgo.cc"

#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEventAuxInfo.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/TriggerReader.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/TopMassReco.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/JECFWLite.h"

//for IVF
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include <DataFormats/GeometrySurface/interface/Surface.h>
#include "Math/SMatrix.h"

//for LHE info
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

//Move class definition to Ntupler.h ?
//#include "VHbbAnalysis/VHbbDataFormats/interface/Ntupler.h"

//btagging
#include "VHbbAnalysis/VHbbDataFormats/interface/BTagWeight.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/BTagReshaping.h"
//trigger
#include "VHbbAnalysis/VHbbDataFormats/interface/TriggerWeight.h"


#include "VHbbAnalysis/HbbAnalyzer/interface/HbbAnalyzerAlgo.h"

#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>       // std::vector


#define MAXPDF 100
#define MAXJ 130
#define MAXL 110
#define MAXB 110
#define MAXT 160
#define nMetUnc 24 

int main(int argc, char* argv[]) {
  gROOT->Reset();
  
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");
  AutoLibraryLoader::enable();

  int ievt=0;  
  int totalCount=0;
 

  //FILE LOOP

  //  TFile* inFile = new TFile(inputFile.c_str(), "read");
  std::vector<std::string> inputFiles_;
  inputFiles_.push_back("TT_Tune4C_13TeV-pythia8-tauola_PAT.root");

  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile) {
    std::cout << iFile << std::endl;
    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
    if(inFile==0) { std::cout << "FAILED " << inputFiles_[iFile] << std::endl; continue; }

    // loop the events
      
      fwlite::Event ev(inFile);

      //
      // here I have to insert all the machinery to prepare a VHBBEvent and auxiliaries
      //
      
      
      HbbAnalyzerAlgo theHbbAlgo;
      
      VHbbEventAuxInfo vhbbAux; 

      //EVENT LOOP

      for(ev.toBegin(); !ev.atEnd() ; ++ev, ++ievt)        {
	totalCount++;
	// call the produce from HbbAnalyzerNew
	VHbbEvent hbbInfo;
	VHbbEventAuxInfo auxInfo;
	
	theHbbAlgo.produce(ev,hbbInfo, auxInfo);
	
	const VHbbEventAuxInfo & aux = (auxInfo);
	
	std::cout <<" VHBBEVENT: "<<totalCount<<" " <<hbbInfo.muInfo.size()<<std::endl;

	
	//
	// fake test, I stop here
	//


      }
  }
  exit(0);
}
