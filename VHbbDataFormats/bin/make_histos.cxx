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
#include "DataFormats/Math/interface/deltaR.h"
/*
std::vector<float> runBDTWe(VHbbProxy* proxy){
  //
  // idea is to use out of the proxy just the candidates
  const std::vector<VHbbCandidate> *pcand = proxy.getVHbbCandidate();
  // check whether the candidate is indeed a We
  if (pcand->candidateType !=VHbbCandidate::Wen){
    std::cout <<" ERROR - Using Wen BDT on a candidate which is not a Wen"<<std::endl;
    abort();
  }
  //
  static TMVA::Reader *readerElec =0;
  static float var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11,var12,var13,var14;
  if (readerElec ==0 ){
    readerElec = new TMVA::Reader( "!Color:!Silent" );

    readerElec->BookMVA("BDT method","/tigress-hsm/dlopes/root528/root/tmva/test/weightsTestJune16Elec/TMVA_test2_BDT.weights.xml");
    readerElec->BookMVA("MLPBNN method","/tigress-hsm/dlopes/root528/root/tmva/test/weightsTestJune16Elec/TMVA_test2_MLPBNN.weights.xml");
    readerElec->BookMVA("BDTG method","/tigress-hsm/dlopes/root528/root/tmva/test/weightsTestJune16Elec/TMVA_test2_BDTG.weights.xml");
    
    readerElec->AddVariable( "bbMass_dij", &var1 );
    readerElec->AddVariable( "bbPt_dij", &var2 );
    readerElec->AddVariable( "btag1_dij", &var3 );
    readerElec->AddVariable( "btag2_dij", &var4 );
    readerElec->AddVariable( "NAddJet1_dij", &var5 );
    readerElec->AddVariable( "DeltaRbb_dij", &var6 );
    readerElec->AddVariable( "hely1_dij", &var7 );
    readerElec->AddVariable( "DeltaPhiWElec_dij", &var8 );
    readerElec->AddVariable( "ptWElec_dij", &var9 );
    readerElec->AddVariable(  "abs(DeltaThetaFlow_dij)", &var10 );
    readerElec->AddSpectator( "ElecFlag_dij", &var11 );
    readerElec->AddSpectator( "NleptonsMu_dij", &var12 );
    readerElec->AddSpectator( "NleptonsElec_dij", &var13 );
    readerElec->AddSpectator( "MET_dij", &var14 );
  }
  //
  // fill variables, please correct
  //  

  var1 = pcand->H.p4.Mass();
  var2 = pcand->H.p4.Pt();
  var3 = pcand->H.jets[0].csv;
  var4 = pcand->H.jets[1].csv;
  var5 = pcand->additionalJets.size();
  var6 = deltaR(pcand->H.jets[0].Eta(),pcand->H.jets[0].Phi(),pcand->H.jets[1].Eta(),pcand->H.jets[1].Phi());
  var7 = pcand->H.helicities[0];
  var8 = Geom::deltaPhi(pcand->H.jets[1].p4.Phi(),V.p4.Phi()) ;
  var9 = TMath::Sqrt((pcand->V.mets[0].p4.Px()+pcand->V.electrons[0].p4.Px())**2
		     +(pcand->V.mets[0].p4.Py()+pcand->V.electrons[0].p4.Py())**2);
  var10 = TMath::Abs( pcand->deltaTheta);
  var11 = pcand->V.electrons[0].id95; //not clear!!
  var12 = 0; // no add lepton allowed
  var13 = 0; // no add lepton allowed
  var14 = pcand->V.mets[0].p4.Pt();
  double mvaoutputBDT = -1;
  double mvaoutputMLPBNN =-1;
  double mvaoutputBDTG = -1;
  
  std::vector<float> results;
  results.push_back(mvaoutputBDT);
  results.push_back(mvaoutputMLPBNN);
  result.push_back(mvaoutputBDTG);
  return results;
}
*/


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

    VHbbProxy iProxy(0,0, &iCand);

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

