#ifndef TRIGGERWEIGHT_H
#define TRIGGERWEIGHT_H

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>

class TriggerWeight
{
public:  
  TriggerWeight(const edm::ParameterSet& ana) : tscaleHLTmu(0), tscaleIDmu(0) 
  {
   TFile *hltMuFile = new TFile (ana.getParameter<std::string> ("hltMuFileName").c_str(),"read");
   if(hltMuFile)    tscaleHLTmu = (TTree*) hltMuFile->Get("tree");
   TFile *idMuFile = new TFile (ana.getParameter<std::string> ("idMuFileName").c_str(),"read");
   if(idMuFile)   tscaleIDmu = (TTree*) idMuFile->Get("tree");
   
   if(tscaleHLTmu == 0 || tscaleIDmu == 0) 
    {
      std::cout << "ERROR: cannot load Muon Trigger efficiencies" << std::endl;
    } 


  }
   
 
  float scaleMuIsoHLT(float pt1, float eta1)
  {

    if(tscaleHLTmu) return 1;
    float ptMin,ptMax,etaMin,etaMax,scale,error;
    float s1 = 0;
    int count = 0;
    tscaleHLTmu->SetBranchAddress("ptMin",&ptMin);
    tscaleHLTmu->SetBranchAddress("ptMax",&ptMax);
    tscaleHLTmu->SetBranchAddress("etaMin",&etaMin);
    tscaleHLTmu->SetBranchAddress("etaMax",&etaMax);
    tscaleHLTmu->SetBranchAddress("scale",&scale);
    tscaleHLTmu->SetBranchAddress("error",&error);
    
    for(int jentry = 0; jentry < tscaleHLTmu->GetEntries(); jentry++)
      {
	tscaleHLTmu->GetEntry(jentry);
	if((pt1 > ptMin) && (pt1 < ptMax) && (eta1 > etaMin) && (eta1 < etaMax))
	  {
	    s1 = scale;
	    count++;   
	  }
      }
    
    if(count == 0 || s1 == 0) 
      {
	//caleFile->Close();
	return 1;
      }
    
    
    //aleFile->Close();
    return (s1);
  }
  
  
  
  float scaleMuID(float pt1, float eta1)
  {
    
    if(tscaleIDmu) return 1;

    float ptMin,ptMax,etaMin,etaMax,scale,error;
    float s1 = 0;
    int count = 0;
    tscaleIDmu->SetBranchAddress("ptMin",&ptMin);
    tscaleIDmu->SetBranchAddress("ptMax",&ptMax);
    tscaleIDmu->SetBranchAddress("etaMin",&etaMin);
    tscaleIDmu->SetBranchAddress("etaMax",&etaMax);
    tscaleIDmu->SetBranchAddress("scale",&scale);
    tscaleIDmu->SetBranchAddress("error",&error);
    
    for(int jentry = 0; jentry < tscaleIDmu->GetEntries(); jentry++)
      {
	
	tscaleIDmu->GetEntry(jentry);
	if((pt1 > ptMin) && (pt1 < ptMax) && (eta1 > etaMin) && (eta1 < etaMax))
	  {
	    s1 = scale;
	    count++;   
	  }
      }
    
    if(count == 0 || s1 == 0) 
      {
	//caleFile->Close();
	return 1;
      }
    
    //aleFile->Close();
    return (s1);
    
  }

private:
  TTree * tscaleHLTmu;
  TTree * tscaleIDmu;
};

#endif
