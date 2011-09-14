#ifndef TRIGGERWEIGHT_H
#define TRIGGERWEIGHT_H

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/TriggerZnunuCurve.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/MultiThresholdEfficiency.h"
#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>

class TriggerWeight
{
public:  
  TriggerWeight(const edm::ParameterSet& ana) : tscaleHLTele1(0),tscaleHLTele2(0),tscaleHLTeleJet1(0),tscaleHLTeleJet2(0),tscaleIDele(0), tscaleHLTmu(0), tscaleIDmu(0),combiner2Thr(2)
  {
   tscaleHLTmu=openFile(ana,"hltMuFileName");
   tscaleIDmu=openFile(ana,"idMuFileName");
   tscaleHLTele1=openFile(ana,"hltEle1FileName");
   tscaleHLTele2=openFile(ana,"hltEle2FileName");
   tscaleIDele=openFile(ana,"idEleFileName");
   tscaleHLTeleJet1=openFile(ana,"hltJetEle1FileName");
   tscaleHLTeleJet2=openFile(ana,"hltJetEle2FileName");

   if(tscaleHLTmu == 0 || tscaleIDmu == 0) 
    {
      std::cout << "ERROR: cannot load Muon Trigger efficiencies" << std::endl;
    } 
   


  }
 
  TTree * openFile(const edm::ParameterSet& ana, const char * name)
  {
   TFile *hltMuFile = new TFile (ana.getParameter<std::string> ("hltMuFileName").c_str(),"read");
   if(hltMuFile)   return (TTree*) hltMuFile->Get("tree");
   else return 0;
  }

  std::pair<float,float> efficiencyFromPtEta(float pt1, float eta1, TTree *t)
  {
     float s1 = 1.,err=1.;
     std::pair<float,float> r(s1,err);
     if(!t) return r;
     float ptMin,ptMax,etaMin,etaMax,scale,error;
     int count = 0;
     t->SetBranchAddress("ptMin",&ptMin);
     t->SetBranchAddress("ptMax",&ptMax);
     t->SetBranchAddress("etaMin",&etaMin);
     t->SetBranchAddress("etaMax",&etaMax);
     t->SetBranchAddress("scale",&scale);
     t->SetBranchAddress("error",&error);

    for(int jentry = 0; jentry < t->GetEntries(); jentry++)
      {
        t->GetEntry(jentry);
        if((pt1 > ptMin) && (pt1 < ptMax) && (eta1 > etaMin) && (eta1 < etaMax))
          {
            s1 = scale;
            err=error;
            count++;
          }
      }

    if(count == 0 || s1 == 0)
      {
        return r;
      }

     r.first=s1;
     r.second = err;
    return (r);
  }
 
  float scaleMuIsoHLT(float pt1, float eta1)
  {
    return efficiencyFromPtEta(pt1,eta1,tscaleHLTmu).first;
/*    if(!tscaleHLTmu) return 1;
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
    return (s1);*/
  }
  
  
  
  float scaleMuID(float pt1, float eta1)
  {
    return efficiencyFromPtEta(pt1,eta1,tscaleIDmu).first;

    // changed to !tscale...
    if(!tscaleIDmu) return 1;

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

double  scaleMetHLT( double met){

    float s1 = 1;
    TF1 * f = new TF1 ("f",TriggerZnunuCurve::trigMet, 0,99999, 0, "triggerZnunuCurve"  );
    
    s1 = f->Eval(met);
    
    return (s1);
    
  }
  

double scaleJet30Jet25( std::vector<float> eta, std::vector<float> pt)
{
return 1;
//     if(!tscaleEleJet1) return 1;
//     if(!tscaleEleJet2) return 1;
     vector<vector<float > > bah;
     combiner2Thr.weight<Trigger1High2Loose>(bah);
}



private:
  TTree * tscaleHLTele1;
  TTree * tscaleHLTele2;
  TTree * tscaleHLTeleJet1;
  TTree * tscaleHLTeleJet2;
  TTree * tscaleIDele;

  TTree * tscaleHLTmu;
  TTree * tscaleIDmu;
  MultiThresholdEfficiency combiner2Thr; 
};

#endif
