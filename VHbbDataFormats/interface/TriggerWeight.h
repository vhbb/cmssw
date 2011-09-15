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
  TriggerWeight(const edm::ParameterSet& ana) : combiner2Thr(2)
  {
   tscaleHLTmu=openFile(ana,"hltMuFileName");
   tscaleIDmu=openFile(ana,"idMuFileName");
   tscaleHLTele1=openFile(ana,"hltEle1FileName");
   tscaleHLTele2=openFile(ana,"hltEle2FileName");
   tscaleID80Ele=openFile(ana,"idEle80FileName");
   tscaleID95Ele=openFile(ana,"idEle95FileName");
   tscaleHLTeleJet1=openFile(ana,"hltJetEle1FileName");
   tscaleHLTeleJet2=openFile(ana,"hltJetEle2FileName");
   tscaleRecoEle=openFile(ana,"recoEleFileName");
//   tscalePFMHTele=openFile(ana,"hltPFMHTEleFileName");
   tscaleSingleEleMay=openFile(ana,"hltSingleEleMayFileName");
   tscaleSingleEleV4=openFile(ana,"hltSingleEleV4FileName");

   if(tscaleHLTmu == 0 || tscaleIDmu == 0) 
    {
      std::cout << "ERROR: cannot load Muon Trigger efficiencies" << std::endl;
    } 
   


  }
 
  TTree * openFile(const edm::ParameterSet& ana, const char * name)
  {
   TFile *hltMuFile = new TFile (ana.getParameter<std::string> (name).c_str(),"read");
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
  }
  
  
  
  float scaleMuID(float pt1, float eta1)
  {
    return efficiencyFromPtEta(pt1,eta1,tscaleIDmu).first;
  }

double  scaleMetHLT( double met){

    float s1 = 1;
    TF1 * f = new TF1 ("f",TriggerZnunuCurve::trigMet, 0,99999, 0, "triggerZnunuCurve"  );
    
    s1 = f->Eval(met);
    
    return (s1);
    
  }
  

double scaleDoubleEle17Ele8( std::vector<float> pt, std::vector<float> eta )
{
   std::vector< std::vector<float> > allEleWithEffs;
for(unsigned int j=0; j< pt.size(); j++)
 {
  std::vector<float> thisEleEffs;
  thisEleEffs.push_back(efficiencyFromPtEta(pt[j],eta[j],tscaleHLTele1).first);
  thisEleEffs.push_back(efficiencyFromPtEta(pt[j],eta[j],tscaleHLTele2).first);
  allEleWithEffs.push_back(thisEleEffs);
 }

  return   combiner2Thr.weight<Trigger1High2Loose>(allEleWithEffs);

}

double scaleSingleEleMay( float pt, float eta){    return efficiencyFromPtEta(pt,eta,tscaleSingleEleMay).first;}
double scaleSingleEleV4( float pt, float eta){   return efficiencyFromPtEta(pt,eta,tscaleSingleEleV4).first; }
double scaleID80Ele( float pt, float eta) {      return efficiencyFromPtEta(pt,eta,tscaleID80Ele).first; }
double scaleID95Ele( float pt, float eta) {      return efficiencyFromPtEta(pt,eta,tscaleID95Ele).first; }
double scaleRecoEle( float pt, float eta){     return efficiencyFromPtEta(pt,eta,tscaleRecoEle).first; }
double scalePFMHTEle( float MetPFPt){
    double weightPFMHTrigger=0.;

    //FIXME: read from file
    if(MetPFPt>0. && MetPFPt<5.) weightPFMHTrigger=0.3834;
    if(MetPFPt>5. && MetPFPt<10.) weightPFMHTrigger=0.4493;
    if(MetPFPt>10. && MetPFPt<15.) weightPFMHTrigger=0.5676;
    if(MetPFPt>15. && MetPFPt<20.) weightPFMHTrigger=0.6474;
    if(MetPFPt>20. && MetPFPt<25.) weightPFMHTrigger=0.7695;
    if(MetPFPt>25. && MetPFPt<30.) weightPFMHTrigger=0.8936;
    if(MetPFPt>30. && MetPFPt<35.) weightPFMHTrigger=0.9304;
    if(MetPFPt>35. && MetPFPt<40.) weightPFMHTrigger=0.9620;
    if(MetPFPt>40. && MetPFPt<45.) weightPFMHTrigger=0.9894;
    if(MetPFPt>45. && MetPFPt<50.) weightPFMHTrigger=0.9863;
    if(MetPFPt>50. && MetPFPt<60.) weightPFMHTrigger=0.9978;
    if(MetPFPt>60.) weightPFMHTrigger=1;
    return weightPFMHTrigger;
}

double scaleJet30Jet25( std::vector<float> pt, std::vector<float> eta)
{

  std::vector< std::vector<float> > allJetsWithEffs;
for(unsigned int j=0; j< pt.size(); j++)
 {
  std::vector<float> thisJetEffs;
  thisJetEffs.push_back(efficiencyFromPtEta(pt[j],eta[j],tscaleHLTeleJet1).first);
  thisJetEffs.push_back(efficiencyFromPtEta(pt[j],eta[j],tscaleHLTeleJet2).first);
  allJetsWithEffs.push_back(thisJetEffs);
 }

  return   combiner2Thr.weight<Trigger1High2Loose>(allJetsWithEffs);
}



private:
  TTree * tscaleHLTele1;
  TTree * tscaleHLTele2;
  TTree * tscaleHLTeleJet1;
  TTree * tscaleHLTeleJet2;
  TTree * tscaleID80Ele;
  TTree * tscaleID95Ele;
  TTree * tscaleRecoEle;
//  TTree * tscalePFMHTele;
  TTree * tscaleSingleEleMay;
  TTree * tscaleSingleEleV4;

  TTree * tscaleHLTmu;
  TTree * tscaleIDmu;
  MultiThresholdEfficiency combiner2Thr; 
};

#endif
