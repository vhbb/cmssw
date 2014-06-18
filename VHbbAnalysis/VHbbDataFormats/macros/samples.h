#include <string>
#include <vector>
#include <TFile.h>
#include <TH1F.h>
#include <iostream>


struct Sample {
  Sample(float xs,std::string n,std::string f, int c, bool isdata,float datalumi=-1.)
   : xsec(xs),luminosity(datalumi),name(n),filename(f),color(c),data(isdata),f(0),nevents(-1) {}

  float lumi() {  if(data) return luminosity; else return numberOfEvents()/xsec; }
  float scale(float l) { return l/lumi();}
  TFile * file() { if(f) return f; else return f=TFile::Open(filename.c_str());}
  float numberOfEvents() 
  {
   if(nevents !=-1) return nevents;
   else
   {
    return ((TH1F*)file()->Get("NoCut/CountNoCut"))->GetEntries();
   }
 }  

   void dump(float l)
   {
    std::cout << name <<  "\t& " << xsec << "\t& " <<  lumi()/1000 << "/fb \t& " << scale(l) << std::endl;
   }

  float nevents;
  float xsec;
  float luminosity;
  std::string name;
  std::string filename;
  int color;
  bool data;
  TFile * f;
}; 

/*
DoubleElectron_HBB_EDMNtupleV1_ProcV2_CandV1_may_histos.root
DoubleElectron_HBB_EDMNtupleV1_ProcV2_CandV1_prompt_histos.root
DoubleMu_HBB_EDMNtupleV1_ProcV2_CandV1_prompt_histos.root
DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root
METBTag_HBB_EDMNtupleV1_ProcV2_CandV1_may_histos.root
MET_HBB_EDMNtupleV1_ProcV2_CandV1_prompt_histos.root
SingleMu_HBB_EDMNtupleV1_ProcV2_CandV1_may_histos.root
SingleMu_HBB_EDMNtupleV1_ProcV2_CandV1_prompt_histos.root
Tbar_TuneZ2_s-channel_7TeV-powheg-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root
Tbar_TuneZ2_t-channel_7TeV-powheg-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root
Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root
Tbar_TuneZ2_tW-channel-DS_7TeV-powheg-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root
TTJets_TuneZ2_7TeV-madgraph-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root
T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root
WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root
WW_TuneZ2_7TeV_pythia6_tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root
WZ_TuneZ2_7TeV_pythia6_tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root
ZZ_TuneZ2_7TeV_pythia6_tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root

*/

std::vector<Sample> samples()
{
 std::vector<Sample> s;

/*
235224077.721 Multi_DT/DoubleElectron_May10Rereco/res/lumiSummary.json
814780712.247 Multi_DT/DoubleElectron_PromptReco/res/lumiSummary.json
833435176.929 Multi_DT/DoubleMu_PromptReco/res/lumiSummary.json
235558534.699 Multi_DT/METBTag_May10Rereco/res/lumiSummary.json
787760254.134 Multi_DT/MET_PromptReco/res/lumiSummary.json
231959363.356 Multi_DT/SingleElectron_May10Rereco/res/lumiSummary.json
811788598.096 Multi_DT/SingleElectron_PromptReco/res/lumiSummary.json
219100833.649 Multi_DT/SingleMu_May10Rereco/res/lumiSummary.json
788877165.417 Multi_DT/SingleMu_PromptReco/res/lumiSummary.json


*/

// s.push_back(Sample(1000,"data","SingleMu_HBB_EDMNtupleV1_ProcV2_CandV1_may_histos.root",0 , true,219));
// s.push_back(Sample(1000,"data","SingleMu_HBB_EDMNtupleV1_ProcV2_CandV1_prompt_histos.root",0 , true,789));
 s.push_back(Sample(1000,"data","SingleMu_HBB_EDMNtupleV1_ProcV2_CandV1_merge_histos.root",0 , true,788.8+219));

// s.push_back(Sample(1000,"data","DoubleElectron_HBB_EDMNtupleV1_ProcV2_CandV1_may_histos.root", 0, true ,235.22));
// s.push_back(Sample(1000,"data","DoubleElectron_HBB_EDMNtupleV1_ProcV2_CandV1_prompt_histos.root",1 , true,814.5 ));
// s.push_back(Sample(1000,"data","DoubleElectron_HBB_EDMNtupleV1_ProcV2_CandV1_merge_histos.root",1 , true,235.22 + 814.8));

// s.push_back(Sample(1000,"data","DoubleMu_HBB_EDMNtupleV1_ProcV2_CandV1_prompt_histos.root", 1, true, 500.159));

// s.push_back(Sample(1000,"data","METBTag_HBB_EDMNtupleV1_ProcV2_CandV1_may_histos.root", 1, true,235));
// s.push_back(Sample(1000,"data","MET_HBB_EDMNtupleV1_ProcV2_CandV1_prompt_histos.root", 1, true,784.12));


  s.push_back(Sample(165,"TTbar","TTJets_TuneZ2_7TeV-madgraph-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root", kBlue , false ));

/*
histMassWjetLF->SetFillColor(kSpring-6);
histMassWjetHF->SetFillColor(kSpring);
histMassTTbar->SetFillColor(kBlue);
histMassQCD->SetFillColor(kMagenta);
histMassWW->SetFillColor(kOrange+10);
histMassWZ->SetFillColor(kOrange+10);
histMassSingleToptW->SetFillColor(kTeal);
*/

int stcolor=kTeal;

// s.push_back(Sample(1.44,"Single Top","Tbar_TuneZ2_s-channel_7TeV-powheg-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root", stcolor, false ));
 s.push_back(Sample(22.65,"Single Top","Tbar_TuneZ2_t-channel_7TeV-powheg-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root", stcolor, false ));
 s.push_back(Sample(7.87,"Single Top","Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root", stcolor, false));
//s.push_back(Sample(7.87,"Single Top","Tbar_TuneZ2_tW-channel-DS_7TeV-powheg-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root", stcolor, false));
 s.push_back(Sample(7.87,"Single Top","T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root", stcolor, false));

 float wxsec= 31314.;
 float wxsec100= 31314./27770.*194.6;
//TOT: 18904365 b: 363441 c: 6264682 l: 12276242
 float t=18904365;
 float b=363441;
 float c=6264682;
 float l=12276242;
 
 s.push_back(Sample(wxsec*b/t,"Wb","WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histosB.root", kSpring, false ));
 s.push_back(Sample(wxsec*c/t,"Wc","WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histosC.root", kSpring-4, false ));
 s.push_back(Sample(wxsec*l/t,"Wl","WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histosL.root", kSpring-6, false ));

 float zxsecMG=2475;
 s.push_back(Sample(3048*0.0441,"Zb","DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histosB.root",9 ,false ));
 s.push_back(Sample(3048*0.244,"Zc","DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histosC.root",11 ,false ));
 s.push_back(Sample(3048*0.711,"Zl","DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histosL.root",12 ,false ));

int VVcolor=5;
 s.push_back(Sample(42.9,"VV","WW_TuneZ2_7TeV_pythia6_tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root",kOrange+10 , false ));
 s.push_back(Sample(18.3,"VV","WZ_TuneZ2_7TeV_pythia6_tauola_HBB_EDMNtupleV1_ProcV2_histos.root",kOrange+10 , false ));
 s.push_back(Sample(5.9,"VV","ZZ_TuneZ2_7TeV_pythia6_tauola_HBB_EDMNtupleV1_ProcV2_CandV1_histos.root",kOrange+10 , false ));

 return s;
}
