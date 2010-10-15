#include <vector>
#include "TFile.h"
#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <sys/stat.h> 
#include "TKey.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TEventList.h"
#include "bbbarcorr.h"
#define MAXJETS 100

////////////////////////
// USAGE
//
// .L bbbarcorr_new.C+
// analyze(SET, HLT), e.g. analyze("Pythia6","HLT_Jet30U")
//
// to analyze a file in particular: analyze("Pythia6","HLT_Jet30U",2)
//
// to sum histos e.g.: sumHisto("Pythia6","HLT_Jet30U") 
// to sum data histos with prescales: sumHistoPrescaled("Pythia6","HLT_Jet30U")
////////////////////////

///Configuration
//general selection cuts
float BJET_PT = 30;
float BJET_ETA = 2;
float JET1_ETA = 3;
float JET_PT = BJET_PT;
float JET_ETA = 2;

//btag cuts
string BALGO = "ssvhp";
string BALGOWP = "ssvhpt";
float BDISCR = 2;

//normalization in pb
double LUMI = 1;
bool applyPrescale = false;

//label for the outfile
string LABEL = "";

//where to find the ttrees
string DIR="/scratch/leo/";

//file containing btag purity fit
string jetCorrFile = "histo-Pythia6-flavourComp.root";
TFile *f_jetPurity;
TF1 *f_jet_pur_pt;

//file containing the IVF corrections:
string ivfCorrFile = "/scratch/leo/efficCurves.root";

///variables
map <string, TH1F*> histos;
map <string, TH2F*> histos2D;

std::vector<TFile*> files;

std::vector<double> xsec; //pb
//std::vector<double> nEvents; 
std::vector<double> lumi; 
std::vector<double> lumiFactor; 
std::vector<int> upperpthatcut; 
double totLumiData; 

//trigger windows
float jet1_minPt, jet1_maxPt;
string jet1_minPt_lab, jet1_maxPt_lab;

//pthat window
std::vector<double> maxPthat;
std::vector<double> minPthat;

//variables
map<string,float> varFloat;
map<string,double> varDouble;
map<string,float*> varFloatArr;
map<string, int> varInt;
map<string, int*> varIntArr;
TH1F *hcorrIVF; 
double xMax_IVFCorr;
double xMin_IVFCorr;
double binW_IVFCorr;

//categories labels
string lab_fl[] = {"2b","1b1c","1b1l","1b1o",
	       "1c1c","1c1l","1c1o",
	       "1l1l","1l1o",
	       "2o"};
string lab_eff[] = {"efficnum","efficden","puritynum","purityden"};
string lab_cat[] = {"cat2bmat","cat2bnotmat","cat2bnotkin","cat0b","catmanyb","catnocat"};
string lab_proc[] = {"fcr","fex","gsp","otherproc"};


//functions
void createCommonHistos(string SET);
vector<TFile*> getFiles(string SET, string HLT);
void setTrigger(string HLT);
void setBranches(TTree* myTree);
void createJetHistos(string SET);
void createIVFHistos(string SET);
int analyzeJetEvent(string SET, double W);
int analyzeIVFEvent(event thisEvent, string SET, double W);
void sumHisto(string SET, string HLT);
void sumHistoPrescaled(string SET, string HLT);

void corrFunctionIVF(string strigger);
float delta_phi(float phi1, float phi2);
float delta_R(float phi1, float phi2, float eta1, float eta2);

/// the main analisys routine
void analyze(string SET="Data", string HLT="HLT_Jet15U", int SET_n=-1){
  files = getFiles(SET,HLT);
  setTrigger(HLT);

  //selecting files
  int init_i=0;
  int end_i = files.size();
  if(SET_n!=-1){
    init_i=SET_n;
    end_i = init_i+1;
  }

  //files loop
  for(int f_i=init_i; f_i< end_i; f_i++){
    //the output file
    string fname = files[f_i]->GetName();
    cout << "Analyzing " << fname << endl;
    fname = fname.replace( fname.begin(), fname.begin()+fname.rfind("/")+1,"" );
    struct stat stDirInfo;
    if(stat(SET.c_str(),&stDirInfo) != 0) system(("mkdir "+SET).c_str());
    TFile *f_out = TFile::Open((SET+"/histo-"+HLT+"-"+fname ).c_str(),"RECREATE");
    cout << "Output file: " << SET+"/histo-"+HLT+"-"+fname << endl;
    //f_out->cd();
    //f_out->mkdir(HLT.c_str());
    //f_out->cd(HLT.c_str());

    //booking histos
    createCommonHistos(SET);
    createJetHistos(SET);
    createIVFHistos(SET);

    //Selecting a reduced tree
    TTree *myTree = (TTree*)files[f_i]->Get("bcanalyzer/tEvents");
    int events_n = myTree->GetEntries();
    string mySel = HLT+"==1 &&";
    if(SET!="Data") mySel=""; //temporary fix for HLT collections issue
    mySel += " PFJet_JetIDLoose[0]==1 && PFJet_pt[0]>= "+jet1_minPt_lab;
    myTree->Draw(">>selList", mySel.c_str()); 
    TEventList *mySelList = (TEventList*)gDirectory->Get("selList"); 
    myTree->SetEventList(mySelList);

    //set weights
    double W = 1;
    if(SET!="Data") W  = (xsec[f_i]*LUMI)/((double) events_n);
    else if (SET=="Data" && applyPrescale) W = lumiFactor[f_i];

    //setting branches
    setBranches(myTree);
    event thisEv; 
    myTree->SetBranchAddress("BEvents",&thisEv);

    //IVF correction
    struct stat stFileInfo; 
    if(stat(ivfCorrFile.c_str(),&stFileInfo) == 0){
      if(HLT=="HLT_Jet30U") corrFunctionIVF("jet30");
      else if(HLT=="HLT_Jet50U") corrFunctionIVF("jet50");
      else if(HLT=="HLT_Jet15U") corrFunctionIVF("jet15");
    }else{
      cout << "Warning: IVF correction file not found" << endl;
    }
    //BTag eff correction
    if(stat(jetCorrFile.c_str(),&stFileInfo) == 0){
      f_jetPurity = TFile::Open("histo-Pythia6-flavourComp.root");
      f_jetPurity->GetObject("f_purity_pt", f_jet_pur_pt);
    }else{
      cout << "Warning: Jet correction file not found" << endl;
    }

    //event loop
    int counter = 0;
    while(mySelList->GetEntry(counter++) !=-1){
      int entry = mySelList->GetEntry(counter-1);
      myTree->GetEntry(entry);
      bool keepEvent=true;

      if(SET!="Data"){
	if( maxPthat[f_i] == -1 && minPthat[f_i]==-1 ) keepEvent = true;
        else if( varDouble["pthat"]>maxPthat[f_i] || varDouble["pthat"] <=minPthat[f_i] ) keepEvent = false;
      }
      if(!keepEvent) continue;

      //checking for leading jet
      if( varIntArr["Jet_IDLoose"][0]!=1) continue;
      if( varFloatArr["Jet_pt"][0]<jet1_minPt || varFloatArr["Jet_pt"][0] > jet1_maxPt ) continue;
      if( fabs(  varFloatArr["Jet_eta"][0]) > JET1_ETA) continue;

      //event analyzers
      int ret_code=0;
      ret_code = analyzeJetEvent(SET, W);
      ret_code =  analyzeIVFEvent(thisEv, SET, W);

    }//end event loop

    f_out->Write();
    f_out->Close();
    delete myTree;
    histos.clear();
    histos2D.clear();
  }//end file loop

  
}





//actual analyzer for jet-based analysis
int analyzeJetEvent(string SET, double W){
  //  cout << "JET" << endl;
  int bJets_i[100];
  int bJetsMC_i[100];
  int n_bJets = 0;
  int n_bJetsMC = 0;


  //jet loop
  for(int j=0; j< varInt["Jet_n"];j++){
    if( varIntArr["Jet_IDLoose"][j]!=1) continue;
    if( varFloatArr["Jet_pt"][j] < BJET_PT ) continue;
    if( fabs(varFloatArr["Jet_eta"][j]) > BJET_ETA) continue;

    histos["Jet_pt"]->Fill( varFloatArr["Jet_pt"][j],W);
    histos["Jet_eta"]->Fill( varFloatArr["Jet_eta"][j],W);
    if( abs(varIntArr["Jet_flavour"][j])==1 || 
	abs(varIntArr["Jet_flavour"][j])==2 ||
	abs(varIntArr["Jet_flavour"][j])==3 ) { histos["Jet_lJetsMC_pt"]->Fill(varFloatArr["Jet_pt"][j] ,W); histos["Jet_lJetsMC_eta"]->Fill(varFloatArr["Jet_eta"][j] ,W);}
    else if( abs(varIntArr["Jet_flavour"][j])==4) {histos["Jet_cJetsMC_pt"]->Fill(varFloatArr["Jet_pt"][j] ,W);  histos["Jet_cJetsMC_eta"]->Fill(varFloatArr["Jet_eta"][j] ,W);}
    else if( abs(varIntArr["Jet_flavour"][j])==5) { bJetsMC_i[n_bJetsMC++] = j;  histos["Jet_bJetsMC_pt"]->Fill(varFloatArr["Jet_pt"][j] ,W);  histos["Jet_bJetsMC_eta"]->Fill(varFloatArr["Jet_eta"][j] ,W);}


    histos["Jet_discr"]->Fill( varFloatArr["Jet_bDiscr"][j], W);
    if(  varFloatArr["Jet_bDiscr"][j] > BDISCR){
      bJets_i[n_bJets++] = j;
      histos["Jet_bJets_pt"]->Fill(varFloatArr["Jet_pt"][j] ,W);
      histos["Jet_bJets_eta"]->Fill(varFloatArr["Jet_eta"][j] ,W);
    }
  }//end jet loop

  histos["Jet_bJets_n"]->Fill( n_bJets ,W);

  //True MC bjet content
  if(SET!="Data" && n_bJetsMC==2){
    int bj1 = bJetsMC_i[0]; int bj2 = bJetsMC_i[1];
    float dR = delta_R( varFloatArr["Jet_phi"][ bj1 ],varFloatArr["Jet_phi"][ bj2 ],  varFloatArr["Jet_eta"][ bj1 ],varFloatArr["Jet_eta"][ bj2 ]  );
    float dPhi = delta_phi(  varFloatArr["Jet_phi"][ bj1 ],varFloatArr["Jet_phi"][ bj2 ] );
    float dEta = fabs( varFloatArr["Jet_eta"][ bj1 ] - varFloatArr["Jet_eta"][ bj2 ]);

    histos["Jet_bJetsMC_dR"]->Fill(dR,W);
    histos["Jet_bJetsMC_dEta"]->Fill(dEta,W);
    histos["Jet_bJetsMC_dPhi"]->Fill(dPhi,W);

  }
  
  //2 bjets part
  if(n_bJets!=2) return 1;
  int bj1 = bJets_i[0]; int bj2 = bJets_i[1];
  float dR = delta_R( varFloatArr["Jet_phi"][ bj1 ],varFloatArr["Jet_phi"][ bj2 ],  varFloatArr["Jet_eta"][ bj1 ],varFloatArr["Jet_eta"][ bj2 ]  );
  histos["Jet_dR_2b"]->Fill(dR,W);
  float dPhi = delta_phi(  varFloatArr["Jet_phi"][ bj1 ],varFloatArr["Jet_phi"][ bj2 ] );
  histos["Jet_dPhi_2b"]->Fill(dPhi,W);
  float dEta = fabs( varFloatArr["Jet_eta"][ bj1 ] - varFloatArr["Jet_eta"][ bj2 ]);
  histos["Jet_dEta_2b"]->Fill( dEta ,W);
  
  histos["Jet_pt_2b"]->Fill(  varFloatArr["Jet_pt"][bj1],W);
  histos["Jet_pt_2b"]->Fill(  varFloatArr["Jet_pt"][bj2],W);
  histos["Jet_eta_2b"]->Fill(  varFloatArr["Jet_eta"][bj1],W);
  histos["Jet_eta_2b"]->Fill(  varFloatArr["Jet_eta"][bj2],W);
  
  histos["Jet_ptAsymm_2b"]->Fill( fabs(varFloatArr["Jet_pt"][bj1]-varFloatArr["Jet_pt"][bj2])/( varFloatArr["Jet_pt"][bj1]+varFloatArr["Jet_pt"][bj2]),W);

  float newW = varFloatArr["Jet_bEff"][0]*varFloatArr["Jet_bEff"][1];
  float meanPur = f_jet_pur_pt->Eval( varFloatArr["Jet_pt"][bj1])*f_jet_pur_pt->Eval( varFloatArr["Jet_pt"][bj2])  ;
  newW = meanPur/newW;
  histos["Jet_dR_2b_CORR"]->Fill(dR,newW*W);
  histos["Jet_dPhi_2b_CORR"]->Fill(dPhi,newW*W);
  histos["Jet_dEta_2b_CORR"]->Fill(dEta,newW*W);


  //invariant mass
  float invM = (varFloatArr["Jet_E"][bj1]+varFloatArr["Jet_E"][bj2])*(varFloatArr["Jet_E"][bj1]+varFloatArr["Jet_E"][bj2]);
  invM -= (varFloatArr["Jet_px"][bj1]+varFloatArr["Jet_px"][bj2])*(varFloatArr["Jet_px"][bj1]+varFloatArr["Jet_px"][bj2]);
  invM -= (varFloatArr["Jet_py"][bj1]+varFloatArr["Jet_py"][bj2])*(varFloatArr["Jet_py"][bj1]+varFloatArr["Jet_py"][bj2]);
  invM -= (varFloatArr["Jet_pz"][bj1]+varFloatArr["Jet_pz"][bj2])*(varFloatArr["Jet_pz"][bj1]+varFloatArr["Jet_pz"][bj2]);
  invM = sqrt(invM);
  histos["Jet_2bJets_Mass"]->Fill(invM,W);
  
  //flavour dependent an
  if(SET!="Data"){
    histos2D["Jet_dEtadPhi_2b"]->Fill( dEta,dPhi, W);
    histos2D["Jet_pt1pt2_2b"]->Fill(varFloatArr["Jet_pt"][bj1], varFloatArr["Jet_pt"][bj2],W);

    int f1 = abs( varIntArr["Jet_flavour"][bj1]); int f2 = abs( varIntArr["Jet_flavour"][bj2]);
    if( f2==5 ) {int t_f2=f1; f1=f2; f2=t_f2;} //f1 is always a b, if possible
    string diFl = "";
    if(f1==5){
      if(f2==5) diFl ="2b";
      else if(f2==4) diFl ="1b1c";
      else if(f2==1 || f2==2 || f2==3) diFl ="1b1l";
      else diFl ="1b1o";
    }
    else if(f1==4){
      if(f2==4) diFl="1c1c";
      else if(f2==1 || f2==2 || f2==3)  diFl="1c1l";
      else diFl="1c1o";
    }
    else if(f1==1 || f1==2 || f1==3 ){
      if(f2==1 || f2==2 || f2==3)  diFl="1l1l";
      else diFl="1l1o";
    }
    else diFl ="2o";
    
    histos["Jet_pt_2b_"+diFl]->Fill( varFloatArr["Jet_pt"][bj1],W);
    histos["Jet_pt_2b_"+diFl]->Fill( varFloatArr["Jet_pt"][bj2],W);
    histos["Jet_dR_2b_"+diFl]->Fill( dR,W);
    histos["Jet_dPhi_2b_"+diFl]->Fill( dPhi,W );
    histos["Jet_dEta_2b_"+diFl]->Fill( dEta,W );
    histos["Jet_2bJets_Mass_"+diFl]->Fill(invM,W);
    histos2D["Jet_dEtadPhi_2b_"+diFl]->Fill( dEta,dPhi, W);
    histos2D["Jet_pt1pt2_2b_"+diFl]->Fill(varFloatArr["Jet_pt"][bj1], varFloatArr["Jet_pt"][bj2],W);


  }
  
  return 0;
}





//IVF analyzer
int analyzeIVFEvent(event thisEv, string SET, double W){
  //selecting events with only 2 b ??? Lukas, is this necessary?
  //yes, we want to have events with exactly two B since otherwise 
  //it is not clear which are the B's produced in same process
  if( thisEv.nV!=2 && thisEv.nB!=2) return 1;
  histos["bVert_N"]->Fill(thisEv.nV, W);

  int flavCat = -1, matCat=-1; 
  if(thisEv.dRvv>-700){
    //mass sum cut
    if((thisEv.massV1 + thisEv.massV2) > 4.5){
      histos["Vert_dR_2b"]->Fill(thisEv.dRvv,W);
      histos["Vert_dPhi_2b"]->Fill(fabs(thisEv.dPhivv),W);
      histos["Vert_dEta_2b"]->Fill(fabs(thisEv.dEtavv),W);
      histos["Vert_goodEv"]->Fill(thisEv.eventNo,W);
      histos["Vert_dRingoodEv"]->Fill(thisEv.eventNo, thisEv.dRvv*W); //for Lukas: is this a TH1F or a TH2F?
      histos["Vert_pt_selvert"]->Fill(thisEv.ptV1,W); 
      histos["Vert_pt_selvert"]->Fill(thisEv.ptV2,W); 
      
      //corrected plots
      if(hcorrIVF!=0){
	double weight = hcorrIVF->GetBinContent((int)((thisEv.dRvv-xMin_IVFCorr)/binW_IVFCorr+1.)); 
	histos["Vert_dR_2b_CORR"]->Fill(thisEv.dRvv,W*weight);
	histos["Vert_dPhi_2b_CORR"]->Fill(fabs(thisEv.dPhivv),W*weight);
	histos["Vert_dEta_2b_CORR"]->Fill(fabs(thisEv.dEtavv),W*weight);
      }
      //else std::cout << "NO FILE FOR CORRECTED IVF PLOTS LOADED\n";
      
      if(SET!="Data"){
	string diFl="";
	if(thisEv.flavors==55){ diFl = "2b"; flavCat=0; }
	else if(thisEv.flavors==54){ diFl = "1b1c"; flavCat=0; }
	else if(thisEv.flavors==51){ diFl = "1b1l"; flavCat=1; }
	else if(thisEv.flavors==44){ diFl = "1c1c"; flavCat=2; }
	else if(thisEv.flavors==41){ diFl = "1c1l"; flavCat=3; }
	else if(thisEv.flavors==11){ diFl = "1l1l"; flavCat=4; }
	else cout << "No flavour!" << endl;
	
	histos["Vert_dR_2b_"+diFl]->Fill(thisEv.dRvv,W);
	histos["Vert_dPhi_2b_"+diFl]->Fill(fabs(thisEv.dPhivv),W);
	histos["Vert_dEta_2b_"+diFl]->Fill(fabs(thisEv.dEtavv),W);
      
	//categories
	//2B
	string myCat="catnotcat"; 
	if(thisEv.nB==2){
	  //B kinematics
	  if(thisEv.ptB1>15 && thisEv.ptB2>15 && fabs(thisEv.etaB1)<2.4 && fabs(thisEv.etaB2)<2.4){
	    if(thisEv.nMat==2){ myCat="cat2bmat"; matCat=0;}
	    else { myCat="cat2bnotmat"; matCat=1;}
	  }
	  else{myCat="cat2bnotkin"; matCat=2;}
	}
	else if(thisEv.nB==0){ myCat="cat0b"; matCat=3;}
	else if(thisEv.nB>2){ myCat="catmanyb"; matCat=4;}

	histos["Vert_dR_2b_"+myCat]->Fill(thisEv.dRvv,W);
	histos["Vert_dPhi_2b_"+myCat]->Fill(thisEv.dPhivv,W);
	histos["Vert_dEta_2b_"+myCat]->Fill(thisEv.dEtavv,W);

	//process
	string myProc="";
	if(thisEv.process==1) myProc="fcr";
	else if(thisEv.process==2) myProc="fex";
	else if(thisEv.process==3) myProc="gsp";
	else myProc="otherproc";
	histos["Vert_dR_2b_"+myProc]->Fill(thisEv.dRvv,W);
	histos["Vert_dPhi_2b_"+myProc]->Fill(thisEv.dPhivv,W);
	histos["Vert_dEta_2b_"+myProc]->Fill(thisEv.dEtavv,W);

      
      //V kinematics
      if(thisEv.ptV1>8 && thisEv.ptV2>8){
	histos["Vert_dR_purityden"]->Fill(thisEv.dRvv,W);
	histos["Vert_dPhi_purityden"]->Fill(fabs(thisEv.dPhivv),W);
	histos["Vert_dEta_purityden"]->Fill(fabs(thisEv.dEtavv),W);
	
	//2B
	if(thisEv.nB==2){
	  //B kinematics
	  if(thisEv.ptB1>15 && thisEv.ptB2>15 && fabs(thisEv.etaB1)<2.4 && fabs(thisEv.etaB2)<2.4){
	    histos["Vert_dR_efficnum"]->Fill(thisEv.dRbb,W);
	    histos["Vert_dR_puritynum"]->Fill(thisEv.dRbb,W);
	    histos["Vert_dPhi_efficnum"]->Fill(fabs(thisEv.dPhibb),W);
	    histos["Vert_dPhi_puritynum"]->Fill(fabs(thisEv.dPhibb),W);
	    histos["Vert_dEta_efficnum"]->Fill(fabs(thisEv.dEtabb),W);
	    histos["Vert_dEta_puritynum"]->Fill(fabs(thisEv.dEtabb),W);
	  }
	}
      }

      if(flavCat!=-1 && matCat!=-1)
	histos2D["Vert_flavCat_vs_matCat"]->Fill(flavCat,matCat,W);
      }
    }
  }
  
  //2B, Bkinematics
  if(thisEv.nB==2 && thisEv.ptB1>15 && thisEv.ptB2>15 && fabs(thisEv.etaB1)<2.4 && fabs(thisEv.etaB2)<2.4) {
    histos["Vert_dR_efficden"]->Fill(thisEv.dRbb,W);
    histos["Vert_dPhi_efficden"]->Fill(fabs(thisEv.dPhibb),W);
    histos["Vert_dEta_efficden"]->Fill(fabs(thisEv.dEtabb),W);
  }
  return 0;
}


//a simple call to hadd
void sumHisto(string SET, string HLT){
  files = getFiles(SET,HLT);
  
  //selecting files
  int init_i=0;
  int end_i = files.size();

  //files loop
  string summedFileName = "histo-"+SET+"-"+HLT+"-total.root";
  string command = "hadd -f "+summedFileName+" ";
  for(int f_i=init_i; f_i< end_i; f_i++){
    //the output file
    string fname = files[f_i]->GetName();
    fname = fname.replace( fname.begin(), fname.begin()+fname.rfind("/")+1,"" );
    fname=SET+"/histo-"+HLT+"-"+fname;
    command += fname+" ";
  }

  system( command.c_str() );
}

void sumHistoPrescaled(string SET, string HLT){
  files = getFiles(SET,HLT);

  //selecting files
  int init_i=0;
  int end_i = files.size();

  string summedFileName = "histo-"+SET+"-"+HLT+"-total-wPrescales.root";

  string histoNames[100];
  TH1F * resultHistos[100];
  TH1F * histos[100][100];
  TF1 *weights[10];
  int TOTHIST=0;
  TFile *hfiles[100];

  for(int f_i=init_i; f_i< end_i; f_i++){
    TOTHIST=0;
    string fname = files[f_i]->GetName();
    fname = fname.replace( fname.begin(), fname.begin()+fname.rfind("/")+1,"" );
    fname=SET+"/histo-"+HLT+"-"+fname;
    hfiles[f_i] = TFile::Open(fname.c_str());
    char*ws = new char[30]; sprintf(ws, "%.4f", lumiFactor[f_i] );    
    char*is = new char[30]; sprintf(is, "w%d", f_i );
    //cout << fname << " " << ws << " " << is << endl;
    weights[f_i] = new TF1(is,ws,-10000,10000); 

    TIter next( hfiles[f_i]->GetListOfKeys());
    TKey *key;
    while ((key=(TKey*)next())) histoNames[TOTHIST++] = key->GetName();

    for(int h=0;h<TOTHIST;h++){
      if( histoNames[h].find("Jet")==string::npos) continue;
      histos[f_i][h] = (TH1F*) hfiles[f_i]->FindObjectAny(histoNames[h].c_str());
      histos[f_i][h]->Multiply(weights[f_i]);

      if(f_i==0) resultHistos[h] = (TH1F*)histos[f_i][h]->Clone();
      else resultHistos[h]->Add(histos[f_i][h]);
    }

  }

  TFile outf(summedFileName.c_str(),"recreate");
  for(int h=0;h<TOTHIST;h++){
    if( histoNames[h].find("Jet")==string::npos) continue;
    resultHistos[h]->Write();
  }
  outf.Write();
  outf.Close();

  cout << "Result file is: " << summedFileName << endl;
  //system( command.c_str() );
}




void setBranches(TTree *myTree){
  //event
  varDouble["pthat"] = 0;
  myTree->SetBranchAddress("pthat", &(varDouble["pthat"]) );

  //jets
  varFloatArr["Jet_px"]  =new float[MAXJETS]; varFloatArr["Jet_py"]  =new float[MAXJETS]; varFloatArr["Jet_pz"]  =new float[MAXJETS]; 
  varFloatArr["Jet_E"]  =new float[MAXJETS]; 
  varFloatArr["Jet_pt"]  =new float[MAXJETS];   varFloatArr["Jet_eta"]  =new float[MAXJETS];
  varFloatArr["Jet_phi"]  =new float[MAXJETS];   varFloatArr["Jet_bDiscr"]  =new float[MAXJETS];
  varIntArr["Jet_flavour"]  =new int[MAXJETS];   varIntArr["Jet_IDLoose"]  =new int[MAXJETS]; 
  varFloatArr["Jet_bEff"] = new float[MAXJETS];
  varInt["Jet_n"] = 0;

  myTree->SetBranchAddress("PFJet_nJets", &(varInt["Jet_n"]) );
  myTree->SetBranchAddress("PFJet_pt", varFloatArr["Jet_pt"]);
  myTree->SetBranchAddress("PFJet_px", varFloatArr["Jet_px"]);
  myTree->SetBranchAddress("PFJet_py", varFloatArr["Jet_py"]);
  myTree->SetBranchAddress("PFJet_pz", varFloatArr["Jet_pz"]);
  myTree->SetBranchAddress("PFJet_E", varFloatArr["Jet_E"]);

  myTree->SetBranchAddress("PFJet_phi", varFloatArr["Jet_phi"]);
  myTree->SetBranchAddress("PFJet_eta", varFloatArr["Jet_eta"]);
  myTree->SetBranchAddress(("PFJet_btag_"+BALGO).c_str(), varFloatArr["Jet_bDiscr"]);
  
  myTree->SetBranchAddress(("PFJet_btag_"+BALGOWP+"_beff").c_str(), varFloatArr["Jet_bEff"]);

  myTree->SetBranchAddress("PFJet_JetIDLoose", varIntArr["Jet_IDLoose"]);
  myTree->SetBranchAddress("PFJet_JetFlavour", varIntArr["Jet_flavour"]);
}



//remember to remove nEvents (useless if samples are not filtered)
vector<TFile*> getFiles(string SET, string HLT){
  // xsecs are in nb
  std::vector<TFile*> files;
  float totLumi = 0, grandTotal=0;
  std::vector<float> totalLumis;
  lumiFactor.clear();
  minPthat.clear(); maxPthat.clear(); xsec.clear();
  
  if(SET=="Data"){
    /*
    totLumi = 6900; totalLumis.push_back(totLumi); grandTotal+=totLumi;
    files.push_back( TFile::Open( (DIR+"anV3-MinimumBias-Commissioning10-GOODCOLL-Jun14thSkim_v1-V9.root").c_str() ) ); lumi.push_back(totLumi);
    if(HLT=="HLT_Jet15U") lumiFactor.push_back(1); else if(HLT=="HLT_Jet30U") lumiFactor.push_back(1);
    else if(HLT=="HLT_Jet50U") lumiFactor.push_back(1); else if(HLT=="HLT_Jet70U") lumiFactor.push_back(1);  else if(HLT=="HLT_Jet100U") lumiFactor.push_back(1);
    */
    totLumi = 8000; totalLumis.push_back(totLumi); grandTotal+=totLumi;
    files.push_back( TFile::Open( (DIR+"anV3-MinimumBias-Commissioning10-SD_JetMETTau-Jun14thSkim_v1-V9.root").c_str() ) ); lumi.push_back(totLumi);
    if(HLT=="HLT_Jet15U") lumiFactor.push_back(1); else if(HLT=="HLT_Jet30U") lumiFactor.push_back(1);
    else if(HLT=="HLT_Jet50U") lumiFactor.push_back(1); else if(HLT=="HLT_Jet70U") lumiFactor.push_back(1);  else if(HLT=="HLT_Jet100U") lumiFactor.push_back(1);
    
    totLumi = 4900;totalLumis.push_back(totLumi);grandTotal+=totLumi;
    files.push_back( TFile::Open( (DIR+"anV3-JetMETTau-Run2010A-Jun14thReReco_v2-V9.root").c_str() ) ); lumi.push_back(totLumi);
    if(HLT=="HLT_Jet15U") lumiFactor.push_back(1); else if(HLT=="HLT_Jet30U") lumiFactor.push_back(1);
    else if(HLT=="HLT_Jet50U") lumiFactor.push_back(1); else if(HLT=="HLT_Jet70U") lumiFactor.push_back(1);  else if(HLT=="HLT_Jet100U") lumiFactor.push_back(1);

    totLumi = 58600;totalLumis.push_back(totLumi);grandTotal+=totLumi;
    files.push_back( TFile::Open( (DIR+"anV3-JetMETTau-PromptReco-v4_Runs137437-139558-V9.root").c_str() ) ); lumi.push_back(totLumi);
    if(HLT=="HLT_Jet15U") lumiFactor.push_back(totLumi/4100); else if(HLT=="HLT_Jet30U") lumiFactor.push_back(1);
    else if(HLT=="HLT_Jet50U") lumiFactor.push_back(1); else if(HLT=="HLT_Jet70U") lumiFactor.push_back(1);  else if(HLT=="HLT_Jet100U") lumiFactor.push_back(1);

    totLumi = 82100;totalLumis.push_back(totLumi);grandTotal+=totLumi;
    files.push_back( TFile::Open( (DIR+"anV3-JetMETTau-Run2010A-Jul16thReReco-v1-V9.root").c_str() ) ); lumi.push_back(totLumi);
    if(HLT=="HLT_Jet15U") lumiFactor.push_back(totLumi/1960); else if(HLT=="HLT_Jet30U") lumiFactor.push_back(totLumi/24600);
    else if(HLT=="HLT_Jet50U") lumiFactor.push_back(1); else if(HLT=="HLT_Jet70U") lumiFactor.push_back(1);  else if(HLT=="HLT_Jet100U") lumiFactor.push_back(1);

    totLumi = 64600;totalLumis.push_back(totLumi);grandTotal+=totLumi;
    files.push_back( TFile::Open( (DIR+"anV3-JetMETTau-PromptReco-v4_Runs140160-140388-V9.root").c_str() ) ); lumi.push_back(totLumi);
    if(HLT=="HLT_Jet15U") lumiFactor.push_back(totLumi/1100); else if(HLT=="HLT_Jet30U") lumiFactor.push_back(totLumi/11000);
    else if(HLT=="HLT_Jet50U") lumiFactor.push_back(1); else if(HLT=="HLT_Jet70U") lumiFactor.push_back(1);  else if(HLT=="HLT_Jet100U") lumiFactor.push_back(1);

    totLumi = 41100;totalLumis.push_back(totLumi);grandTotal+=totLumi;
    files.push_back( TFile::Open( (DIR+"anV3-JetMETTau-PromptReco-v4_Runs140389-141887-V9.root").c_str() ) ); lumi.push_back(totLumi);
    if(HLT=="HLT_Jet15U") lumiFactor.push_back(totLumi/760); else if(HLT=="HLT_Jet30U") lumiFactor.push_back(totLumi/7600);
    else if(HLT=="HLT_Jet50U") lumiFactor.push_back(1); else if(HLT=="HLT_Jet70U") lumiFactor.push_back(1);  else if(HLT=="HLT_Jet100U") lumiFactor.push_back(1);

    totLumi = 1450000; totalLumis.push_back(totLumi);grandTotal+=totLumi;
    files.push_back( TFile::Open( (DIR+"anV3-JetMET-Run2010A-PromptReco-v4_Runs141950-143731-V9B.root").c_str() ) ); lumi.push_back(totLumi);
    if(HLT=="HLT_Jet15U") lumiFactor.push_back(totLumi/7000); else if(HLT=="HLT_Jet30U") lumiFactor.push_back(totLumi/140600);
    else if(HLT=="HLT_Jet50U") lumiFactor.push_back(1); else if(HLT=="HLT_Jet70U") lumiFactor.push_back(1);  else if(HLT=="HLT_Jet100U") lumiFactor.push_back(1);

    totLumi = 1342100; totalLumis.push_back(totLumi);grandTotal+=totLumi;
    files.push_back( TFile::Open( (DIR+"anV3-JetMET-Run2010A-PromptReco-v4_Runs143732-144114-V9B.root").c_str() ) ); lumi.push_back(totLumi);
    if(HLT=="HLT_Jet15U") lumiFactor.push_back(totLumi/2500); else if(HLT=="HLT_Jet30U") lumiFactor.push_back(totLumi/49200);
    else if(HLT=="HLT_Jet50U") lumiFactor.push_back(1); else if(HLT=="HLT_Jet70U") lumiFactor.push_back(1);  else if(HLT=="HLT_Jet100U") lumiFactor.push_back(1);

    float grandScaledTotal=0;
    cout <<"--- "<<HLT << " Total Lumi: "<< grandTotal/1000 << "/nb" << endl;
    for (int l=0;l<(int)lumiFactor.size();l++) {
      cout << totalLumis[l] << " " <<lumiFactor[l] << endl;
      grandScaledTotal += totalLumis[l]/lumiFactor[l];
    }
    cout << "--- Effective Lumi: " << grandScaledTotal/1000 << "/nb" << endl;
  }

  else if(SET=="Pythia6"){
    minPthat.push_back(15); maxPthat.push_back(30);
    files.push_back( TFile::Open( (DIR+"anV3-QCD_Pt15_Spring10-V8b.root").c_str() ) ); 
    xsec.push_back(876215000.0); //nEvents.push_back(6090500);

    minPthat.push_back(30); maxPthat.push_back(80);
    files.push_back( TFile::Open( (DIR+"anV3-QCD_Pt30_Spring10-V8b.root").c_str() ) ); 
    xsec.push_back(60411000); //nEvents.push_back(4989664);
  
    minPthat.push_back(80); maxPthat.push_back(170);
    files.push_back( TFile::Open( (DIR+"anV3-QCD_Pt80_Spring10-V8b.root").c_str() ) ); 
    xsec.push_back(923821); //nEvents.push_back(2971800);

    minPthat.push_back(170); maxPthat.push_back(300);
    files.push_back( TFile::Open( (DIR+"anV3-QCD_Pt170_Spring10-V8b.root").c_str() ) ); 
    xsec.push_back(25474.9); //nEvents.push_back(3091950);

    minPthat.push_back(300); maxPthat.push_back(10000000);
    files.push_back( TFile::Open( (DIR+"anV3-QCD_Pt300_Spring10-V8b.root").c_str() ) );
    xsec.push_back(1256); 
  }

  else if(SET=="Herwig"){
    minPthat.push_back(15); maxPthat.push_back(30);
    files.push_back( TFile::Open( (DIR+"anV3-QCD_Pt15-herwig-Summer10-V8b.root").c_str() ) );
    xsec.push_back(893300000); //nEvents.push_back(6090500);

    minPthat.push_back(30); maxPthat.push_back(80);
    files.push_back( TFile::Open( (DIR+"anV3-QCD_Pt30-herwig-Summer10-V8b.root").c_str() ) );
    xsec.push_back(62290000); //nEvents.push_back(4989664);

    minPthat.push_back(80); maxPthat.push_back(170);
    files.push_back( TFile::Open( (DIR+"anV3-QCD_Pt80-herwig-Summer10-V8b.root").c_str() ) );
    xsec.push_back(988700); //nEvents.push_back(2971800);

    minPthat.push_back(170); maxPthat.push_back(300);
    files.push_back( TFile::Open( (DIR+"anV3-QCD_Pt170-herwig-Summer10-V8b.root").c_str() ) );
    xsec.push_back(28030); //nEvents.push_back(3091950);

    minPthat.push_back(300); maxPthat.push_back(10000000);
    files.push_back( TFile::Open( (DIR+"anV3-QCD_Pt300-herwig-Summer10-V8b.root").c_str() ) );
    xsec.push_back(1401);
  }

  else if(SET=="Madgraph"){
    //statistics too low
    //minPthat.push_back(-1); maxPthat.push_back(-1);
    //files.push_back( TFile::Open( (DIR+"anV3-QCD_Pt50to100-madgraph-Summer10-v8b.root").c_str() ) );
    //xsec.push_back(30000000); //nEvents.push_back(6090500);

    minPthat.push_back(-1); maxPthat.push_back(-1);
    files.push_back( TFile::Open( (DIR+"anV3-QCD_Pt100to250-madgraph-Summer10-v8b.root").c_str() ) );
    xsec.push_back(7000000); //nEvents.push_back(6090500);

    minPthat.push_back(-1); maxPthat.push_back(-1);
    files.push_back( TFile::Open( (DIR+"anV3-QCD_Pt250to500-madgraph-Summer10-v8b.root").c_str() ) );
    xsec.push_back(171000); //nEvents.push_back(6090500);

  }


  else cout << "WARNING: "<< SET << " is no valid set " << endl;

  return files;
}


void setTrigger(string HLT){
  jet1_minPt = 0;
  jet1_maxPt = 0;

  if(HLT=="HLT_Jet15U"){ jet1_minPt =56;  jet1_maxPt =84;jet1_minPt_lab ="56";  jet1_maxPt_lab ="84"; }
  else if(HLT=="HLT_Jet30U"){ jet1_minPt =84;  jet1_maxPt =120; jet1_minPt_lab ="84";  jet1_maxPt_lab ="120"; }
  else if(HLT=="HLT_Jet50U"){ jet1_minPt =120;  jet1_maxPt =10000000; jet1_minPt_lab ="120";  jet1_maxPt_lab ="";}

}


void createCommonHistos(string SET){

  string labels[2] = {"Vert","Jet"};
  string labelsC[2] = {"","_CORR"};
  //string fl[] = {"2b","1b1c","1b1l","1c1c","1c1l","1l1l"};

  int nLabels = (int) ((float)sizeof(labels)/(float)sizeof(string));
  int nLabelsC = (int) ((float)sizeof(labelsC)/(float)sizeof(string));

  for(int h=0; h< nLabels; h++){
    for(int hC=0; hC< nLabelsC; hC++){
      histos[labels[h]+"_dR_2b"+labelsC[hC]] = new TH1F((labels[h]+"_dR_2b"+labelsC[hC]).c_str() ,"",60,0,6);  
      histos[labels[h]+"_dR_2b"+labelsC[hC]]->Sumw2();
      histos[labels[h]+"_dPhi_2b"+labelsC[hC]] = new TH1F((labels[h]+"_dPhi_2b"+labelsC[hC]).c_str() ,"",32,0,3.2);  
      histos[labels[h]+"_dPhi_2b"+labelsC[hC]]->Sumw2();
      histos[labels[h]+"_dEta_2b"+labelsC[hC]] = new TH1F((labels[h]+"_dEta_2b"+labelsC[hC]).c_str() ,"",100,0,10);  
      histos[labels[h]+"_dEta_2b"+labelsC[hC]]->Sumw2();
    }
    
    if(SET!="Data"){
      int lab_fl_n= (int) ( (float)sizeof(lab_fl))/((float) sizeof(string));

      for(int f=0;f< lab_fl_n;f++){
	histos[labels[h]+"_dR_2b_"+lab_fl[f]] = new TH1F( (labels[h]+"_dR_2b_"+lab_fl[f]).c_str(),"",60,0,6); histos[labels[h]+"_dR_2b_"+lab_fl[f] ]->Sumw2();
	histos[labels[h]+"_dEta_2b_"+lab_fl[f]] = new TH1F( (labels[h]+"_dEta_2b_"+lab_fl[f]).c_str(),"",100,0,10); histos[labels[h]+"_dEta_2b_"+lab_fl[f] ]->Sumw2();
	histos[labels[h]+"_dPhi_2b_"+lab_fl[f]] = new TH1F( (labels[h]+"_dPhi_2b_"+lab_fl[f]).c_str(),"",32,0,3.2); histos[labels[h]+"_dPhi_2b_"+lab_fl[f] ]->Sumw2();
      }  
    }
  }
  
}



void createJetHistos(string SET){

  histos["Jet_pt"] = new TH1F("Jet_pt" ,"",200,0,2000);  histos["Jet_pt"]->Sumw2();
  histos["Jet_eta"] = new TH1F("Jet_eta" ,"",100,-3,3);  histos["Jet_eta"]->Sumw2();
  histos["Jet_discr"] = new TH1F("Jet_discr","",200,-10,10); histos["Jet_discr"]->Sumw2();
  
  histos["Jet_pt_2b"] = new TH1F("Jet_pt_2b","",200,0,2000);  histos["Jet_pt_2b"]->Sumw2();
  histos["Jet_eta_2b"] = new TH1F("Jet_eta_2b","",100,-3,3);  histos["Jet_eta_2b"]->Sumw2();
  histos["Jet_ptAsymm_2b"] = new TH1F("Jet_ptAsymm_2b","",100,0,1); histos["Jet_ptAsymm_2b"]->Sumw2();

  histos["Jet_bJets_pt"] = new TH1F("Jet_bJets_pt","",200,0,2000); histos["Jet_bJets_pt"]->Sumw2();
  histos["Jet_bJets_eta"] = new TH1F("Jet_bJets_eta","",100,-3,3); histos["Jet_bJets_eta"]->Sumw2();
  histos["Jet_bJets_n"] = new TH1F("Jet_bJets_n","",10,0,10); histos["Jet_bJets_n"]->Sumw2();

  histos["Jet_2bJets_Mass"] = new TH1F("Jet_2bJets_Mass",";mass(bjet,bjet;events)",200,0,2000); histos["Jet_2bJets_Mass"]->Sumw2();

  histos2D["Jet_dEtadPhi_2b"] = new TH2F("Jet_dEtadPhi_2b",";#Delta #eta;#Delta #phi",100,0,10,32,0,3.2); histos2D["Jet_dEtadPhi_2b"]->Sumw2();
  histos2D["Jet_pt1pt2_2b"] = new TH2F("Jet_pt1pt2_2b",";p_{T}^{1};p_{T}^{2}",200,0,2000,200,0,2000); histos2D["Jet_pt1pt2_2b"]->Sumw2();

  //string fl[] = {"2b","1b1c","1b1l","1c1c","1c1l","1l1l"};
  if(SET!="Data"){
    int lab_fl_n= (int) ( (float)sizeof(lab_fl))/((float) sizeof(string));
    for(int f=0;f<lab_fl_n;f++){
      histos["Jet_pt_2b_"+lab_fl[f]] = new TH1F( ("Jet_pt_2b_"+lab_fl[f]).c_str(),"",200,0,2000); histos["Jet_pt_2b_"+lab_fl[f] ]->Sumw2();
      histos["Jet_2bJets_Mass_"+lab_fl[f]] = new TH1F( ("Jet_2bJets_Mass_"+lab_fl[f]).c_str(),";mass(bjet,bjet;events)",200,0,2000); histos["Jet_2bJets_Mass_"+lab_fl[f]]->Sumw2();
      histos2D["Jet_dEtadPhi_2b_"+lab_fl[f]] = new TH2F(("Jet_dEtadPhi_2b_"+lab_fl[f]).c_str(),";#Delta #eta;#Delta #phi",100,0,10,32,0,3.2); histos2D["Jet_dEtadPhi_2b_"+lab_fl[f]]->Sumw2();
      histos2D["Jet_pt1pt2_2b_"+lab_fl[f]] = new TH2F(("Jet_pt1pt2_2b_"+lab_fl[f]).c_str(),";p_{T}^{1};p_{T}^{2}",200,0,2000,200,0,2000); histos2D["Jet_pt1pt2_2b_"+lab_fl[f]]->Sumw2();

    }  
    
    histos["Jet_bJetsMC_dR"] = new TH1F("Jet_bJetsMC_dR","",60,0,6); histos["Jet_bJetsMC_dR"]->Sumw2();
    histos["Jet_bJetsMC_dPhi"] = new TH1F("Jet_bJetsMC_dPhi","",100,0,10); histos["Jet_bJetsMC_dPhi"]->Sumw2();
    histos["Jet_bJetsMC_dEta"] = new TH1F("Jet_bJetsMC_dEta","",32,0,3.2); histos["Jet_bJetsMC_dEta"]->Sumw2();

    
    histos["Jet_bJetsMC_pt"] = new TH1F("Jet_bJetsMC_pt","",400,0,2000); histos["Jet_bJetsMC_pt"]->Sumw2();
    histos["Jet_cJetsMC_pt"] = new TH1F("Jet_cJetsMC_pt","",400,0,2000); histos["Jet_cJetsMC_pt"]->Sumw2();
    histos["Jet_lJetsMC_pt"] = new TH1F("Jet_lJetsMC_pt","",400,0,2000); histos["Jet_lJetsMC_pt"]->Sumw2();
    histos["Jet_bJetsMC_eta"] = new TH1F("Jet_bJetsMC_eta","",100,-3,3); histos["Jet_bJetsMC_eta"]->Sumw2();
    histos["Jet_cJetsMC_eta"] = new TH1F("Jet_cJetsMC_eta","",100,-3,3); histos["Jet_cJetsMC_eta"]->Sumw2();
    histos["Jet_lJetsMC_eta"] = new TH1F("Jet_lJetsMC_eta","",100,-3,3); histos["Jet_lJetsMC_eta"]->Sumw2();
  }


}

void createIVFHistos(string SET){

  histos["bVert_N"] = new TH1F("Vert_bVert_N","for events with 2B or 2V",10,0,10);  histos["bVert_N"]->Sumw2();
  histos["Vert_goodEv"] = new TH1F("Vert_goodEv","",61000,0,61000); histos["Vert_goodEv"]->Sumw2();
  histos["Vert_dRingoodEv"] = new TH1F("Vert_dRingoodEv","",61000,0,61000); histos["Vert_dRingoodEv"]->Sumw2();
  histos["Vert_pt_selvert"] = new TH1F("Vert_pt_goodVert","",200,0,200);  histos["Vert_pt_selvert"]->Sumw2();
  
  histos["trigger_ratio_50_30"] = new TH1F("Trigger_ratio_50_30","",20000,130000,150000);  histos["trigger_ratio_50_30"]->Sumw2();
  histos["trigger_temp_50"] = new TH1F("Trigger_temp_50","",20000,130000,150000); histos["trigger_temp_50"]->Sumw2();
  histos["trigger_temp_5030"] = new TH1F("Trigger_temp_30","",20000,130000,150000); histos["trigger_temp_5030"]->Sumw2();
  histos["trigger_ratio_50_15"] = new TH1F("Trigger_ratio_50_15","",20000,130000,150000); histos["trigger_ratio_50_15"]->Sumw2();
  histos["trigger_temp_5015"] = new TH1F("Trigger_temp_15","",20000,130000,150000); histos["trigger_temp_5015"]->Sumw2();
  histos["trigger_ratio_50_6"] = new TH1F("Trigger_ratio_50_6","",20000,130000,150000); histos["trigger_ratio_50_6"]->Sumw2();
  histos["trigger_temp_506"] = new TH1F("Trigger_temp_6","",20000,130000,150000);  histos["trigger_temp_506"]->Sumw2();
  
  if(SET!="Data"){
    histos2D["Vert_flavCat_vs_matCat"] = new TH2F("Vert_flavCat_vs_matCat","",6,0,6,5,0,5); histos2D["Vert_flavCat_vs_matCat"]->Sumw2(); 

    int lab_eff_n= (int) ( (float)sizeof(lab_eff))/((float) sizeof(string));
    for(int f=0;f<lab_eff_n;f++){
      string l = lab_eff[f];
      histos["Vert_dR_"+l] =  new TH1F( ("Vert_dR_"+l).c_str(),"",60,0,6); histos["Vert_dR_"+l]->Sumw2();
      histos["Vert_dPhi_"+l] =  new TH1F( ("Vert_dPhi_"+l).c_str(),"",32,0,3.2); histos["Vert_dPhi_"+l]->Sumw2();
      histos["Vert_dEta_"+l] =  new TH1F( ("Vert_dEta_"+l).c_str(),"",100,0,10); histos["Vert_dEta_"+l]->Sumw2();
    }

    int lab_cat_n= (int) ( (float)sizeof(lab_cat))/((float) sizeof(string));
    for(int f=0;f<lab_cat_n;f++){
      string l = lab_cat[f];
      histos["Vert_dR_2b_"+l] =  new TH1F( ("Vert_dR_2b_"+l).c_str(),"",60,0,6); histos["Vert_dR_2b_"+l]->Sumw2();
      histos["Vert_dPhi_2b_"+l] =  new TH1F( ("Vert_dPhi_2b_"+l).c_str(),"",32,0,3.2); histos["Vert_dPhi_2b_"+l]->Sumw2();
      histos["Vert_dEta_2b_"+l] =  new TH1F( ("Vert_dEta_2b_"+l).c_str(),"",100,0,10); histos["Vert_dEta_2b_"+l]->Sumw2();
    }
     
    int lab_proc_n= (int) ( (float)sizeof(lab_proc))/((float) sizeof(string));
    for(int f=0;f<lab_proc_n;f++){
      string l = lab_proc[f];
      histos["Vert_dR_2b_"+l] =  new TH1F( ("Vert_dR_2b_"+l).c_str(),"",60,0,6); histos["Vert_dR_2b_"+l]->Sumw2();
      histos["Vert_dPhi_2b_"+l] =  new TH1F( ("Vert_dPhi_2b_"+l).c_str(),"",32,0,3.2); histos["Vert_dPhi_2b_"+l]->Sumw2();
      histos["Vert_dEta_2b_"+l] =  new TH1F( ("Vert_dEta_2b_"+l).c_str(),"",100,0,10); histos["Vert_dEta_2b_"+l]->Sumw2();
    }
  }

}

void corrFunctionIVF(string strigger){
  TFile *fCorr = TFile::Open(ivfCorrFile.c_str(),"UPDATE");
  fCorr->GetObject((strigger+"/corrfunc").c_str(),hcorrIVF);
  xMax_IVFCorr = hcorrIVF->GetXaxis()->GetXmax();
  xMin_IVFCorr = hcorrIVF->GetXaxis()->GetXmin();
  binW_IVFCorr = (xMax_IVFCorr-xMin_IVFCorr)/hcorrIVF->GetXaxis()->GetNbins();
  
}

float delta_phi(float phi1, float phi2){

  Float_t pi_greco = 3.1415;
  Float_t adphi = fabs(phi1 -phi2);
  Float_t dphi = (phi1 -phi2);
  if (adphi > fabs(2*pi_greco-adphi)){
    if(dphi > 0){ dphi -= 2*pi_greco;}
    if(dphi < 0){ dphi += 2*pi_greco;}
    dphi = fabs(2*pi_greco-adphi);
  }

  return fabs(dphi);
}

float delta_R(float phi1, float phi2, float eta1, float eta2){
  float dPhi = delta_phi(phi1,phi2);
  float dEta = eta1 - eta2;
  float dR = dEta*dEta + dPhi*dPhi;

  return sqrt(dR);
}
