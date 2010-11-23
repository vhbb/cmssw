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
#include "TF2.h"
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
////general selection cuts
//Jets
float BJET_PT = 30;
float BJET_ETA = 2;
float JET1_ETA = 3;
float JET_PT = BJET_PT;
float JET_ETA = 2;
//Bhads
float BHAD_ETA = 2;
float BHAD_PT = 15;
//Bcands
float BCAND_ETA = 2;
float BCAND_PT = 8;

//btag cuts
string BALGO = "ssvhp";
string BALGOWP = "ssvhpt";
float BDISCR = 2;
float BTAGEFF_Scaling = 0.9;//0.9; //btag eff scaling (andrea docet)
float BTAGEFF_Error = 0.1;//0.1; //btag eff relative flat error (andrea docet)


//normalization in pb
double LUMI = 1;
bool applyPrescale = false;

//label for the outfile
string LABEL = "-BBCorr_scaledEff"; //add a trailing - please, eg -mylabel

//where to find the ttrees
string DIR="/scratch/leo/";

//file containing btag purity fit
string jetCorrFile = "histo-btagPurity-Pythia6.root";
TFile *f_jetPurity;
//TF1 *f_jet_pur_pt;
TH1F * Jet_BPurity_pt;
// histo with jet/bhadron efficiency (defined with getBHadJetEff())
TH1F * Jet_BHadrJetEff_pt;
//file with dR corrections: the actual file is called e.g. BTag_dR_correction_HLT_Jet30U.root
//string jetDrCorrFile = "BTag_dR_corr_BB";
TH1F * Jet_BHadrJetEff_dR;


//file containing the IVF corrections:
string ivfCorrFile = "histo-bvertEffic-Pythia6.root";

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
string lab_eff[] = {"efficnum","efficden","puritynum","puritydenNONBB","puritydenBB"};
string lab_cat[] = {"cat2bmat","cat2bnotmat","cat2bnotkin","cat0b","catmanyb","catnocat"};
string lab_proc[] = {"fcr","fex","gsp","otherproc"};


//functions
void createCommonHistos(string SET, string HLT);
vector<TFile*> getFiles(string SET, string HLT);
void setTrigger(string HLT);
void setBranches(TTree* myTree);
void createJetHistos(string SET);
void createIVFHistos(string SET);
int analyzeJetEvent(event thisEvent,string SET, double W);
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

    //Selecting a reduced tree
    TTree *myTree = (TTree*)files[f_i]->Get("bcanalyzer/tEvents");
    int events_n = myTree->GetEntries();
    string mySel = "";
    if(HLT!="") mySel += HLT+"==1 ";
    //if(SET!="Data") mySel=""; //temporary fix for HLT collections issue
    //mySel += " && PFJet_JetIDLoose[0]==1 && PFJet_pt[0]>= "+jet1_minPt_lab;
    myTree->Draw(">>selList", mySel.c_str()); 
    TEventList *mySelList = (TEventList*)gDirectory->Get("selList"); 
    myTree->SetEventList(mySelList);

    fname = fname.replace( fname.begin(), fname.begin()+fname.rfind("/")+1,"" );
    struct stat stDirInfo;
    if(stat(SET.c_str(),&stDirInfo) != 0) system(("mkdir "+SET).c_str());
    TFile *f_out = TFile::Open((SET+"/histo-"+HLT+LABEL+"-"+fname ).c_str(),"RECREATE");
    cout << "Output file: " << SET+"/histo-"+HLT+LABEL+"-"+fname << endl;

    //booking histos
    createCommonHistos(SET,HLT);
    createJetHistos(SET);
    createIVFHistos(SET);

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
      else if(HLT=="HLT_L1Jet6U") corrFunctionIVF("jet6");
      else if(HLT=="HLT_Jet70U") corrFunctionIVF("jet70");
      else if(HLT=="HLT_Jet100U") corrFunctionIVF("jet100");
    }else{
      cout << "Warning: IVF correction file not found" << endl;
    }
    //BTag eff correction
    if(stat(jetCorrFile.c_str(),&stFileInfo) == 0){
      f_jetPurity = TFile::Open(jetCorrFile.c_str());
      f_jetPurity->GetObject( ("BHadrEff_dR_"+HLT).c_str(), Jet_BHadrJetEff_dR);
      f_jetPurity->GetObject("Jet_BPurity_pt", Jet_BPurity_pt);

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

      if( fabs(  varFloatArr["Jet_eta"][0]) > JET1_ETA) continue;

      //ljet hlt plot
      histos["Jet_Jet1ptHLT"]->Fill( varFloatArr["Jet_pt"][0],W);
      
      if( varFloatArr["Jet_pt"][0]<jet1_minPt || varFloatArr["Jet_pt"][0] > jet1_maxPt ) continue;

      //event analyzers
      int ret_code=0;
      ret_code = analyzeJetEvent(thisEv, SET, W);
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
int analyzeJetEvent(event thisEv, string SET, double W){
  //  cout << "JET" << endl;
  int bJets_i[100];
  int bJetsMC_i[100];
  int n_bJets = 0;
  int n_bJetsMC = 0;

  float bestBH1_dR=1000000, bestBH2_dR=1000000;
  int bestBH1_i = -1, bestBH2_i=-1;
  float measure = 0.5;

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

    //bhadron matching
    float dR1 = delta_R( varFloatArr["Jet_phi"][j],thisEv.phiB1,  varFloatArr["Jet_eta"][j],thisEv.etaB1);
    float dR2 = delta_R( varFloatArr["Jet_phi"][j],thisEv.phiB2,  varFloatArr["Jet_eta"][j],thisEv.etaB2);

    if(dR1 < dR2 && dR1<bestBH1_dR){ bestBH1_dR= dR1; bestBH1_i = j;}
    else if(dR1 > dR2 && dR2<bestBH2_dR) { bestBH2_dR= dR2; bestBH2_i = j;}

  }//end jet loop

  histos["Jet_bJets_n"]->Fill( n_bJets ,W);

  //True MC bjet content
  if(SET!="Data" ){
    
    if(n_bJetsMC==2){
    int bj1 = bJetsMC_i[0]; int bj2 = bJetsMC_i[1];
    float dR = delta_R( varFloatArr["Jet_phi"][ bj1 ],varFloatArr["Jet_phi"][ bj2 ],  varFloatArr["Jet_eta"][ bj1 ],varFloatArr["Jet_eta"][ bj2 ]  );
    float dPhi = delta_phi(  varFloatArr["Jet_phi"][ bj1 ],varFloatArr["Jet_phi"][ bj2 ] );
    float dEta = fabs( varFloatArr["Jet_eta"][ bj1 ] - varFloatArr["Jet_eta"][ bj2 ]);

    histos["Jet_2bJetsMC_dR"]->Fill(dR,W);
    histos["Jet_2bJetsMC_dEta"]->Fill(dEta,W);
    histos["Jet_2bJetsMC_dPhi"]->Fill(dPhi,W);
    histos["Jet_2bJetsMC_pt"]->Fill(varFloatArr["Jet_pt"][ bj2 ],W);
    histos["Jet_2bJetsMC_pt"]->Fill(varFloatArr["Jet_pt"][ bj1 ],W);
    }
    
  //jet-MC Bhadr eff
    if( thisEv.nB==2 ) {
      //BHadron - BJet eff: matching
      if(thisEv.ptB1>BHAD_PT && thisEv.ptB2>BHAD_PT && fabs(thisEv.etaB1)<BHAD_ETA && fabs(thisEv.etaB2)<BHAD_ETA){
	if(bestBH1_dR < measure){
	  histos["BHadr_BHadrJetMatched_pt"]->Fill(thisEv.ptB1,W);
	  histos["BHadr_BHadrJetMatched_eta"]->Fill(thisEv.etaB1,W);
	  histos["Jet_BHadrJetMatched_pt"]->Fill( varFloatArr["Jet_pt"][bestBH1_i] ,W);
	}
	if(bestBH2_dR < measure){
	  histos["BHadr_BHadrJetMatched_pt"]->Fill(thisEv.ptB2,W);
	  histos["BHadr_BHadrJetMatched_eta"]->Fill(thisEv.etaB2,W);
	  histos["Jet_BHadrJetMatched_pt"]->Fill( varFloatArr["Jet_pt"][bestBH2_i] ,W);
	}
	if(bestBH1_dR < measure && bestBH2_dR < measure){
	  histos["BHadr_BHadrJetMatched_dR"]->Fill(thisEv.dRbb,W);
	  histos["BHadr_BHadrJetMatched_dPhi"]->Fill(thisEv.dPhibb,W);
	  histos["BHadr_BHadrJetMatched_dEta"]->Fill(thisEv.dEtabb,W);

	  float dR = delta_R( varFloatArr["Jet_phi"][ bestBH1_i ],varFloatArr["Jet_phi"][ bestBH2_i ],  varFloatArr["Jet_eta"][ bestBH1_i ],varFloatArr["Jet_eta"][ bestBH2_i ]  );
	  histos["Jet_BHadrJetMatched_dR"]->Fill(dR,W);
	}

	histos["BHadr_pt"]->Fill(thisEv.ptB1,W);
	histos["BHadr_eta"]->Fill(thisEv.etaB1,W);
	histos["BHadr_pt"]->Fill(thisEv.ptB2,W);
	histos["BHadr_eta"]->Fill(thisEv.etaB2,W);
	
	histos["BHadr_dR"]->Fill(thisEv.dRbb,W);
	histos["BHadr_dPhi"]->Fill(thisEv.dPhibb,W);
	histos["BHadr_dEta"]->Fill(thisEv.dEtabb,W);
      }
    }
  }
  

  //2 bjets part
  if(n_bJets!=2) return 1;
  int bj1 = bJets_i[0]; int bj2 = bJets_i[1];
  float dR = delta_R( varFloatArr["Jet_phi"][ bj1 ],varFloatArr["Jet_phi"][ bj2 ],  varFloatArr["Jet_eta"][ bj1 ],varFloatArr["Jet_eta"][ bj2 ]  );
  histos["Jet_dR_2b"]->Fill(dR,W);
  histos2D["Jet_dRvsJet1pt_2b"]->Fill(dR,varFloatArr["Jet_pt"][ 0],W);
  float dPhi = delta_phi(  varFloatArr["Jet_phi"][ bj1 ],varFloatArr["Jet_phi"][ bj2 ] );
  histos["Jet_dPhi_2b"]->Fill(dPhi,W);
  float dEta = fabs( varFloatArr["Jet_eta"][ bj1 ] - varFloatArr["Jet_eta"][ bj2 ]);
  histos["Jet_dEta_2b"]->Fill( dEta ,W);
  
  histos["Jet_pt_2b"]->Fill(  varFloatArr["Jet_pt"][bj1],W);
  histos["Jet_pt_2b"]->Fill(  varFloatArr["Jet_pt"][bj2],W);
  histos["Jet_eta_2b"]->Fill(  varFloatArr["Jet_eta"][bj1],W);
  histos["Jet_eta_2b"]->Fill(  varFloatArr["Jet_eta"][bj2],W);
  
  float ptAsymm = fabs(varFloatArr["Jet_pt"][bj1]-varFloatArr["Jet_pt"][bj2])/( varFloatArr["Jet_pt"][bj1]+varFloatArr["Jet_pt"][bj2]);
  histos["Jet_ptAsymm_2b"]->Fill( fabs(varFloatArr["Jet_pt"][bj1]-varFloatArr["Jet_pt"][bj2])/( varFloatArr["Jet_pt"][bj1]+varFloatArr["Jet_pt"][bj2]),W);

  //corrections
  float btagEff_corr = varFloatArr["Jet_bEff"][0]*varFloatArr["Jet_bEff"][1];
  if(SET=="Data") btagEff_corr *= BTAGEFF_Scaling; 
  int purBin1 = Jet_BPurity_pt->FindBin(varFloatArr["Jet_pt"][bj1]);
  int purBin2 = Jet_BPurity_pt->FindBin(varFloatArr["Jet_pt"][bj2]);
  float meanPur = Jet_BPurity_pt->GetBinContent(purBin1)*Jet_BPurity_pt->GetBinContent(purBin2);
  float purErr = (Jet_BPurity_pt->GetBinError(purBin1)/Jet_BPurity_pt->GetBinContent(purBin1))*(Jet_BPurity_pt->GetBinError(purBin1)/Jet_BPurity_pt->GetBinContent(purBin1)) +
    (Jet_BPurity_pt->GetBinError(purBin2)/Jet_BPurity_pt->GetBinContent(purBin2))*(Jet_BPurity_pt->GetBinError(purBin2)/Jet_BPurity_pt->GetBinContent(purBin2));
  //corr errors (relative, quadratic)
  float corrErr = (varFloatArr["Jet_bEffErr"][0]*varFloatArr["Jet_bEffErr"][0]) +
    (varFloatArr["Jet_bEffErr"][1]*varFloatArr["Jet_bEffErr"][1]);
  if(SET=="Data" && BTAGEFF_Error!=0) corrErr = BTAGEFF_Error*BTAGEFF_Error;
  corrErr += purErr;

  float newW = meanPur/btagEff_corr;

  histos["Jet_dR_2b_BTAGCORR"]->Fill(dR,newW*W);
  histos["Jet_dPhi_2b_BTAGCORR"]->Fill(dPhi,newW*W);
  histos["Jet_dEta_2b_BTAGCORR"]->Fill(dEta,newW*W);
  histos["Jet_ptAsymm_2b_BTAGCORR"]->Fill( ptAsymm,newW*W );

  int effBin1 = Jet_BHadrJetEff_dR->FindBin(dR);
  float eff1 = Jet_BHadrJetEff_dR->GetBinContent(effBin1);   
  //float eff1=0; int effBin1=-1;

  if(eff1!=0){
    float effErr1 = Jet_BHadrJetEff_dR->GetBinError(effBin1)/eff1;
    corrErr += effErr1*effErr1;

    histos["Jet_dR_2b_ErrCORR"]->Fill(dR,corrErr);
    histos["Jet_dEta_2b_ErrCORR"]->Fill(dEta,corrErr);
    histos["Jet_dPhi_2b_ErrCORR"]->Fill(dPhi,corrErr);

    histos["Jet_dR_2b_CORR"]->Fill(dR,newW*W/eff1);
    histos["Jet_dPhi_2b_CORR"]->Fill(dPhi,newW*W/eff1);
    histos["Jet_dEta_2b_CORR"]->Fill(dEta,newW*W/eff1);

    histos2D["Jet_dRvsJet1pt_2b_CORR"]->Fill(dR,varFloatArr["Jet_pt"][ 0],W);
  }

  //invariant mass
  float invM = (varFloatArr["Jet_E"][bj1]+varFloatArr["Jet_E"][bj2])*(varFloatArr["Jet_E"][bj1]+varFloatArr["Jet_E"][bj2]);
  invM -= (varFloatArr["Jet_px"][bj1]+varFloatArr["Jet_px"][bj2])*(varFloatArr["Jet_px"][bj1]+varFloatArr["Jet_px"][bj2]);
  invM -= (varFloatArr["Jet_py"][bj1]+varFloatArr["Jet_py"][bj2])*(varFloatArr["Jet_py"][bj1]+varFloatArr["Jet_py"][bj2]);
  invM -= (varFloatArr["Jet_pz"][bj1]+varFloatArr["Jet_pz"][bj2])*(varFloatArr["Jet_pz"][bj1]+varFloatArr["Jet_pz"][bj2]);
  invM = sqrt(invM);
  histos["Jet_2bJets_Mass"]->Fill(invM,W);
  
  //flavour dependent an
  if(SET!="Data"){
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
    histos["Jet_ptAsymm_2b_"+diFl]->Fill( ptAsymm,W );

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

  //////////////
  //NO! THERE ARE HISTOS WITH ONLY B OR ONLY VERTICES REQUIRED!!!
//   if( thisEv.nV!=2 ) return 1;
  
//   if(SET!="Data" && thisEv.nB!=2) return 1;
//   if(SET!="Data" && fabs(thisEv.etaB1)>=BHAD_ETA && fabs(thisEv.etaB2)>=BHAD_ETA) return 1;
  //////////////

  //not 2V AND not 2B -->skip event
  if(thisEv.nB!=2 && thisEv.nV!=2) return 1; 


  histos["bVert_N"]->Fill(thisEv.nV, W);

  int flavCat = -1, matCat=-1; 
//   if(thisEv.dRvv>-700 && fabs(thisEv.etaV1)<BCAND_ETA && fabs(thisEv.etaV2)<BCAND_ETA){
  if(thisEv.dRvv>-700){
    //mass sum cut
    if((thisEv.massV1 + thisEv.massV2) > 4.5 && thisEv.ptV1>BCAND_PT && thisEv.ptV2>BCAND_PT){//ADD PT CUT ON VERTEX HERE and/or ETA CUT ON VERTEX
      histos["Vert_dR_2b"]->Fill(thisEv.dRvv,W);
      histos2D["Vert_dRvsJet1pt_2b"]->Fill(thisEv.dRvv,varFloatArr["Jet_pt"][ 0],W);

      histos["Vert_dPhi_2b"]->Fill(fabs(thisEv.dPhivv),W);
      histos["Vert_dEta_2b"]->Fill(fabs(thisEv.dEtavv),W);
      histos["Vert_goodEv"]->Fill(thisEv.eventNo,W);
      histos["Vert_dRingoodEv"]->Fill(thisEv.eventNo, thisEv.dRvv*W); //for Lukas: is this a TH1F or a TH2F?
      histos["Vert_pt_selvert"]->Fill(thisEv.ptV1,W); 
      histos["Vert_pt_selvert"]->Fill(thisEv.ptV2,W); 

      //corrected plots
      if(hcorrIVF!=0){
// 	double weight = hcorrIVF->GetBinContent((int)((thisEv.dRvv-xMin_IVFCorr)/binW_IVFCorr+1.)); 
	double weight = hcorrIVF->GetBinContent(hcorrIVF->FindBin(thisEv.dRvv)); 

	histos["Vert_dR_2b_CORR"]->Fill(thisEv.dRvv,W*weight);
	histos2D["Vert_dRvsJet1pt_2b_CORR"]->Fill(thisEv.dRvv,varFloatArr["Jet_pt"][ 0],W*weight);

	histos["Vert_dPhi_2b_CORR"]->Fill(fabs(thisEv.dPhivv),W*weight);
	histos["Vert_dEta_2b_CORR"]->Fill(fabs(thisEv.dEtavv),W*weight);

// 	double weightError = hcorrIVF->GetBinError((int)((thisEv.dRvv-xMin_IVFCorr)/binW_IVFCorr+1.) )/weight;
	double weightError = hcorrIVF->GetBinError(hcorrIVF->FindBin(thisEv.dRvv));
// 	double we2 = weightError*weightError; //-->Add the errors linearly due to correlation
	histos["Vert_dR_2b_ErrCORR"]->Fill(thisEv.dRvv,weightError);
	histos["Vert_dEta_2b_ErrCORR"]->Fill(thisEv.dEtavv,weightError);
        histos["Vert_dPhi_2b_ErrCORR"]->Fill(thisEv.dPhivv,weightError);


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
	  if(thisEv.ptB1>BHAD_PT && thisEv.ptB2>BHAD_PT && fabs(thisEv.etaB1)<BHAD_ETA && fabs(thisEv.etaB2)<BHAD_ETA){
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
	if(thisEv.ptV1>BCAND_PT && thisEv.ptV2>BCAND_PT){
	  if(thisEv.nB==0){
	    histos["Vert_dR_puritydenNONBB"]->Fill(thisEv.dRvv,W);
	    histos["Vert_dPhi_puritydenNONBB"]->Fill(fabs(thisEv.dPhivv),W);
	    histos["Vert_dEta_puritydenNONBB"]->Fill(fabs(thisEv.dEtavv),W);
	  }
	  else{
	    histos["Vert_dR_puritydenBB"]->Fill(thisEv.dRvv,W);
	    histos["Vert_dPhi_puritydenBB"]->Fill(fabs(thisEv.dPhivv),W);
	    histos["Vert_dEta_puritydenBB"]->Fill(fabs(thisEv.dEtavv),W);
	  }

	  //2B
	  if(thisEv.nB==2){
	    //B kinematics
	    if(thisEv.ptB1>BHAD_PT && thisEv.ptB2>BHAD_PT && fabs(thisEv.etaB1)<BHAD_ETA && fabs(thisEv.etaB2)<BHAD_ETA){
	      if(thisEv.nMat==2){
		histos["Vert_dR_efficnum"]->Fill(thisEv.dRbb,W);
		histos["Vert_dR_puritynum"]->Fill(thisEv.dRbb,W);
		histos["Vert_dPhi_efficnum"]->Fill(fabs(thisEv.dPhibb),W);
		histos["Vert_dPhi_puritynum"]->Fill(fabs(thisEv.dPhibb),W);
		histos["Vert_dEta_efficnum"]->Fill(fabs(thisEv.dEtabb),W);
		histos["Vert_dEta_puritynum"]->Fill(fabs(thisEv.dEtabb),W);
	      }
	    }
	  }
	}
      }//end only for MC


//       if(flavCat!=-1 && matCat!=-1)
// 	histos2D["Vert_flavCat_vs_matCat"]->Fill(flavCat,matCat,W);
//       }
    }//end mass sum and pt cut
  }//end if dRvv>-777 --> 2 vertex

  //2B, Bkinematics
  if(SET!="Data" && thisEv.nB==2 && thisEv.ptB1>BHAD_PT && thisEv.ptB2>BHAD_PT && fabs(thisEv.etaB1)<BHAD_ETA && fabs(thisEv.etaB2)<BHAD_ETA) {
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
  string summedFileName = "histo-"+SET+"-"+HLT+LABEL+"-total.root";
  string command = "hadd -f "+summedFileName+" ";
  for(int f_i=init_i; f_i< end_i; f_i++){
    //the output file
    string fname = files[f_i]->GetName();
    fname = fname.replace( fname.begin(), fname.begin()+fname.rfind("/")+1,"" );
    fname=SET+"/histo-"+HLT+LABEL+"-"+fname;
    command += fname+" ";
  }

  system( command.c_str() );
}

/// summing data histos with loaded prescales
void sumHistoPrescaled(string SET, string HLT){
  files = getFiles(SET,HLT);

  //selecting files
  int init_i=0;
  int end_i = files.size();

  string summedFileName = "histo-"+SET+"-"+HLT+LABEL+"-total-wPrescales.root";

  string histoNames[100];
  TH1F * resultHistos[100];
  TH1F * histos[100][100];

  string histo2DNames[100];
  TH2F * resultHistos2D[100];
  TH2F * histos2D[100][100];


  TF1 *weights[10];
  TF2 *weights2D[10];

  int TOTHIST=0;
  int TOTHIST2D=0;

  TFile *hfiles[100];

  for(int f_i=init_i; f_i< end_i; f_i++){
    //TOTHIST=0;
    string fname = files[f_i]->GetName();
    fname = fname.replace( fname.begin(), fname.begin()+fname.rfind("/")+1,"" );
    fname=SET+"/histo-"+HLT+LABEL+"-"+fname;
    hfiles[f_i] = TFile::Open(fname.c_str());
    
    char*ws = new char[30]; sprintf(ws, "%.4f", lumiFactor[f_i] );    
    char*ws2D = new char[30]; sprintf(ws2D, "x-x+y-y+%.4f", lumiFactor[f_i] );
    char*is = new char[30]; sprintf(is, "w%d", f_i );
    char*is2D = new char[30]; sprintf(is2D, "w%d", f_i );

    //cout << fname << " " << ws << " " << is << endl;
    weights[f_i] = new TF1(is,ws,-10000,10000); 
    weights2D[f_i] = new TF2(is2D,ws2D,-10000,10000,-10000,10000);

    if(f_i==init_i){
      TIter next( hfiles[f_i]->GetListOfKeys());

      TKey *key;
      while ((key=(TKey*)next())) {
	if((string)key->GetClassName()=="TH1F")
	  histoNames[TOTHIST++] = key->GetName();
	else if((string)key->GetClassName()=="TH2F"){
	  histo2DNames[TOTHIST2D++] = key->GetName();
	}
      }
    }
    
    for(int h=0;h<TOTHIST;h++){
      histos[f_i][h] = (TH1F*) hfiles[f_i]->FindObjectAny(histoNames[h].c_str());
      histos[f_i][h]->Multiply(weights[f_i]);

      if(f_i==init_i) resultHistos[h] = (TH1F*)histos[f_i][h]->Clone();
      else resultHistos[h]->Add(histos[f_i][h]);
    }

    for(int h=0;h<TOTHIST2D;h++){
      histos2D[f_i][h] = (TH2F*) hfiles[f_i]->FindObjectAny(histo2DNames[h].c_str());
      histos2D[f_i][h]->Multiply(weights2D[f_i]);
      if(f_i==init_i) resultHistos2D[h] = (TH2F*)histos2D[f_i][h]->Clone();
      else resultHistos2D[h]->Add(histos2D[f_i][h]);
    }
  }

  TFile outf(summedFileName.c_str(),"recreate");
  for(int h=0;h<TOTHIST;h++){
    //if( histoNames[h].find("Jet")==string::npos) continue;
    resultHistos[h]->Write();
  }

  for(int h=0;h<TOTHIST2D;h++){
    //if( histo2DNames[h].find("Jet")==string::npos) continue;
    resultHistos2D[h]->Write();
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
  varFloatArr["Jet_bEff"] = new float[MAXJETS];   varFloatArr["Jet_bEffErr"] = new float[MAXJETS];
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
  myTree->SetBranchAddress(("PFJet_btag_"+BALGOWP+"_befferr").c_str(), varFloatArr["Jet_bEffErr"]);


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
    if(HLT=="HLT_L1Jet6U"){
      cout << "please insert here the proper files" << endl;
      //JetMETTauMonitor/Run2010A-PromptReco-v4/RECO Runs 140160-140388
      totLumi = 64610.239; totalLumis.push_back(totLumi);
      grandTotal+=totLumi;
      files.push_back( TFile::Open( (DIR+"anV3-JetMETTauMonitor-Run2010A_140160-140388_v1-v9.root").c_str() ) );
      lumi.push_back(totLumi);
      lumiFactor.push_back(totLumi/65.742);

      //JetMETTauMonitor/Run2010A-Jun14thReReco_v1/RECO Runs 132440-137028
      //Edit totLumi Value (inverse micro barn)
      totLumi = 12949.814; totalLumis.push_back(totLumi);
      grandTotal+=totLumi;
      //Edit FileName, i.e. change the *.root to the actual file
      files.push_back( TFile::Open( (DIR+"anV3-JetMETTauMonitor-Jun14thReReco-Run2010A_132440-137028_v1-v9.root").c_str()
				   ) ); lumi.push_back(totLumi);
      //Fill in totLumi/effLumiForTheHLT
      lumiFactor.push_back(totLumi/2838.119);

      //JetMETTauMonitor/Run2010A-PromptReco-v4/RECO Runs 137437-139558
      totLumi = 61278.016; totalLumis.push_back(totLumi);
      grandTotal+=totLumi;
      files.push_back( TFile::Open( (DIR+"anV3-JetMETTauMonitor-Run2010A_137437-139558_v1-V9_Inc.root").c_str()
				   ) ); lumi.push_back(totLumi);
      lumiFactor.push_back(totLumi/194.437);

      //MinimumBias/Run2010A-Jul16thReReco-v1/RECO Runs 139779-140159
      totLumi = 119107.887; totalLumis.push_back(totLumi);
      grandTotal+=totLumi;
      files.push_back( TFile::Open(  (DIR+"anV3-MinimumBias-Jul16thReReco-Run2010A_139779-140159_v1-v9.root").c_str()
				   ) ); lumi.push_back(totLumi);
      lumiFactor.push_back(totLumi/159.036);

      //MinimumBias/Commissioning10-GOODCOLL-Jun14thSkim_v1
      totLumi = 8025.534; totalLumis.push_back(totLumi);
      grandTotal+=totLumi;
      files.push_back( TFile::Open(  (DIR+"anV3-MinimumBias-Comm10-GoodColl-Jun14thSkim_132440-135735_v1-v9.root").c_str()
				   ) ); lumi.push_back(totLumi);
      lumiFactor.push_back(totLumi/2771.957);

      //MinimumBias/Commissioning10-SD_JetMETTauMonitor-Jun14thSkim_v1
      totLumi = 12949.814; totalLumis.push_back(totLumi);
      grandTotal+=totLumi;
      files.push_back( TFile::Open(  (DIR+"anV3-MinimumBias-JetMETTauMonitor-Jun14th_132440-137028_v1-V9_Inc.root").c_str()
				   ) ); lumi.push_back(totLumi);
      lumiFactor.push_back(totLumi/2838.119);

    }
    
    else{
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
    }  
    float grandScaledTotal=0;
    cout <<"--- "<<HLT << " Total Lumi: "<< grandTotal/1000 << "/nb" << endl;
    for (int l=0;l<(int)lumiFactor.size();l++) {
      cout << "\t Total Lumi: "<< totalLumis[l] << "\t LumiFactor: " <<lumiFactor[l] << "\t File: "<< files[l]->GetName()  << endl;
      grandScaledTotal += totalLumis[l]/lumiFactor[l];
    }
    cout << "--- Effective Lumi: " << grandScaledTotal/1000 << "/nb" << endl;
  }
  
  else if(SET=="Pythia6"){
    //low stat
    if(HLT=="HLT_L1Jet6U"){
    minPthat.push_back(15); maxPthat.push_back(30);
    files.push_back( TFile::Open( (DIR+"anV3-QCD_Pt15_Spring10-V8b.root").c_str() ) ); 
    xsec.push_back(876215000.0); //nEvents.push_back(6090500);
    }
    //if(HLT!="HLT_Jet50U"){//low stat
    minPthat.push_back(30); maxPthat.push_back(80);
     files.push_back( TFile::Open( (DIR+"anV3-QCD_Pt30_Spring10-V8b.root").c_str() ) ); 
     xsec.push_back(60411000); //nEvents.push_back(4989664);
      //}
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

  else if(SET=="Pythia6BB"){

    minPthat.push_back(0); maxPthat.push_back(1000000000);
    files.push_back( TFile::Open( (DIR+"anV3-InclusiveBB_Pt30_Spring10-V8b.root").c_str() ) );
    xsec.push_back(0.069*60000000); //nEvents.push_back(2971800);
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

  //values for Jet70, Jet100 missing, waiting for turnon curves
  if(HLT=="HLT_Jet15U"){ jet1_minPt =56;  jet1_maxPt =84;jet1_minPt_lab ="56";  jet1_maxPt_lab ="84"; }
  else if(HLT=="HLT_Jet30U"){ jet1_minPt =84;  jet1_maxPt =120; jet1_minPt_lab ="84";  jet1_maxPt_lab ="120"; }
  else if(HLT=="HLT_Jet50U"){ jet1_minPt =120;  jet1_maxPt =10000000; jet1_minPt_lab ="120";  jet1_maxPt_lab ="";}
  else if(HLT=="HLT_L1Jet6U"){ jet1_minPt =37;  jet1_maxPt =56; jet1_minPt_lab ="37";  jet1_maxPt_lab ="56";}

  else if(HLT=="HLT_Jet15UOpen"){ jet1_minPt =56;  jet1_maxPt =10000000;jet1_minPt_lab ="56";  jet1_maxPt_lab =""; }
  else if(HLT=="HLT_Jet30UOpen"){ jet1_minPt =84;  jet1_maxPt =10000000; jet1_minPt_lab ="84";  jet1_maxPt_lab =""; }
  else if(HLT=="HLT_Jet50UOpen"){ jet1_minPt =120;  jet1_maxPt =10000000; jet1_minPt_lab ="120";  jet1_maxPt_lab ="";}


//   /////////////////////////////
//   //TEMP: 
//   jet1_maxPt =10000000; 
//   ////////////////////////////////
}


void createCommonHistos(string SET, string HLT){

  string labels[2] = {"Vert","Jet"};
  string labelsC[3] = {"","_CORR","_ErrCORR"};
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
      float lowBin=0, hiBin=2000;
      if(HLT=="HLT_Jet15U") {lowBin=56; hiBin=84;}
      else if(HLT=="HLT_Jet30U") {lowBin=84; hiBin=120;}
      else if(HLT=="HLT_Jet50U") {lowBin=120; hiBin=500;}
      else if(HLT=="HLT_L1Jet6U") {lowBin=37; hiBin=56;}

      histos2D[labels[h]+"_dRvsJet1pt_2b"+labelsC[hC]] = new TH2F((labels[h]+"_dRvsJet1pt_2b"+labelsC[hC]).c_str() ,";#Delta R; Jet1 p_{T}",60,0,6,30,lowBin,hiBin);
      histos2D[labels[h]+"_dRvsJet1pt_2b"+labelsC[hC]]->Sumw2();

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

  histos["Jet_Jet1ptHLT"] = new TH1F("Jet_Jet1ptHLT" ,"",2000,0,2000);  histos["Jet_Jet1ptHLT"]->Sumw2();

  histos["Jet_pt"] = new TH1F("Jet_pt" ,"",200,0,2000);  histos["Jet_pt"]->Sumw2();
  histos["Jet_eta"] = new TH1F("Jet_eta" ,"",100,-3,3);  histos["Jet_eta"]->Sumw2();
  histos["Jet_discr"] = new TH1F("Jet_discr","",200,-10,10); histos["Jet_discr"]->Sumw2();
  
  histos["Jet_pt_2b"] = new TH1F("Jet_pt_2b","",200,0,2000);  histos["Jet_pt_2b"]->Sumw2();
  histos["Jet_eta_2b"] = new TH1F("Jet_eta_2b","",100,-3,3);  histos["Jet_eta_2b"]->Sumw2();
  histos["Jet_ptAsymm_2b"] = new TH1F("Jet_ptAsymm_2b","",100,0,1); histos["Jet_ptAsymm_2b"]->Sumw2();
  histos["Jet_ptAsymm_2b_CORR"] = new TH1F("Jet_ptAsymm_2b_CORR","",100,0,1); histos["Jet_ptAsymm_2b_CORR"]->Sumw2();

  histos["Jet_dR_2b_BTAGCORR"] = new TH1F("Jet_dR_2b_BTAGCORR" ,";MC BHadron #Delta R",60,0,6);  histos["Jet_dR_2b_BTAGCORR"]->Sumw2();
  histos["Jet_dPhi_2b_BTAGCORR"] = new TH1F("Jet_dPhi_2b_BTAGCORR" ,";MC BHadron #Delta #phi",32,0,3.2);  histos["Jet_dPhi_2b_BTAGCORR"]->Sumw2();
  histos["Jet_dEta_2b_BTAGCORR"] = new TH1F("Jet_dEta_2b_BTAGCORR" ,";MC BHadron #Delta #eta",100,0,10);  histos["Jet_dEta_2b_BTAGCORR"]->Sumw2();
  histos["Jet_ptAsymm_2b_BTAGCORR"] = new TH1F("Jet_ptAsymm_2b_BTAGCORR","",100,0,1); histos["Jet_ptAsymm_2b_BTAGCORR"]->Sumw2();

  histos["Jet_bJets_pt"] = new TH1F("Jet_bJets_pt","",200,0,2000); histos["Jet_bJets_pt"]->Sumw2();
  histos["Jet_bJets_eta"] = new TH1F("Jet_bJets_eta","",100,-3,3); histos["Jet_bJets_eta"]->Sumw2();
  histos["Jet_bJets_n"] = new TH1F("Jet_bJets_n","",10,0,10); histos["Jet_bJets_n"]->Sumw2();

  histos["Jet_2bJets_Mass"] = new TH1F("Jet_2bJets_Mass",";mass(bjet,bjet;events)",200,0,2000); histos["Jet_2bJets_Mass"]->Sumw2();

  histos2D["Jet_dEtadPhi_2b"] = new TH2F("Jet_dEtadPhi_2b",";#Delta #eta;#Delta #phi",100,0,10,32,0,3.2); histos2D["Jet_dEtadPhi_2b"]->Sumw2();
  histos2D["Jet_pt1pt2_2b"] = new TH2F("Jet_pt1pt2_2b",";p_{T}^{1};p_{T}^{2}",200,0,2000,200,0,2000); histos2D["Jet_pt1pt2_2b"]->Sumw2();

  histos["Jet_BHadrJetMatched_pt"] = new TH1F("Jet_BHadrJetMatched_pt" ,";MC BHadrJetMatchedon p_{T}",200,0,2000);  histos["Jet_BHadrJetMatched_pt"]->Sumw2();;
  histos["Jet_BHadrJetMatched_eta"]= new TH1F("Jet_BHadrJetMatched_eta" ,";MC BHadrJetMatchedon #eta",100,-3,3);  histos["Jet_BHadrJetMatched_eta"]->Sumw2();
  histos["Jet_BHadrJetMatched_dR"] = new TH1F("Jet_BHadrJetMatched_dR" ,";MC BHadrJetMatchedon #Delta R",60,0,6);  histos["Jet_BHadrJetMatched_dR"]->Sumw2();
  histos["Jet_BHadrJetMatched_dPhi"] = new TH1F("Jet_BHadrJetMatched_dPhi" ,";MC BHadrJetMatchedon #Delta #phi",32,0,3.2);  histos["Jet_BHadrJetMatched_dPhi"]->Sumw2();
  histos["Jet_BHadrJetMatched_dEta"] = new TH1F("Jet_BHadrJetMatched_dEta" ,";MC BHadrJetMatchedon #Delta #eta",100,0,10);  histos["Jet_BHadrJetMatched_dEta"]->Sumw2();

  histos["BHadr_pt"] = new TH1F("BHadr_pt" ,";MC BHadron p_{T}",200,0,2000);  histos["BHadr_pt"]->Sumw2();;
  histos["BHadr_eta"]= new TH1F("BHadr_eta" ,";MC BHadron #eta",100,-3,3);  histos["BHadr_eta"]->Sumw2();
  histos["BHadr_dR"] = new TH1F("BHadr_dR" ,";MC BHadron #Delta R",60,0,6);  histos["BHadr_dR"]->Sumw2();
  histos["BHadr_dPhi"] = new TH1F("BHadr_dPhi" ,";MC BHadron #Delta #phi",32,0,3.2);  histos["BHadr_dPhi"]->Sumw2();
  histos["BHadr_dEta"] = new TH1F("BHadr_dEta" ,";MC BHadron #Delta #eta",100,0,10);  histos["BHadr_dEta"]->Sumw2();

  histos["BHadr_BHadrJetMatched_pt"] = new TH1F("BHadr_BHadrJetMatched_pt" ,";MC BHadrJetMatchedon p_{T}",200,0,2000);  histos["BHadr_BHadrJetMatched_pt"]->Sumw2();;
  histos["BHadr_BHadrJetMatched_eta"]= new TH1F("BHadr_BHadrJetMatched_eta" ,";MC BHadrJetMatchedon #eta",100,-3,3);  histos["BHadr_BHadrJetMatched_eta"]->Sumw2();
  histos["BHadr_BHadrJetMatched_dR"] = new TH1F("BHadr_BHadrJetMatched_dR" ,";MC BHadrJetMatchedon #Delta R",60,0,6);  histos["BHadr_BHadrJetMatched_dR"]->Sumw2();
  histos["BHadr_BHadrJetMatched_dPhi"] = new TH1F("BHadr_BHadrJetMatched_dPhi" ,";MC BHadrJetMatchedon #Delta #phi",32,0,3.2);  histos["BHadr_BHadrJetMatched_dPhi"]->Sumw2();
  histos["BHadr_BHadrJetMatched_dEta"] = new TH1F("BHadr_BHadrJetMatched_dEta" ,";MC BHadrJetMatchedon #Delta #eta",100,0,10);  histos["BHadr_BHadrJetMatched_dEta"]->Sumw2();




  //string fl[] = {"2b","1b1c","1b1l","1c1c","1c1l","1l1l"};
  if(SET!="Data"){
    int lab_fl_n= (int) ( (float)sizeof(lab_fl))/((float) sizeof(string));
    for(int f=0;f<lab_fl_n;f++){
      histos["Jet_ptAsymm_2b_"+lab_fl[f]] = new TH1F( ("Jet_ptAsymm_2b_"+lab_fl[f]).c_str(),"",100,0,1); histos["Jet_ptAsymm_2b_"+lab_fl[f] ]->Sumw2();

      histos["Jet_pt_2b_"+lab_fl[f]] = new TH1F( ("Jet_pt_2b_"+lab_fl[f]).c_str(),"",200,0,2000); histos["Jet_pt_2b_"+lab_fl[f] ]->Sumw2();
      histos["Jet_2bJets_Mass_"+lab_fl[f]] = new TH1F( ("Jet_2bJets_Mass_"+lab_fl[f]).c_str(),";mass(bjet,bjet;events)",200,0,2000); histos["Jet_2bJets_Mass_"+lab_fl[f]]->Sumw2();
      histos2D["Jet_dEtadPhi_2b_"+lab_fl[f]] = new TH2F(("Jet_dEtadPhi_2b_"+lab_fl[f]).c_str(),";#Delta #eta;#Delta #phi",100,0,10,32,0,3.2); histos2D["Jet_dEtadPhi_2b_"+lab_fl[f]]->Sumw2();
      histos2D["Jet_pt1pt2_2b_"+lab_fl[f]] = new TH2F(("Jet_pt1pt2_2b_"+lab_fl[f]).c_str(),";p_{T}^{1};p_{T}^{2}",200,0,2000,200,0,2000); histos2D["Jet_pt1pt2_2b_"+lab_fl[f]]->Sumw2();

    }  
    
    histos["Jet_2bJetsMC_dR"] = new TH1F("Jet_2bJetsMC_dR","",60,0,6); histos["Jet_2bJetsMC_dR"]->Sumw2();
    histos["Jet_2bJetsMC_dPhi"] = new TH1F("Jet_2bJetsMC_dPhi","",32,0,3.2); histos["Jet_2bJetsMC_dPhi"]->Sumw2();
    histos["Jet_2bJetsMC_dEta"] = new TH1F("Jet_2bJetsMC_dEta","",100,0,10); histos["Jet_2bJetsMC_dEta"]->Sumw2();
    
    histos["Jet_2bJetsMC_pt"] = new TH1F("Jet_2bJetsMC_pt","",200,0,2000); histos["Jet_2bJetsMC_pt"]->Sumw2();

    histos["Jet_bJetsMC_pt"] = new TH1F("Jet_bJetsMC_pt","",200,0,2000); histos["Jet_bJetsMC_pt"]->Sumw2();
    histos["Jet_cJetsMC_pt"] = new TH1F("Jet_cJetsMC_pt","",200,0,2000); histos["Jet_cJetsMC_pt"]->Sumw2();
    histos["Jet_lJetsMC_pt"] = new TH1F("Jet_lJetsMC_pt","",200,0,2000); histos["Jet_lJetsMC_pt"]->Sumw2();
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

