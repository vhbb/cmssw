// @(#)root/tmva $Id: TMVAClassification.cxx 37399 2010-12-08 15:22:07Z evt $
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassification                                                 *
 *                                                                                *
 * This executable provides examples for the training and testing of the          *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans.        *
 *                                                                                *
 * Compile and run the example with the following commands                        *
 *                                                                                *
 *    make                                                                        *
 *    ./TMVAClassification <Methods>                                              *
 *                                                                                *
 * where: <Methods> = "method1 method2"                                           *
 *        are the TMVA classifier names                                           *
 *                                                                                *
 * example:                                                                       *
 *    ./TMVAClassification Fisher LikelihoodPCA BDT                               *
 *                                                                                *
 * If no method given, a default set is of classifiers is used                    *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <../macros/macro.C>), which can be conveniently    *
 * invoked through a GUI launched by the command                                  *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 **********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Config.h"

#include "../macros/samples.h"

// read input data file with ascii format (otherwise ROOT) ?
Bool_t ReadDataFromAsciiIFormat = kFALSE;

int main( int argc, char** argv )
{
   //---------------------------------------------------------------
   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 1;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // 
   // --- 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 1;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // --- Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Support Vector Machine 
   Use["SVM"]             = 0;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------

   std::cout << std::endl << "==> Start TMVAClassification" << std::endl;

   bool batchMode(false);
   bool useDefaultMethods(true);

//    // Select methods (don't look at this code - not of interest)
//    for (int i=1; i<argc; i++) {
//       std::string regMethod(argv[i]);
//       if(regMethod=="-b" || regMethod=="--batch") {
//          batchMode=true;
//          continue;
//       }
//       if (Use.find(regMethod) == Use.end()) {
//          std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
//          for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
//          std::cout << std::endl;
//          return 1;
//       }
//       useDefaultMethods = false;
//    }

//    if (!useDefaultMethods) {
//       for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
//       for (int i=1; i<argc; i++) {
//          std::string regMethod(argv[i]);
//          if(regMethod=="-b" || regMethod=="--batch") continue;
//          Use[regMethod] = 1;
//       }
//    }

   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins
   //   string channel="Zmm";

//    std::cout << "argc = " << argc << std::endl;
   std::string channel;
   if( argc > 1 )
     channel=argv[1];
   else{
     std::cerr << "no channel selected" << std::endl;
     return -1;
   }
   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA.root" );
   outfileName=channel+outfileName;
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is 
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
   (TMVA::gConfig().GetIONames()).fWeightFileExtension = channel+"_weight";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

   float BR115 = 0.704;
   float BR120 = 0.648;
   float BR125 = 0.577;
   float BR130 = 0.493;
   float BR135 = 0.403;

   // 115
   float xSecWH = 0.7546*BR115;
   float xSecZH = 0.4107*BR115;
   float backgroundXsecWJets =  31314;
   float backgroundXsecDY = 3048;
   float backgroundXsecWZ = 18.3;     
   float backgroundXsecWW = 42.9;     
   float backgroundXsecZZ = 5.9;
   float backgroundXsecTT = 165;
   float backgroundXsecT_tchannel = 41.92;
   float backgroundXsecT_tWDRchannel = 7.87;
   float backgroundXsecTbar_tchannel = 22.65;
   float backgroundXsecTbar_tWDRchannel = 7.87;
//    float backgroundXsecTbar_tWDSchannel = 25;

   //get Lumi from something serius
   float lumi = 2047.0+8872.0;
   
   float signalXsec;
   std::vector<double> backgroundXsec;
   std::vector<TFile*> myInputFile_signal;
   std::vector<TFile*> myInputFile_background;

   TCut mycuts = "";
   TCut mycutb = "";
   bool origCuts = true;


   //signal for Z
   if( channel == "Zmm" or channel == "Zee" or channel == "Znn" ){
     myInputFile_signal.push_back(TFile::Open("MC_files/file_ZH_ZToLL_HToBB_M-115_7TeV-powheg-herwigppTreeFile.root"));
     signalXsec = xSecZH;
   }
   //signal for W
   if( channel == "Wm" or channel == "We" ){
     myInputFile_signal.push_back(TFile::Open("MC_files/file_ZH_ZToLL_HToBB_M-115_7TeV-powheg-herwigppTreeFile.root"));
     signalXsec = xSecWH;
   }
   //background
   myInputFile_background.push_back(TFile::Open("MC_files/file_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauolaTreeFile.root"));
   backgroundXsec.push_back(backgroundXsecDY);
   myInputFile_background.push_back(TFile::Open("MC_files//file_WJetsToLNu_TuneZ2_7TeV-madgraph-tauolaTreeFile.root"));
   backgroundXsec.push_back(backgroundXsecWJets);
   myInputFile_background.push_back(TFile::Open("MC_files//file_TTJets_TuneZ2_7TeV-madgraph-tauolaTreeFile.root"));
   backgroundXsec.push_back(backgroundXsecTT);
   myInputFile_background.push_back(TFile::Open("MC_files//file_T_TuneZ2_t-channel_7TeV-powheg-tauolaTreeFile.root"));
   backgroundXsec.push_back(backgroundXsecT_tchannel);
   myInputFile_background.push_back(TFile::Open("MC_files//file_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauolaTreeFile.root"));
   backgroundXsec.push_back(backgroundXsecT_tWDRchannel);
   myInputFile_background.push_back(TFile::Open("MC_files//file_Tbar_TuneZ2_t-channel_7TeV-powheg-tauolaTreeFile.root"));
   backgroundXsec.push_back(backgroundXsecTbar_tchannel);
   myInputFile_background.push_back(TFile::Open("MC_files//file_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauolaTreeFile.root"));
   backgroundXsec.push_back(backgroundXsecTbar_tWDRchannel);
   //      myInputFile_background.push_back(TFile::Open("MC_files//file_Tbar_TuneZ2_tW-channel-DS_7TeV-powheg-tauolaTreeFile.root"));
   //      backgroundXsec.push_back(backgroundXsecTbar_tWDSchannel);
   myInputFile_background.push_back(TFile::Open("MC_files///file_ZZ_TuneZ2_7TeV_pythia6_tauolaTreeFile.root"));
   backgroundXsec.push_back(backgroundXsecZZ);
   myInputFile_background.push_back(TFile::Open("MC_files/file_WW_TuneZ2_7TeV_pythia6_tauolaTreeFile.root"));
   backgroundXsec.push_back(backgroundXsecWW);
   //WZ still missing
   //      myInputFile_background.push_back(TFile::Open("MC_files/"));
   //      backgroundXsec.push_back(backgroundXsecWZ);
   myInputFile_background.push_back(TFile::Open("MC_files//file_ZZ_TuneZ2_7TeV_pythia6_tauolaTreeFile.root"));
   backgroundXsec.push_back(backgroundXsecZZ);


   if(channel == "Zmm"){
     factory->AddVariable( "bbMass", "bb mass"      , "GeV/c2", 'F' );
     factory->AddVariable( "VMass", "vector mass" , "GeV/c2", 'F' );
     factory->AddVariable( "bbPt", "bb pt"        , "GeV/c" , 'F' );
     factory->AddVariable( "VPt", "vector pt"  , "GeV/c" , 'F' );
     factory->AddVariable( "btag1", "btag1"        , "csv"   , 'F' );
     factory->AddVariable( "btag2", "btag2"        , "csv"   , 'F' );
     factory->AddVariable( "DeltaPhiVH", "DeltaPhi(V,H)"    , ""  , 'F' );
     factory->AddVariable( "DeltaEtabb", "DeltaEta(b,b)"     , "" , 'F' );

     //PRESELECTION
     mycuts="(bPt1>20) && (bPt2>20) && (btag1>0.5) && (btag2>0.5) && (NaddJet<2) && (DeltaPhiVH>2.4)";
     mycutb="(bPt1>20) && (bPt2>20) && (btag1>0.5) && (btag2>0.5) && (NaddJet<2) && (DeltaPhiVH>2.4)";

   }
   else if(channel == "Zee"){
     factory->AddVariable( "bbMass", "bb mass"      , "GeV/c2", 'F' );
     factory->AddVariable( "VMass", "vector mass" , "GeV/c2", 'F' );
     factory->AddVariable( "bbPt", "bb pt"        , "GeV/c" , 'F' );
     factory->AddVariable( "VPt", "vector pt"  , "GeV/c" , 'F' );
     factory->AddVariable( "btag1", "btag1"        , "csv"   , 'F' );
     factory->AddVariable( "btag2", "btag2"        , "csv"   , 'F' );
     factory->AddVariable( "DeltaPhiVH", "DeltaPhi(V,H)"    , ""  , 'F' );
     factory->AddVariable( "DeltaEtabb", "DeltaEta(b,b)"     , "" , 'F' );

     //PRESELECTION
     mycuts="bPt1>20 && bPt2>20 && btag1>0.5 && btag2>0.5 && NaddJet<2 && DeltaPhiVH>2.4";
     mycutb="bPt1>20 && bPt2>20 && btag1>0.5 && btag2>0.5 && NaddJet<2 && DeltaPhiVH>2.4";

   }
   else if(channel == "Znn"){
     factory->AddVariable( "bbMass", "bb mass"      , "GeV/c2", 'F' );
     factory->AddVariable( "bbPt", "bb pt"        , "GeV/c" , 'F' );
     factory->AddVariable( "pfMET", "met pt"  , "GeV/c" , 'F' );
     factory->AddVariable( "btag1", "btag1"        , "csv"   , 'F' );
     factory->AddVariable( "btag2", "btag2"        , "csv"   , 'F' );
     factory->AddVariable( "DeltaPhiVH", "DeltaPhi(V,H)"    , ""  , 'F' );
     //NaddJet is imposed to be 0
     //     factory->AddVariable( "NaddJet", "NaddJet"      , ""      , 'F' );

     //PRESELECTION
     mycuts="bPt1>80 && bPt2>20 && bbPt>160 && btag1>0.5 && btag2>0.5 && NaddJet<1 && deltaPhipfMETjet1>0.5 && deltaPhipfMETjet2 && pfMETsig>5";
     mycutb="bPt1>80 && bPt2>20 && bbPt>160 && btag1>0.5 && btag2>0.5 && NaddJet<1 && deltaPhipfMETjet1>0.5 && deltaPhipfMETjet2 && pfMETsig>5";

   } 
   else if(channel == "We"){
     factory->AddVariable( "bbMass", "bb mass"      , "GeV/c2", 'F' );
     factory->AddVariable( "bbPt", "bb pt"        , "GeV/c" , 'F' );
     factory->AddVariable( "VPt", "vector pt"  , "GeV/c" , 'F' );
     factory->AddVariable( "btag1", "btag1"        , "csv"   , 'F' );
     factory->AddVariable( "btag2", "btag2"        , "csv"   , 'F' );
     factory->AddVariable( "DeltaPhiVH", "DeltaPhi(V,H)"    , ""  , 'F' );
     factory->AddVariable( "DeltaEtabb", "DeltaEta(b,b)"     , "" , 'F' );
     //NaddJet is imposed to be 0
     //    factory->AddVariable( "NaddJet", "NaddJet"      , ""      , 'F' );

     //PRESELECTION
     mycuts="bPt1>30 && bPt2>30 && bbPt>150 && VPt>150 && btag1>0.4 && btag2>0.4 && NaddJet<1 && pfMETsig>2";
     mycutb="bPt1>30 && bPt2>30 && bbPt>150 && VPt>150 && btag1>0.4 && btag2>0.4 && NaddJet<1 && pfMETsig>2";
     
   }
   else if(channel == "Wm"){
     factory->AddVariable( "bbMass", "bb mass"      , "GeV/c2", 'F' );
     factory->AddVariable( "bbPt", "bb pt"        , "GeV/c" , 'F' );
     factory->AddVariable( "VPt", "vector pt"  , "GeV/c" , 'F' );
     factory->AddVariable( "btag1", "btag1"        , "csv"   , 'F' );
     factory->AddVariable( "btag2", "btag2"        , "csv"   , 'F' );
     factory->AddVariable( "DeltaPhiVH", "DeltaPhi(V,H)"    , ""  , 'F' );
     factory->AddVariable( "DeltaEtabb", "DeltaEta(b,b)"     , "" , 'F' );
     //NaddJet is imposed to be 0
     //     factory->AddVariable( "NaddJet", "NaddJet"      , ""      , 'F' );

     //PRESELECTION
     mycuts="bPt1>30 && bPt2>30 && bbPt>150 && VPt>150 && btag1>0.4 && btag2>0.4 && NaddJet<1";
     mycutb="bPt1>30 && bPt2>30 && bbPt>150 && VPt>150 && btag1>0.4 && btag2>0.4 && NaddJet<1";
     
   }

   std::vector<Sample> samples_signal;
   std::vector<Sample> samples_background;

   std::string str_toReplace("TreeFile");

   std::vector<TTree*> myInputTree_signal;
   std::vector<TTree*> myInputTree_backgound;
   for(unsigned int i = 0; i < myInputFile_signal.size(); ++i ){
     std::string name = myInputFile_signal.at(0)->GetName();
     samples_signal.push_back(Sample(signalXsec,"signal", name.replace(name.find("TreeFile"),str_toReplace.length() ,"_histos" ) , kRed,false));
     TTree * tmpTree = (TTree*)myInputFile_signal.at(i)->Get("treeMVA");
     factory->AddSignalTree( tmpTree , samples_signal.at(i).scale(lumi) );
   }
   for(unsigned int i = 0; i < myInputFile_background.size(); ++i ){
     std::string name = myInputFile_background.at(i)->GetName();
     samples_background.push_back(Sample(backgroundXsec.at(i),"background", name.replace(name.find("TreeFile"),str_toReplace.length() ,"_histos" ) ,i+10,false));
     TTree * tmpbackTree = (TTree*)myInputFile_background.at(i)->Get("treeMVA");
     factory->AddBackgroundTree( tmpbackTree,  samples_background.at(i).scale(lumi) );
   }

   //from Michele and Matt
   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
					"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );


   // ---- Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod( TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   if (Use["CutsD"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsPCA",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

   if (Use["CutsSA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   // Likelihood ("naive Bayes estimator")
   if (Use["Likelihood"])
      factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
                           "H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

   // Decorrelated likelihood
   if (Use["LikelihoodD"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD",
                           "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

   // PCA-transformed likelihood
   if (Use["LikelihoodPCA"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA",
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 

   // Use a kernel density estimator to approximate the PDFs
   if (Use["LikelihoodKDE"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE",
                           "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

   // Use a variable-dependent mix of splines and kernel density estimator
   if (Use["LikelihoodMIX"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX",
                           "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

   // Test the multi-dimensional probability density estimator
   // here are the options strings for the MinMax and RMS methods, respectively:
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
   if (Use["PDERS"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERS",
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

   if (Use["PDERSD"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSD",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

   if (Use["PDERSPCA"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSPCA",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

   // Multi-dimensional likelihood estimator using self-adapting phase-space binning
   if (Use["PDEFoam"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam",
                           "H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0333:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

   if (Use["PDEFoamBoost"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoamBoost",
                           "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod( TMVA::Types::kKNN, "KNN",
                           "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   // H-Matrix (chi2-squared) method
   if (Use["HMatrix"])
      factory->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V" );

   // Linear discriminant (same as Fisher discriminant)
   if (Use["LD"])
      factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher discriminant (same as LD)
   if (Use["Fisher"])
      factory->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher with Gauss-transformed input variables
   if (Use["FisherG"])
      factory->BookMethod( TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

   // Composite classifier: ensemble (tree) of boosted Fisher classifiers
   if (Use["BoostedFisher"])
      factory->BookMethod( TMVA::Types::kFisher, "BoostedFisher", "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2");

   // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

   if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_SA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   if (Use["FDA_MT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   if (Use["FDA_GAMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   if (Use["FDA_MCMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MCMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"])
      factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

   if (Use["MLPBFGS"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   if (Use["MLPBNN"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

   // CF(Clermont-Ferrand)ANN
   if (Use["CFMlpANN"])
      factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

   // Tmlp(Root)ANN
   if (Use["TMlpANN"])
      factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

   // Support Vector Machine
   if (Use["SVM"])
      factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

   if (Use["BDTB"]) // Bagging
      factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

   // For an example of the category classifier, see: TMVAClassificationCategory

   // For an example of the category classifier usage, see: TMVAClassificationCategory

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","GA");

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl
             << "==> TMVAClassification is done!" << std::endl
             << std::endl
             << "==> To view the results, launch the GUI: \"root -l ./TMVAGui.C\"" << std::endl
             << std::endl;

   // Clean up
   delete factory;
}

