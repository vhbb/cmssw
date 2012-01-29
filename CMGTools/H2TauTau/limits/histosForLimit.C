
TH1F* changeBinning(TH1F* hinput, int rebin,float xmax){
  int nbins=(int)(xmax);
  TString hname=hinput->GetName();
  TH1F* houtput=new TH1F("htmp",hinput->GetTitle(),nbins,0,xmax);
  for(Int_t b=1;b<nbins;b++){
    houtput->SetBinContent(b,hinput->GetBinContent(b));
    houtput->SetBinError(b,hinput->GetBinError(b));
  }
  delete hinput;

  houtput->Rebin(rebin);
  houtput->SetName(hname);
  return houtput;
}

std::string fillString(std::string input){
  if(input.size()>15){
    cout<<"input string size too long"<<endl;
    exit(0);
  }
  std::string output=input;
  for(Int_t i=0;i<15;i++)
    if(i>=input.size())output+=" ";
  
  return output;
}

std::string fillStringLong(std::string input){
  if(input.size()>40){
    cout<<"input string size too long"<<endl;
    exit(0);
  }
  std::string output=input;
  for(Int_t i=0;i<40;i++)
    if(i>=input.size())output+=" ";
  
  return output;
}

std::string fillInt(int input){

  char  inchar[100];
  sprintf(inchar,"%d",input);
  std::string newinput(inchar);
  if(newinput.size()>15){
    cout<<"input string size too long"<<endl;
    exit(0);
  }
  std::string output=newinput;
  for(Int_t i=0;i<15;i++)
    if(i>=newinput.size())output+=" ";
  
  return output;
}

std::string fillFloat(float input){

  char  inchar[100];
  sprintf(inchar,"%.3f",input);
  std::string newinput(inchar);
  if(newinput.size()>15){
    cout<<"input string size too long"<<endl;
    exit(0);
  }
  std::string output=newinput;
  for(Int_t i=0;i<15;i++)
    if(i>=newinput.size())output+=" ";
    
 return output;
}

void histosForLimit(long sm, long mass){
  //defines samples
  gROOT->ProcessLine(".x ../workdir/TauMuConfig_limits.C");
  
  //normalizes backgrounds for the data
  if(!analysis.init()){cout<<" could not init"<<endl;return 0;}  
  if(!analysis.scaleHistos())return 0;
  

  /* ***************************************
    -access backgrounds like this (updated Dec 2)
    -everything should be obtained from class interface/TauMuAnalysis.h
    -individual backgrounds obtained using getSample(), except for ZtoTauTau and QCD 
    -methods return ownership of object --> delete objects
    ****************************************** 
    analysis.getTotalDataSS(TString histoname);//sum of SS Data samples 
    analysis.getTotalData(TString histoname);
    analysis.getSample(TString samplename, TString histoname);//can be used to get histo for any MC or Data sample
    analysis.getQCD(TString histoname);
    analysis.getZToTauTau(TString histoname);//Z-->tau tau (either from MC or Embedded)
    analysis.getTotalBackground(TString histoname);
  */


  //histogram used for fit in limit calcualtion
  TString histoname=TString("diTauMassHisto_SM")+sm;//"diTauMassHisto_SM1" "diTauMassHisto_SM2"
  

  TCanvas C("Canvas");
  TString filename=TString("plots_SM")+sm+"_mH"+mass+".ps";
  C.Print(filename+"[");
  
  Int_t rebin=30;
  Int_t xmax=350;
  

  TH1F* htmp = (TH1F*)(analysis.getSample(TString("HiggsVH")+mass,histoname));
  TH1F* VH =(TH1F*)(htmp->Clone(TString("VH")+mass));//rename otherwise default name may conflict below
  VH=changeBinning(VH,rebin,xmax);//use 10GeV wide bins from 0-->300
  delete htmp;
  C.Clear();
  VH->SetTitle(TString("VH")+mass);
  VH->Draw();
  C.Print(filename);
  
 
  TH1F* htmp = (TH1F*)(analysis.getSample(TString("HiggsGG")+mass,histoname));
  TH1F* SM =(TH1F*)(htmp->Clone(TString("SM")+mass));//rename otherwise default name may conflict below
  SM=changeBinning(SM,rebin,xmax);//use 10GeV wide bins from 0-->300
  delete htmp;
  C.Clear();
  SM->SetTitle(TString("SM")+mass);
  SM->Draw();
  C.Print(filename);
 

  TH1F* htmp = (TH1F*)(analysis.getSample(TString("HiggsVBF")+mass,histoname));
  TH1F* VBF =(TH1F*)(htmp->Clone(TString("VBF")+mass));//rename otherwise default name may conflict below
  VBF=changeBinning(VBF,rebin,xmax);
  delete htmp;
  C.Clear();
  VBF->SetTitle(TString("VBF")+mass);
  VBF->Draw();
  C.Print(filename);

  
  TH1F* htmp = (TH1F*)(analysis.getZToTauTau(histoname));
  TH1F* ZTT =(TH1F*)(htmp->Clone("ZTT"));//rename otherwise default name may conflict below
  ZTT=changeBinning(ZTT,rebin,xmax);
  delete htmp;
  C.Clear();
  ZTT->SetTitle("ZTT");
  ZTT->Draw();
  C.Print(filename);

  
  TH1F* htmp = (TH1F*)(analysis.getSample("ZToMuMu",histoname));
  TH1F* ZL =(TH1F*)(htmp->Clone("ZL"));//rename otherwise default name may conflict below
  ZL=changeBinning(ZL,rebin,xmax);
  delete htmp;
  C.Clear();
  ZL->SetTitle("ZL");
  ZL->Draw();
  C.Print(filename);


  TH1F* htmp = (TH1F*)(analysis.getSample("ZToLJet",histoname));
  TH1F* ZJ =(TH1F*)(htmp->Clone("ZJ"));//rename otherwise default name may conflict below
  ZJ=changeBinning(ZJ,rebin,xmax);
  delete htmp;
  C.Clear();
  ZJ->SetTitle("ZJ");
  ZJ->Draw();
  C.Print(filename);

  TH1F* htmp = 0;
  if(sm==0)htmp=(TH1F*)(analysis.getSample("WJetsToLNu",histoname));
  if(sm==1)htmp=(TH1F*)(analysis.getSample("WJetsToLNu",histoname));
  if(sm==2)htmp=(TH1F*)(analysis.getSample("W3JetsToLNu",histoname));
  TH1F* W =(TH1F*)(htmp->Clone("W"));//rename otherwise default name may conflict below
  W=changeBinning(W,rebin,xmax);
  delete htmp;
  C.Clear();
  W->SetTitle("W");
  W->Draw();
  C.Print(filename);

  TH1F* htmp = (TH1F*)(analysis.getSample("TTJets",histoname));
  TH1F* TT =(TH1F*)(htmp->Clone("TT"));//rename otherwise default name may conflict below
  TT=changeBinning(TT,rebin,xmax);
  delete htmp;
  C.Clear();
  TT->SetTitle("TT");
  TT->Draw();
  C.Print(filename);
  
  TH1F* htmp = (TH1F*)(analysis.getSample("WW",histoname));
  TH1F* WW =(TH1F*)(htmp->Clone("WW"));//rename otherwise default name may conflict below
  WW=changeBinning(WW,rebin,xmax);
  delete htmp;
  C.Clear();
  WW->SetTitle("WW");
  WW->Draw();
  C.Print(filename);
  
  TH1F* htmp = (TH1F*)(analysis.getSample("WZ",histoname));
  TH1F* WZ =(TH1F*)(htmp->Clone("WZ"));//rename otherwise default name may conflict below
  WZ=changeBinning(WZ,rebin,xmax);
  delete htmp;
  C.Clear();
  WZ->SetTitle("WZ");
  WZ->Draw();
  C.Print(filename);
  
  TH1F* htmp = (TH1F*)(analysis.getSample("ZZ",histoname));
  TH1F* ZZ =(TH1F*)(htmp->Clone("ZZ"));//rename otherwise default name may conflict below
  ZZ=changeBinning(ZZ,rebin,xmax);
  delete htmp;
  C.Clear();
  ZZ->SetTitle("ZZ");
  ZZ->Draw();
  C.Print(filename);
  
  TH1F* VV =(TH1F*)(WW->Clone("VV"));
  VV->Add(WZ);
  VV->Add(ZZ);
  
  TH1F* htmp=0;
  if(sm==0)htmp=analysis.getQCD(histoname);
  if(sm==1)htmp=analysis.getSMQCD(histoname);
  if(sm==2)htmp=analysis.getSMQCD2(histoname);
  TH1F* QCD=(TH1F*)(htmp->Clone("QCD"));//rename otherwise default name may conflict below
  QCD=changeBinning(QCD,rebin,xmax);
  delete htmp;
  C.Clear();
  QCD->SetTitle("QCD");
  QCD->Draw();
  C.Print(filename);

  
  TH1F* htmp=0;
  if(sm==0)htmp=analysis.getTotalBackground(histoname,0,0);
  if(sm==1)htmp=analysis.getTotalBackground(histoname,1,1);
  if(sm==2)htmp=analysis.getTotalBackground(histoname,2,2);
  TH1F* bkg_obs = (TH1F*)(htmp->Clone("bkg_obs"));//rename otherwise default name may conflict below
  bkg_obs=changeBinning(bkg_obs,rebin,xmax);
  delete htmp;
  
  TH1F* htmp=analysis.getTotalData(histoname);
  TH1F* data_obs = (TH1F*)(htmp->Clone("data_obs"));//rename otherwise default name may conflict below
  data_obs=changeBinning(data_obs,rebin,xmax);
  delete htmp;


  C.Clear();
  TH1F*hSig=(TH1F*)(VH->Clone("hSig"));
  hSig->Add(SM);
  hSig->Add(VBF);
  hSig->Scale(5);
  hSig->Add(bkg_obs);
  hSig->SetLineColor(2);
  hSig->Draw("hist");
  data_obs->Draw("pesame");
  bkg_obs->Draw("histsame");
  C.Print(filename);

  cout<<" Data yield = "<<data_obs->Integral()<<"  Background yield = "<<bkg_obs->Integral()<<"  5xSM="<<hSig->Integral()<<endl;

  delete hSig;

  TH1F* htmp=analysis.getTotalDataSS(histoname);
  TH1F* data_obsSS=(TH1F*)(htmp->Clone("data_obsSS"));//rename otherwise default name may conflict below
  data_obsSS=changeBinning(data_obsSS,rebin,xmax);
  delete htmp;
  C.Clear();
  data_obsSS->SetTitle("data_obsSS");
  data_obsSS->Draw();
  C.Print(filename);

  C.Print(filename+"]");


  ////////////////Histograms file

  char rootfilename[256];
  sprintf(rootfilename,"muTau_SM%d_mH%d.root",sm,mass);  
  TFile output(rootfilename,"recreate");
  TDirectory* dir = output.mkdir(TString("muTau_SM")+sm);  
  dir->cd();

  VH->Write();
  SM->Write();
  VBF->Write();
  ZTT ->Write();
  ZL->Write();
  ZJ->Write();
  W->Write();
  TT->Write();
  VV->Write();
  QCD->Write();
  data_obs->Write();

  output.Close();


  
  //////////////////////////////////////Data Card
  ofstream file;
  file.open((const char*)(TString("muTau_SM")+sm+"_mH"+mass+".txt"));

  file << "imax 1" <<endl;
  file << "jmax *" <<endl;
  file << "kmax *" <<endl;
  file << "shapes *  *  "<< rootfilename << "  $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC" << endl;
  file << endl;
  file << "observation " << data_obs->Integral()<<endl;
  file<<endl;


  file <<fillStringLong("bin")<<fillString(" ")
       <<fillString(" ")
       <<fillString((const char*)(TString("muTau_SM")+sm))
       <<fillString((const char*)(TString("muTau_SM")+sm))
       <<fillString((const char*)(TString("muTau_SM")+sm))
       <<fillString((const char*)(TString("muTau_SM")+sm))
       <<fillString((const char*)(TString("muTau_SM")+sm))
       <<fillString((const char*)(TString("muTau_SM")+sm))
       <<fillString((const char*)(TString("muTau_SM")+sm))
       <<fillString((const char*)(TString("muTau_SM")+sm))
       <<fillString((const char*)(TString("muTau_SM")+sm))
       <<fillString((const char*)(TString("muTau_SM")+sm))<<endl;
  file <<fillStringLong("process")<<fillString(" ")
       <<fillString(" ")
       <<fillString((const char*)(TString("VH")+mass))
       <<fillString((const char*)(TString("SM")+mass))
       <<fillString((const char*)(TString("VBF")+mass))
       <<fillString("ZTT")
       <<fillString("QCD")
       <<fillString("W")
       <<fillString("ZJ")
       <<fillString("ZL")
       <<fillString("TT")
       <<fillString("VV")<<endl;
  file <<fillStringLong("process")<<fillString(" ")
       <<fillString(" ")
       <<fillInt(-2)
       <<fillInt(-1)
       <<fillInt(0)
       <<fillInt(1)
       <<fillInt(2)
       <<fillInt(3)
       <<fillInt(4)
       <<fillInt(5)
       <<fillInt(6)
       <<fillInt(7)<<endl;
  file <<fillStringLong("rate")<<fillString(" ")
       <<fillString(" ")
       <<fillFloat(VH->Integral())
       <<fillFloat(SM->Integral())
       <<fillFloat(VBF->Integral())
       <<fillFloat(ZTT->Integral())
       <<fillFloat(QCD->Integral())
       <<fillFloat(W->Integral())
       <<fillFloat(ZJ->Integral())
       <<fillFloat(ZL->Integral())
       <<fillFloat(TT->Integral())
       <<fillFloat(VV->Integral())<<endl;

  file << endl;
  file <<"-------------------------------"<<endl;


  file <<fillStringLong("lumi")<<fillString("lnN")
       <<fillString(" ")
       <<fillFloat(1.045)
       <<fillFloat(1.045)
       <<fillFloat(1.045)
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<"luminosity"<<endl;

  file <<fillStringLong("CMS_eff_m")<<fillString("lnN")
       <<fillString(" ")
       <<fillFloat(1.020)
       <<fillFloat(1.020)
       <<fillFloat(1.020)
       <<fillFloat(1.020)
       <<fillString("-")
       <<fillString("-")
       <<fillFloat(1.020)
       <<fillFloat(1.020)
       <<fillFloat(1.020)
       <<fillFloat(1.020)
       <<"Muon ID/HLT"<<endl;

  file <<fillStringLong("CMS_eff_t")<<fillString("lnN")
       <<fillString(" ")
       <<fillFloat(1.06)
       <<fillFloat(1.06)
       <<fillFloat(1.06)
       <<fillFloat(1.06)
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillFloat(1.06)
       <<fillFloat(1.06)
       <<"Tau ID"<<endl;

  file <<fillStringLong("CMS_htt_zttNorm")<<fillString("lnN")
       <<fillString(" ")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillFloat(1.033)
       <<fillString("-")
       <<fillString("-")
       <<fillFloat(1.033)
       <<fillFloat(1.033)
       <<fillString("-")
       <<fillString("-")
       <<"Norm Ztt"<<endl;

  file <<fillStringLong("CMS_htt_ttbarNorm")<<fillString("lnN")
       <<fillString(" ")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillFloat(1.075)
       <<fillString("-")
       <<"Norm ttbar"<<endl;

  file <<fillStringLong("CMS_htt_DiBosonNorm")<<fillString("lnN")
       <<fillString(" ")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillFloat(1.30)
       <<"Norm DiBoson"<<endl;

  file <<fillStringLong("CMS_htt_WNorm")<<fillString("lnN")
       <<fillString(" ")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillFloat(1.066)
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<"Norm W+Jets"<<endl;

//  file <<fillStringLong("CMS_htt_muTau_SM0_WNorm")<<fillString("gmN")
file <<fillStringLong("CMS_htt_muTau_WNorm")<<fillString("gmN")
       <<fillInt(W->Integral()/analysis.getWJetsSignalToSBFraction())
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillFloat(analysis.getWJetsSignalToSBFraction())
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<"W Background"<<endl;

//file <<fillStringLong("CMS_htt_muTau_SM0_QCDNorm")<<fillString("gmN")
file <<fillStringLong("CMS_htt_muTau_QCDNorm")<<fillString("gmN")
       <<fillInt(QCD->Integral()/analysis.getQCDOStoSSRatio())
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillFloat(analysis.getQCDOStoSSRatio())
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<"QCD Background"<<endl;

  file <<fillStringLong("CMS_htt_ZJetFakeTau")<<fillString("lnN")
       <<fillString(" ")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillFloat(1.10)
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<"Z(jet->tau) background"<<endl;

  file <<fillStringLong("CMS_htt_ZLeptonFakeTau")<<fillString("lnN")
       <<fillString(" ")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillFloat(1.10)
       <<fillString("-")
       <<fillString("-")
       <<"Z(lepton->tau) background"<<endl;

  file <<fillStringLong("CMS_scale_j")<<fillString("lnN")
       <<fillString(" ")
       <<fillFloat(0.96)
       <<fillFloat(0.99)
       <<fillFloat(0.92)
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillFloat(0.97)
       <<fillFloat(0.94)
       <<"Jet scale"<<endl;

  file <<fillStringLong("CMS_scale_met")<<fillString("lnN")
       <<fillString(" ")
       <<fillFloat(1.05)
       <<fillFloat(1.05)
       <<fillFloat(1.05)
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillFloat(1.05)
       <<fillFloat(1.05)
       <<fillFloat(1.02)
       <<fillFloat(1.06)
       <<"Jet scale"<<endl;

  file <<fillStringLong("pdf_vh")<<fillString("lnN")
       <<fillString(" ")
       <<fillFloat(1.08)
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<"PDF VBF"<<endl;

  file <<fillStringLong("pdf_gg")<<fillString("lnN")
       <<fillString(" ")
       <<fillString("-")
       <<fillFloat(1.03)
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<"PDF VBF"<<endl;

  file <<fillStringLong("pdf_qqbar")<<fillString("lnN")
       <<fillString(" ")
       <<fillString("-")
       <<fillString("-")
       <<fillFloat(1.03)
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<"PDF VBF"<<endl;


  file <<fillStringLong("QCDscale_ggH")<<fillString("lnN")
       <<fillString(" ")
       <<fillString("-")
       <<fillFloat(1.12)
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<"QCD scale GGF"<<endl;

  file <<fillStringLong("QCDscale_qqH")<<fillString("lnN")
       <<fillString(" ")
       <<fillString("-")
       <<fillString("-")
       <<fillFloat(1.035)
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<"QCD scale VBF"<<endl;


  file <<fillStringLong("UEPS")<<fillString("lnN")
       <<fillString(" ")
       <<fillFloat(0.96)
       <<fillFloat(0.96)
       <<fillFloat(0.96)
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<fillString("-")
       <<"Underlying events & PS"<<endl;

  file.close();
  
  gROOT->ProcessLine(".q");
}


//   TH1F* SM115_CMS_scale_tUp = (TH1F*)SM115->Clone("SM115_CMS_scale_tUp");
//   TH1F* SM115_CMS_scale_tDown = (TH1F*)SM115->Clone("SM115_CMS_scale_tDown");
//   TH1F* SM115_CMS_scale_jUp = (TH1F*)SM115->Clone("SM115_CMS_scale_jUp");
//   TH1F* SM115_CMS_scale_jDown = (TH1F*)SM115->Clone("SM115_CMS_scale_jDown");
  
//   TH1F* VBF115_CMS_scale_tUp = (TH1F*)VBF115->Clone("VBF115_CMS_scale_tUp");
//   TH1F* VBF115_CMS_scale_tDown = (TH1F*)VBF115->Clone("VBF115_CMS_scale_tDown");
//   TH1F* VBF115_CMS_scale_jUp = (TH1F*)VBF115->Clone("VBF115_CMS_scale_jUp");
//   TH1F* VBF115_CMS_scale_jDown = (TH1F*)VBF115->Clone("VBF115_CMS_scale_jDown");  

//   TH1F* ZTT_CMS_scale_tUp = (TH1F*)ZTT->Clone("ZTT_CMS_scale_tUp");
//   TH1F* ZTT_CMS_scale_tDown = (TH1F*)ZTT->Clone("ZTT_CMS_scale_tDown");
//   TH1F* ZTT_CMS_scale_jUp = (TH1F*)ZTT->Clone("ZTT_CMS_scale_jUp");
//   TH1F* ZTT_CMS_scale_jDown = (TH1F*)ZTT->Clone("ZTT_CMS_scale_jDown");
  
//   TH1F* ZL_CMS_scale_tUp = (TH1F*)ZL->Clone("ZL_CMS_scale_tUp");
//   TH1F* ZL_CMS_scale_tDown = (TH1F*)ZL->Clone("ZL_CMS_scale_tDown");
//   TH1F* ZL_CMS_scale_jUp = (TH1F*)ZL->Clone("ZL_CMS_scale_jUp");
//   TH1F* ZL_CMS_scale_jDown = (TH1F*)ZL->Clone("ZL_CMS_scale_jDown");
  
//   TH1F* ZJ_CMS_scale_tUp = (TH1F*)ZJ->Clone("ZJ_CMS_scale_tUp");
//   TH1F* ZJ_CMS_scale_tDown = (TH1F*)ZJ->Clone("ZJ_CMS_scale_tDown");
//   TH1F* ZJ_CMS_scale_jUp = (TH1F*)ZJ->Clone("ZJ_CMS_scale_jUp");
//   TH1F* ZJ_CMS_scale_jDown = (TH1F*)ZJ->Clone("ZJ_CMS_scale_jDown");

//   TH1F* W_CMS_scale_tUp = (TH1F*)W->Clone("W_CMS_scale_tUp");
//   TH1F* W_CMS_scale_tDown = (TH1F*)W->Clone("W_CMS_scale_tDown");
//   TH1F* W_CMS_scale_jUp = (TH1F*)W->Clone("W_CMS_scale_jUp");
//   TH1F* W_CMS_scale_jDown = (TH1F*)W->Clone("W_CMS_scale_jDown");

//   TH1F* TT_CMS_scale_tUp = (TH1F*)TT->Clone("TT_CMS_scale_tUp");
//   TH1F* TT_CMS_scale_tDown = (TH1F*)TT->Clone("TT_CMS_scale_tDown");
//   TH1F* TT_CMS_scale_jUp = (TH1F*)TT->Clone("TT_CMS_scale_jUp");
//   TH1F* TT_CMS_scale_jDown = (TH1F*)TT->Clone("TT_CMS_scale_jDown");
    
//   SM115_CMS_scale_tUp   ->Write();
//   SM115_CMS_scale_tDown ->Write();
//   SM115_CMS_scale_jUp   ->Write();
//   SM115_CMS_scale_jDown ->Write();

//   VBF115_CMS_scale_tUp   ->Write();
//   VBF115_CMS_scale_tDown ->Write();
//   VBF115_CMS_scale_jUp   ->Write();
//   VBF115_CMS_scale_jDown ->Write();

//   ZTT_CMS_scale_tUp   ->Write();
//   ZTT_CMS_scale_tDown ->Write();
//   ZTT_CMS_scale_jUp   ->Write();
//   ZTT_CMS_scale_jDown ->Write();

//   ZL_CMS_scale_tUp   ->Write();
//   ZL_CMS_scale_tDown ->Write();
//   ZL_CMS_scale_jUp   ->Write();
//   ZL_CMS_scale_jDown ->Write();

//   ZJ_CMS_scale_tUp   ->Write();
//   ZJ_CMS_scale_tDown ->Write();
//   ZJ_CMS_scale_jUp   ->Write();
//   ZJ_CMS_scale_jDown ->Write();

//   W_CMS_scale_tUp   ->Write();
//   W_CMS_scale_tDown ->Write();
//   W_CMS_scale_jUp   ->Write();
//   W_CMS_scale_jDown ->Write();
//   TT_CMS_scale_tUp   ->Write();
//   TT_CMS_scale_tDown ->Write();
//   TT_CMS_scale_jUp   ->Write();
//   TT_CMS_scale_jDown ->Write();
