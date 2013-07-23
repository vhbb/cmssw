void controlPlots() {
  gStyle->SetOptStat(0); 
  gStyle->SetEndErrorSize(3); 
//   TStyle *dataST = gStyle->Clone(); 
//   TStyle *mcST = gStyle->Clone(); 
//   dataST->SetErrorX(0); 

  //systematic errors
  double systerror_effCorr = 0.13; 
  double systerror_various = 0.113; 
  double systerror = sqrt(systerror_effCorr*systerror_effCorr+systerror_various*systerror_various); 

  //old
//   double systerror = 0.113; 


  double maxr = -1.0, minr =1000.0;
  double maxa = -1.0, mina =1000.0;



  string dir="/scratch/wehrlilu/BBCORR/MACROS/SUMFILES/";
//   string dirPYTH="tempBDCOMPplots/";
//   string dirMADG="tempJESminus/";
//   string dirMCAT="tempJESplus/";
  string dirPYTH="MC/";
  string dirDATA="DATA/";
  
  string mcfilename1 = "histo-Pythia6-HLT_";
  string mcfilename2 = "-BBCorr_scaledEff-total_ETACUTOPEN.root"; 
  string datafilename1 = "histo-Data-HLT_";
  string datafilename2 = "-BBCorr_scaledEff-total_ETACUTOPEN.root"; 

  //Jet50U

  string strigger; 
//   string histoname="Vert_dR_efficden";



  std::vector<double> vratio; 
  std::vector<double> vasym; 
  std::vector<double> vratioDATA; 
  std::vector<double> vasymDATA;

  std::vector<double> vratioErr; 
  std::vector<double> vasymErr; 
  std::vector<double> vratioDATAErr; 
  std::vector<double> vasymDATAErr;
  std::vector<double> vratioDATAErr_stat; 
  std::vector<double> vasymDATAErr_stat;

  for(unsigned int t=0; t<3; t++){
//   for(unsigned int t=0; t<1; t++){
    switch(t){
    case 0: strigger="Jet15U"; break;
    case 1: strigger="Jet30U"; break;
    case 2: strigger="Jet50U"; break;

    }

    TFile *fePYTH = TFile::Open((dir+dirPYTH+mcfilename1+strigger+datafilename2).c_str());
    //data with pythia correction
    TFile *feDATA = TFile::Open((dir+dirDATA+datafilename1+strigger+datafilename2).c_str());


    TH1F *hePYTH; fePYTH->GetObject("Vert_dR_efficden",hePYTH);
    TH1F *heDATA; feDATA->GetObject("Vert_dR_2b_CORR",heDATA);

    TCanvas *cp = new TCanvas("cp","cp",100,200,600,400);
    hePYTH->Draw(); 
    TCanvas *cd = new TCanvas("cd","cd",200,200,600,400);
    heDATA->Draw(); 

//     hePYTH->Draw();
    double b1 = 0.0, b4 = 0.0;
    double b1e2 = 0.0, b4e2 = 0.0; //statistical errors quadratically added
    double b1DATA = 0.0, b4DATA = 0.0;
    double b1DATAe2 = 0.0, b4DATAe2 = 0.0; //statistical errors q. a. 

    if(hePYTH->GetNbinsX()!=hePYTH->GetNbinsX()) std::cout << "NOT SAME NUMBER OF BINS...\n";

    for(unsigned int i=1; i<=hePYTH->GetNbinsX(); i++){
      if(hePYTH->GetBinCenter(i)<0.8){
	b1+=hePYTH->GetBinContent(i); 
	b1e2 +=hePYTH->GetBinError(i)*hePYTH->GetBinError(i); 
	b1DATA+=heDATA->GetBinContent(i);
	b1DATAe2 +=heDATA->GetBinError(i)*heDATA->GetBinError(i); 
      }
      else if(hePYTH->GetBinCenter(i)>2.4){
	b4+=hePYTH->GetBinContent(i); 
	b4e2 +=hePYTH->GetBinError(i)*hePYTH->GetBinError(i); 
	b4DATA+=heDATA->GetBinContent(i); 
	b4DATAe2 +=heDATA->GetBinError(i)*heDATA->GetBinError(i); 
      }
    }
    std::cout << "b1, b1DATA " << b1 << " +- " << sqrt(b1e2) << " "  
	      << b1DATA << " +- " << sqrt(b1DATAe2) << std::endl;
    std::cout << "b4, b4DATA " << b4 << " +- " << sqrt(b4e2) << " " 
	      << b4DATA << " +- " << sqrt(b4DATAe2) << std::endl;

//     std::cout << strigger << " ratio b1/b4 (PYTHIA) = " << b1 << "/" << b4  << " = " << b1/b4 << "\n";
//     std::cout << strigger << " ratio b1/b4 (DATA)   = " << b1DATA << "/" << b4DATA  << " = " << b1DATA/b4DATA << "\n";
//     std::cout << strigger << " asym b1,b4 (PYTHIA) = " << b1 << "," << b4  << " = " << (b1-b4)/(b1+b4) << "\n";
//     std::cout << strigger << " asym b1,b4 (DATA)   = " << b1DATA << "," << b4DATA  << " = " << (b1DATA-b4DATA)/(b1DATA+b4DATA) << "\n";

    vratio.push_back(b1/b4);
//     hratio->Fill(t,b1/b4); 
    if(maxr<b1/b4) maxr = b1/b4; 
    if(minr>b1/b4) minr = b1/b4;

    vasym.push_back((b1-b4)/(b1+b4));
//     hasym->Fill(t,(b1-b4)/(b1+b4));  
    if(maxa<(b1-b4)/(b1+b4)) maxa = (b1-b4)/(b1+b4); 
    if(mina>(b1-b4)/(b1+b4)) mina = (b1-b4)/(b1+b4);

    vratioDATA.push_back(b1DATA/b4DATA);
//     hratioDATA->Fill(t,b1DATA/b4DATA);  
//     hratioDATA_stat->Fill(t,b1DATA/b4DATA); 
    if(maxr<b1DATA/b4DATA) maxr = b1DATA/b4DATA; 
    if(minr>b1DATA/b4DATA) minr = b1DATA/b4DATA;

    vasymDATA.push_back((b1DATA-b4DATA)/(b1DATA+b4DATA)); 
//     hasymDATA->Fill(t,(b1DATA-b4DATA)/(b1DATA+b4DATA)); 
//     hasymDATA_stat->Fill(t,(b1DATA-b4DATA)/(b1DATA+b4DATA));
    if(maxa<(b1DATA-b4DATA)/(b1DATA+b4DATA)) maxa = (b1DATA-b4DATA)/(b1DATA+b4DATA); 
    if(mina>(b1DATA-b4DATA)/(b1DATA+b4DATA)) mina = (b1DATA-b4DATA)/(b1DATA+b4DATA);

    std::cout << endl;

    //////////////////
    //systematic error: 
   //affects only data and only GSP (b1DATAe):
   double b1DATAe = sqrt(b1DATAe2); //squared statistical error 
   double b1DATAe2_stat = b1DATAe2; //statistical error for histo with stat. error only
   b1DATAe = sqrt((b1DATA*systerror)*(b1DATA*systerror) + b1DATAe2); //adding quadratically syst. error to statistical one
   b1DATAe2 = b1DATAe*b1DATAe; 
   std::cout << "Syst. error: " << systerror << " abs " << systerror*b1DATA << std::endl;
   std::cout << "b1DATAe " << sqrt(b1DATAe2) << std::endl;
    //////////////////


    //ratio: b1/b4: 
    double rate = 1.0/b4*sqrt(b1e2+b4e2/b4/b4); 
    double ratDATAe = 1.0/b4DATA*sqrt(b1DATAe2+b4DATAe2/b4DATA/b4DATA); 
    double ratDATAe_stat = 1.0/b4DATA*sqrt(b1DATAe2_stat+b4DATAe2/b4DATA/b4DATA); 
    //asym b1-b4/b1+b4
    double asyme = 2.0/(b1+b4)/(b1+b4)*sqrt(b4*b4*b1e2+b1*b1*b4e2); 
    double asymDATAe = 2.0/(b1DATA+b4DATA)/(b1DATA+b4DATA)*sqrt(b4DATA*b4DATA*b1DATAe2+b1DATA*b1DATA*b4DATAe2); 
    double asymDATAe_stat = 2.0/(b1DATA+b4DATA)/(b1DATA+b4DATA)*sqrt(b4DATA*b4DATA*b1DATAe2_stat+b1DATA*b1DATA*b4DATAe2); 
    vratioErr.push_back(rate);
    vratioDATAErr.push_back(ratDATAe);
    vratioDATAErr_stat.push_back(ratDATAe_stat);
    vasymErr.push_back(asyme);
    vasymDATAErr.push_back(asymDATAe);
    vasymDATAErr_stat.push_back(asymDATAe_stat);
//     hratio->SetBinError(t+1,rate); 
//     hratioDATA->SetBinError(t+1,ratDATAe); 
//     hratioDATA_stat->SetBinError(t+1,ratDATAe_stat); 
//     hasym->SetBinError(t+1,asyme); 
//     hasymDATA->SetBinError(t+1,asymDATAe); 
//     hasymDATA_stat->SetBinError(t+1,asymDATAe_stat); 
    std::cout << "FEHLER MC RATIO " << rate << "\n";
    std::cout << "FEHLER MC ASYM  " << asyme << "\n";
  }
  

  int n = 3; 
  double x[3] = {72.0, 106.0, 150.0}; 
  double exl[3] = {16.0,22.0,30.0}; 
  double exh[3] = {12,14,30}; 
  double yratio[3] = {vratio[0],vratio[1],vratio[2]};
  double eyratiol[3] = {vratioErr[0],vratioErr[1],vratioErr[2]}; 
  double eyratioh[3] = {vratioErr[0],vratioErr[1],vratioErr[2]}; 
  double yratioDATA[3] = {vratioDATA[0],vratioDATA[1],vratioDATA[2]};
  double eyratioDATAl[3] = {vratioDATAErr[0],vratioDATAErr[1],vratioDATAErr[2]}; 
  double eyratioDATAh[3] = {vratioDATAErr[0],vratioDATAErr[1],vratioDATAErr[2]}; 
  double yratioDATA_stat[3] = {vratioDATA[0],vratioDATA[1],vratioDATA[2]};
  double eyratioDATA_statl[3] = {vratioDATAErr_stat[0],vratioDATAErr_stat[1],vratioDATAErr_stat[2]}; 
  double eyratioDATA_stath[3] = {vratioDATAErr_stat[0],vratioDATAErr_stat[1],vratioDATAErr_stat[2]}; 
  double yasym[3] = {vasym[0],vasym[1],vasym[2]};
  double eyasyml[3] = {vasymErr[0],vasymErr[1],vasymErr[2]}; 
  double eyasymh[3] = {vasymErr[0],vasymErr[1],vasymErr[2]}; 
  double yasymDATA[3] = {vasymDATA[0],vasymDATA[1],vasymDATA[2]};
  double eyasymDATAl[3] = {vasymDATAErr[0],vasymDATAErr[1],vasymDATAErr[2]}; 
  double eyasymDATAh[3] = {vasymDATAErr[0],vasymDATAErr[1],vasymDATAErr[2]}; 
  double yasymDATA_stat[3] = {vasymDATA[0],vasymDATA[1],vasymDATA[2]};
  double eyasymDATA_statl[3] = {vasymDATAErr_stat[0],vasymDATAErr_stat[1],vasymDATAErr_stat[2]}; 
  double eyasymDATA_stath[3] = {vasymDATAErr_stat[0],vasymDATAErr_stat[1],vasymDATAErr_stat[2]}; 


  TGraphAsymmErrors *hratio = new TGraphAsymmErrors(n,x,yratio,exl,exh,eyratiol,eyratioh); 
  TGraphAsymmErrors *hratioDATA = new TGraphAsymmErrors(n,x,yratioDATA,exl,exh,eyratioDATAl,eyratioDATAh); 
  TGraphAsymmErrors *hratioDATA_stat = new TGraphAsymmErrors(n,x,yratioDATA_stat,exl,exh,eyratioDATA_statl,eyratioDATA_stath); 
  TGraphAsymmErrors *hasym = new TGraphAsymmErrors(n,x,yasym,exl,exh,eyasyml,eyasymh); 
  TGraphAsymmErrors *hasymDATA = new TGraphAsymmErrors(n,x,yasymDATA,exl,exh,eyasymDATAl,eyasymDATAh); 
  TGraphAsymmErrors *hasymDATA_stat = new TGraphAsymmErrors(n,x,yasymDATA_stat,exl,exh,eyasymDATA_statl,eyasymDATA_stath); 
  
//   hratio->GetXaxis()->SetLabelSize(0.08); 
//   hasym->GetXaxis()->SetLabelSize(0.08); 
//   hratio->GetYaxis()->SetLabelSize(0.05); 
//   hasym->GetYaxis()->SetLabelSize(0.05); 
//   hratio->GetYaxis()->SetTitleSize(0.06); 
//   hasym->GetYaxis()->SetTitleSize(0.06); 
//   hratio->GetXaxis()->SetTitleSize(0.06); 
//   hasym->GetXaxis()->SetTitleSize(0.06); 
//   hratio->GetXaxis()->SetTitleOffset(0.9); 
//   hasym->GetXaxis()->SetTitleOffset(0.9); 
//   hratio->GetXaxis()->SetBinLabel(1,"56"); 
//   hratio->GetXaxis()->SetBinLabel(2,"84"); 
//   hratio->GetXaxis()->SetBinLabel(3,"120"); 
//   hasym->GetXaxis()->SetBinLabel(1,"56"); 
//   hasym->GetXaxis()->SetBinLabel(2,"84"); 
//   hasym->GetXaxis()->SetBinLabel(3,"120"); 
  hratio->SetTitle("");
  hasym->SetTitle("");
  hratio->GetYaxis()->SetTitle("GSP/FCR"); 
  hasym->GetYaxis()->SetTitle("GSP-FCR/GSP+FCR"); 
  hratio->GetXaxis()->SetTitle("Leading Jet p_{T} (GeV)"); 
  hasym->GetXaxis()->SetTitle("Leading Jet p_{T} (GeV)"); 
  hratio->GetYaxis()->SetTitleOffset(1.15); 
  hasym->GetYaxis()->SetTitleOffset(1.15); 
  hratio->GetXaxis()->SetTitleOffset(1.15); 
  hasym->GetXaxis()->SetTitleOffset(1.15); 


//   hratio->SetMarkerStyle(20); 
//   hasym->SetMarkerStyle(20); 
  hratio->SetMarkerColor(8); 
  hasym->SetMarkerColor(8); 

  hratioDATA->SetMarkerStyle(20); 
  hasymDATA->SetMarkerStyle(20); 
  hratioDATA->SetLineWidth(2);
  hasymDATA->SetLineWidth(2);
  hratioDATA_stat->SetLineWidth(2);
  hasymDATA_stat->SetLineWidth(2);

  hratio->SetFillColor(8); 
  hasym->SetFillColor(8);
  hratio->SetLineColor(8); 
  hasym->SetLineColor(8);

//   TH1F *hratioDATA_stat = hratioDATA->Clone(); 
//   TH1F *hasymDATA_stat = hasymDATA->Clone(); 

  ///here was for


  hratio->GetYaxis()->SetRangeUser(minr-fabs(minr/4.0),maxr+maxr/5.0); 
  hasym->GetYaxis()->SetRangeUser(mina-fabs(mina/2.0),maxa+maxa/5.0); 

  
  TCanvas *c1 = new TCanvas("c1","c1",300,200,600,600); 
//   c1->SetGridy(); 
  c1->SetBottomMargin(0.14);
  c1->SetLeftMargin(0.13); 

  hratio->Draw("AE2"); 
  hratioDATA->Draw("PE1same"); 
  hratioDATA_stat->Draw("PE1same"); 

  //leg
  TLegend *leg = new TLegend(0.6,0.25,0.9,0.4);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hratio,"PYTHIA","F");
  leg->AddEntry(hratioDATA,"corrected DATA","LP");
  leg->Draw(); 


  TLatex * tex = new TLatex(45,3.4,"CMS #sqrt{s}=7 TeV, L = 3 pb^{-1}");
  tex->SetTextSize(0.040);
  tex->SetLineWidth(2);
  tex->Draw();



  TCanvas *c2 = new TCanvas("c2","c2",500,200,600,600); 
//   c2->SetGridy();
  c2->SetBottomMargin(0.14);
  c2->SetLeftMargin(0.13); 
  hasym->Draw("AE2"); 
  hasymDATA->Draw("PE1same"); 
  hasymDATA_stat->Draw("PE1same"); 

  leg->Draw(); 
  TLatex * tex2 = new TLatex(45,0.57,"CMS #sqrt{s}=7 TeV, L = 3 pb^{-1}");
  tex2->SetTextSize(0.040);
  tex2->SetLineWidth(2);
  tex2->Draw();

  
//   //MIKKO
//   TLatex *tt = new TLatex();
//   tt->SetTextSize(0.05);
//   tt->SetTextFont(42);
//   tt->SetTextAlign(31); // align right

//   const int nylab = 12;
//   const double ylabels[nylab] = {0, 1, 2, 1, 2, 1, 0, 1, 2, 3, 4, 5};
//   for (int ii = 0; ii != nylab; ++ii)
//     tt->DrawTextNDC(0.14, 0.135+0.0665*ii, Form("%1.0f",ylabels[ii]));


}


void getXcoordinates(string strigger="jet30"){


  TFile *fDATA = TFile::Open("/shome/wehrlilu/SecVtx/SVProd/SVVALIDATION/forbatchjobs/temp/BBC/CMSSW_3_7_0_patch2/crab/analyze/rootfiles/dataFiles/withIP/allDATAupTo144114.root","READ");

  TTree *tDATA = (TTree*)fDATA->Get("bcanalyzer/tEvents");
  TH1F *lpjpt = new TH1F("lpjpt","leading jet pt",1000,0,1000); 

  string cut; 
  if(strigger=="jet15") cut="BEvents.jet15==1 && BEvents.ptHardestPJ>56  && BEvents.ptHardestPJ<84";
  else if(strigger=="jet30") cut="BEvents.jet30==1 && BEvents.ptHardestPJ>84 && BEvents.ptHardestPJ<120";
  else if(strigger=="jet50") cut="BEvents.jet50==1 && BEvents.ptHardestPJ>120";
  else {std::cout << "no correct trigger path\n"; break;}

  tDATA->Draw("BEvents.ptHardestPJ>>lpjpt",(cut).c_str());

  lpjpt->Draw(); 

  double mean = 0.0, ntot=0.0;
  for(unsigned int i=1; i<1001; i++){
    mean += lpjpt->GetBinContent(i)*lpjpt->GetBinCenter(i);
    ntot += lpjpt->GetBinContent(i); 
    
  }

  std::cout << strigger << " mean " << lpjpt->GetMean() << " +- " << lpjpt->GetMeanError() << "\n";
}
