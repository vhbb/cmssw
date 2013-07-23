{
  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();

  string jetType = "PF";
  TFile *f0 = TFile::Open("/scratch/leo/anV3-MinimumBias-Commissioning10-GOODCOLL-Jun14thSkim_v1-V9.root");
  TFile *f1 = TFile::Open("/scratch/leo/anV3-MinimumBias-Commissioning10-SD_JetMETTau-Jun14thSkim_v1-V9.root");
  TFile *f2 = TFile::Open("/scratch/leo/anV3-JetMET-Run2010A-PromptReco-v4_Runs141950-143731-V9B.root");

  TTree *t0 = (TTree*) f0->Get("bcanalyzer/tEvents");
  TTree *t1 = (TTree*) f1->Get("bcanalyzer/tEvents");
  TTree *t2 = (TTree*) f2->Get("bcanalyzer/tEvents");

  t0->Draw( "PFJet_pt[0] >>h_HLT_L1Jet6U_den(80,0,200)", "eventRunNumber<135518" ) ;
  t0->Draw( "PFJet_pt[0] >>h_HLT_L1Jet6U_num(80,0,200)", "HLT_L1Jet6U==1 && eventRunNumber<135518") ;
  t0->Draw( "PFJet_pt[0] >>h_HLT_Jet15U_den(80,0,200)", "abs(HLT_L1Jet6U_PFJet_eta)<3 && HLT_L1Jet6U==1" ) ;
  t0->Draw( "PFJet_pt[0] >>h_HLT_Jet15U_num(80,0,200)", "abs(HLT_L1Jet6U_PFJet_eta)<3  && HLT_L1Jet6U==1 && HLT_Jet15U==1") ;
  t1->Draw( "PFJet_pt[0] >>h_HLT_Jet30U_den(80,0,200)", "abs(HLT_Jet15U_PFJet_eta)<3 && HLT_Jet15U==1 " ) ;
  t1->Draw( "PFJet_pt[0] >>h_HLT_Jet30U_num(80,0,200)", "abs(HLT_Jet15U_PFJet_eta)<3  && HLT_Jet15U==1 && HLT_Jet30U==1 ") ;
  t2->Draw( "PFJet_pt[0] >>h_HLT_Jet50U_den(80,0,200)", "abs(HLT_Jet30U_PFJet_eta)<3 && HLT_Jet30U==1") ;
  t2->Draw( "PFJet_pt[0] >>h_HLT_Jet50U_num(80,0,200)", "abs(HLT_Jet30U_PFJet_eta)<3  &&abs(HLT_Jet50U_PFJet_eta)<3 && ( HLT_Jet30U==1 && HLT_Jet50U==1) " ) ;
  //t2->Draw( "PFJet_pt[0] >>h_HLT_Jet70U_den(80,0,200)", "abs(HLT_Jet50U_PFJet_eta)<3 && HLT_Jet50U==1 && HLT_Jet50U_prescale==1") ;
  //t2->Draw( "PFJet_pt[0] >>h_HLT_Jet70U_num(80,0,200)", "abs(HLT_Jet50U_PFJet_eta)<3  &&( HLT_Jet50U==1 && HLT_Jet70U==1) && HLT_Jet50U_prescale==1" ) ;

  //h_HLT_Jet30U_den->Draw();

  TCanvas *myc = new TCanvas( ("TurnOn_"+jetType).c_str(), ("TurnOn_"+jetType).c_str() );
  TGraphAsymmErrors *g[10];
  TMultiGraph *mg = new TMultiGraph();
  TLegend *l = new TLegend(0.1,0.2,0.3,0.4);
  

  //g[0] = new TGraphAsymmErrors(h_HLT_Jet15U_num, h_HLT_Jet15U_den);
  //g[0]->SetLineColor(kBlack); l->AddEntry(g[0],"HLT_Jet15U","l"); mg->Add(g[0]);

  int g_i=0;
  g[g_i] = new TGraphAsymmErrors(h_HLT_L1Jet6U_num, h_HLT_L1Jet6U_den);
  g[g_i]->SetLineColor(kBlack); l->AddEntry(g[g_i],"HLT_L1Jet6U","l"); mg->Add(g[g_i++]);
  g[g_i] = new TGraphAsymmErrors(h_HLT_Jet15U_num, h_HLT_Jet15U_den);
  g[g_i]->SetLineColor(kGreen); l->AddEntry(g[g_i],"HLT_Jet15U","l"); mg->Add(g[g_i++]);
  g[g_i] = new TGraphAsymmErrors(h_HLT_Jet30U_num, h_HLT_Jet30U_den);
  g[g_i]->SetLineColor(kBlue); l->AddEntry(g[g_i],"HLT_Jet30U","l"); mg->Add(g[g_i++]);
  g[g_i] = new TGraphAsymmErrors(h_HLT_Jet50U_num, h_HLT_Jet50U_den);
  g[g_i]->SetLineColor(kRed); l->AddEntry(g[g_i],"HLT_Jet50U","l"); mg->Add(g[g_i++]);
  

/*
  g[1] = new TGraphAsymmErrors(h_pt_HLT_Jet15U, h_pt1);
  g[1]->SetLineColor(kBlue); l->AddEntry(g[1],"HLT_Jet15U","l");mg->Add(g[1]);
  g[2] = new TGraphAsymmErrors(h_pt_HLT_Jet30U, h_pt2);
  g[2]->SetLineColor(kRed); l->AddEntry(g[2],"HLT_Jet30U","l");mg->Add(g[2]);
  g[3] = new TGraphAsymmErrors(h_pt_HLT_Jet50U, h_pt3);
  g[3]->SetLineColor(kGreen); l->AddEntry(g[3],"HLT_Jet50U","l");mg->Add(g[3]);
  */

  mg->SetTitle(jetType.c_str());
  
  mg->Draw("apl");
  mg->GetXaxis()->SetTitle("jet p_{T}");
  mg->GetYaxis()->SetTitle("turn on");

  mg->Draw("apl");
  l->Draw();

}
