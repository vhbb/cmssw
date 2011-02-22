int script(void)
{


  gROOT->ProcessLine(".L ~/tdrstyle.C");
  setTDRStyle();

  //style modifications:
  gStyle->SetErrorX(0.5);
  gStyle->SetPadTopMargin(0.08);//.005
  gStyle->SetPadBottomMargin(0.10);//0.13
  gStyle->SetMarkerStyle(1);
//   gStyle->SetStatFont(62); 
//   gStyle->SetTitleFont(42);
//   gStyle->SetTitleFont(62);
//   gStyle->SetTitleFont(62, "XYZ");
//   gStyle->SetLabelFont(62, "XYZ");
//   gStyle->SetTitleFontSize(0.02);



  //TFile *file = new TFile ("sv_out_30um_10000_corrected1_no_cuts.root");
	//TFile *file = new TFile ("sv_out_30um_MC_full_sample.root");
  
	
	//TFile *file = new TFile ("sv_out_30um_MC_new_code.root");
  // TFile *file = new TFile ("sv_out_30um_Z_10000_combinations.root");
	//TFile *file = new TFile ("sv_out_30um_Z_MC_PVcheck.root");
	//TFile *file = new TFile ("sv_out_30um_Z_Data_PVcheck.root");
//	TFile *file = new TFile ("mcall2.root");
//	TFile *file1 = new TFile ("mcall2.root");
	TFile *file = new TFile ("evMIX_sv_out_20um_3D_MC_pt_50-80.root");
	TFile *file1 = new TFile ("evMIX_sv_out_20um_3D_MC_pt_30-inf-all.root");
//sv_out_20um_3D_MC_pt_30-inf.root");
	//TFile *file2 = new TFile ("sv_out_30um_really3D_data.root");
	TFile *file2 = new TFile ("evMIX_all-pv.root");
	//TFile *file = new TFile ("sv_out_30um_Z_MC_pt_50-80.root");
	//TFile *file = new TFile ("sv_out_20um_3D_MC_pt_50-80.root");
  
	TH1F::SetDefaultSumw2(kTRUE);
double a[15] = {0.0,0.2,0.4,0.8,2.0,2.8,3.2,3.6,4.4,6.0};
int nbins = 9;
//double a[15] = {0.0,0.8,2.4,4.0,6.0};
//int nbins = 4;
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
/*  TH1F *h1 = new TH1F("h1","h1",nbins,a);
  TH1F *h2 = new TH1F("h2","h2",nbins,a);
  TH1F *h11 = new TH1F("h11","h11",nbins,a);
  TH1F *h12 = new TH1F("h12","h12",nbins,a);
  TH1F *h21 = new TH1F("h21","h21",nbins,a);
  TH1F *h22 = new TH1F("h22","h22",nbins,a);
*/

  TH1F *h1 = new TH1F("h1","h1",10,0.,4.);
  TH1F *h2 = new TH1F("h2","h2",10,0.,4.);
  TH1F *h11 = new TH1F("h11","h11",10,0.,4.);
  TH1F *h12 = new TH1F("h12","h12",10,0.,4.);
  TH1F *h21 = new TH1F("h21","h22",10,0.,4.);
  TH1F *h22 = new TH1F("h22","h22",10,0.,4.);

/*  TH1F *h1 = new TH1F("h1","h1",60,0.,6.);
  TH1F *h2 = new TH1F("h2","h2",60,0.,6.);
  TH1F *h11 = new TH1F("h11","h11",60,0.,6.);
  TH1F *h12 = new TH1F("h12","h12",60,0.,6.);
  TH1F *h21 = new TH1F("h21","h22",60,0.,6.);
  TH1F *h22 = new TH1F("h22","h22",60,0.,6.);
*/
TCanvas *c1=new TCanvas("c1","c1",800,600);
  c1->Divide(3,3);
  c1->cd(1);
c1->SetObjectStat(false);

  TTree * Delta_r =  (TTree *)file->Get("svAnalysis/events/Delta_r");
  TTree * Delta_R =  (TTree *)file->Get("svAnalysis/events/Delta_R");

  Delta_r->Draw("delta_r>>h1","delta_r!=0 && jetpt > 50","");
  c1->cd(2);
  Delta_R->Draw("Delta_R>>h2","Delta_R!=0 && jetpt > 50","");

  Delta_r =  (TTree *)file1->Get("svAnalysis/events/Delta_r");
  Delta_R =  (TTree *)file1->Get("svAnalysis/events/Delta_R");

  c1->cd(4);
  Delta_r->Draw("delta_r>>h11","delta_r!=0 && jetpt > 50","");
  c1->cd(5);
  Delta_R->Draw("Delta_R>>h12","Delta_R!=0 && jetpt > 50","");

  Delta_r =  (TTree *)file2->Get("svAnalysis/events/Delta_r");
  Delta_R =  (TTree *)file2->Get("svAnalysis/events/Delta_R");

  c1->cd(7);
  Delta_r->Draw("delta_r>>h21","delta_r!=0 && maxdist < 0.002 && jetpt > 50","");
  c1->cd(8);
  Delta_R->Draw("Delta_R>>h22","Delta_R!=0 && maxdist < 0.002 && jetpt > 50","");

  TH1F *rat1 =(TH1F *) h2->Clone("Ratio_5080");
  rat1->Divide(h2,h1,1.,1.,"B"); 
  TH1F *rat2 =(TH1F *) h21->Clone("Ratio_30inf");
  rat2->Divide(h12,h11,1.,1.,"B"); 
  TH1F *rat3 =(TH1F *) h22->Clone("Data");
  rat3->Divide(h22,h21,1.,1.,"B"); 
  rat3->SetMarkerStyle(20); 

  rat2->SetFillColor(4);
  rat2->SetLineColor(0); 
  rat2->SetFillStyle(3001);
  rat1->SetFillStyle(3005);
  rat1->SetFillColor(2);
  rat1->SetLineColor(0); 


  c1->cd(3);
  rat1->Draw("E3"); 
  c1->cd(6);
  rat2->Draw("E3"); 
  c1->cd(9);
  rat3->Draw("E1"); 
//  c1->SaveAs("c1.root");
//  c1->SaveAs("c1.png");

  TCanvas *c2=new TCanvas("c2","c2",600,600);
  c2->SetObjectStat(false);
  TLegend *le =  new TLegend(0.2,0.15,0.8,0.40);
  le->AddEntry(rat2,"Simulation (#hat{p}_{T}>30 GeV)","F");
  le->AddEntry(rat1,"Simulation (50 GeV < #hat{p}_{T} < 80 GeV)","F");
  le->AddEntry(rat3,"Data","P");
  le->SetFillColor(0);
  le->SetBorderSize(0); 
//  rat2->Draw("E3same"); 
  rat2->Draw("E2"); 
  rat2->GetXaxis()->SetTitle("\\Delta R");
  rat2->GetYaxis()->SetTitle("relative efficiency");
  rat2->GetYaxis()->SetTitleOffset(1.1);
  rat2->GetYaxis()->SetRangeUser(0.0,1.0); 
  rat2->GetYaxis()->SetTitleSize(0.045); //0.035
  rat2->GetYaxis()->SetLabelSize(0.038); //0.03
  rat2->GetXaxis()->SetTitleSize(0.045); //0.035
  rat2->GetXaxis()->SetLabelSize(0.038); //0.03

  TLatex *   tex1 = new TLatex(0.0,1.05,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
  tex1->SetTextSize(0.040); //0.044
  tex1->SetLineWidth(2);
  tex1->Draw();

  rat1->Draw("E2same"); 
  rat3->Draw("E1same"); 
  rat2->SetFillColor(2);
  rat2->SetMarkerColor(2); 
  rat2->SetMarkerSize(0); 
  rat2->SetMarkerStyle(7); 	
  rat1->SetFillColor(4);
  rat1->SetMarkerColor(4); 
  rat1->SetMarkerSize(0); 	
  rat1->SetMarkerStyle(7); 	

  le->Draw();
  c2->SetTickx(); 
  c2->SetTicky(); 
  c2->RedrawAxis(); 
//  c2->SaveAs("c2.root");
//  c2->SaveAs("c2.png");
  TCanvas *c3=new TCanvas("c3","c3",600,600);
  c3->SetObjectStat(false);
  c3->cd();
  TH1F *rat11 = (TH1F *) rat3->Clone("dataoverMC5080");
  TH1F *rat12 = (TH1F *) rat3->Clone("dataoverMC30");
  rat11->Divide(rat1);
  rat12->Divide(rat2);
  rat11->SetFillColor(4);
  rat12->SetFillColor(2);
  rat11->SetFillStyle(3001);
  rat12->SetFillStyle(3005);
  rat11->SetLineColor(4);
  rat12->SetLineColor(2);

  rat12->Draw("E2");
  rat12->GetXaxis()->SetTitle("\\Delta R");
  rat11->Draw("E2same");

  TH1F *rat21 = (TH1F *) rat11->Clone("dataoverMC5080_scaled");
  TH1F *rat22 = (TH1F *) rat12->Clone("dataoverMC30_scaled");

  rat11->Fit("pol0","","",0.8,4);
  rat11->GetFunction("pol0")->SetLineColor(4);
  double scale1 =  rat11->GetFunction("pol0")->GetParameter(0);
  rat12->Fit("pol0","","",0.8,4);
  rat12->GetFunction("pol0")->SetLineColor(2);
  double scale2 =  rat12->GetFunction("pol0")->GetParameter(0);

  rat12->Draw("E2");
  rat12->GetXaxis()->SetTitle("\\Delta R");
  rat11->Draw("E2same");
  //c3->SaveAs("c3.root");
  //c3->SaveAs("c3.png");

  TCanvas *c4=new TCanvas("c4","c4",600,600);
  c4->SetObjectStat(false);
  rat21->Scale(1./scale1);
  rat22->Scale(1./scale2);
  rat21->SetFillStyle(3001);
  rat22->SetFillStyle(3005);
  rat22->SetMarkerStyle(20);
  rat22->SetMarkerColor(1);
  rat22->SetLineColor(1);
  rat22->Draw("E1");
  rat22->GetXaxis()->SetTitle("\\Delta R");
  rat22->GetYaxis()->SetTitle("data/simulation scaled efficiencies");
  rat22->GetYaxis()->SetTitleSize(0.045); //0.035
  rat22->GetYaxis()->SetLabelSize(0.035); //0.03
  rat22->GetXaxis()->SetTitleSize(0.045); //0.035
  rat22->GetXaxis()->SetLabelSize(0.035); //0.03

  rat22->GetYaxis()->SetRangeUser(0.9,1.1);
// //   rat22->GetYaxis()->SetTitleOffset(1.5);
  gStyle->SetGridColor(kGray+2); 
  c4->SetGridy(1); 
  TLatex *   tex = new TLatex(0.0,1.11,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
  tex->SetTextSize(0.040); //0.044
  tex->SetLineWidth(2);
  tex->Draw();
  c4->SetTickx(); 
  c4->SetTicky(); 
  c4->RedrawAxis(); 

//  c4->SaveAs("c4.root");
//  c4->SaveAs("c4.png");

//rat21->Draw("E2same");


  //c2->Divide(2,1);
  //c2->cd(1);

  //  h1->SetLineColor(1);
  // h2->SetLineColor(2);

 // h1->Sumw2(); 
  //h2->Sumw2();
  
//TH1F *rat = h1->Clone(); rat->SetName("Ratio");
     // Clone one of the histograms
  //h1->Draw();
  //   //   hSelected->SetLineColor(2);
  //   //   hSelected->Draw("same");
  //c1->cd(2);
  //h2->Draw();
    // rat->Divide(h2,h1,1.,1.,"B");
    //rat->Draw();
//   setTDRStyle();

}
