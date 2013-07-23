#include <string>
#include <vector>
#include <map>
#include <TH1F.h>
#include <TProfile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <math.h>
#include <TVector.h>
#include <iostream>
double a[15] = {0.0,0.8,2.4,4.0};
int nbins = 3;
//double a[15] = {0.0,0.2,0.4,0.8,2.4,3.2,4.0};
//int nbins = 6;
//double a[15] = {0.0,0.2,0.4,0.8,2.0,2.8,3.2,3.6,4.4,6.0};
//int nbins = 9;

std::map<std::string,std::string> translate;

struct DrawStyle
{
 std::string def;
 std::string legend;
 int marker; 
 int fill; 
};

struct Sample
{
 Sample(){}
 Sample(const char *n,const char *l,int c, DrawStyle d,std::string hn="",bool drs=false): drawSystematics(drs), histoName(hn),name(n),label(l),prefix("/scratch/leo/finalOpenBins36X_No6U_NewDPhiBins/histo"),postfix("BBCorr_HadrJetMatchPt_OpenHLT-total.root"),color(c),ds(d) {}
 bool drawSystematics;
 std::string histoName;
 std::string name;
 std::string label;
 std::string prefix;
 std::string postfix;
 int color;
 std::map<std::string,TFile *> tfile;
 std::string fileName(std::string trigger)
   {
     return prefix+"-"+name+"-"+trigger+"-"+postfix;
   }
 TFile * file(std::string trigger)
   {
     if(tfile.find(trigger) == tfile.end()) tfile[trigger] = TFile::Open(fileName(trigger).c_str());
     return tfile[trigger];
   } 
 TObject * get(std::string trigger, std::string name)
   {
    if(file(trigger) == 0) return 0;
    return  file(trigger)->Get(name.c_str());
   }

 DrawStyle ds;
};


struct PtBin
{
 PtBin() {}
 PtBin(std::string tr,double x,double xl,double xh, double systematics1, double systematics2):
  trigger(tr),centralValue(x),exl(xl),exh(xh),sys1(systematics1),sys2(systematics2) {}

 std::string trigger;
 double centralValue;
 double exl,exh;
 double sys1; // linearly added
 double sys2; // quadratically added
};

struct GraphPoint
{
 float y;
 float eyl;
 float eyh;
 float esysyl;
 float esysyh;
 bool invalid;
};

struct PlotterVsPtBin
{
 std::vector<Sample> samples;
 std::vector<PtBin>  bins;
 virtual std::string titleX() {return  "leading jet p_{T} (GeV)";}
 virtual std::string titleY() {return  "";}
 virtual TLegend * legend() {
//   TLegend * l = new TLegend(0.6,0.1,0.9,0.28);
   TLegend * l = new TLegend(0.45,0.15,0.8,0.30);
   l->SetFillColor(0);
   l->SetBorderSize(0);
   return l;
 }

 virtual GraphPoint computePoint(Sample s,PtBin b) =0 ; 
 virtual TGraphAsymmErrors * draw(std::string options,Sample s) 
 {
    std::cerr << "DRAW" << std::endl;
    //workaround for madgraph having only 3 points
    int n=0; //bins.size();
    for(int i=0;i<bins.size();i++) {
     if(!computePoint(s,bins[i]).invalid) n++;
    }
    TVectorD x(n),y(n),exl(n),exh(n),eyl(n),eyh(n),esysyl(n),esysyh(n);
    int i=0;
    for(int j=0;j<bins.size();j++) {
     GraphPoint p = computePoint(s,bins[j]);
     if(p.invalid) continue;
     x[i]=bins[j].centralValue;
     exl[i]=bins[j].exl;
     exh[i]=bins[j].exh;
     y[i]=p.y;
     eyl[i]=p.eyl;
     eyh[i]=p.eyh;
     esysyl[i]=p.esysyl;
     esysyh[i]=p.esysyh;
     std::cout << "TEXT " << s.name << " "  << x[i] << " " << y[i] << " " << eyl[i] << " " << esysyl[i] << std::endl; 
     i++;
    }
    TGraphAsymmErrors * gr = new TGraphAsymmErrors(x,y,exl,exh,eyl,eyh);
    TGraphAsymmErrors * grSys = new TGraphAsymmErrors(x,y,exl,exh,esysyl,esysyh);
    gr->SetLineColor(s.color);
    gr->SetFillColor(s.color);
    gr->SetMarkerColor(s.color);
    gr->SetMarkerStyle(s.ds.marker);
    std::cout<< options.c_str() << std::endl;
    gr->Draw(options.c_str());
    if(s.drawSystematics)   grSys->Draw((options+"same").c_str());

    return gr;
 }


 void doAll(std::string name)
 {
   TCanvas * c1 = new TCanvas(name.c_str(),name.c_str(),100,300,500,500);
   if(name=="testasymmetry"){
//      c1->Range(44.8,-0.2,161.8,0.8);
     c1->SetLeftMargin(0.1653226);
     c1->SetRightMargin(0.03427419);
   }
   c1->SetTicky();
   c1->SetTickx();
   gStyle->SetOptStat(0);
   c1->SetObjectStat(false);
   TLegend * l = legend();
   std::string options="A";
   float mx=-1e99,mi=1e99;
   TAxis * y_axis=0;
   TAxis * x_axis=0;
   std::cout << "Loop on samples" << samples.size() << std::endl;
   for(unsigned int i=0;i< samples.size();i++)
   {
    std::cerr << "sample" << i  << std::endl;
     TGraphAsymmErrors & d = * draw(options+samples[i].ds.def,samples[i]);
     std::cout << &d << std::endl;
     l->AddEntry(&d,samples[i].label.c_str(),samples[i].ds.legend.c_str());
     if(d.GetMaximum() > mx) mx=d.GetMaximum();
     if(d.GetMaximum() < mi) mi=d.GetMinimum();

     if(options == "A")
      {
        y_axis=d.GetYaxis();
        x_axis=d.GetXaxis();
        options="same";
      }
     if(name=="trendplot"){
       d.SetFillStyle(samples[i].ds.fill); 
     }
   }
   mx*=1.2;
   if(mi>0) mi*=0.8; 
   if(mi<0) mi*=1.2; 
   std::cout << mi << " " <<mx <<std::endl; 
   mx = 4; mi = 0; 
   if(name=="testratio"){
     mx=4;
     mi=0;
   }
   else if(name=="testasymmetry"){
     mx=0.7;
     mi=-0.1;
   }
   else if(name=="trendplot"){
     mx=6.0; 
     mi=0.0; 
   }
   x_axis->SetTitle(titleX().c_str());
   y_axis->SetTitle(titleY().c_str());
   y_axis->SetRangeUser(mi,mx);
   y_axis->SetTitleOffset(1.19);
   if(name=="testasymmetry") y_axis->SetTitleOffset(1.36); 

//   x_axis->SetRangeUser(0,3.9);
   if(name=="trendplot"){
     l->SetX1NDC(0.5947581); 
     l->SetY1NDC(0.1737288); 
     l->SetX2NDC(0.9435484); 
     l->SetY2NDC(0.3241525); 

// 0.5947581,0.1737288,0.9435484,0.3241525
   }
   l->Draw("same");
   TLatex * tex, * tex2; 
   tex = new TLatex(bins[0].centralValue,(mx-mi)*0.03+mx,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
   if(name=="testratio" || name=="testasymmetry") tex = new TLatex(64,(mx-mi)*0.03+mx,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
   if(name=="trendplot")  tex = new TLatex(43.1,(mx-mi)*0.03+mx,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
   tex->SetTextSize(0.040); //0.044
   tex->SetLineWidth(2);
   tex->Draw();

   tex2 = new TLatex(bins[0].centralValue,3.5,"#splitline{p_{T}^{B} > 15 GeV, |\\eta^{B}| < 2.0}{|\\eta^{Jet}| < 3.0} ");
   if(name=="testasymmetry") tex2 = new TLatex(bins[0].centralValue,0.6,"#splitline{p_{T}^{B} > 15 GeV, |\\eta^{B}| < 2.0}{|\\eta^{Jet}| < 3.0} ");
   if(name=="trendplot") tex2 = new TLatex(120.2976,5.128152,"#splitline{p_{T}^{B} > 15 GeV, |\\eta^{B}| < 2.0}{|\\eta^{Jet}| < 3.0} "); 
   tex2->SetTextSize(0.04);
   tex2->SetLineWidth(2);
   tex2->Draw();

   c1->SaveAs((name+".root").c_str());
   c1->SaveAs((name+".png").c_str());
   c1->SaveAs((name+".pdf").c_str());
 }

};

struct PlotterRatio : public PlotterVsPtBin
{
 virtual std::string titleY() {return  "\\sigma_{\\DeltaR < 0.8} / \\sigma_{\\DeltaR > 2.4}";}
 virtual GraphPoint computePoint(Sample s,PtBin b) {
    GraphPoint r;
    TH1F * h =(TH1F*) s.get(b.trigger,s.histoName.c_str());
    if(h==0) {
       std::cout<<" ERROR " <<  b.trigger << " " << s.histoName << std::endl;   
       r.invalid=true;
       return r;
     }
    r.invalid=false;
    double b1=0,b4=0,b1e2=0,b4e2=0; // values and stat error
    double b1e2_sy2=0,b4e2_sy2=0; // quad added syst error
    for(unsigned int i=1; i<=h->GetNbinsX(); i++){
      if(h->GetBinCenter(i)<0.8){
        b1+=h->GetBinContent(i);
        b1e2 +=h->GetBinError(i)*h->GetBinError(i);
        b1e2_sy2+=b.sys2*b.sys2*h->GetBinContent(i)*h->GetBinContent(i);
      }
      if(h->GetBinCenter(i)>2.4){
        b4+=h->GetBinContent(i);
        b4e2 +=h->GetBinError(i)*h->GetBinError(i);
        b4e2_sy2+=b.sys2*b.sys2*h->GetBinContent(i)*h->GetBinContent(i);
      }
    }
   std::cout << "TEXT2 " << s.name << " " << b.centralValue <<  " "  << b1 << " " << b4 << " " << sqrt(b1e2+b1e2_sy2) << " " << sqrt(b4e2+b4e2_sy2)  << std::endl;

    cout << b1e2 << " " <<  b1e2_sy2 << endl;
    r.y = b1/b4;
    if(b4==0) r.y=0;
//  r.eyl = 1.0/b4*sqrt(b1e2+b4e2/b4/b4*b1*b1);
    r.eyl = r.y*sqrt(b1e2/b1/b1+b4e2/b4/b4);
    r.eyh = r.eyl;
    double sys2 = r.y*r.y*(b4e2_sy2/b4/b4+b1e2_sy2/b1/b1);
    r.esysyh = sqrt(r.eyl*r.eyl+b.sys1*b.sys1*r.y*r.y+sys2);
    r.esysyl = r.esysyh;
    return r;
 }

};

struct PlotterAsym : public PlotterVsPtBin
{
 virtual std::string titleY() {return  "#frac{\\sigma_{\\DeltaR < 0.8} - \\sigma_{\\DeltaR > 2.4}}{\\sigma_{\\DeltaR < 0.8} + \\sigma_{\\DeltaR > 2.4}}";}
 virtual GraphPoint computePoint(Sample s,PtBin b) {
    GraphPoint r;
    TH1F * h =(TH1F*) s.get(b.trigger,s.histoName.c_str());
    if(h==0) {
       std::cout<<" ERROR " <<  b.trigger << " " << s.histoName << std::endl;
       r.invalid=true;
       return r;
     }
    r.invalid=false;
    double b1=0,b4=0,b1e2=0,b4e2=0; // values and stat error
    double b1e2_sy2=0,b4e2_sy2=0; // quad added syst error
    for(unsigned int i=1; i<=h->GetNbinsX(); i++){
      if(h->GetBinCenter(i)<0.8){
        b1+=h->GetBinContent(i);
        b1e2 +=h->GetBinError(i)*h->GetBinError(i);
        b1e2_sy2+=b.sys2*b.sys2*h->GetBinContent(i)*h->GetBinContent(i);
      }
      if(h->GetBinCenter(i)>2.4){
        b4+=h->GetBinContent(i);
        b4e2 +=h->GetBinError(i)*h->GetBinError(i);
        b4e2_sy2+=b.sys2*b.sys2*h->GetBinContent(i)*h->GetBinContent(i);
      }
    }
    r.y = (b1-b4)/(b1+b4);
//    if(b4==0) r.y=0;
    cout << b1e2 << " " <<  b1e2_sy2 << endl;
    r.eyl = 2.0/(b1+b4)/(b1+b4)*sqrt(b4*b4*b1e2+b1*b1*b4e2);
    r.eyh = r.eyl;
    double R=b1/b4;
    double coef=2.*R/(R+1.)/(R+1.);
    double sys = 2.0/(b1+b4)/(b1+b4)*sqrt(b4*b4*b1e2_sy2+b1*b1*b4e2_sy2);
    r.esysyh = sqrt(r.eyl*r.eyl+b.sys1*b.sys1*coef*coef+sys*sys);
    r.esysyl = r.esysyh;
    return r;
 }

};


struct PlotterByTriggers
{
 std::vector<std::string> triggers;
 std::vector<Sample> samples;
 virtual TH1 * draw(std::string options,Sample s,std::string trigger) {return 0; }
 virtual std::string titleX() { return "";} 
 virtual std::string titleY() { return "";} 
 virtual std::string filename() {return "tmp";}
 virtual TLegend * legend() { 
//   TLegend * l = new TLegend(0.6,0.1,0.9,0.28);
//   TLegend * l = new TLegend(0.45,0.25,0.8,0.40); // no same
//    TLegend * l = new TLegend(0.25,0.43,0.5,0.68); // same
//above for others than pt?
   TLegend * l = new TLegend(0.32,0.50,0.5,0.75); // same


   l->SetBorderSize(0);
   l->SetFillColor(0);
   return l;
 }
 void doAll(std::string name, bool same=false)
 {
   TPad * c1;
   TPad * c2;
  if(same) 
   {
      c2  = new TCanvas("all","all",100,300,500,500);
      TLatex *   t1 = new TLatex(0.88,0.04,"#Delta R");
      t1->SetTextSize(0.038); //0.044
      t1->SetLineWidth(2);
      t1->Draw();

      TPad *pad = new TPad("pad1", "The pad 20% of the height",0.0,0.90,1.0,1.0);
      TPad *pad2 = new TPad("pad2", "side",0.0,0.0,0.1,1.0);
      pad2->SetFillStyle(0);

      c1 = new TPad("pad2", "The pad 80% of the height",0.0,0.07,1.0,0.9303);
      c1->SetFillStyle(0);
      pad->Draw();
      pad2->Draw();

      c1->Draw();
      c1->Divide(1,triggers.size());
      pad->cd();
//       TLatex *   tex = new TLatex(0.11,0.15,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
      TLatex *   tex = new TLatex(0.15,0.42,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");

      tex->SetTextSize(0.4); //0.044
      tex->SetLineWidth(2);
      tex->Draw();
      pad2->cd();
//       TLatex *   tex2 = new TLatex(0.38,0.265,"average p_{T}^{rec} softer B candidate (GeV)");
      TLatex *   tex2 = new TLatex(0.8,0.265,"average p_{T}^{rec} softer B candidate (GeV)");
      tex2->SetTextAngle(90);
      tex2->SetTextSize(0.38); //0.044
      tex2->SetLineWidth(2);
      tex2->Draw();

   } 
   for(int i=0;i< triggers.size();i++) {
   std::string tr = triggers[i];
   if(same) 
   {
     cout << "cd " << i << endl; 
    TVirtualPad * p = c1->cd(i+1);
     p->SetTicky();
     p->SetTickx();
     p->SetFillStyle(0);

   }
     else
   {  c1  = new TCanvas(tr.c_str(),tr.c_str(),100,300,500,500);}
    
   c1->SetTicky();
   c1->SetTickx();
   gStyle->SetOptStat(0);
   c1->SetObjectStat(false);

   doTrigger(triggers[i]);
   TLatex *   tex = new TLatex(0.8,35, translate[tr].c_str());
   tex->SetTextSize(0.090); //0.044
   tex->SetLineWidth(2);
   tex->Draw();

if(!same)
{
   TLatex *   tex = new TLatex(0.0,0.03*(30)+40,"CMS    #sqrt{s} = 7 TeV, L = 3.1 pb^{-1}");
   tex->SetTextSize(0.040); //0.044
   tex->SetLineWidth(2);
   tex->Draw();

   c1->SaveAs((filename()+"_"+tr+".root").c_str());
   c1->SaveAs((filename()+"_"+tr+".pdf").c_str());
   c1->SaveAs((filename()+"_"+tr+".png").c_str());
} 

   }
if(same)
 {
   c2->cd();
   c2->SaveAs((filename()+"_all.root").c_str());
   c2->SaveAs((filename()+"_all.pdf").c_str());
   c2->SaveAs((filename()+"_all.png").c_str());
}
//   f->Close();
 }

 void doTrigger(std::string tr)
 {
/* TCanvas * c1 = new TCanvas(tr.c_str(),tr.c_str(),100,300,500,500);
   c1->SetTicky();
   c1->SetTickx();
   gStyle->SetOptStat(0);

   c1->SetObjectStat(false);*/
   TLegend * l = legend();
//   l->AddEntry((TObject *)0,tr.c_str());
   std::string options="";
   float mx=-1e99,mi=1e99;
   TAxis * y_axis=0;
   TAxis * x_axis=0;
   for(int i=0;i< samples.size();i++)
   {
     TH1 & d = * draw(options+samples[i].ds.def,samples[i],tr);
     l->AddEntry(&d,samples[i].label.c_str(),samples[i].ds.legend.c_str());

     if(d.GetMaximum() > mx) mx=d.GetMaximum();
     if(d.GetMaximum() < mi) mi=d.GetMinimum();

     if(options == "") 
      {
        y_axis=d.GetYaxis();
        x_axis=d.GetXaxis();
        options="same"; 
      }
   }
   mx*=1.2;
   x_axis->SetTitle( (titleX()+"      ").c_str());
   y_axis->SetTitle(titleY().c_str());
   x_axis->SetTitleOffset(0.4);
   x_axis->SetLabelSize(0.1);//0.1
   x_axis->SetTitleSize(0.1);
   y_axis->SetTitleOffset(0.9);
   y_axis->SetLabelOffset(0.010);
   y_axis->SetLabelSize(0.10);
//    y_axis->SetTitleSize(0.08);


//    y_axis->SetTitleSize(0.071); //added for GeV

   mi=10;
   mx=40;
   y_axis->SetRangeUser(10,mx);
   x_axis->SetRangeUser(0,3.9);
   l->Draw("same");
/* c1->SaveAs((filename()+"_"+tr+".root").c_str());
   c1->SaveAs((filename()+"_"+tr+".pdf").c_str());
   c1->SaveAs((filename()+"_"+tr+".png").c_str());*/
 }


 
};

struct PlotterPtVsDr : public PlotterByTriggers
{
  PlotterPtVsDr(const char * h = "Vert_ptrecsoftVSdR_2b",std::string ytitle=""/*p_{T}^{rec} softer B candidate (GeV)"*/): histoName(h), ty(ytitle){}
 virtual TH1 * draw(std::string options,Sample s,std::string trigger)
 {
   TH2F * h =(TH2F*) s.get(trigger,histoName.c_str());
   TProfile * p = h->ProfileX((s.name+"_"+trigger+"_"+histoName).c_str());
   p->SetLineColor(s.color);
   p->SetFillColor(s.color);
   p->SetMarkerColor(s.color);
   p->SetMarkerStyle(s.ds.marker);
// p->Rebin(10);
   std::string na("rebin_");
   TH1 * p2 = p->Rebin(nbins,(na+p->GetName()).c_str(),a);
   p2->Draw(options.c_str());
   return p2;
 }
// virtual std::string titleX() {return  "\\DeltaR";}
 virtual std::string titleY() {return  ty;}
 virtual std::string filename() {return  "Profile_"+histoName;}
 std::string histoName;
 std::string ty;
};

class PlotterRatioXY : public PlotterRatio
{

  virtual TGraphAsymmErrors * draw(std::string options,Sample s)
 {
    std::cerr << "DRAW" << std::endl;
    //workaround for madgraph having only 3 points
    int n=0; //bins.size();
    for(int i=0;i<bins.size();i++) {
     if(!computePoint(s,bins[i]).invalid) n++;
    }
    TVectorD x(n),y(n),exl(n),exh(n),eyl(n),eyh(n),esysyl(n),esysyh(n);
    int i=0;
    for(int j=0;j<bins.size();j++) {
     GraphPoint p = computePoint(s,bins[j]);
     if(p.invalid) continue;
     TH1F *hx=0;
     hx=((TH1F *)s.get(bins[j].trigger,"Vert_Jet1pt_2b"));
     if(hx==0)  hx=((TH1F *)s.get(bins[j].trigger,"Jet_Jet1pt"));
     x[i]=hx->GetMean();
     exl[i]=hx->GetRMS()/sqrt(hx->GetEntries()/2);
     cout << "DEBUG " << hx->GetEntries() << endl;
     exh[i]=exl[i];
     y[i]=p.y;
     eyl[i]=p.eyl;
     eyh[i]=p.eyh;
     esysyl[i]=p.esysyl;
     esysyh[i]=p.esysyh;
     std::cout << "TEXT " << s.name << " "  << x[i] << " " << y[i] << " " << eyl[i] << " " << esysyl[i] << std::endl;
     i++;
    }
    TGraphAsymmErrors * gr = new TGraphAsymmErrors(x,y,exl,exh,eyl,eyh);
    TGraphAsymmErrors * grSys = new TGraphAsymmErrors(x,y,exl,exh,esysyl,esysyh);
    gr->SetLineColor(s.color);
    gr->SetFillColor(s.color);
    gr->SetMarkerColor(s.color);
    gr->SetFillStyle(s.ds.fill+s.color-1);
    gr->SetMarkerStyle(s.ds.marker);
    std::cout<< options.c_str() << std::endl;
    gr->Draw(options.c_str());
    if(s.drawSystematics)   grSys->Draw((options+"same").c_str());

    return gr;
 }


};

void pt(bool same=false,bool sim=false)
{

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
//   tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.05, "XYZ"); //0.05
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.04, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

  //style modifications:
  gStyle->SetErrorX(0.5);
  gStyle->SetPadTopMargin(0.08);//.005
  gStyle->SetPadBottomMargin(0.10);//0.13
  gStyle->SetMarkerStyle(1); 
  gStyle->SetPalette(1);
  gStyle->SetFuncWidth(2);




translate["HLT_Jet15U"]="p^{Jet}_{T} > 56 GeV";
translate["HLT_Jet30U"]="p^{Jet}_{T} > 84 GeV";
translate["HLT_Jet50U"]="p^{Jet}_{T} > 120 GeV";
DrawStyle mcStyle;
mcStyle.def="E2";
mcStyle.legend="F";
mcStyle.marker=1;

DrawStyle dataStyle;
dataStyle.def="PE1";
dataStyle.legend="P";
dataStyle.marker=20;
std::vector<std::string> triggers;
//triggers.push_back("HLT_L1Jet6U");
triggers.push_back("HLT_Jet15U");
triggers.push_back("HLT_Jet30U");
triggers.push_back("HLT_Jet50U");

std::vector<Sample> samples;
samples.push_back(Sample("Pythia6","PYTHIA",8,mcStyle));
//samples.push_back(Sample("MadgraphBB","MadGraph",kBlue -7,mcStyle));
samples.push_back(Sample("Data","Data ",1,dataStyle));

PlotterPtVsDr pl1;
pl1.triggers = triggers;
pl1.samples = samples;
pl1.doAll("outtest.root",same);
if(sim){
PlotterPtVsDr pl2("Vert_ptsimsoftVSdR_2b","p_{T}^{sim} softer B hadron");
pl2.triggers = triggers;
//pl2.triggers.push_back("HLT_L1Jet6U");
pl2.samples.push_back(Sample("MadgraphBB","MadGraph",kBlue+2,mcStyle));
pl2.samples.push_back(Sample("Pythia6","PYTHIA",8,mcStyle));
pl2.doAll("outtest.root",same);
}


}



void ratio()
{

//   gROOT->ProcessLine(".L ~/tdrstyle.C");
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
//   tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.05, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.04, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

  //style modifications:
  gStyle->SetErrorX(0.5);
  gStyle->SetPadTopMargin(0.08);//.005
  gStyle->SetPadBottomMargin(0.10);//0.13
  gStyle->SetMarkerStyle(1); 
  gStyle->SetPalette(1);
  gStyle->SetFuncWidth(2);



DrawStyle mcStyle;
mcStyle.def="E3";
mcStyle.legend="F";
mcStyle.marker=1;

DrawStyle dataStyle;
dataStyle.def="PE1";
dataStyle.legend="P";
dataStyle.marker=20;
std::vector<PtBin> ptbins;
//ptbins.push_back(PtBin("HLT_L1Jet6U",36,0,0,0.21));
ptbins.push_back(PtBin("HLT_Jet15U",72,0,0,0.106,0.13));
ptbins.push_back(PtBin("HLT_Jet30U",106,0,0,0.099,0.13));
ptbins.push_back(PtBin("HLT_Jet50U",150,0,0,0.083,0.13));

std::vector<Sample> samples;
samples.push_back(Sample("Pythia6","PYTHIA",8,mcStyle,"Vert_dR_efficden",false));
samples.push_back(Sample("MadgraphBBMC","MadGraph",kBlue+2,mcStyle,"Vert_dR_efficden",false));
samples.push_back(Sample("Data","Data ",1,dataStyle,"Vert_dR_2b_CORR",true));

///scratch/wehrlilu/BBCORR/MACROS/SUMFILES/DATA/histo-Data-HLT_Jet15U-BBCorr_scaledEff-total_ETACUTOPEN.root
//samples[2].prefix=("/scratch/wehrlilu/BBCORR/MACROS/SUMFILES/DATA/histo");
//samples[2].postfix=("BBCorr_scaledEff-total_ETACUTOPEN.root");

PlotterRatio pl1;
pl1.bins = ptbins;
pl1.samples = samples;
pl1.doAll("testratio");

PlotterAsym pl2;
pl2.bins = ptbins;
pl2.samples = samples;
pl2.doAll("testasymmetry");

}

void trendRatioPlot()
{

//   gROOT->ProcessLine(".L ~/tdrstyle.C");
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
//   tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.05, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.04, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

  //style modifications:
  gStyle->SetErrorX(0.5);
  gStyle->SetPadTopMargin(0.08);//.005
  gStyle->SetPadBottomMargin(0.10);//0.13
  gStyle->SetMarkerStyle(1); 
  gStyle->SetPalette(1);
  gStyle->SetFuncWidth(2);



// DrawStyle mcStyle;
// mcStyle.def="E3";
// mcStyle.legend="F";
// mcStyle.marker=1;

DrawStyle dataStyle1;
dataStyle1.def="PE2";
dataStyle1.legend="F";
dataStyle1.marker=20;
dataStyle1.fill=3002; 

DrawStyle dataStyle2;
dataStyle2.def="PE2";
dataStyle2.legend="F";
dataStyle2.marker=20;
dataStyle2.fill=3003; 

DrawStyle dataStyle3;
dataStyle3.def="PE2";
dataStyle3.legend="F";
dataStyle3.marker=4;
dataStyle3.fill=3004; 

std::vector<PtBin> ptbins;
//ptbins.push_back(PtBin("HLT_L1Jet6U",36,0,0,0.21));
ptbins.push_back(PtBin("HLT_Jet15U",72,16,12,0.106,0.13));
ptbins.push_back(PtBin("HLT_Jet30U",106,22,14,0.099,0.13));
ptbins.push_back(PtBin("HLT_Jet50U",150,30,30,0.083,0.13));

std::vector<Sample> samples;
// samples.push_back(Sample("Pythia6","PYTHIA",8,mcStyle,"Vert_dR_efficden",false));
// samples.push_back(Sample("MadgraphBBMC","MadGraph",kBlue+2,mcStyle,"Vert_dR_efficden",false));
 samples.push_back(Sample("Data","Runs > 143731",8,dataStyle1,"Vert_dR_2b_CORR",false)); //after
 samples.push_back(Sample("Data","Runs < 141950",kBlue+2,dataStyle2,"Vert_dR_2b_CORR",false)); //before
 samples.push_back(Sample("Data","Runs 141950-143731",1,dataStyle3,"Vert_dR_2b_CORR",false)); //between

///scratch/wehrlilu/BBCORR/MACROS/SUMFILES/DATA/histo-Data-HLT_Jet15U-BBCorr_scaledEff-total_ETACUTOPEN.root
//samples[2].prefix=("/scratch/wehrlilu/BBCORR/MACROS/SUMFILES/DATA/histo");
//samples[2].postfix=("BBCorr_scaledEff-total_ETACUTOPEN.root");
 samples[0].prefix=("/shome/wehrlilu/BBCORR/forUdatingMacrosInCVS/CMSSW_3_8_6/src/UserCode/BbCorrelation/Macros/trendPlot/AFTER143731/histo");
 samples[1].prefix=("/shome/wehrlilu/BBCORR/forUdatingMacrosInCVS/CMSSW_3_8_6/src/UserCode/BbCorrelation/Macros/trendPlot/BEF141950/histo");
 samples[2].prefix=("/shome/wehrlilu/BBCORR/forUdatingMacrosInCVS/CMSSW_3_8_6/src/UserCode/BbCorrelation/Macros/trendPlot/BETWEEN/histo");

PlotterRatio pl1;
pl1.bins = ptbins;
pl1.samples = samples;
pl1.doAll("trendplot");

// PlotterAsym pl2;
// pl2.bins = ptbins;
// pl2.samples = samples;
// pl2.doAll("testasymmetry");

}



void systemPhaseSpace()
{
DrawStyle mcStyle;
mcStyle.def="E3";
mcStyle.legend="F";
mcStyle.marker=1;

DrawStyle dataStyle;
dataStyle.def="PE1";
dataStyle.legend="P";
dataStyle.marker=20;
std::vector<PtBin> ptbins;
//ptbins.push_back(PtBin("HLT_L1Jet6U",36,0,0,0.21));
ptbins.push_back(PtBin("HLT_Jet15U",72,0,0,0.106,0.13));
ptbins.push_back(PtBin("HLT_Jet30U",106,0,0,0.099,0.13));
ptbins.push_back(PtBin("HLT_Jet50U",150,0,0,0.083,0.13));

std::vector<Sample> samples;
samples.push_back(Sample("Data","Data (8GeV cut)",1,dataStyle,"Vert_dR_2b_CORR",true));
samples.push_back(Sample("Data","Data (10GeV cut)",2,dataStyle,"Vert_dR_2b_CORR",true));

samples[1].postfix=("BBCorr_HadrJetMatchPt_OpenHLT_BCand10GeV-total.root");

PlotterRatio pl1;
pl1.bins = ptbins;
pl1.samples = samples;
pl1.doAll("testratio");


}

void systemPV()
{
DrawStyle mcStyle;
mcStyle.def="E3";
mcStyle.legend="F";
mcStyle.marker=1;

DrawStyle dataStyle;
dataStyle.def="PE3L";
dataStyle.legend="F";
dataStyle.marker=20;
dataStyle.fill=3002;
std::vector<PtBin> ptbins;
//ptbins.push_back(PtBin("HLT_L1Jet6U",36,0,0,0.21));
ptbins.push_back(PtBin("HLT_Jet15U",72,0,0,0.,0.));
ptbins.push_back(PtBin("HLT_Jet30U",106,0,0,0.,0.));
ptbins.push_back(PtBin("HLT_Jet50U",150,0,0,0.,0.));

std::vector<Sample> samples;
samples.push_back(Sample("Data","Data (standard 36X)",1,dataStyle,"Vert_dR_2b_CORR",false));
samples.push_back(Sample("Data386","Data (2PV <1cm 38X)",2,dataStyle,"Vert_dR_2b_CORR",false));
samples.push_back(Sample("Data386","Data (2PV >1cm 38X)",4,dataStyle,"Vert_dR_2b_CORR",false));
samples.push_back(Sample("Data386","Data (1PV  38X)",6,dataStyle,"Vert_dR_2b_CORR",false));

samples[1].postfix="BBCorr_HadrJetMatchPt_OpenHLT_no2VCheck_2VL1cm-total.root";
samples[2].postfix="BBCorr_HadrJetMatchPt_OpenHLT_no2VCheck_2VG1cm-total.root";
samples[3].postfix="BBCorr_HadrJetMatchPt_OpenHLT_no2VCheck_1V-total.root";

/*PlotterAsym pl;
pl.bins = ptbins;
pl.samples = samples;
pl.doAll("testasym");
*/
//PlotterRatio pl1;
PlotterRatioXY pl1;
pl1.bins = ptbins;
pl1.samples = samples;
pl1.doAll("testratio");


}

