#include<samples.h>
#include <iostream> 
#include <TCanvas.h>
#include <TRegexp.h>
#include <TLegend.h>
#include <THStack.h>
#include "setTDRStyle.C"
#include <TROOT.h>
#include "customize.h"

void plottingmacro()
{
 setTDRStyle();
 gROOT->ForceStyle();
 initOptions();

 std::vector<Sample> s = samples();
 Sample data(1,"fake data","S1.root",0,true,1000);

 for(size_t i=0;i< s.size();i++) if(s[i].data) {data=s[i];break;}
 data.file()->ls(); 
 for(size_t i=0;i< s.size();i++) s[i].dump(data.lumi());

 std::vector<std::string> names;

 TList * subs = data.file()->GetListOfKeys();
 for(size_t i=0;i< subs->GetSize();i++)
  {
    TList * objs = ((TDirectoryFile *) data.file()->Get(subs->At(i)->GetName()))->GetListOfKeys();
     for(size_t j=0;j< objs->GetSize();j++)
     {
         names.push_back(subs->At(i)->GetName()+std::string("/")  + objs->At(j)->GetName());
 //      std::cout << subs->At(i)->GetName() << "/"  << objs->At(j)->GetName() << std::endl;
         //TODO: select plots via regexp
     }
    
  }

 for(size_t i = 0 ; i < names.size() ; i++) 
  {
   std::map<std::string,TH1F *> grouped;
   TString n=names[i];
   if(!n.Contains(TRegexp("V.*RegionH.*mu.*HiggsMass"))) continue;
   TCanvas *c = new TCanvas();
   c->SetLogy(true);
   c->SetTitle(names[i].c_str());
   TH1F *hd = ((TH1F*)data.file()->Get(names[i].c_str()));
   Options o=options[names[i]];
   hd->Rebin(o.rebin);
   hd->SetMarkerStyle(21);
   hd->Draw("E1");
   hd->SetYTitle(o.yaxis.c_str());
   THStack * sta = new THStack("sta",hd->GetTitle());
   TLegend * l = new TLegend(o.legendx1,o.legendy1,o.legendx2,o.legendy2); //0.7,0.1,0.9,0.6);
  
   l->AddEntry(hd, "Data","LP");

   for(size_t j=0;j< s.size() ;j++) 
   { 
       if(!s[j].data) 
      {
       TH1F * h = ((TH1F*)s[j].file()->Get(names[i].c_str()));
       h->Scale(s[j].scale(data.lumi()));
       h->SetLineColor(s[j].color);
       h->SetFillColor(s[j].color);
       h->Rebin(options[names[i]].rebin);
       if(grouped.find(s[j].name)==grouped.end()) {
          grouped[s[j].name]=(TH1F *)h->Clone(("_"+names[i]).c_str());
          l->AddEntry(h,s[j].name.c_str(),"F");
       }
       else
       {
        grouped[s[j].name]->Add(h);
       }
       sta->Add(h);
//     h->Draw("same");
      }
   }
   sta->Draw("same");
   hd->Draw("E1same");
   hd->GetYaxis()->SetRangeUser(options[names[i]].min,options[names[i]].max);
   l->Draw();


   std::cout << names[i] << " d: " <<  hd->Integral() << " ";
   THStack * sta2 = new THStack("sta2",hd->GetTitle());
   float tot=0;
   float toterr2=0;

   for(std::map<std::string,TH1F *>::iterator it = grouped.begin(); it != grouped.end();it++)
   {
             std::cout << it->first << " " << it->second->Integral() << " | " ;
             if(it->second->GetEntries() > 0) {
             float er=1.*sqrt(it->second->GetEntries())/it->second->GetEntries()*it->second->Integral();
             toterr2+=er*er;
             }
	     tot+=it->second->Integral();
             sta2->Add(it->second);
   }   
    std::cout << " Tot: " << tot << "+-" << sqrt(toterr2) <<  " SF: " << hd->Integral()/tot << std::endl;
    c = new TCanvas();
    sta2->Draw();
    hd->Draw("E1,same");
    sta2->GetYaxis()->SetRangeUser(options[names[i]].min,options[names[i]].max);
    hd->GetYaxis()->SetRangeUser(options[names[i]].min,options[names[i]].max);
    l->Draw();

   

  }

}
