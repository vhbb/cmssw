#ifndef CUTANDHISTOS_H
#define CUTANDHISTOS_H
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbProxy.h"
#include <string>
#include <sstream>
#include <vector>
#include <TROOT.h>
#include <TFile.h>

class VHbbEvent;

class Histos {
public:
   virtual void book(TFile &f, std::string suffix) = 0;
   virtual void fill(VHbbProxy &, float w) = 0;
};

class Cut {
 public:
   virtual bool pass(VHbbProxy &) = 0;
   virtual std::string name() = 0;
   virtual bool operator()(VHbbProxy &iProxy) {return pass(iProxy); }
};

class CutSet : public Cut {
 public:
 void add(Cut *c) {cuts.push_back(c);}  
 bool pass(VHbbProxy &iProxy) {
  bool result=true;
  for(size_t i=0; i< cuts.size(); i++) 
    if( ! (cuts.at(i)->pass(iProxy)) ) 
      result=false;
  return result;
 } 
 std::string name() {
   std::stringstream s;
   for(size_t i=0; i< cuts.size(); i++) {
     s << "_" << cuts.at(i)->name();
   }
 return s.str();
 }

private:
 std::vector<Cut *> cuts;
 
};

class CutsAndHistos {
public:
  CutsAndHistos() {}
  CutsAndHistos(CutSet & c, std::vector<Histos *> & h):
    cutset(c),
    histos(h) {}
  
  void book(TFile &f) {
    std::string suffix=cutset.name();
    for(size_t i=0; i< histos.size(); i++) 
      histos.at(i)->book(f,suffix);
  }
  
  void process(VHbbProxy & iProxy,float w)
  {
    if(cutset.pass(iProxy))
      for(size_t i=0; i< histos.size(); i++) 
	histos.at(i)->fill(iProxy,w);
  }

  
  CutSet cutset;
  std::vector<Histos *> histos;
    
}; 





#endif
