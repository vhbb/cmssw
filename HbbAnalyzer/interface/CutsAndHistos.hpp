#ifndef CUTANDHISTOS_H
#define CUTANDHISTOS_H
#include "../interface/VHbbEvent.h"
#include <string>
#include <sstream>
#include <vector>
#include <TROOT.h>
#include <TFile.h>

class VHbbEvent;
//class TFolder;

/*class VHbbEventWithHypothesis : public VHbbEvent
{
 TLorentzVector H;
 TLorentzVector V;
};
*/

class Histos {
public:
   virtual void book(TFile &f, std::string suffix) = 0;
   virtual void fill(const VHbbEvent &, float w) = 0;
};

class Cut {
 public:
   virtual bool pass(const VHbbEvent &) = 0;
   virtual std::string name() = 0;
   virtual bool operator()(const VHbbEvent &e) {return pass(e); }
};


class CutSet : public Cut {
 public:
 void add(Cut *c) {cuts.push_back(c);}  
 bool pass(const VHbbEvent &e) {
  bool result=true;
  for(size_t i=0; i< cuts.size(); i++) 
    if( ! (cuts.at(i)->pass(e)) ) 
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
  
  void process(const VHbbEvent &  e,float w) 
  {
    if(cutset.pass(e))
      for(size_t i=0; i< histos.size(); i++) 
	histos.at(i)->fill(e,w);
  }

  
  CutSet cutset;
  std::vector<Histos *> histos;
    
}; 





#endif
