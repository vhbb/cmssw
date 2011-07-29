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

class NoCut : public Cut {
 public:
   virtual bool pass(VHbbProxy &)  { return true;}
   virtual std::string name() {return "NoCut"; }
};

/// One parameter cut class
class PCut : public Cut
{
 public:
 PCut(float cut) : m_cut(cut) {}
 void setCut(float cut) {m_cut=cut;}

 std::string cutValueString()
  {
     std::stringstream s;
     s << m_cut;
     return s.str();
  } 

 private:
 float m_cut;
 
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


///CutSet of PCut, with scanning functions
class PCutSet : public Cut {
 public:
 void add(PCut *c) {cuts.push_back(c);}
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
 std::vector<PCut *> cuts;

};



class CutsAndHistos {
public:
  CutsAndHistos() {}
  CutsAndHistos(Cut * c, std::vector<Histos *> & h):
    cut(c),
    histos(h) {}
  CutsAndHistos(Cut * c, Histos * h):
    cut(c) {
       histos.push_back(h);
   }
  //TODO: implement destructor for all pointers received
  
  void book(TFile &f) {
    std::string suffix=cut->name();
    for(size_t i=0; i< histos.size(); i++) 
      histos.at(i)->book(f,suffix);
  }
  
  void process(VHbbProxy & iProxy,float w)
  {
    if(cut->pass(iProxy))
      for(size_t i=0; i< histos.size(); i++) 
	histos.at(i)->fill(iProxy,w);
  }

  
  Cut * cut;
  std::vector<Histos *> histos;
    
}; 




#endif
