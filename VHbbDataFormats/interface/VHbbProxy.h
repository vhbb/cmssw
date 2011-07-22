#ifndef VHBBPROXY_H
#define VHBBPROXY_H

#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEventAuxInfo.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"

class VHbbProxy {

 public:

  VHbbProxy(){}
 VHbbProxy(const VHbbEvent *ev, const VHbbEventAuxInfo* aux, const std::vector<VHbbCandidate> *cand):
    iEvent(ev),
      iAuxInfo(aux),  
    iCandidate(cand) {}

    const VHbbEvent *getVHbbEvent() { return iEvent; };
    const VHbbEventAuxInfo *getVHbbEventAuxInfo() { return iAuxInfo; };
    const std::vector<VHbbCandidate> *getVHbbCandidate() { return iCandidate; };
    
 private:
    
    const VHbbEvent *iEvent;
    const VHbbEventAuxInfo *iAuxInfo;
    const std::vector<VHbbCandidate> *iCandidate;

};

#endif
