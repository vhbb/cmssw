#ifndef VHBBPROXY_H
#define VHBBPROXY_H

#include "VHbbAnalysis/HbbAnalyzer/interface/VHbbEvent.h"
#include "VHbbAnalysis/HbbAnalyzer/interface/VHbbCandidate.h"

class VHbbProxy {

 public:

  VHbbProxy(){}
  VHbbProxy(const VHbbEvent *ev, const std::vector<VHbbCandidate> *cand):
    iEvent(ev),
    iCandidate(cand) {}

    const VHbbEvent *getVHbbEvent() { return iEvent; };
    const std::vector<VHbbCandidate> *getVHbbCandidate() { return iCandidate; };
    
 private:
    
    const VHbbEvent *iEvent;
    const std::vector<VHbbCandidate> *iCandidate;

};

#endif
