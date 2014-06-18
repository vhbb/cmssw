#ifndef VHBBPROXY_H
#define VHBBPROXY_H

#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEventAuxInfo.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/TriggerReader.h"


class VHbbProxy {

 public:

  VHbbProxy(){}
 VHbbProxy(const VHbbEvent *vhbbev, const VHbbEventAuxInfo* aux, const std::vector<VHbbCandidate> *cand, TriggerReader * tr):
    iEvent(vhbbev),
      iAuxInfo(aux),  
    iCandidate(cand),iTrigger(tr) {}

    const VHbbEvent *getVHbbEvent() { return iEvent; };
    const VHbbEventAuxInfo *getVHbbEventAuxInfo() { return iAuxInfo; };
    const std::vector<VHbbCandidate> *getVHbbCandidate() { return iCandidate; };
    TriggerReader * trigger() {return iTrigger; }
    
 private:
    
    const VHbbEvent *iEvent;
    const VHbbEventAuxInfo *iAuxInfo;
    const std::vector<VHbbCandidate> *iCandidate;
    TriggerReader * iTrigger;

};

#endif
