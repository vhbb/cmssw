#ifndef TRIGGERREADER_H
#define TRIGGERREADER_H
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include <string>
#include <map>
#include <TString.h>
#include <TRegexp.h>
#include <regex.h>
#include <stdio.h>

class TriggerReader {
 public:
    TriggerReader(bool passAllEvents=false) : passAll(passAllEvents) {}
    
    void setEvent( fwlite::Event * e) {  ev=e;}
  
    bool accept(const std::string & triggername)
    {
     if(passAll) return true;
     fwlite::Handle<edm::TriggerResults> hTriggerResults;
     hTriggerResults.getByLabel(*ev,"TriggerResults","","HLT");

     
     regex_t regex;
     int reti;


     reti = regcomp(&regex, triggername.c_str(), 0);
//     reti = regcomp(&regex, "HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v.*", 0);
     if( reti ){ std::cerr << "Could not compile regex" << std::endl;}


     std::map<std::string,size_t>::iterator nit;
     if(ev->getRun().run()  != cacheRun[triggername] || nameMap.find(triggername) == nameMap.end())
     {
      std::cout << "new run" << ev->getRun().run() << std::endl;
      cacheRun[triggername]=ev->getRun().run();
      edm::TriggerNames const&  triggerNames = ev->triggerNames(*hTriggerResults);
      std::string oldiname="whatever";
      nit=nameMap.find(triggername);
      if(nit!=nameMap.end() && nit->second < triggerNames.size()) 
        {
           oldiname=triggerNames.triggerName(nit->second);
        }
//      std::cout << "whatever3" << std::endl;
//      std::cout << "reg1 " << regexec(&regex, oldiname.c_str(), 0, NULL, 0)  << std::endl;
      if(nit==nameMap.end() || ( regexec(&regex, oldiname.c_str(), 0, NULL, 0) != 0 ))   //; ! oldiname.Contains( TRegexp(triggername)) ) // changed, search it
        {
         std::cout << "searching " << triggername << std::endl;
         for (unsigned i = 0; i < triggerNames.size(); ++i) {
           std::string iname(triggerNames.triggerName(i));
       
           if( regexec(&regex, iname.c_str(), 0, NULL, 0) == 0 )    //iname.Contains(TRegexp(triggername)) )
             {
               nameMap[triggername]= i;
               std::cout << "FOUND: " << triggerNames.triggerName(i) << " is bit  " << i << std::endl;
             }
         }
       
      }
   } 
    nit=nameMap.find(triggername);
   if(nit==nameMap.end()) 
    {
     std::cout << "ERROR: trigger name not found" << std::endl;
     edm::TriggerNames const&  triggerNames = ev->triggerNames(*hTriggerResults);
    // for (unsigned i = 0; i < triggerNames.size(); ++i)  std::cout << triggerNames.triggerName(i) << " is bit  " << i << "looking for: "<< triggername <<  std::endl;
     nameMap[triggername]=100000000; // meaning not found in this run 
       
     return false; 
    }
   if(nit->second==100000000) return false;
   return hTriggerResults->accept(nit->second);
}

private:
  std::map<std::string,size_t> nameMap;
  std::map<std::string,size_t> cacheRun;
//  unsigned int cacheRun;
  fwlite::Event * ev;
  bool passAll;

};

#endif
