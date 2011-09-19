#ifndef HLT_FROM_SOUVIK_H
#define HLT_FROM_SOUVIK_H
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;


class HLTInfoDumperGeneral: public edm::EDProducer 
{
  public:
  HLTInfoDumperGeneral(const edm::ParameterSet&);
   
  private:
  void produce(edm::Event&, const edm::EventSetup&);
  std::string hltPath_;
  edm::InputTag trigTag_, trigEvt_;
  int nFilters_;
  typedef std::vector<edm::InputTag> VInputTag;
  VInputTag filterNames_;
};

HLTInfoDumperGeneral::HLTInfoDumperGeneral(const ParameterSet& cfg)
{
  hltPath_=cfg.getUntrackedParameter<std::string>("HLTPath");
  trigTag_=cfg.getUntrackedParameter<edm::InputTag>("TriggerResults");
  trigEvt_=cfg.getUntrackedParameter<edm::InputTag>("TriggerEvent");
  filterNames_=cfg.getUntrackedParameter<VInputTag>("FilterNames");
  
  produces<unsigned int>("Passed").setBranchAlias("Passed");
  nFilters_=filterNames_.size();
  for (int i=0; i<nFilters_; ++i)
  {
    std::string filterName=filterNames_.at(i).label();
    std::string pT_string=filterName+"PT";
    std::string phi_string=filterName+"Phi";
    std::string eta_string=filterName+"Eta";
    std::string mass_string=filterName+"Mass";
    std::string id_string=filterName+"ID";
    produces<std::vector<float> >(pT_string).setBranchAlias(pT_string);
    produces<std::vector<float> >(phi_string).setBranchAlias(phi_string);
    produces<std::vector<float> >(eta_string).setBranchAlias(eta_string);
    produces<std::vector<float> >(mass_string).setBranchAlias(mass_string);
    produces<std::vector<int> >(id_string).setBranchAlias(id_string);
  }
}

void HLTInfoDumperGeneral::produce(Event &evt, const EventSetup &iSetup) 
{
  auto_ptr<unsigned int> passed(new unsigned int); *passed=0;
  std::vector<float> objPT;
  std::vector<float> objPhi;
  std::vector<float> objEta;
  std::vector<float> objMass;
  std::vector<int> objID;
  std::map<std::string, std::vector<float> > objPT_filter;
  std::map<std::string, std::vector<float> > objPhi_filter;
  std::map<std::string, std::vector<float> > objEta_filter;
  std::map<std::string, std::vector<float> > objMass_filter;
  std::map<std::string, std::vector<int> > objID_filter;
  
  // Extract trigger result information
  Handle<TriggerResults> triggerResults;
  if (!evt.getByLabel(trigTag_, triggerResults)) 
  {
    std::cout<<"TriggerResults does not exist"<<std::endl;
    return;
  }
  evt.getByLabel(trigTag_, triggerResults);
  const edm::TriggerNames &trigNames = evt.triggerNames(*triggerResults);
  for (size_t i=0; i<triggerResults->size(); ++i) 
  {
    std::string trigName = trigNames.triggerName(i);
    if (trigName.find(hltPath_)==0)
    {
      if (triggerResults->accept(i))
      {
        *passed=1;
      }
      break;
    }
  }
  
  objPT_filter.clear();
  objPhi_filter.clear();
  objEta_filter.clear();
  objMass_filter.clear();
  objID_filter.clear();
  if (*passed==1)
  {
    // Now extract triggger object information
    Handle<trigger::TriggerEvent> triggerEvent;
    if (!evt.getByLabel(trigEvt_, triggerEvent))
    {
      std::cout<<"TriggerEvent does not exist"<<std::endl;
      return;
    }
    evt.getByLabel(trigEvt_, triggerEvent);
    const trigger::TriggerObjectCollection &toc=triggerEvent->getObjects();
    
    // Now extract 4-vectors for all the filter objects
    for (int i=0; i<nFilters_; ++i)
    {
      edm::InputTag filterTag=filterNames_.at(i);
      int filter_index=triggerEvent->filterIndex(filterTag);
      if (filter_index!=triggerEvent->sizeFilters())
      {
        const trigger::Keys& keys=triggerEvent->filterKeys(filter_index);
        objPT.clear();
        objPhi.clear();
        objEta.clear();
        objMass.clear();
        objID.clear();
        for (unsigned int k=0; k<keys.size(); k++)
        {
          objPT.push_back(toc[keys[k]].pt());
          objPhi.push_back(toc[keys[k]].phi());
          objEta.push_back(toc[keys[k]].eta());
          objMass.push_back(toc[keys[k]].mass());
          objID.push_back(toc[keys[k]].id());          
        }
        std::string filterName=filterTag.label();
        std::string pT_string=filterName+"PT";
        std::string phi_string=filterName+"Phi";
        std::string eta_string=filterName+"Eta";
        std::string mass_string=filterName+"Mass";
        std::string id_string=filterName+"ID";
        objPT_filter.insert(std::pair<std::string, std::vector<float> >(pT_string, objPT));
        objPhi_filter.insert(std::pair<std::string, std::vector<float> >(phi_string, objPhi));
        objEta_filter.insert(std::pair<std::string, std::vector<float> >(eta_string, objEta));
        objMass_filter.insert(std::pair<std::string, std::vector<float> >(mass_string, objMass));
        objID_filter.insert(std::pair<std::string, std::vector<int> >(id_string, objID));
      }
      else
      {
        std::cout<<"Filter name "<<filterTag<<" not found"<<std::endl;
      }
    }
  }
  
  evt.put(passed, "Passed");
  for (std::map<std::string, std::vector<float> >::iterator i=objPT_filter.begin(); i!=objPT_filter.end(); ++i)
  {
    auto_ptr<std::vector<float> >objPT_ptr(new std::vector<float>);
    for (unsigned j=0; j<(i->second).size(); ++j) objPT_ptr->push_back((i->second).at(j));
    evt.put(objPT_ptr, i->first);
  }
  for (std::map<std::string, std::vector<float> >::iterator i=objPhi_filter.begin(); i!=objPhi_filter.end(); ++i)
  {
    auto_ptr<std::vector<float> >objPhi_ptr(new std::vector<float>);
    for (unsigned j=0; j<(i->second).size(); ++j) objPhi_ptr->push_back((i->second).at(j));
    evt.put(objPhi_ptr, i->first);
  }
  for (std::map<std::string, std::vector<float> >::iterator i=objEta_filter.begin(); i!=objEta_filter.end(); ++i)
  {
    auto_ptr<std::vector<float> >objEta_ptr(new std::vector<float>);
    for (unsigned j=0; j<(i->second).size(); ++j) objEta_ptr->push_back((i->second).at(j));
    evt.put(objEta_ptr, i->first);
  }
  for (std::map<std::string, std::vector<float> >::iterator i=objMass_filter.begin(); i!=objMass_filter.end(); ++i)
  {
    auto_ptr<std::vector<float> >objMass_ptr(new std::vector<float>);
    for (unsigned j=0; j<(i->second).size(); ++j) objMass_ptr->push_back((i->second).at(j));
    evt.put(objMass_ptr, i->first);
  }
  for (std::map<std::string, std::vector<int> >::iterator i=objID_filter.begin(); i!=objID_filter.end(); ++i)
  {
    auto_ptr<std::vector<int> >objID_ptr(new std::vector<int>);
    for (unsigned j=0; j<(i->second).size(); ++j) objID_ptr->push_back((i->second).at(j));
    evt.put(objID_ptr, i->first);
  }
 
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( HLTInfoDumperGeneral );

#endif
