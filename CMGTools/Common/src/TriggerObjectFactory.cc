#include "CMGTools/Common/interface/TriggerObjectFactory.h"

#include "DataFormats/Provenance/interface/ProcessHistory.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/Registry.h"

#include <algorithm>
#include <map>
#include <set>

cmg::TriggerObjectFactory::event_ptr cmg::TriggerObjectFactory::create(const edm::Event& iEvent, const edm::EventSetup& iSetup) const{

    const std::string isRealDataString("isRealData");
    const bool isRealData = iEvent.isRealData();

    // Just return empty if no trigger info found.
    cmg::TriggerObjectFactory::event_ptr result(new cmg::TriggerObjectFactory::collection);
    if ( processName_ == "*" || processName_ == "@" )
      return result;

    edm::Handle<edm::TriggerResults> triggerResults;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

    // Use the auto-process name.
    const edm::InputTag triggerResultsLabel = edm::InputTag( triggerResultsLabel_.label(), triggerResultsLabel_.instance(), processName_);

    iEvent.getByLabel(triggerResultsLabel, triggerResults);
    iEvent.getByLabel(triggerObjectsLabel_, triggerObjects);

    const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);
    edm::TriggerNames::Strings const& names = triggerNames.triggerNames();
    
    unsigned int pset = hlt_.prescaleSet(iEvent,iSetup);

    // Store which triggers passed and failed.
    std::map<std::string,bool> triggerMap;
    for(edm::TriggerNames::Strings::const_iterator it = names.begin(); it != names.end(); ++it) {
        const unsigned int i = triggerNames.triggerIndex(*it);
        if(i >= triggerNames.size())
          continue;
        triggerMap[*it] = triggerResults->wasrun(i) && triggerResults->accept(i);
    }

    if(useTriggerObjects_) {
        typedef std::multimap<cmg::TriggerObject,std::string> pathMap;
        pathMap triggerObjectPaths;
        for(std::map<std::string, bool>::const_iterator it = triggerMap.begin();
            it != triggerMap.end(); ++it) {
            // Store the trigger objects.
            if(it->second) {
              // Current trigger passed for this object.
                for(unsigned int index = 0; index < triggerObjects->size(); index++) {
                    pat::TriggerObjectStandAlone sa = triggerObjects->at(index);
                    if(sa.path(it->first, true, false)) {
                        pat::TriggerObjectPtr o(triggerObjects, index);
                        cmg::TriggerObject to(o, it->first);
                        result->push_back(to);
                        triggerObjectPaths.insert(std::make_pair(to, it->first));
                    }
                }
            }
        }
        std::sort(result->begin(), result->end());
        std::reverse(result->begin(), result->end());

        // Filter out the same trigger object from different trigger
        // paths.
        std::set<cmg::TriggerObject> triggers;
        for(collection::iterator it = result->begin(); it != result->end(); ++it) {
            std::pair<std::set<cmg::TriggerObject>::iterator, bool> set_it = triggers.insert(*it);

            // This is not unique so we erase it.
            if(!set_it.second) {
                it = result->erase(it);
                --it;
            } else {

                // Now reset the trigger match.
                for(std::map<std::string, bool>::const_iterator jt = triggerMap.begin();
                    jt != triggerMap.end(); ++jt) {
                    triggerMap[jt->first] = false;
                }

                // This is unique so we set the selection.
                std::pair<pathMap::iterator, pathMap::iterator> ret;
                ret = triggerObjectPaths.equal_range(*it);

                // All of the trigger objects that are equal.
                for(pathMap::iterator kt= ret.first; kt!= ret.second; ++kt) {
                    triggerMap[(*kt).second] = true;
                }

                // Store *all* of the trigger results in each object,
                // but most will be false as not relevant to this
                // object.
                for(std::map<std::string,bool>::const_iterator jt = triggerMap.begin();
                    jt != triggerMap.end(); ++jt) {
                  const int prescale = (pset > 0) ? hlt_.prescaleValue(pset,jt->first) : pset;      
                  it->addSelectionWithPrescale(jt->first, jt->second,prescale);
                }
                it->addSelection(isRealDataString, isRealData);
            }
        }

    } else {
        // Just store a single trigger object and set which triggers
        // passed.
        cmg::TriggerObject to;
        for(std::map<std::string,bool>::const_iterator jt = triggerMap.begin(); jt != triggerMap.end(); ++jt) {
            const int prescale = (pset > 0) ? hlt_.prescaleValue(pset,jt->first) : pset;
            to.addSelectionWithPrescale(jt->first,jt->second,prescale);
            //std::cout << "Prescale: " << jt->first << ":" << prescale << "\t - " << jt->second << " " << to.getPrescale(jt->first) << std::endl;
        }
        to.addSelection(isRealDataString,isRealData);
        result->push_back(to);
    }
    return result;
}

std::string cmg::TriggerObjectFactory::getProcessName(edm::Run& iRun, const std::string& tag) {

    //code taken from PatTriggerEventProducer
    //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatAlgos/plugins/PATTriggerEventProducer.cc?revision=1.15&view=markup

    std::string nameProcess = "*";
    // Determine process name from last run module "tag" in process
    // history of input.
    const edm::ProcessHistory& processHistory( iRun.processHistory() );
    edm::ProcessConfiguration processConfiguration;
    edm::ParameterSet processPSet;
    // Unbroken loop, which relies on time ordering (accepts the last
    // found entry).
    for ( edm::ProcessHistory::const_iterator iHist = processHistory.begin(); iHist != processHistory.end(); ++iHist ) {
      if ( processHistory.getConfigurationForProcess( iHist->processName(), processConfiguration )     &&
           edm::pset::Registry::instance()->getMapped( processConfiguration.parameterSetID(), processPSet ) &&
           processPSet.exists( tag )
         ) {
        nameProcess = iHist->processName();
      }
    }
    // If no name was found default to nothing and hope for the best.
    if(nameProcess == "*") {
         edm::LogError( "cmg::TriggerObjectFactory" ) << "trigger::TriggerEvent product with label '"
            << tag << "' not produced according to process history of input data\n"
            << "No trigger information produced.";
        nameProcess = "@";//mark as error state
    }
    return nameProcess;
}

bool cmg::TriggerObjectFactory::beginRun(edm::Run& iRun, edm::EventSetup const& iSetup) {

    if( (processName_ != "") && ( (processName_ == "*") || (processName_ == "@") ) ) {
        processName_ = getProcessName(iRun,tagTriggerEvent_);
    }
    bool changed = false;
    return hlt_.init(iRun,iSetup,processName_,changed);
}
