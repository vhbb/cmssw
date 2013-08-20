#ifndef SimBJet_H
#define SimBJet_H

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
typedef  std::vector<reco::Jet> JetCollection;
typedef  edm::Ref<SimBHadronCollection> SimBHadronRef; 
typedef  edm::RefProd<SimBHadronCollection> SimBHadronRefProd; 
typedef  edm::RefVector<SimBHadronCollection> SimBHadronRefVector; 

#endif

