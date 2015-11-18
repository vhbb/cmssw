triggerObjectCollectionsFull = {
    "caloMet":[["hltMet","","HLT"],"hltMET90","HLT_PFMET90_PFMHT90_IDTight*"]
}
triggerObjectCollectionsOnlyPt = {
    "l1Met":[["hltL1extraParticles","MET"]]
}
triggerObjectCollectionsOnlySize = {
    "caloJets":[["hltAK4CaloJetsCorrectedIDPassed"]],
}

triggerObjectCollectionsOnlySize.update(triggerObjectCollectionsFull)
triggerObjectCollectionsOnlySize.update(triggerObjectCollectionsOnlyPt)

