import itertools

import ROOT

from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.physicsobjects.PhysicsObject import PhysicsObject

class AdditionalBoost( Analyzer ):

    def declareHandles(self):
        super(AdditionalBoost, self).declareHandles()
        
        self.handles['ungroomedFatjets'] = AutoHandle( ("ca15PFJetsCHS","","EX"), "std::vector<reco::PFJet>")
        self.handles['trimmedFatjets'] = AutoHandle(   ("ca15PFTrimmedJetsCHS","","EX"), "std::vector<reco::PFJet>")
        
        self.handles['tau1'] = AutoHandle( ("ca15PFJetsCHSNSubjettiness","tau1","EX"), "edm::ValueMap<float>")
        self.handles['tau2'] = AutoHandle( ("ca15PFJetsCHSNSubjettiness","tau2","EX"), "edm::ValueMap<float>")
        self.handles['tau3'] = AutoHandle( ("ca15PFJetsCHSNSubjettiness","tau3","EX"), "edm::ValueMap<float>")

        self.handles['httCandJets'] = AutoHandle(   ("looseMultiRHTT","","EX"), "std::vector<reco::BasicJet>")
        self.handles['httCandInfos'] = AutoHandle(   ("looseMultiRHTT","","EX"), "vector<reco::HTTTopJetTagInfo>")

    def process(self, event):

        self.readCollections( event.input )

        ######## 
        # Trimmed Fatjets
        ########
        event.trimmedFatjets = map(PhysicsObject, self.handles['trimmedFatjets'].product()) 

        ######## 
        # Ungroomed Fatjets + NSubjettiness
        ########
        event.ungroomedFatjets = map(PhysicsObject, self.handles['ungroomedFatjets'].product()) 

        tau1 = self.handles['tau1'].product()
        tau2 = self.handles['tau2'].product()
        tau3 = self.handles['tau3'].product()
    
        for i in xrange(0, len(tau1)):

            event.ungroomedFatjets[i].tau1 = tau1.get(i)
            event.ungroomedFatjets[i].tau2 = tau2.get(i)
            event.ungroomedFatjets[i].tau3 = tau3.get(i)


        ######## 
        # HEPTopTagger
        ########

        candJets = self.handles['httCandJets'].product()
        candInfos = self.handles['httCandInfos'].product()

        event.httCandidates = map(PhysicsObject, candJets) 

        for i in xrange(0, len(candJets)):            
            event.httCandidates[i].fW = candInfos[i].properties().fW
            event.httCandidates[i].Rmin = candInfos[i].properties().Rmin
            event.httCandidates[i].RminExpected = candInfos[i].properties().RminExpected


        return True


