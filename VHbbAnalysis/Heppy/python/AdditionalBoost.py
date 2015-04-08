import itertools

import ROOT
import sys

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

        self.handles['httCandJets'] = AutoHandle(   ("looseOptRHTT","","EX"), "std::vector<reco::BasicJet>")
        self.handles['httCandInfos'] = AutoHandle(   ("looseOptRHTT","","EX"), "vector<reco::HTTTopJetTagInfo>")

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

            event.httCandidates[i].fRec = candInfos[i].properties().fRec
            event.httCandidates[i].Ropt = candInfos[i].properties().Ropt
            event.httCandidates[i].RoptCalc = candInfos[i].properties().RoptCalc
            event.httCandidates[i].ptForRoptCalc = candInfos[i].properties().ptForRoptCalc
                                 
            # HTT return the subjet-pair closest to the W-mass as W-subjets
            # Could be improved by b-tagging if we run into a problem
            
            [sj_w1, sj_w2, sj_nonw] = [con.__deref__() for con in candJets[i].getJetConstituents() if not con.isNull()]
            
            event.httCandidates[i].sjW1pt   = sj_w1.pt()
            event.httCandidates[i].sjW1eta  = sj_w1.eta()
            event.httCandidates[i].sjW1phi  = sj_w1.phi()
            event.httCandidates[i].sjW1mass = sj_w1.mass()
                
            event.httCandidates[i].sjW2pt   = sj_w2.pt()  
            event.httCandidates[i].sjW2eta  = sj_w2.eta() 
            event.httCandidates[i].sjW2phi  = sj_w2.phi() 
            event.httCandidates[i].sjW2mass = sj_w2.mass()

            event.httCandidates[i].sjNonWpt   = sj_nonw.pt()  
            event.httCandidates[i].sjNonWeta  = sj_nonw.eta() 
            event.httCandidates[i].sjNonWphi  = sj_nonw.phi() 
            event.httCandidates[i].sjNonWmass = sj_nonw.mass()

        return True


