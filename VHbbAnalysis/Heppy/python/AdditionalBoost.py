import itertools

import ROOT
import sys

from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.physicsobjects.PhysicsObject import PhysicsObject

class AdditionalBoost( Analyzer ):

    skip_ca15 = False

    def declareHandles(self):
        super(AdditionalBoost, self).declareHandles()
        
        self.handles['ak08ungroomed']     = AutoHandle( ("ak08PFJetsCHS","","EX"), "std::vector<reco::PFJet>")
        self.handles['ak08softdrop']      = AutoHandle( ("ak08PFSoftdropJetsCHS","","EX"), "std::vector<reco::PFJet>")
        self.handles['ak08pruned']        = AutoHandle( ("ak08PFPrunedJetsCHS","","EX"), "std::vector<reco::BasicJet>")
        self.handles['ak08prunedsubjets'] = AutoHandle( ("ak08PFPrunedJetsCHS","SubJets","EX"), "std::vector<reco::PFJet>")

        self.handles['ak08tau1'] = AutoHandle( ("ak08PFJetsCHSNSubjettiness","tau1","EX"), "edm::ValueMap<float>")
        self.handles['ak08tau2'] = AutoHandle( ("ak08PFJetsCHSNSubjettiness","tau2","EX"), "edm::ValueMap<float>")
        self.handles['ak08tau3'] = AutoHandle( ("ak08PFJetsCHSNSubjettiness","tau3","EX"), "edm::ValueMap<float>")

        self.handles['ak08bbtag'] = AutoHandle( ("ak08PFJetsCHSpfBoostedDoubleSecondaryVertexBJetTags","","EX"), 
                                                "edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>")

        self.handles['ak08prunedsubjetbtag'] = AutoHandle( ("ak08PFPrunedJetsCHSpfCombinedInclusiveSecondaryVertexV2BJetTags","","EX"), 
                                                           "edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>")

        if not AdditionalBoost.skip_ca15:
            self.handles['ca15ungroomed']     = AutoHandle( ("ca15PFJetsCHS","","EX"), "std::vector<reco::PFJet>")
            self.handles['ca15trimmed']       = AutoHandle( ("ca15PFTrimmedJetsCHS","","EX"), "std::vector<reco::PFJet>")
            self.handles['ca15softdrop']      = AutoHandle( ("ca15PFSoftdropJetsCHS","","EX"), "std::vector<reco::PFJet>")
            self.handles['ca15pruned']        = AutoHandle( ("ca15PFPrunedJetsCHS","","EX"), "std::vector<reco::BasicJet>")
            self.handles['ca15prunedsubjets'] = AutoHandle( ("ca15PFPrunedJetsCHS","SubJets","EX"), "std::vector<reco::PFJet>")

            self.handles['ca15tau1'] = AutoHandle( ("ca15PFJetsCHSNSubjettiness","tau1","EX"), "edm::ValueMap<float>")
            self.handles['ca15tau2'] = AutoHandle( ("ca15PFJetsCHSNSubjettiness","tau2","EX"), "edm::ValueMap<float>")
            self.handles['ca15tau3'] = AutoHandle( ("ca15PFJetsCHSNSubjettiness","tau3","EX"), "edm::ValueMap<float>")

            self.handles['httCandJets']  = AutoHandle( ("looseOptRHTT","","EX"), "std::vector<reco::BasicJet>")
            self.handles['httCandInfos'] = AutoHandle( ("looseOptRHTT","","EX"), "vector<reco::HTTTopJetTagInfo>")

            self.handles['ca15bbtag'] = AutoHandle( ("ca15PFJetsCHSpfBoostedDoubleSecondaryVertexBJetTags","","EX"), 
                                                    "edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>")

            self.handles['ca15prunedsubjetbtag'] = AutoHandle( ("ca15PFPrunedJetsCHSpfCombinedInclusiveSecondaryVertexV2BJetTags","","EX"), 
                                                               "edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>")


    def process(self, event):

        self.readCollections( event.input )

        ######## 
        # Ungroomed Fatjets + NSubjettiness
        ########

        for prefix in ["ak08", "ca15"]:

            if AdditionalBoost.skip_ca15 and ("ca15" in prefix):
                continue

            # Four Vector
            setattr(event, prefix+"ungroomed", map(PhysicsObject, self.handles[prefix+'ungroomed'].product()))

            # N-Subjettiness
            tau1 = self.handles[prefix+'tau1'].product()
            tau2 = self.handles[prefix+'tau2'].product()
            tau3 = self.handles[prefix+'tau3'].product()
    
            for i in xrange(0, len(tau1)):
                getattr(event, prefix+"ungroomed")[i].tau1 = tau1.get(i)
                getattr(event, prefix+"ungroomed")[i].tau2 = tau2.get(i)
                getattr(event, prefix+"ungroomed")[i].tau3 = tau3.get(i)

            # Double b-tag
            newtags =  self.handles[prefix+'bbtag'].product()
            for i in xrange(0,len(newtags)) :
                for j in getattr(event, prefix+"ungroomed"):
                    if  j.physObj == newtags.key(i).get():
                        j.bbtag = newtags.value(i)

                                                                
        ######## 
        # Groomed Fatjets
        ########

        for fj_name in ['ak08softdrop', 'ak08pruned', 'ca15trimmed', 'ca15softdrop', 'ca15pruned']:
            
            if AdditionalBoost.skip_ca15 and ("ca15" in fj_name):
                continue

            setattr(event, fj_name, map(PhysicsObject, self.handles[fj_name].product()))


        ######## 
        # Subjets 
        ########

        for fj_name in ['ak08pruned','ca15pruned']:

            if AdditionalBoost.skip_ca15 and ("ca15" in fj_name):
                continue

            setattr(event, fj_name + "subjets", map(PhysicsObject, self.handles[fj_name+"subjets"].product()))
            
            newtags =  self.handles[fj_name+'subjetbtag'].product()
            for i in xrange(0,len(newtags)) :
                for j in getattr(event, fj_name+"subjets"):
                    if  j.physObj == newtags.key(i).get():
                        j.btag = newtags.value(i)


        ######## 
        # HEPTopTagger
        ########

        if not AdditionalBoost.skip_ca15:
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


