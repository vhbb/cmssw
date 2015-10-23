import itertools
import math

import ROOT
import sys
import os

from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.physicsobjects.PhysicsObject import PhysicsObject


class AdditionalBoost( Analyzer ):

    skip_ca15 = False

    def declareHandles(self):
        super(AdditionalBoost, self).declareHandles()
        
        self.handles['ak08']     = AutoHandle( ("slimmedJetsAK8",""), "std::vector<pat::Jet>")

        self.handles['ak08pruned']        = AutoHandle( ("ak08PFPrunedJetsCHS","","EX"), "std::vector<reco::BasicJet>")
        self.handles['ak08prunedsubjets'] = AutoHandle( ("ak08PFPrunedJetsCHS","SubJets","EX"), "std::vector<reco::PFJet>")

        self.handles['ak08softdropsubjets'] = AutoHandle( ("slimmedJetsAK8PFCHSSoftDropPacked","SubJets"), "std::vector<pat::Jet>")

        self.handles['ak08bbtag'] = AutoHandle( ("slimmedJetsAK8pfBoostedDoubleSecondaryVertexBJetTags","","EX"), 
                                                "edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>")

        self.handles['ak08prunedsubjetbtag'] = AutoHandle( ("ak08PFPrunedJetsCHSpfCombinedInclusiveSecondaryVertexV2BJetTags","","EX"), 
                                                           "edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>")


        self.handles['ak08ipTagInfos']     = AutoHandle( ("slimmedJetsAK8ImpactParameterTagInfos","","EX"), "vector<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo> >")

        self.handles['ak08svTagInfos']     = AutoHandle( ("slimmedJetsAK8pfInclusiveSecondaryVertexFinderTagInfos", "","EX"), "vector<reco::TemplatedSecondaryVertexTagInfo<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo>,reco::VertexCompositePtrCandidate> >")


        if not AdditionalBoost.skip_ca15:
        
            self.handles['ca15ipTagInfos']     = AutoHandle( ("ca15PFJetsCHSImpactParameterTagInfos","","EX"), "vector<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo> >")

            self.handles['ca15svTagInfos']     = AutoHandle( ("ca15PFJetsCHSpfInclusiveSecondaryVertexFinderTagInfos", "","EX"), "vector<reco::TemplatedSecondaryVertexTagInfo<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo>,reco::VertexCompositePtrCandidate> >")

            self.handles['ca15ungroomed']     = AutoHandle( ("ca15PFJetsCHS","","EX"), "std::vector<reco::PFJet>")
            self.handles['ca15trimmed']       = AutoHandle( ("ca15PFTrimmedJetsCHS","","EX"), "std::vector<reco::PFJet>")
            self.handles['ca15softdrop']      = AutoHandle( ("ca15PFSoftdropJetsCHS","","EX"), "std::vector<reco::PFJet>")
            self.handles['ca15softdropz2b1']  = AutoHandle( ("ca15PFSoftdropZ2B1JetsCHS","","EX"), "std::vector<reco::PFJet>")
            self.handles['ca15pruned']        = AutoHandle( ("ca15PFPrunedJetsCHS","","EX"), "std::vector<reco::BasicJet>")
            self.handles['ca15prunedsubjets'] = AutoHandle( ("ca15PFPrunedJetsCHS","SubJets","EX"), "std::vector<reco::PFJet>")

            self.handles['ca15tau1'] = AutoHandle( ("ca15PFJetsCHSNSubjettiness","tau1","EX"), "edm::ValueMap<float>")
            self.handles['ca15tau2'] = AutoHandle( ("ca15PFJetsCHSNSubjettiness","tau2","EX"), "edm::ValueMap<float>")
            self.handles['ca15tau3'] = AutoHandle( ("ca15PFJetsCHSNSubjettiness","tau3","EX"), "edm::ValueMap<float>")

            self.handles['ca15softdropz2b1tau1'] = AutoHandle( ("ca15PFSoftdropZ2B1JetsCHSNSubjettiness","tau1","EX"), "edm::ValueMap<float>")
            self.handles['ca15softdropz2b1tau2'] = AutoHandle( ("ca15PFSoftdropZ2B1JetsCHSNSubjettiness","tau2","EX"), "edm::ValueMap<float>")
            self.handles['ca15softdropz2b1tau3'] = AutoHandle( ("ca15PFSoftdropZ2B1JetsCHSNSubjettiness","tau3","EX"), "edm::ValueMap<float>")

            self.handles['httCandJets']  = AutoHandle( ("looseOptRHTT","","EX"), "std::vector<reco::BasicJet>")
            self.handles['httCandInfos'] = AutoHandle( ("looseOptRHTT","","EX"), "vector<reco::HTTTopJetTagInfo>")


            self.handles['httSubjetBtags'] = AutoHandle( ("looseOptRHTTpfCombinedInclusiveSecondaryVertexV2BJetTags","","EX"), 
                                                         "edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>")
            

            self.handles['ca15bbtag'] = AutoHandle( ("ca15PFJetsCHSpfBoostedDoubleSecondaryVertexBJetTags","","EX"), 
                                                    "edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>")

            self.handles['ca15prunedsubjetbtag'] = AutoHandle( ("ca15PFPrunedJetsCHSpfCombinedInclusiveSecondaryVertexV2BJetTags","","EX"), 
                                                               "edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>")

    def process(self, event):
                
        run = event.input.eventAuxiliary().id().run()
        lumi = event.input.eventAuxiliary().id().luminosityBlock()
        eventId = event.input.eventAuxiliary().id().event()
        
        self.readCollections( event.input )

        ######## 
        # AK8 Jets from MiniAOD + Subjet btags
        ########

        setattr(event, "ak08", map(PhysicsObject, self.handles["ak08"].product()))

        setattr(event, "ak0softdropsubjets", map(PhysicsObject, self.handles["ak08softdropsubjets"].product()))

        # bb-tag Output
        newtags =  self.handles['ak08bbtag'].product()

        # Loop over jets                        
        for ij, jet in enumerate(getattr(event, "ak08")):

            # Fill bb-tag
            for i in xrange(len(newtags)) :
                if jet.physObj == newtags.key(i).get():
                    jet.bbtag = newtags.value(i)

        # end of loop over jets



        ######## 
        # Ungroomed Fatjets + NSubjettiness + Hbb Tagging
        ########

        for prefix in ["ca15"]:

            if AdditionalBoost.skip_ca15 and ("ca15" in prefix):
                continue

            # N-Subjettiness
            tau1 = self.handles[prefix+'tau1'].product()
            tau2 = self.handles[prefix+'tau2'].product()
            tau3 = self.handles[prefix+'tau3'].product()

            # bb-tag Output
            newtags =  self.handles[prefix+'bbtag'].product()
                
            # Four Vector
            setattr(event, prefix+"ungroomed", map(PhysicsObject, self.handles[prefix+'ungroomed'].product()))

            # Loop over jets                        
            for ij, jet in enumerate(getattr(event, prefix+"ungroomed")):

                # Fill N-Subjettiness
                jet.tau1 = tau1.get(ij)
                jet.tau2 = tau2.get(ij)
                jet.tau3 = tau3.get(ij)

                # Fill bb-tag
                for i in xrange(len(newtags)) :
                    if jet.physObj == newtags.key(i).get():
                        jet.bbtag = newtags.value(i)
                                    
            # end of loop over jets

        ######## 
        # Softdrop Fatjets + NSubjettiness
        ########

        for fj_name in ["ca15softdropz2b1"]:

            if AdditionalBoost.skip_ca15 and ("ca15" in fj_name):
                continue
                
            # Set the four-vector
            setattr(event, fj_name, map(PhysicsObject, self.handles[fj_name].product()))

            # N-Subjettiness
            tau1 = self.handles[fj_name+'tau1'].product()
            tau2 = self.handles[fj_name+'tau2'].product()
            tau3 = self.handles[fj_name+'tau3'].product()

            # Loop over jets                        
            for ij, jet in enumerate(getattr(event, fj_name)):

                # Fill N-Subjettiness
                jet.tau1 = tau1.get(ij)
                jet.tau2 = tau2.get(ij)
                jet.tau3 = tau3.get(ij)
                                    
            # end of loop over jets

                                                                
        ######## 
        # Groomed Fatjets
        ########

        for fj_name in ['ak08pruned', 'ca15trimmed', 'ca15softdrop', 'ca15pruned']:
            
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
            sjbtags = self.handles['httSubjetBtags'].product()

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

                # Get the correct b-tag
                for ib in xrange(0, len(sjbtags)) :
                    if  sj_w1 == sjbtags.key(ib).get():
                        event.httCandidates[i].sjW1btag = sjbtags.value(ib)


                event.httCandidates[i].sjW2pt   = sj_w2.pt()  
                event.httCandidates[i].sjW2eta  = sj_w2.eta() 
                event.httCandidates[i].sjW2phi  = sj_w2.phi() 
                event.httCandidates[i].sjW2mass = sj_w2.mass()


                # Get the correct b-tag
                for ib in xrange(0, len(sjbtags)) :
                    if  sj_w2 == sjbtags.key(ib).get():
                        event.httCandidates[i].sjW2btag = sjbtags.value(ib)

                event.httCandidates[i].sjNonWpt   = sj_nonw.pt()  
                event.httCandidates[i].sjNonWeta  = sj_nonw.eta() 
                event.httCandidates[i].sjNonWphi  = sj_nonw.phi() 
                event.httCandidates[i].sjNonWmass = sj_nonw.mass()

                # Get the correct b-tag
                for ib in xrange(0, len(sjbtags)) :
                    if  sj_nonw == sjbtags.key(ib).get():
                        event.httCandidates[i].sjNonWbtag = sjbtags.value(ib)

        return True


