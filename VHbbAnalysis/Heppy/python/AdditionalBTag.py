from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar import deltaR
from copy import deepcopy
from math import *
import itertools
import ROOT
class AdditionalBTag( Analyzer ):

    def declareHandles(self):
        super(AdditionalBTag, self).declareHandles()
        self.handles['btag'] = AutoHandle( ("combinedInclusiveSecondaryVertexV2BJetTags","","EX"), "edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>")
        self.handles['btagcsv'] = AutoHandle( ("combinedSecondaryVertexBJetTags","","EX"), "edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>")

    def addNewBTag(self,event):
        newtags =  self.handles['btag'].product()
        for i in xrange(0,len(newtags)) :
             for j in event.cleanJets :
                if j.physObj == newtags.key(i).get() :
                    j.btagnew=newtags.value(i)
        newtags =  self.handles['btagcsv'].product()
        for i in xrange(0,len(newtags)) :
             for j in event.cleanJets :
                if j.physObj == newtags.key(i).get() :
                    j.btagcsv=newtags.value(i)


    def process(self, event):
        self.readCollections( event.input )
        self.addNewBTag(event)
        return True


