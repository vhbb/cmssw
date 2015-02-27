from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar import deltaR
from copy import deepcopy
from math import *
import itertools
import ROOT
import array
class AdditionalBTag( Analyzer ):

    def declareHandles(self):
        super(AdditionalBTag, self).declareHandles()
        self.handles['btag'] = AutoHandle( ("combinedInclusiveSecondaryVertexV2BJetTags","","EX"), "edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>")
        self.handles['btagcsv'] = AutoHandle( ("combinedSecondaryVertexBJetTags","","EX"), "edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>")
        reader = ROOT.TMVA.Reader()
        self.Jet_btagCSV = array.array('f',[0]) 
        self.Jet_btagCMVA = array.array('f',[0]) 
        self.Jet_btagProb = array.array('f',[0]) 
        self.Jet_btagBProb = array.array('f',[0]) 
        self.Jet_pt = array.array('f',[0]) 
        self.Jet_eta = array.array('f',[0]) 
        reader.AddVariable("Jet_btagCSV",self.Jet_btagCSV)
        reader.AddVariable("Jet_btagCMVA",self.Jet_btagCMVA)
        reader.AddVariable("Jet_btagProb",self.Jet_btagProb)
        reader.AddVariable("Jet_btagBProb",self.Jet_btagBProb)
        reader.AddSpectator("Jet_pt",self.Jet_pt)
        reader.AddSpectator("Jet_eta",self.Jet_eta)
        reader.BookMVA("BDT","TMVAClassification_BDT.weights.xml")
        self.reader=reader
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
        for j in event.cleanJets :
            self.Jet_btagCSV[0]=j.btag('combinedInclusiveSecondaryVertexV2BJetTags') 
            self.Jet_btagCMVA[0]=j.btag('combinedMVABJetTags') 
            self.Jet_btagProb[0]=j.btag('jetProbabilityBJetTags') 
            self.Jet_btagBProb[0]=j.btag('jetBProbabilityBJetTags') 
            self.Jet_pt[0] = j.pt()
            self.Jet_eta[0] = j.eta()         
            j.btagBDT = self.reader.EvaluateMVA("BDT")
#            print self.Jet_btagCSV, self.Jet_btagCMVA, self.Jet_btagProb, self.Jet_btagBProb, self.Jet_pt, self.Jet_eta,j.btagBDT
        return True


