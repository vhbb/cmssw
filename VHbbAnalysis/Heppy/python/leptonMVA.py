from array import array
from math import *

import ROOT
#import os
#if "/smearer_cc.so" not in ROOT.gSystem.GetLibraries(): 
#    ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/TTHAnalysis/python/plotter/smearer.cc+" % os.environ['CMSSW_BASE']);
#if "/mcCorrections_cc.so" not in ROOT.gSystem.GetLibraries(): 
#    ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/TTHAnalysis/python/plotter/mcCorrections.cc+" % os.environ['CMSSW_BASE']);

def jetLepAwareJEC(lep): # use only if jetAna.calculateSeparateCorrections==True
    p4l = lep.p4()
    l = ROOT.TLorentzVector(p4l.Px(),p4l.Py(),p4l.Pz(),p4l.E())
    if not hasattr(lep.jet,'rawFactor'): return l # if lep==jet (matched to lepton object itself)
    p4j = lep.jet.p4()
    j = ROOT.TLorentzVector(p4j.Px(),p4j.Py(),p4j.Pz(),p4j.E())
    if ((j*lep.jet.rawFactor()-l).Rho()<1e-4): return l # matched to jet containing only the lepton
    j = (j*lep.jet.rawFactor()-l*(1.0/lep.jet.l1corrFactor()))*lep.jet.corrFactor()+l
    return j
def ptRelv2(lep): # use only if jetAna.calculateSeparateCorrections==True
    m = jetLepAwareJEC(lep)
    p4l = lep.p4()
    l = ROOT.TLorentzVector(p4l.Px(),p4l.Py(),p4l.Pz(),p4l.E())
    if ((m-l).Rho()<1e-4): return 0 # lep.jet==lep (no match) or jet containing only the lepton
    return l.Perp((m-l).Vect())

class MVAVar:
    def __init__(self,name,func,corrfunc=None):
        self.name = name
        self.var  = array('f',[0.])
        self.func = func
        self.corrfunc = corrfunc
    def set(self,lep,ncorr): ## apply correction ncorr times
        self.var[0] = self.func(lep)
        if self.corrfunc:
            for i in range(ncorr):
                self.var[0] = self.corrfunc(self.var[0], lep.pdgId(),lep.pt(),lep.eta(),lep.mcMatchId,lep.mcMatchAny)
class MVATool:
    def __init__(self,name,xml,specs,vars):
        self.name = name
        self.reader = ROOT.TMVA.Reader("Silent")
        self.specs = specs
        self.vars  = vars
        for s in specs: self.reader.AddSpectator(s.name,s.var)
        for v in vars:  self.reader.AddVariable(v.name,v.var)
        #print "Would like to load %s from %s! " % (name,xml)
        self.reader.BookMVA(name,xml)
    def __call__(self,lep,ncorr): ## apply correction ncorr times
        for s in self.specs: s.set(lep,ncorr)
        for s in self.vars:  s.set(lep,ncorr)
        return self.reader.EvaluateMVA(self.name)   
class CategorizedMVA:
    def __init__(self,catMvaPairs):
        self.catMvaPairs = catMvaPairs
    def __call__(self,lep,ncorr):
        for c,m in self.catMvaPairs:
            if c(lep): return m(lep,ncorr)
        return -99.

_CommonSpect = [ 
]
_CommonVars = {
 'WithPtV2':[ 
    MVAVar("LepGood_pt",lambda x: x.pt()),
    MVAVar("LepGood_miniRelIsoCharged",lambda x: getattr(x,'miniAbsIsoCharged',-99)/x.pt()), 
    MVAVar("LepGood_miniRelIsoNeutral",lambda x: getattr(x,'miniAbsIsoNeutral',-99)/x.pt()), 
    MVAVar("LepGood_jetPtRelv2", lambda x : ptRelv2(x) if hasattr(x,'jet') else -1),
    MVAVar("LepGood_jetPtRatio := min(LepGood_jetPtRatio_LepAwareJECv2,1.5)", lambda x : min((x.pt()/jetLepAwareJEC(x).Pt() if hasattr(x,'jet') else -1), 1.5)),
    MVAVar("LepGood_jetBTagCSV := max(LepGood_jetBTagCSV,0)", lambda x : max( (x.jet.btag('pfCombinedInclusiveSecondaryVertexV2BJetTags') if hasattr(x.jet, 'btag') else -99) ,0.)),
    MVAVar("LepGood_sip3d",lambda x: x.sip3D() if abs(x.edB(x.PV3D))>0. else 0.),
    MVAVar("LepGood_dxy := log(abs(LepGood_dxy))",lambda x: log(abs(x.dxy()))),
    MVAVar("LepGood_dz  := log(abs(LepGood_dz))", lambda x: log(abs(x.dz()))),
 ],
}

_MuonVars = {
 'WithPtV2': [
    MVAVar("LepGood_segmentCompatibility",lambda x: x.segmentCompatibility()), 
 ],
}

_ElectronVars = {
 'WithPtV2': [
    MVAVar("LepGood_mvaIdPhys14",lambda x: x.mvaRun2("NonTrigPhys14")),
 ]
}


class LeptonMVA:
    def __init__(self, kind, basepath, isMC):
        global _CommonVars, _CommonSpect, _ElectronVars
        #print "Creating LeptonMVA of kind %s, base path %s" % (kind, basepath)
        self._isMC = isMC
        self._kind = kind
        muVars = _CommonVars[kind] + _MuonVars[kind]
        elVars = _CommonVars[kind] + _ElectronVars[kind]
        self.mu = CategorizedMVA([
            ( lambda x: abs(x.eta()) <  1.5, MVATool("BDTG",basepath%"mu_eta_b",_CommonSpect,muVars) ),
            ( lambda x: abs(x.eta()) >= 1.5, MVATool("BDTG",basepath%"mu_eta_e",_CommonSpect,muVars) ),
        ])
        self.el = CategorizedMVA([
            ( lambda x: abs(x.eta()) <  0.8                          , MVATool("BDTG",basepath%"el_eta_cb",_CommonSpect,elVars) ),
            ( lambda x: abs(x.eta()) >= 0.8 and abs(x.eta()) <  1.479, MVATool("BDTG",basepath%"el_eta_fb",_CommonSpect,elVars) ),
            ( lambda x: abs(x.eta()) >= 1.479                        , MVATool("BDTG",basepath%"el_eta_ec",_CommonSpect,elVars) ),
        ])
    def __call__(self,lep,ncorr="auto"):
        if ncorr == "auto": ncorr = 0 # (1 if self._isMC else 0)
        retVal = -99
        if abs(lep.pdgId()) == 11:
            retVal = self.el(lep,ncorr)
        elif abs(lep.pdgId()) == 13:
            retVal = self.mu(lep,ncorr)
        return retVal
    

