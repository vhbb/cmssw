from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar import deltaR
from copy import deepcopy
from math import *
import itertools
import ROOT
def Boost(self,boost):
   bx=boost.X()
   by=boost.Y()
   bz=boost.Z()
   b2 = bx*bx + by*by + bz*bz; 
   gamma = 1.0 / sqrt(1.0 - b2);
   bp = bx*self.X() + by*self.Y() + bz*self.Z();
   gamma2 =  (gamma - 1.0)/b2 if b2 >0 else  0.0;
   #print gamma2,gamma,bp,self.X(),self.Y(),self.Z(),self.T()
   self.SetXYZT(self.X() + gamma2*bp*bx + gamma*bx*self.T(),
   self.Y() + gamma2*bp*by + gamma*by*self.T(),
   self.Z() + gamma2*bp*bz + gamma*bz*self.T(),
   (gamma*(self.T() + bp)))
   return self
class VHbbAnalyzer( Analyzer ):
    '''Analyze VH events
    '''

    def declareHandles(self):
        super(VHbbAnalyzer, self).declareHandles()
        if getattr(self.cfg_ana,"doSoftActivity", False) :
            self.handles['pfCands'] =  AutoHandle( 'packedPFCandidates', 'std::vector<pat::PackedCandidate>' )
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


    def beginLoop(self,setup):
        super(VHbbAnalyzer,self).beginLoop(setup)
        if "outputfile" in setup.services :
            setup.services["outputfile"].file.cd()
            self.inputCounter = ROOT.TH1F("Count","Count",1,0,2)

    def doSoftActivity(self,event) :
        event.jetsForVBF = [x for x in event.cleanJetsAll if self.cfg_ana.higgsJetsPreSelection(x) ]
        #compute SoftActivity only for events passing VBF selection
        if len(event.jetsForVBF) < 4 or  event.jetsForVBF[0] < 70 or  event.jetsForVBF[1] < 55 or  event.jetsForVBF[2] < 35 or  event.jetsForVBF[3] < 20 :
            return 
        event.jetsForVBF.sort(key=lambda x:x.pt(),reverse=True)
        event.jetsForVBF=event.jetsForVBF[:4]
        event.bJetsForVBF=sorted(event.jetsForVBF,key = lambda jet : jet.btag('combinedInclusiveSecondaryVertexV2BJetTags'), reverse=True)[:2]
        j1=event.bJetsForVBF[0]
        j2=event.bJetsForVBF[1]
        dRbb = deltaR(j1.eta(),j1.phi(),j2.eta(),j2.phi())
        event.pfCands = list(self.handles['pfCands'].product())
        inputs=ROOT.std.vector(ROOT.heppy.ReclusterJets.LorentzVector)() 
        used=[]
        for j in event.jetsForVBF :
            used.extend(j.daughterPtrVector())
        remainingPF = [x for x in event.pfCands if x.charge() != 0 and abs(x.eta()) < 2.5 and  x.pt() > 0.3 and x.fromPV() >=2 and x not in used] 
        etaMin = min(event.jetsForVBF,key=lambda x:x.eta()).eta()+0.4
        etaMax = max(event.jetsForVBF,key=lambda x:x.eta()).eta()-0.4
        dR0=0.4
        for pf in remainingPF :
#             if pf.eta() > etaMin and pf.eta() < etaMax:
                dr1=deltaR(j1.eta(),j1.phi(),pf.eta(),pf.phi())               
                dr2=deltaR(j2.eta(),j2.phi(),pf.eta(),pf.phi())               
                if dr1+dr2 > dRbb + 2*dR0: 
                   inputs.push_back(pf.p4())
        clusterizer=ROOT.heppy.ReclusterJets(inputs,-1,0.4)
        jets=clusterizer.getGrouping(1)
        event.softActivityJets =  [ ROOT.reco.Particle.LorentzVector(p4) for p4 in jets ]
        event.softActivityJets.sort(key=lambda x:x.pt(), reverse=True)


    def makeJets(self,event,b):
	inputs=ROOT.std.vector(ROOT.heppy.ReclusterJets.LorentzVector)()
        event.pfCands = list(self.handles['pfCands'].product())
#        print "original jets pt,eta,phi \t\t",map(lambda x:"%s,%s,%s --"%(x.pt(),x.eta(),x.phi()),event.cleanJets)
#        print "BquarksFromH pt,eta,phi \t\t",map(lambda x:"%s,%s,%s -- "%(x.pt(),x.eta(),x.phi()),event.genbquarksFromH)
#        print "Inv mass",(event.genbquarksFromH[0].p4()+event.genbquarksFromH[1].p4()).M()
#        print "pt from sum ",(event.genbquarksFromH[0].p4()+event.genbquarksFromH[1].p4()).Pt()
#        print "pt from h",event.genHiggsBoson[0].pt()

        copyhb=map(lambda x: deepcopy(x.p4()),event.genbquarksFromH)
        map(lambda x: Boost(x,b),copyhb)
#        print "BquarksFromH(boost) pt,eta,phi,p \t\t",map(lambda x:"%s,%s,%s,%s -- "%(x.pt(),x.eta(),x.phi(),x.P()),copyhb)
#        print "Inv mass (boost)",(copyhb[0]+copyhb[1]).M()
#        print "pt  (boost)",(copyhb[0]+copyhb[1]).Pt()
      # print "boost",b.X(),b.Y(),b.Z(),"phi",b.Phi()
	for pf in event.pfCands :
	     if pf.fromPV() or pf.charge()==0 :
                p4copy=ROOT.heppy.ReclusterJets.LorentzVector(pf.p4())
                bst=Boost(p4copy,b)
       #        if bst.pt() > 20 : 
       #          print "   candidate orig,boost",pf.pt(),bst.pt(),pf.phi(),bst.phi()
		inputs.push_back(bst)
	clusterizer=ROOT.heppy.ReclusterJets(inputs,-1,0.4)
	jets = clusterizer.getGrouping(10)
        #event.jee = list(self.handles['jee'].product())
#        print "Boosted jets:",
	for j in list(jets):
		oldpt=j.pt()
		oldeta=j.eta()
		oldphi=j.phi()
                Boost(j,-b)
 #               if j.pt() > 20 :
 #                   print oldpt,oldeta,oldphi," -> ",j.pt(),j.eta(),j.phi(),"|",
#	print " "
    def doFakeMET(self,event):
	#fake MET from Zmumu
	event.fakeMET = ROOT.reco.Particle.LorentzVector(0.,0.,0.,0.)
	event.fakeMET.sumet = 0
	if event.Vtype == 0 :
		event.fakeMET=event.met.p4() + event.V
                event.fakeMET.sumet = event.met.sumEt() - event.V.pt()

    def doHiggsHighCSV(self,event) :
        #leading csv interpretation
        event.hJetsCSV=sorted(event.jetsForHiggs,key = lambda jet : jet.btag('combinedInclusiveSecondaryVertexV2BJetTags'), reverse=True)[0:2]
        event.aJetsCSV = [x for x in event.cleanJets if x not in event.hJetsCSV]
        event.hjidxCSV=[event.cleanJets.index(x) for x in event.hJetsCSV ]
        event.ajidxCSV=[event.cleanJets.index(x) for x in event.aJetsCSV ]
        event.aJetsCSV+=event.cleanJetsFwd
        event.HCSV = event.hJetsCSV[0].p4()+event.hJetsCSV[1].p4()

    def doHiggsHighPt(self,event) :
        #highest pair interpretations
        event.hJets=list(max(itertools.combinations(event.jetsForHiggs,2), key = lambda x : (x[0].p4()+x[1].p4()).pt() ))
        event.aJets = [x for x in event.cleanJets if x not in event.hJets]
        event.hjidx=[event.cleanJets.index(x) for x in event.hJets ]
        event.ajidx=[event.cleanJets.index(x) for x in event.aJets ]
        event.aJets+=event.cleanJetsFwd
        hJetsByCSV = sorted(event.hJets , key =  lambda jet : jet.btag('combinedInclusiveSecondaryVertexV2BJetTags'), reverse=True)
        event.hjidxDiJetPtByCSV = [event.cleanJets.index(x) for x in hJetsByCSV]
        event.H = event.hJets[0].p4()+event.hJets[1].p4()


    def doHiggs3cj(self,event) :
        event.H3cj = ROOT.reco.Particle.LorentzVector(0.,0.,0.,0.)
        event.hJets3cj = []
        event.aJets3cj = []
	event.hjidx3cj = []         
	event.ajidx3cj = []         
        event.minDr3 = 0         
        #3 jets interpretations, for the closest 3 central jets
        if (len(event.jetsForHiggs) > 2 and event.jetsForHiggs[2] > 15.): 
           aJetsForHiggs = [x for x in event.jetsForHiggs if x not in event.hJets]
           j3=min(aJetsForHiggs, key = lambda x : min( deltaR( x.eta(), x.phi(), event.hJets[0].eta(), event.hJets[0].phi()),deltaR( x.eta(), x.phi(), event.hJets[1].eta(), event.hJets[1].phi() ) ) )
           event.hJets3cj=event.hJets
           event.hJets3cj.append(j3)
#          event.hJets3cj=list(min(itertools.combinations(event.jetsForHiggs,3), key = lambda x : ( x[2]>15 and (deltaR( x[0].eta(), x[0].phi(), x[2].eta(), x[2].phi()) +  deltaR( x[1].eta(), x[1].phi(), x[2].eta(), x[2].phi()) ) ) ))
           event.aJets3cj = [x for x in event.cleanJets if x not in event.hJets3cj]
           event.hjidx3cj=[event.cleanJets.index(x) for x in event.hJets3cj ]
           event.ajidx3cj=[event.cleanJets.index(x) for x in event.aJets3cj ]
           event.aJets3cj+=event.cleanJetsFwd
           event.H3cj = event.hJets3cj[0].p4()+event.hJets3cj[1].p4()+event.hJets3cj[2].p4()
           event.minDr3 = min(deltaR( event.hJets3cj[0].eta(), event.hJets3cj[0].phi(), event.hJets3cj[2].eta(), event.hJets3cj[2].phi()) ,  deltaR( event.hJets3cj[1].eta(), event.hJets3cj[1].phi(), event.hJets3cj[2].eta(), event.hJets3cj[2].phi()))  


    def classifyMCEvent(self,event):
        if self.cfg_comp.isMC:
		event.VtypeSim = -1
		if len(event.genvbosons) == 1:
			# ZtoLL events, same flavour leptons
			if event.genvbosons[0].pdgId()==23 and event.genvbosons[0].numberOfDaughters()>1 and abs(event.genvbosons[0].daughter(0).pdgId()) == abs(event.genvbosons[0].daughter(1).pdgId()):
				if abs(event.genvbosons[0].daughter(0).pdgId()) == 11:
 					#Ztoee
					event.VtypeSim = 1
				if abs(event.genvbosons[0].daughter(0).pdgId()) == 13:
 					#ZtoMuMu
					event.VtypeSim = 0
				if abs(event.genvbosons[0].daughter(0).pdgId()) == 15:
 					#ZtoTauTau
					event.VtypeSim = 5
				if abs(event.genvbosons[0].daughter(0).pdgId()) in [12,14,16]:
					event.VtypeSim = 4
			#WtoLNu events	
			if abs(event.genvbosons[0].pdgId())==24 and event.genvbosons[0].numberOfDaughters()==2:
				if abs(event.genvbosons[0].daughter(0).pdgId()) == 11 and abs(event.genvbosons[0].daughter(1).pdgId()) == 12:
					#WtoEleNu_e
					event.VtypeSim = 3
			if abs(event.genvbosons[0].daughter(0).pdgId()) == 13 and abs(event.genvbosons[0].daughter(1).pdgId()) == 14:
					#WtoMuNu_mu
					event.VtypeSim = 2
			#to be added: WtoTauNu

		if len(event.genvbosons)>1:
			#print 'more than one W/Zbosons?'
			event.VtypeSim = -2
#		if event.VtypeSim == -1:
#			print '===================================='
#			print ' --------- Debug VtypeSim -1 --------'
#			print '# genVbosons: ',len(event.genvbosons), '| #daughters ', event.genvbosons[0].numberOfDaughters()
#			for i in xrange (0, event.genvbosons[0].numberOfDaughters() ) :
#				print 'daughter ',i ,'| pdgId', event.genvbosons[0].daughter(i).pdgId()	

    def classifyEvent(self,event):
	#assign events to analysis (Vtype)
	#enum CandidateType{Zmumu, Zee, Wmun, Wen, Znn,  Zemu, Ztaumu, Ztaue, Wtaun, Ztautau, Zbb, UNKNOWN};
	event.Vtype=-1
        nLep=len(event.selectedLeptons)	
	event.vLeptons=[]
	#WH requires exactly one selected lepton
	wElectrons=[x for x in event.selectedElectrons if self.cfg_ana.wEleSelection(x) ]
	wMuons=[x for x in event.selectedMuons if self.cfg_ana.wMuSelection(x) ]
	zElectrons=[x for x in event.selectedElectrons if self.cfg_ana.zEleSelection(x) ]
	zMuons=[x for x in event.selectedMuons if self.cfg_ana.zMuSelection(x) ]

        zMuons.sort(key=lambda x:x.pt(),reverse=True)
        zElectrons.sort(key=lambda x:x.pt(),reverse=True)
	if len(zMuons) >=  2 :
#              print  zMuons[0].pt()
              if zMuons[0].pt() > self.cfg_ana.zLeadingMuPt :
		    for i in xrange(1,len(zMuons)):
			if zMuons[0].charge()*zMuons[i].charge()<0 :
	                      event.Vtype = 0
			      event.vLeptons =[zMuons[0],zMuons[i]]
			      break
	elif len(zElectrons) >=  2 :
#	    for i in zElectrons[0].electronIDs()  :
#			print i.first,i.second
              if zElectrons[0].pt() > self.cfg_ana.zLeadingElePt :
		    for i in xrange(1,len(zElectrons)):
			if zElectrons[0].charge()*zElectrons[i].charge()<0 :
	                      event.Vtype = 1
			      event.vLeptons =[zElectrons[0],zElectrons[i]]
			      break
	elif len(wElectrons) + len(wMuons) == 1: 
		if abs(event.selectedLeptons[0].pdgId())==13 :
			event.Vtype = 2
			event.vLeptons =event.selectedLeptons
		if abs(event.selectedLeptons[0].pdgId())==11 :
			event.Vtype = 3
			event.vLeptons =event.selectedLeptons
        elif len(zElectrons) + len(zMuons) > 0 :
                event.Vtype = 5 #there are some loose (Z selection) leptons but not matching the W/Z above requirements
	else :
		event.Vtype = 4	#no leptons at all, apply MET cut
		#apply MET cut
		if  event.met.pt() < 80 :
                         event.Vtype = -1


	event.V=sum(map(lambda x:x.p4(), event.vLeptons),ROOT.reco.Particle.LorentzVector(0.,0.,0.,0.))
	if event.Vtype > 1 :	
		event.V+=ROOT.reco.Particle.LorentzVector(event.met.p4().x(),event.met.p4().y(),0,event.met.p4().pt())
        if event.V.Et() > event.V.pt() :
           event.V.goodMt = sqrt(event.V.Et()**2-event.V.pt()**2)
        else :
           event.V.goodMt = -sqrt(-event.V.Et()**2+event.V.pt()**2)

	event.aLeptons = [x for x in event.inclusiveLeptons if x not in event.vLeptons]

	return True

    def fillTauIndices(self,event) :
        for j in event.cleanJetsAll :
            j.tauIdxs = [event.inclusiveTaus.index(x) for x in j.taus if j.taus in event.inclusiveTaus]
        for t in event.inclusiveTaus :
            dRmin = 1.e+3
            t.jetIdx = -1
            for jIdx, j in enumerate(event.cleanJetsAll) :
                dR = None
                if t.isPFTau():
                   dR = deltaR(t.p4Jet().eta(),t.p4Jet().phi(),j.eta(),j.phi())
                else:
                   dR = deltaR(t.pfEssential().p4CorrJet_.eta(),t.pfEssential().p4CorrJet_.phi(),j.eta(),j.phi())
                if dR < 0.3 and dR < dRmin :
                    t.jetIdx = jIdx
                    dRmin = dR

    def initOutputs (self,event) : 
        event.hJets = []
        event.aJets = []
        event.hjidx = []
        event.ajidx = []
        event.hJetsCSV = []
        event.aJetsCSV = []
        event.hjidxCSV = []
        event.ajidxCSV = []
        event.hJets3cj = []
        event.aJets3cj = []
        event.hjidx3cj = []
        event.ajidx3cj = []
        event.aLeptons = []
        event.vLeptons = []
        event.H = ROOT.reco.Particle.LorentzVector(0.,0.,0.,0.)
        event.HCSV = ROOT.reco.Particle.LorentzVector(0.,0.,0.,0.)
        event.H3cj = ROOT.reco.Particle.LorentzVector(0.,0.,0.,0.)
        event.V = ROOT.reco.Particle.LorentzVector(0.,0.,0.,0.)
        event.minDr3=-1
        event.V.goodMt=0
        event.hjidxDiJetPtByCSV = []
        event.softActivityJets=[]
    def process(self, event):
        self.readCollections( event.input )
        self.inputCounter.Fill(1)
        self.initOutputs(event)
#	event.pfCands = self.handles['pfCands'].product()
# 	met = event.met
        #self.addNewBTag(event)
   	self.classifyMCEvent(event)
	self.classifyEvent(event) 
	self.doFakeMET(event)

	#substructure threshold, make configurable
	ssTrheshold = 200.
	# filter events with less than 2 jets with pt 20
        event.jetsForHiggs = [x for x in event.cleanJets if self.cfg_ana.higgsJetsPreSelection(x) ]
	if not   len(event.jetsForHiggs) >= 2 : # and event.jetsForHiggs[1] > 20.) : # or(len(event.cleanJets) == 1 and event.cleanJets[0] > ssThreshold ) ) :
		return self.cfg_ana.passall
        if event.Vtype < 0 and not ( sum(x.pt() > 30 for x in event.jetsForHiggs) >= 4 or sum(x.pt() for x in event.jetsForHiggs[:4]) > 160 ):
                return self.cfg_ana.passall

	self.doHiggsHighCSV(event)
	self.doHiggsHighPt(event)
        self.doHiggs3cj(event)
        self.fillTauIndices(event)
        if getattr(self.cfg_ana,"doSoftActivity", False) :
            self.doSoftActivity(event)

    #    event.jee = list(self.handles['jee'].product())
	#for j in list(jets)[0:3]:
	#	print j.pt(),
	#print " "
	
	#to implement
	#if some threshold: 
	#   self.computeSubStructuresStuff()  
   	
	#self.doIVFHiggs()
	#self.computePullAngle()
	#

	#perhaps in different producers:	
	# LHE weights
	# Trigger weights
	# gen level VH specific info
	# add soft jet info to jets
	# PU weights
	# SIM B hadrons information
	# MET corrections (in MET analyzer)
  	# trigger flags
		
       # Hbalance = ROOT.TLorentzVector()
       # Hbalance.SetPtEtaPhiM(event.V.pt(),0,event.V.phi(),125.)
       # hh=event.genHiggsBoson[0]
       # Hbalance2 = ROOT.TLorentzVector()
       # Hbalance2.SetPtEtaPhiM(hh.pt(),hh.eta(),hh.phi(),125.)
       # print "Hbalance pt,e,phi",Hbalance.Pt(),Hbalance.E(),Hbalance.Phi()
       # print "gen H    pt,e,phi",hh.pt(),hh.p4().E(),hh.phi()
       # print "Hbalance2 pt,e,phi",Hbalance2.Pt(),Hbalance2.E(),Hbalance2.Phi()
       # print "boostVector:",Hbalance.BoostVector().X(),Hbalance.BoostVector().Y(),Hbalance.BoostVector().Z()
       # print "boostVector2:",Hbalance2.BoostVector().X(),Hbalance2.BoostVector().Y(),Hbalance2.BoostVector().Z()
        #print "Unboosted Hbalance",Boost(Hbalance,-Hbalance.BoostVector()).Pt()
#       print "hbalance boost vector",Hbalance.BoostVector()
       # hhh=deepcopy(Hbalance2)
       # Boost(hhh,-Hbalance2.BoostVector())
       # print "unboosted pt",hhh.Pt()
#	self.makeJets(event,-Hbalance2.BoostVector())
        return True


