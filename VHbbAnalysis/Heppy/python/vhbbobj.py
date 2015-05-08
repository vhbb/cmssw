
#!/bin/env python
from math import *
import ROOT
#from CMGTools.TTHAnalysis.signedSip import *
from PhysicsTools.Heppy.analyzers.objects.autophobj import *
import copy

leptonTypeVHbb = NTupleObjectType("leptonTypeVHbb", baseObjectTypes = [ leptonType ], variables = [
    # Loose id 
    NTupleVariable("looseIdSusy", lambda x : x.looseIdSusy if hasattr(x, 'looseIdSusy') else -1, int, help="Loose ID for Susy ntuples (always true on selected leptons)"),
    NTupleVariable("looseIdPOG", lambda x : x.muonID("POG_ID_Loose") if abs(x.pdgId()) == 13 else -1, int, help="Loose ID for Susy ntuples (always true on selected leptons)"),
    # Isolations with the two radia
    NTupleVariable("chargedHadRelIso03",  lambda x : x.chargedHadronIsoR(0.3)/x.pt(), help="PF Rel Iso, R=0.3, charged hadrons only"),
    NTupleVariable("chargedHadRelIso04",  lambda x : x.chargedHadronIsoR(0.4)/x.pt(), help="PF Rel Iso, R=0.4, charged hadrons only"),
    NTupleVariable("eleSieie",    lambda x : x.full5x5_sigmaIetaIeta() if abs(x.pdgId())==11 else -1., help="sigma IEtaIEta for electrons"),
    NTupleVariable("eleDEta",     lambda x : x.deltaEtaSuperClusterTrackAtVtx() if abs(x.pdgId())==11 else -1., help="delta eta for electrons"),
    NTupleVariable("eleDPhi",     lambda x : x.deltaPhiSuperClusterTrackAtVtx() if abs(x.pdgId())==11 else -1., help="delta phi for electrons"),
    NTupleVariable("eleHoE",      lambda x : x.hadronicOverEm() if abs(x.pdgId())==11 else -1., help="H/E for electrons"),
    NTupleVariable("eleMissingHits",      lambda x : x.lostInner() if abs(x.pdgId())==11 else -1., help="Missing hits for electrons"),
    NTupleVariable("eleChi2",      lambda x : x.gsfTrack().normalizedChi2() if abs(x.pdgId())==11 else -1., help="Track chi squared for electrons' gsf tracks"),
    # Extra electron id variables
#    NTupleVariable("convVetoFull", lambda x : (x.passConversionVeto() and x.gsfTrack().trackerExpectedHitsInner().numberOfLostHits() == 0) if abs(x.pdgId())==11 else 1, int, help="Conv veto + no missing hits for electrons, always true for muons."),
    #NTupleVariable("eleMVAId",     lambda x : (x.electronID("POG_MVA_ID_NonTrig") + 2*x.electronID("POG_MVA_ID_Trig")) if abs(x.pdgId()) == 11 else -1, int, help="Electron mva id working point: 0=none, 1=non-trig, 2=trig, 3=both"),
#    NTupleVariable("eleMVAId",     lambda x : (x.electronID("POG_MVA_ID_NonTrig_full5x5") + 2*x.electronID("POG_MVA_ID_Trig_full5x5")) if abs(x.pdgId()) == 11 else -1, int, help="Electron mva id working point (2012, full5x5 shapes): 0=none, 1=non-trig, 2=trig, 3=both"),
#    NTupleVariable("tightCharge",  lambda lepton : ( lepton.isGsfCtfScPixChargeConsistent() + lepton.isGsfScPixChargeConsistent() ) if abs(lepton.pdgId()) == 11 else 2*(lepton.innerTrack().ptError()/lepton.innerTrack().pt() < 0.2), int, help="Tight charge criteria"),
    #NTupleVariable("mvaId",        lambda lepton : lepton.mvaNonTrigV0() if abs(lepton.pdgId()) == 11 else 1, help="EGamma POG MVA ID for non-triggering electrons (as HZZ); 1 for muons"),
#    NTupleVariable("mvaId",         lambda lepton : lepton.mvaNonTrigV0(full5x5=True) if abs(lepton.pdgId()) == 11 else 1, help="EGamma POG MVA ID for non-triggering electrons (as HZZ); 1 for muons"),
#    NTupleVariable("mvaIdTrig",     lambda lepton : lepton.mvaTrigV0(full5x5=True)    if abs(lepton.pdgId()) == 11 else 1, help="EGamma POG MVA ID for triggering electrons; 1 for muons"),
    # Muon-speficic info
    NTupleVariable("nStations",    lambda lepton : lepton.numberOfMatchedStations() if abs(lepton.pdgId()) == 13 else 4, help="Number of matched muons stations (4 for electrons)"),
    NTupleVariable("trkKink",      lambda lepton : lepton.combinedQuality().trkKink if abs(lepton.pdgId()) == 13 else 0, help="Tracker kink-finder"), 
    NTupleVariable("caloCompatibility",      lambda lepton : lepton.caloCompatibility() if abs(lepton.pdgId()) == 13 else 0, help="Calorimetric compatibility"), 
    NTupleVariable("globalTrackChi2",      lambda lepton : lepton.globalTrack().normalizedChi2() if abs(lepton.pdgId()) == 13 and lepton.globalTrack().isNonnull() else 0, help="Global track normalized chi2"), 
    NTupleVariable("nChamberHits", lambda lepton: lepton.globalTrack().hitPattern().numberOfValidMuonHits() if abs(lepton.pdgId()) == 13 and lepton.globalTrack().isNonnull() else -1, help="Number of muon chamber hits (-1 for electrons)"),
    NTupleVariable("isPFMuon", lambda lepton: lepton.isPFMuon() if abs(lepton.pdgId()) == 13 else 0, help="1 if muon passes particle flow ID"),
    NTupleVariable("isGlobalMuon", lambda lepton: lepton.isGlobalMuon() if abs(lepton.pdgId()) == 13 else 0, help="1 if muon is global muon"),
    NTupleVariable("isTrackerMuon", lambda lepton: lepton.isTrackerMuon() if abs(lepton.pdgId()) == 13 else 0, help="1 if muon is tracker muon"),
    NTupleVariable("pixelHits", lambda lepton : lepton.innerTrack().hitPattern().numberOfValidPixelHits() if abs(lepton.pdgId()) == 13 and lepton.innerTrack().isNonnull() else -1, help="Number of pixel hits (-1 for electrons)"),
    # Extra tracker-related id variables
    NTupleVariable("trackerLayers", lambda x : (x.track() if abs(x.pdgId())==13 else x.gsfTrack()).hitPattern().trackerLayersWithMeasurement(), int, help="Tracker Layers"),
    NTupleVariable("pixelLayers", lambda x : (x.track() if abs(x.pdgId())==13 else x.gsfTrack()).hitPattern().pixelLayersWithMeasurement(), int, help="Pixel Layers"),
    # TTH-id related variables
    NTupleVariable("mvaTTH",     lambda lepton : lepton.mvaValue if hasattr(lepton,'mvaValue') else -1, help="Lepton MVA (ttH version)"),
    NTupleVariable("jetPtRatio", lambda lepton : lepton.pt()/lepton.jet.pt() if hasattr(lepton,'jet') else -1, help="pt(lepton)/pt(nearest jet)"),
    NTupleVariable("jetBTagCSV", lambda lepton : lepton.jet.btag('combinedSecondaryVertexBJetTags') if hasattr(lepton,'jet') and hasattr(lepton.jet, 'btag') else -99, help="btag of nearest jet"),
    NTupleVariable("jetDR",      lambda lepton : deltaR(lepton.eta(),lepton.phi(),lepton.jet.eta(),lepton.jet.phi()) if hasattr(lepton,'jet') else -1, help="deltaR(lepton, nearest jet)"),
    # MC-match info
#    NTupleVariable("mcMatchId",  lambda x : x.mcMatchId, int, mcOnly=True, help="Match to source from hard scatter (25 for H, 6 for t, 23/24 for W/Z)"),
#    NTupleVariable("mcMatchAny",  lambda x : x.mcMatchAny, int, mcOnly=True, help="Match to any final state leptons: -mcMatchId if prompt, 0 if unmatched, 1 if light flavour, 2 if heavy flavour (b)"),
#    NTupleVariable("mcMatchTau",  lambda x : x.mcMatchTau, int, mcOnly=True, help="True if the leptons comes from a tau"),
])

##------------------------------------------  
## TAU
##------------------------------------------  

tauTypeVHbb = NTupleObjectType("tauTypeVHbb", baseObjectTypes = [ tauType ], variables = [
    NTupleVariable("idxJetMatch", lambda x : x.jetIdx, int, help="index of the matching jet"),
    NTupleVariable("genMatchType", lambda x : x.genMatchType, int,mcOnly=True, help="..FILLME PLEASE..")
])

##------------------------------------------  
## JET
##------------------------------------------  

jetTypeVHbb = NTupleObjectType("jet",  baseObjectTypes = [ jetType ], variables = [
    NTupleVariable("idxFirstTauMatch", lambda x : x.tauIdxs[0] if len(x.tauIdxs) > 0 else -1, int,help='index of the first matching tau'),
    NTupleVariable("hadronFlavour", lambda x : x.mcFlavour, int,     mcOnly=True, help="match to heavy hadrons"),
    NTupleVariable("btagBDT", lambda x : getattr(x,"btagBDT",-99), help="btag"),
    NTupleVariable("btagProb", lambda x : x.btag('jetProbabilityBJetTags') , help="btag"),
    NTupleVariable("btagBProb", lambda x : x.btag('jetBProbabilityBJetTags') , help="btag"),
    NTupleVariable("btagnew",   lambda x : getattr(x,"btagnew",-2), help="newest btag discriminator"),
    NTupleVariable("btagCSVV0",   lambda x : getattr(x,"btagcsv",-2), help="should be the old CSV discriminator"),
   # NTupleVariable("mcMatchId",    lambda x : x.mcMatchId,   int, mcOnly=True, help="Match to source from hard scatter (25 for H, 6 for t, 23/24 for W/Z)"),
   # NTupleVariable("puId", lambda x : x.puJetIdPassed, int,     mcOnly=False, help="puId (full MVA, loose WP, 5.3.X training on AK5PFchs: the only thing that is available now)"),
   # NTupleVariable("id",    lambda x : x.jetID("POG_PFID") , int, mcOnly=False,help="POG Loose jet ID"),
    NTupleVariable("chHEF", lambda x : x.chargedHadronEnergyFraction(), float, mcOnly = False, help="chargedHadronEnergyFraction (relative to uncorrected jet energy)"),
    NTupleVariable("neHEF", lambda x : x.neutralHadronEnergyFraction(), float, mcOnly = False,help="neutralHadronEnergyFraction (relative to uncorrected jet energy)"),
    NTupleVariable("chEmEF", lambda x : x.chargedEmEnergyFraction(), float, mcOnly = False,help="chargedEmEnergyFraction (relative to uncorrected jet energy)"),
    NTupleVariable("neEmEF", lambda x : x.neutralEmEnergyFraction(), float, mcOnly = False,help="neutralEmEnergyFraction (relative to uncorrected jet energy)"),
    NTupleVariable("chMult", lambda x : x.chargedMultiplicity(), int, mcOnly = False,help="chargedMultiplicity from PFJet.h"),
    NTupleVariable("leadTrackPt", lambda x : x.leadTrackPt() , float, mcOnly = False, help="pt of the leading track in the jet"),   
    NTupleVariable("mcEta",   lambda x : x.mcJet.eta() if getattr(x,"mcJet",None) else 0., mcOnly=True, help="eta of associated gen jet"),
    NTupleVariable("mcPhi",   lambda x : x.mcJet.phi() if getattr(x,"mcJet",None) else 0., mcOnly=True, help="phi of associated gen jet"),
    NTupleVariable("mcM",   lambda x : x.mcJet.p4().M() if getattr(x,"mcJet",None) else 0., mcOnly=True, help="mass of associated gen jet"),
    NTupleVariable("leptonPdgId",   lambda x : x.leptons[0].pdgId() if len(x.leptons) > 0 else -99, mcOnly=False, help="pdg id of the first associated lepton"),
    NTupleVariable("leptonPt",   lambda x : x.leptons[0].pt() if len(x.leptons) > 0 else -99, mcOnly=False, help="pt of the first associated lepton"),
    NTupleVariable("leptonPtRel",   lambda x : ptRel(x.leptons[0].p4(),x.p4()) if len(x.leptons) > 0 else -99, mcOnly=False, help="ptrel of the first associated lepton"),
    NTupleVariable("leptonPtRelInv",   lambda x : ptRel(x.p4(),x.leptons[0].p4()) if len(x.leptons) > 0 else -99, mcOnly=False, help="ptrel Run1 definition of the first associated lepton"),
    NTupleVariable("leptonDeltaR",   lambda x : deltaR(x.leptons[0].p4().eta(),x.leptons[0].p4().phi(),x.p4().eta(),x.p4().phi()) if len(x.leptons) > 0 else -99, mcOnly=False, help="deltaR of the first associated lepton"),
    NTupleVariable("vtxMass",   lambda x : x.userFloat("vtxMass"), mcOnly=False, help="vtxMass from btag"),
    NTupleVariable("vtxNtracks",   lambda x : x.userFloat("vtxNtracks"), mcOnly=False, help="number of tracks at vertex from btag"),
    NTupleVariable("vtxPt",   lambda x : sqrt(x.userFloat("vtxPx")**2 + x.userFloat("vtxPy")**2), mcOnly=False, help="pt of vertex from btag"),
    NTupleVariable("vtx3DSig",   lambda x : x.userFloat("vtx3DSig"), mcOnly=False, help="decay len significance of vertex from btag"),
    NTupleVariable("vtx3DVal",   lambda x : x.userFloat("vtx3DVal"), mcOnly=False, help="decay len of vertex from btag"),
    NTupleVariable("vtxPosX",   lambda x : x.userFloat("vtxPosX"), mcOnly=False, help="X coord of vertex from btag"),
    NTupleVariable("vtxPosY",   lambda x : x.userFloat("vtxPosY"), mcOnly=False, help="Y coord of vertex from btag"), 
    NTupleVariable("vtxPosZ",   lambda x : x.userFloat("vtxPosZ"), mcOnly=False, help="Z coord of vertex from btag"),
   # QG variables:
    NTupleVariable("qgl",   lambda x : getattr(x,'qgl', 0) , float, mcOnly=False,help="QG Likelihood"),
    NTupleVariable("ptd",   lambda x : getattr(x,'ptd', 0), float, mcOnly=False,help="QG input variable: ptD"),
    NTupleVariable("axis2",   lambda x : getattr(x,'axis2', 0) , float, mcOnly=False,help="QG input variable: axis2"),
    NTupleVariable("mult",   lambda x : getattr(x,'mult', 0) , int, mcOnly=False,help="QG input variable: total multiplicity"),
 ])
 
##------------------------------------------  
## FAT JET
##------------------------------------------  

# Four Vector + Nsubjettiness

fatjetType = NTupleObjectType("fatjet",  baseObjectTypes = [ fourVectorType ], variables = [
    NTupleVariable("tau1",  lambda x : x.tau1, help="Nsubjettiness (1 axis)"),
    NTupleVariable("tau2",  lambda x : x.tau2, help="Nsubjettiness (2 axes)"),
    NTupleVariable("tau3",  lambda x : x.tau3, help="Nsubjettiness (3 axes)"),
    ])


##------------------------------------------  
## HEPTopTagger Candidate
##------------------------------------------  

# Four Vector + fW + Rmin + RminExp

httType = NTupleObjectType("htt",  baseObjectTypes = [ fourVectorType ], variables = [
    NTupleVariable("fW",  lambda x : x.fW, help="relative W width"),
    NTupleVariable("Rmin",  lambda x : x.Rmin, help="optimal value of R"),
    NTupleVariable("RminExpected",  lambda x : x.RminExpected, help="expected value of R"),
    ])
   

##------------------------------------------  
## SECONDARY VERTEX CANDIDATE
##------------------------------------------  
  
svType = NTupleObjectType("sv", baseObjectTypes = [ fourVectorType ], variables = [
    NTupleVariable("charge",   lambda x : x.charge(), int),
    NTupleVariable("ntracks", lambda x : x.numberOfDaughters(), int, help="Number of tracks (with weight > 0.5)"),
    NTupleVariable("chi2", lambda x : x.vertexChi2(), help="Chi2 of the vertex fit"),
    NTupleVariable("ndof", lambda x : x.vertexNdof(), help="Degrees of freedom of the fit, ndof = (2*ntracks - 3)" ),
    NTupleVariable("dxy",  lambda x : x.dxy.value(), help="Transverse distance from the PV [cm]"),
    NTupleVariable("edxy", lambda x : x.dxy.error(), help="Uncertainty on the transverse distance from the PV [cm]"),
    NTupleVariable("ip3d",  lambda x : x.d3d.value(), help="3D distance from the PV [cm]"),
    NTupleVariable("eip3d", lambda x : x.d3d.error(), help="Uncertainty on the 3D distance from the PV [cm]"),
    NTupleVariable("sip3d", lambda x : x.d3d.significance(), help="S_{ip3d} with respect to PV (absolute value)"),
    NTupleVariable("cosTheta", lambda x : x.cosTheta, help="Cosine of the angle between the 3D displacement and the momentum"),
    NTupleVariable("jetPt",  lambda x : x.jet.pt() if x.jet != None else 0, help="pT of associated jet"),
    NTupleVariable("jetBTag",  lambda x : x.jet.btag('combinedSecondaryVertexBJetTags') if x.jet != None else -99, help="CSV b-tag of associated jet"),
    NTupleVariable("mcMatchNTracks", lambda x : x.mcMatchNTracks, int, mcOnly=True, help="Number of mc-matched tracks in SV"),
    NTupleVariable("mcMatchNTracksHF", lambda x : x.mcMatchNTracksHF, int, mcOnly=True, help="Number of mc-matched tracks from b/c in SV"),
    NTupleVariable("mcMatchFraction", lambda x : x.mcMatchFraction, mcOnly=True, help="Fraction of mc-matched tracks from b/c matched to a single hadron (or -1 if mcMatchNTracksHF < 2)"),
    NTupleVariable("mcFlavFirst", lambda x : x.mcFlavFirst, int, mcOnly=True, help="Flavour of last ancestor with maximum number of matched daughters"),
    NTupleVariable("mcFlavHeaviest", lambda x : x.mcFlavHeaviest, int, mcOnly=True, help="Flavour of heaviest hadron with maximum number of matched daughters"),
])

heavyFlavourHadronType = NTupleObjectType("heavyFlavourHadron", baseObjectTypes = [ genParticleType ], variables = [
    NTupleVariable("flav", lambda x : x.flav, int, mcOnly=True, help="Flavour"),
    NTupleVariable("sourceId", lambda x : x.sourceId, int, mcOnly=True, help="pdgId of heaviest mother particle (stopping at the first one heaviest than 175 GeV)"),
    NTupleVariable("svMass",   lambda x : x.sv.mass() if x.sv else 0, help="SV: mass"),
    NTupleVariable("svPt",   lambda x : x.sv.pt() if x.sv else 0, help="SV: pt"),
    NTupleVariable("svCharge",   lambda x : x.sv.charge() if x.sv else -99., int, help="SV: charge"),
    NTupleVariable("svNtracks", lambda x : x.sv.numberOfDaughters() if x.sv else 0, int, help="SV: Number of tracks (with weight > 0.5)"),
    NTupleVariable("svChi2", lambda x : x.sv.vertexChi2() if x.sv else -99., help="SV: Chi2 of the vertex fit"),
    NTupleVariable("svNdof", lambda x : x.sv.vertexNdof() if x.sv else -99., help="SV: Degrees of freedom of the fit, ndof = (2*ntracks - 3)" ),
    NTupleVariable("svDxy",  lambda x : x.sv.dxy.value() if x.sv else -99., help="SV: Transverse distance from the PV [cm]"),
    NTupleVariable("svEdxy", lambda x : x.sv.dxy.error() if x.sv else -99., help="SV: Uncertainty on the transverse distance from the PV [cm]"),
    NTupleVariable("svIp3d",  lambda x : x.sv.d3d.value() if x.sv else -99., help="SV: 3D distance from the PV [cm]"),
    NTupleVariable("svEip3d", lambda x : x.sv.d3d.error() if x.sv else -99., help="SV: Uncertainty on the 3D distance from the PV [cm]"),
    NTupleVariable("svSip3d", lambda x : x.sv.d3d.significance() if x.sv else -99., help="SV: S_{ip3d} with respect to PV (absolute value)"),
    NTupleVariable("svCosTheta", lambda x : x.sv.cosTheta if x.sv else -99., help="SV: Cosine of the angle between the 3D displacement and the momentum"),
    NTupleVariable("jetPt",  lambda x : x.jet.pt() if x.jet != None else 0, help="Jet: pT"),
    NTupleVariable("jetBTag",  lambda x : x.jet.btag('combinedSecondaryVertexBJetTags') if x.jet != None else -99, help="CSV b-tag of associated jet"),
])
shiftedMetType= NTupleObjectType("shiftedMetType", baseObjectTypes=[twoVectorType], variables=[
    NTupleVariable("sumEt", lambda x : x.sumEt() ),
])

primaryVertexType = NTupleObjectType("primaryVertex", variables = [
    NTupleVariable("x",    lambda x : x.x()),
    NTupleVariable("y",   lambda x : x.y()),
    NTupleVariable("z",   lambda x : x.z()),
#    NTupleVariable("score",  lambda x : x.mass()), # to be added for 74X
])

genTauJetType = NTupleObjectType("genTauJet", baseObjectTypes = [ genParticleType ], variables = [
    NTupleVariable("decayMode", lambda x : x.decayMode, int, mcOnly=True, help="Generator level tau decay mode"),
])

def ptRel(p4,axis):
    a=ROOT.TVector3(axis.Vect().X(),axis.Vect().Y(),axis.Vect().Z())
    o=ROOT.TLorentzVector(p4.Px(),p4.Py(),p4.Pz(),p4.E())
    return o.Perp(a)

