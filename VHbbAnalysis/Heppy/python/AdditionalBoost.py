import itertools
import math

import ROOT
import sys
import os

from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.physicsobjects.PhysicsObject import PhysicsObject

# Fastjet-Contrib is not in the path per default.
# We need it for n-subjettiness recalculation
os.environ['ROOT_INCLUDE_PATH'] = "/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/fastjet-contrib/1.014-odfocd/include:" + os.environ['ROOT_INCLUDE_PATH']

ROOT.gSystem.Load("libfastjet")
ROOT.gSystem.Load("libfastjetcontribfragile")
ROOT.gSystem.Load("libRecoBTagSecondaryVertex")

ROOT.gInterpreter.ProcessLine('#include "fastjet/contrib/Njettiness.hh"')

#ROOT.gInterpreter.ProcessLine('#include "RecoBTag/SecondaryVertex/interface/MvaBoostedDoubleSecondaryVertexEstimator.h"')

# Helper function as there seem to be no python bindings for
# DataFormats/Math/interface/deltaR.h
def deltaR2(a, b):    
    phi1 = float(a.phi())
    eta1 = float(a.eta())

    phi2 = float(b.phi())
    eta2 = float(b.eta())
    
    return pow(abs(abs(phi1 - phi2) - math.pi) - math.pi, 2) + pow(eta1 - eta2, 2)

njettiness_08 = ROOT.fastjet.contrib.Njettiness(ROOT.fastjet.contrib.OnePass_KT_Axes(), 
                                                  ROOT.fastjet.contrib.NormalizedMeasure(1.0, 0.8))

njettiness_15 = ROOT.fastjet.contrib.Njettiness(ROOT.fastjet.contrib.OnePass_KT_Axes(), 
                                                  ROOT.fastjet.contrib.NormalizedMeasure(1.0, 1.5))

# Only needed if we want the re-calculate the bb-tag for cross-checking
#mvaID_ca15 = ROOT.MvaBoostedDoubleSecondaryVertexEstimator("/shome/gregor/VHBB-743/CMSSW_7_4_3_patch1/src/RecoBTag/SecondaryVertex/data/BoostedDoubleSV_CA15_BDT.weights.xml.gz")
#mvaID_ak08 = ROOT.MvaBoostedDoubleSecondaryVertexEstimator("/shome/gregor/VHBB-743/CMSSW_7_4_3_patch1/src/RecoBTag/SecondaryVertex/data/BoostedDoubleSV_AK8_BDT.weights.xml.gz")

# Helper function to calculate Hbb tagging input variables
def calcBBTagVariables(jet, 
                       muonTagInfos, 
                       elecTagInfos, 
                       ipTagInfo, 
                       svTagInfo,
                       njettiness,
                       maxSVDeltaRToJet):
    DEBUG = False

    z_ratio          = -1.
    tau_dot          = -1.
    tau_21           = -1.
    SV_pt_0          = -1.
    SV_mass_0        = -1.
    SV_EnergyRatio_0 = -1.
    SV_EnergyRatio_1 = -1.
    tau_21           = -1.
    contSV           = 0
    vertexNTracks    = 0

    # Re-calculate N-subjettiness independently and get axees
    fjParticles = ROOT.std.vector("fastjet::PseudoJet")()
    for dau in jet.daughterPtrVector():
        if dau.isNonnull() and dau.isAvailable():
            fjParticles.push_back( ROOT.fastjet.PseudoJet(dau.px(), dau.py(), dau.pz(), dau.energy()))

    tau1_new = njettiness.getTau(1, fjParticles)
    tau2_new = njettiness.getTau(2, fjParticles)
    currentAxes = njettiness.currentAxes()

    if not tau1_new == 0:
        tau_21 = tau2_new/tau1_new

    selectedTracks = ipTagInfo.selectedTracks()
    trackSize = len(selectedTracks)
    vertexRef = ipTagInfo.primaryVertex()

    allKinematics = ROOT.reco.TrackKinematics()

    if DEBUG:
        print "trackSize=", trackSize

    for itt in range(trackSize):

        ptrack = ROOT.reco.btag.toTrack(selectedTracks[itt])
        ptrackRef = selectedTracks[itt]

        track_PVweight = 0.
        pcand = ptrackRef.get()
        if (pcand.fromPV() == ROOT.pat.PackedCandidate.PVUsedInFit):
            track_PVweight = 1.

        if (track_PVweight>0.):
            if DEBUG:
                print "Adding track", itt, "with weight", track_PVweight
            allKinematics.add(ptrack, track_PVweight)
    # End of loop over tracks

    jetDir = jet.momentum().Unit()

    if DEBUG:
        print  "JetDir = ", jetDir.X(), jetDir.Y(), jetDir.Z()

    #std::map<double, size_t> VTXmass;
    # Just use a Python dict insteas of map

    VTXmass = {}

    for vtx in range(svTagInfo.nVertices()):

        vertexNTracks += (svTagInfo.secondaryVertex(vtx)).numberOfSourceCandidatePtrs()
        flightDir = svTagInfo.flightDirection(vtx);
              
        if DEBUG:
            print "vtx = ", vtx, "flightDir = ", flightDir.x(), flightDir.y(), flightDir.z()

        if ( deltaR2(flightDir, jetDir)<(maxSVDeltaRToJet*maxSVDeltaRToJet)):
            if DEBUG:
                print "pass"
            contSV += 1
            VTXmass[svTagInfo.secondaryVertex(vtx).p4().mass()]=vtx

    # End of loop over vertices

    cont=0

    # GlobalVector is defined DataFormats/GeometryVector/interface/GlobalVector.h
    # typedef Global3DVector GlobalVector;
    # typedef Vector3DBase< float, GlobalTag>    Global3DVector;
    # Since there is no dict for GlobalVector build it directly from the template
    #
    # Commented out constructors to prevent mem-leak
    # flightDir_0 = ROOT.Vector3DBase(float, ROOT.GlobalTag)()
    # flightDir_1 = ROOT.Vector3DBase(float, ROOT.GlobalTag)()
    # SV_p4_0 = ROOT.reco.Candidate.LorentzVector()
    # SV_p4_1 = ROOT.reco.Candidate.LorentzVector()

    for mass in reversed(sorted(VTXmass.keys())):
        cont += 1

        index = VTXmass[mass]

        vertex = svTagInfo.secondaryVertex(index)

        if DEBUG:
            print mass, index

        vtxKinematics = ROOT.reco.TrackKinematics()

        # vertexKinematics
        tracks = vertex.daughterPtrVector()
        for track in tracks:
            mytrack = track.bestTrack()
            vtxKinematics.add(mytrack, 1.0)

        allSum = allKinematics.weightedVectorSum()
        vertexSum = vtxKinematics.weightedVectorSum()

        if cont==1:
            SV_mass_0 = vertex.p4().mass()
            if not allSum.E() == 0:
                SV_EnergyRatio_0 = vertexSum.E() / allSum.E()
            SV_pt_0 = vertex.p4().pt()
            flightDir_0 = svTagInfo.flightDirection(index)
            SV_p4_0 = vertex.p4()

            if deltaR2(flightDir_0, currentAxes[1]) < deltaR2(flightDir_0, currentAxes[0]):
                tau_dot = (currentAxes[1].px()*flightDir_0.x()+currentAxes[1].py()*flightDir_0.y()+currentAxes[1].pz()*flightDir_0.z())/(math.sqrt(currentAxes[1].modp2())*flightDir_0.mag())
            else:
                tau_dot = (currentAxes[0].px()*flightDir_0.x()+currentAxes[0].py()*flightDir_0.y()+currentAxes[0].pz()*flightDir_0.z())/(math.sqrt(currentAxes[0].modp2())*flightDir_0.mag())                            
        # end cont==1

        if cont==2:         
            if not allSum.E() == 0:
                SV_EnergyRatio_1 = vertexSum.E() / allSum.E()
            flightDir_1 = svTagInfo.flightDirection(index)
            SV_p4_1 = vertex.p4()
            z_ratio = math.sqrt(deltaR2(flightDir_0,flightDir_1))*SV_pt_0/(SV_p4_0+SV_p4_1).mass()
            break
        # end cont==2
    # end of loop over vertex/mass

    nSM = muonTagInfos.leptons()
    nSE = elecTagInfos.leptons()
    nSL = nSM + nSE

    if DEBUG:
        print "SL SM SE", nSL, nSM, nSE

    # Lepton pT-rel and IP2D
    PFLepton_ptrel = -1.
    PFLepton_IP2D  = -1.

    # PFMuon information
    for i_lepton in range(muonTagInfos.leptons()):

        PFMuon_ptrel  = muonTagInfos.properties(i_lepton).ptRel

        if (PFMuon_ptrel > PFLepton_ptrel ):
            PFLepton_ptrel = PFMuon_ptrel
            PFLepton_IP2D  = muonTagInfos.properties(i_lepton).sip2d
            # End of loop over leptons

    # PFElectron information
    for i_lepton in range(elecTagInfos.leptons()):

        PFElectron_ptrel  = elecTagInfos.properties(i_lepton).ptRel

        if (PFElectron_ptrel > PFLepton_ptrel ):
            PFLepton_ptrel = PFElectron_ptrel
            PFLepton_IP2D  = elecTagInfos.properties(i_lepton).sip2d
    # End of loop over leptons

    #jet.bbtag_recalc = mvaID.mvaValue(PFLepton_ptrel, 
    #                                  z_ratio, 
    #                                  tau_dot, 
    #                                  SV_mass_0, 
    #                                  SV_EnergyRatio_0, 
    #                                  SV_EnergyRatio_1, 
    #                                  PFLepton_IP2D, 
    #                                  tau_21, 
    #                                  nSL, 
    #                                  vertexNTracks)      
    
    # Fill bb tagging input variables
    jet.PFLepton_ptrel       = PFLepton_ptrel     
    jet.z_ratio              = z_ratio            
    jet.tau_dot              = tau_dot            
    jet.SV_mass_0            = SV_mass_0          
    jet.SV_EnergyRatio_0     = SV_EnergyRatio_0   
    jet.SV_EnergyRatio_1     = SV_EnergyRatio_1   
    jet.PFLepton_IP2D        = PFLepton_IP2D      
    jet.tau_21               = tau_21             
    jet.nSL                  = nSL                
    jet.vertexNTracks        = vertexNTracks       
# End of calcBBTagVariables

class AdditionalBoost( Analyzer ):

    skip_ca15 = False

    def declareHandles(self):
        super(AdditionalBoost, self).declareHandles()
        
        self.handles['ak08']     = AutoHandle( ("slimmedJetsAK8","","PAT"), "std::vector<pat::Jet>")

        self.handles['ak08pruned']        = AutoHandle( ("ak08PFPrunedJetsCHS","","EX"), "std::vector<reco::BasicJet>")
        self.handles['ak08prunedsubjets'] = AutoHandle( ("ak08PFPrunedJetsCHS","SubJets","EX"), "std::vector<reco::PFJet>")

        self.handles['ak08softdropsubjets'] = AutoHandle( ("slimmedJetsAK8PFCHSSoftDropPacked","SubJets","PAT"), "std::vector<pat::Jet>")

        self.handles['ak08bbtag'] = AutoHandle( ("slimmedJetsAK8pfBoostedDoubleSecondaryVertexBJetTags","","EX"), 
                                                "edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>")

        self.handles['ak08prunedsubjetbtag'] = AutoHandle( ("ak08PFPrunedJetsCHSpfCombinedInclusiveSecondaryVertexV2BJetTags","","EX"), 
                                                           "edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>")


        self.handles['ak08ipTagInfos']     = AutoHandle( ("slimmedJetsAK8ImpactParameterTagInfos","","EX"), "vector<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo> >")

        self.handles['ak08svTagInfos']     = AutoHandle( ("slimmedJetsAK8pfInclusiveSecondaryVertexFinderTagInfos", "","EX"), "vector<reco::TemplatedSecondaryVertexTagInfo<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo>,reco::VertexCompositePtrCandidate> >")

        self.handles['ak08muonTagInfos']     = AutoHandle( ("slimmedJetsAK8softPFMuonsTagInfos", "","EX"), "vector<reco::TemplatedSoftLeptonTagInfo<edm::Ptr<reco::Candidate> > >")

        self.handles['ak08elecTagInfos']     = AutoHandle( ("slimmedJetsAK8softPFElectronsTagInfos", "","EX"), "vector<reco::TemplatedSoftLeptonTagInfo<edm::Ptr<reco::Candidate> > >")


        if not AdditionalBoost.skip_ca15:
        
            self.handles['ca15ipTagInfos']     = AutoHandle( ("ca15PFJetsCHSImpactParameterTagInfos","","EX"), "vector<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo> >")

            self.handles['ca15svTagInfos']     = AutoHandle( ("ca15PFJetsCHSpfInclusiveSecondaryVertexFinderTagInfos", "","EX"), "vector<reco::TemplatedSecondaryVertexTagInfo<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo>,reco::VertexCompositePtrCandidate> >")

            self.handles['ca15muonTagInfos']     = AutoHandle( ("ca15PFJetsCHSsoftPFMuonsTagInfos", "","EX"), "vector<reco::TemplatedSoftLeptonTagInfo<edm::Ptr<reco::Candidate> > >")

            self.handles['ca15elecTagInfos']     = AutoHandle( ("ca15PFJetsCHSsoftPFElectronsTagInfos", "","EX"), "vector<reco::TemplatedSoftLeptonTagInfo<edm::Ptr<reco::Candidate> > >")

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


            # bb-tag Inputs
            muonTagInfos = self.handles['ak08muonTagInfos'].product()[ij]
            elecTagInfos = self.handles['ak08elecTagInfos'].product()[ij]
            ipTagInfo    = self.handles['ak08ipTagInfos'].product()[ij]
            svTagInfo    = self.handles['ak08svTagInfos'].product()[ij]

            calcBBTagVariables(jet, 
                               muonTagInfos, 
                               elecTagInfos, 
                               ipTagInfo, 
                               svTagInfo,
                               njettiness_08,
                               maxSVDeltaRToJet = 0.7,)

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

                # bb-tag Inputs
                muonTagInfos = self.handles['ca15muonTagInfos'].product()[ij]
                elecTagInfos = self.handles['ca15elecTagInfos'].product()[ij]
                ipTagInfo    = self.handles['ca15ipTagInfos'].product()[ij]
                svTagInfo    = self.handles['ca15svTagInfos'].product()[ij]

                calcBBTagVariables(jet, 
                                   muonTagInfos, 
                                   elecTagInfos, 
                                   ipTagInfo, 
                                   svTagInfo,
                                   njettiness_15,
                                   maxSVDeltaRToJet = 1.3)
                                    
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


