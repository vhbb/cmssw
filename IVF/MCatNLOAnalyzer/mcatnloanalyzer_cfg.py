import FWCore.ParameterSet.Config as cms

process = cms.Process("MCATNLO")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/shome/arizzi/CMSSW_3_5_6/src/out0GeV_light.root'
        #'file:/shome/arizzi/CMSSW_3_5_6/src/out30GeV_light.root'
        #'file:/shome/arizzi/CMSSW_3_5_6/src/out80GeV_200M.root'
    )
)
process.TFileService = cms.Service("TFileService",
                                  fileName = cms.string("bbdist_mcatnlo_pthat0_lj120_eta2.root")
                                  #fileName = cms.string("bbdist_pythia_pthattemp.root")
                                  #fileName = cms.string("bbdist_madgraph_pthattemp.root")
                                  #fileName = cms.string("bbdist_mcatnlo_pthat80_20M_pthat.root")
                                  #fileName = cms.string("bbdist_mcatnlo_pthat0_20M_pthat_v2.root")
)
process.bhadrons = cms.EDProducer('MCBHadronProducer',
                                      quarkId = cms.uint32(5)
                                  )

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True) 
    )

process.bbdistmcatnlo = cms.EDAnalyzer('MCatNLOAnalyzer',
                                    simbsrc = cms.InputTag('bhadrons','',''),
                                    minPt = cms.untracked.double(15),   #for b, not for LJ
					#madgraph:
                                    #minPt = cms.untracked.double(30),   #for b, not for LJ
                                    maxEta = cms.untracked.double(2.0), #for b, not for LJ
                                    minPThat = cms.untracked.double(0), #pthat is pt of hardest b quark (id 5) 
                                    maxPThat = cms.untracked.double(30), #pthat is pt of hardest b quark (id 5)
				    ptcutLeadingJet = cms.untracked.double(120) #pt cut for leading jet (84 default)
)

#########################
process.load("RecoJets.JetProducers.ak5GenJets_cfi")
process.genParticlesForJets = cms.EDProducer("InputGenJetsParticleSelector",
    src = cms.InputTag("genParticles"),
    ignoreParticleIDs = cms.vuint32(
         1000022,
         1000012, 1000014, 1000016,
         2000012, 2000014, 2000016,
         1000039, 5100039,
         4000012, 4000014, 4000016,
         9900012, 9900014, 9900016,
         39),
    partonicFinalState = cms.bool(False),
    excludeResonances = cms.bool(True),
    excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
    tausAsJets = cms.bool(False)
)
#########################

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printList = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint = cms.untracked.int32(-1)
)

#MCATNLO
#process.p = cms.Path(process.bhadrons * process.bbdistmcatnlo)
#MadGraph
#process.p = cms.Path(process.printList * process.bhadrons * process.bbdistmcatnlo)
process.p = cms.Path(process.bhadrons * process.bbdistmcatnlo)
#PYTHIA
#process.p = cms.Path(process.genParticlesForJets * process.ak5GenJets * process.bhadrons * process.bbdistmcatnlo)
