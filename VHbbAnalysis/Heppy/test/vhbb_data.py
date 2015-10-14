from vhbb import *

sample.isMC=False
sample.isData=True
sample.files=[
 "root://xrootd.unl.edu//store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/254/212/00000/620024C4-3D45-E511-B851-02163E014328.root"
#/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/162/00000/160C08A3-4227-E511-B829-02163E01259F.root"
]
sample.json="json.txt"


FlagsAna.processName='RECO'
TrigAna.triggerBits = triggerTableData

# and the following runs the process directly 
if __name__ == '__main__':
    from PhysicsTools.HeppyCore.framework.looper import Looper 
    looper = Looper( 'Loop', config, nPrint = 1, nEvents = 1000)

    import time
    import cProfile
    p = cProfile.Profile(time.clock)
    p.runcall(looper.loop)
    p.print_stats()
    looper.write()
