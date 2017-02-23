#! /usr/bin/env python

# First import basic vhbb so we have the sample available
# Need to set sample.isMC/isData correctly before importing vhbb_combined
# so vhbb_combined knows which modules to schedule 

from vhbb import *

sample.isMC=False
sample.isData=True

from vhbb_combined import *
sample.json="json.txt"
sample.files=[
    #"root://xrootd.ba.infn.it//store/data/Run2015D/BTagCSV/MINIAOD/16Dec2015-v1/50000/00AF8EB4-70AB-E511-9271-00266CFAE7AC.root"
    "F0B09550-7DEA-E611-A445-B8CA3A70A5E8.root"
    #"root://131.169.191.218:1094//store/data/Run2016B/SingleMuon//MINIAOD/PromptReco-v2/000/273/537/00000/6EDA3436-A91F-E611-83EA-02163E014736.root"
    ]

TriggerObjectsAna.triggerObjectInputTag = ('selectedPatTrigger','','PAT')
FlagsAna.processName='RECO'
TrigAna.triggerBits = triggerTableData
L1TriggerAna.processName = 'RECO'

#import code
#code.interact(local=locals())


# and the following runs the process directly 
if __name__ == '__main__':
    from PhysicsTools.HeppyCore.framework.looper import Looper 
    looper = Looper( 'Loop', config, nPrint = 1, nEvents = 20000)

    import time
    import cProfile
    p = cProfile.Profile(time.clock)
    p.runcall(looper.loop)
    p.print_stats()
    looper.write()
