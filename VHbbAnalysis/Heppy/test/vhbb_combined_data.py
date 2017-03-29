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
    "root://cmsxrootd.hep.wisc.edu//store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver3-v1/80000/F0B09550-7DEA-E611-A445-B8CA3A70A5E8.root"
    ]

TrigAna.triggerBits = triggerTableData

#import code
#code.interact(local=locals())


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
