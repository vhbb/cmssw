from vhbb import *
from VHbbAnalysis.Heppy.AdditionalBTag import AdditionalBTag
btagana=cfg.Analyzer(
    verbose=False,
    class_object=AdditionalBTag,
)
sequence.insert(sequence.index(VHbb),btagana)
from PhysicsTools.Heppy.utils.cmsswPreprocessor import CmsswPreprocessor
preprocessor = CmsswPreprocessor("newbtag.py")
config.preprocessor=preprocessor
if __name__ == '__main__':
    from PhysicsTools.HeppyCore.framework.looper import Looper 
    looper = Looper( 'Loop', config, nPrint = 20, nEvents = 3000)
    import time
    import cProfile
    p = cProfile.Profile(time.clock)
    p.runcall(looper.loop)
    p.print_stats()
#    looper.loop()
    looper.write()
