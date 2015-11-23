from vhbb import *
from PhysicsTools.HeppyCore.framework.looper import Looper
print "making looper"
looper = Looper('Loop', config, nPrint=0, nEvents=0)
of = open("tree.py", "w")
of.write(looper.analyzers[-1].getPythonWrapper())
of.close()
