from vhbb import *
from PhysicsTools.HeppyCore.framework.looper import Looper
looper = Looper('Loop', config)
of = open("tree.py", "w")
of.write(looper.analyzers[-1].getPythonWrapper())
of.close()
