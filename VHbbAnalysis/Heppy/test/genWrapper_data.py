from vhbb import *


from vhbb_combined_data import *
sample.json="json.txt"
sample.isMC=False
sample.isData=True
from PhysicsTools.HeppyCore.framework.looper import Looper
looper = Looper('Loop', config, nPrint=0, nEvents=0)
of = open("tree_data.py", "w")
of.write(looper.analyzers[-1].getPythonWrapper())
of.close()
