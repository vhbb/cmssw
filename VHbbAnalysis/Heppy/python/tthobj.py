
#!/bin/env python
from math import *
import ROOT
from PhysicsTools.Heppy.analyzers.objects.autophobj import *

##------------------------------------------  
## TAU
##------------------------------------------  

tauTypeVHbb = NTupleObjectType("tauTypeVHbb", baseObjectTypes = [ tauType ], variables = [
    NTupleVariable("genMatchType", lambda x : x.genMatchType, int)
])

