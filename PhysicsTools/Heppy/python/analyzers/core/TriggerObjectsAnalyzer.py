import ROOT

from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.analyzers.core.AutoFillTreeProducer  import NTupleVariable
from PhysicsTools.HeppyCore.utils.deltar import matchObjectCollection, matchObjectCollection3
import PhysicsTools.HeppyCore.framework.config as cfg

""" Check if triggerName is in names using wildcard"""
def FindTrigger(triggerName,names):
    if triggerName=="": return True
    triggerName = triggerName.replace("*","")
    triggerName = triggerName.replace("$","")
    for name in names:
        if triggerName in name: return True
    return False

"""
Example
trColl = triggerCollection(("hltMet","","HLT"),"hltMET90","HLT_PFMET90_PFMHT90_IDTight*")

It gives:
trColl.collectionText ="hltMet::HLT"
trColl.collectionName ="hltMet"
trColl.filterName ="hltMET90"
trColl.path ="HLT_PFMET90_PFMHT90_IDTight*"
"""
class triggerCollection(object):
    def __init__(self, triggerObjectsCfg):
        collectionInfo = triggerObjectsCfg[0]
        collectionLabel = collectionInfo[0] if len(collectionInfo)>0 else ""
        collectionInstance = collectionInfo[1] if len(collectionInfo)>1 else ""
        collectionProcess = collectionInfo[2] if len(collectionInfo)>2 else "HLT"
        self.collectionText = collectionLabel + ":" + collectionInstance + ":" + collectionProcess

        self.collectionName = collectionLabel 
        if collectionInstance != "": self.collectionName = self.collectionName + "_" + collectionInstance
        if collectionProcess != "HLT": self.collectionName = self.collectionName + "_" + collectionProcess

        self.filterName = triggerObjectsCfg[1] if len(triggerObjectsCfg)>1 else ""
        self.path = triggerObjectsCfg[2] if len(triggerObjectsCfg)>2 else ""

class TriggerObjectsAnalyzer( Analyzer ):
    def __init__(self, cfg_ana, cfg_comp, looperName ):
        super(TriggerObjectsAnalyzer,self).__init__(cfg_ana,cfg_comp,looperName)
        triggerObjectsCfgs = getattr(self.cfg_ana,"triggerObjectsCfgs",[])
        self.triggerObjectInputTag = getattr(self.cfg_ana,"triggerObjectInputTag",("","",""))
        self.triggerBitsInputTag = getattr(self.cfg_ana,"triggerBitsInputTag",("","",""))
        self.triggerObjectsInfos = {}
        self.names = None
        for collectionName in triggerObjectsCfgs.keys():
            self.triggerObjectsInfos[collectionName] = triggerCollection(triggerObjectsCfgs[collectionName])
        
    def declareHandles(self):
        super(TriggerObjectsAnalyzer, self).declareHandles()
        self.handles['TriggerBits'] = AutoHandle( self.triggerBitsInputTag, 'edm::TriggerResults' )
        self.handles['TriggerObjects'] = AutoHandle( self.triggerObjectInputTag, 'std::vector<pat::TriggerObjectStandAlone>' )

    def beginLoop(self, setup):
        super(TriggerObjectsAnalyzer,self).beginLoop(setup)

    def process(self, event):
        self.readCollections( event.input )
        # get the trigger names (only in the first event)
        if self.names is None:
            triggerBits = self.handles['TriggerBits'].product()
            self.names = event.input.object().triggerNames(triggerBits)
        # get the trigger object
        allTriggerObjects = self.handles['TriggerObjects'].product()
        # for each collection name save the trigger object matching with the triggerObjectsInfo (eg. )
        for collectionName in self.triggerObjectsInfos.keys():
            triggerObjectsInfo = self.triggerObjectsInfos[collectionName]
            objects = []
            for ob in allTriggerObjects:
                if (triggerObjectsInfo.collectionText!="::HLT") and triggerObjectsInfo.collectionText!=ob.collection(): continue
                if (triggerObjectsInfo.filterName!="") and not (triggerObjectsInfo.filterName in ob.filterLabels()): continue
                if (triggerObjectsInfo.path!=""):
                    ob.unpackPathNames(self.names)
                    if not FindTrigger(triggerObjectsInfo.path, ob.pathNames()): continue
                objects.append(ob)
            setattr(event,'trgObjects_'+collectionName,objects)

setattr(TriggerObjectsAnalyzer,"defaultConfig",cfg.Analyzer(
    TriggerObjectsAnalyzer, name="TriggerObjectsAnalyzerDefault",
    triggerObjectsCfgs = {"caloJets":(("hltAK4CaloJetsCorrectedIDPassed")),"caloMet":(("hltMet","","HLT"),"hltMET90","HLT_PFMET90_PFMHT90_IDTight*")},
    triggerObjectInputTag = ('selectedPatTrigger','','PAT'),
    triggerBitsInputTag = ('TriggerResults','','HLT')
)
)


