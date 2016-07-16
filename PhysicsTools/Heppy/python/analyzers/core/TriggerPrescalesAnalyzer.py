import ROOT

from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.analyzers.core.AutoFillTreeProducer import NTupleVariable
import PhysicsTools.HeppyCore.framework.config as cfg

""" triggerNameDict: return a dictionary with the trigger number for each trigger name"""
def updateTriggerNumberDict(triggerNameDict,names):
    for i in range(names.size()):
        triggerName = names.triggerName(i)
        for name in triggerNameDict:
            if name in triggerName:
                triggerNameDict[name] = i
    
    return triggerNameDict

"""
TriggerPrescalesAnalyzer is a class that ..........
"""
class TriggerPrescalesAnalyzer( Analyzer ):
    def __init__(self, cfg_ana, cfg_comp, looperName ):
        super(TriggerPrescalesAnalyzer,self).__init__(cfg_ana,cfg_comp,looperName)
        self.processName                = getattr(self.cfg_ana,"processName","PAT")
        self.triggerBitsInputTag        = getattr(self.cfg_ana,"triggerBitsInputTag",('TriggerResults','','HLT'))
        self.l1Prescales                = getattr(self.cfg_ana,"l1Prescales",False)
        self.triggerNameDict            = {}
        for triggerName in getattr(self.cfg_ana,"triggerList",[]):
            self.triggerNameDict[triggerName] = -1
        self.triggerPrescalesInfos      = {}
        self.runNumber                  = -1
        self.names                      = None
        
    def declareHandles(self):
        super(TriggerPrescalesAnalyzer, self).declareHandles()
        self.handles['triggerBits']         = AutoHandle( self.triggerBitsInputTag, 'edm::TriggerResults' )
        self.handles['triggerPrescales']    = AutoHandle( ('patTrigger','',self.processName), 'pat::PackedTriggerPrescales' )
        if self.l1Prescales:
            self.handles['triggerPrescalesL1Min']    = AutoHandle( ('patTrigger','l1min',self.processName), 'pat::PackedTriggerPrescales' )
            self.handles['triggerPrescalesL1Max']    = AutoHandle( ('patTrigger','l1max',self.processName), 'pat::PackedTriggerPrescales' )

    def beginLoop(self, setup):
        super(TriggerPrescalesAnalyzer,self).beginLoop(setup)
        for triggerName in self.triggerNameDict:
            setup.globalVariables.append( NTupleVariable('prescalesTest_%s'%triggerName, eval("lambda ev : ev.prescalesTest_%s"%triggerName), int, mcOnly=False, help="trigger prescale of %s"%triggerName ))
            if self.l1Prescales:
                setup.globalVariables.append( NTupleVariable('HLTprescaleL1Min_%s'%triggerName, eval("lambda ev : ev.HLTprescaleL1Min_%s"%triggerName), int, mcOnly=False,    help="Min L1 seed prescale for %s"%triggerName ))
                setup.globalVariables.append( NTupleVariable('HLTprescaleL1Max_%s'%triggerName, eval("lambda ev : ev.HLTprescaleL1Max_%s"%triggerName), int, mcOnly=False,    help="Max L1 seed prescale for %s"%triggerName ))

    def process(self, event):
        self.readCollections( event.input )
        run = event.input.eventAuxiliary().id().run()
        
        # get the trigger names (only for the first event of each run)
        if self.runNumber!= run:
            triggerBits             = self.handles['triggerBits'].product()
            self.names              = event.input.object().triggerNames(triggerBits)
            self.runNumber          = run
            self.triggerNameDict    = updateTriggerNumberDict(self.triggerNameDict,self.names)
        
        # get the trigger prescales
        triggerPrescales = self.handles['triggerPrescales'].product()
        if self.l1Prescales:
            triggerPrescalesL1Min = self.handles['triggerPrescalesL1Min'].product()
            triggerPrescalesL1Max = self.handles['triggerPrescalesL1Max'].product()
        
        for triggerName,triggerNumber in self.triggerNameDict.items():
            prescale = -1
            if triggerNumber>=0:
                prescale = triggerPrescales.getPrescaleForIndex(triggerNumber)
            setattr(event,'prescalesTest_'+triggerName,prescale)
            if self.l1Prescales:
                prescaleL1Max = -1
                prescaleL1Min = -1
                if triggerNumber>=0:
                    prescaleL1Max = triggerPrescalesL1Min.getPrescaleForIndex(triggerNumber)
                    prescaleL1Min = triggerPrescalesL1Max.getPrescaleForIndex(triggerNumber)
                setattr(event,'HLTprescaleL1Min_'+triggerName,prescaleL1Max)
                setattr(event,'HLTprescaleL1Max_'+triggerName,prescaleL1Min)

setattr(TriggerPrescalesAnalyzer,"defaultConfig",cfg.Analyzer(
    TriggerPrescalesAnalyzer,
    name="TriggerPrescalesAnalyzerDefault",
    triggerList = ["HLT_PFMET90_PFMHT90_IDTight","HLT_PFHT"],
    processName = 'PAT',
    triggerBitsInputTag = ('TriggerResults','','HLT'),
    l1Prescales = True
)
)


