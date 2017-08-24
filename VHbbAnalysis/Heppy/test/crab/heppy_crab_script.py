#!/usr/bin/env python
import os
import PhysicsTools.HeppyCore.framework.config as cfg
cfg.Analyzer.nosubdir=True

import PSet
import sys
import re
print "ARGV:",sys.argv
JobNumber=sys.argv[1]
crabFiles=PSet.process.source.fileNames
print crabFiles
firstInput = crabFiles[0]
tested=False
forceaaa=False
print "--------------- using edmFileUtil to convert PFN to LFN -------------------------"
for i in xrange(0,len(crabFiles)) :
     if os.getenv("GLIDECLIENT_Group","") != "overflow" and  os.getenv("GLIDECLIENT_Group","") != "overflow_conservative" and not forceaaa:
       print "Data is local"
       pfn=os.popen("edmFileUtil -d %s"%(crabFiles[i])).read() 
       pfn=re.sub("\n","",pfn)
       print crabFiles[i],"->",pfn
       if not tested:
         print "Testing file open"
         import ROOT
         testfile=ROOT.TFile.Open(pfn)
         if testfile and testfile.IsOpen() :
            print "Test OK"
       	    crabFiles[i]=pfn
            testfile.Close()
            #tested=True
         else :
            print "Test open failed, forcing AAA"
            crabFiles[i]="root://cms-xrd-global.cern.ch/"+crabFiles[i]
            forceaaa=True
       else :
	    crabFiles[i]=pfn


     else:
       print "Data is not local, using AAA/xrootd"
       crabFiles[i]="root://cms-xrd-global.cern.ch/"+crabFiles[i]

import imp
handle = open("heppy_config.py", 'r')
cfo = imp.load_source("heppy_config", "heppy_config.py", handle)
config = cfo.config
handle.close()

#replace files with crab ones
config.components[0].files=crabFiles

#adjust global tag for DATA
print "Setting Data GT if needed"
mm=re.match('.*(Run2016.).*',crabFiles[0])
gtmap={}
gtmap["Run2016B"]='Summer16_23Sep2016BCDV4_DATA'
gtmap["Run2016C"]='Summer16_23Sep2016BCDV4_DATA'
gtmap["Run2016D"]='Summer16_23Sep2016BCDV4_DATA'
gtmap["Run2016E"]='Summer16_23Sep2016EFV4_DATA'
gtmap["Run2016F"]='Summer16_23Sep2016EFV4_DATA'
gtmap["Run2016G"]='Summer16_23Sep2016GV4_DATA'
gtmap["Run2016H"]='Summer16_23Sep2016HV4_DATA'

for x in config.sequence :
  if x.name == "PhysicsTools.Heppy.analyzers.objects.JetAnalyzer.JetAnalyzer_1" :
    JetAna=x
if mm :
  JetAna.dataGT=gtmap[mm.group(1)]
  print "Updated data GT: ",   JetAna.dataGT

if hasattr(PSet.process.source, "lumisToProcess"):
    config.preprocessor.options["lumisToProcess"] = PSet.process.source.lumisToProcess

os.system("ps aux |grep heppy")
from PhysicsTools.HeppyCore.framework.looper import Looper
looper = Looper( 'Output', config, nPrint = 1)
looper.loop()
looper.write()

#print PSet.process.output.fileName
#os.system("ls -lR")
os.rename("Output/tree.root", "tree.root")
#os.system("ls -lR")

import ROOT
f=ROOT.TFile.Open('tree.root')
entries=f.Get('tree').GetEntries()

fwkreport='''<FrameworkJobReport>
<ReadBranches>
</ReadBranches>
<PerformanceReport>
  <PerformanceSummary Metric="StorageStatistics">
    <Metric Name="Parameter-untracked-bool-enabled" Value="true"/>
    <Metric Name="Parameter-untracked-bool-stats" Value="true"/>
    <Metric Name="Parameter-untracked-string-cacheHint" Value="application-only"/>
    <Metric Name="Parameter-untracked-string-readHint" Value="auto-detect"/>
    <Metric Name="ROOT-tfile-read-totalMegabytes" Value="0"/>
    <Metric Name="ROOT-tfile-write-totalMegabytes" Value="0"/>
  </PerformanceSummary>
</PerformanceReport>

<GeneratorInfo>
</GeneratorInfo>

<InputFile>
<LFN>%s</LFN>
<PFN></PFN>
<Catalog></Catalog>
<InputType>primaryFiles</InputType>
<ModuleLabel>source</ModuleLabel>
<GUID></GUID>
<InputSourceClass>PoolSource</InputSourceClass>
<EventsRead>1</EventsRead>

</InputFile>

<File>
<LFN></LFN>
<PFN>tree.root</PFN>
<Catalog></Catalog>
<ModuleLabel>HEPPY</ModuleLabel>
<GUID></GUID>
<OutputModuleClass>PoolOutputModule</OutputModuleClass>
<TotalEvents>%s</TotalEvents>
<BranchHash>dc90308e392b2fa1e0eff46acbfa24bc</BranchHash>
</File>

</FrameworkJobReport>''' % (firstInput,entries)


os.system("perl -pe 's/.*cmsswPreProcessing.root.*/\<PFN\>tree.root\<\/PFN\>/' -i ./FrameworkJobReport.xml")

#f1=open('./FrameworkJobReport.xml', 'w+')
#f1.write(fwkreport)
