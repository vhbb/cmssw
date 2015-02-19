from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'VHBB_HEPPY_V8_003'
config.General.workArea = 'crab_projects_V8_003'
config.General.transferLogs=True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'heppy_crab_fake_pset.py'
config.JobType.scriptExe = 'heppy_crab_script.sh'
import os
os.system("tar czf python.tar.gz --dereference --directory $CMSSW_BASE python")
config.JobType.inputFiles = ['heppy_config.py','heppy_crab_script.py','python.tar.gz','MVAJetTags_620SLHCX_Phase1And2Upgrade.db','newbtag.py']
#config.JobType.outputFiles = ['tree.root']

config.section_("Data")
config.Data.inputDataset = '/WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp/Phys14DR-PU40bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
config.Data.outLFN = '/store/user/arizzi/VHBBHeppyV8/'
config.Data.publication = True
config.Data.publishDataName = 'VHBB_HEPPY_V8'

config.section_("Site")
config.Site.storageSite = "T2_IT_Pisa"

#config.Data.ignoreLocality = True
