from WMCore.Configuration import Configuration
config = Configuration()

version="V25b"
submission="_004"

config.section_("General")
config.General.requestName = 'VHBB_HEPPY_'+version+submission
config.General.workArea = '/scratch/arizzi/crab_sub/crab_projects_'+version+submission
config.General.transferLogs=True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'heppy_crab_fake_pset.py'
config.JobType.scriptExe = 'heppy_crab_script.sh'
import os
os.system("tar czf python.tar.gz --directory $CMSSW_BASE python `find $CMSSW_BASE/src -name python | perl -pe s#$CMSSW_BASE/## `")

#os.system("tar czf python.tar.gz --dereference --directory $CMSSW_BASE python")
config.JobType.maxMemoryMB = 2400
config.JobType.maxJobRuntimeMin = 2000

config.JobType.inputFiles = ['heppy_config.py',
                             'heppy_crab_script.py',
                             'python.tar.gz',
#                             'MVAJetTags_620SLHCX_Phase1And2Upgrade.db',
                             'combined_cmssw.py',
                             '../vhbb.py',
                             '../vhbb_combined.py',
                             '../TMVA_blikelihood_vbf_cmssw76.weights.xml',
                             'TMVAClassification_BDT.weights.xml',
                             'puData.root',
                             'puMC.root',
                             'puDataMinus.root',
                             'puDataPlus.root',
                              'json.txt',
                              '../silver.txt',
                              #"../Zll-spring15.weights.xml",
                              #"../Wln-spring15.weights.xml",
                              #"../Znn-spring15.weights.xml",
                              #"../VBF-spring15.weights.xml",
                              #"../ttbar-spring15.weights.xml",
                              #"../ttbar-fall15.weights.xml",
                              #"../ttbar-fall15_TargetGenOverPt_GenPtCut0.weights.xml",
                              '../ttbar-G25-500k-13d-300t.weights.xml',
                              '../triggerEmulation.root',
			      #'../ttbar-spring16-80X.weights.xml',
                              '../TMVA_blikelihood_vbf_cmssw76_h21trained.weights.xml',
]
#config.JobType.outputFiles = ['tree.root']

config.section_("Data")
config.Data.inputDataset = '/ZH_HToBB_ZToLL_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 50000 #30
#onfig.Data.totalUnits = 2
config.Data.outLFNDirBase = '/store/user/arizzi/VHBBHeppy%s/'%version
config.Data.publication = True
config.Data.outputDatasetTag = 'VHBB_HEPPY_%s'%version
#only for data
#config.Data.lumiMask = 'json.txt'

config.section_("Site")
config.Site.storageSite = "T2_IT_Pisa"
config.Site.blacklist = ['T2_FR_CCIN2P3','T1_FR_CCIN2P3','T2_BR_SPRACE','T2_ES_CIEMAT']

#config.Data.ignoreLocality = True
