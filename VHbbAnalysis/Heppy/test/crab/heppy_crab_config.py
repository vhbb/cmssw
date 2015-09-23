from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'vhbb_heppy_blike'
config.General.workArea = 'crab_projects_blike'
config.General.transferLogs=True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'heppy_crab_fake_pset.py'
#config.JobType.scriptExe = 'heppy_crab_script.sh'
import os
os.system("tar czf python.tar.gz --dereference --directory $CMSSW_BASE python")
config.JobType.inputFiles = ['heppy_config.py',
                             'heppy_crab_script.py',
                             'python.tar.gz',
                             'MVAJetTags_620SLHCX_Phase1And2Upgrade.db',
                             'combined_cmssw.py',
                             '../vhbb.py',
                             'TMVAClassification_BDT.weights.xml',
                             'pdfQG_AK4chs_antib_13TeV_v1.root',
#                             '../jec/PHYS14_V4_MC_L1FastJet_AK4PFchs.txt',  
#                             '../jec/PHYS14_V4_MC_L2Relative_AK4PFchs.txt',  
#                             '../jec/PHYS14_V4_MC_L3Absolute_AK4PFchs.txt',
#                             '../jec/Uncertainty_FAKE.txt',
#                             '../csv/csv_rwt_hf_IT_FlatSF_2015_07_27.root',
#                             '../csv/csv_rwt_lf_IT_FlatSF_2015_07_27.root',
                             'Wln_weights_phys14.xml',
                             'Zll_weights_phys14.xml',
                             'Znn_weights_phys14.xml',
										'TMVA_blikelihood_vbf_singlebtag.xml'
]
#config.JobType.outputFiles = ['tree.root']

config.section_("Data")
config.Data.inputDataset = '/VBFHToBB_M-125_13TeV_powheg_pythia8_weightfix/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
#config.Data.totalUnits = 1
config.Data.outLFNDirBase = '/store/user/nchernya/vhbb_blike/'
config.Data.publication = False
config.Data.publishDataName = 'VHBB_heppy_blike'

config.section_("Site")
config.Site.storageSite = "T2_AT_Vienna"

#config.Data.ignoreLocality = True
