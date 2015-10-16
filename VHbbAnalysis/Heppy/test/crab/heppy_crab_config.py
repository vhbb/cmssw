from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'vhbb_heppy_blike_final_v13_v2b_id_sig'
config.General.workArea = 'crab_projects_blike_final_v13_v2b_id_sig'
config.General.transferLogs=True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'heppy_crab_fake_pset.py'
config.JobType.scriptExe = 'heppy_crab_script.sh'
import os
os.system("tar czf python.tar.gz --dereference --directory $CMSSW_BASE python")
config.JobType.inputFiles = ['heppy_config.py',
                             'heppy_crab_script.py',
                             'python.tar.gz',
                             'MVAJetTags_620SLHCX_Phase1And2Upgrade.db',
                             'combined_cmssw.py',
                             '../vhbb.py',
                              '../vhbb_combined.py',
                             'TMVAClassification_BDT.weights.xml',
                             '../pdfQG_AK4chs_13TeV_v2b.root',
                             'puData.root',
                             'puMC.root',
                              '../json.txt',
                              "../Zll-spring15.weights.xml",
                              "../Wln-spring15.weights.xml",
                              "../Znn-spring15.weights.xml",
                              "../VBF-spring15.weights.xml",
										'../TMVA_blikelihood_vbf_singlebtag_v13_id.xml'
]
#config.JobType.outputFiles = ['tree.root']

config.section_("Data")
config.Data.inputDataset = '/VBFHToBB_M-125_13TeV_powheg_pythia8_weightfix/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset = '/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits = 6
config.Data.outLFNDirBase = '/store/user/nchernya/vhbb_final_v13_0_v2b_id_sig/'
config.Data.publication = False
config.Data.publishDataName = 'VHBB_heppy_blike_v2bsig'

config.section_("Site")
config.Site.storageSite = "T2_IT_Rome"

#config.Data.ignoreLocality = True
