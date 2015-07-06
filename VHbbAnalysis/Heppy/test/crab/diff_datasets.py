import sys,os

# REMOVE DATASET NAMES CONTAINING:

remove = ['http','prime','Grav','SUSY','QstarTo','RSGluonTo','WRTo','TstarTstar','Unpart','LQTo','Radion','BstarTo',
          'HToZZ','HToWW','HToG','HToT','/ADD','GJet','GluGluToZZ','TTbarDM','HToInvisible','WToENu_M','WToMuNu_M','WToTauNu_M',
          'ttHJetToGG','ttHJetToNonbb','ttHJetToTT','Muminus_Pt','/Muplus','Photon','SinglePion','ZZTo4L','DoubleElectron',
          'SingleEta','tGamma','JPsiToMuMu','JpsiToMuMu','mtop1','BdToJpsiKs','tZq_','GG_M',
          'DYJetsToLL_M-1000to1500','DYJetsToLL_M-100to200','DYJetsToLL_M-10to50','DYJetsToLL_M-1500to2000',
          'DYJetsToLL_M-2000to3000','DYJetsToLL_M-400to500','DYJetsToLL_M-500to700','DYJetsToLL_M-500to700',
          'DYJetsToLL_M-200to400','DYJetsToLL_M-700to800','BuToJpsiK','GluGluHToZG']


# FILELIST OF AVAILABLE DATASETS ON DAS AS VALID

# das_valid = [line.rstrip('\n').rstrip('\r') for line in open('all_datasets_MCRUN2_25ns.txt')]
das_valid = os.popen('python ./das_client.py --limit=0 --query="dataset=/*/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v*/MINIAODSIM"').read().split('\n')
das_valid = filter(None, das_valid)

# FILELIST OF AVAILABLE DATASETS ON DAS AS PRODUCTION

das_production = os.popen('python ./das_client.py --limit=0 --query="dataset dataset=/*/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v*/MINIAODSIM status=PRODUCTION"').read().split('\n')
das_production = filter(None, das_production)

# FILTER DATASETS WITH REMOVE LIST
print '\nretrieved DAS lists contain',len(das_valid),'datasets in VALID state,',len(das_production),'datasets in PRODUCTION state'

for i in remove:
  das_valid = [ x for x in das_valid if i not in x ]
  das_production = [ x for x in das_production if i not in x ]

print '\nWARNING: filtering datasets containing',remove
print '\nfiltered lists contain',len(das_valid),'datasets in VALID state,',len(das_production),'datasets in PRODUCTION state'

# CHECK EXISTING LIST OF DATASETS TO BE PROCESSED

filename = 'datasets_MCRUN2_25ns.txt'
vhbb_all = open(filename).read()

print 'HBB filelist ',filename,'contains',len(filter(None, vhbb_all.split('\n'))),'datasets'

print '\nDATASETS AVAILABLE ON DAS AS VALID AND NOT (YET) INCLUDED IN THE HBB LIST\n'
for line in das_valid:
  if line not in vhbb_all:
    print line
    
print '\nDATASETS AVAILABLE ON DAS AS PRODUCTION AND NOT (YET) INCLUDED IN THE HBB LIST\n'
for line in das_production:
  if line not in vhbb_all:
    print line

print '\n'