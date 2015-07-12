import sys,os

# REMOVE DATASET NAMES CONTAINING:

remove = [
          'http','prime','GravToZZ','GravitonToGluonGluon','GravitonToQuarkQuark','GravToGG','GravToWW','ToHHTo2B',
          'SUSY','QstarTo','RSGluonTo','WRTo','TstarTstar','Unpart','LQTo','BstarTo','WpWpJJ','WZTo3LNu',
          'HToZZ','HToWW','HToG','HToT','/ADD','/GJet','GluGluToZZ','TTbarDM','HToInvisible','WToENu_M','WToMuNu_M','WToTauNu_M',
          'ttHJetToGG','ttHJetToTT','Muminus_Pt','/Muplus','Photon','SinglePion','ZZTo4L','DoubleElectron',
          'SingleEta','tGamma','JPsiToMuMu','JpsiToMuMu','mtop1','BdToJpsiKs','tZq_','GG_M','HToNonbb',
          'DYJetsToLL_M-1000to1500','DYJetsToLL_M-100to200','DYJetsToLL_M-10to50','DYJetsToLL_M-1500to2000',
          'DYJetsToLL_M-2000to3000','DYJetsToLL_M-400to500','DYJetsToLL_M-500to700','DYJetsToLL_M-500to700',
          'DYJetsToLL_M-200to400','DYJetsToLL_M-700to800','DYJetsToLL_M-800to1000','BuToJpsiK','GluGluHToZG',
          'GGJets','Monotop_S','TTJets_Mtt-','TT_Mtt-',
          '/ttHJetToNonbb_M120_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM',
          '/ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM',
          '/ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext1-v3/MINIAODSIM',
          '/ttHJetToNonbb_M130_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM',

          ]

# FILELIST OF AVAILABLE DATASETS ON DAS AS VALID

# das_valid = [line.rstrip('\n').rstrip('\r') for line in open('all_datasets_MCRUN2_25ns.txt')]
das_valid = os.popen('python ./das_client.py --limit=0 --query="dataset=/*/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9*/MINIAODSIM"').read().split('\n')
das_valid = filter(None, das_valid)

# FILELIST OF AVAILABLE DATASETS ON DAS AS PRODUCTION

das_production = os.popen('python ./das_client.py --limit=0 --query="dataset dataset=/*/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9*/MINIAODSIM status=PRODUCTION"').read().split('\n')
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
    
vhbb_prod = []
print '\nDATASETS AVAILABLE ON DAS AS PRODUCTION AND NOT (YET) INCLUDED IN THE HBB LIST\n'
for line in das_production:
  if line not in vhbb_all:
    print line
  if line in vhbb_all:
    vhbb_prod.append(line)

print '\nDATASETS INCLUDED IN THE HBB LIST AND STILL IN PRODUCTION\n'
for line in vhbb_prod:
    print line

print '\nDATASETS INCLUDED IN THE HBB LIST NOT IN PRODUCTION NOR IN VALID state\n'
for line in vhbb_prod:
  if (line not in das_production) and (line not in das_valid):
    print line

print '\n'