from ROOT import *
from array import array

#lumiA = 0.709
#lumiB = 2.844
#lumiA = 0.7
#lumiB = 2.8
#lumiT = 1.5
lumiA = 0.001*(809.379 + 82.136)
lumiB = 4.404 # 4.398 for ele
lumiC = 0.495+4.655+1.656
verboseOutput = True

def makeAvg(inputFiles,inputWeights,outputFile,verbose=False):
  if verbose:
    print "Input files: "," ".join(inputFiles)
    print "Output file: %s" % outputFile
    print "Input weights: "," ".join(["%.4f" % w for w in inputWeights])
  if len(inputFiles) != len(inputWeights):
    raise Exception,"Need a weight for each file"
  wtot = sum(inputWeights)
  if verbose: print "Total weight: %.4f" % wtot
  tfs = []
  tts = []
  for file in inputFiles:
    tfs.append(TFile(file))
    tfs[-1].cd()
    tts.append(tfs[-1].Get("tree"))
  nentries = tts[0].GetEntries()
  for tree in tts:
    if tree.GetEntries() != nentries:
      raise Exception,"Trees have different numbers of entries"
  tf = TFile(outputFile,"RECREATE")
  tt = TTree("tree","tree")
  ptMin = array('f',[0])
  branch_ptMin = tt.Branch("ptMin",ptMin,"ptMin/F")
  ptMax = array('f',[0])
  branch_ptMax = tt.Branch("ptMax",ptMax,"ptMax/F")
  etaMin = array('f',[0])
  branch_etaMin = tt.Branch("etaMin",etaMin,"etaMin/F")
  etaMax = array('f',[0])
  branch_etaMax = tt.Branch("etaMax",etaMax,"etaMax/F")
  scale = array('f',[0])
  branch_scale = tt.Branch("scale",scale,"scale/F")
  error = array('f',[0])
  branch_error = tt.Branch("error",error,"error/F")
  for i in range(nentries):
    for tree in tts:
      tree.GetEntry(i)
    ptMin[0] = tts[-1].ptMin
    ptMax[0] = tts[-1].ptMax
    etaMin[0] = tts[-1].etaMin
    etaMax[0] = tts[-1].etaMax
    scale[0] = 0.
    error[0] = 0.
    if verbose: print "Bin %.3f %.3f %.3f %.3f:" % (ptMin[0],ptMax[0],etaMin[0],etaMax[0])
    for i in range(len(tts)):
      if tts[i].ptMin != ptMin[0] or tts[i].ptMax != ptMax[0] or tts[i].etaMin != etaMin[0] or tts[i].etaMax != etaMax[0]:
        raise Exception,"Trees have different binning!"
      scale[0] += (inputWeights[i]/wtot)*tts[i].scale
      if tts[i].scale > 0.:
        error[0] += ((inputWeights[i]/wtot)*(tts[i].error/tts[i].scale))**2
      if verbose: print "    %i: %.3f %.3f" % (i,tts[i].scale,tts[i].error)
    error[0] = scale[0]*(error[0])**0.5
    if verbose: print "  Total: %.3f %.3f" % (scale[0],error[0])
    tt.Fill()
  tf.cd()
  tt.Write()
  tf.Close()
  if verbose: print

# example
# makeAvg(["DoubleMuDz.TrigEff.2012A.root","DoubleMuDz.TrigEff.2012B.root"],[0.5,1.7],"DoubleMuDz.TrigEff.2012AB.root",True)

# do everything in directory
from os import listdir

#for fileA in listdir("."):
#  if (fileA.count("TrigEff") or fileA.count("MuRecoId") or fileA.count("EleRecoId")) and fileA.count("2012A.root") and not fileA.count("2012AB"):
#    fileB = fileA.replace("2012A","2012B")
#    if not listdir(".").count(fileB):
#      raise Exception,"%s exists but %s does not" % (fileA,fileB)
#    makeAvg([fileA,fileB],[lumiA,lumiB],fileA.replace("2012A","2012AB"),verboseOutput)

for fileA in listdir("."):
  if (fileA.count("TrigEff") or fileA.count("MuRecoId") or fileA.count("EleRecoId")) and (fileA.count("2012A.root") or fileA.count("2012A.PUAB.root")) and not fileA.count("2012AB"):
    if fileA.count("PUAB"):
      fileB = fileA.replace("2012A","2012B")
      if not listdir(".").count(fileB):
        raise Exception,"%s exists but %s does not" % (fileA,fileB)
      if fileA.count("Ele"):
        makeAvg([fileA,fileB],[lumiA,lumiB-0.006],fileA.replace("2012A","2012AB").replace(".PUAB",""),verboseOutput)
      else:
        makeAvg([fileA,fileB],[lumiA,lumiB-0.006],fileA.replace("2012A","2012AB").replace(".PUAB",""),verboseOutput)
    else:  
      fileB = fileA.replace("2012A","2012B")
      fileC = fileA.replace("2012A","2012C")
      if not listdir(".").count(fileB):
        raise Exception,"%s exists but %s does not" % (fileA,fileB)
      if not listdir(".").count(fileC):
        raise Exception,"%s exists but %s does not" % (fileA,fileC)
      if fileA.count("Ele"):
        makeAvg([fileA,fileB,fileC],[lumiA,lumiB-0.006,lumiC],fileA.replace("2012A","2012ABC"),verboseOutput)
      else:
        makeAvg([fileA,fileB,fileC],[lumiA,lumiB-0.006,lumiC],fileA.replace("2012A","2012ABC"),verboseOutput)
      if fileA.count("TrigEff"):
        fileB = fileA.replace("2012A","2012B")
        if not listdir(".").count(fileB):
          raise Exception,"%s exists but %s does not" % (fileA,fileB)
        if fileA.count("Ele"):
          makeAvg([fileA,fileB],[lumiA,lumiB-0.006],fileA.replace("2012A","2012AB"),verboseOutput)
        else:
          makeAvg([fileA,fileB],[lumiA,lumiB-0.006],fileA.replace("2012A","2012AB"),verboseOutput)
                            
