#! /usr/bin/env python


import sys
import os
import commands
import string
import time

files = [
"110/vhbb_WS_BDT_H110_WenHighPt_8TeV.root",
"110/vhbb_WS_BDT_H110_WenLowPt_8TeV.root",
"110/vhbb_WS_BDT_H110_WenMidPt_8TeV.root",
"110/vhbb_WS_BDT_H110_WmnHighPt_8TeV.root",
"110/vhbb_WS_BDT_H110_WmnLowPt_8TeV.root",
"110/vhbb_WS_BDT_H110_WmnMidPt_8TeV.root",
"110/vhbb_WS_BDT_M110_ZeeHighPt_8TeV.root",
"110/vhbb_WS_BDT_M110_ZeeLowPt_8TeV.root",
"110/vhbb_WS_BDT_M110_ZmmHighPt_8TeV.root",
"110/vhbb_WS_BDT_M110_ZmmLowPt_8TeV.root",
"110/vhbb_Znn_J12_8TeV.root",
"110/Wtn_BDT_newBinning_110.root",
"115/vhbb_WS_BDT_H115_WenHighPt_8TeV.root",
"115/vhbb_WS_BDT_H115_WenLowPt_8TeV.root",
"115/vhbb_WS_BDT_H115_WenMidPt_8TeV.root",
"115/vhbb_WS_BDT_H115_WmnHighPt_8TeV.root",
"115/vhbb_WS_BDT_H115_WmnLowPt_8TeV.root",
"115/vhbb_WS_BDT_H115_WmnMidPt_8TeV.root",
"115/vhbb_WS_BDT_M115_ZeeHighPt_8TeV.root",
"115/vhbb_WS_BDT_M115_ZeeLowPt_8TeV.root",
"115/vhbb_WS_BDT_M115_ZmmHighPt_8TeV.root",
"115/vhbb_WS_BDT_M115_ZmmLowPt_8TeV.root",
"115/vhbb_Znn_J12_8TeV.root",
"115/Wtn_BDT_newBinning_115.root",
"120/vhbb_WS_BDT_H120_WenHighPt_8TeV.root",
"120/vhbb_WS_BDT_H120_WenLowPt_8TeV.root",
"120/vhbb_WS_BDT_H120_WenMidPt_8TeV.root",
"120/vhbb_WS_BDT_H120_WmnHighPt_8TeV.root",
"120/vhbb_WS_BDT_H120_WmnLowPt_8TeV.root",
"120/vhbb_WS_BDT_H120_WmnMidPt_8TeV.root",
"120/vhbb_WS_BDT_M120_ZeeHighPt_8TeV.root",
"120/vhbb_WS_BDT_M120_ZeeLowPt_8TeV.root",
"120/vhbb_WS_BDT_M120_ZmmHighPt_8TeV.root",
"120/vhbb_WS_BDT_M120_ZmmLowPt_8TeV.root",
"120/vhbb_Znn_J12_8TeV.root",
"120/Wtn_BDT_newBinning_120.root",
"125/vhbb_WS_BDT_H125_WenHighPt_8TeV.root",
"125/vhbb_WS_BDT_H125_WenLowPt_8TeV.root",
"125/vhbb_WS_BDT_H125_WenMidPt_8TeV.root",
"125/vhbb_WS_BDT_H125_WmnHighPt_8TeV.root",
"125/vhbb_WS_BDT_H125_WmnLowPt_8TeV.root",
"125/vhbb_WS_BDT_H125_WmnMidPt_8TeV.root",
"125/vhbb_WS_BDT_M125_ZeeHighPt_8TeV.root",
"125/vhbb_WS_BDT_M125_ZeeLowPt_8TeV.root",
"125/vhbb_WS_BDT_M125_ZmmHighPt_8TeV.root",
"125/vhbb_WS_BDT_M125_ZmmLowPt_8TeV.root",
"125/vhbb_Znn_J12_8TeV.root",
"125/Wtn_BDT_newBinning_125.root",
"130/vhbb_WS_BDT_H130_WenHighPt_8TeV.root",
"130/vhbb_WS_BDT_H130_WenLowPt_8TeV.root",
"130/vhbb_WS_BDT_H130_WenMidPt_8TeV.root",
"130/vhbb_WS_BDT_H130_WmnHighPt_8TeV.root",
"130/vhbb_WS_BDT_H130_WmnLowPt_8TeV.root",
"130/vhbb_WS_BDT_H130_WmnMidPt_8TeV.root",
"130/vhbb_WS_BDT_M130_ZeeHighPt_8TeV.root",
"130/vhbb_WS_BDT_M130_ZeeLowPt_8TeV.root",
"130/vhbb_WS_BDT_M130_ZmmHighPt_8TeV.root",
"130/vhbb_WS_BDT_M130_ZmmLowPt_8TeV.root",
"130/vhbb_Znn_J12_8TeV.root",
"130/Wtn_BDT_newBinning_130.root",
"135/vhbb_WS_BDT_H135_WenHighPt_8TeV.root",
"135/vhbb_WS_BDT_H135_WenLowPt_8TeV.root",
"135/vhbb_WS_BDT_H135_WenMidPt_8TeV.root",
"135/vhbb_WS_BDT_H135_WmnHighPt_8TeV.root",
"135/vhbb_WS_BDT_H135_WmnLowPt_8TeV.root",
"135/vhbb_WS_BDT_H135_WmnMidPt_8TeV.root",
"135/vhbb_WS_BDT_M135_ZeeHighPt_8TeV.root",
"135/vhbb_WS_BDT_M135_ZeeLowPt_8TeV.root",
"135/vhbb_WS_BDT_M135_ZmmHighPt_8TeV.root",
"135/vhbb_WS_BDT_M135_ZmmLowPt_8TeV.root",
"135/vhbb_Znn_J12_8TeV.root",
"135/Wtn_BDT_newBinning_135.root",
"140/vhbb_WS_BDT_H140_WenHighPt_8TeV.root",
"140/vhbb_WS_BDT_H140_WenLowPt_8TeV.root",
"140/vhbb_WS_BDT_H140_WenMidPt_8TeV.root",
"140/vhbb_WS_BDT_H140_WmnHighPt_8TeV.root",
"140/vhbb_WS_BDT_H140_WmnLowPt_8TeV.root",
"140/vhbb_WS_BDT_H140_WmnMidPt_8TeV.root",
"140/vhbb_WS_BDT_M140_ZeeHighPt_8TeV.root",
"140/vhbb_WS_BDT_M140_ZeeLowPt_8TeV.root",
"140/vhbb_WS_BDT_M140_ZmmHighPt_8TeV.root",
"140/vhbb_WS_BDT_M140_ZmmLowPt_8TeV.root",
"140/vhbb_Znn_J12_8TeV.root",
"140/Wtn_BDT_newBinning_140.root",
"145/vhbb_WS_BDT_H145_WenHighPt_8TeV.root",
"145/vhbb_WS_BDT_H145_WenLowPt_8TeV.root",
"145/vhbb_WS_BDT_H145_WenMidPt_8TeV.root",
"145/vhbb_WS_BDT_H145_WmnHighPt_8TeV.root",
"145/vhbb_WS_BDT_H145_WmnLowPt_8TeV.root",
"145/vhbb_WS_BDT_H145_WmnMidPt_8TeV.root",
"145/vhbb_WS_BDT_M145_ZeeHighPt_8TeV.root",
"145/vhbb_WS_BDT_M145_ZeeLowPt_8TeV.root",
"145/vhbb_WS_BDT_M145_ZmmHighPt_8TeV.root",
"145/vhbb_WS_BDT_M145_ZmmLowPt_8TeV.root",
"145/vhbb_Znn_J12_8TeV.root",
"145/Wtn_BDT_newBinning_145.root",
"150/vhbb_WS_BDT_H150_WenHighPt_8TeV.root",
"150/vhbb_WS_BDT_H150_WenLowPt_8TeV.root",
"150/vhbb_WS_BDT_H150_WenMidPt_8TeV.root",
"150/vhbb_WS_BDT_H150_WmnHighPt_8TeV.root",
"150/vhbb_WS_BDT_H150_WmnLowPt_8TeV.root",
"150/vhbb_WS_BDT_H150_WmnMidPt_8TeV.root",
"150/vhbb_WS_BDT_M150_ZeeHighPt_8TeV.root",
"150/vhbb_WS_BDT_M150_ZeeLowPt_8TeV.root",
"150/vhbb_WS_BDT_M150_ZmmHighPt_8TeV.root",
"150/vhbb_WS_BDT_M150_ZmmLowPt_8TeV.root",
"150/vhbb_Znn_J12_8TeV.root",
"150/Wtn_BDT_newBinning_150.root"
]



###intialize
os.system("more produceSystHeader.C | grep \"FileName =\" > fn1")
os.system("more IntermediateMassMaker.C | grep \"IFILE = \" > fn2")
fn1 = open("fn1", "r")
fn2 = open("fn2","r")
sn1 = fn1.read().rstrip('\n')
sn2 = fn2.read().rstrip('\n')
os.system("sed -i 's@%s@TString FileName = \"XXX\";@g' produceSystHeader.C" % sn1)
os.system("sed -i 's@%s@  IFILE = \"XXX\";@g' IntermediateMassMaker.C" % sn2)





mass=["110","110.5","111","111.5","112","112.5","113","113.5","114","114.5", "115","115.5","116","116.5","117","117.5","118","118.5","119","119.5", "120","120.5","121","121.5","122","122.5","123","123.5","124","124.5", "125","125.5","126","126.5","127","127.5","128","128.5","129","129.5", "130","130.5","131","131.5","132","132.5","133","133.5","134","134.5", "135","135.5","136","136.5","137","137.5","138","138.5","139","139.5", "140","140.5","141","141.5","142","142.5","143","143.5","144","144.5", "145","145.5","146","146.5","147","147.5","148","148.5","149","149.5", "150", "124.6","124.7","124.8","124.9","125.1","125.2","125.3","125.4", "125.6","125.7","125.8","125.9","126.1","126.2","126.3","126.4"]

for item in mass: 
  os.system("mkdir -p newcards/%s" % item)




i = 0
BDT = "BDT"
Type = -1
#Type 0 = Zll, 1 = Wln, 2 = Znn, 3 = Wtn
max = 1
Bin = 0

for file in files:
  Bin = 0
  if "Zee" in file:
    BDT = "CMS_vhbb_BDT_Zll_8TeV"
    Type = 0
  if "Zmm" in file:
    BDT = "CMS_vhbb_BDT_Zll_8TeV"
    Type = 0
  if "We" in file: 
    BDT = "CMS_vhbb_BDT_Wln_8TeV"  
    Type = 1
  if "Wm" in file: 
    BDT = "CMS_vhbb_BDT_Wln_8TeV"  
    Type = 1
  if "Zn" in file: 
    Type = 2
    max = 3
  if "Wt" in file: 
    BDT = "BDT"  
    Type = 3
  
  while Bin < max:
    if Type == 2:
      if Bin == 0:
        BDT = "CMS_vhbb_BDT_ZnunuLowPt_8TeV"  
      if Bin == 1:
        BDT = "CMS_vhbb_BDT_ZnunuMedPt_8TeV"  
      if Bin == 2:
        BDT = "CMS_vhbb_BDT_ZnunuHighPt_8TeV"  
    os.system("rm -rf SYSTS.h count")
    os.system("sed -i 's@XXX@%s@g' produceSystHeader.C IntermediateMassMaker.C" % file)
    time.sleep(1)
    os.system("root -b -l -q produceSystHeader.C\(%i\) | grep RooDataHist | wc -l > count" % (Bin))
    time.sleep(1)
    fCount = open("count", "r")
    sCount = fCount.read().rstrip('\n')
    time.sleep(1)
    os.system("echo 'const int NS = %s; ' > SYSTS.h" % sCount)
    time.sleep(1)
    os.system("echo 'TString systs[] = {' >> SYSTS.h")
    time.sleep(1)
    os.system("root -b -l -q produceSystHeader.C\(%i\) | grep RooDataHist | sed 's/RooDataHist:://g' | sed 's/(%s)//g' | sed 's/^/\"/g' | sed 's/$/\",/g' | sed '$s/,/};/' >> SYSTS.h" % (Bin, BDT))
    time.sleep(1)
    os.system("root -b -l -q IntermediateMassMaker.C++\(%i\)" % Bin)
    time.sleep(1)
    os.system("sed -i 's@%s@XXX@g' produceSystHeader.C IntermediateMassMaker.C" % file)
    Bin = Bin + 1

