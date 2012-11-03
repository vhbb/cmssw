#! /usr/bin/env python


import sys
import os
import commands
import string

xSec7ZH = [0.4721,0.4655,0.4589,0.4525,0.4462,0.4400,0.4340,0.4280,0.4221,0.4164, 0.4107, 0.4052, 0.3998, 0.3945, 0.3893, 0.3842, 0.3791, 0.3742, 0.3693, 0.3645, 0.3598, 0.3551, 0.3505, 0.3459, 0.3414, 0.3370, 0.3326, 0.3283, 0.3241, 0.3199, 0.3158, 0.3117, 0.3077, 0.3038, 0.2999, 0.2961, 0.2923, 0.2886, 0.2849, 0.2813, 0.2778, 0.2743, 0.2709, 0.2675, 0.2642, 0.2609, 0.2577, 0.2545, 0.2514, 0.2483, 0.2453, 0.31908, 0.31826, 0.31744, 0.31662, 0.314980, 0.314160, 0.313340, 0.312520, 0.310900, 0.310100, 0.309300, 0.308500, 0.306920, 0.306140, 0.305360, 0.304580]



xSec7WH  = [0.8754,0.8623,0.8495,0.8368,0.8244,	0.8122,	0.8003,	0.7885,	0.7770,	0.7657,	0.7546,	0.7439,	0.7333,	0.7230,	0.7129,	0.7030,	0.6933,	0.6837,	0.6744,	0.6651,	0.6561,	0.6472,	0.6384,	0.6297,	0.6212,	0.6129,	0.6046,	0.5965,	0.5885,	0.5806,	0.5729,	0.5652,	0.5576,	0.5501,	0.5428,	0.5355,	0.5284,	0.5213,	0.5144,	0.5075,	0.5008,	0.4942,	0.4877,	0.4813,	0.4749,	0.4687,	0.4626,	0.4566,	0.4506,	0.4448,	0.4390,0.579060,0.577520,0.575980,0.574440,0.571360,0.569820,0.568280,0.566740,0.563680,0.562160,0.560640,0.559120,0.556100,0.554600,0.553100,0.551600]

xSec8ZH = [0.5869,0.5788,0.5708,0.5629,0.5552,0.5476,0.5402,0.5329,0.5258,0.5187,0.5117,0.5049,0.4981,0.4916,0.4850,0.4787,0.4724,0.4662,0.4602,0.4542,0.4483,0.4426,0.4368,0.4312,0.4257,0.4203,0.4150,0.4096,0.4044,0.3993,0.3943,0.3893,0.3843,0.3794,0.3746,0.3699,0.3652,0.3606,0.3561,0.3516,0.3473,0.3430,0.3388,0.3347,0.3306,0.3266,0.3226,0.3188,0.3149,0.3112,0.3074,0.39830,0.39730,0.39630,0.39530,0.39330,0.39230,0.39130,0.39030,0.38830,0.38730,0.38630,0.38530,0.38332,0.38234,0.38136,0.38038]
          
xSec8WH = [1.060 ,1.045 ,1.030 ,1.015 ,	0.9998,	0.9852,	0.9709,	0.9570,	0.9432,	0.9297,	0.9165,	0.9035,	0.8907,	0.8782,	0.8659,	0.8538,	0.8420,	0.8303,	0.8187,	0.8075,	0.7966,	0.7859,	0.7753,	0.7649,	0.7547,	0.7446,	0.7347,	0.7249,	0.7154,	0.7060,	0.6966,	0.6873,	0.6782,	0.6691,	0.6602,	0.6515,	0.6429,	0.6344,	0.6260,	0.6177,	0.6095,	0.6015,	0.5936,	0.5859,	0.5783,	0.5708,	0.5634,	0.5562,	0.5491,	0.5420,	0.5351, 0.70412,0.70224,0.70036,0.69848,0.69474,0.69288,0.69102,0.68916,0.68548,0.68366,0.68184,0.68002,0.67638,0.67456,0.67274,0.67092]

mass=["110","110.5","111","111.5","112","112.5","113","113.5","114","114.5","115","115.5","116","116.5","117","117.5","118","118.5","119","119.5","120","120.5","121","121.5","122","122.5","123","123.5","124","124.5","125","125.5","126","126.5","127","127.5","128","128.5","129","129.5","130","130.5","131","131.5","132","132.5","133","133.5","134","134.5","135","124.6","124.7","124.8","124.9","125.1","125.2","125.3","125.4","125.6","125.7","125.8","125.9","126.1","126.2","126.3","126.4"];


xSec = xSec8ZH


files=[
"vhbb_DC_BDT_M110_Zll_8TeV.txt",
"vhbb_DC_BDT_M115_Zll_8TeV.txt",
"vhbb_DC_BDT_M120_Zll_8TeV.txt",
"vhbb_DC_BDT_M125_Zll_8TeV.txt",
"vhbb_DC_BDT_M130_Zll_8TeV.txt",
"vhbb_DC_BDT_M135_Zll_8TeV.txt"
]



def Process(file, toMass, fromMass):
       newcard = "./%s" % mass[toMass]
       newcard += "/vhbb_Zll_8TeV.txt"
       os.system('more %s | grep rate > a' % file)
       os.system("more a |  awk '{print $2}'  > aVH1")
       os.system("more a |  awk '{print $8}'  > aVH2")
       os.system("more a |  awk '{print $14}'  > aVH3")
       os.system("more a |  awk '{print $20}'  > aVH4")
 
       fVH1 = open("aVH1", "r")
       fVH2 = open("aVH2", "r")
       fVH3 = open("aVH3", "r")
       fVH4 = open("aVH4", "r")


       sVH1 = fVH1.read().rstrip('\n')
       sVH2 = fVH2.read().rstrip('\n')
       sVH3 = fVH3.read().rstrip('\n')
       sVH4 = fVH4.read().rstrip('\n')


       VH1 = float(sVH1)
       VH2 = float(sVH2)
       VH3 = float(sVH3)
       VH4 = float(sVH4)


       bVH1 = (xSec[toMass]/xSec[fromMass])*VH1
       bVH2 = (xSec[toMass]/xSec[fromMass])*VH2
       bVH3 = (xSec[toMass]/xSec[fromMass])*VH3
       bVH4 = (xSec[toMass]/xSec[fromMass])*VH4

       os.system("sed 's/%s/%f/g' a > b" % (sVH1,bVH1))
       os.system("sed -i 's/%s/%f/g' b" % (sVH2,bVH2))
       os.system("sed -i 's/%s/%f/g' b" % (sVH3,bVH3))
       os.system("sed -i 's/%s/%f/g' b" % (sVH4,bVH4))

       fRO = open("a", "r")
       sRO = fRO.read().rstrip('\n')
       fRF = open("b", "r")
       sRF = fRF.read().rstrip('\n')
       os.system("sed 's/%s/%s/g' %s > %s" % (sRO,sRF,file,newcard))

       os.system("sed -i 's/vhbb_WS_BDT_M110_ZeeHighPt_8TeV.root/vhbb_ZeeHighPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M115_ZeeHighPt_8TeV.root/vhbb_ZeeHighPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M120_ZeeHighPt_8TeV.root/vhbb_ZeeHighPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M125_ZeeHighPt_8TeV.root/vhbb_ZeeHighPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M130_ZeeHighPt_8TeV.root/vhbb_ZeeHighPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M135_ZeeHighPt_8TeV.root/vhbb_ZeeHighPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M110_ZmmHighPt_8TeV.root/vhbb_ZmmHighPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M115_ZmmHighPt_8TeV.root/vhbb_ZmmHighPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M120_ZmmHighPt_8TeV.root/vhbb_ZmmHighPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M125_ZmmHighPt_8TeV.root/vhbb_ZmmHighPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M130_ZmmHighPt_8TeV.root/vhbb_ZmmHighPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M135_ZmmHighPt_8TeV.root/vhbb_ZmmHighPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M110_ZeeLowPt_8TeV.root/vhbb_ZeeLowPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M115_ZeeLowPt_8TeV.root/vhbb_ZeeLowPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M120_ZeeLowPt_8TeV.root/vhbb_ZeeLowPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M125_ZeeLowPt_8TeV.root/vhbb_ZeeLowPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M130_ZeeLowPt_8TeV.root/vhbb_ZeeLowPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M135_ZeeLowPt_8TeV.root/vhbb_ZeeLowPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M110_ZmmLowPt_8TeV.root/vhbb_ZmmLowPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M115_ZmmLowPt_8TeV.root/vhbb_ZmmLowPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M120_ZmmLowPt_8TeV.root/vhbb_ZmmLowPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M125_ZmmLowPt_8TeV.root/vhbb_ZmmLowPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M130_ZmmLowPt_8TeV.root/vhbb_ZmmLowPt_8TeV.root/g' %s " %newcard)
       os.system("sed -i 's/vhbb_WS_BDT_M135_ZmmLowPt_8TeV.root/vhbb_ZmmLowPt_8TeV.root/g' %s " %newcard)

       os.system("echo %s written" % (newcard))

for item in files:

 if "110" in item:
  Process(item,0,0);
  Process(item,1,0);
  Process(item,2,0);
  Process(item,3,0);
  Process(item,4,0);
  
 if "115" in item:
  for x in range(5,15):
    Process(item, x , 10);
 
 if "120" in item:
  for x in range(15,25):
    Process(item, x , 20);

 
 if "125" in item:
  for x in range(25,35):
    Process(item, x , 30);
  for x in range(51,67):
    Process(item, x , 30);
 
 if "130" in item:
  for x in range(35,45):
    Process(item, x , 40);
 
 if "135" in item:
  for x in range(45,51):
    Process(item, x , 50);
 
