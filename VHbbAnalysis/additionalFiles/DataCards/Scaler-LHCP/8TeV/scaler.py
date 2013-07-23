#! /usr/bin/env python


import sys
import os
import commands
import string


files=[
#"110/vhbb_Wln_8TeV.txt",
#"110/vhbb_Wtn_8TeV.txt",
"110/vhbb_Zll_8TeV.txt",
#"110/vhbb_Znn_8TeV.txt",
#"115/vhbb_Wln_8TeV.txt",
#"115/vhbb_Wtn_8TeV.txt",
"115/vhbb_Zll_8TeV.txt",
#"115/vhbb_Znn_8TeV.txt",
#"120/vhbb_Wln_8TeV.txt",
#"120/vhbb_Wtn_8TeV.txt",
"120/vhbb_Zll_8TeV.txt",
#"120/vhbb_Znn_8TeV.txt",
#"125/vhbb_Wln_8TeV.txt",
#"125/vhbb_Wtn_8TeV.txt",
"125/vhbb_Zll_8TeV.txt",
#"125/vhbb_Znn_8TeV.txt",
#"130/vhbb_Wln_8TeV.txt",
#"130/vhbb_Wtn_8TeV.txt",
"130/vhbb_Zll_8TeV.txt",
#"130/vhbb_Znn_8TeV.txt",
#"135/vhbb_Wln_8TeV.txt",
#"135/vhbb_Wtn_8TeV.txt",
"135/vhbb_Zll_8TeV.txt",
#"135/vhbb_Znn_8TeV.txt",
#"140/vhbb_Wln_8TeV.txt",
#"140/vhbb_Wtn_8TeV.txt",
"140/vhbb_Zll_8TeV.txt",
#"140/vhbb_Znn_8TeV.txt",
#"145/vhbb_Wln_8TeV.txt",
#"145/vhbb_Wtn_8TeV.txt",
"145/vhbb_Zll_8TeV.txt",
#"145/vhbb_Znn_8TeV.txt",
#"150/vhbb_Wln_8TeV.txt",
#"150/vhbb_Wtn_8TeV.txt",
"150/vhbb_Zll_8TeV.txt",
#"150/vhbb_Znn_8TeV.txt"
]



xSec7ZH = [0.4721,0.4655,0.4589,0.4525,0.4462,0.4400,0.4340,0.4280,0.4221,0.4164, 0.4107, 0.4052, 0.3998, 0.3945, 0.3893, 0.3842, 0.3791, 0.3742, 0.3693, 0.3645, 0.3598, 0.3551, 0.3505, 0.3459, 0.3414, 0.3370, 0.3326, 0.3283, 0.3241, 0.3199, 0.3158, 0.3117, 0.3077, 0.3038, 0.2999, 0.2961, 0.2923, 0.2886, 0.2849, 0.2813, 0.2778, 0.2743, 0.2709, 0.2675, 0.2642, 0.2609, 0.2577, 0.2545, 0.2514, 0.2483, 0.2453, 0.31908, 0.31826, 0.31744, 0.31662, 0.314980, 0.314160, 0.313340, 0.312520, 0.310900, 0.310100, 0.309300, 0.308500, 0.306920, 0.306140, 0.305360, 0.304580]

xSec7WH  = [0.8754,0.8623,0.8495,0.8368,0.8244,	0.8122,	0.8003,	0.7885,	0.7770,	0.7657,	0.7546,	0.7439,	0.7333,	0.7230,	0.7129,	0.7030,	0.6933,	0.6837,	0.6744,	0.6651,	0.6561,	0.6472,	0.6384,	0.6297,	0.6212,	0.6129,	0.6046,	0.5965,	0.5885,	0.5806,	0.5729,	0.5652,	0.5576,	0.5501,	0.5428,	0.5355,	0.5284,	0.5213,	0.5144,	0.5075,	0.5008,	0.4942,	0.4877,	0.4813,	0.4749,	0.4687,	0.4626,	0.4566,	0.4506,	0.4448,	0.4390,0.579060,0.577520,0.575980,0.574440,0.571360,0.569820,0.568280,0.566740,0.563680,0.562160,0.560640,0.559120,0.556100,0.554600,0.553100,0.551600]

xSec8ZH = [0.5869, 0.5788, 0.5708, 0.5629, 0.5552, 0.5476, 0.5402, 0.5329, 0.5258, 0.5187, 0.5117, 0.5049, 0.4981, 0.4916, 0.4850, 0.4787, 0.4724, 0.4662, 0.4602, 0.4542, 0.4483, 0.4426, 0.4368, 0.4312, 0.4257, 0.4203, 0.4150, 0.4096, 0.4044, 0.3993, 0.3943, 0.3893, 0.3843, 0.3794, 0.3746, 0.3699, 0.3652, 0.3606, 0.3561, 0.3516, 0.3473, 0.3430, 0.3388, 0.3347, 0.3306, 0.3266, 0.3226, 0.3188, 0.3149, 0.3112, 0.3074, 0.3038, 0.3001, 0.2966, 0.2930, 0.2895, 0.2861, 0.2827, 0.2793, 0.2760, 0.2728, 0.2696, 0.2664, 0.2632, 0.2601, 0.2571, 0.2541, 0.2511, 0.2482, 0.2453, 0.2424, 0.2396, 0.2368, 0.2341, 0.2314, 0.2287, 0.2261, 0.2235, 0.2209, 0.2184, 0.2159, 0.39830, 0.39730, 0.39630, 0.39530, 0.39330, 0.39230, 0.39130, 0.39030, 0.38830, 0.38730, 0.38630, 0.38530, 0.38332, 0.38234, 0.38136, 0.38038]

xSec8WH = [1.060, 1.045, 1.030, 1.015, 0.9998, 0.9852, 0.9709, 0.9570, 0.9432, 0.9297, 0.9165, 0.9035, 0.8907, 0.8782, 0.8659, 0.8538, 0.8420, 0.8303, 0.8187, 0.8075, 0.7966, 0.7859, 0.7753, 0.7649, 0.7547, 0.7446, 0.7347, 0.7249, 0.7154, 0.7060, 0.6966, 0.6873, 0.6782, 0.6691, 0.6602, 0.6515, 0.6429, 0.6344, 0.6260, 0.6177, 0.6095, 0.6015, 0.5936, 0.5859, 0.5783, 0.5708, 0.5634, 0.5562, 0.5491, 0.5420, 0.5351, 0.5283, 0.5215, 0.5149, 0.5084, 0.5020, 0.4956, 0.4894, 0.4833, 0.4772, 0.4713, 0.4655, 0.4597, 0.4541, 0.4484, 0.4430, 0.4375, 0.4322, 0.4268, 0.4216, 0.4164, 0.4113, 0.4062, 0.4012, 0.3963, 0.3915, 0.3867, 0.3820, 0.3773, 0.3727, 0.3681, 0.70412, 0.70224, 0.70036, 0.69848, 0.69474, 0.69288, 0.69102, 0.68916, 0.68548, 0.68366, 0.68184, 0.68002, 0.67638, 0.67456, 0.67274, 0.67092 ]

mass=["110","110.5","111","111.5","112","112.5","113","113.5","114","114.5", "115","115.5","116","116.5","117","117.5","118","118.5","119","119.5", "120","120.5","121","121.5","122","122.5","123","123.5","124","124.5", "125","125.5","126","126.5","127","127.5","128","128.5","129","129.5", "130","130.5","131","131.5","132","132.5","133","133.5","134","134.5", "135","135.5","136","136.5","137","137.5","138","138.5","139","139.5", "140","140.5","141","141.5","142","142.5","143","143.5","144","144.5", "145","145.5","146","146.5","147","147.5","148","148.5","149","149.5", "150", "124.6","124.7","124.8","124.9","125.1","125.2","125.3","125.4", "125.6","125.7","125.8","125.9","126.1","126.2","126.3","126.4"]

def ProcessZll(file, toMass, fromMass):
       newcard = "newcards/" + mass[toMass] + "/vhbb_Zll_8TeV.txt"
       os.system('more %s | grep rate > a' % file)
       os.system("more a |  awk '{print $2}'  > aVH1")
       os.system("more a |  awk '{print $10}'  > aVH2")
       os.system("more a |  awk '{print $18}'  > aVH3")
       os.system("more a |  awk '{print $26}'  > aVH4")
 
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


       bVH1 = (xSec8ZH[toMass]/xSec8ZH[fromMass])*VH1
       bVH2 = (xSec8ZH[toMass]/xSec8ZH[fromMass])*VH2
       bVH3 = (xSec8ZH[toMass]/xSec8ZH[fromMass])*VH3
       bVH4 = (xSec8ZH[toMass]/xSec8ZH[fromMass])*VH4
    
       os.system("sed 's/%s/%f/g' a > b" % (sVH1,bVH1))
       os.system("sed -i 's/%s/%f/g' b" % (sVH2,bVH2))
       os.system("sed -i 's/%s/%f/g' b" % (sVH3,bVH3))
       os.system("sed -i 's/%s/%f/g' b" % (sVH4,bVH4))

       fRO = open("a", "r")
       sRO = fRO.read().rstrip('\n')
       fRF = open("b", "r")
       sRF = fRF.read().rstrip('\n')
       os.system("sed 's/%s/%s/g' %s > %s" % (sRO,sRF,file,newcard))

       os.system("sed -i 's/_WS_BDT_//g' %s " %newcard)
       os.system("sed -i 's/M110//g' %s " %newcard)
       os.system("sed -i 's/M115//g' %s " %newcard)
       os.system("sed -i 's/M120//g' %s " %newcard)
       os.system("sed -i 's/M125//g' %s " %newcard)
       os.system("sed -i 's/M130//g' %s " %newcard)
       os.system("sed -i 's/M135//g' %s " %newcard)
       os.system("sed -i 's/M140//g' %s " %newcard)
       os.system("sed -i 's/M145//g' %s " %newcard)
       os.system("sed -i 's/M150//g' %s " %newcard)

       os.system("echo %s written" % (newcard))



def ProcessWln(file, toMass, fromMass):
       newcard = "newcards/" + mass[toMass] + "/vhbb_Wln_8TeV.txt"
       os.system('more %s | grep rate > a' % file)
       os.system("more a |  awk '{print $2} '  > aZH1")
       os.system("more a |  awk '{print $3}'   > aWH1")
       os.system("more a |  awk '{print $15}'  > aZH2")
       os.system("more a |  awk '{print $16}'  > aWH2")
       os.system("more a |  awk '{print $28}'  > aZH3")
       os.system("more a |  awk '{print $29}'  > aWH3")
       os.system("more a |  awk '{print $41}'  > aZH4")
       os.system("more a |  awk '{print $42}'  > aWH4")
       os.system("more a |  awk '{print $54}'  > aZH5")
       os.system("more a |  awk '{print $55}'  > aWH5")
       os.system("more a |  awk '{print $67}'  > aZH6")
       os.system("more a |  awk '{print $68}'  > aWH6")
 
       fZH1 = open("aZH1", "r")
       fZH2 = open("aZH2", "r")
       fZH3 = open("aZH3", "r")
       fZH4 = open("aZH4", "r")
       fZH5 = open("aZH5", "r")
       fZH6 = open("aZH6", "r")
       fWH1 = open("aWH1", "r")
       fWH2 = open("aWH2", "r")
       fWH3 = open("aWH3", "r")
       fWH4 = open("aWH4", "r")
       fWH5 = open("aWH5", "r")
       fWH6 = open("aWH6", "r")

       sZH1 = fZH1.read().rstrip('\n')
       sZH2 = fZH2.read().rstrip('\n')
       sZH3 = fZH3.read().rstrip('\n')
       sZH4 = fZH4.read().rstrip('\n')
       sZH5 = fZH5.read().rstrip('\n')
       sZH6 = fZH6.read().rstrip('\n')
       sWH1 = fWH1.read().rstrip('\n')
       sWH2 = fWH2.read().rstrip('\n')
       sWH3 = fWH3.read().rstrip('\n')
       sWH4 = fWH4.read().rstrip('\n')
       sWH5 = fWH5.read().rstrip('\n')
       sWH6 = fWH6.read().rstrip('\n')


       ZH1 = float(sZH1)
       ZH2 = float(sZH2)
       ZH3 = float(sZH3)
       ZH4 = float(sZH4)
       ZH5 = float(sZH5)
       ZH6 = float(sZH6)

       WH1 = float(sWH1)
       WH2 = float(sWH2)
       WH3 = float(sWH3)
       WH4 = float(sWH4)
       WH5 = float(sWH5)
       WH6 = float(sWH6)


       bZH1 = (xSec8ZH[toMass]/xSec8ZH[fromMass])*ZH1
       bZH2 = (xSec8ZH[toMass]/xSec8ZH[fromMass])*ZH2
       bZH3 = (xSec8ZH[toMass]/xSec8ZH[fromMass])*ZH3
       bZH4 = (xSec8ZH[toMass]/xSec8ZH[fromMass])*ZH4
       bZH5 = (xSec8ZH[toMass]/xSec8ZH[fromMass])*ZH5
       bZH6 = (xSec8ZH[toMass]/xSec8ZH[fromMass])*ZH6

       bWH1 = (xSec8WH[toMass]/xSec8WH[fromMass])*WH1
       bWH2 = (xSec8WH[toMass]/xSec8WH[fromMass])*WH2
       bWH3 = (xSec8WH[toMass]/xSec8WH[fromMass])*WH3
       bWH4 = (xSec8WH[toMass]/xSec8WH[fromMass])*WH4
       bWH5 = (xSec8WH[toMass]/xSec8WH[fromMass])*WH5
       bWH6 = (xSec8WH[toMass]/xSec8WH[fromMass])*WH6



       os.system("sed 's/%s/%f/g' a > b" %(sZH1,bZH1))
       os.system("sed -i 's/%s/%f/g' b" % (sZH2,bZH2))
       os.system("sed -i 's/%s/%f/g' b" % (sZH3,bZH3))
       os.system("sed -i 's/%s/%f/g' b" % (sZH4,bZH4))
       os.system("sed -i 's/%s/%f/g' b" % (sZH5,bZH5))
       os.system("sed -i 's/%s/%f/g' b" % (sZH6,bZH6))

       os.system("sed -i 's/%s/%f/g' b" % (sWH1,bWH1))
       os.system("sed -i 's/%s/%f/g' b" % (sWH2,bWH2))
       os.system("sed -i 's/%s/%f/g' b" % (sWH3,bWH3))
       os.system("sed -i 's/%s/%f/g' b" % (sWH4,bWH4))
       os.system("sed -i 's/%s/%f/g' b" % (sWH5,bWH5))
       os.system("sed -i 's/%s/%f/g' b" % (sWH6,bWH6))


       fRO = open("a", "r")
       sRO = fRO.read().rstrip('\n')
       fRF = open("b", "r")
       sRF = fRF.read().rstrip('\n')
       os.system("sed 's/%s/%s/g' %s > %s" % (sRO,sRF,file,newcard))

       os.system("sed -i 's/_WS_BDT_//g' %s " %newcard)
       os.system("sed -i 's/H110//g' %s " %newcard)
       os.system("sed -i 's/H115//g' %s " %newcard)
       os.system("sed -i 's/H120//g' %s " %newcard)
       os.system("sed -i 's/H125//g' %s " %newcard)
       os.system("sed -i 's/H130//g' %s " %newcard)
       os.system("sed -i 's/H135//g' %s " %newcard)
       os.system("sed -i 's/H140//g' %s " %newcard)
       os.system("sed -i 's/H145//g' %s " %newcard)
       os.system("sed -i 's/H150//g' %s " %newcard)

       os.system("echo %s written" % (newcard))


def ProcessZnn(file, toMass, fromMass):
       newcard = "newcards/" + mass[toMass] + "/vhbb_Znn_8TeV.txt"
       os.system('more %s | grep rate > a' % file)
       os.system("more a |  awk '{print $2}'  > aZH1")
       os.system("more a |  awk '{print $3}'  > aWH1")
       os.system("more a |  awk '{print $15}'  > aZH2")
       os.system("more a |  awk '{print $16}'  > aWH2")
       os.system("more a |  awk '{print $28}'  > aZH3")
       os.system("more a |  awk '{print $29}'  > aWH3")
 
       fZH1 = open("aZH1", "r")
       fZH2 = open("aZH2", "r")
       fZH3 = open("aZH3", "r")
       fWH1 = open("aWH1", "r")
       fWH2 = open("aWH2", "r")
       fWH3 = open("aWH3", "r")

       sZH1 = fZH1.read().rstrip('\n')
       sZH2 = fZH2.read().rstrip('\n')
       sZH3 = fZH3.read().rstrip('\n')
       sWH1 = fWH1.read().rstrip('\n')
       sWH2 = fWH2.read().rstrip('\n')
       sWH3 = fWH3.read().rstrip('\n')



       ZH1 = float(sZH1)
       ZH2 = float(sZH2)
       ZH3 = float(sZH3)
       WH1 = float(sWH1)
       WH2 = float(sWH2)
       WH3 = float(sWH3)

       bZH1 = (xSec8ZH[toMass]/xSec8ZH[fromMass])*ZH1
       bZH2 = (xSec8ZH[toMass]/xSec8ZH[fromMass])*ZH2
       bZH3 = (xSec8ZH[toMass]/xSec8ZH[fromMass])*ZH3
       bWH1 = (xSec8WH[toMass]/xSec8WH[fromMass])*WH1
       bWH2 = (xSec8WH[toMass]/xSec8WH[fromMass])*WH2
       bWH3 = (xSec8WH[toMass]/xSec8WH[fromMass])*WH3


       os.system("sed 's/%s/%f/g' a > b" % (sWH1,bWH1))
       os.system("sed -i 's/%s/%f/g' b" % (sWH2,bWH2))
       os.system("sed -i 's/%s/%f/g' b" % (sWH3,bWH3))
       os.system("sed -i 's/%s/%f/g' b" % (sZH1,bZH1))
       os.system("sed -i 's/%s/%f/g' b" % (sZH2,bZH2))
       os.system("sed -i 's/%s/%f/g' b" % (sZH3,bZH3))

       fRO = open("a", "r")
       sRO = fRO.read().rstrip('\n')
       fRF = open("b", "r")
       sRF = fRF.read().rstrip('\n')
       os.system("sed 's/%s/%s/g' %s > %s" % (sRO,sRF,file,newcard))
       os.system("sed -i 's/vhbb_Znn_J12_8TeV.root/vhbb_Znn_8TeV.root/g' %s " % (newcard))
       os.system("echo %s written" % (newcard))

def ProcessWtn(file, toMass, fromMass):
       newcard = "newcards/" + mass[toMass] + "/vhbb_Wtn_8TeV.txt"
       os.system('more %s | grep rate > a' % file)
       os.system("more a |  awk '{print $2} '  > aZH1")
       os.system("more a |  awk '{print $3}'   > aWH1")
 
       fZH1 = open("aZH1", "r")
       fWH1 = open("aWH1", "r")

       sZH1 = fZH1.read().rstrip('\n')
       sWH1 = fWH1.read().rstrip('\n')


       ZH1 = float(sZH1)
       WH1 = float(sWH1)


       bZH1 = (xSec8ZH[toMass]/xSec8ZH[fromMass])*ZH1
       bWH1 = (xSec8WH[toMass]/xSec8WH[fromMass])*WH1

       os.system("sed 's/%s/%f/g' a > b" %(sZH1,bZH1))
       os.system("sed -i 's/%s/%f/g' b" % (sWH1,bWH1))

       fRO = open("a", "r")
       sRO = fRO.read().rstrip('\n')
       fRF = open("b", "r")
       sRF = fRF.read().rstrip('\n')
       os.system("sed 's/%s/%s/g' %s > %s" % (sRO,sRF,file,newcard))

       os.system("sed -i 's/Wtn_BDT_newBinning_110.root/vhbb_Wtn_8TeV.root/g' %s " % (newcard))
       os.system("sed -i 's/Wtn_BDT_newBinning_115.root/vhbb_Wtn_8TeV.root/g' %s " % (newcard))
       os.system("sed -i 's/Wtn_BDT_newBinning_120.root/vhbb_Wtn_8TeV.root/g' %s " % (newcard))
       os.system("sed -i 's/Wtn_BDT_newBinning_125.root/vhbb_Wtn_8TeV.root/g' %s " % (newcard))
       os.system("sed -i 's/Wtn_BDT_newBinning_130.root/vhbb_Wtn_8TeV.root/g' %s " % (newcard))
       os.system("sed -i 's/Wtn_BDT_newBinning_135.root/vhbb_Wtn_8TeV.root/g' %s " % (newcard))
       os.system("sed -i 's/Wtn_BDT_newBinning_140.root/vhbb_Wtn_8TeV.root/g' %s " % (newcard))
       os.system("sed -i 's/Wtn_BDT_newBinning_145.root/vhbb_Wtn_8TeV.root/g' %s " % (newcard))
       os.system("sed -i 's/Wtn_BDT_newBinning_150.root/vhbb_Wtn_8TeV.root/g' %s " % (newcard))

       os.system("echo %s written" % (newcard))



for item in mass: 
  os.system("mkdir -p newcards/%s" % item)


for item in files:

 if "Wln" in item:

  if "110" in item:
   ProcessWln(item,0,0);
   ProcessWln(item,1,0);
   ProcessWln(item,2,0);
   ProcessWln(item,3,0);
   ProcessWln(item,4,0);
  
  if "115" in item:
   for x in range(5,15):
     ProcessWln(item, x , 10);
 
  if "120" in item:
   for x in range(15,25):
     ProcessWln(item, x , 20);
 
  if "125" in item:
   for x in range(25,35):
     ProcessWln(item, x , 30);
   for x in range(81,97):
     ProcessWln(item, x , 30);
 
  if "130" in item:
   for x in range(35,45):
     ProcessWln(item, x , 40);
  
  if "135" in item:
   for x in range(45,55):
     ProcessWln(item, x , 50);
 
  if "140" in item:
   for x in range(55,65):
     ProcessWln(item, x , 60);
 
  if "145" in item:
   for x in range(65,75):
     ProcessWln(item, x , 70);

  if "150" in item:
   for x in range(75,81):
     ProcessWln(item, x , 80);

 if "Zll" in item:

   if "110" in item:
    ProcessZll(item,0,0);
    ProcessZll(item,1,0);
    ProcessZll(item,2,0);
    ProcessZll(item,3,0);
    ProcessZll(item,4,0);
    
   if "115" in item:
    for x in range(5,15):
      ProcessZll(item, x , 10);
   
   if "120" in item:
    for x in range(15,25):
      ProcessZll(item, x , 20);
   
   if "125" in item:
    for x in range(25,35):
      ProcessZll(item, x , 30);
    for x in range(81,97):
      ProcessZll(item, x , 30);
   
   if "130" in item:
    for x in range(35,45):
      ProcessZll(item, x , 40);
   
   if "135" in item:
    for x in range(45,55):
      ProcessZll(item, x , 50);
   
   if "140" in item:
    for x in range(55,65):
      ProcessZll(item, x , 60);
   
   if "145" in item:
    for x in range(65,75):
      ProcessZll(item, x , 70);
  
   if "150" in item:
    for x in range(75,81):
      ProcessZll(item, x , 80);
  
 if "Znn" in item:

   if "110" in item:
    ProcessZnn(item,0,0);
    ProcessZnn(item,1,0);
    ProcessZnn(item,2,0);
    ProcessZnn(item,3,0);
    ProcessZnn(item,4,0);
    
   if "115" in item:
    for x in range(5,15):
      ProcessZnn(item, x , 10);
   
   if "120" in item:
    for x in range(15,25):
      ProcessZnn(item, x , 20);
   
   if "125" in item:
    for x in range(25,35):
      ProcessZnn(item, x , 30);
    for x in range(81,97):
      ProcessZnn(item, x , 30);
   
   if "130" in item:
    for x in range(35,45):
      ProcessZnn(item, x , 40);
   
   if "135" in item:
    for x in range(45,55):
      ProcessZnn(item, x , 50);
   
   if "140" in item:
    for x in range(55,65):
      ProcessZnn(item, x , 60);
   
   if "145" in item:
    for x in range(65,75):
      ProcessZnn(item, x , 70);
  
   if "150" in item:
    for x in range(75,81):
      ProcessZnn(item, x , 80);
  
 if "Wtn" in item:

   if "110" in item:
    ProcessWtn(item,0,0);
    ProcessWtn(item,1,0);
    ProcessWtn(item,2,0);
    ProcessWtn(item,3,0);
    ProcessWtn(item,4,0);
    
   if "115" in item:
    for x in range(5,15):
      ProcessWtn(item, x , 10);
   
   if "120" in item:
    for x in range(15,25):
      ProcessWtn(item, x , 20);
   
   if "125" in item:
    for x in range(25,35):
      ProcessWtn(item, x , 30);
    for x in range(81,97):
      ProcessWtn(item, x , 30);
   
   if "130" in item:
    for x in range(35,45):
      ProcessWtn(item, x , 40);
   
   if "135" in item:
    for x in range(45,55):
      ProcessWtn(item, x , 50);
   
   if "140" in item:
    for x in range(55,65):
      ProcessWtn(item, x , 60);
   
   if "145" in item:
    for x in range(65,75):
      ProcessWtn(item, x , 70);
  
   if "150" in item:
    for x in range(75,81):
      ProcessWtn(item, x , 80);
  
             
            
             
            
             
            
             
            
  

