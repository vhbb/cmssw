#! /usr/bin/env python

import sys
import os
import commands
import string


mass=["110.root","110_5.root","111.root","111_5.root","112.root","112_5.root","113.root","113_5.root","114.root","114_5.root","115.root","115_5.root","116.root","116_5.root","117.root","117_5.root","118.root","118_5.root","119.root","119_5.root","120.root","120_5.root","121.root","121_5.root","122.root","122_5.root","123.root","123_5.root","124.root","124_5.root","125.root","125_5.root","126.root","126_5.root","127.root","127_5.root","128.root","128_5.root","129.root","129_5.root","130.root","130_5.root","131.root","131_5.root","132.root","132_5.root","133.root","133_5.root","134.root","134_5.root","135.root"]


os.system("rm *Tight*txt")
os.system("rm *Med*txt")
path="./"  # insert the path to the directory of interest
dirList=os.listdir(path)
for fname in dirList:
   if ".root" in fname:
     if "Zmm" in fname:
      if "Tight" in fname:
        os.system("mv %s %s" % (fname,"vhbb_Zmm2_7TeV.input.root"))
      if "Med" in fname:
        os.system("mv %s %s" % (fname,"vhbb_Zmm1_7TeV.input.root"))
     if "Zee" in fname:
      if "Tight" in fname:
        os.system("mv %s %s" % (fname,"vhbb_Zee2_7TeV.input.root"))
      if "Med" in fname:
        os.system("mv %s %s" % (fname,"vhbb_Zee1_7TeV.input.root"))
 
   if ".txt" in fname:
        for item in mass:
          os.system("sed -i 's@vhbb_WS_BDTRTight_M110_Zmm_7TeV%s@vhbb_Zmm2_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRTight_M115_Zmm_7TeV%s@vhbb_Zmm2_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRTight_M120_Zmm_7TeV%s@vhbb_Zmm2_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRTight_M125_Zmm_7TeV%s@vhbb_Zmm2_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRTight_M130_Zmm_7TeV%s@vhbb_Zmm2_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRTight_M135_Zmm_7TeV%s@vhbb_Zmm2_7TeV.input.root@g' %s" % (item,fname))

          os.system("sed -i 's@vhbb_WS_BDTRMed_M110_Zmm_7TeV%s@vhbb_Zmm1_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRMed_M115_Zmm_7TeV%s@vhbb_Zmm1_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRMed_M120_Zmm_7TeV%s@vhbb_Zmm1_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRMed_M125_Zmm_7TeV%s@vhbb_Zmm1_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRMed_M130_Zmm_7TeV%s@vhbb_Zmm1_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRMed_M135_Zmm_7TeV%s@vhbb_Zmm1_7TeV.input.root@g' %s" % (item,fname))

          os.system("sed -i 's@vhbb_WS_BDTRTight_M110_Zee_7TeV%s@vhbb_Zee2_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRTight_M115_Zee_7TeV%s@vhbb_Zee2_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRTight_M120_Zee_7TeV%s@vhbb_Zee2_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRTight_M125_Zee_7TeV%s@vhbb_Zee2_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRTight_M130_Zee_7TeV%s@vhbb_Zee2_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRTight_M135_Zee_7TeV%s@vhbb_Zee2_7TeV.input.root@g' %s" % (item,fname))

          os.system("sed -i 's@vhbb_WS_BDTRMed_M110_Zee_7TeV%s@vhbb_Zee1_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRMed_M115_Zee_7TeV%s@vhbb_Zee1_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRMed_M120_Zee_7TeV%s@vhbb_Zee1_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRMed_M125_Zee_7TeV%s@vhbb_Zee1_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRMed_M130_Zee_7TeV%s@vhbb_Zee1_7TeV.input.root@g' %s" % (item,fname))
          os.system("sed -i 's@vhbb_WS_BDTRMed_M135_Zee_7TeV%s@vhbb_Zee1_7TeV.input.root@g' %s" % (item,fname))
        os.system("mv %s vhbb_Zll_7TeV.txt" % fname) 
          
