#! /usr/bin/env python

#setenv SVNGRP svn+ssh://svn.cern.ch/reps/cmshcg/trunk/ichep2012
#svn co $SVNGRP


import sys
import os
import commands
import string


mass=["110","110_5","111","111_5","112","112_5","113","113_5","114","114_5","115","115_5","116","116_5","117","117_5","118","118_5","119","119_5","120","120_5","121","121_5","122","122_5","123","123_5","124","124_5","125","125_5","126","126_5","127","127_5","128","128_5","129","129_5","130","130_5","131","131_5","132","132_5","133","133_5","134","134_5","135"]
massSVN=["110","110.5","111","111.5","112","112.5","113","113.5","114","114.5","115","115.5","116","116.5","117","117.5","118","118.5","119","119.5","120","120.5","121","121.5","122","122.5","123","123.5","124","124.5","125","125.5","126","126.5","127","127.5","128","128.5","129","129.5","130","130.5","131","131.5","132","132.5","133","133.5","134","134.5","135"]

def remove(n):
	os.system("rm %s/%s"%(path,n))
	os.system("svn rm %s/%s"%(path,n))
	
for x in range(0,51):
    print massSVN[x] 
    
    
    path=massSVN[x]  # insert the path to the directory of interest
    dirList=os.listdir(path)
    
    
    
  # for fname in dirList:
  #     if "Wmunu" in fname:
  #      remove(fname)
  #     if "Wenu" in fname:
  #      remove(fname)
  #     if "vhbb" in fname:
  #      remove(fname)
  #     if "WHbb" in fname:
  #      remove(fname)
  #     if "ZnnH" in fname:
  #      remove(fname)
    
   # os.system('svn commit -m "removed files"')
    
    
    os.system("cp ../../tmpArea/Zll/7TeV/%s/vhbb* %s" % (mass[x],massSVN[x]))
    os.system("cp ../../tmpArea/Zll/8TeV/%s/vhbb* %s" % (mass[x],massSVN[x]))
    
    os.system("cp ../../tmpArea/WH/7TeV/%s/vhbb* %s" % (mass[x],massSVN[x]))
    os.system("cp ../../tmpArea/WH/8TeV/%s/vhbb* %s" % (mass[x],massSVN[x]))
    
    os.system("cp ../../tmpArea/Znunu/7TeV/%s/vhbb* %s" % (mass[x],massSVN[x]))
    os.system("cp ../../tmpArea/Znunu/8TeV/%s/vhbb* %s" % (mass[x],massSVN[x]))
    
    os.system('svn add %s/vhbb_*' % massSVN[x])
    
    os.system('svn commit -m "added new vhbb files Mass %s"' % massSVN[x])
    
