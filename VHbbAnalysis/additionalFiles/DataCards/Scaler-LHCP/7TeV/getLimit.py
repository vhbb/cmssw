#! /usr/bin/env python


import sys
import os
import commands
import string


#print sys.argv[1] 
files = [
"110/vhbb_Wln_8TeV.txt",
"115/vhbb_Wln_8TeV.txt",
"120/vhbb_Wln_8TeV.txt",
"125/vhbb_Wln_8TeV.txt",
"130/vhbb_Wln_8TeV.txt",
"135/vhbb_Wln_8TeV.txt"
]

i = 0;

for item in files:  
   os.system("/tmp/madfish/temp/Intermediate/CMSSW_5_2_4_patch4/bin/slc5_amd64_gcc462/combine -M Asymptotic %s |  grep \"r < \" | awk '{print $5}' > a%i" % (item,i))
   i = 1 + i

os.system("echo 'double y_observed[n_points]=        {'> left")
os.system("echo 'double y_down_points2[n_points]=    {'>> left")
os.system("echo 'double y_down_points1[n_points]=    {'>> left")
os.system("echo 'double y_vals[n_points]=            {'>> left")
os.system("echo 'double y_up_points1[n_points]=      {'>> left")
os.system("echo 'double y_up_points2[n_points]=      {'>> left")

os.system("echo ','>commas")
os.system("echo ','>>commas")
os.system("echo ','>>commas")
os.system("echo ','>>commas")
os.system("echo ','>>commas")
os.system("echo ','>>commas")

os.system("echo '};'>right")
os.system("echo '};'>>right")
os.system("echo '};'>>right")
os.system("echo '};'>>right")
os.system("echo '};'>>right")
os.system("echo '};'>>right")


os.system("paste left a0 commas a1 commas a2 commas a3 commas a4 commas a5 right > Result")
os.system("more Result")
os.system("rm -rf left right commas a0 a1 a2 a3 a4 a5")
