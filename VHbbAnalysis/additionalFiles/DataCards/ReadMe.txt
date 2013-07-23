0) copy ALL files from UserCode/VHbbAnalysis/additionalFiles/DataCards/ into a
folder 

1) mkdir 110 110.5 111 111.5 112 112.5 113 113.5 114 114.5 115 115.5 116 116.5 117 117.5 118 118.5 119 119.5 120 120.5 121 121.5 122 122.5 123 123.5 124 124.5 125 125.5 126 126.5 127 127.5 128 128.5 129 129.5 130 130.5 131 131.5 132 132.5 133 133.5 134 134.5 135 124.6 124.7 124.8 124.9 125.1 125.2 125.3 125.4 125.6 125.7 125.8 125.9 126.1 126.2 126.3 126.4
   
2) IntermediateMassMaker - 
   Change:  const int bins = X
   to the number of bins in your histos
   YOU SHOULD PLOT YOUR HISTO IN THE ACTUAL DATACARD TO VERIFY THIS VALUE!


3) Edit Systematic names properly!
   s < x (x is number of systs)
   c < y (y is the number of channels)

4) Copy all .root files into the same directory as the script (should contain
mass points in the filenames) Edit filesXXX.h (XXX = Zll, Znn, Wln) to point
to these files and correct the array size


4) run IntermediateMassMaker - 

   source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.02/x86_64-slc5-gcc43-opt/root/bin/thisroot.csh
   source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.csh

   Be careful not to run all 67 mass point loops at once, you will likely crash the computer. I'd say to run 1 or 2 loops (mass pt DCs) max  at a time.
   
   RUN IT COMPILED, IT IS 10X FASTER IE: root IntermediateMassMakerXX.C++   

5) edit scalerXXX.py:
   a)   set the proper energy and channel all xsecs are available, is xSec8ZH, xSec7ZH, xSec7WH, or xSec8WH
   b) edit files =[ ] to datacards in order of masses
   c) correct the sed (replacement) commands to edit cards .root WS to new official formats
      (the ones from IntermediateMassMaker that do not have Masses in the
filenames)
   d) check that signal columns are correct (awk '{print $x}' commands have x
 as signal columns!

6) run scalerXXX.py

7) edit/run tester.py to test all DCs

8) To submit:
 
first edit submitter.py to point to the correct source folders
(kinit)

mkdir HCP2012
cp submitter.py HCP2012/
cd HCP2012
setenv SVNGRP svn+ssh://svn.cern.ch/reps/cmshcg
svn co $SVNGRP/trunk/hcp2012 
svn update (IMPORTANT - ALWAYS DO THIS BEFORE COMMIT)
python submitter.py
svn commit
 
