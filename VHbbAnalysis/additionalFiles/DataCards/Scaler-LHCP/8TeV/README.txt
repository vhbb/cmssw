To run WS .root scaler:

add .root files to IntermediateMassMaker.py (files = [] list), the filenames (or paths) should contain mass value (ie 110,115,120..) & Channel (Wen,Wmn,Zll,Znn,Wtn).

then: python IntermediateMassMaker.py

you can use adder.py to add Low/High/Medium separated cards into 1 WS

============================================================================

To run DC .txt scaler:

same procedure as before except editing scaler.py


then: python scaler.py

=============================================================================


tester.py can test all datacards at once 

To use for 7TeV change all Xsec usages to use 7TeV arrays (provided in scaler.py
and xsec7TeV.h).
