To run WS .root scaler:

add .root files to IntermediateMassMaker.py (files = [] list), the filenames (or paths) should contain mass value (ie 110,115,120..) & Channel (Wen,Wmn,Zll,Znn,Wtn).

then: python IntermediateMassMaker.py



To run DC .txt scaler:

same procedure as before except editing scaler.py


then: python scaler.py


To use for 7TeV change all Xsec usages to use 7TeV arrays (provided in scaler.py
and xsec7TeV.h).
