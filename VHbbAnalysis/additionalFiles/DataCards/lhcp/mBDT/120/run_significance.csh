combine -M ProfileLikelihood -m 120 --signif --pvalue -t -1 --expectSignal=1 vhbb_VH_7p8TeV.txt  > signPre.txt
combine -M ProfileLikelihood -m 120 --signif --pvalue -t -1 --toysFreq --expectSignal=1 vhbb_VH_7p8TeV.txt  >  signPost.txt
combine -M ProfileLikelihood -m 120 --signif --pvalue  vhbb_VH_7p8TeV.txt  >  signObs.txt
