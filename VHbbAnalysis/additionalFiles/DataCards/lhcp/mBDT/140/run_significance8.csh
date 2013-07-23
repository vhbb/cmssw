combine -M ProfileLikelihood -m 140 --signif --pvalue -t -1 --expectSignal=1 vhbb_VH_8TeV.txt  > signPre8.txt
combine -M ProfileLikelihood -m 140 --signif --pvalue -t -1 --toysFreq --expectSignal=1 vhbb_VH_8TeV.txt  >  signPost8.txt
combine -M ProfileLikelihood -m 140 --signif --pvalue  vhbb_VH_8TeV.txt  >  signObs8.txt
