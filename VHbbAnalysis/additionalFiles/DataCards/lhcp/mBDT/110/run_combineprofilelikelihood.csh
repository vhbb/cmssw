echo "*** Zll ***" >! p1.log
combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --expectSignal=1 vhbb_Zll_8TeV.txt >> p1.log
echo "*** Znn ***" >! p2.log
combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --expectSignal=1 vhbb_Znn_8TeV.txt >> p2.log
echo "*** Wln ***" >! p3.log
combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --expectSignal=1 vhbb_Wln_8TeV.txt >> p3.log
echo "*** Wtn ***" >! p4.log
combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --expectSignal=1 vhbb_Wtn_8TeV.txt >> p4.log
echo "*** combo ***" >! p5.log
combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --expectSignal=1 vhbb_VH_8TeV.txt >> p5.log
echo "*** ZH ***" >! p6.log
combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --expectSignal=1 vhbb_ZH_8TeV.txt >> p6.log
echo "*** WH ***" >! p7.log
combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --expectSignal=1 vhbb_WH_8TeV.txt >> p7.log
echo "*** combo 7+8 ***" >! p8.log
combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --expectSignal=1 vhbb_VH_7p8TeV.txt >> p8.log

cat p1.log p2.log p3.log p4.log p5.log p6.log p7.log p8.log
