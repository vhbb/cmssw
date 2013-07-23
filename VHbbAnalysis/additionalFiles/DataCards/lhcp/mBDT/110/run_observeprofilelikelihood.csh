echo "*** Zll ***" >! P1.log
combine -M ProfileLikelihood -m 125 --signif --pvalue vhbb_Zll_8TeV.txt >> P1.log
echo "*** Znn ***" >! P2.log
combine -M ProfileLikelihood -m 125 --signif --pvalue vhbb_Znn_8TeV.txt >> P2.log
echo "*** Wln ***" >! P3.log
combine -M ProfileLikelihood -m 125 --signif --pvalue vhbb_Wln_8TeV.txt >> P3.log
echo "*** Wtn ***" >! P4.log
combine -M ProfileLikelihood -m 125 --signif --pvalue vhbb_Wtn_8TeV.txt >> P4.log
echo "*** combo ***" >! P5.log
combine -M ProfileLikelihood -m 125 --signif --pvalue vhbb_VH_8TeV.txt >> P5.log
echo "*** ZH ***" >! P6.log
combine -M ProfileLikelihood -m 125 --signif --pvalue vhbb_ZH_8TeV.txt >> P6.log
echo "*** WH ***" >! P7.log
combine -M ProfileLikelihood -m 125 --signif --pvalue vhbb_WH_8TeV.txt >> P7.log
echo "*** combo 7+8 ***" >! P8.log
combine -M ProfileLikelihood -m 125 --signif --pvalue vhbb_VH_7p8TeV.txt >> P8.log

cat P1.log P2.log P3.log P4.log P5.log P6.log P7.log P8.log
