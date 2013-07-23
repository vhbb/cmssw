echo "*** Zll ***" >! a1.log
combine -M Asymptotic -t -1 vhbb_Zll_8TeV.txt >> a1.log
echo "*** Znn ***" >! a2.log
combine -M Asymptotic -t -1 vhbb_Znn_8TeV.txt >> a2.log
echo "*** Wln ***" >! a3.log
combine -M Asymptotic -t -1 vhbb_Wln_8TeV.txt >> a3.log
echo "*** Wtn ***" >! a4.log
combine -M Asymptotic -t -1 vhbb_Wtn_8TeV.txt >> a4.log
echo "*** combo ***" >! a5.log
combine -M Asymptotic -t -1 vhbb_VH_8TeV.txt >> a5.log
echo "*** ZH ***" >! a6.log
combine -M Asymptotic -t -1 vhbb_ZH_8TeV.txt >> a6.log
echo "*** WH ***" >! a7.log
combine -M Asymptotic -t -1 vhbb_WH_8TeV.txt >> a7.log
echo "*** combo 7+8 ***" >! a8.log
combine -M Asymptotic -t -1 vhbb_VH_7p8TeV.txt >> a8.log

cat a1.log a2.log a3.log a4.log a5.log a6.log a7.log a8.log
