echo "*** Zll ***" >! A1.log
combine -M Asymptotic vhbb_Zll_8TeV.txt >> A1.log
echo "*** Znn ***" >! A2.log
combine -M Asymptotic vhbb_Znn_8TeV.txt >> A2.log
echo "*** Wln ***" >! A3.log
combine -M Asymptotic vhbb_Wln_8TeV.txt >> A3.log
echo "*** Wtn ***" >! A4.log
combine -M Asymptotic vhbb_Wtn_8TeV.txt >> A4.log
echo "*** combo ***" >! A5.log
combine -M Asymptotic vhbb_VH_8TeV.txt >> A5.log
echo "*** ZH ***" >! A6.log
combine -M Asymptotic vhbb_ZH_8TeV.txt >> A6.log
echo "*** WH ***" >! A7.log
combine -M Asymptotic vhbb_WH_8TeV.txt >> A7.log
echo "*** combo 7+8 ***" >! A8.log
combine -M Asymptotic vhbb_VH_7p8TeV.txt >> A8.log

cat A1.log A2.log A3.log A4.log A5.log A6.log A7.log A8.log
