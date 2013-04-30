echo "*** Zll ***" >! m1.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_Zll_8TeV.txt > /dev/null
cp mlfit.root mlfit_Zll.root
echo "*** Znn ***" >! m2.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_Znn_8TeV.txt > /dev/null
cp mlfit.root mlfit_Znn.root
echo "*** Wln ***" >! m3.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_Wln_8TeV.txt > /dev/null
cp mlfit.root mlfit_Wln.root
echo "*** Wtn ***" >! m4.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_Wtn_8TeV.txt > /dev/null
cp mlfit.root mlfit_Wtn.root
echo "*** combo ***" >! m5.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_VH_8TeV.txt > /dev/null
cp mlfit.root mlfit_VH.root
echo "*** ZH ***" >! m6.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_ZH_8TeV.txt > /dev/null
cp mlfit.root mlfit_ZH.root
echo "*** WH ***" >! m7.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_WH_8TeV.txt > /dev/null
cp mlfit.root mlfit_WH.root
echo "*** combo 7+8 ***" >! m8.log
#combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_VH_7p8TeV.txt > /dev/null
#cp mlfit.root mlfit_VH78.root

cat m1.log m2.log m3.log m4.log m5.log m6.log m7.log m8.log

python diffNuisances.py mlfit_Zll.root -f html > ~/www/vhbb_LHCP_preapproval/nui_20130422/mlfit_mBDT_Zll_8TeV.html
python diffNuisances.py mlfit_Znn.root -f html > ~/www/vhbb_LHCP_preapproval/nui_20130422/mlfit_mBDT_Znn_8TeV.html
python diffNuisances.py mlfit_Wln.root -f html > ~/www/vhbb_LHCP_preapproval/nui_20130422/mlfit_mBDT_Wln_8TeV.html
python diffNuisances.py mlfit_Wtn.root -f html > ~/www/vhbb_LHCP_preapproval/nui_20130422/mlfit_mBDT_Wtn_8TeV.html
python diffNuisances.py mlfit_VH.root  -f html > ~/www/vhbb_LHCP_preapproval/nui_20130422/mlfit_mBDT_VH_8TeV.html
python diffNuisances.py mlfit_ZH.root  -f html > ~/www/vhbb_LHCP_preapproval/nui_20130422/mlfit_mBDT_ZH_8TeV.html
python diffNuisances.py mlfit_WH.root  -f html > ~/www/vhbb_LHCP_preapproval/nui_20130422/mlfit_mBDT_WH_8TeV.html

