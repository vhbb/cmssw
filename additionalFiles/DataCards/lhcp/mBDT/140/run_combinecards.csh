# Datacards
cp ../../Zll/BDT/vhbb_DC_BDT_M140_Zll_8TeV.txt          vhbb_Zll_8TeV.txt
cp ../../Znn/res_20130425/mBDT/140/vhbb_Znn_J12_combo_8TeV.txt       vhbb_Znn_8TeV.txt
cp ../../Wln/LHCp_preApproval6/3bins/VH_MultiBDT/combined_vhbb_DC_BDT_H140_8TeV.txt vhbb_Wln_8TeV.txt
cp ../../Wtn/LHCpTauLimits_preApproval_3/VH_BDT/Wtn_BDT_newBinning_140.txt          vhbb_Wtn_8TeV.txt

# Workspaces
cp ../../Zll/BDT/vhbb_WS_BDT_M140_Z*Pt_8TeV.root .
cp ../../Znn/res_20130425/mBDT/140/vhbb_Znn_J12_8TeV.root .
cp ../../Wln/LHCp_preApproval6/3bins/VH_MultiBDT/vhbb_WS_BDT_H140_W*Pt_8TeV.root .
cp ../../Wtn/LHCpTauLimits_preApproval_3/VH_BDT/Wtn_BDT_newBinning_140.root .

# Datacards, Workspaces 7TeV
#cp /afs/cern.ch/user/j/jiafulow/public/vhbb_LHCP_20140422/combo7TeV/BDT/140/vhbb_*txt .
#cp /afs/cern.ch/user/j/jiafulow/public/vhbb_LHCP_20140422/combo7TeV/BDT/140/vhbb_*root .

# Combine cards
combineCards.py Zll=vhbb_Zll_8TeV.txt Znn=vhbb_Znn_8TeV.txt >! vhbb_ZH_8TeV.txt
combineCards.py Wln=vhbb_Wln_8TeV.txt Wtn=vhbb_Wtn_8TeV.txt >! vhbb_WH_8TeV.txt
combineCards.py ZH=vhbb_ZH_8TeV.txt WH=vhbb_WH_8TeV.txt >! vhbb_VH_8TeV.txt
#combineCards.py VH_8TeV=vhbb_VH_8TeV.txt VH_7TeV=vhbb_VH_7TeV.txt >! vhbb_VH_7p8TeV.txt
