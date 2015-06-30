#s -lR
tar xvzf python.tar.gz --directory $CMSSW_BASE 
ls -lR .
echo "ENV..................................."
env 
echo "VOMS"
voms-proxy-info -all
echo "CMSSW BASE, python path, pwd"
echo $CMSSW_BASE 
echo $PYTHON_PATH
echo $PWD 
cp lib/slc*/* $CMSSW_BASE/lib/slc*
cp lib/slc*/.* $CMSSW_BASE/lib/slc*
echo "AFTER COPY content of $CMSSW_BASE/lib/slc*"
ls -lR  $CMSSW_BASE/lib/slc*

cp -r interface/* $CMSSW_BASE/interface/
echo "AFTER COPY content of $CMSSW_BASE/interface"
ls -lR  $CMSSW_BASE/interface/

cp -r src/* $CMSSW_BASE/src/
echo "AFTER COPY content of $CMSSW_BASE/src"
ls -lR  $CMSSW_BASE/src/

PROXYFILE=`grep "BEGIN CERTIFICATE" * | perl -pe 's/:.*//'  | grep -v heppy | tail -n 1`
export X509_USER_PROXY=$PWD/$PROXYFILE
echo Found Proxy in: $X509_USER_PROXY
MD5SUM=`cat python.tar.gz heppy_config.py | md5sum | awk '{print $1}'`

cat <<EOF > fakeprov.txt
Processing History:
 HEPPY '' '"CMSSW_X_y_Z"' [1]  ($MD5SUM)
EOF

cat <<EOF > $CMSSW_BASE/bin/$SCRAM_ARCH/edmProvDump
#!/bin/sh
cat fakeprov.txt
EOF

chmod +x $CMSSW_BASE/bin/$SCRAM_ARCH/edmProvDump

echo "Which edmProvDump"
which edmProvDump
edmProvDump

#echo PLUGIN DUMP 
#edmPluginDump -f  | grep -A2 -B2 HTTTopJetProducer
# Move JEC files into flace
mkdir jec
mv PHYS14_V4_MC_L1FastJet_AK4PFchs.txt jec/
mv PHYS14_V4_MC_L2Relative_AK4PFchs.txt jec/
mv PHYS14_V4_MC_L3Absolute_AK4PFchs.txt jec/
mkdir csv
mv csv*root csv/
mv Uncertainty_FAKE.txt jec/
python heppy_crab_script.py $1
echo "======================== CMSRUN LOG ============================"
tail -n 200 Output/cmsRun.log 
