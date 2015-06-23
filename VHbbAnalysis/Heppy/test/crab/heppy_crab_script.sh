ls -lR
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

# Update library path
# Needed so recompiled modules are found 
# the architecture is now (CMSSW 743) gcc491
export LD_LIBRARY_PATH=./lib/slc6_amd64_gcc491:$LD_LIBRARY_PATH 
export ROOT_INCLUDE_PATH=.:./src:$ROOT_INCLUDE_PATH

# Move JEC files into flace
mkdir jec
mv PHYS14_V4_MC_L1FastJet_AK4PFchs.txt jec/
mv PHYS14_V4_MC_L2Relative_AK4PFchs.txt jec/
mv PHYS14_V4_MC_L3Absolute_AK4PFchs.txt jec/
mv Uncertainty_FAKE.txt jec/

mkdir csv
mv csv_rwt_hf_IT_FlatSF.root csv/
mv csv_rwt_lf_IT_FlatSF.root csv/

# Verify plugins
echo "------ Begin: edmPluginDump | grep -i CandidateBoostedDoubleSecondaryVertexESProducer"
edmPluginDump | grep -i CandidateBoostedDoubleSecondaryVertexESProduce
echo "------ End edmPluginDump | grep -i CandidateBoostedDoubleSecondaryVertexESProducer"

python heppy_crab_script.py $1
