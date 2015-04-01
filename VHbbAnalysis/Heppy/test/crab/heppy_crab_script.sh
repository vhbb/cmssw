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



python heppy_crab_script.py $1
