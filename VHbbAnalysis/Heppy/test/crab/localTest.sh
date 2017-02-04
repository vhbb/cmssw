#!/bin/sh
FILES=`echo "print ' '.join(config.JobType.inputFiles)" | python -i heppy_crab_config.py `
DIR=/tmp/localtest$2
mkdir $DIR

echo -e "process.source.fileNames = ['$1']\nprint process.dumpPython()" | python -i heppy_crab_fake_pset.py > $DIR/PSet.py
for i in $FILES ; do
cp $i $DIR
done
cp heppy_crab_script.sh $DIR
cd $DIR
scramv1 project CMSSW  $CMSSW_VERSION
#cd  $CMSSW_VERSION
#eval `scramv1 runtime -sh`
#cd -
#MSSW_BASE=$DIR 
./heppy_crab_script.sh $1
