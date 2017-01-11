#!/bin/sh
FILES=`echo "print ' '.join(config.JobType.inputFiles)" | python -i heppy_crab_config.py `
mkdir localtest$2
echo -e "process.source.fileNames = ['$1']\nprint process.dumpPython()" | python -i heppy_crab_fake_pset.py > localtest$2/PSet.py
for i in $FILES ; do
cp $i localtest$2
done
cd localtest$2
../heppy_crab_script.sh $1
