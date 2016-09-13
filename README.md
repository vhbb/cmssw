# cmssw
cmsrel CMSSW_8_0_19

cd CMSSW_8_0_19/src

cmsenv

git cms-addpkg GeneratorInterface/RivetInterface

rivet-buildplugin GeneratorInterface/RivetInterface/src/HiggsTemplateCrossSections.cc; 

scram b -j 8; 

cmsRun testRivet.py
