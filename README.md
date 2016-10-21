// cmssw
cmsrel CMSSW_8_0_19

cd CMSSW_8_0_19/src

cmsenv

git cms-merge-topic -u perrozzi:rivet_hepmc

// not needed anymore
// rivet-buildplugin GeneratorInterface/RivetInterface/src/HiggsTemplateCrossSections.cc; 

scram b -j 8; 

cd GeneratorInterface/RivetInterface/test 

cmsRun testRivet.py
