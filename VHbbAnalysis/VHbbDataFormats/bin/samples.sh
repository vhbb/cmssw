#!/bin/sh
for i in  `dbs  --search --url "http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet" --query "find dataset where dataset like *HBB_EDM* and dataset.status like VALID*" | grep USER` ; do 

#/ZZ_TuneZ2_7TeV_pythia6_tauola/tboccali-HBB_EDMNtupleV1_ProcV2-424d8284eb2215daef76cdae9e8db75a/USER

#Tbar_TuneZ2_s-channel_7TeV-powheg-tauola_HBB_EDMNtupleV1_ProcV2
echo $i
NAME=`echo $i | perl -pe 's#/##; s#/tboccali-#_#; s#-[0-9a-z]*/USER##'` 

echo $NAME
#ls ../samples3/$NAME.txt
dbs  --search --url "http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet" --query "find file where dataset = "$i" and dataset.status like VALID*"  | grep root | perl -pe 's#/store#dcap://cmsdcache/pnfs/pi.infn.it/data/cms/store#' > $NAME.txt 


done
