#ifndef JECFWLITE
#define JECFWLITE

#include <iostream>
#include <string>
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"


class JECFWLite
{
public:

   JECFWLite(std::string base)
   {
        std::string prefix = base + "/Summer12V3MC";
        parMC.push_back( JetCorrectorParameters((prefix+"_L1FastJet_AK5PFchs.txt").c_str()));
        parMC.push_back( JetCorrectorParameters((prefix+"_L2Relative_AK5PFchs.txt").c_str()));
        parMC.push_back( JetCorrectorParameters((prefix+"_L3Absolute_AK5PFchs.txt").c_str()));
        jetCorrectorMC= new FactorizedJetCorrector(parMC);

/*        prefix = base + "/ReferenceMC";
        parMCRef.push_back( JetCorrectorParameters((prefix+"_L1FastJet_AK5PF.txt").c_str()));
        parMCRef.push_back( JetCorrectorParameters((prefix+"_L2Relative_AK5PF.txt").c_str()));
        parMCRef.push_back( JetCorrectorParameters((prefix+"_L3Absolute_AK5PF.txt").c_str()));
        jetCorrectorMCRef= new FactorizedJetCorrector(parMCRef);
*/

        prefix = base + "/ReferenceMC";
        parMCRef.push_back( JetCorrectorParameters((prefix+"_L1FastJet_AK5PFchs.txt").c_str()));
        parMCRef.push_back( JetCorrectorParameters((prefix+"_L2Relative_AK5PFchs.txt").c_str()));
        parMCRef.push_back( JetCorrectorParameters((prefix+"_L3Absolute_AK5PFchs.txt").c_str()));
        jetCorrectorMCRef= new FactorizedJetCorrector(parMCRef);


        prefix = base + "/Summer12V3DATA";
        parData.push_back( JetCorrectorParameters((prefix+"_L1FastJet_AK5PFchs.txt").c_str()));
        parData.push_back( JetCorrectorParameters((prefix+"_L2Relative_AK5PFchs.txt").c_str()));
        parData.push_back( JetCorrectorParameters((prefix+"_L3Absolute_AK5PFchs.txt").c_str()));
        parData.push_back( JetCorrectorParameters((prefix+"_L2L3Residual_AK5PFchs.txt").c_str()));
        jetCorrectorData= new FactorizedJetCorrector(parData);

        prefix = base + "/Reference";
        parDataRef.push_back( JetCorrectorParameters((prefix+"_L1FastJet_AK5PFchs.txt").c_str()));
        parDataRef.push_back( JetCorrectorParameters((prefix+"_L2Relative_AK5PFchs.txt").c_str()));
        parDataRef.push_back( JetCorrectorParameters((prefix+"_L3Absolute_AK5PFchs.txt").c_str()));
        parDataRef.push_back( JetCorrectorParameters((prefix+"_L2L3Residual_AK5PFchs.txt").c_str()));
        jetCorrectorDataRef= new FactorizedJetCorrector(parDataRef);

        
   }

   VHbbEvent::SimpleJet correct(const VHbbEvent::SimpleJet &j, float rho, bool isMC, bool checkRef = false)
   {
     VHbbEvent::SimpleJet c=j;
     FactorizedJetCorrector * corr=0;
     if(checkRef && isMC) corr=jetCorrectorMCRef;
     if(!checkRef && isMC) corr=jetCorrectorMC;
     if(checkRef && !isMC) corr=jetCorrectorDataRef;
     if(!checkRef && !isMC) corr=jetCorrectorData;
     corr->setJetEta(j.p4.Eta());
     corr->setJetPt(j.ptRaw);
     corr->setJetA(j.jetArea);
     corr->setRho(rho); 
     c.p4 = j.ptRaw*corr->getCorrection();
     corr->setJetEta(j.p4.Eta());
     corr->setJetPt(j.ptRaw);
     corr->setJetA(j.jetArea);
     corr->setRho(rho);

     if(checkRef)
       {
          if(fabs(c.p4.Pt()-j.p4.Pt()) > 0.01 * (c.p4.Pt()+j.p4.Pt())  )
           {
		std::cout << "ERROR CORRECTIONS ARE NOT CLOSING: " << c.p4.Pt() << " vs " <<  j.p4.Pt() << " raw "  << j.ptRaw << " new corr " << corr->getCorrection() << " old " <<  j.p4.Pt()/j.ptRaw <<  std::endl;
           }
 
       }
     return c;
   }   

   std::vector<JetCorrectorParameters> parMC;
   std::vector<JetCorrectorParameters> parMCRef;
   std::vector<JetCorrectorParameters> parData;
   std::vector<JetCorrectorParameters> parDataRef;
   FactorizedJetCorrector * jetCorrectorMC;
   FactorizedJetCorrector * jetCorrectorMCRef;
   FactorizedJetCorrector * jetCorrectorData;
   FactorizedJetCorrector * jetCorrectorDataRef;

};

#endif
