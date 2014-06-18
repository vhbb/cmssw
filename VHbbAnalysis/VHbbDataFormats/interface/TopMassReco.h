#ifndef TOPMASSRECO_H
#define TOPMASSRECO_H

#include <cmath>
#include <vector>

#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"

struct TopHypo {
  TLorentzVector p4;
  TLorentzVector p4W;
};

namespace TopMassReco {
      
       TopHypo topMass(const TLorentzVector & lepton, const TLorentzVector & bJet,const TLorentzVector &  met){
            TopHypo top;
            const double mW2 = 6464.16; //80.4**2;
            double pzNu=0;

            //Did not recalculate any formulas, but seems to check the neutrino solutions to get met z component
        
            const double a = mW2 - lepton.M2() + 2.0*lepton.Px()*met.Px()+2.0*lepton.Py()*met.Py();
            const double A = 4.0*(lepton.E()*lepton.E() - lepton.Pz()*lepton.Pz());
            const double B = -4.0*a*lepton.Pz();
            const double C = 4.0*lepton.E()*lepton.E()*(met.Px()*met.Px()+met.Py()*met.Py()) - a*a;
            const double tmpRoot = B*B - 4.0*A*C;
            if (tmpRoot<0) pzNu = - B/(2*A);
            else{
                const double tmpSol1 = (-B + TMath::Sqrt(tmpRoot))/(2.0*A);
                const double tmpSol2 = (-B - TMath::Sqrt(tmpRoot))/(2.0*A);
                if (fabs(tmpSol2-lepton.Pz()) < fabs(tmpSol1-lepton.Pz())){ 
                    pzNu = tmpSol2;
                } 
                else pzNu = tmpSol1;
            }
            const double nuE = TMath::Sqrt(met.Px()*met.Px()+met.Py()*met.Py()+pzNu*pzNu);

            //new met vector
            TLorentzVector metHyp = TLorentzVector(met.Px(),met.Py(),pzNu,nuE); 

            //W = lepton+met
            top.p4W = (lepton+metHyp);
            //Top = lepton+jet+met
            top.p4 = (lepton+bJet+metHyp);

            return top;
        }
          
}

#endif
