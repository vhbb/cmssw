#ifndef TOPMASSRECO_H
#define TOPMASSRECO_H

#include <cmath>
#include <vector>

#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"

struct topHypo {
  TLorentzVector p4;
  TLorentzVector p4W;
} topQuark;

struct topMassReco {
      
        topHypo operator()(const VHbbCandidate vhCand) const {
            //Look for additional jet that is most b-tagged
            topQuark.p4 = TLorentzVector(0,0,0,0);
            topQuark.p4W = TLorentzVector(0,0,0,0);
            int topJet=-99;
            double minBtag=-9999.;
            for(unsigned int j=0; j < vhCand.additionalJets.size(); j++ ){
                if (vhCand.additionalJets[j].csv > minBtag) topJet = j;
            }
            if(topJet < 0) return topQuark; //If no additional jet, do no computation
            const TLorentzVector bJet = vhCand.additionalJets[topJet].p4;
            const TLorentzVector met = vhCand.V.mets.at(0).p4;
            // just support semiletonic top
            if( vhCand.candidateType == VHbbCandidate::Zee || vhCand.candidateType == VHbbCandidate::Zmumu || vhCand.candidateType == VHbbCandidate::Znn ) return topQuark;
            else if( vhCand.candidateType == VHbbCandidate::Wmun ){
                const double leptonMass = 0.105;
                const TLorentzVector lepton = vhCand.V.muons[0].p4;
                return calcMass(lepton,bJet,met,leptonMass);
            }
            else if( vhCand.candidateType == VHbbCandidate::Wen ){
                const double leptonMass = 0.005;
                const TLorentzVector lepton = vhCand.V.electrons[0].p4;
                return calcMass(lepton,bJet,met,leptonMass);
            }
            else return topQuark;
        }

        static topHypo calcMass(const TLorentzVector lepton, const TLorentzVector bJet,const TLorentzVector met, const double leptonMass){
            const double mW = 80.4;
            double pzNu=0;
        
            //Did not recalculate any formulas, but seems to check the neutrino solutions to get met z component
            const double a = mW*mW - leptonMass*leptonMass + 2.0*lepton.Px()*met.Px()+2.0*lepton.Py()*met.Py();
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
            topQuark.p4W = (lepton+metHyp);
            //Top = lepton+jet+met
            topQuark.p4 = (lepton+bJet+metHyp);

            return topQuark;
        }
          
};

#endif
