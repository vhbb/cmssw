#include "VHbbAnalysis/VHbbDataFormats/interface/Cuts200X.h"

/*
class SignalPreSelectionWen : public Cut {
class SignalPreSelectionWmun : public Cut {
class SignalPreSelectionZee : public Cut {
class SignalPreSelectionZmumu : public Cut {
class SignalPreSelectionZnn : public Cut {
class HPtCut : public PCut
class VPtCut : public PCut
class DoubleBTagCut : public PCut
class SingleBTagCut : public PCut
class VHDeltaPhiCut : public PCut
class AdditionalJetsCut : public PCut
class AdditionalLeptonsCut : public PCut
class METCut : public PCut
class METSignificanceCut : public PCut
class JetMETDeltaPhiCut : public PCut
class DiJetMassMinCut : public PCut
class DiJetMassMaxCut : public PCut

Variable W(`n)H Z(``)H Z(nn)H
pT(b1) > 30 > 20 > 80
pT(b2) > 30 > 20 > 30
pT(jj) > 165 > 100 > 160
pT(V) > 160 > 100 –
CSV1 CSVT CSVT CSVT
CSV2 > 0.5 > 0.5 > 0.5
Df(V,H) > 2.95 > 2.90 > 2.90
Naj = 0 < 2 –
Nal = 0 – = 0
pfMET > 30(W(en)H) – > 160
pfMETsig – – > 5
Df(pfMET, J) – – > 1.5
M(jj)(110) 95–125 90–120 95–125
M(jj)(115) 100–130 95–125 100–130
M(jj)(120) 105–135 100–130 105–135
M(jj)(125) 110–140 105–135 110–140
M(jj)(130) 115–145 110–140 115–145
M(jj)(135) 120–150 115–145 120–150
*/
CutSet buildSignalSelectionZee(float mass)
{
  CutSet result;
  result.add(new SignalPreSelectionZee);
  result.add(new HPtCut(100));
  result.add(new VPtCut(100));
  result.add(new DoubleBTagCut(0.5));
  result.add(new SingleBTagCut(CSVT));
  result.add(new AdditionalJetsCut(2)); // < 2 
//  result.add(new AdditionalLeptonsCut(1));  // < 1
//  result.add(new METCut(30));
  result.add(new DiJetMassMinCut(mass-15. -5.));
  result.add(new DiJetMassMaxCut(mass+15. -5.));
 return result;
}

CutSet buildSignalSelectionZmumu(float mass)
{
  CutSet result;
  result.add(new SignalPreSelectionZee);
  result.add(new HPtCut(100));
  result.add(new VPtCut(100));
  result.add(new DoubleBTagCut(0.5));
  result.add(new SingleBTagCut(CSVT));
  result.add(new AdditionalJetsCut(2)); // < 2 
//  result.add(new AdditionalLeptonsCut(1));  // < 1
//  result.add(new METCut(30));
  result.add(new DiJetMassMinCut(mass-15. -5.));
  result.add(new DiJetMassMaxCut(mass+15. -5.));
 return result;
}


CutSet buildSignalSelectionWen(float mass)
{
  CutSet result;
  result.add(new SignalPreSelectionWen); 
  result.add(new HPtCut(165)); 
  result.add(new VPtCut(160)); 
  result.add(new DoubleBTagCut(0.5)); 
  result.add(new SingleBTagCut(CSVT)); 
  result.add(new AdditionalJetsCut(1)); // < 1 
  result.add(new AdditionalLeptonsCut(1));  // < 1
  result.add(new METCut(30));  
  result.add(new DiJetMassMinCut(mass-15.));  
  result.add(new DiJetMassMaxCut(mass+15.));  
 return result;
}

CutSet buildSignalSelectionWmun(float mass)
{
  CutSet result;
  result.add(new SignalPreSelectionWmun);
  result.add(new HPtCut(165));
  result.add(new VPtCut(160));
  result.add(new DoubleBTagCut(0.5));
  result.add(new SingleBTagCut(CSVT));
  result.add(new AdditionalJetsCut(1)); // < 1 
  result.add(new AdditionalLeptonsCut(1));  // < 1
//  result.add(new METCut(30));
  result.add(new DiJetMassMinCut(mass-15.));
  result.add(new DiJetMassMaxCut(mass+15.));
 return result;
}




