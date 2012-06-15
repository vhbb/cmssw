#include "VHbbAnalysis/VHbbDataFormats/interface/HbbCandidateFinderAlgo.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidateTools.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEventAuxInfo.h"

#include "DataFormats/Math/interface/deltaR.h"
#include <algorithm>

#include <iostream>
#include<cstdlib>

struct CompareJetPt {
  bool operator()( const VHbbEvent::SimpleJet& j1, const  VHbbEvent::SimpleJet& j2 ) const {
    return j1.p4.Pt() > j2.p4.Pt();
  }
};

struct CompareBTag {
  bool operator()(const  VHbbEvent::SimpleJet& j1, const  VHbbEvent::SimpleJet& j2 ) const {
    return j1.csv > j2.csv;
  }
};

VHbbCandidate HbbCandidateFinderAlgo::changeHiggs(bool useHighestPtHiggs , const VHbbCandidate & old)
{
  VHbbCandidateTools selector(verbose_);

  VHbbCandidate temp(old);
  VHbbEvent::SimpleJet j1,j2;
  std::vector<VHbbEvent::SimpleJet> addJets;
  std::vector<VHbbEvent::SimpleJet> jets;

  for(size_t i=0; i < old.H.jets.size();i++) jets.push_back(old.H.jets[i]);
  for(size_t i=0; i < old.additionalJets.size();i++) jets.push_back(old.additionalJets[i]);

  bool foundJets;

  if (useHighestPtHiggs == false){
    foundJets = findDiJets(jets,j1,j2,addJets) ;
  }else{
    foundJets= findDiJetsHighestPt(jets,j1,j2,addJets) ;
  }

  temp.H.jets.clear();
  temp.H.jets.push_back(j1);
  temp.H.jets.push_back(j2);
  temp.H.p4 = (j1).p4+(j2).p4;
  TVector3 higgsBoost;
  higgsBoost = (temp.H.p4).BoostVector();
  temp.H.helicities.clear();
  temp.H.helicities.push_back(selector.getHelicity(j1,higgsBoost));
  temp.H.helicities.push_back(selector.getHelicity(j2,higgsBoost));
  temp.H.deltaTheta = selector.getDeltaTheta(j1,j2);
  temp.additionalJets = addJets;
  return temp;

}



void HbbCandidateFinderAlgo::run (const VHbbEvent* event, std::vector<VHbbCandidate>  & candidates,const VHbbEventAuxInfo & aux){
  //
  // search for leptons
  //
  std::vector<VHbbEvent::MuonInfo> mu;
  std::vector<VHbbEvent::ElectronInfo> ele;
//  std::vector<VHbbEvent::MuonInfo> mu;
  std::vector<unsigned int> muPos;
  findMuons(event->muInfo,mu, muPos,aux);
//  std::vector<VHbbEvent::ElectronInfo> ele;
  std::vector<unsigned int> elePos;
  findElectrons(event->eleInfo,ele, elePos,aux);

  std::vector<VHbbEvent::TauInfo> tau;	 
  std::vector<unsigned int> tauPos;	 
  findTaus(event->tauInfo,tau, tauPos);

//FIXME: cleaning here

  //
  // first find the jets
  //

  VHbbCandidateTools selector(verbose_);
  std::vector<VHbbEvent::SimpleJet> noOverlap;
  for(size_t j=0; j< event->simpleJets2.size(); j++)
  {
     float overlap=false;
     for(size_t i=0; i< mu.size(); i++) {
       if(deltaR(mu[i].p4.Eta(),mu[i].p4.Phi(),event->simpleJets2[j].p4.Eta(),event->simpleJets2[j].p4.Phi()) < 0.5) overlap=true; 
     }
     for(size_t i=0; i< ele.size(); i++) {
       if(deltaR(ele[i].p4.Eta(),ele[i].p4.Phi(),event->simpleJets2[j].p4.Eta(),event->simpleJets2[j].p4.Phi()) < 0.5) overlap=true; 
     }
   
     if(!overlap) noOverlap.push_back(event->simpleJets2[j]);
     else 
     {
    //    std::cout << "jet removed in cleaning" << std::endl;
     }   
  }


  VHbbEvent::SimpleJet j1,j2;
  std::vector<VHbbEvent::SimpleJet> addJets;
  bool foundJets;
  if (useHighestPtHiggs_ == false){
    foundJets = findDiJets(noOverlap,j1,j2,addJets) ;
  }else{
    foundJets= findDiJetsHighestPt(noOverlap,j1,j2,addJets) ;
  }

  if (verbose_){
    std::cout <<" Found Dijets: "<<foundJets<< " Additional: "<<addJets.size()<< std::endl;
  }
  
//  if (foundJets == false) return;


  bool foundHardJets;
  VHbbEvent::HardJet fatj1;
  std::vector<VHbbEvent::SimpleJet> subJetsout;
  std::vector<VHbbEvent::SimpleJet> addJetsFat;
  // foundHardJets= findFatJet(event->hardJets,event->subJets,event->filterJets,fatj1,subJetsout, event->simpleJets2, addJetsFat, mu, ele) ;
  foundHardJets= findFatJet(event->hardJets,event->subJets,event->filterJets,fatj1,subJetsout, event->simpleJets2, addJetsFat, event->muInfo, event->eleInfo) ;

  if (foundJets == false && foundHardJets == false) return;


  
  std::vector<VHbbEvent::METInfo> met;

  //  findMET(event->pfmet, met);    
   findMET(event->pfmetType1corr, met);

  if (verbose_){
    std::cout <<" Electrons: "<< ele.size()<<std::endl;
    std::cout <<" Muons    : "<< mu.size()<<std::endl;
    std::cout <<" Taus    : "<< tau.size()<<std::endl;	
    std::cout <<" MET      : "<< met.size()<<std::endl;
  }
  if (ele.size()<1 && mu.size() < 1 && met.size()<1 && tau.size()<1) return;
  //
  // fill!
  //
  VHbbCandidate temp;
  temp.H.HiggsFlag = foundJets;
  if(foundJets){
  temp.H.jets.push_back(j1);
  temp.H.jets.push_back(j2);
  temp.H.p4 = (j1).p4+(j2).p4;
  TVector3 higgsBoost;
  higgsBoost = (temp.H.p4).BoostVector();
  temp.H.helicities.push_back(selector.getHelicity(j1,higgsBoost));
  temp.H.helicities.push_back(selector.getHelicity(j2,higgsBoost));
  temp.H.deltaTheta = selector.getDeltaTheta(j1,j2);
  }

  temp.FatH.FatHiggsFlag= foundHardJets;
  if(foundHardJets){
  temp.FatH.p4 = fatj1.p4;
  temp.FatH.subjetsSize=subJetsout.size();
  if(subJetsout.size()==2) {temp.FatH.jets.push_back(subJetsout[0]);temp.FatH.jets.push_back(subJetsout[1]);}
  if(subJetsout.size()==3) {temp.FatH.jets.push_back(subJetsout[0]);
  temp.FatH.jets.push_back(subJetsout[1]);temp.FatH.jets.push_back(subJetsout[2]);}
  }

  std::vector<VHbbEvent::TauInfo> tauNoCandidateJetOverlap;	 
  std::vector<unsigned int> tauPosNoCandidateJetOverlap;	 
  removeTauOverlapWithJets(tau,temp.H.jets,tauNoCandidateJetOverlap,tauPos,tauPosNoCandidateJetOverlap);

  temp.additionalJets = addJets;
  temp.additionalJetsFat = addJetsFat;
  temp.V.mets = met;
  temp.V.muons = mu;
  temp.V.electrons = ele;
  temp.V.taus = tauNoCandidateJetOverlap;	
  
  //
  // now see which kind of andidate this can be
  // 

  VHbbCandidate result;
  bool ok = false;
  //
  // first: hZmumu
  //
  
  if (verbose_){
    std::cout <<" START SELECTION "<<std::endl;
  }
  //
  //
  // change must allow a candidate to be both Zee and Zmumu

  // Zmumu
  result = selector.getHZmumuCandidate(temp,ok,muPos);
  if ( ok == true ){
    result.setCandidateType(VHbbCandidate::Zmumu);
    candidates.push_back(result);
  }
  //      HZee
  result = selector. getHZeeCandidate(temp,ok,elePos);
  if ( ok == true ){
    result.setCandidateType(VHbbCandidate::Zee);
    candidates.push_back(result);
  }
  //HWmunu
  result = selector. getHWmunCandidate(temp,ok,muPos);
  if ( ok == true ){
    result.setCandidateType(VHbbCandidate::Wmun);
    candidates.push_back(result);
  }
  // HWenu
  result = selector. getHWenCandidate(temp,ok,elePos);
  if ( ok == true ){
    result.setCandidateType(VHbbCandidate::Wen);
    candidates.push_back(result);
  }

  // New tau categorizations currently commented out	 
   /*	 
   result = selector.getHWtaunCandidate(temp,ok,tauPosNoCandidateJetOverlap);	 
   if ( ok == true ){	 
     if (verbose_) std::cout << "We have a taun candidate" << std::endl;	 
     result.setCandidateType(VHbbCandidate::Wtaun);	 
     candidates.push_back(result);	 
   }	 
   result = selector.getHZtaumuCandidate(temp,ok,muPos,tauPosNoCandidateJetOverlap);	 
   if ( ok == true ){	 
     if (verbose_) std::cout << "We have a HZtaumu candidate" << std::endl;	 
     result.setCandidateType(VHbbCandidate::Ztaumu);	 
     candidates.push_back(result);	 
   }	 
   */

  if (candidates.size()!=0 ) return;

  // HZnn - look at it only if nothing found up to now
  result = selector. getHZnnCandidate(temp,ok);
  if ( ok == true ){
    result.setCandidateType(VHbbCandidate::Znn);
    candidates.push_back(result);
  }
  return;
}

void HbbCandidateFinderAlgo::findMET(const VHbbEvent::METInfo & met, std::vector<VHbbEvent::METInfo>& out){
  //

  //  just preselection: met significance > 2 
  // removed on request by A. Rizzi - no filter at all
  //    if (met.metSig >2 ) out.push_back(met);
  out.push_back(met);
    if (verbose_){
      std::cout <<" CandidateFinder: Input MET = "<<met.metSig<<" Output MET = "<<out.size()<<std::endl;
    }
    
}
bool HbbCandidateFinderAlgo::jetID(const VHbbEvent::SimpleJet & j)
{
    if(j.neutralHadronEFraction > 0.99) return false;
    if(j.neutralEmEFraction > 0.99) return false;
    if(fabs(j.p4.Eta())<2.4 && j.chargedEmEFraction > 0.99) return false;
    if(fabs(j.p4.Eta())<2.4 && j.chargedHadronEFraction == 0) return false;
    if(fabs(j.p4.Eta())<2.4 && j.ntracks == 0) return false;
    if(j.nConstituents <= 1) return false;
return true;
}


bool HbbCandidateFinderAlgo::findDiJets (const std::vector<VHbbEvent::SimpleJet>& jetsin, VHbbEvent::SimpleJet& j1, VHbbEvent::SimpleJet& j2,std::vector<VHbbEvent::SimpleJet>& addJets){
  


if (verbose_){
   std::cout <<" CandidateFinder: Input Jets = "<<jetsin.size()<<std::endl;
 }


 CompareBTag  bTagComparator;
 CompareJetPt ptComparator;
 float etaThr = 2.5;


 if (jetsin.size()<2) return false;

 std::vector<VHbbEvent::SimpleJet> jets = jetsin;
 
 std::sort(jets.begin(), jets.end(), bTagComparator);

 //
 // now I need at least 2 with pt > threshold
 //
 unsigned int index1=999999, index2=999999;
 for (unsigned int i =0; i< jets.size(); ++i){
   if (jets[i].p4.Pt()> jetPtThreshold && fabs(jets[i].p4.Eta()) < etaThr /*&& jetID(jets[i])*/){
     if (index1 == 999999) {
       index1=i;
     }else if (index2 == 999999){
       index2=i;
       break;
     }
   }
}
 if (index1==999999 || index2== 999999) return false;

if (jets[index1].p4.Pt()<(jets[index2].p4.Pt())){
  std::swap (index1,index2);
 }
 j1 = jets[index1];
 j2 = jets[index2];
 
 //
 // additional jets
 //
 if (jets.size()>2){
   for (unsigned int i=0 ; i< jets.size(); ++i){
     if (i != index1 && i != index2 )
       addJets.push_back(jets[i]);
   }
 }
 
 if (verbose_){
   std::cout <<" CandidateFinder: Output Jets = "<<2<<" Additional = "<<addJets.size()<<std::endl;
 }
 
 
 std::sort(addJets.begin(), addJets.end(), ptComparator);
 return true;
 
 



}


bool HbbCandidateFinderAlgo::findDiJetsHighestPt (const std::vector<VHbbEvent::SimpleJet>& jetsin, VHbbEvent::SimpleJet& j1, VHbbEvent::SimpleJet& j2,std::vector<VHbbEvent::SimpleJet>& addJets){
  

  if (verbose_){
    std::cout <<" CandidateFinder: Input Jets = "<<jetsin.size()<<std::endl;
  }
  if (jetsin.size()<2) return false;
  
  float etaThr = 2.5;
  std::vector<VHbbEvent::SimpleJet> jets = jetsin;
  //loop over the dijets and save the one with highest Pt
  
  CompareJetPt ptComparator;
  std::sort(jets.begin(), jets.end(), ptComparator);
  //
  // so if i<j, pt(i)>pt(j)
  //
  
  float highestPt = -1000;
  unsigned  int highesti=999999,highestj=999999;
  for (unsigned int i =0; i< jets.size()-1; ++i){
    for (unsigned int j =i+1; j< jets.size(); ++j){
      float pt = (jets[i].p4+jets[j].p4).Pt();
      if (pt>highestPt && jets[j].p4.Pt()> jetPtThreshold && jets[i].p4.Pt()> jetPtThreshold && fabs(jets[i].p4.Eta()) < etaThr &&  fabs(jets[j].p4.Eta()) < etaThr && jetID(jets[i]) && jetID(jets[j]) ){
	highestPt = pt;
	highesti=i;
	highestj=j;
      }
    }
  }
  
  if (highesti == 999999 || highestj == 999999) return false;
  j1 = jets[highesti]; 
  j2 = jets[highestj]; 
  
  for (unsigned int i=0; i<jets.size(); ++i){
    if (i!= highesti && i!= highestj)
      addJets.push_back(jets[i]);
  }
  
  if (verbose_){
    std::cout <<" CandidateFinder: Output Jets = "<<2<<" Additional = "<<addJets.size()<<std::endl;
  }
  std::sort(addJets.begin(), addJets.end(), ptComparator);
  return true;
}


bool HbbCandidateFinderAlgo::findFatJet (const std::vector<VHbbEvent::HardJet>& jetsin, const std::vector<VHbbEvent::SimpleJet>& subjetsin, const std::vector<VHbbEvent::SimpleJet>& filterjetsin, VHbbEvent::HardJet& fatj1,std::vector<VHbbEvent::SimpleJet>& subJetsout, const std::vector<VHbbEvent::SimpleJet>& ak5jetsin, std::vector<VHbbEvent::SimpleJet>& addJetsFat, const std::vector<VHbbEvent::MuonInfo>& muons, const std::vector<VHbbEvent::ElectronInfo>& electrons){

  if (verbose_){
    std::cout <<" CandidateFinder: Input Jets = "<<jetsin.size()<<std::endl;
  }
  if (jetsin.size()<1) return false;
 
  float etaThr = 2.5;
  std::vector<VHbbEvent::HardJet> hardjets = jetsin;
  std::vector<VHbbEvent::SimpleJet> subjets = subjetsin;
  std::vector<VHbbEvent::SimpleJet> filterjets = filterjetsin;

  fatj1=hardjets[0];  // to avoid warning if subjet below fail selection


/*  TMatrixD *pointerEta = new TMatrixD(90,80);
  TMatrixD* pointerPhi = new TMatrixD(90,80);
  for (unsigned int i =0; i< hardjets.size(); i++){
   for (unsigned int j=0; j< i.constituents; j++){
   TMatrixDRow(*pointerEta,i)(j)=i.etaSub[j];
   TMatrixDRow(*pointerPhi,i)(j)=i.phiSub[j];
   }
  }
*/

// debug
/*    std::cout << "hardjet size: " << hardjets.size() << "\n";
    std::cout << "subjet size: " << subjets.size() << "\n";
    std::cout << "filterjet size: " << filterjets.size() << "\n";
       for(unsigned int kk=0;kk<subjets.size();kk++){
    std::cout << "subjet  pt: " << subjets[kk].p4.Pt() << " eta,phi " <<
    subjets[kk].p4.Eta()  << " , " << subjets[kk].p4.Phi()  << "\n";
      }
       for(unsigned int kk=0;kk<filterjets.size();kk++){
    std::cout << "filterjet  pt: " << filterjets[kk].p4.Pt() << " eta,phi " <<
    filterjets[kk].p4.Eta()  << " , " << filterjets[kk].p4.Phi()  << "\n";
      }*/
// debug

  double minBtag1=-9999.;
  for (unsigned int i =0; i< hardjets.size(); i++){
     int subJetIn[300]; 
     for(int k=0;k<300;k++)subJetIn[k]=-99;
     int subJetIn1st[2]; 
     for(int k=0;k<2;k++)subJetIn1st[k]=-99;
//    TMatrixDRow* roweta=new TMatrixDRow(*pointerEta,i);	 
//    TMatrixDRow* rowphi=new TMatrixDRow(*pointerPhi,i);	 
 

// debug
//    std::cout << "HardJet pt: " << hardjets[i].p4.Pt() << " # daughters " << 
//    hardjets[i].constituents << "\n";
// debug

   if(hardjets[i].constituents<4) continue;

// first get ja and jb  from first decomposition
     int In1st=0; 
     for(int j=0;j<2;j++){
// debug
//    std::cout << "hardJet constituent pt: " << hardjets[i].subFourMomentum[j].Pt() << " eta,phi " <<
//    hardjets[i].subFourMomentum[j].Eta()  << " , " << hardjets[i].subFourMomentum[j].Phi()  << "\n";
// debug
       for(unsigned int kk=0;kk<subjets.size();kk++){
//       if(subjets.eta[kk]==roweta->GetPtr()[j] && subjets.phi[kk]==rowphi->GetPtr()[j])
//       subJetIn1st[In1st]=kk;} 
//        if(subjets[kk].p4.Eta()==hardjets[i].etaSub[j] && subjets[kk].p4.Phi()==hardjets[i].phiSub[j])
//        subJetIn1st[In1st]=kk;} 
//       if(subjets[kk].p4.Eta()==hardjets[i].subFourMomentum[j].Eta() && subjets[kk].p4.Phi()==hardjets[i].subFourMomentum[j].Phi())
//        subJetIn1st[In1st]=kk;}
       double deltaR_1=deltaR(subjets[kk].p4.Eta(),subjets[kk].p4.Phi(),hardjets[i].subFourMomentum[j].Eta(),hardjets[i].subFourMomentum[j].Phi()); if(deltaR_1<0.01) subJetIn1st[In1st]=kk;} 
       In1st++;
     }

// then get all subjets from decomposition with
     for(int j=2;j<hardjets[i].constituents;j++){
// debug
//    std::cout << "hardJet constituent pt: " << hardjets[i].subFourMomentum[j].Pt() << " eta,phi " <<
//    hardjets[i].subFourMomentum[j].Eta()  << " , " << hardjets[i].subFourMomentum[j].Phi()  << "\n";
// debug
//    cout << "n: " << i << "," << j << " hardjetsub eta: " << roweta->GetPtr()[j] << " phi: " << 
//    rowphi->GetPtr()[j] << "\n";
       for(unsigned int kk=0;kk<filterjets.size();kk++){
//       if(subjets.eta[kk]==roweta->GetPtr()[j] && subjets.phi[kk]==rowphi->GetPtr()[j])
//       subJetIn[j]=kk;} 
//        if(subjets[kk].p4.Eta()==hardjets[i].etaSub[j] && subjets[kk].p4.Phi()==hardjets[i].phiSub[j])
//        subJetIn[j]=kk;} 
//       if(subjets[kk].p4.Eta()==hardjets[i].subFourMomentum[j].Eta() && subjets[kk].p4.Phi()==hardjets[i].subFourMomentum[j].Phi()) subJetIn[j]=kk;
       double deltaR_1=deltaR(filterjets[kk].p4.Eta(),filterjets[kk].p4.Phi(),hardjets[i].subFourMomentum[j].Eta(),hardjets[i].subFourMomentum[j].Phi()); if(deltaR_1<0.01) subJetIn[j-2]=kk;}
           } 
 

//debug
//  std::cout << "index in subjetTag: " << subJetIn1st[0] << "," << subJetIn1st[1] << "\n";
//  std::cout << "index in subfilterTag: " << subJetIn[0] << "," << subJetIn[1] << "," << subJetIn[2] << "\n";
// debug
    
     if(subJetIn1st[0]==-99 || subJetIn1st[1]==-99) continue;
     if(subJetIn[0]==-99 || subJetIn[1]==-99) continue;

     int nBtag=0;
     for(int j=0;j<2;j++){
     if(subjets[subJetIn1st[j]].csv>0.) nBtag++;}

     int nPt=0;
     for(int j=0;j<2;j++){
////     if(subjets[subJetIn1st[j]].p4.Pt()>30. && fabs(subjets[subJetIn1st[j]].p4.Eta())<etaThr && jetID(subjets[subJetIn1st[j]])) nPt++;}
     if(filterjets[subJetIn[j]].p4.Pt()>20. && fabs(filterjets[subJetIn[j]].p4.Eta())<etaThr && jetID(filterjets[subJetIn[j]])) nPt++;}

//     if(nBtag<2 || nPt<2) continue;
     if(nPt<2) continue;

///// lepton overlap
     int muOverlap=0;
  for (unsigned int it = 0; it < muons.size();  ++it){
    if (
        muons[it]. globChi2<10 &&
        muons[it].nPixelHits>= 1 &&
        muons[it].globNHits != 0 &&
        muons[it].nHits > 10 &&
        //tracker
        (muons[it].cat & 0x1) &&
        //global
        (muons[it].cat & 0x2) &&
        muons[it].nMatches >=2 &&
        muons[it].ipDb<.2 &&
        (muons[it].pfChaIso+muons[it].pfPhoIso+muons[it].pfNeuIso)/muons[it].p4.Pt()<.15  &&
        fabs(muons[it].p4.Eta())<2.4 &&
        muons[it].p4.Pt()>20 ) {
     for(int j=0;j<2;j++){
    if(deltaR(muons[it].p4.Eta(),muons[it].p4.Phi(),filterjets[subJetIn[j]].p4.Eta(),filterjets[subJetIn[j]].p4.Phi())<0.3) muOverlap++;}
  }
  }

    int elecOverlap=0;
  for (unsigned int  it = 0; it< electrons.size(); ++it){
    if (
        // fake
        (fabs(electrons[it].id95 - 7)) < 0.1  &&
        fabs(electrons[it].p4.Eta()) < 2.5 &&
//Remove this workaround as now we have the proper flags
//      !( fabs(electrons[it].p4.Eta()) < 1.57 && fabs(electrons[it].p4.Eta()) > 1.44) &&
        electrons[it].p4.Pt()>15 //  I use the minimum ok for both Z and W
         && (electrons[it].pfChaIso+electrons[it].pfPhoIso+electrons[it].pfNeuIso)/electrons[it].p4.Pt()<.15
        ){
     for(int j=0;j<2;j++){
    if(deltaR(electrons[it].p4.Eta(),electrons[it].p4.Phi(),filterjets[subJetIn[j]].p4.Eta(),filterjets[subJetIn[j]].p4.Phi())<0.3) elecOverlap++;}
    }
  }


   if(muOverlap>0) continue;
   if(elecOverlap>0) continue;

//      if(subjets[subJetIn1st[0]].csv+subjets[subJetIn1st[1]].csv>minBtag1){
//      minBtag1=subjets[subJetIn1st[0]].csv+subjets[subJetIn1st[1]].csv;       
       if((subjets[subJetIn1st[0]].p4+subjets[subJetIn1st[1]].p4).Pt()>minBtag1){
       minBtag1=(subjets[subJetIn1st[0]].p4+subjets[subJetIn1st[1]].p4).Pt();
/*       double filtpt=0;
       if(subJetIn[0]!=-99 && subJetIn[1]!=-99) filtpt=(filterjets[subJetIn[0]].p4+filterjets[subJetIn[1]].p4).Pt(); 
       if(subJetIn[0]!=-99 && subJetIn[1]!=-99 && subJetIn[2]!=-99) filtpt=(filterjets[subJetIn[0]].p4+filterjets[subJetIn[1]].p4+filterjets[subJetIn[2]].p4).Pt();
       if(filtpt>minBtag1){
       minBtag1=filtpt;
*/
//       if(hardjets[i].p4.Pt()>minBtag1){
//       minBtag1=hardjets[i].p4.Pt();  
       fatj1=hardjets[i];
       subJetsout.clear();
       if(subJetIn[0]!=-99) subJetsout.push_back(filterjets[subJetIn[0]]);
       if(subJetIn[1]!=-99) subJetsout.push_back(filterjets[subJetIn[1]]);
       if(subJetIn[2]!=-99) subJetsout.push_back(filterjets[subJetIn[2]]);
       }  

  } // loop hard jet


 //
 // additional jets
 //
   std::vector<VHbbEvent::SimpleJet> ak5jets = ak5jetsin;
 
   addJetsFat.clear(); 
   for (unsigned int i=0 ; i< ak5jets.size(); ++i){
     int overlap=0;
     for (unsigned int j=0 ; j< subJetsout.size(); ++j){ 
      if(subJetsout[j].p4.Pt() <20.) continue; 
      if (deltaR(ak5jets[i].p4.Eta(),ak5jets[i].p4.Phi(),subJetsout[j].p4.Eta(),subJetsout[j].p4.Phi())<0.3) overlap++;
     }   
      if(overlap==0) addJetsFat.push_back(ak5jets[i]);
   }
 

  
  return true;
}

void HbbCandidateFinderAlgo::removeTauOverlapWithJets(const std::vector<VHbbEvent::TauInfo>& taus, const std::vector<VHbbEvent::SimpleJet>& jets, std::vector<VHbbEvent::TauInfo>& out, const std::vector<unsigned int>& oldPositions,std::vector<unsigned int>& positions) {
  for (unsigned int it = 0; it < taus.size();  ++it){
    bool overlap = false;
    for (unsigned int jit = 0 ; jit < jets.size() ; ++jit) {
      if (taus[it].p4.DeltaR(jets[jit].p4) < 0.2) {
	overlap = true;
	if (verbose_) {
	  std::cout << "Found overlap of tau (pt,eta,phi)=(" << taus[it].p4.Pt() <<","<< taus[it].p4.Eta() <<","<< taus[it].p4.Phi() << ")"
		    << "with candidate jet (pt,eta,phi)=(" << jets[jit].p4.Pt() <<","<< jets[jit].p4.Eta() <<","<< jets[jit].p4.Phi() << ")"
		    << std::endl;
	}
      }
    }
    if (!overlap) {
      if (verbose_) {
	std::cout << "No overlap with tau (pt,eta,phi)=(" << taus[it].p4.Pt() <<","<< taus[it].p4.Eta() <<","<< taus[it].p4.Phi() << "); keeping it " << std::endl;
      }
      out.push_back(taus[it]);
	             positions.push_back(oldPositions[it]);
    }
  }
}
 
void HbbCandidateFinderAlgo::findTaus(const std::vector<VHbbEvent::TauInfo>& taus, std::vector<VHbbEvent::TauInfo>& out, std::vector<unsigned int>& positions){
  if (verbose_) std::cout << "[SCZ] Tau size=" << taus.size() << std::endl;
  for (unsigned int it = 0; it < taus.size();  ++it){
    /*
      myPatTau.tauID("decayModeFinding"); // cuts on tau invariant mass, etc
      myPatTau.tauID("byLooseCombinedIsolationDeltaBetaCorr") // isolated
      taus, corrected for PU use DB technique
      myPatTau.tauID("againstMuonTight");  // remove muons faking hadronic taus
      myPatTau.tauID("againstElectronLoose/againstElectronMedium/againstElectronMVA");
      // remove electrons faking hadronic taus, choose based on your fake - e background.
      */
 
    if (verbose_) {
      std::cout << "(pt,decayModeFinding,byLooseCombinedIsolationDeltaBetaCorr,againstMuonTight,againstElectronLoose,againstElectronMedium,againstElectronMVA)=("
		<< taus[it].p4.Pt() << ","
		<< (taus[it].decayModeFinding>0.5) << ","
		<< (taus[it].byLooseCombinedIsolationDeltaBetaCorr>0.5) << ","
		<< (taus[it].againstMuonTight>0.5) << ","
		<< (taus[it].againstElectronLoose>0.5) <<","
		<< (taus[it].againstElectronMedium>0.5) <<","
		<< (taus[it].againstElectronMVA>0.5) << ")" << std::endl;
    }
    
    if (taus[it].decayModeFinding>0.5&&taus[it].byLooseCombinedIsolationDeltaBetaCorr>0.5&&taus[it].againstMuonTight>0.5&&taus[it].againstElectronLoose>0.5&&taus[it].p4.Pt()>20.) {
      out.push_back(taus[it]);
      positions.push_back(it);
    }
  }
  if (verbose_){
    std::cout <<" CandidateFinder: Input Taus = "<<taus.size()<<" Output Taus = "<<out.size()<<std::endl;
  }
}


void HbbCandidateFinderAlgo::findMuons(const std::vector<VHbbEvent::MuonInfo>& muons, std::vector<VHbbEvent::MuonInfo>& out, std::vector<unsigned int>& positions,const  VHbbEventAuxInfo & aux){
  /* Use:
For both W -> mu nu and Z -> mu mu, we adopt the standard VBTF muon selection described in VbtfWmunuBaselineSelection. The explicit cuts are reproduced here:

    We use RECO (pf?) Muons that are both Global and Tracker
    chi2/ndof < 10 for the global muon fit
    The track associated to the muon must have
        >= 1 pixel hits
        >= 10 pixel + strip hits
        >= 1 valid hit in the muon chambers
        >= 2 muon stations
        |dxy| < 0.2
        |eta| < 2.4 
    PF Relative combined isolation (R) is required to be < 0.15
        R = [pfChaIso + pfNeuIso + pfPhoIso] / pT(mu) computed in a cone of radius 0.3 in eta-phi 
    pT(mu) > 20 GeV 
  */
  //  for (std::vector<VHbbEvent::MuonInfo>::const_iterator it = muons.begin(); it!= muons.end(); ++it){

/* New iso:

Alternatively, for analysis using rho correction, Effective Areas are also provided using following prescription:
Correction to be done as PFIsoCorr = PF(PFNoPU) – Max ((PF(Nh+Ph) - ρ’EACombined),0.0)) where ρ’=max(ρ,0.0) and with a 0.5 GeV threshold on neutrals
Rho is neutral rho, defined in full tracker acceptance, with a 0.5 GeV threshold on neutrals. This can be taken, starting from 50X, from the event directly (double_kt6PFJetsCentralNeutral_rho_RECO.obj) For its exact definition see [2].
Values of Effective Areas EACombined are provided in this link, page 9 (see PF Combined Column, DeltaR >0/4 )
[2] kt6PFJetsCentralNeutral = kt6PFJets.clone( src = cms.InputTag("pfAllNeutralHadronsAndPhotons"), Ghost_EtaMax = cms.double(3.1), Rho_EtaMax = cms.double(2.5), inputEtMin = cms.double(0.5) )

Effective area, last column matters for us:
           |                   PF Neutral                  |                  PF  Photons                 |                PF  Combined                 |
     Muon Eta       |    DR < 0.30       |    DR < 0.40      |    DR < 0.30       |    DR < 0.40      |    DR < 0.30       |    DR < 0.40       |
0.0 < |eta| < 1.0  | 0.107 +/- 0.010 | 0.166 +/- 0.013 | 0.274 +/- 0.017 | 0.504 +/- 0.020 | 0.382 +/- 0.034 | 0.674 +/- 0.020 |
1.0 < |eta| < 1.5  | 0.141 +/- 0.016 | 0.259 +/- 0.023 | 0.161 +/- 0.019 | 0.306 +/- 0.027 | 0.317 +/- 0.045 | 0.565 +/- 0.027 |
1.5 < |eta| < 2.0  | 0.159 +/- 0.017 | 0.247 +/- 0.025 | 0.079 +/- 0.015 | 0.198 +/- 0.026 | 0.242 +/- 0.040 | 0.442 +/- 0.026 |
2.0 < |eta| < 2.2  | 0.102 +/- 0.022 | 0.220 +/- 0.036 | 0.168 +/- 0.035 | 0.287 +/- 0.048 | 0.326 +/- 0.076 | 0.515 +/- 0.048 |
2.2 < |eta| < 2.3  | 0.096 +/- 0.030 | 0.340 +/- 0.072 | 0.359 +/- 0.072 | 0.525 +/- 0.074 | 0.462 +/- 0.159 | 0.821 +/- 0.074 |
2.3 < |eta| < 2.4  | 0.104 +/- 0.036 | 0.216 +/- 0.056 | 0.294 +/- 0.064 | 0.488 +/- 0.083 | 0.372 +/- 0.141 | 0.660 +/- 0.083 |

*/

  for (unsigned int it = 0; it < muons.size();  ++it){
float mincor=0.0;
float minrho=0.0;
float rhoN = std::max(aux.puInfo.rhoNeutral,minrho);
float eta=muons[it].p4.Eta();
float area=0.5;
if(fabs(eta)>0.0 && fabs(eta) <= 1.0) {area=0.674;}
if(fabs(eta)>1.0 && fabs(eta) <= 1.5) {area=0.565;}
if(fabs(eta)>1.5 && fabs(eta) <= 2.0) {area=0.442;}
if(fabs(eta)>2.0 && fabs(eta) <= 2.2) {area=0.515;}
if(fabs(eta)>2.2 && fabs(eta) <= 2.3) {area=0.821;}
if(fabs(eta)>2.3 && fabs(eta) <= 2.4) {area=0.660;}
float pfCorrIso = (muons[it].pfChaIso+ std::max(muons[it].pfPhoIso+muons[it].pfNeuIso-rhoN*area,mincor))/muons[it].p4.Pt();
    if (
        muons[it].isPF &&
	muons[it]. globChi2<10 &&
	muons[it].nPixelHits>= 1 &&
        muons[it].globNHits != 0 &&
	muons[it].nValidLayers > 5 &&
//	muons[it].nHits > 10 &&
        //tracker
//	(muons[it].cat & 0x2) && 
	//global
	(muons[it].cat & 0x1) && 
	muons[it].nMatches >=2 &&
	muons[it].ipDb<.2 &&
	muons[it].zPVPt<0.5 &&
	//	(muons[it].hIso+muons[it].eIso+muons[it].tIso)/muons[it].p4.Pt()<.15 &&
        pfCorrIso < 0.12 &&
	fabs(muons[it].p4.Eta())<2.4 &&
	muons[it].p4.Pt()>20 ) {
      out.push_back(muons[it]);
      positions.push_back(it);
  }
  }
    if (verbose_){
      std::cout <<" CandidateFinder: Input Muons = "<<muons.size()<<" Output Muons = "<<out.size()<<std::endl;
    }


}


void HbbCandidateFinderAlgo::findElectrons(const std::vector<VHbbEvent::ElectronInfo>& electrons, std::vector<VHbbEvent::ElectronInfo>& out,  std::vector<unsigned int>& positions,const  VHbbEventAuxInfo & aux){
 /*
We adopt the standard cut-based selection from VBTF described in detail here.

    Z -> ee
        gsf (pf?) electrons
        VBTF WP95
        |eta|<2.5, excluding the gap 1.44 < |eta| < 1.57
        pT(e) > 20 

    W -> e nu
        gsf (pf?) electrons
        VBTF WP80
        |eta|<2.5, excluding the gap 1.44 < |eta| < 1.57
        pT(e) > 30 
  */
  //  for (std::vector<VHbbEvent::ElectronInfo>::const_iterator it = electrons.begin(); it!= electrons.end(); ++it){
/*
      isocorr = chargediso + max(PFIso(&gamma;) - rho * Aeff(&gamma), 0.) +  max(PFIso(NH) - rho * Aeff(NH), 0.) 
0<abs(eta)<1.0   Aeff(NH) = 0.024 +/- 0.001      Aeff(γ) = 0.081 +/- 0.001       Aeff(γ+NH) = 0.10 +/- 0.002
1.0<abs(eta)<1.479       Aeff(NH) = 0.037 +/- 0.001      Aeff(γ) = 0.084 +/- 0.003       Aeff(γ+NH) = 0.12 +/- 0.003
1.479<abs(eta)<2.0       Aeff(NH) = 0.037 +/- 0.001      Aeff(γ) = 0.048 +/- 0.001       Aeff(γ+NH) = 0.085 +/- 0.002
2.0<abs(eta)<2.2         Aeff(NH) = 0.023 +/- 0.001      Aeff(γ) = 0.089 +/- 0.002       Aeff(γ+NH) = 0.11 +/- 0.003
2.2<abs(eta)<2.3         Aeff(NH) = 0.023 +/- 0.002      Aeff(γ) = 0.092 +/- 0.004       Aeff(γ+NH) = 0.12 +/- 0.004
2.3<abs(eta)<2.4         Aeff(NH) = 0.021 +/- 0.002      Aeff(γ) = 0.097 +/- 0.004       Aeff(γ+NH) = 0.12 +/- 0.005
2.4< abs(eta)<123213     Aeff(NH) = 0.021 +/- 0.003      Aeff(γ) = 0.11 +/- 0.004        Aeff(γ+NH) = 0.13 +/- 0.006

*/
  for (unsigned int  it = 0; it< electrons.size(); ++it){
float mincor=0.0;
float minrho=0.0;
float rho = std::max(aux.puInfo.rho25Iso,minrho);
float eta=electrons[it].p4.Eta();
float areagamma=0.5;
float areaNH=0.5;
float areaComb=0.5;

if(fabs(eta) <= 1.0 ) {areagamma=0.14; areaNH=0.044; areaComb=0.18;}
if(fabs(eta) > 1.0 &&  fabs(eta) <= 1.479 ) {areagamma=0.13; areaNH=0.065; areaComb=0.20;}
if(fabs(eta) > 1.479 &&  fabs(eta) <= 2.0 ) {areagamma=0.079; areaNH=0.068; areaComb=0.15;}
if(fabs(eta) > 2.0 &&  fabs(eta) <= 2.2 ) {areagamma=0.13; areaNH=0.057; areaComb=0.19;}
if(fabs(eta) > 2.2 &&  fabs(eta) <= 2.3 ) {areagamma=0.15; areaNH=0.058; areaComb=0.21;}
if(fabs(eta) > 2.3 &&  fabs(eta) <= 2.4 ) {areagamma=0.16; areaNH=0.061; areaComb=0.22;}
if(fabs(eta) > 2.4  ) {areagamma=0.18; areaNH=0.11; areaComb=0.29;}
/*
if(fabs(eta) <= 1.0 ) {areagamma=0.081; areaNH=0.024; areaComb=0.10;}
if(fabs(eta) > 1.0 &&  fabs(eta) <= 1.479 ) {areagamma=0.084; areaNH=0.037; areaComb=0.12;}
if(fabs(eta) > 1.479 &&  fabs(eta) <= 2.0 ) {areagamma=0.048; areaNH=0.037; areaComb=0.085;}
if(fabs(eta) > 2.0 &&  fabs(eta) <= 2.2 ) {areagamma=0.089; areaNH=0.023; areaComb=0.11;}
if(fabs(eta) > 2.2 &&  fabs(eta) <= 2.3 ) {areagamma=0.092; areaNH=0.023; areaComb=0.12;}
if(fabs(eta) > 2.3 &&  fabs(eta) <= 2.4 ) {areagamma=0.097; areaNH=0.021; areaComb=0.12;}
if(fabs(eta) > 2.4  ) {areagamma=0.11; areaNH=0.021; areaComb=0.13;}
*/

//Correct electron photon double count
float pho=electrons[it].pfPhoIso;
if(electrons[it].innerHits>0) 
{ 
 pho-=electrons[it].pfPhoIsoDoubleCounted;
}

float pfCorrIso = (electrons[it].pfChaIso+ std::max(pho-rho*areagamma,mincor )+std::max(electrons[it].pfNeuIso-rho*areaNH,mincor))/electrons[it].p4.Pt();
float iso=pfCorrIso;
float id=electrons[it].mvaOutTrig;
bool wp70=((fabs(eta) < 0.8 && id>0.977 && iso < 0.093) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.956 && iso < 0.095) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.966 && iso < 0.171));
bool wp80=((fabs(eta) < 0.8 && id>0.913 && iso < 0.105) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.964 && iso < 0.178) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.899 && iso < 0.150));
bool wp85=((fabs(eta) < 0.8 && id>0.929 && iso < 0.135) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.931 && iso < 0.159) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.805 && iso < 0.155));
bool wp90=((fabs(eta) < 0.8 && id>0.877 && iso < 0.177) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.794 && iso < 0.180) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.846 && iso < 0.244));
bool wp95=((fabs(eta) < 0.8 && id>0.858 && iso < 0.253) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.425 && iso < 0.225) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.759 && iso < 0.308));
bool wpHWW=((fabs(eta) < 0.8 && id>0.94 && iso < 0.15) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.85 && iso < 0.15) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.92 && iso < 0.15));


    if (

	// fake
/*   "(isEE || isEB) && !isEBEEGap &&"
  " (chargedHadronIso + neutralHadronIso + photonIso)/pt <0.10 &&"
   "dB < 0.02 && "  #dB is computed wrt PV but is transverse only, no info about dZ(vertex) 
   "( "
   "(isEE && ("
   "abs(deltaEtaSuperClusterTrackAtVtx) < 0.005 &&  abs(deltaPhiSuperClusterTrackAtVtx) < 0.02 && sigmaIetaIeta < 0.03 && hadronicOverEm < 0.10 &&  abs(1./ecalEnergy*(1.-eSuperClusterOverP)) < 0.05 "
   ")) || " 
    "(isEB && (  "
    "abs(deltaEtaSuperClusterTrackAtVtx) < 0.004 &&  abs(deltaPhiSuperClusterTrackAtVtx) < 0.03 && sigmaIetaIeta < 0.01 && hadronicOverEm < 0.12 && abs(1./ecalEnergy*(1.-eSuperClusterOverP)) < 0.05"
    "))"
#or use mvaNonTrigV0 and mvaTrigV0
    ")" */
//	(fabs(electrons[it].id95 - 7)) < 0.1  &&
        wp95 &&
	 (
         (electrons[it].isEE  &&
        //fabs(electrons[it].Deta) < 0.009 &&
        //fabs(electrons[it].Dphi) < 0.1 &&
        electrons[it].sihih < 0.03  &&
        electrons[it].HoE < 0.10  &&
        electrons[it].innerHits == 0  &&
        (electrons[it].tIso/electrons[it].p4.Pt()) < 0.2 && 
        (electrons[it].eIso/electrons[it].p4.Pt()) < 0.2 &&
        (electrons[it].hIso/electrons[it].p4.Pt()) < 0.2) 
	    || 
        (electrons[it].isEB &&
        //fabs(electrons[it].Deta) < 0.007 &&
        //fabs(electrons[it].Dphi) < 0.015 &&
        electrons[it].sihih < 0.01  &&
        electrons[it].HoE < 0.12  &&
        electrons[it].innerHits == 0  &&
        (electrons[it].tIso/electrons[it].p4.Pt()) < 0.2 && 
        (electrons[it].eIso/electrons[it].p4.Pt()) < 0.2 &&
        (electrons[it].hIso/electrons[it].p4.Pt()) < 0.2) 
        ) && 
//2012 cut based ELE ID
  /*      fabs(electrons[it].dxy) < 0.02  &&
        fabs(electrons[it].dz) < 0.1  &&
        ((electrons[it].isEE  &&
        fabs(electrons[it].Deta) < 0.005 &&
        fabs(electrons[it].Dphi) < 0.02 &&
        electrons[it].sihih < 0.03  &&
        electrons[it].HoE < 0.10  &&
        fabs(electrons[it].fMVAVar_IoEmIoP) < 0.05  
	) || 
        (electrons[it].isEB &&
       fabs(electrons[it].Deta) < 0.004 &&
        fabs(electrons[it].Dphi) < 0.03 &&
        electrons[it].sihih < 0.01  &&
        electrons[it].HoE < 0.12  &&
        fabs(electrons[it].fMVAVar_IoEmIoP) < 0.05 
         )
        ) &&
        pfCorrIso < 0.1 &&
	fabs(electrons[it].p4.Eta()) < 2.5 &&
*/
//Remove this workaround as now we have the proper flags
//	!( fabs(electrons[it].p4.Eta()) < 1.57 && fabs(electrons[it].p4.Eta()) > 1.44) &&
	electrons[it].p4.Pt()>20 //  I use the minimum ok for both Z and W
	){
/*      std::cout << "dxy .dz .isEE  .Deta  .Dphi  .sihih  .HoE  .fMVAVar_IoEmIoP  .isEB " <<std::endl;
      std::cout << fabs(electrons[it].dxy)  << " " <<  fabs(electrons[it].dz) << " " << electrons[it].isEE << " " << 
       fabs(electrons[it].Deta) 
<< " " << fabs(electrons[it].Dphi)
<< " " <<electrons[it].sihih
<< " " <<electrons[it].HoE
<< " " <<fabs(electrons[it].fMVAVar_IoEmIoP)
<< " " <<electrons[it].isEB << std::endl;
*/

      out.push_back(electrons[it]);
      positions.push_back(it);
    }  
  }
    if (verbose_){
      std::cout <<" CandidateFinder: Input Electrons = "<<electrons.size()<<" Output Electrons = "<<out.size()<<std::endl;
    }


}

