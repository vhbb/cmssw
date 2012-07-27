#include "VHbbAnalysis/VHbbDataFormats/interface/HbbCandidateFinderAlgo.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidateTools.h"
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



void HbbCandidateFinderAlgo::run (const VHbbEvent* event, std::vector<VHbbCandidate>  & candidates){
  //
  // first find the jets
  //

  VHbbCandidateTools selector(verbose_);

  VHbbEvent::SimpleJet j1,j2;
  std::vector<VHbbEvent::SimpleJet> addJets;
  bool foundJets;
  if (useHighestPtHiggs_ == false){
    foundJets = findDiJets(event->simpleJets2,j1,j2,addJets) ;
  }else{
    foundJets= findDiJetsHighestPt(event->simpleJets2,j1,j2,addJets) ;
  }

  if (verbose_){
    std::cout <<" Found Dijets: "<<foundJets<< " Additional: "<<addJets.size()<< std::endl;
  }
  
  if (foundJets == false) return;

  bool foundHardJets;
  VHbbEvent::HardJet fatj1;
  std::vector<VHbbEvent::SimpleJet> subJetsout;
  foundHardJets= findFatJet(event->hardJets,event->subJets,event->filterJets,fatj1,subJetsout) ;

//  if (foundHardJets == false) return;


  //
  // search for leptons
  //
  std::vector<VHbbEvent::MuonInfo> mu;
  std::vector<unsigned int> muPos;
  findMuons(event->muInfo,mu, muPos);
  std::vector<VHbbEvent::ElectronInfo> ele;
  std::vector<unsigned int> elePos;
  findElectrons(event->eleInfo,ele, elePos);
  std::vector<VHbbEvent::TauInfo> tau;
  std::vector<unsigned int> tauPos;
  findTaus(event->tauInfo,tau, tauPos);
  
  std::vector<VHbbEvent::METInfo> met;
  findMET(event->pfmet, met);

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
  temp.H.jets.push_back(j1);
  temp.H.jets.push_back(j2);
  temp.H.p4 = (j1).p4+(j2).p4;
  TVector3 higgsBoost;
  higgsBoost = (temp.H.p4).BoostVector();
  temp.H.helicities.push_back(selector.getHelicity(j1,higgsBoost));
  temp.H.helicities.push_back(selector.getHelicity(j2,higgsBoost));
  temp.H.deltaTheta = selector.getDeltaTheta(j1,j2);

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
  	//Zemu
	result = selector. getHZemuCandidate(temp,ok,elePos);
	if ( ok == true ){
		result.setCandidateType(VHbbCandidate::Zemu);
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
    if(j.chargedEmEFraction > 0.99) return false;
    if(j.chargedHadronEFraction == 0) return false;
    if(j.ntracks == 0) return false;
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


bool HbbCandidateFinderAlgo::findFatJet (const std::vector<VHbbEvent::HardJet>& jetsin, const std::vector<VHbbEvent::SimpleJet>& subjetsin, const std::vector<VHbbEvent::SimpleJet>& filterjetsin, VHbbEvent::HardJet& fatj1,std::vector<VHbbEvent::SimpleJet>& subJetsout){

  if (verbose_){
    std::cout <<" CandidateFinder: Input Jets = "<<jetsin.size()<<std::endl;
  }
  if (jetsin.size()<1) return false;
  
//  float etaThr = 2.5;
  std::vector<VHbbEvent::HardJet> hardjets = jetsin;
  std::vector<VHbbEvent::SimpleJet> subjets = subjetsin;
  std::vector<VHbbEvent::SimpleJet> filterjets = filterjetsin;

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
     if(subjets[subJetIn1st[j]].p4.Pt()>30.) nPt++;}

//     if(nBtag<2 || nPt<2) continue;
     if(nPt<2) continue;

//      if(subjets[subJetIn1st[0]].csv+subjets[subJetIn1st[1]].csv>minBtag1){
//      minBtag1=subjets[subJetIn1st[0]].csv+subjets[subJetIn1st[1]].csv;       
       if((subjets[subJetIn1st[0]].p4+subjets[subJetIn1st[1]].p4).Pt()>minBtag1){
       minBtag1=(subjets[subJetIn1st[0]].p4+subjets[subJetIn1st[1]].p4).Pt();
       fatj1=hardjets[i];
       subJetsout.clear();
       if(subJetIn[0]!=-99) subJetsout.push_back(filterjets[subJetIn[0]]);
       if(subJetIn[1]!=-99) subJetsout.push_back(filterjets[subJetIn[1]]);
       if(subJetIn[2]!=-99) subJetsout.push_back(filterjets[subJetIn[2]]);
       }  

  } // loop hard jet

  
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


void HbbCandidateFinderAlgo::findMuons(const std::vector<VHbbEvent::MuonInfo>& muons, std::vector<VHbbEvent::MuonInfo>& out, std::vector<unsigned int>& positions){
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
	//	(muons[it].hIso+muons[it].eIso+muons[it].tIso)/muons[it].p4.Pt()<.15 &&
	(muons[it].pfChaIso+muons[it].pfPhoIso+muons[it].pfNeuIso)/muons[it].p4.Pt()<.15  &&
	fabs(muons[it].p4.Eta())<2.4 &&
	muons[it].p4.Pt()>8 ) {
      out.push_back(muons[it]);
      positions.push_back(it);
  }
  }
    if (verbose_){
      std::cout <<" CandidateFinder: Input Muons = "<<muons.size()<<" Output Muons = "<<out.size()<<std::endl;
    }


}


void HbbCandidateFinderAlgo::findElectrons(const std::vector<VHbbEvent::ElectronInfo>& electrons, std::vector<VHbbEvent::ElectronInfo>& out,  std::vector<unsigned int>& positions){
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
  for (unsigned int  it = 0; it< electrons.size(); ++it){
    if (
	// fake
	(fabs(electrons[it].id95 - 7)) < 0.1  &&
	fabs(electrons[it].p4.Eta()) < 2.5 &&
//Remove this workaround as now we have the proper flags
//	!( fabs(electrons[it].p4.Eta()) < 1.57 && fabs(electrons[it].p4.Eta()) > 1.44) &&
	electrons[it].p4.Pt()>8 //  I use the minimum ok for both Z and W
	){
      out.push_back(electrons[it]);
      positions.push_back(it);
    }  
  }
    if (verbose_){
      std::cout <<" CandidateFinder: Input Electrons = "<<electrons.size()<<" Output Electrons = "<<out.size()<<std::endl;
    }


}

