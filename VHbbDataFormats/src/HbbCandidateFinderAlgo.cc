#include "VHbbAnalysis/VHbbDataFormats/interface/HbbCandidateFinderAlgo.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidateTools.h"
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
  //
  // search for leptons
  //
  std::vector<VHbbEvent::MuonInfo> mu;
  std::vector<unsigned int> muPos;
  findMuons(event->muInfo,mu, muPos);
  std::vector<VHbbEvent::ElectronInfo> ele;
  std::vector<unsigned int> elePos;
  findElectrons(event->eleInfo,ele, elePos);
  
  std::vector<VHbbEvent::METInfo> met;
  findMET(event->pfmet, met);

  if (verbose_){
    std::cout <<" Electrons: "<< ele.size()<<std::endl;
    std::cout <<" Muons    : "<< mu.size()<<std::endl;
    std::cout <<" MET      : "<< met.size()<<std::endl;
  }
  if (ele.size()<1 && mu.size() < 1 && met.size()<1) return;

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
  temp.additionalJets = addJets;
  temp.V.mets = met;
  temp.V.muons = mu;
  temp.V.electrons = ele;
  
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
   if (jets[i].p4.Pt()> jetPtThreshold && fabs(jets[i].p4.Eta()) < etaThr && jetID(jets[i])){
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
	muons[it].nHits >= 10 &&
        //tracker
	(muons[it].cat & 0x1) && 
	//global
	(muons[it].cat & 0x2) && 
//TODO ?	muons[it].validMuStations >=2 &&
	muons[it].ipDb<.2 &&
	//	(muons[it].hIso+muons[it].eIso+muons[it].tIso)/muons[it].p4.Pt()<.15 &&
	(muons[it].pfChaIso+muons[it].pfPhoIso+muons[it].pfNeuIso)/muons[it].p4.Pt()<.15 &&
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
	(fabs(electrons[it].id95r - 7)) < 0.1  &&
	fabs(electrons[it].p4.Eta()) < 2.5 &&
	!( fabs(electrons[it].p4.Eta()) < 1.57 && fabs(electrons[it].p4.Eta()) > 1.44) &&
	electrons[it].p4.Pt()>20 //  I use the minimum ok for both Z and W
	){
      out.push_back(electrons[it]);
      positions.push_back(it);
    }  
  }
    if (verbose_){
      std::cout <<" CandidateFinder: Input Electrons = "<<electrons.size()<<" Output Electrons = "<<out.size()<<std::endl;
    }


}

