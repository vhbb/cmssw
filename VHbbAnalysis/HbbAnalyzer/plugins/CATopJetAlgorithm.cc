// Original author: Brock Tweedie (JHU)
// Ported to CMSSW by: Sal Rappoccio (JHU)
// $Id: CATopJetAlgorithm.cc,v 1.1 2011/06/08 17:25:32 tboccali Exp $

#include "VHbbAnalysis/HbbAnalyzer/interface/CATopJetAlgorithm.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/angle.h"


using namespace std;
using namespace reco;
using namespace edm;


//  Run the algorithm
//  ------------------
void CATopJetAlgorithm::run( const vector<fastjet::PseudoJet> & cell_particles, 
							vector<CompoundPseudoJet> & hardjetsOutput,
							edm::EventSetup const & c )  
{
	if ( verbose_ ) cout << "Welcome to CATopSubJetAlgorithm::run" << endl;
	
	// Sum Et of the event
	double sumEt = 0.;
	
	//make a list of input objects ordered by ET and calculate sum et
	// list of fastjet pseudojet constituents
	for (unsigned i = 0; i < cell_particles.size(); ++i) {
		sumEt += cell_particles[i].perp();
	}
	
	// Determine which bin we are in for et clustering
	int sumEtBinId = -1;
	for ( unsigned int i = 0; i < sumEtBins_.size(); ++i ) {
		if ( sumEt > sumEtBins_[i] ) sumEtBinId = i;
	}
	if ( verbose_ ) cout << "Using sumEt = " << sumEt << ", bin = " << sumEtBinId << endl;
	
	// If the sum et is too low, exit
	if ( sumEtBinId < 0 ) {    
		return;
	}
	
	// empty 4-vector
	fastjet::PseudoJet blankJetA(0,0,0,0);
	blankJetA.set_user_index(-1);
	const fastjet::PseudoJet blankJet = blankJetA;
	
	// Define adjacency variables which depend on which sumEtBin we are in
	double deltarcut = deltarBins_[sumEtBinId];
	double nCellMin = nCellBins_[sumEtBinId];
	
	if ( verbose_ )cout<<"useAdjacency_ = "<<useAdjacency_<<endl;
	if ( verbose_ && useAdjacency_==0)cout<<"No Adjacency"<<endl;
	if ( verbose_ && useAdjacency_==1)cout<<"using deltar adjacency"<<endl;
	if ( verbose_ && useAdjacency_==2)cout<<"using modified adjacency"<<endl;
	if ( verbose_ && useAdjacency_==3)cout<<"using calorimeter nearest neighbor based adjacency"<<endl;
	if ( verbose_ && useAdjacency_==1)cout << "Using deltarcut = " << deltarcut << endl;
	if ( verbose_ && useAdjacency_==3)cout << "Using nCellMin = " << nCellMin << endl;
	
	
	// Define strategy, recombination scheme, and jet definition
	fastjet::Strategy strategy = fastjet::Best;
	fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
	
	// pick algorithm
	fastjet::JetAlgorithm algorithm = static_cast<fastjet::JetAlgorithm>( algorithm_ );
	
	fastjet::JetDefinition jetDef( algorithm, 
								  rBins_[sumEtBinId], recombScheme, strategy);
	
	if ( verbose_ ) cout << "About to do jet clustering in CA" << endl;
	// run the jet clustering
	fastjet::ClusterSequence clusterSeq(cell_particles, jetDef);
	
	if ( verbose_ ) cout << "Getting inclusive jets" << endl;
	// Get the transient inclusive jets
//	vector<fastjet::PseudoJet> inclusiveJets = clusterSeq.inclusive_jets(ptMin_);
        vector<fastjet::PseudoJet> inclusiveJets = fastjet::sorted_by_pt(clusterSeq.inclusive_jets(ptMin_));
	
	if ( verbose_ ) cout << "Getting central jets" << endl;
	// Find the transient central jets
	vector<fastjet::PseudoJet> centralJets;
	for (unsigned int i = 0; i < inclusiveJets.size(); i++) {
		
//		if (inclusiveJets[i].perp() > ptMin_ && fabs(inclusiveJets[i].rapidity()) < centralEtaCut_) {
                if (inclusiveJets[i].perp() > ptMin_ && fabs(inclusiveJets[i].eta()) < centralEtaCut_) {
			centralJets.push_back(inclusiveJets[i]);
		}
	}
	// Sort the transient central jets in Et
////	GreaterByEtPseudoJet compEt;
////	sort( centralJets.begin(), centralJets.end(), compEt );
	
	// These will store the 4-vectors of each hard jet
	vector<math::XYZTLorentzVector> p4_hardJets;
	
	// These will store the indices of each subjet that 
	// are present in each jet
	vector<vector<int> > indices( centralJets.size() );
	
	// Loop over central jets, attempt to find substructure
	vector<fastjet::PseudoJet>::iterator jetIt = centralJets.begin(),
    centralJetsBegin = centralJets.begin(),
    centralJetsEnd = centralJets.end();
	if ( verbose_ )cout<<"Loop over jets"<<endl;
	int i=0;
	for ( ; jetIt != centralJetsEnd; ++jetIt ) {
		if ( verbose_ )cout<<"\nJet "<<i<<endl;
		i++;
		fastjet::PseudoJet localJet = *jetIt;
		
		// Get the 4-vector for this jet
		p4_hardJets.push_back( math::XYZTLorentzVector(localJet.px(), localJet.py(), localJet.pz(), localJet.e() ));
		
		// jet decomposition.  try to find 3 or 4 hard, well-localized subjets, characteristic of a boosted top.
		if ( verbose_ )cout<<"local jet pt = "<<localJet.perp()<<endl;
		if ( verbose_ )cout<<"deltap = "<<ptFracBins_[sumEtBinId]<<endl;
		
		double ptHard = ptFracBins_[sumEtBinId]*localJet.perp();
		vector<fastjet::PseudoJet> leftoversAll;
		
		// stage 1:  primary decomposition.  look for when the jet declusters into two hard subjets
		if ( verbose_ ) cout << "Doing decomposition 1" << endl;
		fastjet::PseudoJet ja, jb;
		vector<fastjet::PseudoJet> leftovers1;
		bool hardBreak1 = decomposeJet(localJet,clusterSeq,cell_particles,ptHard,nCellMin,deltarcut,ja,jb,leftovers1);
		leftoversAll.insert(leftoversAll.end(),leftovers1.begin(),leftovers1.end());
		
/*		// stage 2:  secondary decomposition.  look for when the hard subjets found above further decluster into two hard sub-subjets
		//
		// ja -> jaa+jab ?
		if ( verbose_ ) cout << "Doing decomposition 2. ja->jaa+jab?" << endl;
		fastjet::PseudoJet jaa, jab;
		vector<fastjet::PseudoJet> leftovers2a;
		bool hardBreak2a = false;
		if (hardBreak1)  hardBreak2a = decomposeJet(ja,clusterSeq,cell_particles,ptHard,nCellMin,deltarcut,jaa,jab,leftovers2a);
		leftoversAll.insert(leftoversAll.end(),leftovers2a.begin(),leftovers2a.end());
		// jb -> jba+jbb ?
		if ( verbose_ ) cout << "Doing decomposition 2. ja->jba+jbb?" << endl;
		fastjet::PseudoJet jba, jbb;
		vector<fastjet::PseudoJet> leftovers2b;
		bool hardBreak2b = false;
		if (hardBreak1)  hardBreak2b = decomposeJet(jb,clusterSeq,cell_particles,ptHard,nCellMin,deltarcut,jba,jbb,leftovers2b);
		leftoversAll.insert(leftoversAll.end(),leftovers2b.begin(),leftovers2b.end());
*/		
////                double DeltaR12=TMath::Sqrt((ja.eta()-jb.eta())*(ja.eta()-jb.eta())+(ja.phi()-jb.phi())*(ja.phi()-jb.phi()));
////                double y=TMath::Min(ja.perp2(),jb.perp2())*DeltaR12*DeltaR12/(localJet.m()*localJet.m());
 
                double DeltaR12=TMath::Sqrt(ja.squared_distance(jb));
                double y=ja.kt_distance(jb)/localJet.m2();  
      
                if(y>0.09 && hardBreak1){
 
                fastjet::PseudoJet hardSubA = ja; 
                fastjet::PseudoJet hardSubB = jb; 
 
                double Rfilt=TMath::Min(0.3,DeltaR12/2.);
                fastjet::JetDefinition jetDef2( algorithm,
                                  Rfilt, recombScheme, strategy);
 
                if ( verbose_ ) cout << "About to do sub jet clustering in CA with R=Rfilt" << endl;
//              vector<fastjet::PseudoJet>  hcand; 
//              hcand.push_back(localJet);
                vector<fastjet::PseudoJet> hcandConstituents = clusterSeq.constituents( localJet );
// run the jet clustering
                fastjet::ClusterSequence clusterSubSeq(hcandConstituents, jetDef2);
 
               if ( verbose_ ) cout << "Getting inclusive sub jets" << endl;
// Get the transient inclusive jets
//               vector<fastjet::PseudoJet> inclusiveSubJets = clusterSubSeq.inclusive_jets(ptMin_);
               vector<fastjet::PseudoJet> inclusiveSubJets = fastjet::sorted_by_pt(clusterSubSeq.inclusive_jets(ptMin_));
////               vector<fastjet::PseudoJet> inclusiveSubJets = fastjet::sorted_by_pt(clusterSubSeq.inclusive_jets());

               if ( verbose_ ) cout << "Getting central sub jets" << endl;
// Find the transient central jets
               vector<fastjet::PseudoJet> centralSubJets;
               for (unsigned int i = 0; i < inclusiveSubJets.size(); i++) {
//                if (inclusiveSubJets[i].perp() > ptMin_ && fabs(inclusiveSubJets[i].rapidity()) < centralEtaCut_) {
                if (inclusiveSubJets[i].perp() > ptMin_ && fabs(inclusiveSubJets[i].eta()) < centralEtaCut_) {
//v2                if (inclusiveSubJets[i].perp() > ptMin_) {
                centralSubJets.push_back(inclusiveSubJets[i]);
                }
               }     
      
//////              sort(centralSubJets.begin(), centralSubJets.end(), compEt );
 
              if ( hardSubA.user_index() >= 0 ) centralSubJets.push_back(hardSubA);
              if ( hardSubB.user_index() >= 0 ) centralSubJets.push_back(hardSubB);

		
		// create the subjets objects to put into the "output" objects
		vector<CompoundPseudoSubJet>  subjetsOutput;
		std::vector<fastjet::PseudoJet>::const_iterator itSubJetBegin = centralSubJets.begin(),
		itSubJet = itSubJetBegin, itSubJetEnd = centralSubJets.end();
		for (; itSubJet != itSubJetEnd; ++itSubJet ){
			//       if ( verbose_ ) cout << "Adding input collection element " << (*itSubJet).user_index() << endl;
			//       if ( (*itSubJet).user_index() >= 0 && (*itSubJet).user_index() < cell_particles.size() )
			
			// Get the transient subjet constituents from fastjet
			vector<fastjet::PseudoJet> subjetFastjetConstituents = clusterSeq.constituents( *itSubJet );
			
			// Get the indices of the constituents
			vector<int> constituents;
			
			// Loop over the constituents and get the indices
			vector<fastjet::PseudoJet>::const_iterator fastSubIt = subjetFastjetConstituents.begin(),
			transConstBegin = subjetFastjetConstituents.begin(),
			transConstEnd = subjetFastjetConstituents.end();
			for ( ; fastSubIt != transConstEnd; ++fastSubIt ) {
				if ( fastSubIt->user_index() >= 0 && static_cast<unsigned int>(fastSubIt->user_index()) < cell_particles.size() ) {
					constituents.push_back( fastSubIt->user_index() );
				}
			}
			
			// Make a CompoundPseudoSubJet object to hold this subjet and the indices of its constituents
			subjetsOutput.push_back( CompoundPseudoSubJet( *itSubJet, constituents ) );
		}
		
		
		// Make a CompoundPseudoJet object to hold this hard jet, and the subjets that make it up
		hardjetsOutput.push_back( CompoundPseudoJet( *jetIt, subjetsOutput ) );
		
	 }  // y cut
	}  // loop central jets
}




//-----------------------------------------------------------------------
// determine whether two clusters (made of calorimeter towers) are living on "adjacent" cells.  if they are, then
// we probably shouldn't consider them to be independent objects!
//
// From Sal: Ignoring genjet case
//
bool CATopJetAlgorithm::adjacentCells(const fastjet::PseudoJet & jet1, const fastjet::PseudoJet & jet2, 
									  const vector<fastjet::PseudoJet> & cell_particles,
									  const fastjet::ClusterSequence & theClusterSequence,
									  double nCellMin ) const {
	
	
	double eta1 = jet1.rapidity();
	double phi1 = jet1.phi();
	double eta2 = jet2.rapidity();
	double phi2 = jet2.phi();
	
	double deta = abs(eta2 - eta1) / 0.087;
	double dphi = fabs( reco::deltaPhi(phi2,phi1) ) / 0.087;
	
	return ( ( deta + dphi ) <= nCellMin );
}



//-------------------------------------------------------------------------
// attempt to decompose a jet into "hard" subjets, where hardness is set by ptHard
//
bool CATopJetAlgorithm::decomposeJet(const fastjet::PseudoJet & theJet, 
									 const fastjet::ClusterSequence & theClusterSequence, 
									 const vector<fastjet::PseudoJet> & cell_particles,
									 double ptHard, double nCellMin, double deltarcut,
									 fastjet::PseudoJet & ja, fastjet::PseudoJet & jb, 
									 vector<fastjet::PseudoJet> & leftovers) const {
	
	bool goodBreak;
	fastjet::PseudoJet j = theJet;
	double InputObjectPt = j.perp();
	if ( verbose_ )cout<<"Input Object Pt = "<<InputObjectPt<<endl;
	if ( verbose_ )cout<<"ptHard = "<<ptHard<<endl;
	leftovers.clear();
	if ( verbose_ )cout<<"start while loop"<<endl;
	
	while (1) {                                                      // watch out for infinite loop!
		goodBreak = theClusterSequence.has_parents(j,ja,jb);
		if (!goodBreak){
			if ( verbose_ )cout<<"bad break. this is one cell. can't decluster anymore."<<endl;
			break;         // this is one cell, can't decluster anymore
		}
		
		if ( verbose_ )cout<<"good break. ja Pt = "<<ja.perp()<<" jb Pt = "<<jb.perp()<<endl;
		
		/// Adjacency Requirement ///
		
		// check if clusters are adjacent using a constant deltar adjacency.
		double clusters_deltar=fabs(ja.eta()-jb.eta())+fabs(deltaPhi(ja.phi(),jb.phi()));
		
		if ( verbose_  && useAdjacency_ ==1)cout<<"clusters_deltar = "<<clusters_deltar<<endl;
		if ( verbose_  && useAdjacency_ ==1)cout<<"deltar cut = "<<deltarcut<<endl;
		
		if ( useAdjacency_==1 && clusters_deltar < deltarcut){
			if ( verbose_ )cout<<"clusters too close. consant adj. break."<<endl;
			break;
		} 
		
		// Check if clusters are adjacent using a DeltaR adjacency which is a function of pT.
		double clusters_deltaR=deltaR( ja.eta(), ja.phi(), jb.eta(), jb.phi() );
		
		if ( verbose_  && useAdjacency_ ==2)cout<<"clusters_deltaR = "<<clusters_deltaR<<endl;
		if ( verbose_  && useAdjacency_ ==2)cout<<"0.4-0.0004*InputObjectPt = "<<0.4-0.0004*InputObjectPt<<endl;
		
		if ( useAdjacency_==2 && clusters_deltaR < 0.4-0.0004*InputObjectPt)
		{
			if ( verbose_ )cout<<"clusters too close. modified adj. break."<<endl;
			break;
		} 

		// Check if clusters are adjacent in the calorimeter. 
		if ( useAdjacency_==3 &&  adjacentCells(ja,jb,cell_particles,theClusterSequence,nCellMin) ){                  
			if ( verbose_ )cout<<"clusters too close in the calorimeter. calorimeter adj. break."<<endl;
			break;         // the clusters are "adjacent" in the calorimeter => shouldn't have decomposed
		}
				
		if ( verbose_ )cout<<"clusters pass distance cut"<<endl;
		
		/// Pt Fraction Requirement ///
		
		if ( verbose_ )cout<<"ptHard = "<<ptHard<<endl;
 
                if(ja.m() > jb.m() && ja.m()< 0.67*j.m() ) return true;
                else if(jb.m() > ja.m() && jb.m()< 0.67*j.m() ) return true;

                else if (ja.m() > jb.m()) {                              // broke into one hard and one soft, ditch the soft one and try again
                j = ja;
                vector<fastjet::PseudoJet> particles = theClusterSequence.constituents(jb);
                leftovers.insert(leftovers.end(),particles.begin(),particles.end());
                }
               else {
               j = jb;
               vector<fastjet::PseudoJet> particles = theClusterSequence.constituents(ja);
               leftovers.insert(leftovers.end(),particles.begin(),particles.end());
               }

		
	}
	
	if ( verbose_ )cout<<"did not decluster."<<endl;  // did not decluster into hard subjets
	
	ja.reset(0,0,0,0);
	jb.reset(0,0,0,0);
	leftovers.clear();
	return false;
}
