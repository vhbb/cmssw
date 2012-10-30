#include <TH1F.h>
#include <TH3F.h>
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "VHbbAnalysis/VHbbDataFormats/interface/HbbCandidateFinderAlgo.h" 
#include "VHbbAnalysis/VHbbDataFormats/src/HbbCandidateFinderAlgo.cc"

#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEventAuxInfo.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/TriggerReader.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/TopMassReco.h"

//for IVF
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include <DataFormats/GeometrySurface/interface/Surface.h>
#include "Math/SMatrix.h"

#include "AnalysisDataFormats/TauAnalysis/interface/PFMEtSignCovMatrix.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"
#include <TMatrixD.h>

//Move class definition to Ntupler.h ?
//#include "VHbbAnalysis/VHbbDataFormats/interface/Ntupler.h"

//for SimBhad: FIXME
//#include "VHbbAnalysis/BAnalysis/interface/SimBHadron.h"
//btagging
#include "VHbbAnalysis/VHbbDataFormats/interface/BTagWeight.h"
//trigger
#include "VHbbAnalysis/VHbbDataFormats/interface/TriggerWeight.h"

#include <sstream>
#include <string>

#define MAXJ 30
#define MAXL 50
#define MAXB 10
#define MAXT 60
#define nMetUnc 24 

#define GENPTOLOR(a) TLorentzVector((a).px(), (a).py(), (a).pz(), (a).energy())
#define GENPTOLORP(a) TLorentzVector((a)->px(), (a)->py(), (a)->pz(),(a)->energy())

//eleEnDown/Up, muEn, tauEn, JES, JER, Unclustered for type1 MET [0-11] and than type1p2 MET [12-23]

TLorentzVector VisibleP(const reco::Candidate & p) {
	TLorentzVector mcp4 = GENPTOLOR(p);
	if (p.numberOfDaughters() == 0) {
		if (abs(p.pdgId()) == 12 || abs(p.pdgId()) == 14 || abs(p.pdgId()) == 16) {
			TLorentzVector zeroVector(0.,0.,0.,0.);
			return zeroVector;
		} else {
			return mcp4;
		}
	} else { // has daughters
		TLorentzVector sum(0.,0.,0.,0.);
		for (unsigned int i = 0 ; i < p.numberOfDaughters() ; i++) {
			const reco::Candidate *d = p.daughter(i);
			sum += VisibleP(*d);
		}
		return sum;
	}
}// VisibleP

double evalEt( double pt, double eta, double phi, double e){
	TLorentzVector j;
	j.SetPtEtaPhiE(pt,eta,phi, e );
	return j.Et(); 
	
}


double evalMt( double pt, double eta, double phi, double e){
	TLorentzVector j;
	j.SetPtEtaPhiE(pt,eta,phi, e );
	return j.Mt(); 
	
}


struct CompareDeltaR {
	CompareDeltaR(TLorentzVector p4dir_): p4dir(p4dir_) {}
	bool operator()( const VHbbEvent::SimpleJet& j1, const  VHbbEvent::SimpleJet& j2 ) const {
		return j1.p4.DeltaR(p4dir) > j2.p4.DeltaR(p4dir);
	}
	TLorentzVector p4dir;
};

const GlobalVector flightDirection(const TVector3 pv, const reco::Vertex &sv){
	GlobalVector fdir(sv.position().X() - pv.X(),
					  sv.position().Y() - pv.Y(),
					  sv.position().Z() - pv.Z());
	return fdir;
}

bool jsonContainsEvent (const std::vector< edm::LuminosityBlockRange > &jsonVec,
                        const edm::EventBase &event)
{
	// if the jsonVec is empty, then no JSON file was provided so all
	// events should pass
	if (jsonVec.empty())
		{
			return true;
		}
	bool (* funcPtr) (edm::LuminosityBlockRange const &,
					  edm::LuminosityBlockID const &) = &edm::contains;
	edm::LuminosityBlockID lumiID (event.id().run(), 
								   event.id().luminosityBlock());
	std::vector< edm::LuminosityBlockRange >::const_iterator iter = 
    std::find_if (jsonVec.begin(), jsonVec.end(),
				  boost::bind(funcPtr, _1, lumiID) );
	return jsonVec.end() != iter;
	
}

float resolutionBias(float eta)
{
	// return 0;//Nominal!
	if(eta< 1.1) return 0.05;
	if(eta< 2.5) return 0.10;
	if(eta< 5) return 0.30;
	return 0;
}


//FIXME : need to update EDM ntuple with SimBhadron infos
// typedef struct
// {
//   void set(const SimBHadron & sbhc, int i){
//     mass[i] = sbhc.mass();
//     pt[i] = sbhc.pt();
//     eta[i] = sbhc.eta();
//     phi[i] = sbhc.phi();
//     vtx_x[i] = sbhc.decPosition.x();
//     vtx_y[i] = sbhc.decPosition.y();
//     vtx_z[i] = sbhc.decPosition.z();
//     pdgId[i] = sbhc.pdgId();
//     status[i] = sbhc.status();
//   };
//   void reset(){
//     for(int i=0; i < MAXB; ++i){
//       mass[i] = -99; pt[i] = -99; eta[i] = -99; phi[i] = -99; vtx_x[i] = -99; vtx_y[i] = -99; vtx_z[i] = -99; pdgId[i] = -99; status[i] = -99;
//     }
//   };
//   float mass[MAXB];
//   float pt[MAXB];
//   float eta[MAXB];
//   float phi[MAXB];
//   float vtx_x[MAXB];
//   float vtx_y[MAXB];
//   float vtx_z[MAXB];
//   int pdgId[MAXB];
//   int status[MAXB];
// //   int quarkStatus[MAXB];
// //   int brotherStatus[MAXB];
// //   int otherId[MAXB];
// //   bool etaOk[MAXB];
// //   bool simOk[MAXB];
// //   bool trackOk[MAXB];
// //   bool cutOk[MAXB];
// //   bool cutNewOk[MAXB];
// //   bool mcMatchOk[MAXB];
// //   bool matchOk[MAXB];
// } SimBHadronInfo;


typedef struct
	{
		void set( const reco::SecondaryVertex & recoSv, const TVector3 recoPv, int isv){
			pt[isv]   = recoSv.p4().Pt();
			eta[isv]  = flightDirection(recoPv,recoSv).eta();
			phi[isv]  = flightDirection(recoPv,recoSv).phi();
			massBcand[isv] = recoSv.p4().M();
			massSv[isv] = recoSv.p4().M();
			dist3D[isv] = recoSv.dist3d().value();  
			distSig3D[isv] = recoSv.dist3d().significance(); 
			dist2D[isv] = recoSv.dist2d().value();  
			distSig2D[isv] = recoSv.dist2d().significance(); 
			dist3D_norm[isv] = recoSv.dist3d().value()/recoSv.p4().Gamma(); 
		};
		void reset(){
			for(int i = 0; i < MAXB; ++i){
				massBcand[i] = -99; massSv[i]= -99; pt[i] = -99; eta[i] = -99; phi[i] = -99; dist3D[i] = -99; distSig3D[i] = -99; dist2D[i] = -99; distSig2D[i] = -99; dist3D_norm[i] = -99;
			}
		};
		float massBcand[MAXB];
		float massSv[MAXB];
		float pt[MAXB]; 
		float eta[MAXB]; 
		float phi[MAXB];
		float dist3D[MAXB]; 
		float distSig3D[MAXB];
		float dist2D[MAXB]; 
		float distSig2D[MAXB];
		float dist3D_norm[MAXB];
	} IVFInfo;


typedef struct 
	{
		float mass; 
		float pt;
		float eta;
		float phi;
		float dR;
		float dPhi;
		float dEta;
	} HiggsInfo;

typedef struct 
	{
		float mass; 
		float pt;
		float eta;
		float phi;
		float status;
		float charge;
		float momid;
	} genParticleInfo;



typedef struct 
	
	{
		float bmass; 
		float bpt;
		float beta;
		float bphi;
		float bstatus;
		float wdau1mass; 
		float wdau1pt;
		float wdau1eta;
		float wdau1phi;
		float wdau1id;
		float wdau2mass; 
		float wdau2pt;
		float wdau2eta;
		float wdau2phi;
		float wdau2id;
		
		
	} genTopInfo;



typedef struct
	{
		bool FatHiggsFlag; 
		float mass;
		float pt;
		float eta;
		float phi;
		float filteredmass;
		float filteredpt; 
		//  float dR;
		//  float dPhi;
		//  float dEta;
	} FatHiggsInfo;

typedef struct
	{
		bool HiggsCSVtkSharing;
		bool HiggsIPtkSharing;
		bool HiggsSVtkSharing;
		bool FatHiggsCSVtkSharing;
		bool FatHiggsIPtkSharing;
		bool FatHiggsSVtkSharing;
	} TrackSharingInfo;

typedef struct 
	{
		float mass;  //MT in case of W
		float pt;
		float eta;
		float phi;
	} TrackInfo;


struct  LeptonInfo
{
	void reset()
	{
		for(int i =0; i < MAXL;i++){ 
			mass[i]=-99; pt[i]=-99; eta[i]=-99; phi[i]=-99; aodCombRelIso[i]=-99; pfCombRelIso[i]=-99; photonIso[i]=-99; neutralHadIso[i]=-99; chargedHadIso[i]=-99; chargedPUIso[i]=-99; particleIso[i]=-99; dxy[i]=-99; dz[i]=-99; type[i]=-99;  genPt[i]=-99; genEta[i]=-99; genPhi[i]=-99;  
			id80[i]=-99; id95[i]=-99; vbtf[i]=-99; id80NoIso[i]=-99;
			acop[i]=-99; 
			charge[i]=-99;wp70[i]=-99; wp80[i]=-99;wp85[i]=-99;wp90[i]=-99;wp95[i]=-99;wpHWW[i]=-99;
			pfCorrIso[i]=-99.; id2012tight[i]=-99; idMVAnotrig[i]=-99; idMVAtrig[i]=-99;innerHits[i]=-99.;
			decayModeFinding[i] =-99.;
			byLooseCombinedIsolationDeltaBetaCorr[i] =-99.;
			againstMuonTight[i] =-99.;
			againstElectronLoose[i] =-99.;
			againstElectronMedium[i] =-99.;
			againstElectronMVA[i] =-99.;
			NsignalPFChargedHadrCands[i] =-99;
			NsignalPFGammaCands[i] =-99;
			leadPFChargedHadrCandPt[i] = -99.;
			byLooseIsolation[i]  = -99.;
			byMediumIsolation[i] = -99.;
			byTightIsolation[i] = -99.;
		}
	}
	
	template <class Input> void set(const Input & i, int j,int t)
	{
		type[j]=t;
		pt[j]=i.p4.Pt(); 
		mass[j]=i.p4.M();
		eta[j]=i.p4.Eta();
		phi[j]=i.p4.Phi();
		aodCombRelIso[j]=(i.hIso+i.eIso+i.tIso)/i.p4.Pt();
		pfCombRelIso[j]=(i.pfChaIso+i.pfPhoIso+i.pfNeuIso)/i.p4.Pt();
		photonIso[j]=i.pfPhoIso;
		neutralHadIso[j]=i.pfNeuIso;
		chargedHadIso[j]=i.pfChaIso;
		chargedPUIso[j]=i.pfChaPUIso;
		charge[j]=i.charge;
		acop[j]=i.acop;
		if(i.mcFourMomentum.Pt() > 0)
			{ 
				genPt[j]=i.mcFourMomentum.Pt();
				genEta[j]=i.mcFourMomentum.Eta();
				genPhi[j]=i.mcFourMomentum.Phi();
			}
		setSpecific(i,j);
	}
	template <class Input> void setSpecific(const Input & i, int j)
	{
	}      
	
	
	
	
	float mass[MAXL];  //MT in case of W
	float pt[MAXL];
	float eta[MAXL];
	float phi[MAXL];
	float aodCombRelIso[MAXL];
	float pfCombRelIso[MAXL];
	float photonIso[MAXL];
	float neutralHadIso[MAXL];
	float chargedHadIso[MAXL];
	float chargedPUIso[MAXL];
	float particleIso[MAXL];
	float genPt[MAXL];
	float genEta[MAXL];
	float genPhi[MAXL];
	float dxy[MAXL];
	float dz[MAXL];
	int type[MAXL];
	float id80[MAXL];
	float id95[MAXL];
	float vbtf[MAXL];
	float id80NoIso[MAXL];
	float charge[MAXL];
	float acop[MAXL];
	
	float pfCorrIso[MAXL];
	float id2012tight[MAXL];
	float idMVAnotrig[MAXL];
	float idMVAtrig[MAXL];
	float idMVApresel[MAXL];
	float wp70[MAXL]; 
	float wp80[MAXL]; 
	float wp85[MAXL]; 
	float wp90[MAXL]; 
	float wp95[MAXL]; 
	float wpHWW[MAXL];
	float innerHits[MAXL]; 
	float photonIsoDoubleCount[MAXL];
	float decayModeFinding[MAXL],byLooseCombinedIsolationDeltaBetaCorr[MAXL],againstMuonTight[MAXL],againstElectronLoose[MAXL],againstElectronMedium[MAXL],againstElectronMVA[MAXL];
	int NsignalPFChargedHadrCands[MAXL], NsignalPFGammaCands[MAXL];
	float leadPFChargedHadrCandPt[MAXL], byLooseIsolation[MAXL], byMediumIsolation[MAXL], byTightIsolation[MAXL]; 	
};

template <> void LeptonInfo::setSpecific<VHbbEvent::ElectronInfo>(const VHbbEvent::ElectronInfo & i, int j){
	id80[j]=i.id80;
	id95[j]=i.id95;
	id80NoIso[j]=(i.innerHits ==0 && !(fabs(i.convDist)<0.02 && fabs(i.convDcot)<0.02) &&
				  ((i.isEB && i.sihih<0.01 && fabs(i.Dphi)<0.06 && fabs(i.Deta)<0.004) || (i.isEE && i.sihih<0.03 && fabs(i.Dphi)<0.03  && fabs(i.Deta)<0.007)));
}
template <> void LeptonInfo::setSpecific<VHbbEvent::MuonInfo>(const VHbbEvent::MuonInfo & i, int j){
	dxy[j]=i.ipDb;
	dz[j]=i.zPVPt;
	vbtf[j]=( i.globChi2<10 && i.nPixelHits>= 1 && i.globNHits != 0 && i.nHits > 10 && i.cat & 0x1 && i.cat & 0x2 && i.nMatches >=2 && i.ipDb<.2 &&
			 (i.pfChaIso+i.pfPhoIso+i.pfNeuIso)/i.p4.Pt()<.15  && fabs(i.p4.Eta())<2.4 && i.p4.Pt()>20 ) ;
	
}
template <> void LeptonInfo::setSpecific<VHbbEvent::TauInfo>(const VHbbEvent::TauInfo &i, int j) {
	//  std::cout << "In tau set specific" << std::endl;
	
	decayModeFinding[j] = i.decayModeFinding;
	byLooseCombinedIsolationDeltaBetaCorr[j] =i.byLooseCombinedIsolationDeltaBetaCorr;
	againstMuonTight[j] =i.againstMuonTight;
	againstElectronLoose[j] =i.againstElectronLoose;
	againstElectronMedium[j] =i.againstElectronMedium;
	againstElectronMVA[j] =i.againstElectronMVA;
	NsignalPFChargedHadrCands[j] =i.NsignalPFChargedHadrCands;
	NsignalPFGammaCands[j] =i.NsignalPFGammaCands;
	leadPFChargedHadrCandPt[j] = i.leadPFChargedHadrCandPt;
	byLooseIsolation[j] = i.byLooseIsolation;
	byMediumIsolation[j] = i.byMediumIsolation;
	byTightIsolation[j] = i.byTightIsolation;
}


typedef struct 
	{
		float et; 
		float sumet;   
		float sig;
		float phi;
	} METInfo;


typedef struct 
	{
		
		float  et[nMetUnc]; float phi[nMetUnc]; float sumet[nMetUnc]; 
		template <class Input> void set(const Input & i, int j)
		{
			et[j]=i.p4.Pt(); 
			phi[j]=i.p4.Phi();
			sumet[j]=i.sumEt;
			
		}      
	} METUncInfo ;


typedef struct 
	{
		float mht;
		float ht;  
		float sig;
		float phi;
	} MHTInfo;

typedef struct 
	{
		float mass;
		float pt;
		float wMass;
	} TopInfo;

typedef struct 
	{
		int run;
		int lumi;
		int event;
		int json;
	} EventInfo;

typedef struct 
	{
		void set(const VHbbEvent::SimpleJet & j, int i) 
		{
			pt[i]=j.p4.Pt();
			eta[i]=j.p4.Eta();
			phi[i]=j.p4.Phi();
			e[i]=j.p4.E();
			Mt[i]=evalEt(j.p4.Pt(),j.p4.Eta(),j.p4.Phi(),j.p4.E());
			Et[i]=evalMt(j.p4.Pt(),j.p4.Eta(),j.p4.Phi(),j.p4.E());
			csv[i]=j.csv;
			csvivf[i]=j.csvivf;
			cmva[i]=j.cmva;
			numTracksSV[i] = j.vtxNTracks;
			vtxMass[i]= j.vtxMass;
			vtxPt[i]= j.vtxP4.Pt();
			vtxEta[i]= j.vtxP4.Eta();
			vtxPhi[i]= j.vtxP4.Phi();
			vtxE[i]= j.vtxP4.E();			
			vtx3dL[i] = j.vtx3dL;
			vtx3deL[i] = j.vtx3deL;
			chf[i]=j.chargedHadronEFraction;
			nhf[i]  =j.neutralHadronEFraction;
			cef[i]  =j.chargedEmEFraction;
			nef[i]  =j.neutralEmEFraction;
			nconstituents[i]  =j.nConstituents;
			nch[i]=j.ntracks;
			SF_CSVL[i]=j.SF_CSVL;
			SF_CSVM[i]=j.SF_CSVM;
			SF_CSVT[i]=j.SF_CSVT;
			SF_CSVLerr[i]=j.SF_CSVLerr;
			SF_CSVMerr[i]=j.SF_CSVMerr;
			SF_CSVTerr[i]=j.SF_CSVTerr;
			
			flavour[i]=j.flavour;
			isSemiLeptMCtruth[i]=j.isSemiLeptMCtruth;
			isSemiLept[i]=j.isSemiLept;
			SoftLeptpdgId[i] = j.SoftLeptpdgId;
			SoftLeptIdlooseMu[i] = j.SoftLeptIdlooseMu;
			SoftLeptId95[i] = j.SoftLeptId95;
			SoftLeptPt[i] = j.SoftLeptPt;
			SoftLeptdR[i] = j.SoftLeptdR;  
			SoftLeptptRel[i] = j.SoftLeptptRel;
			SoftLeptRelCombIso[i] = j.SoftLeptRelCombIso;
			if(j.bestMCp4.Pt() > 0)
				{
					genPt[i]=j.bestMCp4.Pt();
					genEta[i]=j.bestMCp4.Eta();
					genPhi[i]=j.bestMCp4.Phi();
				}
			JECUnc[i]=j.jecunc;
			id[i]=jetId(i);
			ptRaw[i]=j.ptRaw;
			ptLeadTrack[i]=j.ptLeadTrack; 
			charge[i]=j.charge;     
			jetArea[i]=j.jetArea;
			
		}
		bool jetId(int i)
		{
			if(nhf[i] > 0.99) return false;
			if(nef[i] > 0.99) return false;
			if(nconstituents[i]  <= 1) return false;
			if(fabs(eta[i])<2.5) {
				if(cef[i] > 0.99) return false;
				if(chf[i] == 0) return false;
				if(nch[i]== 0) return false;
			}
			return true;
		}
		void reset()
		{
			for(int i=0;i<MAXJ;i++) {
				pt[i]=-99; eta[i]=-99; phi[i]=-99;e[i]=-99;csv[i]=-99,Et[i]=-99,Mt[i]=-99;
				csvivf[i]=-99; cmva[i]=-99;
				cosTheta[i]=-99; numTracksSV[i]=-99; chf[i]=-99; nhf[i]=-99; cef[i]=-99; nef[i]=-99; nch[i]=-99; nconstituents[i]=-99; flavour[i]=-99; isSemiLeptMCtruth[i]=-99; isSemiLept[i]=-99;      
				SoftLeptpdgId[i] = -99; SoftLeptIdlooseMu[i] = -99;  SoftLeptId95[i] =  -99;   SoftLeptPt[i] = -99;  SoftLeptdR[i] = -99;   SoftLeptptRel[i] = -99; SoftLeptRelCombIso[i] = -99;  
				genPt[i]=-99; genEta[i]=-99; genPhi[i]=-99; JECUnc[i]=-99; ptRaw[i]=-99.; ptLeadTrack[i]=-99.; selectedTauDR[i] = -99.; charge[i]=-99; jetArea[i]=-99;
			}
		}
		float pt[MAXJ];
		float eta[MAXJ];
		float phi[MAXJ];
		float e[MAXJ];
		float Mt[MAXJ];
		float Et[MAXJ];
		float csv[MAXJ];
		float csvivf[MAXJ];
		float cmva[MAXJ];
		float cosTheta[MAXJ];
		int numTracksSV[MAXJ];
		float chf[MAXJ];
		float nhf[MAXJ];
		float cef[MAXJ];
		float nef[MAXJ];
		float nch[MAXJ];
		float nconstituents[MAXJ];
		float flavour[MAXJ];
		int isSemiLept[MAXJ];
		int isSemiLeptMCtruth[MAXJ];
		int SoftLeptpdgId[MAXJ] ;
		int SoftLeptIdlooseMu[MAXJ] ;  
		int SoftLeptId95[MAXJ]  ;   
		float SoftLeptPt[MAXJ] ;  
		float SoftLeptdR[MAXJ] ;   
		float SoftLeptptRel[MAXJ] ; 
		float SoftLeptRelCombIso[MAXJ];
		float genPt[MAXJ];
		float genEta[MAXJ];
		float genPhi[MAXJ];
		float JECUnc[MAXJ];
		float vtxMass[MAXJ];
		float vtxPt[MAXJ];
		float vtxEta[MAXJ];
		float vtxPhi[MAXJ];
		float vtxE[MAXJ];		
		float vtx3dL [MAXJ];
		float vtx3deL[MAXJ];
		bool id[MAXJ];
		float SF_CSVL[MAXJ];
		float SF_CSVM[MAXJ];
		float SF_CSVT[MAXJ]; 
		float SF_CSVLerr[MAXJ];
		float SF_CSVMerr[MAXJ];
		float SF_CSVTerr[MAXJ];  
		float ptRaw[MAXJ];
		float ptLeadTrack[MAXJ];
		float selectedTauDR[MAXJ];
		float charge[MAXJ];
		float jetArea[MAXJ];
		
	} JetInfo;

int main(int argc, char* argv[]) 
{
	gROOT->Reset();
	
	
	
	TTree *_outTree;
	IVFInfo IVF;
	//FIXME
	//  SimBHadronInfo SimBs;
	float rho,rho25;
	int nPVs;
	METInfo MET;
	METInfo fakeMET;
	METInfo METnoPU;
	METInfo METnoPUCh;
	METInfo METtype1corr;
	METInfo METtype1p2corr;
	METInfo METnoPUtype1corr;
	METInfo METnoPUtype1p2corr;
	
	MHTInfo MHT;
	
	METUncInfo  metUnc;
	
	
	TopInfo top;
	EventInfo EVENT;
	//  JetInfo jet1,jet2, addJet1, addJet2;
	// lepton1,lepton2;
	JetInfo hJets, aJets, fathFilterJets;
	LeptonInfo vLeptons, aLeptons;
	int naJets=0, nhJets=0, nfathFilterJets=0;
	HiggsInfo H,SVH; //add here the fatjet higgs
	FatHiggsInfo FatH;
	genParticleInfo genZ, genZstar, genWstar, genW,  genH, genB, genBbar, genTau, genTauMu; //add here the fatjet higgs
	genTopInfo genTop, genTbar;
	TrackInfo V;
	int nvlep=0,nalep=0; 
	
	TrackSharingInfo TkSharing; // track sharing info;
	
	float HVdPhi,HVMass,HMETdPhi,VMt,VMte,VMtmu,deltaPullAngle, deltaPullAngleAK7, deltaPullAngle2;
	float deltaPullAngle2AK7, gendrcc, gendrbb, genZpt, genWpt, genHpt, weightTrig, weightTrigMay, weightTrigV4;
	float weightTrigMET, weightTrigOrMu30, minDeltaPhijetMET,  jetPt_minDeltaPhijetMET, PUweight, PUweight2011B;
	float PU0,PUp1,PUm1;
	float CovMartixfirst, CovMartixsecond, CovMatrixCorrelation, minDeltaPhijetMETZtau, svmass;
	float PzetaVis, PzetaMiss, Pzeta25, Pzeta85, Pzeta45, Pzeta65;
	bool genZtaumu, genWmunu, genWenu, genZmumu, genZee, genZnunu;
	
	float tau_decayModeFinding[200], tau_byLooseCombinedIsolationDeltaBetaCorr[200], tau_againstMuonTight[200], tau_againstElectronLoose[200], tau_againstElectronMedium[200], tau_againstElectronMVA[200];
	for (int t=0; t!=200; t++) tau_decayModeFinding[t]=tau_byLooseCombinedIsolationDeltaBetaCorr[t]=tau_againstMuonTight[t]=tau_againstElectronLoose[t]=tau_againstElectronMedium[t]=tau_againstElectronMVA[t]=-1.;
	int nRecoTaus;
	
	float weightEleRecoAndId,weightEleTrigJetMETPart, weightEleTrigElePart,weightEleTrigEleAugPart;
	float  weightTrigMET80, weightTrigMET100,    weightTrig2CJet20 , weightTrigMET150  , weightTrigMET802CJet, weightTrigMET1002CJet, weightTrigMETLP ;
	
	//FIXME
	//  float SimBs_dr = -99, SimBs_dEta= -99, SimBs_dPhi = -99, SimBs_Hmass = -99, SimBs_Hpt = -99, SimBs_Heta = -99, SimBs_Hphi = -99;
	int WplusMode,WminusMode;
	int Vtype,nSvs,numJets,numBJets,eventFlav;
	//   bool isMET80_CJ80, ispfMHT150, isMET80_2CJ20,isMET65_2CJ20, isJETID,isIsoMu17;
	bool triggerFlags[500],hbhe,ecalFlag,totalKinematics,  cschaloFlag,  hcallaserFlag,   trackingfailureFlag ;
	
	
	float btag1T2CSF=1.,btag2TSF=1.,btag1TSF=1.,btagA0CSF=1., btagA0TSF=1., btag2CSF=1., btag1TA1C=1.;
	// ----------------------------------------------------------------------
	// First Part: 
	//
	//  * enable the AutoLibraryLoader 
	//  * book the histograms of interest 
	//  * open the input file
	// ----------------------------------------------------------------------
	
	// load framework libraries
	gSystem->Load("libFWCoreFWLite");
	gSystem->Load("libDataFormatsFWLite");
	AutoLibraryLoader::enable();
	
	// parse arguments
	if ( argc < 2 ) {
		return 0;
	}
	
	std::vector<VHbbCandidate> * candZlocal = new std::vector<VHbbCandidate>;
	std::vector<VHbbCandidate> * candWlocal = new std::vector<VHbbCandidate>;
	
	// get the python configuration
	PythonProcessDesc builder(argv[1]);
	const edm::ParameterSet& in  = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput" );
	const edm::ParameterSet& out = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteOutput");
	const edm::ParameterSet& ana = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("Analyzer");
	
	std::vector<edm::LuminosityBlockRange> jsonVector;
	if ( in.exists("lumisToProcess") ) 
		{
			std::vector<edm::LuminosityBlockRange> const & lumisTemp =
			in.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> > ("lumisToProcess");
			jsonVector.resize( lumisTemp.size() );
			copy( lumisTemp.begin(), lumisTemp.end(), jsonVector.begin() );
		}
	
	// now get each parameter
	int maxEvents_( in.getParameter<int>("maxEvents") );
	int skipEvents_( in.getParameter<int>("skipEvents") );
	int runMin_( in.getParameter<int>("runMin") );
	int runMax_( in.getParameter<int>("runMax") );
	bool Verbage = ana.getParameter<bool>("verbose");
	std::string outputFile_( out.getParameter<std::string>("fileName" ) );
	std::vector<std::string> triggers( ana.getParameter<std::vector<std::string> >("triggers") );
	double btagThr =  ana.getParameter<double>("bJetCountThreshold" );
	bool fromCandidate = ana.getParameter<bool>("readFromCandidates");
	bool useHighestPtHiggsZ = ana.getParameter<bool>("useHighestPtHiggsZ");
	bool useHighestPtHiggsW = ana.getParameter<bool>("useHighestPtHiggsW");
	HbbCandidateFinderAlgo * algoZ = new HbbCandidateFinderAlgo(ana.getParameter<bool>("verbose"), ana.getParameter<double>("jetPtThresholdZ"),useHighestPtHiggsZ);
	HbbCandidateFinderAlgo * algoW = new HbbCandidateFinderAlgo(ana.getParameter<bool>("verbose"), ana.getParameter<double>("jetPtThresholdW"),useHighestPtHiggsW );
	HbbCandidateFinderAlgo * algoRecoverLowPt = new HbbCandidateFinderAlgo(ana.getParameter<bool>("verbose"), 15, true);
	
	TriggerWeight triggerWeight(ana);
	BTagWeight btag(2); // 2 operating points "Custom" = 0.5 and "Tight = 0.898"
	BTagSampleEfficiency btagEff( ana.getParameter<std::string>("btagEffFileName" ).c_str() ); 
	
	std::vector<std::string> inputFiles_( in.getParameter<std::vector<std::string> >("fileNames") );
	//  std::string inputFile( in.getParameter<std::string> ("fileName") );
	
	std::string PUmcfileName_ = in.getParameter<std::string> ("PUmcfileName") ;
	std::string PUmcfileName2011B_ = in.getParameter<std::string> ("PUmcfileName2011B") ;
	
	std::string PUdatafileName_ = in.getParameter<std::string> ("PUdatafileName") ;
	std::string PUdatafileName2011B_ = in.getParameter<std::string> ("PUdatafileName2011B") ;
	std::string Weight3DfileName_ = in.getParameter<std::string> ("Weight3DfileName") ;
	
	
	
	bool isMC_( ana.getParameter<bool>("isMC") );  
	TriggerReader trigger(isMC_);
	TriggerReader patFilters(false);
	
	edm::LumiReWeighting   lumiWeights;
	edm::Lumi3DReWeighting   lumiWeights2011B;
	if(isMC_)
		{
			lumiWeights = edm::LumiReWeighting(PUmcfileName_,PUdatafileName_ , "pileup", "pileup");
			
			lumiWeights2011B = edm::Lumi3DReWeighting(PUmcfileName2011B_,PUdatafileName2011B_ , "pileup", "pileup");
			if(Weight3DfileName_!="")
				{ lumiWeights2011B.weight3D_init(Weight3DfileName_.c_str()); }
			else
				{
					lumiWeights2011B.weight3D_init(1.0); // generate the weights the fisrt time;
				}
			
		}
	
	//   TFile *_outPUFile	= new TFile((outputFile_+"_PU").c_str(), "recreate");	
	TFile *_outFile	= new TFile(outputFile_.c_str(), "recreate");	
	TH1F *  count = new TH1F("Count","Count", 1,0,2 );
	TH1F *  countWithPU = new TH1F("CountWithPU","CountWithPU", 1,0,2 );
	TH1F *  countWithPU2011B = new TH1F("CountWithPU2011B","CountWithPU2011B", 1,0,2 );
	TH3F *  input3DPU = new TH3F("Input3DPU","Input3DPU", 36,-0.5,35.5,36,-0.5,35.5, 36,-0.5,35.5 );
	
	TH1F * pu = new TH1F("pileup","",51,-0.5,50.5);
	_outTree = new TTree("tautree", "myTree");
	
	_outTree->Branch("H"		,  &H	            ,  "mass/F:pt/F:eta:phi/F:dR/F:dPhi/F:dEta/F");
	_outTree->Branch("V"		,  &V	            ,  "mass/F:pt/F:eta:phi/F");
	_outTree->Branch("FatH"          ,  &FatH               ,  "FatHiggsFlag/I:mass/F:pt/F:eta:phi/F:filteredmass/F:filteredpt/F");
	_outTree->Branch("genZ"		,  &genZ	            ,  "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
	_outTree->Branch("genZstar"		,  &genZstar	            ,  "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
	_outTree->Branch("genW"		,  &genW	            ,  "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
	_outTree->Branch("genWstar"		,  &genWstar	            ,  "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
	_outTree->Branch("genH"		,  &genH	            ,  "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
	_outTree->Branch("genB"		,  &genB	            ,  "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
	_outTree->Branch("genBbar"		,  &genBbar	            ,  "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
	
	_outTree->Branch("genTop"		,  &genTop	            ,  "bmass/F:bpt/F:beta:bphi/F:bstatus/F:wdau1mass/F:wdau1pt/F:wdau1eta:wdau1phi/F:wdau1id/F:wdau2mass/F:wdau2pt/F:wdau2eta:wdau2phi/F:wdau2id/F");
	_outTree->Branch("genTbar"		,  &genTbar	            ,  "bmass/F:bpt/F:beta:bphi/F:bstatus/F:wdau1mass/F:wdau1pt/F:wdau1eta:wdau1phi/F:wdau1id/F:wdau2mass/F:wdau2pt/F:wdau2eta:wdau2phi/F:wdau2id/F");
	
	_outTree->Branch("genTau"		,  &genTau	            ,  "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
	
	_outTree->Branch("genTauMu"		,  &genTauMu	            ,  "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
	
	_outTree->Branch("TkSharing", &TkSharing, "HiggsCSVtkSharing/b:HiggsIPtkSharing:HiggsSVtkSharing:FatHiggsCSVtkSharing:FatHiggsIPtkSharing:FatHiggsSVtkSharing");
	
	// only for testing tau reco:
	_outTree->Branch("nRecoTaus", &nRecoTaus, "nRecoTaus/I");
	_outTree->Branch("tau_decayModeFinding", &tau_decayModeFinding, "tau_decayModeFinding[nRecoTaus]/F");
	_outTree->Branch("tau_byLooseCombinedIsolationDeltaBetaCorr", &tau_byLooseCombinedIsolationDeltaBetaCorr, "tau_byLooseCombinedIsolationDeltaBetaCorr[nRecoTaus]/F");
	_outTree->Branch("tau_againstMuonTight", &tau_againstMuonTight, "tau_againstMuonTight[nRecoTaus]/F");
	_outTree->Branch("tau_againstElectronLoose", &tau_againstElectronLoose, "tau_againstElectronLoose[nRecoTaus]/F");
	_outTree->Branch("tau_againstElectronMedium", &tau_againstElectronMedium, "tau_againstElectronMedium[nRecoTaus]/F");
	_outTree->Branch("tau_againstElectronMVA", &tau_againstElectronMVA, "tau_againstElectronMVA[nRecoTaus]/F");
	
	
	_outTree->Branch("genZtaumu", &genZtaumu, "genZtaumu/O");
	_outTree->Branch("genWmunu", &genWmunu, "genWmunu/O");
	_outTree->Branch("genWenu", &genWenu, "genWenu/O");
	_outTree->Branch("genZmumu", &genZmumu, "genZmumu/O");
	_outTree->Branch("genZee", &genZee, "genZee/O");
	_outTree->Branch("genZnunu", &genZnunu, "genZnunu/O");
	
	_outTree->Branch("nhJets"		,  &nhJets	            ,  "nhJets/I");
	_outTree->Branch("nfathFilterJets",   &nfathFilterJets,   "nfathFilterJets/I");
	_outTree->Branch("naJets"		,  &naJets	            ,  "naJets/I");
	
	
	
	_outTree->Branch("hJet_pt",hJets.pt ,"pt[nhJets]/F");
	_outTree->Branch("hJet_eta",hJets.eta ,"eta[nhJets]/F");
	_outTree->Branch("hJet_phi",hJets.phi ,"phi[nhJets]/F");
	_outTree->Branch("hJet_e",hJets.e ,"e[nhJets]/F");
	_outTree->Branch("hJet_Et",hJets.Et ,"Et[nhJets]/F");
	_outTree->Branch("hJet_Mt",hJets.Mt ,"Mt[nhJets]/F");
	_outTree->Branch("hJet_csv",hJets.csv ,"csv[nhJets]/F");
	_outTree->Branch("hJet_csvivf",hJets.csvivf ,"csvivf[nhJets]/F");
	_outTree->Branch("hJet_cmva",hJets.cmva ,"cmva[nhJets]/F");
	_outTree->Branch("hJet_cosTheta",hJets.cosTheta ,"cosTheta[nhJets]/F");
	_outTree->Branch("hJet_numTracksSV",hJets.numTracksSV ,"numTracksSV[nhJets]/I");
	_outTree->Branch("hJet_chf",hJets.chf ,"chf[nhJets]/F");
	_outTree->Branch("hJet_nhf",hJets.nhf ,"nhf[nhJets]/F");
	_outTree->Branch("hJet_cef",hJets.cef ,"cef[nhJets]/F");
	_outTree->Branch("hJet_nef",hJets.nef ,"nef[nhJets]/F");
	_outTree->Branch("hJet_nch",hJets.nch ,"nch[nhJets]/F");
	_outTree->Branch("hJet_nconstituents",hJets.nconstituents ,"nconstituents[nhJets]");
	_outTree->Branch("hJet_flavour",hJets.flavour ,"flavour[nhJets]/F");
	_outTree->Branch("hJet_isSemiLept",hJets.isSemiLept ,"isSemiLept[nhJets]/I");
	_outTree->Branch("hJet_isSemiLeptMCtruth",hJets.isSemiLeptMCtruth ,"isSemiLeptMCtruth[nhJets]/I");
	_outTree->Branch("hJet_SoftLeptpdgId", hJets.SoftLeptpdgId , "SoftLeptpdgId[nhJets]/I");
	_outTree->Branch("hJet_SoftLeptIdlooseMu", hJets.SoftLeptIdlooseMu , "SoftLeptIdlooseMu[nhJets]/I");
	_outTree->Branch("hJet_SoftLeptId95", hJets.SoftLeptId95 , "SoftLeptId95[nhJets]/I");
	_outTree->Branch("hJet_SoftLeptPt", hJets.SoftLeptPt , "SoftLeptPt[nhJets]/F");
	_outTree->Branch("hJet_SoftLeptdR", hJets.SoftLeptdR , "SoftLeptdR[nhJets]/F");
	_outTree->Branch("hJet_SoftLeptptRel", hJets.SoftLeptptRel , "SoftLeptptRel[nhJets]/F");
	_outTree->Branch("hJet_SoftLeptRelCombIso", hJets.SoftLeptRelCombIso , "SoftLeptRelCombIso[nhJets]/F");
	_outTree->Branch("hJet_genPt",hJets.genPt ,"genPt[nhJets]/F");
	_outTree->Branch("hJet_genEta",hJets.genEta ,"genEta[nhJets]/F");
	_outTree->Branch("hJet_genPhi",hJets.genPhi ,"genPhi[nhJets]/F");
	_outTree->Branch("hJet_JECUnc",hJets.JECUnc ,"JECUnc[nhJets]/F");
	_outTree->Branch("hJet_vtxMass",hJets.vtxMass ,"vtxMass[nhJets]/F");
	_outTree->Branch("hJet_vtxPt",hJets.vtxPt ,"vtxPt[nhJets]/F");
	_outTree->Branch("hJet_vtxEta",hJets.vtxEta ,"vtxEta[nhJets]/F");
	_outTree->Branch("hJet_vtxPhi",hJets.vtxPhi ,"vtxPhi[nhJets]/F");
	_outTree->Branch("hJet_vtxE",hJets.vtxE ,"vtxE[nhJets]/F");	
	_outTree->Branch("hJet_vtx3dL",hJets.vtx3dL ,"vtx3dL[nhJets]/F");
	_outTree->Branch("hJet_vtx3deL",hJets.vtx3deL ,"vtx3deL[nhJets]/F");
	_outTree->Branch("hJet_id",hJets.id ,"id[nhJets]/b");
	_outTree->Branch("hJet_SF_CSVL",hJets.SF_CSVL ,"SF_CSVL[nhJets]/b");
	_outTree->Branch("hJet_SF_CSVM",hJets.SF_CSVM ,"SF_CSVM[nhJets]/b");
	_outTree->Branch("hJet_SF_CSVT",hJets.SF_CSVT ,"SF_CSVT[nhJets]/b");
	_outTree->Branch("hJet_SF_CSVLerr",hJets.SF_CSVLerr ,"SF_CSVLerr[nhJets]/b");
	_outTree->Branch("hJet_SF_CSVMerr",hJets.SF_CSVMerr ,"SF_CSVMerr[nhJets]/b");
	_outTree->Branch("hJet_SF_CSVTerr",hJets.SF_CSVTerr ,"SF_CSVTerr[nhJets]/b");
	_outTree->Branch("hJet_ptRaw",hJets.ptRaw ,"ptRaw[nhJets]/F");
	_outTree->Branch("hJet_ptLeadTrack",hJets.ptLeadTrack ,"ptLeadTrack[nhJets]/F");
	_outTree->Branch("hJet_charge",hJets.charge ,"charge[nhJets]/F");
	_outTree->Branch("hJet_jetArea",hJets.jetArea ,"jetArea[nhJets]/F");

	
	_outTree->Branch("fathFilterJets_pt",fathFilterJets.pt ,"pt[nfathFilterJets]/F");
	_outTree->Branch("fathFilterJets_eta",fathFilterJets.eta ,"eta[nfathFilterJets]/F");
	_outTree->Branch("fathFilterJets_phi",fathFilterJets.phi ,"phi[nfathFilterJets]/F");
	_outTree->Branch("fathFilterJets_e",fathFilterJets.e ,"e[nfathFilterJets]/F");
	_outTree->Branch("fathFilterJets_csv",fathFilterJets.csv ,"csv[nfathFilterJets]/F");
	_outTree->Branch("fathFilterJets_csvivf",fathFilterJets.csvivf ,"csvivf[nfathFilterJets]/F");
	_outTree->Branch("fathFilterJets_cmva",fathFilterJets.cmva ,"cmva[nfathFilterJets]/F");
	_outTree->Branch("fathFilterJets_flavour",fathFilterJets.flavour ,"flavour[nfathFilterJets]/F");
	
	
	_outTree->Branch("aJet_pt",aJets.pt ,"pt[naJets]/F");
	_outTree->Branch("aJet_eta",aJets.eta ,"eta[naJets]/F");
	_outTree->Branch("aJet_phi",aJets.phi ,"phi[naJets]/F");
	_outTree->Branch("aJet_e",aJets.e ,"e[naJets]/F");
	_outTree->Branch("aJet_Et",aJets.Et ,"Et[naJets]/F");
	_outTree->Branch("aJet_Mt",aJets.Mt ,"Mt[naJets]/F");
	_outTree->Branch("aJet_csv",aJets.csv ,"csv[naJets]/F");
	_outTree->Branch("aJet_csvivf",aJets.csvivf ,"csvivf[naJets]/F");
	_outTree->Branch("aJet_cmva",aJets.cmva ,"cmva[naJets]/F");
	_outTree->Branch("aJet_cosTheta",aJets.cosTheta ,"cosTheta[naJets]/F");
	_outTree->Branch("aJet_numTracksSV",aJets.numTracksSV ,"numTracksSV[naJets]/I");
	_outTree->Branch("aJet_chf",aJets.chf ,"chf[naJets]/F");
	_outTree->Branch("aJet_nhf",aJets.nhf ,"nhf[naJets]/F");
	_outTree->Branch("aJet_cef",aJets.cef ,"cef[naJets]/F");
	_outTree->Branch("aJet_nef",aJets.nef ,"nef[naJets]/F");
	_outTree->Branch("aJet_nch",aJets.nch ,"nch[naJets]/F");
	_outTree->Branch("aJet_nconstituents",aJets.nconstituents ,"nconstituents[naJets]");
	_outTree->Branch("aJet_flavour",aJets.flavour ,"flavour[naJets]/F");
	_outTree->Branch("aJet_isSemiLept",aJets.isSemiLept ,"isSemiLept[naJets]/I");
	_outTree->Branch("aJet_isSemiLeptMCtruth",aJets.isSemiLeptMCtruth ,"isSemiLeptMCtruth[naJets]/I");
	_outTree->Branch("aJet_SoftLeptpdgId",aJets.SoftLeptpdgId , "SoftLeptpdgId[naJets]/I");
	_outTree->Branch("aJet_SoftLeptIdlooseMu", aJets.SoftLeptIdlooseMu , "SoftLeptIdlooseMu[naJets]/I");
	_outTree->Branch("aJet_SoftLeptId95", aJets.SoftLeptId95 , "SoftLeptId95[naJets]/I");
	_outTree->Branch("aJet_SoftLeptPt", aJets.SoftLeptPt , "SoftLeptPt[naJets]/F");
	_outTree->Branch("aJet_SoftLeptdR", aJets.SoftLeptdR , "SoftLeptdR[naJets]/F");
	_outTree->Branch("aJet_SoftLeptptRel", aJets.SoftLeptptRel , "SoftLeptptRel[naJets]/F");
	_outTree->Branch("aJet_SoftLeptRelCombIso", aJets.SoftLeptRelCombIso , "SoftLeptRelCombIso[naJets]/F");
	_outTree->Branch("aJet_jetArea",aJets.jetArea ,"jetArea[naJets]/F");
	_outTree->Branch("aJet_charge",aJets.charge ,"charge[naJets]/F");
	_outTree->Branch("aJet_selectedTauDR",aJets.selectedTauDR ,"selectedTauDR[naJets]/F");

	_outTree->Branch("aJet_genPt",aJets.genPt ,"genPt[naJets]/F");
	_outTree->Branch("aJet_genEta",aJets.genEta ,"genEta[naJets]/F");
	_outTree->Branch("aJet_genPhi",aJets.genPhi ,"genPhi[naJets]/F");
	_outTree->Branch("aJet_JECUnc",aJets.JECUnc ,"JECUnc[naJets]/F");
	_outTree->Branch("aJet_vtxMass",aJets.vtxMass ,"vtxMass[naJets]/F");
	_outTree->Branch("aJet_vtx3dL",aJets.vtx3dL ,"vtx3dL[naJets]/F");
	_outTree->Branch("aJet_vtx3deL",aJets.vtx3deL ,"vtx3deL[naJets]/F");
	_outTree->Branch("aJet_id",aJets.id ,"id[naJets]/b");
	_outTree->Branch("aJet_SF_CSVL",aJets.SF_CSVL ,"SF_CSVL[naJets]/b");
	_outTree->Branch("aJet_SF_CSVM",aJets.SF_CSVM ,"SF_CSVM[naJets]/b");
	_outTree->Branch("aJet_SF_CSVT",aJets.SF_CSVT ,"SF_CSVT[naJets]/b");
	_outTree->Branch("aJet_SF_CSVLerr",aJets.SF_CSVLerr ,"SF_CSVLerr[naJets]/b");
	_outTree->Branch("aJet_SF_CSVMerr",aJets.SF_CSVMerr ,"SF_CSVMerr[naJets]/b");
	_outTree->Branch("aJet_SF_CSVTerr",aJets.SF_CSVTerr ,"SF_CSVTerr[naJets]/b");
	
	_outTree->Branch("numJets"      ,  &numJets         ,  "numJets/I"       );                
	_outTree->Branch("numBJets"      ,  &numBJets         ,  "numBJets/I"       );                
	_outTree->Branch("deltaPullAngle", &deltaPullAngle  ,  "deltaPullAngle/F");
	_outTree->Branch("deltaPullAngle2", &deltaPullAngle2  ,  "deltaPullAngle2/F");
	_outTree->Branch("gendrcc"    , &gendrcc      ,  "gendrcc/F");
	_outTree->Branch("gendrbb"    , &gendrbb      ,  "gendrbb/F");
	_outTree->Branch("genZpt"    , &genZpt      ,  "genZpt/F");
	_outTree->Branch("genWpt"    , &genWpt      ,  "genWpt/F");
	_outTree->Branch("genHpt"    , &genHpt      ,  "genHpt/F");
	_outTree->Branch("weightTrig"        , &weightTrig          ,  "weightTrig/F");
	_outTree->Branch("weightTrigMay"        , &weightTrigMay          ,  "weightTrigMay/F");
	_outTree->Branch("weightTrigV4"        , &weightTrigV4          ,  "weightTrigV4/F");
	_outTree->Branch("weightTrigMET"        , &weightTrigMET          ,  "weightTrigMET/F");
	_outTree->Branch("weightTrigOrMu30"  , &weightTrigOrMu30    , "weightTrigOrMu30/F");
	_outTree->Branch("weightEleRecoAndId"        , &weightEleRecoAndId     ,  "weightEleRecoAndId/F");
	_outTree->Branch("weightEleTrigJetMETPart"        , &weightEleTrigJetMETPart          ,  "weightEleTrigJetMETPart/F");
	_outTree->Branch("weightEleTrigElePart"        , &weightEleTrigElePart          ,  "weightEleTrigElePart/F");
	_outTree->Branch("weightEleTrigEleAugPart"        , &weightEleTrigEleAugPart          ,  "weightEleTrigEleAugPart/F");
	
	_outTree->Branch("weightTrigMET80"        , &weightTrigMET80          , "weightTrigMET80/F");
	_outTree->Branch("weightTrigMET100"        , &weightTrigMET100          , "weightTrigMET100/F");
	_outTree->Branch("weightTrig2CJet20"        , &weightTrig2CJet20          , "weightTrig2CJet20/F");
	_outTree->Branch("weightTrigMET150"        , &weightTrigMET150          , "weightTrigMET150/F");
	_outTree->Branch("weightTrigMET802CJet"        , &weightTrigMET802CJet          , "weightTrigMET802CJet/F");
	_outTree->Branch("weightTrigMET1002CJet"        , &weightTrigMET1002CJet          , "weightTrigMET1002CJet/F");
	_outTree->Branch("weightTrigMETLP"        , &weightTrigMETLP          , "weightTrigMETLP/F");
	
	
	_outTree->Branch("deltaPullAngleAK7", &deltaPullAngleAK7  ,  "deltaPullAngleAK7/F");
	_outTree->Branch("deltaPullAngle2AK7", &deltaPullAngle2AK7  ,  "deltaPullAngle2AK7/F");
	_outTree->Branch("PU0",       &PU0  ,  "PU0/F");
	_outTree->Branch("PUm1",       &PUm1  ,  "PUm1/F");
	_outTree->Branch("PUp1",       &PUp1  ,  "PUp1/F");
	_outTree->Branch("PUweight",       &PUweight  ,  "PUweight/F");
	_outTree->Branch("PUweight2011B",       &PUweight2011B  ,  "PUweight2011B/F");
	_outTree->Branch("eventFlav",       &eventFlav  ,  "eventFlav/I");
	
	
	
    
	_outTree->Branch("Vtype"     ,  &Vtype   ,   "Vtype/I" );                
	_outTree->Branch("HVdPhi"     ,  &HVdPhi   ,   "HVdPhi/F" );                
	_outTree->Branch("HVMass"     ,  &HVMass   ,   "HVMass/F" );                
	_outTree->Branch("HMETdPhi"     ,  &HMETdPhi   ,   "HMETdPhi/F" );                
	_outTree->Branch("VMt"  	,  &VMt      ,   "VMt/F"    );   
  	_outTree->Branch("VMte"  	,  &VMte      ,   "VMte/F"    );             	
	_outTree->Branch("VMtmu"  	,  &VMtmu      ,   "VMtmu/F"    );             	
	
	_outTree->Branch("nvlep"	,  &nvlep    ,   "nvlep/I");
	_outTree->Branch("nalep"	,  &nalep    ,   "nalep/I");
	
	_outTree->Branch("vLepton_mass",vLeptons.mass ,"mass[nvlep]/F");
	_outTree->Branch("vLepton_pt",vLeptons.pt ,"pt[nvlep]/F");
	_outTree->Branch("vLepton_eta",vLeptons.eta ,"eta[nvlep]");
	_outTree->Branch("vLepton_phi",vLeptons.phi ,"phi[nvlep]/F");
	_outTree->Branch("vLepton_aodCombRelIso",vLeptons.aodCombRelIso ,"aodCombRelIso[nvlep]/F");
	_outTree->Branch("vLepton_pfCombRelIso",vLeptons.pfCombRelIso ,"pfCombRelIso[nvlep]/F");
	_outTree->Branch("vLepton_photonIso",vLeptons.photonIso ,"photonIso[nvlep]/F");
	_outTree->Branch("vLepton_neutralHadIso",vLeptons.neutralHadIso ,"neutralHadIso[nvlep]/F");
	_outTree->Branch("vLepton_chargedHadIso",vLeptons.chargedHadIso ,"chargedHadIso[nvlep]/F");
	_outTree->Branch("vLepton_chargedPUIso",vLeptons.chargedPUIso ,"chargedPUIso[nvlep]/F");
	_outTree->Branch("vLepton_particleIso",vLeptons.particleIso ,"particleIso[nvlep]/F");
	_outTree->Branch("vLepton_dxy",vLeptons.dxy ,"dxy[nvlep]/F");
	_outTree->Branch("vLepton_dz",vLeptons.dz ,"dz[nvlep]/F");
	_outTree->Branch("vLepton_type",vLeptons.type ,"type[nvlep]/I");
	_outTree->Branch("vLepton_id80",vLeptons.id80 ,"id80[nvlep]/F");
	_outTree->Branch("vLepton_id95",vLeptons.id95 ,"id95[nvlep]/F");
	_outTree->Branch("vLepton_vbtf",vLeptons.vbtf ,"vbtf[nvlep]/F");
	_outTree->Branch("vLepton_id80NoIso",vLeptons.id80NoIso ,"id80NoIso[nvlep]/F");
	_outTree->Branch("vLepton_genPt",vLeptons.genPt ,"genPt[nvlep]/F");
	_outTree->Branch("vLepton_genEta",vLeptons.genEta ,"genEta[nvlep]");
	_outTree->Branch("vLepton_genPhi",vLeptons.genPhi ,"genPhi[nvlep]/F");
	_outTree->Branch("vLepton_charge",vLeptons.charge ,"charge[nvlep]/F");
	_outTree->Branch("vLepton_acop",vLeptons.acop ,"acop[nvlep]/F");
	
	_outTree->Branch("aLepton_mass",aLeptons.mass ,"mass[nalep]/F");
	_outTree->Branch("aLepton_pt",aLeptons.pt ,"pt[nalep]/F");
	_outTree->Branch("aLepton_eta",aLeptons.eta ,"eta[nalep]");
	_outTree->Branch("aLepton_phi",aLeptons.phi ,"phi[nalep]/F");
	_outTree->Branch("aLepton_aodCombRelIso",aLeptons.aodCombRelIso ,"aodCombRelIso[nalep]/F");
	_outTree->Branch("aLepton_pfCombRelIso",aLeptons.pfCombRelIso ,"pfCombRelIso[nalep]/F");
	_outTree->Branch("aLepton_photonIso",aLeptons.photonIso ,"photonIso[nalep]/F");
	_outTree->Branch("aLepton_neutralHadIso",aLeptons.neutralHadIso ,"neutralHadIso[nalep]/F");
	_outTree->Branch("aLepton_chargedHadIso",aLeptons.chargedHadIso ,"chargedHadIso[nalep]/F");
	_outTree->Branch("aLepton_chargedPUIso",aLeptons.chargedPUIso ,"chargedPUIso[nalep]/F");
	_outTree->Branch("aLepton_particleIso",aLeptons.particleIso ,"particleIso[nalep]/F");
	_outTree->Branch("aLepton_dxy",aLeptons.dxy ,"dxy[nalep]/F");
	_outTree->Branch("aLepton_dz",aLeptons.dz ,"dz[nalep]/F");
	_outTree->Branch("aLepton_type",aLeptons.type ,"type[nalep]/I");
	_outTree->Branch("aLepton_id80",aLeptons.id80 ,"id80[nalep]/F");
	_outTree->Branch("aLepton_id95",aLeptons.id95 ,"id95[nalep]/F");
	_outTree->Branch("aLepton_vbtf",aLeptons.vbtf ,"vbtf[nalep]/F");
	_outTree->Branch("aLepton_id80NoIso",aLeptons.id80NoIso ,"id80NoIso[nalep]/F");
	_outTree->Branch("aLepton_genPt",aLeptons.genPt ,"genPt[nalep]/F");
	_outTree->Branch("aLepton_genEta",aLeptons.genEta ,"genEta[nalep]");
	_outTree->Branch("aLepton_genPhi",aLeptons.genPhi ,"genPhi[nalep]/F");
	_outTree->Branch("aLepton_charge",aLeptons.charge ,"charge[nalep]/F");
	_outTree->Branch("aLepton_acop",aLeptons.acop ,"acop[nalep]/F");
	
	_outTree->Branch("top"		,  &top	         ,   "mass/F:pt/F:wMass/F");
	_outTree->Branch("WplusMode"		,  &WplusMode	 ,   "WplusMode/I");
	_outTree->Branch("WminusMode"		,  &WminusMode	 ,   "WminusMode/I");
	
	//IVF
	_outTree->Branch("nSvs",&nSvs ,"nSvs/I");
	_outTree->Branch("Sv_massBCand", &IVF.massBcand,"massBcand[nSvs]/F");
	_outTree->Branch("Sv_massSv", &IVF.massSv,"massSv[nSvs]/F");
	_outTree->Branch("Sv_pt", &IVF.pt,"pt[nSvs]/F");
	_outTree->Branch("Sv_eta", &IVF.eta,"eta[nSvs]/F");
	_outTree->Branch("Sv_phi", &IVF.phi,"phi[nSvs]/F");
	_outTree->Branch("Sv_dist3D", &IVF.dist3D,"dist3D[nSvs]/F");
	_outTree->Branch("Sv_dist2D", &IVF.dist2D,"dist2D[nSvs]/F");
	_outTree->Branch("Sv_distSim2D", &IVF.distSig2D,"distSig2D[nSvs]/F");
	_outTree->Branch("Sv_distSig3D", &IVF.distSig3D,"distSig3D[nSvs]/F");
	_outTree->Branch("Sv_dist3D_norm", &IVF.dist3D_norm,"dist3D_norm[nSvs]/F");
	//IVF higgs candidate
	_outTree->Branch("SVH"          ,  &SVH               ,  "mass/F:pt/F:eta:phi/F:dR/F:dPhi/F:dEta/F");
	
	_outTree->Branch("rho"		,  &rho	         ,   "rho/F");
	_outTree->Branch("rho25"		,  &rho25	         ,   "rho25/F");
	_outTree->Branch("nPVs"		,  &nPVs	         ,   "nPVs/I");
	_outTree->Branch("METnoPU"		,  &METnoPU	         ,   "et/F:sumet:sig/F:phi/F");
	_outTree->Branch("METnoPUCh"		,  &METnoPUCh	         ,   "et/F:sumet:sig/F:phi/F");
	_outTree->Branch("MET"		,  &MET	         ,   "et/F:sumet:sig/F:phi/F");
	_outTree->Branch("METtype1corr"		,  &METtype1corr	         ,   "et/F:sumet:sig/F:phi/F");
	_outTree->Branch("METtype1p2corr"		,  &METtype1p2corr	         ,   "et/F:sumet:sig/F:phi/F");
	_outTree->Branch("METnoPUtype1corr"		,  &METnoPUtype1corr	         ,   "et/F:sumet:sig/F:phi/F");
	_outTree->Branch("METnoPUtype1p2corr"		,  &METnoPUtype1p2corr	         ,   "et/F:sumet:sig/F:phi/F");
	
	_outTree->Branch("metUnc_et",&metUnc.et ,"et[24]/F");
	_outTree->Branch("metUnc_phi",&metUnc.phi ,"phi[24]/F");
	_outTree->Branch("metUnc_sumet",&metUnc.sumet ,"sumet[24]/F");
	
	
	
	_outTree->Branch("fakeMET"		,  &fakeMET	         ,   "et/F:sumet:sig/F:phi/F");
	_outTree->Branch("MHT"		,  &MHT	         ,   "mht/F:ht:sig/F:phi/F");
	_outTree->Branch("minDeltaPhijetMET"		,  &minDeltaPhijetMET	         ,   "minDeltaPhijetMET/F");
	_outTree->Branch("minDeltaPhijetMETZtau"		,  &minDeltaPhijetMETZtau	         ,   "minDeltaPhijetMETZtau/F");
	_outTree->Branch("jetPt_minDeltaPhijetMET"		,  &jetPt_minDeltaPhijetMET	         ,   "jetPt_minDeltaPhijetMET/F");
	
	_outTree->Branch("CovMartixfirst"		,  &CovMartixfirst	         ,   "CovMartixfirst/F");
	_outTree->Branch("CovMartixsecond"		,  &CovMartixsecond	         ,   "CovMartixsecond/F");
	_outTree->Branch("CovMatrixCorrelation"		,  &CovMatrixCorrelation	         ,   "CovMatrixCorrelation/F");
	_outTree->Branch("svmass"		,  &svmass	         ,   "svmass/F");
	_outTree->Branch("PzetaVis"		,  &PzetaVis	         ,   "PzetaVis/F");
	_outTree->Branch("PzetaMiss"		,  &PzetaMiss	         ,   "PzetaMiss/F");
	_outTree->Branch("Pzeta25"		,  &Pzeta25	         ,   "Pzeta25/F");
	_outTree->Branch("Pzeta45"		,  &Pzeta45	         ,   "Pzeta45/F");
	_outTree->Branch("Pzeta65"		,  &Pzeta65	         ,   "Pzeta65/F");
	_outTree->Branch("Pzeta85"		,  &Pzeta85	         ,   "Pzeta85/F");
	
	
	std::stringstream s;
	s << "triggerFlags[" << triggers.size() << "]/b";
	_outTree->Branch("triggerFlags", triggerFlags, s.str().c_str()); 
	
	_outTree->Branch("EVENT"		,  &EVENT	         ,   "run/I:lumi/I:event/I:json/I");
	_outTree->Branch("hbhe"		,  &hbhe	         ,   "hbhe/b");
	_outTree->Branch("totalKinematics"		,  &totalKinematics	         ,   "totalKinematics/b");
	_outTree->Branch("ecalFlag"		,  &ecalFlag	         ,   "ecalFlag/b");
	_outTree->Branch("cschaloFlag"		,  &cschaloFlag	         ,   "cschaloFlag/b");
	_outTree->Branch("hcallaserFlag"		,  &hcallaserFlag	         ,   "hcallaserFlag/b");
	_outTree->Branch("trackingfailureFlag"		,  &trackingfailureFlag	         ,   "trackingfailureFlag/b");
	_outTree->Branch("btag1TSF"		,  &btag1TSF	         ,   "btag1TSF/F");
	_outTree->Branch("btag2TSF"		,  &btag2TSF	         ,   "btag2TSF/F");
	_outTree->Branch("btag1T2CSF"	,  &btag1T2CSF	         ,   "btag1T2CSF/F");
	_outTree->Branch("btag2CSF"	,  &btag2CSF	         ,   "btag2CSF/F");
	_outTree->Branch("btagA0CSF"	,  &btagA0CSF	         ,   "btagA0CSF/F");
	_outTree->Branch("btagA0TSF"	,  &btagA0TSF	         ,   "btagA0TSF/F");
	_outTree->Branch("btag1TA1C"	,  &btag1TA1C	         ,   "btag1TA1C/F");
	
	_outTree->Branch("vLepton_decayModeFinding",vLeptons.decayModeFinding,"decayModeFinding[nvlep]/F");
	_outTree->Branch("vLepton_byLooseCombinedIsolationDeltaBetaCorr",vLeptons.byLooseCombinedIsolationDeltaBetaCorr,"byLooseCombinedIsolationDeltaBetaCorr[nvlep]/F");
	_outTree->Branch("vLepton_againstMuonTight",vLeptons.againstMuonTight,"againstMuonTight[nvlep]/F");
	_outTree->Branch("vLepton_againstElectronLoose",vLeptons.againstElectronLoose,"againstElectronLoose[nvlep]/F");
	_outTree->Branch("vLepton_againstElectronMedium",vLeptons.againstElectronMedium,"againstElectronMedium[nvlep]/F");
	_outTree->Branch("vLepton_againstElectronMVA",vLeptons.againstElectronMVA,"againstElectronMVA[nvlep]/F");
	_outTree->Branch("vLepton_NsignalPFChargedHadrCands",vLeptons.NsignalPFChargedHadrCands,"NsignalPFChargedHadrCands[nvlep]/I");
	_outTree->Branch("vLepton_NsignalPFGammaCands",vLeptons.NsignalPFGammaCands,"NsignalPFGammaCands[nvlep]/I");
	_outTree->Branch("vLepton_leadPFChargedHadrCandPt",vLeptons.leadPFChargedHadrCandPt,"leadPFChargedHadrCandPt[nvlep]/F");
	_outTree->Branch("vLepton_byLooseIsolation",vLeptons.byLooseIsolation,"byLooseIsolation[nvlep]/F");
	_outTree->Branch("vLepton_byMediumIsolation",vLeptons.byMediumIsolation,"byMediumIsolation[nvlep]/F");
	_outTree->Branch("vLepton_byTightIsolation",vLeptons.byTightIsolation,"byTightIsolation[nvlep]/F");
	
	_outTree->Branch("aLepton_decayModeFinding",aLeptons.decayModeFinding,"decayModeFinding[nvlep]/F");
	_outTree->Branch("aLepton_byLooseCombinedIsolationDeltaBetaCorr",aLeptons.byLooseCombinedIsolationDeltaBetaCorr,"byLooseCombinedIsolationDeltaBetaCorr[nvlep]/F");
	_outTree->Branch("aLepton_againstMuonTight",aLeptons.againstMuonTight,"againstMuonTight[nvlep]/F");
	_outTree->Branch("aLepton_againstElectronLoose",aLeptons.againstElectronLoose,"againstElectronLoose[nvlep]/F");
	_outTree->Branch("aLepton_againstElectronMedium",aLeptons.againstElectronMedium,"againstElectronMedium[nvlep]/F");
	_outTree->Branch("aLepton_againstElectronMVA",aLeptons.againstElectronMVA,"againstElectronMVA[nvlep]/F");
	_outTree->Branch("aLepton_NsignalPFChargedHadrCands",aLeptons.NsignalPFChargedHadrCands,"NsignalPFChargedHadrCands[nvlep]/I");
	_outTree->Branch("aLepton_NsignalPFGammaCands",aLeptons.NsignalPFGammaCands,"NsignalPFGammaCands[nvlep]/I");
	_outTree->Branch("aLepton_leadPFChargedHadrCandPt",aLeptons.leadPFChargedHadrCandPt,"leadPFChargedHadrCandPt[nvlep]/F");
	_outTree->Branch("aLepton_byLooseIsolation",aLeptons.byLooseIsolation,"byLooseIsolation[nvlep]/F");
	_outTree->Branch("aLepton_byMediumIsolation",aLeptons.byMediumIsolation,"byMediumIsolation[nvlep]/F");
	_outTree->Branch("aLepton_byTightIsolation",aLeptons.byTightIsolation,"byTightIsolation[nvlep]/F");
		
	int ievt=0;  
  	int Nevents = 0;
	int Nevents_Zmumu = 0, Nevents_Zee = 0, Nevents_Wmun = 0, Nevents_Wen = 0, Nevents_Znunu = 0, Nmulti = 0;
	int Nevents_Zemu = 0, Nevents_Ztaumu = 0, Nevents_Ztaue = 0, Nevents_Wtaun = 0, Nevents_Ztautau = 0, Nevents_Zbb = 0;
	int Nevents_LowJetPt = 0;
	//	unsigned int LumiNumber = -99;
	//	unsigned int RunNumber = 0;
	//	unsigned int trueLumiNumber = 999999;
	//	unsigned int trueRunNumber = 0;
	
	//  TFile* inFile = new TFile(inputFile.c_str(), "read");
	for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile) {
		std::cout << iFile << std::endl;
		TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
		if(inFile==0) continue;
		
		// loop the events
		
		fwlite::Event ev(inFile);
		for(ev.toBegin(); !ev.atEnd() ; ++ev, ++ievt)
			{
				if(Verbage) std::cout << "Event number  " << ievt << std::endl;
				
				if (ievt <= skipEvents_) continue;
				if (maxEvents_ >= 0){
					if (ievt > maxEvents_ + skipEvents_) break;
				};
				
				fwlite::Handle< VHbbEventAuxInfo > vhbbAuxHandle; 
				vhbbAuxHandle.getByLabel(ev,"HbbAnalyzerNew");
				const VHbbEventAuxInfo & aux = *vhbbAuxHandle.product();
				
				if(EVENT.run < runMin_ && runMin_ > 0) continue;
				if(EVENT.run > runMax_ && runMax_ > 0) continue;
				
				count->Fill(1.);
				
				
				
				
				
				
				PUweight=1.;
				PUweight2011B=1.;
				if(isMC_){
					
					// PU weights // Run2011A
					std::map<int, unsigned int>::const_iterator puit = aux.puInfo.pus.find(0);
					int npu =puit->second ;
					PUweight =  lumiWeights.weight( npu );        
					pu->Fill(puit->second);
					// PU weight Run2011B
					// PU weight Run2011B
					std::map<int, unsigned int>::const_iterator puit0 =  aux.puInfo.pus.find(0);
					std::map<int, unsigned int>::const_iterator puitm1 = aux.puInfo.pus.find(-1);
					std::map<int, unsigned int>::const_iterator puitp1 = aux.puInfo.pus.find(+1);
					PU0=puit0->second;
					PUp1=puitp1->second;
					PUm1=puitm1->second;
					input3DPU->Fill(PUm1,PU0,PUp1);	  
					PUweight2011B = lumiWeights2011B.weight3D( puitm1->second, puit0->second,puitp1->second); 
					
				}
				countWithPU->Fill(1,PUweight);
				countWithPU2011B->Fill(1,PUweight2011B);
				
				//Write event info 
				EVENT.run = ev.id().run();
				EVENT.lumi = ev.id().luminosityBlock();
				EVENT.event = ev.id().event();
				EVENT.json = jsonContainsEvent (jsonVector, ev);
				//FIXME : need to update EDM ntuple with BHadron infos
				// 	// simBHadrons
				// 	fwlite::Handle<SimBHadronCollection> SBHC;
				// 	SBHC.getByLabel(ev, "bhadrons");
				// 	const SimBHadronCollection sbhc = *(SBHC.product());
				
				const std::vector<VHbbCandidate> * candZ ;
				const std::vector<VHbbCandidate> * candW ;
				VHbbEvent modifiedEvent;;
				const VHbbEvent *  iEvent =0;
				if(fromCandidate)
					{
						fwlite::Handle< std::vector<VHbbCandidate> > vhbbCandHandleZ;
						vhbbCandHandleZ.getByLabel(ev,"hbbBestCSVPt20Candidates");
						candZ = vhbbCandHandleZ.product();
						
						fwlite::Handle< std::vector<VHbbCandidate> > vhbbCandHandle;
						vhbbCandHandle.getByLabel(ev,"hbbHighestPtHiggsPt30Candidates");
						candW = vhbbCandHandle.product();
					}
				else
					{
						candZlocal->clear();
						candWlocal->clear();
						fwlite::Handle< VHbbEvent > vhbbHandle; 
						vhbbHandle.getByLabel(ev,"HbbAnalyzerNew");
						modifiedEvent = *vhbbHandle.product();
						if(isMC_)
							{
								iEvent= &modifiedEvent;
								
								for(size_t j=0; j< modifiedEvent.simpleJets2.size() ; j++)
									{
										TLorentzVector & p4 = modifiedEvent.simpleJets2[j].p4; 
										TLorentzVector & mcp4 = modifiedEvent.simpleJets2[j].bestMCp4;
										if ((fabs(p4.Pt() - mcp4.Pt())/ p4.Pt())<0.5) { //Limit the effect to the core 
											float cor = (p4.Pt()+resolutionBias(fabs(p4.Eta()))*(p4.Pt()-mcp4.Pt()))/p4.Pt();
											p4.SetPtEtaPhiE(p4.Pt()*cor,p4.Eta(), p4.Phi(), p4.E()*cor);
										}
									}
							} else
								{
									iEvent = vhbbHandle.product();
								}  
						
						Nevents++;
						if (Verbage) std::cout << " event number " << Nevents << std::endl;		  
						
						algoZ->run(iEvent,*candZlocal);
						algoW->run(iEvent,*candWlocal);
						
						if(candZlocal->size() == 0 or candZlocal->at(0).H.jets.size() < 2)  //recover low pt 
							{
								candZlocal->clear();
								candWlocal->clear();
								algoRecoverLowPt->run(iEvent,*candZlocal);
								algoRecoverLowPt->run(iEvent,*candWlocal);
								if(candZlocal->size() != 0 && candZlocal->at(0).H.jets.size() > 1) {
									Nevents_LowJetPt++;
									if (Verbage) std::cout << "Event " << ievt<< " has a jet pt 20-15 " << endl;
								}				  
							}
						
						candZ= candZlocal; 
						candW= candWlocal; 
						/*     for(size_t m=0;m<iEvent->muInfo.size();m++)
						 { 
						 
						 if( fabs(iEvent->muInfo[m].p4.Pt()-28.118684) < 0.0001 ||
						 fabs(iEvent->muInfo[m].p4.Pt()-34.853199) < 0.0001 )  
						 {
						 std::cout << "FOUND " << iEvent->muInfo[m].p4.Pt() <<  " " << EVENT.event << " " << candW->size() << " " << candZ->size() << std::endl;
						 }
						 }
						 */
						
					}
				
				const std::vector<VHbbCandidate> * cand = candZ;
				
				
				
				/*  fwlite::Handle< VHbbEvent > vhbbHandle; 
				 vhbbHandle.getByLabel(ev,"HbbAnalyzerNew");
				 const VHbbEvent iEvent = *vhbbHandle.product();
				 */
				
				//      std::clog << "Filling tree "<< std::endl;
				bool isW=false;
				
				// to check how much we gain with jets subtraction 
				
				
				
				genHpt=aux.mcH.size() > 0 ? aux.mcH[0].p4.Pt():-99;
				
				if(cand->size() == 0 or cand->at(0).H.jets.size() < 2) continue;
				//if(cand->size() == 0 ) continue;
				if(Verbage)std::cout << "cand->size() " << cand->size() << std::endl;
				//std::cout << "cand->at(0).H.jets.size() " << cand->at(0).H.jets.size() << std::endl;
				
				unsigned int tauCandidate = 0;
				if(cand->size() > 1 ) 
					{
						std::cout << "MULTIPLE CANDIDATES: " << cand->size() << std::endl;
						Nmulti++;
						for (unsigned int i=0; i!=cand->size(); i++){
							if (cand->at(i).candidateType > 4 && (tauCandidate!=7||tauCandidate!=6)) tauCandidate = i;
							cout << "cand(" << i << ").Vtype: " << cand->at(i).candidateType << endl; 
						}// loop on candidates
					}
				if(cand->at(0).candidateType == VHbbCandidate::Wmun || cand->at(0).candidateType == VHbbCandidate::Wen ) { cand=candW; isW=true; }
				if((cand->at(0).candidateType == VHbbCandidate::Wmun || cand->at(0).candidateType == VHbbCandidate::Wen || cand->at(0).candidateType == VHbbCandidate::Wtaun) && cand->size() < 2) { cand=candW; isW=true;}
 
				if(cand->at(0).H.jets.size() < 2) continue;
				
				// secondary vtxs
				fwlite::Handle<std::vector<reco::Vertex> > SVC;
				SVC.getByLabel(ev,"bcandidates");
				const std::vector<reco::Vertex> svc = *(SVC.product());
				
				//CovMatrix
				TMatrixD metCov(2, 2);
				fwlite::Handle<PFMEtSignCovMatrix> CovM;
				CovM.getByLabel(ev,"pfMEtSignCovMatrix");
				metCov = (*CovM);
				if(Verbage) cout << "first element " << metCov(0,0) << endl;
				CovMartixfirst =  metCov(0,0);
				CovMartixsecond =  metCov(1,1);
				CovMatrixCorrelation =  metCov(0,1);
				
				
				const VHbbCandidate & vhCand =  cand->at(0);
				patFilters.setEvent(&ev,"VH");
				hbhe = patFilters.accept("hbhe");
				ecalFlag = patFilters.accept("ecalFilter");
				totalKinematics = patFilters.accept("totalKinematics");
				cschaloFlag = patFilters.accept("cschaloFilter");   
				hcallaserFlag = patFilters.accept("hcallaserFilter");   
				trackingfailureFlag = patFilters.accept("trackingfailureFilter");   
				
				trigger.setEvent(&ev);
				for(size_t j=0;j < triggers.size();j++)
					triggerFlags[j]=trigger.accept(triggers[j]);
				
				eventFlav=0;
				
				
				
				if(aux.mcBbar.size() > 0 || aux.mcB.size() > 0) eventFlav=5;
				else if(aux.mcC.size() > 0) eventFlav=4;
				
				
				Vtype = vhCand.candidateType;
				
				
				/*nRecoTaus = iEvent->tauInfo.size();
				 for(size_t j=0;j< iEvent->tauInfo.size();j++){
				 tau_decayModeFinding[j] = iEvent->tauInfo[j].decayModeFinding;
				 tau_byLooseCombinedIsolationDeltaBetaCorr[j] = iEvent->tauInfo[j].byLooseCombinedIsolationDeltaBetaCorr;
				 tau_againstMuonTight[j] = iEvent->tauInfo[j].againstMuonTight;
				 tau_againstElectronLoose[j] = iEvent->tauInfo[j].againstElectronLoose;
				 tau_againstElectronMedium[j] = iEvent->tauInfo[j].againstElectronMedium;
				 tau_againstElectronMVA[j] = iEvent->tauInfo[j].againstElectronMVA;
				 }// loop on taus
				 */
				nRecoTaus = vhCand.V.taus.size();
				for(size_t j=0;j< vhCand.V.taus.size();j++){
					tau_decayModeFinding[j] = vhCand.V.taus[j].decayModeFinding;
					tau_byLooseCombinedIsolationDeltaBetaCorr[j] = vhCand.V.taus[j].byLooseCombinedIsolationDeltaBetaCorr;
					tau_againstMuonTight[j] = vhCand.V.taus[j].againstMuonTight;
					tau_againstElectronLoose[j] = vhCand.V.taus[j].againstElectronLoose;
					tau_againstElectronMedium[j] = vhCand.V.taus[j].againstElectronMedium;
					tau_againstElectronMVA[j] = vhCand.V.taus[j].againstElectronMVA;
				}// loop on taus
				
				
				H.mass = vhCand.H.p4.M();
				H.pt = vhCand.H.p4.Pt();
				
				H.eta = vhCand.H.p4.Eta();
				H.phi = vhCand.H.p4.Phi();
				
				
				FatH.FatHiggsFlag = vhCand.FatH.FatHiggsFlag;
				if(vhCand.FatH.FatHiggsFlag){ FatH.mass= vhCand.FatH.p4.M(); 
					FatH.pt = vhCand.FatH.p4.Pt();
					FatH.eta = vhCand.FatH.p4.Eta();
					FatH.phi = vhCand.FatH.p4.Phi();
					
					//        if(vhCand.FatH.FatHiggsFlag)  vhCand.FatH.subjetsSize; 
					nfathFilterJets=vhCand.FatH.subjetsSize;  
					for( int j=0; j < nfathFilterJets; j++ ){
						fathFilterJets.set(vhCand.FatH.jets[j],j);
					}
					
					if(nfathFilterJets==2){
						FatH.filteredmass=(vhCand.FatH.jets[0].p4+vhCand.FatH.jets[1].p4).M();
					FatH.filteredpt=(vhCand.FatH.jets[0].p4+vhCand.FatH.jets[1].p4).Pt();}
					else if(nfathFilterJets==3){
						FatH.filteredmass=(vhCand.FatH.jets[0].p4+vhCand.FatH.jets[1].p4+vhCand.FatH.jets[2].p4).M();
					FatH.filteredpt=(vhCand.FatH.jets[0].p4+vhCand.FatH.jets[1].p4+vhCand.FatH.jets[2].p4).Pt();}
					
				}
				
				
				V.mass = vhCand.V.p4.M();
				if(isW) V.mass = vhCand.Mt(); 
				V.pt = vhCand.V.p4.Pt();
				V.eta = vhCand.V.p4.Eta();
				V.phi = vhCand.V.p4.Phi();
				nhJets=2;
				hJets.set(vhCand.H.jets[0],0);
				hJets.set(vhCand.H.jets[1],1);
				
				
				aJets.reset();
				naJets=vhCand.additionalJets.size();
				numBJets=0;
				if(vhCand.H.jets[0].csv> btagThr) numBJets++;
				if(vhCand.H.jets[1].csv> btagThr) numBJets++;
				for( int j=0; j < naJets && j < MAXJ; j++ ) 
					{
						aJets.set(vhCand.additionalJets[j],j);
						if(vhCand.additionalJets[j].csv> btagThr) numBJets++;
						if (Vtype==VHbbCandidate::Wtaun||Vtype==VHbbCandidate::Ztaumu||Vtype==VHbbCandidate::Ztaue||Vtype==VHbbCandidate::Ztautau) {
							aJets.selectedTauDR[j] = vhCand.V.taus[0].p4.DeltaR(vhCand.additionalJets[j].p4);
						}						
					}   
				numJets = vhCand.additionalJets.size()+2;
				H.dR = deltaR(vhCand.H.jets[0].p4.Eta(),vhCand.H.jets[0].p4.Phi(),vhCand.H.jets[1].p4.Eta(),vhCand.H.jets[1].p4.Phi());
				H.dPhi = deltaPhi(vhCand.H.jets[0].p4.Phi(),vhCand.H.jets[1].p4.Phi());
				H.dEta= TMath::Abs( vhCand.H.jets[0].p4.Eta() - vhCand.H.jets[1].p4.Eta() );
				HVdPhi = fabs( deltaPhi(vhCand.H.p4.Phi(),vhCand.V.p4.Phi()) ) ;
				HVMass = (vhCand.H.p4 + vhCand.V.p4).M() ;
				HMETdPhi = fabs( deltaPhi(vhCand.H.p4.Phi(),vhCand.V.mets.at(0).p4.Phi()) ) ;
				VMt = vhCand.Mt() ;
				VMte = vhCand.Mte();
				VMtmu = vhCand.Mtmu();
				//eltaPullAngle = vhCand.H.deltaTheta;
				deltaPullAngle  =  VHbbCandidateTools::getDeltaTheta(vhCand.H.jets[0],vhCand.H.jets[1]);
				deltaPullAngle2 =  VHbbCandidateTools::getDeltaTheta(vhCand.H.jets[1],vhCand.H.jets[0]);
				
				if(Verbage) cout << "at(0) used for HMETdPhi " << HMETdPhi << endl;

				
				hJets.cosTheta[0]=  vhCand.H.helicities[0];
				hJets.cosTheta[1]=  vhCand.H.helicities[1];
				// METInfo calomet;  METInfo tcmet;  METInfo pfmet;  METInfo mht;  METInfo metNoPU
				MET.et = vhCand.V.mets.at(0).p4.Pt();
				MET.phi = vhCand.V.mets.at(0).p4.Phi();
				MET.sumet = vhCand.V.mets.at(0).sumEt;
				MET.sig = vhCand.V.mets.at(0).metSig;
				
				
				if(Verbage) cout << "middle of MET info" << endl;

				
				fakeMET.sumet = 0;
				fakeMET.sig = 0;
				fakeMET.et = 0;
				fakeMET.phi = 0;
				if( Vtype == VHbbCandidate::Zmumu) {
					TVector3 mu1 = vhCand.V.muons[0].p4.Vect();
					TVector3 mu2 = vhCand.V.muons[1].p4.Vect();
					// Not needed with PFMET
					//		mu1.SetMag( mu1.Mag() - vhCand.V.muons[0].emEnergy - vhCand.V.muons[0].hadEnergy);
					//		mu2.SetMag( mu2.Mag() - vhCand.V.muons[1].emEnergy - vhCand.V.muons[1].hadEnergy);
					TVector3 sum = vhCand.V.mets.at(0).p4.Vect() + mu1 + mu2;
					fakeMET.et = sum.Pt();
					fakeMET.phi = sum.Phi();
					fakeMET.sumet = vhCand.V.mets.at(0).sumEt - mu1.Pt() - mu2.Pt();
				}
				
				
				
				METnoPU.et = iEvent->metNoPU.p4.Pt();
				METnoPU.phi = iEvent->metNoPU.p4.Phi();
				METnoPU.sumet = iEvent->metNoPU.sumEt;
				METnoPU.sig = iEvent->metNoPU.metSig;
				METnoPUCh.et = iEvent->metCh.p4.Pt();
				METnoPUCh.phi = iEvent->metCh.p4.Phi();
				METnoPUCh.sumet = iEvent->metCh.sumEt;
				METnoPUCh.sig = iEvent->metCh.metSig;
				
				METnoPUCh.et = iEvent->metCh.p4.Pt();
				METnoPUCh.phi = iEvent->metCh.p4.Phi();
				METnoPUCh.sumet = iEvent->metCh.sumEt;
				METnoPUCh.sig = iEvent->metCh.metSig;
				
				
				METtype1corr.et = iEvent->pfmetType1corr.p4.Pt();
				METtype1corr.phi = iEvent->pfmetType1corr.p4.Phi();
				METtype1corr.sumet = iEvent->pfmetType1corr.sumEt;
				METtype1corr.sig = iEvent->pfmetType1corr.metSig;
				
				
				METtype1p2corr.et = iEvent->pfmetType1p2corr.p4.Pt();
				METtype1p2corr.phi = iEvent->pfmetType1p2corr.p4.Phi();
				METtype1p2corr.sumet = iEvent->pfmetType1p2corr.sumEt;
				METtype1p2corr.sig = iEvent->pfmetType1p2corr.metSig;
				
				
				METnoPUtype1corr.et = iEvent->pfmetNoPUType1corr.p4.Pt();
				METnoPUtype1corr.phi = iEvent->pfmetNoPUType1corr.p4.Phi();
				METnoPUtype1corr.sumet = iEvent->pfmetNoPUType1corr.sumEt;
				METnoPUtype1corr.sig = iEvent->pfmetNoPUType1corr.metSig;
				
				
				METnoPUtype1p2corr.et = iEvent->pfmetNoPUType1p2corr.p4.Pt();
				METnoPUtype1p2corr.phi = iEvent->pfmetNoPUType1p2corr.p4.Phi();
				METnoPUtype1p2corr.sumet = iEvent->pfmetNoPUType1p2corr.sumEt;
				METnoPUtype1p2corr.sig = iEvent->pfmetNoPUType1p2corr.metSig;
				
				if(Verbage)	std::cout << " iEvent->metUncInfo.size() " << iEvent->metUncInfo.size() << std::endl;
				for(size_t m=0;m<iEvent->metUncInfo.size();m++)
					{ 
						metUnc.set(iEvent->metUncInfo[m], m );
						if(Verbage)	     std::cout << "metUncInfo[" << m <<" ].et = " << metUnc.et[m] << std::endl; 
					}
				
				
				rho = aux.puInfo.rho;
				rho25 = aux.puInfo.rho25;
				nPVs=aux.pvInfo.nVertices; 
				
				if(!fromCandidate) {
					MHT.mht = iEvent->mht.p4.Pt(); 
					MHT.phi = iEvent->mht.p4.Phi(); 
					MHT.ht = iEvent->mht.sumEt; 
					MHT.sig = iEvent->mht.metSig; 
				}
				
				/////////
				// track sharing flags:
				////////
				TkSharing.HiggsCSVtkSharing = TkSharing.HiggsIPtkSharing = TkSharing.HiggsSVtkSharing = TkSharing.FatHiggsCSVtkSharing = TkSharing.FatHiggsIPtkSharing = TkSharing.FatHiggsSVtkSharing = false;
				
				// csv tracks
				if (vhCand.H.jets[0].csvNTracks > 0 && vhCand.H.jets[1].csvNTracks > 0){
					for (int t=0;t!=vhCand.H.jets[0].csvNTracks;t++){
						for (int ti=0;ti!=vhCand.H.jets[1].csvNTracks;ti++){
							if ((int)vhCand.H.jets[0].csvTrackIds[t] == (int)vhCand.H.jets[1].csvTrackIds[ti]){
								TkSharing.HiggsCSVtkSharing = true;
							}// same trackID
						}// loop tracks in second hjet
					}// loop tracks in first hjet
				}// if tracks in jet
				
				// ip tracks
				if (vhCand.H.jets[0].btagNTracks > 0 && vhCand.H.jets[1].btagNTracks > 0){
					for (int t=0;t!=vhCand.H.jets[0].btagNTracks;t++){
						for (int ti=0;ti!=vhCand.H.jets[1].btagNTracks;ti++){
							if ((int)vhCand.H.jets[0].btagTrackIds[t] == (int)vhCand.H.jets[1].btagTrackIds[ti]){
								TkSharing.HiggsIPtkSharing = true;
							}// same trackID
						}// loop tracks in second hjet
					}// loop tracks in first hjet
				}// if tracks in jet
				
				// sv tracks
				if (vhCand.H.jets[0].vtxNTracks > 0 && vhCand.H.jets[1].vtxNTracks > 0){
					for (int t=0;t!=vhCand.H.jets[0].vtxNTracks;t++){
						for (int ti=0;ti!=vhCand.H.jets[1].vtxNTracks;ti++){
							if ((int)vhCand.H.jets[0].vtxTrackIds[t] == (int)vhCand.H.jets[1].vtxTrackIds[ti]){
								TkSharing.HiggsSVtkSharing = true;
							}// same trackID
						}// loop tracks in second hjet
					}// loop tracks in first hjet
				}// if tracks in jet
				
				// tracksharing for Filtered jets:
				if(vhCand.FatH.FatHiggsFlag && nfathFilterJets > 1){
					
					// csv tracks
					if (vhCand.FatH.jets[0].csvNTracks > 0 && vhCand.FatH.jets[1].csvNTracks > 0){
						for (int t=0;t!=vhCand.FatH.jets[0].csvNTracks;t++){
							for (int ti=0;ti!=vhCand.FatH.jets[1].csvNTracks;ti++){
								if ((int)vhCand.FatH.jets[0].csvTrackIds[t] == (int)vhCand.FatH.jets[1].csvTrackIds[ti]){
									TkSharing.FatHiggsCSVtkSharing = true;
								}// same trackID
							}// loop tracks in second hjet
						}// loop tracks in first hjet
					}// if tracks in jet
					
					// ip tracks
					if (vhCand.FatH.jets[0].btagNTracks > 0 && vhCand.FatH.jets[1].btagNTracks > 0){
						for (int t=0;t!=vhCand.FatH.jets[0].btagNTracks;t++){
							for (int ti=0;ti!=vhCand.FatH.jets[1].btagNTracks;ti++){
								if ((int)vhCand.FatH.jets[0].btagTrackIds[t] == (int)vhCand.FatH.jets[1].btagTrackIds[ti]){
									TkSharing.FatHiggsIPtkSharing = true;
								}// same trackID
							}// loop tracks in second hjet
						}// loop tracks in first hjet
					}// if tracks in jet
					
					// sv tracks
					if (vhCand.FatH.jets[0].vtxNTracks > 0 && vhCand.FatH.jets[1].vtxNTracks > 0){
						for (int t=0;t!=vhCand.FatH.jets[0].vtxNTracks;t++){
							for (int ti=0;ti!=vhCand.FatH.jets[1].vtxNTracks;ti++){
								if ((int)vhCand.FatH.jets[0].vtxTrackIds[t] == (int)vhCand.FatH.jets[1].vtxTrackIds[ti]){
									TkSharing.FatHiggsSVtkSharing = true;
								}// same trackID
							}// loop tracks in second hjet
						}// loop tracks in first hjet
					}// if tracks in jet
					
				}// fatH
				
				////////
				////////
				
				
				
				//Secondary Vertices
				IVF.reset();
				nSvs = svc.size();
				const TVector3 recoPv = aux.pvInfo.firstPVInPT2;
				const math::XYZPoint myPv(recoPv);
				//look here for Matrix filling info http://project-mathlibs.web.cern.ch/project-mathlibs/sw/html/SMatrixDoc.html
				std::vector<double> fillMatrix(6);
				for (int i = 0; i<6; ++i) fillMatrix[i] = 0.;
				fillMatrix[0] = TMath::Power(0.002,2);
				fillMatrix[2] = TMath::Power(0.002,2);
				fillMatrix[5] = TMath::Power(0.002,2);
				const ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > myFakeMatrixError(fillMatrix.begin(),fillMatrix.end());
				const reco::Vertex recoVtxPv(myPv, myFakeMatrixError);
				for( int j=0; j < nSvs && j < MAXB; ++j ) {
					const GlobalVector flightDir = flightDirection(recoPv,svc[j]);
					reco::SecondaryVertex recoSv(recoVtxPv, svc[j], flightDir ,true);
					IVF.set( recoSv, recoPv ,j);
				}
				if(nSvs > 2){
					TLorentzVector BCands_H1, BCands_H2, BCands_H;
					BCands_H1.SetPtEtaPhiM(IVF.pt[0], IVF.eta[0], IVF.phi[0], IVF.massBcand[0]);
					BCands_H2.SetPtEtaPhiM(IVF.pt[1], IVF.eta[1], IVF.phi[1], IVF.massBcand[1]);
					BCands_H = BCands_H1 + BCands_H2;
					SVH.dR = deltaR(IVF.eta[0], IVF.phi[0], IVF.eta[1], IVF.phi[1] );
					SVH.dPhi = deltaPhi(IVF.phi[0], IVF.phi[1] );
					SVH.dEta = TMath::Abs(IVF.eta[0] - IVF.eta[1] );
					SVH.mass = BCands_H.M();
					SVH.pt = BCands_H.Pt();
					SVH.eta = BCands_H.Eta();
					SVH.phi = BCands_H.Phi();
				}
				
				//Loop on jets
				double maxBtag=-99999;
				minDeltaPhijetMET = 999;
				minDeltaPhijetMETZtau = 999;
				svmass = 999;
				PzetaVis = 999, PzetaMiss = 999, Pzeta25 = 999, Pzeta45 = 999, Pzeta65 = 999, Pzeta85 = 999;
				TLorentzVector bJet;
				TLorentzVector METZtau;
				std::vector<std::vector<BTagWeight::JetInfo> > btagJetInfos;
				std::vector<float> jet10eta;
				std::vector<float> jet10pt;
				std::vector<float> jet30eta;
				std::vector<float> jet30pt;
				if(fromCandidate)
					{
						//Loop on Higgs Jets
						for(unsigned int j=0; j < vhCand.H.jets.size(); j++ ){
							if (vhCand.H.jets[j].csv > maxBtag) { bJet=vhCand.H.jets[j].p4 ; maxBtag =vhCand.H.jets[j].csv; }
							if (fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), vhCand.H.jets[j].p4.Phi())) < minDeltaPhijetMET) 
								{
									minDeltaPhijetMET=fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), vhCand.H.jets[j].p4.Phi())); 
									METZtau = vhCand.V.mets.at(0).p4 + vhCand.V.p4;
									minDeltaPhijetMETZtau=fabs(deltaPhi( METZtau.Phi(), vhCand.H.jets[j].p4.Phi()));
									jetPt_minDeltaPhijetMET=vhCand.H.jets[j].p4.Pt();
								}
							btagJetInfos.push_back(btagEff.jetInfo(vhCand.H.jets[j]));
						}
						//Loop on Additional Jets
						for(unsigned int j=0; j < vhCand.additionalJets.size(); j++ ){
							if (vhCand.additionalJets[j].csv > maxBtag) { bJet=vhCand.additionalJets[j].p4 ; maxBtag =vhCand.additionalJets[j].csv; }
							/*	      if (fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), vhCand.additionalJets[j].p4.Phi())) < minDeltaPhijetMET) 
							 {
							 minDeltaPhijetMET=fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), vhCand.additionalJets[j].p4.Phi()));
							 jetPt_minDeltaPhijetMET=vhCand.additionalJets[j].p4.Pt();
							 }*/
							if( ( isW && ! useHighestPtHiggsW ) ||  ( ! isW && ! useHighestPtHiggsZ )  )  // btag SF computed using only H-jets if best-H made with dijetPt rather than best CSV
								{
									if(vhCand.additionalJets[j].p4.Pt() > 20)
										btagJetInfos.push_back(btagEff.jetInfo(vhCand.additionalJets[j]));
								}
						}
					} 
				else
					{
						//Loop on all jets
						for(unsigned int j=0; j < iEvent->simpleJets2.size(); j++ ){
							if (iEvent->simpleJets2[j].csv > maxBtag) { bJet=iEvent->simpleJets2[j].p4 ; maxBtag =iEvent->simpleJets2[j].csv; }
							if ( iEvent->simpleJets2[j].p4.Pt() > 20 &&  fabs(iEvent->simpleJets2[j].p4.Eta()) < 2.5&& fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), iEvent->simpleJets2[j].p4.Phi())) < minDeltaPhijetMET)
								{
									minDeltaPhijetMET=fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), iEvent->simpleJets2[j].p4.Phi()));
									minDeltaPhijetMETZtau=fabs(deltaPhi( METZtau.Phi(), vhCand.H.jets[j].p4.Phi()));
									jetPt_minDeltaPhijetMET=iEvent->simpleJets2[j].p4.Pt();
								}
							if(iEvent->simpleJets2[j].p4.Pt() > 10)
								{
									jet10eta.push_back(iEvent->simpleJets2[j].p4.Eta());
									jet10pt.push_back(iEvent->simpleJets2[j].p4.Pt());
								}
							if(iEvent->simpleJets2[j].p4.Pt() > 30)
								{
									jet30eta.push_back(iEvent->simpleJets2[j].p4.Eta());
									jet30pt.push_back(iEvent->simpleJets2[j].p4.Pt());
								}
							
							//For events made with highest CSV, all jets in the event should be taken into account for "tagging" SF (anti tagging is a mess)
							// because for example a light jet not used for the Higgs can have in reality a higher CSV due to SF > 1 and become a higgs jet
							if( ( isW && ! useHighestPtHiggsW ) ||  ( ! isW && ! useHighestPtHiggsZ )  ) 
								{ 
									if(iEvent->simpleJets2[j].p4.Pt() > 20 && fabs(iEvent->simpleJets2[j].p4.Eta()) < 2.5)
										btagJetInfos.push_back(btagEff.jetInfo(iEvent->simpleJets2[j]));
								}
						}
						
						//Loop on Higgs jets
						
						for(unsigned int j=0; j < vhCand.H.jets.size(); j++ ) {
							
							//if we use the highest pt pair, only the two higgs jet should be used to compute the SF because the other jets are excluded 
							// by a criteria (pt of the dijet) that is not btag SF dependent 
							if(!( ( isW && ! useHighestPtHiggsW ) ||  ( ! isW && ! useHighestPtHiggsZ ) )) {
								btagJetInfos.push_back(btagEff.jetInfo(vhCand.H.jets[j]));
							}
						}
						
					}
				vLeptons.reset();
				weightTrig = 1.; // better to default to 1 
				weightTrigMay = -1.;
				weightTrigV4 = -1.; 
				weightTrigOrMu30 = 1.;
				TLorentzVector leptonForTop;
				size_t firstAddMu=0;
				size_t firstAddEle=0;
				size_t firstAddTau=0;
				if(Vtype == VHbbCandidate::Zmumu ){
				Nevents_Zmumu++;
					vLeptons.set(vhCand.V.muons[0],0,13); 
					vLeptons.set(vhCand.V.muons[1],1,13);
					float weightMuonRecoAndID = triggerWeight.scaleMuID(vLeptons.pt[0],vLeptons.eta[0]) * triggerWeight.scaleMuID(vLeptons.pt[1],vLeptons.eta[1]) ;
					float weightTrig1 = triggerWeight.scaleMuIsoHLT(vLeptons.pt[0],vLeptons.eta[0]);
					float weightTrig2 = triggerWeight.scaleMuIsoHLT(vLeptons.pt[1],vLeptons.eta[1]);
					float cweightTrig = weightTrig1 + weightTrig2 - weightTrig1*weightTrig2;
					weightTrig = weightMuonRecoAndID * cweightTrig;
					nvlep=2;
					firstAddMu=2;
				}
				if( Vtype == VHbbCandidate::Zee ){
				Nevents_Zee++;
					vLeptons.set(vhCand.V.electrons[0],0,11);
					vLeptons.set(vhCand.V.electrons[1],1,11);
					nvlep=2;
					firstAddEle=2;
					std::vector<float> pt,eta;
					pt.push_back(vLeptons.pt[0]); eta.push_back(vLeptons.eta[0]);
					pt.push_back(vLeptons.pt[1]); eta.push_back(vLeptons.eta[1]);
					weightEleRecoAndId=triggerWeight.scaleID95Ele(vLeptons.pt[0],vLeptons.eta[0]) * triggerWeight.scaleRecoEle(vLeptons.pt[0],vLeptons.eta[0]) *
					triggerWeight.scaleID95Ele(vLeptons.pt[1],vLeptons.eta[1]) * triggerWeight.scaleRecoEle(vLeptons.pt[1],vLeptons.eta[1]);
					weightEleTrigElePart = triggerWeight.scaleDoubleEle17Ele8(pt,eta); 
					weightEleTrigEleAugPart = triggerWeight.scaleDoubleEle17Ele8Aug(pt,eta); 
					weightTrig = (weightEleTrigElePart*1.14+weightEleTrigEleAugPart*0.98 )/2.12 * weightEleRecoAndId;
				}
				
				
				TMatrixD* metsig = new TMatrixD(2,2);
				(*metsig)[0][0] = metCov(0,0);
				(*metsig)[0][1] = metCov(0,1);
				(*metsig)[1][0] = metCov(1,0);
				(*metsig)[1][1] = metCov(1,1);
				TVector3 m,e,tau,metv;
				
				if(Vtype == VHbbCandidate::Zemu ){
				Nevents_Zemu++;
					leptonForTop=vhCand.V.muons[0].p4;
					vLeptons.set(vhCand.V.electrons[0],0,11); 
					vLeptons.set(vhCand.V.muons[0],1,13);
					float weightMuonHLTRecoAndID = 1.0;
					float weightEleHLTRecoAndId = 1.0;
					//cout <<"electron pt " << vLeptons.pt[0] << " eta " << vLeptons.eta[0] << endl;
					
					if (vLeptons.pt[1]>10&&vLeptons.pt[0]>10){
						weightEleHLTRecoAndId = triggerWeight.scaleEleCrossTrig(vLeptons.pt[0],vLeptons.eta[0]);
						weightMuonHLTRecoAndID =triggerWeight.scaleMuCrossTrig(vLeptons.pt[1],vLeptons.eta[1]);
						//Which trigger more likely to fire if both leptons above upper leg threshold
						if( vLeptons.pt[1]>20&&vLeptons.pt[0]>20){
							if (weightEleHLTRecoAndId>weightMuonHLTRecoAndID)weightMuonHLTRecoAndID = 1.0;
							if (weightEleHLTRecoAndId<weightMuonHLTRecoAndID)weightEleHLTRecoAndId = 1.0;
						}
						if( vLeptons.pt[1]>10&&vLeptons.pt[1]<20){ // Muon bottom leg
							weightEleHLTRecoAndId = 1.0;
						}
						if( vLeptons.pt[0]>10&&vLeptons.pt[0]<20){ // Electron bottom leg
							weightMuonHLTRecoAndID = 1.0;
						}
					}
					
					if (vLeptons.pt[1]>20&&vLeptons.pt[0]>20){
						weightEleHLTRecoAndId = weightEleHLTRecoAndId*triggerWeight.scaleID95Ele(vLeptons.pt[0],vLeptons.eta[0]) * triggerWeight.scaleRecoEle(vLeptons.pt[0],vLeptons.eta[0]);
						weightMuonHLTRecoAndID = weightMuonHLTRecoAndID*triggerWeight.scaleMuID(vLeptons.pt[1],vLeptons.eta[1]);
					} else {
						float MuonID = triggerWeight.scaleMuID(21,vLeptons.eta[1]);
						float EleID = triggerWeight.scaleID95Ele(21,vLeptons.eta[0]) * triggerWeight.scaleRecoEle(21,vLeptons.eta[0]);
						if ( vLeptons.pt[1] <= 15 && vLeptons.pt[1] >10) {
							if (vLeptons.eta[1] >= 1.0) weightMuonHLTRecoAndID = 0.98*weightMuonHLTRecoAndID*MuonID;
							if (vLeptons.eta[1] < 1.0) weightMuonHLTRecoAndID = 0.92*weightMuonHLTRecoAndID*MuonID;
						}
						if ( vLeptons.pt[1] <= 20 && vLeptons.pt[1] >15) {
							if (vLeptons.eta[1] >= 1.0) weightMuonHLTRecoAndID = 0.962*weightMuonHLTRecoAndID*MuonID;
							if (vLeptons.eta[1] < 1.0) weightMuonHLTRecoAndID = 0.948*weightMuonHLTRecoAndID*MuonID;
						}
						if ( vLeptons.pt[0] <= 15 && vLeptons.pt[0] >10) {
							if (vLeptons.eta[0] >= 1.2) weightEleHLTRecoAndId = 1.19*weightEleHLTRecoAndId*EleID;
							if (vLeptons.eta[0] < 1.2) weightEleHLTRecoAndId = 1.11*weightEleHLTRecoAndId*EleID;
						}
						if ( vLeptons.pt[0] <= 20 && vLeptons.pt[0] >15) {
							if (vLeptons.eta[0] >= 1.2) weightEleHLTRecoAndId = 1.06*weightEleHLTRecoAndId*EleID;
							if (vLeptons.eta[0] < 1.2) weightEleHLTRecoAndId = 0.98*weightEleHLTRecoAndId*EleID;
						}
					}					
					
					weightTrig = weightEleHLTRecoAndId * weightMuonHLTRecoAndID;
					nvlep=2;
					firstAddMu=1;
					firstAddEle=1;
					if (Verbage) std::cout <<"event " <<ievt << " is a Zemu event" << endl;
					
					LorentzVector tauelectron = LorentzVector(vhCand.V.electrons[0].p4.Px(),vhCand.V.electrons[0].p4.Py(),vhCand.V.electrons[0].p4.Pz(),vhCand.V.electrons[0].p4.P());
					LorentzVector taumuon = LorentzVector(vhCand.V.muons[0].p4.Px(),vhCand.V.muons[0].p4.Py(),vhCand.V.muons[0].p4.Pz(),vhCand.V.muons[0].p4.P());
					NSVfitStandalone::Vector measuredMET( vhCand.V.mets.at(0).p4.Px(), vhCand.V.mets.at(0).p4.Py(), 0);
					std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
					NSVfitStandalone::LorentzVector p1 = tauelectron;
					measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay,p1));
					NSVfitStandalone::LorentzVector p2 = taumuon;
					measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay,p2));
					unsigned int flagSA = 0;
					NSVfitStandaloneAlgorithm algo(measuredTauLeptons,measuredMET,*metsig,flagSA);
					algo.addLogM(false);
					algo.integrate();
					svmass = algo.getMass(); // mass uncertainty not implemented yet
					//cout << "original mass " << vhCand.V.p4.M() << " SVfit mass " << svmass << " %diff " << (vhCand.V.p4.M()-svmass)/vhCand.V.p4.M() << endl;

					m.SetPtEtaPhi(vhCand.V.muons[0].p4.Pt(),0,vhCand.V.muons[0].p4.Phi());
					e.SetPtEtaPhi(vhCand.V.electrons[0].p4.Pt(),0,vhCand.V.electrons[0].p4.Phi());
					metv.SetPtEtaPhi(vhCand.V.mets.at(0).p4.Pt(),0,vhCand.V.mets.at(0).p4.Phi());
					TVector3 bisector(m.Unit() + e.Unit());
					bisector = bisector.Unit();
					PzetaVis  = (m+e).Dot(bisector);
					PzetaMiss  =  metv.Dot(bisector);
					
					
				}
				
				if(Vtype == VHbbCandidate::Ztaumu ){
				Nevents_Ztaumu++;
					//cout << "Ztaumu candidate" << endl;
					leptonForTop=vhCand.V.muons[0].p4;
					vLeptons.set(vhCand.V.muons[0],0,13); 
					vLeptons.set(vhCand.V.taus[0],1,15);
					nvlep=2;
					firstAddMu=1;
					firstAddTau=1;
				
					LorentzVector hadronictau = LorentzVector(vhCand.V.taus[0].p4.Px(),vhCand.V.taus[0].p4.Py(),vhCand.V.taus[0].p4.Pz(),vhCand.V.taus[0].p4.P());
					LorentzVector taumuon = LorentzVector(vhCand.V.muons[0].p4.Px(),vhCand.V.muons[0].p4.Py(),vhCand.V.muons[0].p4.Pz(),vhCand.V.muons[0].p4.P());
					NSVfitStandalone::Vector measuredMET( vhCand.V.mets.at(0).p4.Px(), vhCand.V.mets.at(0).p4.Py(), 0);
					std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
					NSVfitStandalone::LorentzVector p1 = hadronictau;
					measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay,p1));
					NSVfitStandalone::LorentzVector p2 = taumuon;
					measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay,p2));
					unsigned int flagSA = 0;
					NSVfitStandaloneAlgorithm algo(measuredTauLeptons,measuredMET,*metsig,flagSA);
					algo.addLogM(false);
					algo.integrate();
					svmass = algo.getMass(); // mass uncertainty not implemented yet
					//cout << "original mass " << vhCand.V.p4.M() << " SVfit mass " << svmass << " %diff " << (vhCand.V.p4.M()-svmass)/vhCand.V.p4.M() << endl;
					m.SetPtEtaPhi(vhCand.V.muons[0].p4.Pt(),0,vhCand.V.muons[0].p4.Phi());
					tau.SetPtEtaPhi(vhCand.V.taus[0].p4.Pt(),0,vhCand.V.taus[0].p4.Phi());
					metv.SetPtEtaPhi(vhCand.V.mets.at(0).p4.Pt(),0,vhCand.V.mets.at(0).p4.Phi());
					TVector3 bisector(m.Unit() + tau.Unit());
					bisector = bisector.Unit();
					PzetaVis  = (m+tau).Dot(bisector);
					PzetaMiss  =  metv.Dot(bisector);
					
				}
				if(Vtype == VHbbCandidate::Ztaue ){
				Nevents_Ztaue++;
					//cout << "Ztaue candidate" << endl;
					leptonForTop=vhCand.V.electrons[0].p4;
					vLeptons.set(vhCand.V.electrons[0],0,11); 
					vLeptons.set(vhCand.V.taus[0],1,15);
					nvlep=2;
					firstAddEle=1;
					firstAddTau=1;
					
					LorentzVector hadronictau = LorentzVector(vhCand.V.taus[0].p4.Px(),vhCand.V.taus[0].p4.Py(),vhCand.V.taus[0].p4.Pz(),vhCand.V.taus[0].p4.P());
					LorentzVector tauelectron = LorentzVector(vhCand.V.electrons[0].p4.Px(),vhCand.V.electrons[0].p4.Py(),vhCand.V.electrons[0].p4.Pz(),vhCand.V.electrons[0].p4.P());
					NSVfitStandalone::Vector measuredMET( vhCand.V.mets.at(0).p4.Px(), vhCand.V.mets.at(0).p4.Py(), 0);
					std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
					NSVfitStandalone::LorentzVector p1 = hadronictau;
					measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay,p1));
					NSVfitStandalone::LorentzVector p2 = tauelectron;
					measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay,p2));
					unsigned int flagSA = 0;
					NSVfitStandaloneAlgorithm algo(measuredTauLeptons,measuredMET,*metsig,flagSA);
					algo.addLogM(false);
					algo.integrate();
					svmass = algo.getMass(); // mass uncertainty not implemented yet
					//cout << "original mass " << vhCand.V.p4.M() << " SVfit mass " << svmass << " %diff " << (vhCand.V.p4.M()-svmass)/vhCand.V.p4.M() << endl;
					tau.SetPtEtaPhi(vhCand.V.taus[0].p4.Pt(),0,vhCand.V.taus[0].p4.Phi());
					e.SetPtEtaPhi(vhCand.V.electrons[0].p4.Pt(),0,vhCand.V.electrons[0].p4.Phi());
					metv.SetPtEtaPhi(vhCand.V.mets.at(0).p4.Pt(),0,vhCand.V.mets.at(0).p4.Phi());
					TVector3 bisector(tau.Unit() + e.Unit());
					bisector = bisector.Unit();
					PzetaVis  = (tau+e).Dot(bisector);
					PzetaMiss  =  metv.Dot(bisector);
					
				}
				
				if(Vtype == VHbbCandidate::Wmun ){
				Nevents_Wmun++;
					leptonForTop=vhCand.V.muons[0].p4;
					vLeptons.set(vhCand.V.muons[0],0,13); 
					float weightMuonRecoAndID = triggerWeight.scaleMuID(vLeptons.pt[0],vLeptons.eta[0]);
					float weightTrig1 = triggerWeight.scaleMuIsoHLT(vLeptons.pt[0],vLeptons.eta[0]);
					float cweightTrig = weightTrig1;
					weightTrig = weightMuonRecoAndID * cweightTrig;
					float weightTrig1OrMu30 = triggerWeight.scaleMuOr30IsoHLT(vLeptons.pt[0],vLeptons.eta[0]);
					weightTrigOrMu30 = weightMuonRecoAndID*weightTrig1OrMu30;
					nvlep=1;
					firstAddMu=1;
					
					}
				if( Vtype == VHbbCandidate::Wen ){
				Nevents_Wen++;
					leptonForTop=vhCand.V.electrons[0].p4;
					vLeptons.set(vhCand.V.electrons[0],0,11);
					nvlep=1;
					firstAddEle=1;
					weightTrigMay = triggerWeight.scaleSingleEleMay(vLeptons.pt[0],vLeptons.eta[0]);
					weightTrigV4 = triggerWeight.scaleSingleEleV4(vLeptons.pt[0],vLeptons.eta[0]);
					weightEleRecoAndId=triggerWeight.scaleID80Ele(vLeptons.pt[0],vLeptons.eta[0]) * triggerWeight.scaleRecoEle(vLeptons.pt[0],vLeptons.eta[0]);
					weightEleTrigJetMETPart=triggerWeight.scaleJet30Jet25(jet30pt,jet30eta)*triggerWeight.scalePFMHTEle(MET.et);
					weightEleTrigElePart= weightTrigV4; //this is for debugging only, checking only the V4 part
					
					weightTrigMay*=weightEleRecoAndId;
					weightTrigV4*=weightEleRecoAndId;
					weightTrigV4*=weightEleTrigJetMETPart;
					//	  weightTrig = weightTrigMay * 0.187 + weightTrigV4 * (1.-0.187); //FIXME: use proper lumi if we reload 2.fb
					weightTrig = (weightTrigMay * 0.215 + weightTrigV4 * 1.915)/ 2.13; //FIXME: use proper lumi if we reload 2.fb
				}
				if(Vtype == VHbbCandidate::Wtaun ){
				Nevents_Wtaun++;
					vLeptons.set(vhCand.V.taus[0],0,15); 
					nvlep=1;
					firstAddTau=1;
					
					if (vhCand.V.taus.size()>1){
						LorentzVector hadronictau = LorentzVector(vhCand.V.taus[0].p4.Px(),vhCand.V.taus[0].p4.Py(),vhCand.V.taus[0].p4.Pz(),vhCand.V.taus[0].p4.P());
						LorentzVector hadronictau2 = LorentzVector(vhCand.V.taus[1].p4.Px(),vhCand.V.taus[1].p4.Py(),vhCand.V.taus[1].p4.Pz(),vhCand.V.taus[1].p4.P());
						NSVfitStandalone::Vector measuredMET( vhCand.V.mets.at(0).p4.Px(), vhCand.V.mets.at(0).p4.Py(), 0);
						std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
						NSVfitStandalone::LorentzVector p1 = hadronictau;
						measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay,p1));
						NSVfitStandalone::LorentzVector p2 = hadronictau2;
						measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay,p2));
						unsigned int flagSA = 0;
						NSVfitStandaloneAlgorithm algo(measuredTauLeptons,measuredMET,*metsig,flagSA);
						algo.addLogM(false);
						algo.integrate();
						svmass = algo.getMass(); // mass uncertainty not implemented yet
						//cout << "original mass " << vhCand.V.p4.M() << " SVfit mass " << svmass << " %diff " << (vhCand.V.p4.M()-svmass)/vhCand.V.p4.M() << endl;
						tau.SetPtEtaPhi(vhCand.V.taus[0].p4.Pt(),0,vhCand.V.taus[0].p4.Phi());
						e.SetPtEtaPhi(vhCand.V.taus[1].p4.Pt(),0,vhCand.V.taus[1].p4.Phi());
						metv.SetPtEtaPhi(vhCand.V.mets.at(0).p4.Pt(),0,vhCand.V.mets.at(0).p4.Phi());
						TVector3 bisector(tau.Unit() + e.Unit());
						bisector = bisector.Unit();
						PzetaVis  = (tau+e).Dot(bisector);
						PzetaMiss  =  metv.Dot(bisector);
						
				}
				}
					if(Vtype == VHbbCandidate::Ztautau &&  Vtype != VHbbCandidate::Wtaun ){
					Nevents_Ztautau++;
						vLeptons.set(vhCand.V.taus[0],0,15); 
						vLeptons.set(vhCand.V.taus[1],0,15); 
						nvlep=2;
						firstAddTau=2;
						
						if (vhCand.V.taus.size()>1){
							LorentzVector hadronictau = LorentzVector(vhCand.V.taus[0].p4.Px(),vhCand.V.taus[0].p4.Py(),vhCand.V.taus[0].p4.Pz(),vhCand.V.taus[0].p4.P());
							LorentzVector hadronictau2 = LorentzVector(vhCand.V.taus[1].p4.Px(),vhCand.V.taus[1].p4.Py(),vhCand.V.taus[1].p4.Pz(),vhCand.V.taus[1].p4.P());
							NSVfitStandalone::Vector measuredMET( vhCand.V.mets.at(0).p4.Px(), vhCand.V.mets.at(0).p4.Py(), 0);
							std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
							NSVfitStandalone::LorentzVector p1 = hadronictau;
							measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay,p1));
							NSVfitStandalone::LorentzVector p2 = hadronictau2;
							measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay,p2));
							unsigned int flagSA = 0;
							NSVfitStandaloneAlgorithm algo(measuredTauLeptons,measuredMET,*metsig,flagSA);
							algo.addLogM(false);
							algo.integrate();
							svmass = algo.getMass(); // mass uncertainty not implemented yet
							//cout << "original mass " << vhCand.V.p4.M() << " SVfit mass " << svmass << " %diff " << (vhCand.V.p4.M()-svmass)/vhCand.V.p4.M() << endl;
							tau.SetPtEtaPhi(vhCand.V.taus[0].p4.Pt(),0,vhCand.V.taus[0].p4.Phi());
							e.SetPtEtaPhi(vhCand.V.taus[1].p4.Pt(),0,vhCand.V.taus[1].p4.Phi());
							metv.SetPtEtaPhi(vhCand.V.mets.at(0).p4.Pt(),0,vhCand.V.mets.at(0).p4.Phi());
							TVector3 bisector(tau.Unit() + e.Unit());
							bisector = bisector.Unit();
							PzetaVis  = (tau+e).Dot(bisector);
							PzetaMiss  =  metv.Dot(bisector);
							
						}
					}
				Pzeta25 = PzetaVis - (0.25*PzetaMiss);
				Pzeta45 = PzetaVis - (0.45*PzetaMiss);
				Pzeta65 = PzetaVis - (0.65*PzetaMiss);
				Pzeta85 = PzetaVis - (0.85*PzetaMiss);
				if(Vtype == VHbbCandidate::Zbb )Nevents_Zbb++;

						
				if(isMC_)
					{
						weightTrigMET80 =  triggerWeight.scaleMET80(MET.et);
						weightTrigMET100 =  triggerWeight.scaleMET80(MET.et);
						weightTrig2CJet20 = triggerWeight.scale2CentralJet( jet10pt, jet10eta);
						weightTrigMET150 = triggerWeight.scaleMET150(MET.et);
						weightTrigMET802CJet= weightTrigMET80 * weightTrig2CJet20;
						weightTrigMET1002CJet= weightTrigMET100 * weightTrig2CJet20;
					}
					//cout << "right before Znunu event definition" << endl;
				if( Vtype == VHbbCandidate::Znn ){
				//cout << "have a Znunu event" << endl;
				Nevents_Znunu++;
					nvlep=0;
					float weightTrig1 = triggerWeight.scaleMetHLT(vhCand.V.mets.at(0).p4.Pt());
					weightTrigMETLP = weightTrig1;
					weightTrig = weightTrigMET150 + weightTrigMET802CJet  - weightTrigMET802CJet*weightTrigMET150;
				}
				
				if(weightTrigMay < 0) weightTrigMay=weightTrig;
				if(weightTrigV4 < 0) weightTrigV4=weightTrig;
				if(!isMC_)
					{
						weightTrig = 1.; 
						weightTrigMay = 1.;
						weightTrigV4 = 1.;
						weightEleRecoAndId= 1.;
						weightEleTrigJetMETPart= 1.;
						weightEleTrigElePart= 1.;
						weightEleTrigEleAugPart=1.;
						weightTrigMET80= 1.;
						weightTrigMET100= 1.;
						weightTrig2CJet20= 1.;
						weightTrigMET150= 1.;
						weightTrigMET802CJet= 1.;
						weightTrigMET1002CJet= 1.;
						weightTrigMETLP = 1.;
						
					} 
				aLeptons.reset();
				nalep=0;
				if(fromCandidate)
					{
						for(size_t j=firstAddMu;j< vhCand.V.muons.size();j++) aLeptons.set(vhCand.V.muons[j],nalep++,13);
						for(size_t j=firstAddEle;j< vhCand.V.electrons.size();j++) aLeptons.set(vhCand.V.electrons[j],nalep++,11);
						for(size_t j=firstAddTau;j< vhCand.V.taus.size();j++) aLeptons.set(vhCand.V.taus[j],nalep++,15);
					}
				else
					{
						
						for(size_t j=0;j< iEvent->muInfo.size();j++)
							{ 
								if((j!= vhCand.V.firstLepton && j!= vhCand.V.secondLepton) || ((Vtype != VHbbCandidate::Wmun ) && (Vtype != VHbbCandidate::Zmumu ) && (Vtype != VHbbCandidate::Zemu) && (Vtype != VHbbCandidate::Ztaumu)) )
									aLeptons.set(iEvent->muInfo[j],nalep++,13);
							}
						for(size_t j=0;j< iEvent->eleInfo.size();j++)
							{ 
								if((j!= vhCand.V.firstLepton && j!= vhCand.V.secondLepton) || ((Vtype != VHbbCandidate::Wen ) && (Vtype != VHbbCandidate::Zee ) && (Vtype != VHbbCandidate::Zemu) && (Vtype != VHbbCandidate::Ztaue)))
									aLeptons.set(iEvent->eleInfo[j],nalep++,11);
							}
						for(size_t j=0;j< iEvent->tauInfo.size();j++)
							{ 
								if((j!= vhCand.V.firstLepton && j!= vhCand.V.secondLepton) || ((Vtype != VHbbCandidate::Ztaumu) && (Vtype != VHbbCandidate::Wtaun) && (Vtype != VHbbCandidate::Ztaue) && (Vtype != VHbbCandidate::Ztautau)) ){
									aLeptons.set(iEvent->tauInfo[j],nalep++,15);
									//cout << "tau Pt: " << iEvent->tauInfo[j].p4.Pt() << endl;
								}
							}
					}
				
				
				if(isMC_)
					{
						if(Verbage)std::cout << "BTAGSF " <<  btagJetInfos.size() << " " << btag.weight<BTag1Tight2CustomFilter>(btagJetInfos) << std::endl;
						if ( btagJetInfos.size()< 10) 
							{
								btag1T2CSF = btag.weight<BTag1Tight2CustomFilter>(btagJetInfos);
								btag2TSF = btag.weight<BTag2TightFilter>(btagJetInfos);
								btag1TSF = btag.weight<BTag1TightFilter>(btagJetInfos);
								btagA0CSF = btag.weight<BTagAntiMax0CustomFilter>(btagJetInfos);
								btagA0TSF = btag.weight<BTagAntiMax0TightFilter>(btagJetInfos);
								btag2CSF = btag.weight<BTag2CustomFilter>(btagJetInfos);
								btag1TA1C = btag.weight<BTag1TightAndMax1CustomFilter>(btagJetInfos);
							}
						else
							{
								std::cout << "WARNING:  combinatorics for " << btagJetInfos.size() << " jets is too high (>=10). use SF=1 " << std::endl;
								//TODO: revert to random throw  for this cases
								btag1T2CSF = 1.;
								btag2TSF =  1.;
								btag1TSF =  1.;
								btagA0CSF = 1.;
								btagA0TSF = 1.;
								btag2CSF = 1.;
								btag1TA1C = 1.;
							}
					}
				
				if(maxBtag > -99999)
					{ 
						TopHypo topQuark = TopMassReco::topMass(leptonForTop,bJet,vhCand.V.mets.at(0).p4);
						top.mass = topQuark.p4.M();
						top.pt = topQuark.p4.Pt();
						top.wMass = topQuark.p4W.M();
					} else {
						top.mass = -99;
						top.pt = -99;
						top.wMass = -99;
					}
				
				
				
				//FIXME: too much  warnings... figure out why 
				//         gendrcc=aux.genCCDeltaR(); 
				//         gendrbb=aux.genBBDeltaR(); 
				
				if(isMC_) {
					
					genTau.mass = -1.;
					genTau.pt = 0.;
					genTau.eta = -99.;
					genTau.phi = -99.;
					fwlite::Handle<reco::GenParticleCollection> genParticles;
					genParticles.getByLabel(ev,"savedGenParticles");
					for(size_t i = 0; i < genParticles->size(); ++ i) {
						const reco::GenParticle p = (*genParticles)[i];
						int id = p.pdgId();
						//int st = p.status();
						//int ch = p.charge();
						TLorentzVector mc_p4 = GENPTOLOR(p);
						if (abs(id)==15) {
							TLorentzVector visibleP = VisibleP(p);
							//cout << "found tau" << ", mass: " << visibleP.M() << ", pT: " << visibleP.Pt() << endl;
							genTau.mass = visibleP.M();
							genTau.pt = visibleP.Pt();
							genTau.eta = visibleP.Eta();
							genTau.phi = visibleP.Phi();
							
						}
					}
				}
				
				
				genZpt=aux.mcZ.size() > 0 ? aux.mcZ[0].p4.Pt():-99;
				genWpt=aux.mcW.size() > 0 ? aux.mcW[0].p4.Pt():-99;
				
				bool tauDauFromZ = false;
				genZmumu = false, genZee = false, genZnunu = false;
				// Z* is status=3 and Nmother=2 (q and qbar)  
				// Z is status=2 and Ndau=2 (mup and mum)  
				if(Verbage)std::cout << "Gen Zstar info " << endl;
				for (unsigned int i=0; i<aux.mcZ.size(); i++){
					if (aux.mcZ[i].status==3) {
						genZstar.mass = aux.mcZ[i].p4.M();
						genZstar.pt = aux.mcZ[i].p4.Pt();
						genZstar.eta = aux.mcZ[i].p4.Eta();
						genZstar.phi = aux.mcZ[i].p4.Phi();
						genZstar.status = aux.mcZ[i].status;
						genZstar.charge = aux.mcZ[i].charge;
						if (aux.mcZ[i].momid!=-99) genZstar.momid  =  aux.mcZ[i].momid ;
					}
					
					if (aux.mcZ[i].dauid.size()==2 && aux.mcZ[i].status == 2) {
						if ((abs(aux.mcZ[i].dauid[0])==12 && abs(aux.mcZ[i].dauid[1])==12) || (abs(aux.mcZ[i].dauid[0])==14 && abs(aux.mcZ[i].dauid[1])==14) || (abs(aux.mcZ[i].dauid[0])==16 && abs(aux.mcZ[i].dauid[1])==16)) genZnunu = true;
						if ( abs(aux.mcZ[i].dauid[0])==13 || abs(aux.mcZ[i].dauid[0])==11 || abs(aux.mcZ[i].dauid[0])==15 ) {
							for (unsigned int d = 0; d != aux.mcZ[i].dauid.size(); d++){
								if (abs(aux.mcZ[i].dauid[d])==15) tauDauFromZ = true;
							}// loop Z daughters
							genZ.mass = aux.mcZ[i].p4.M();
							genZ.pt = aux.mcZ[i].p4.Pt();
							genZ.eta = aux.mcZ[i].p4.Eta();
							genZ.phi = aux.mcZ[i].p4.Phi();
							genZ.status = aux.mcZ[i].status;
							genZ.charge = aux.mcZ[i].charge;
							if ( aux.mcZ[i].momid!=-99) genZ.momid = aux.mcZ[i].momid;
							if (abs(aux.mcZ[i].dauid[0])==13 && abs(aux.mcZ[i].dauid[1])==13) genZmumu = true;
							if (abs(aux.mcZ[i].dauid[0])==11 && abs(aux.mcZ[i].dauid[1])==11) genZee = true;
						}
					}
				}
				
				int muTauDau = 0, eTauDau = 0, hasPions = 0, hadTau = 0;
				genTauMu.mass = -1.;
				genTauMu.pt = 0.;
				genTauMu.eta = -99.;
				genTauMu.phi = -99.;
				if (aux.mcTau.size() > 1 && tauDauFromZ==true){
					//cout << "********" << endl;
					//cout << "there are " << aux.mcTau.size() << " taus in event " << ev.id().event() << endl;
					//cout << "********" << endl;
					TLorentzVector hadronicTauDecays, dummy;
					for (unsigned int t=0; t<aux.mcTau.size(); t++){
						bool hasMu = false, hasE = false;
						for (unsigned int dt = 0; dt != aux.mcTau[t].dauid.size(); dt++){
							//cout << "tau[" << t << "].daughter[" << dt << "]: " << aux.mcTau[t].dauid[dt] << endl;
							if (abs(aux.mcTau[t].dauid[dt])==13){
								hasMu=true;
								genTauMu.mass = aux.mcTau[t].dauFourMomentum[dt].M();
								genTauMu.pt = aux.mcTau[t].dauFourMomentum[dt].Pt();
								genTauMu.eta = aux.mcTau[t].dauFourMomentum[dt].Eta();
								genTauMu.phi = aux.mcTau[t].dauFourMomentum[dt].Phi();
							}
							else if (abs(aux.mcTau[t].dauid[dt])==11) hasE=true;
							else{
								dummy.SetPtEtaPhiM(aux.mcTau[t].dauFourMomentum[dt].Pt(), aux.mcTau[t].dauFourMomentum[dt].Eta(), aux.mcTau[t].dauFourMomentum[dt].Phi(), aux.mcTau[t].dauFourMomentum[dt].M());
								hadronicTauDecays = hadronicTauDecays + dummy;
							}// hadro tau
							uint32_t pdgid=abs(aux.mcTau[t].dauid[dt]);
							uint32_t hundreds=pdgid%1000;
							if (hundreds>=111 && hundreds<=215) hasPions=true;
						}// loop on tau daughters
						if (!hasMu && !hasE){
							hadTau++;
							////genTau.mass = hadronicTauDecays.M();
							////genTau.pt = hadronicTauDecays.Pt();
							////genTau.eta = hadronicTauDecays.Eta();
							////genTau.phi = hadronicTauDecays.Phi();
							//genTau.status = hadronicTauDecays.status;
							//genTau.charge = hadronicTauDecays.charge;
							//genTau.momid = hadronicTauDecays.momid;
							//cout << "tau eta: " << aux.mcTau[t].p4.Eta() << " mass: " << aux.mcTau[t].p4.M() << endl;
							//cout << "momid: " << aux.mcTau[t].momid << endl;
						}// hadTau
						//if (hasMu && abs(aux.mcTau[t].momid)==15) muTauDau++;
						//if (hasE && abs(aux.mcTau[t].momid)==15) eTauDau++;
						if (hasMu) muTauDau++;
						if (hasE) eTauDau++;
					}// loop on taus
				}// if 2 gen taus in event
				//cout << "********" << endl;
				//cout << "tau decays: " << hadTau << " hadronic, " << muTauDau << " muon, " << eTauDau << " electron." << endl;
				//cout << "********" << endl;
				
				genZtaumu = false;
				if (muTauDau >= 1 && hadTau >= 1) genZtaumu = true;
				
				for (unsigned int i=0; i<aux.mcW.size(); i++){
					if  ( aux.mcW[i].momid==6 && aux.mcW[i].dauid.size()>1 ){
						genTop.wdau1mass= aux.mcW[i].dauFourMomentum[0].M();
						genTop.wdau1pt= aux.mcW[i].dauFourMomentum[0].Pt();
						genTop.wdau1eta= aux.mcW[i].dauFourMomentum[0].Eta();
						genTop.wdau1phi= aux.mcW[i].dauFourMomentum[0].Phi();
						genTop.wdau1id= aux.mcW[i].dauid[0];
						
						genTop.wdau2mass= aux.mcW[i].dauFourMomentum[1].M();
						genTop.wdau2pt= aux.mcW[i].dauFourMomentum[1].Pt();
						genTop.wdau2eta= aux.mcW[i].dauFourMomentum[1].Eta();
						genTop.wdau2phi= aux.mcW[i].dauFourMomentum[1].Phi();
						genTop.wdau2id= aux.mcW[i].dauid[1];
						
					}
					
					if(Verbage)std::cout << "Gen Tbar info " << endl;
					
					if  ( aux.mcW[i].momid==-6 && aux.mcW[i].dauid.size()>1 ){
						genTbar.wdau1mass= aux.mcW[i].dauFourMomentum[0].M();
						genTbar.wdau1pt= aux.mcW[i].dauFourMomentum[0].Pt();
						genTbar.wdau1eta= aux.mcW[i].dauFourMomentum[0].Eta();
						genTbar.wdau1phi= aux.mcW[i].dauFourMomentum[0].Phi();
						genTbar.wdau1id= aux.mcW[i].dauid[0];
						
						genTbar.wdau2mass= aux.mcW[i].dauFourMomentum[1].M();
						genTbar.wdau2pt= aux.mcW[i].dauFourMomentum[1].Pt();
						genTbar.wdau2eta= aux.mcW[i].dauFourMomentum[1].Eta();
						genTbar.wdau2phi= aux.mcW[i].dauFourMomentum[1].Phi();
						genTbar.wdau2id= aux.mcW[i].dauid[1];
						
					}
					
					genWmunu = false;
					genWenu = false;
					if (aux.mcW[i].status==3 ) {
						genWstar.mass = aux.mcW[i].p4.M();
						genWstar.pt = aux.mcW[i].p4.Pt();
						genWstar.eta = aux.mcW[i].p4.Eta();
						genWstar.phi = aux.mcW[i].p4.Phi();
						genWstar.status = aux.mcW[i].status;
						genWstar.charge = aux.mcW[i].charge;
						if ( aux.mcW[i].momid!=-99) genWstar.momid =  aux.mcW[i].momid;
					}
					
					if(Verbage)std::cout << "Gen W info " << endl;
					
					if (aux.mcW[i].dauid.size()>1 && (abs(aux.mcW[i].dauid[0])==13 || abs(aux.mcW[i].dauid[0])==11 )) {
						genW.mass = aux.mcW[i].p4.M();
						genW.pt = aux.mcW[i].p4.Pt();
						genW.eta = aux.mcW[i].p4.Eta();
						genW.phi = aux.mcW[i].p4.Phi();
						genW.status = aux.mcW[i].status;
						genW.charge = aux.mcW[i].charge;
						if (aux.mcW[i].momid!=-99) genW.momid =  aux.mcW[i].momid;
						if ( (abs(aux.mcW[i].dauid[0])==13 && abs(aux.mcW[i].dauid[1])==14) || (abs(aux.mcW[i].dauid[0])==14 && abs(aux.mcW[i].dauid[1])==13) ){
							genWmunu = true;
						}  
						if ( (abs(aux.mcW[i].dauid[0])==11 && abs(aux.mcW[i].dauid[1])==12) || (abs(aux.mcW[i].dauid[0])==12 && abs(aux.mcW[i].dauid[1])==11) ) genWenu = true;
					}
				}
				if(Verbage)std::cout << "Gen b info " << endl;
				for (unsigned int i=0; i<aux.mcB.size(); i++){
					if (abs(aux.mcB[i].momid)!=5) {
						genB.mass = aux.mcB[i].p4.M();
						genB.pt = aux.mcB[i].p4.Pt();
						genB.eta = aux.mcB[i].p4.Eta();
						genB.phi = aux.mcB[i].p4.Phi();
						genB.status = aux.mcB[i].status;
						genB.charge = aux.mcB[i].charge;
						if (aux.mcB[i].momid!=-99) genB.momid = aux.mcB[i].momid;
					}
					
					if  ( aux.mcB[i].momid==6 ){
						genTop.bmass = aux.mcB[i].p4.M();
						genTop.bpt = aux.mcB[i].p4.Pt();
						genTop.beta = aux.mcB[i].p4.Eta();
						genTop.bphi = aux.mcB[i].p4.Phi();
						genTop.bstatus = aux.mcB[i].status;
						
					}
				}
				if(Verbage)std::cout << "Gen Bbar info " << endl;
				
				for (unsigned int i=0; i<aux.mcBbar.size(); i++){
					if (abs(aux.mcBbar[i].momid)!=5 ) {
						genBbar.mass = aux.mcBbar[i].p4.M();
						genBbar.pt = aux.mcBbar[i].p4.Pt();
						genBbar.eta = aux.mcBbar[i].p4.Eta();
						genBbar.phi = aux.mcBbar[i].p4.Phi();
						genBbar.status = aux.mcBbar[i].status;
						if (aux.mcBbar[i].momid!=-99) genBbar.momid =  aux.mcBbar[i].momid;
					}
					if  ( aux.mcBbar[i].momid==-6 ){
						genTbar.bmass = aux.mcBbar[i].p4.M();
						genTbar.bpt = aux.mcBbar[i].p4.Pt();
						genTbar.beta = aux.mcBbar[i].p4.Eta();
						genTbar.bphi = aux.mcBbar[i].p4.Phi();
						genTbar.bstatus = aux.mcBbar[i].status;
						
						
					}
					
				}
				
				
				if(Verbage)std::cout << "Gen Higgs info " << endl;
				if (aux.mcH.size()>0) {
					genH.mass = aux.mcH[0].p4.M();
					genH.pt = aux.mcH[0].p4.Pt();
					genH.eta = aux.mcH[0].p4.Eta();
					genH.phi = aux.mcH[0].p4.Phi();
					genH.status = aux.mcH[0].status;
					genH.charge = aux.mcH[0].charge;
					if (aux.mcH[0].momid!=-99) genH.momid =  aux.mcH[0].momid;
				}
				
				
				
				
				WminusMode=-99;
				WplusMode=-99;
				for(unsigned int j=0; j< aux.mcW.size();j++)
					{
						for(unsigned int k=0;k< aux.mcW[j].dauid.size();k++)
							{
								int idd=abs(aux.mcW[j].dauid[k]);
								if(idd==11 || idd==13 || idd==15|| (idd<=5 && idd >=1)) 
									{
										if(WminusMode==-99 && aux.mcW[j].charge ==-1) WminusMode = idd;
										if(WplusMode==-99 && aux.mcW[j].charge ==+1) WplusMode = idd;
									}
							}
						/*
						 /// now check if a semileptonic W is also in a bjets....      
						 if ( ( (WminusMode==11 || WminusMode==13 || WminusMode==15  ) || (WplusMode==11 || WplusMode==13 || WplusMode==15  ))  && deltaR(vhCand.H.jets[0].p4.Eta(),vhCand.H.jets[0].p4.Phi(), aux.mcW[j].p4.Eta(), aux.mcW[j].p4.Phi())<0.3 ) hJets.isSemiLeptMCtruth[0]=1;
						 if ( ( (WminusMode==11 || WminusMode==13 || WminusMode==15  ) || (WplusMode==11 || WplusMode==13 || WplusMode==15  ))  && deltaR(vhCand.H.jets[1].p4.Eta(),vhCand.H.jets[1].p4.Phi(), aux.mcW[j].p4.Eta(), aux.mcW[j].p4.Phi())<0.3 ) hJets.isSemiLeptMCtruth[1]=1;
						 
						 for( int j=0; j < naJets && j < MAXJ; j++ ) 
						 {
						 if ((idd==11 || idd==13 || idd==15  )  && deltaR(vhCand.additionalJets[j].p4.Eta(),vhCand.additionalJets[j].p4.Phi(), aux.mcW[j].p4.Eta(), aux.mcW[j].p4.Phi()) <0.3) aJets.isSemiLept[j]=1;
						 
						 }
						 */      
						
					}
				/// Compute pull angle from AK7
				if(!fromCandidate){
					std::vector<VHbbEvent::SimpleJet> ak7wrt1(iEvent->simpleJets3);
					std::vector<VHbbEvent::SimpleJet> ak7wrt2(iEvent->simpleJets3);
					if(ak7wrt1.size() > 1){
						CompareDeltaR deltaRComparatorJ1(vhCand.H.jets[0].p4);
						CompareDeltaR deltaRComparatorJ2(vhCand.H.jets[1].p4);
						std::sort( ak7wrt1.begin(),ak7wrt1.end(),deltaRComparatorJ1 );
						std::sort( ak7wrt2.begin(),ak7wrt2.end(),deltaRComparatorJ2 );
						std::vector<VHbbEvent::SimpleJet> ak7_matched;
						// if the matched are different save them
						if(ak7wrt1[0].p4.DeltaR(ak7wrt2[0].p4) > 0.1) {
							ak7_matched.push_back(ak7wrt1[0]);
							ak7_matched.push_back(ak7wrt2[0]);
						}
						// else look at the second best
						else{
							// ak7wrt1 is best
							if( ak7wrt1[1].p4.DeltaR(vhCand.H.jets[0].p4) < ak7wrt2[1].p4.DeltaR(vhCand.H.jets[1].p4))
								{
									ak7_matched.push_back(ak7wrt1[1]);
									ak7_matched.push_back(ak7wrt2[0]);
								}
							else
								{
									ak7_matched.push_back(ak7wrt1[0]);
									ak7_matched.push_back(ak7wrt2[1]);
								}
						}
						CompareJetPt ptComparator;
						std::sort( ak7_matched.begin(),ak7_matched.end(),ptComparator );
						if(ak7_matched[0].p4.DeltaR(vhCand.H.jets[0].p4) < 0.5
						   and ak7_matched[1].p4.DeltaR(vhCand.H.jets[1].p4) < 0.5)
							{
								deltaPullAngleAK7 =  VHbbCandidateTools::getDeltaTheta(ak7_matched[0],ak7_matched[1]);
								deltaPullAngle2AK7 =  VHbbCandidateTools::getDeltaTheta(ak7_matched[1],ak7_matched[0]);
							}
					}
				}
				
				_outTree->Fill();
				if(Verbage)std::cout << "After Fill outTree" << endl;
				
			}// closed event loop
		
		std::cout << "closing the file: " << inputFiles_[iFile] << std::endl;
		inFile->Close();
		// close input file
	} // loop on files
	
	
	std::cout << "Events: " << ievt <<std::endl;
	std::cout << "Number of Events with 15 GeV Pt threshold for jets " << Nevents_LowJetPt << endl;

	std::cout << "Number of Zmumu Events " << Nevents_Zmumu << endl;
	std::cout << "Number of Zee Events " << Nevents_Zee << endl;
	std::cout << "Number of Wmun Events " << Nevents_Wmun << endl;
	std::cout << "Number of Wen Events " << Nevents_Wen << endl;
	std::cout << "Number of Znunu Events " << Nevents_Znunu << endl;
	std::cout << "Number of Zemu Events " << Nevents_Zemu << endl;
	std::cout << "Number of Ztaumu Events " << Nevents_Ztaumu << endl;
	std::cout << "Number of Ztaue Events " << Nevents_Ztaue << endl;
	std::cout << "Number of Wtaun Events " << Nevents_Wtaun << endl;
	std::cout << "Number of Ztautau Events " << Nevents_Ztautau << endl;
	std::cout << "Number of Zbb Events " << Nevents_Zbb << endl;
	std::cout << "More than one Candidate " << Nmulti << endl;
    
    
	_outFile->cd();
    
	_outTree->Write();
	_outFile->Write();
	_outFile->Close();
	return 0;
}


