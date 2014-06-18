#ifndef SimBEvent_H
#define SimBEvent_H
//#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/JetReco/interface/Jet.h"
#include <vector>
//#include "DataFormats/Math/interface/Vector3D.h"

//#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "DataFormats/Math/interface/LorentzVector.h"
//#include "Math/LorentzVector.h"
//#include "DataFormats/JetReco/interface/CaloJetCollection.h"
//#include "DataFormats/JetReco/interface/CaloJet.h"

//#include <TVector3.h>

#include "ZSV/BAnalysis/interface/SimBHadron.h"
#include "ZSV/BAnalysis/interface/RecoBVertexCand.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class SimBEvent{
public:

  SimBEvent(){}
  
  int eventNo; 
  int runNo; 
  int nSBH; 
  float dRbb; //-1 if != 2 SBH
  float dRvv; //-1 if != 2 unique RSV
  float dEtavv; //-1 if != 2 unique RSV
  float dPhivv; //-1 if != 2 unique RSV
  float dRvvCutsOk; //-1 if != 2 unique RSV WHICH SURVIVED CUTS
  float dEtavvCutsOk; //-1 if != 2 unique RSV WHICH SURVIVED CUTS
  float dPhivvCutsOk; //-1 if != 2 unique RSV WHICH SURVIVED CUTS
  float dRvvCutsOkCentral; //-1 if != 2 unique RSV WHICH SURVIVED CUTS IN GOOD KIN. REGION 
  float dEtavvCutsOkCentral; //-1 if != 2 unique RSV WHICH SURVIVED CUTS IN GOOD KIN. REGION 
  float dPhivvCutsOkCentral; //-1 if != 2 unique RSV WHICH SURVIVED CUTS IN GOOD KIN. REGION 
  float ptB1; //-1 if != 2 SBH
  float ptB2; //-1 if != 2 SBH
  float etaB1; //-1 if != 2 SBH
  float etaB2; //-1 if != 2 SBH
  SimBHadronRefProd SimBHadrons; 

  int nVertices;
  int nUniqueVertices;
  int nBVertices; //not nec. unique
  int nUBVertices; 
  int nDVertices; //not nec. unique
  int nUDVertices; 
  int nSecVertices; //not nec. unique
  int nUSecVertices; 
  int nTertVertices; // not nec. unique
  int nUTertVertices; 
  int nUCutsOkVertices; //unique
  int nUCutsOkVerticesCentral; //unique && in central region
  int nUMatchedVertices; //unique
  int nMCVertices; 
  int flavorV1; // -1 no nUCutsOkVertices==2 event, 
  int flavorV2; // 1 b 2 c 3 l 4 f
  int flavorV1Central; // 1 b 2 c 3 l 4 f
  int flavorV2Central; // 1 b 2 c 3 l 4 f
  int flavorCat; // 1bb 2bc 3bl 4bf 5cc 6other -1not2uselvert
  int flavorCatCentral; // 1bb 2bc 3bl 4bf 5cc 6other -1not2uselvert
  int goodVert1; // nUCutsOkVertices ==2 
  int goodVert2; 
  int goodVert1Central; 
  int goodVert2Central; 
  float scalarMassCutsOkVertexPair; 
  float scalarMassCutsOkCentralVertexPair; 
  float cosDir12Td2; //for cos variables 1 is always the vertex closer to pv	
  float cosTd1Fd1; //all three variables are -1.777 if not two vertices selected
  float cosTd2Fd2; 
  float cosDir12Td2Central; //for cos variables 1 is always the vertex closer to pv	
  float cosTd1Fd1Central; //all three variables are -1.777 if not two vertices selected
  float cosTd2Fd2Central; 
  bool simOk2;   // exactly two B, both |eta|<2.0 && ncSTBeta25pt05 >= 2
  bool trackOk2; // exactly two B, both nRTBeta25pt05 >= 2
  bool trackGoodOk2; // exactly two B, both with two Good Tracks (n valid hits >=7, pt>0.9 GeV)
  bool vertOk2;  // exactly two B, both nVertices>0
  bool mcvertOk2;  // exactly two B, both nVertices>0
  bool cutsOk2;  // exactly two B, both nVerticesSurviveCuts > 0 ((or ==1??))
  bool cutsOk2NEW;  // exactly two B, both nVerticesSurviveCuts > 0 ((or ==1??))
  bool mccutsOk2;  // exactly two B, both nVerticesSurviveCuts > 0 ((or ==1??))
  bool matchOk2; // exactly two B, both nVertices>0 with dR_BSV < 0.1
  bool mcmatchOk2; // exactly two B, both nVertices>0 with dR_BSV < 0.1
  bool cutsmatchOk2; // exactly two B, both nVertices>0 with dR_BSV < 0.1
  bool cutsmatchOk2NEW; // exactly two B, both nVertices>0 with dR_BSV < 0.1
  bool mccutsmatchOk2; // exactly two B, both nVertices>0 with dR_BSV < 0.1
  bool paircutsOk; //not yet defined...
  bool selPoss;  //event selection possible means : simOk2 && vertOk2
  bool matchPoss;//matching possible means : simOk2 && matchOk2 
  bool selected; 
  bool selectedCentral; 
  bool veto4B; 
  bool veto4BCentral;
  float ptrel;  
  float ptrelCentral;  

  RecoBVertexCandRefProd RecoBVCands; 
  reco::VertexRefProd MCBVertices;

  const SimBHadronCollection & sbh() const {return (*SimBHadrons);}
  const reco::VertexCollection & mcbtrackvert() const {return (*MCBVertices);}
  const RecoBVertexCandCollection  & recovert() const {return (*RecoBVCands);}

};

typedef  std::vector<SimBEvent> SimBEventCollection;
typedef  edm::Ref<SimBEventCollection> SimBEventRef;
typedef  edm::RefProd<SimBEventCollection> SimBEventRefProd;
typedef  edm::RefVector<SimBEventCollection> SimBEventRefVector;

#endif

