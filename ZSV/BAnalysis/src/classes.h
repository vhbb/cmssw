#include "ZSV/BAnalysis/interface/SimSecondaryVertex.h"
#include "ZSV/BAnalysis/interface/SimBHadron.h"
#include "ZSV/BAnalysis/interface/SimBJet.h"
#include "ZSV/BAnalysis/interface/SimBEvent.h"
#include "ZSV/BAnalysis/interface/RecoBVertexCand.h"

namespace {
 namespace {

SimPrimaryVertex ap;
edm::Wrapper<SimPrimaryVertex> wrp;
SimSecondaryVertex aa;
SimSecondaryVertexCollection aa1;
SimSecondaryVertexRefProd ssvrp;
SimSecondaryVertexRefVector ssvrv;
edm::Wrapper<SimSecondaryVertexCollection> wr;
edm::Ptr<SimSecondaryVertex> ssv;  
SimBHadron bb;
SimBHadronCollection bb1;
SimBHadronRefProd sbhrp;
SimBHadronRefVector sbhrv;
edm::Wrapper<SimBHadronCollection> wr2;
SimBEvent be;
SimBEventCollection be1;
SimBEventRefProd sberp;
SimBEventRefVector sberv;
edm::Wrapper<SimBEvent> wr4;
edm::Wrapper<SimBEventCollection> wr5;
RecoBVertexCand bvc;
RecoBVertexCandCollection bvc1;
RecoBVertexCandRefProd rbvcrp;
RecoBVertexCandRefVector rbvcrv;
edm::Wrapper<RecoBVertexCandCollection> wr6;



JetCollection jc; 
edm::Wrapper<JetCollection> wr3;

 }
}
