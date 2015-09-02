#include "VHbbAnalysis/Heppy/interface/SignedImpactParameter.h"
#include "MagneticField/UniformEngine/src/UniformMagneticField.h"
#include "MagneticField/ParametrizedEngine/plugins/OAEParametrizedMagneticField.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include <vector>
#include <TMath.h>
#include <Math/SVector.h>

//-------------------------------------------------------------------------------
// CV: temporary copy of MagneticField/ParametrizedEngine/plugins/OAEParametrizedMagneticField.cc ,
//     in order to use this code in CMSSW 7_2_x

#include <FWCore/MessageLogger/interface/MessageLogger.h>

using namespace std;
using namespace magfieldparam;
 

OAEParametrizedMagneticField::OAEParametrizedMagneticField(string T) : 
  theParam(T){}

OAEParametrizedMagneticField::OAEParametrizedMagneticField(const edm::ParameterSet& parameters) : 
  theParam(parameters.getParameter<string>("BValue")) {}

OAEParametrizedMagneticField::~OAEParametrizedMagneticField() {}

GlobalVector
OAEParametrizedMagneticField::inTesla(const GlobalPoint& gp) const {
  if (isDefined(gp)) {
    return inTeslaUnchecked(gp);
  } else {
    edm::LogWarning("MagneticField|FieldOutsideValidity") << " Point " << gp << " is outside the validity region of OAEParametrizedMagneticField";
    return GlobalVector();
  }
}

namespace {
  constexpr float ooh = 1./100;
}

GlobalVector
OAEParametrizedMagneticField::inTeslaUnchecked(const GlobalPoint& gp) const {
  float x[3] = {gp.x()*ooh, gp.y()*ooh, gp.z()*ooh};
  float B[3];
  theParam.getBxyz(x,B);
  return GlobalVector(B[0], B[1], B[2]);
}

bool
OAEParametrizedMagneticField::isDefined(const GlobalPoint& gp) const {
  return (gp.perp2()<(115.f*115.f) && fabs(gp.z())<280.f);
}

// CV: temporary copy of MagneticField/ParametrizedEngine/plugins/TkBfield.cc ,
//     in order to use this code in CMSSW 7_2_x

#include "MagneticField/ParametrizedEngine/plugins/TkBfield.h"
 
#include "FWCore/Utilities/interface/Exception.h"

#include <cmath>
#include<algorithm>
 
namespace {
  BCylParam<float>  fpar1{4.90541f,17.8768f,2.02355f,0.0210538f,0.000321885f,2.37511f,0.00326725f,2.07656f,1.71879f}; // 2.0T-2G
  BCylParam<float>  fpar2{4.41982f,15.7732f,3.02621f,0.0197814f,0.000515759f,2.43385f,0.00584258f,2.11333f,1.76079f}; // 3.0T-2G
  BCylParam<float>  fpar3{4.30161f,15.2586f,3.51926f,0.0183494f,0.000606773f,2.45110f,0.00709986f,2.12161f,1.77038f}; // 3.5T-2G
  BCylParam<float>  fpar4{4.24326f,15.0201f,3.81492f,0.0178712f,0.000656527f,2.45818f,0.00778695f,2.12500f,1.77436f}; // 3.8T-2G
  BCylParam<float>  fpar5{4.21136f,14.8824f,4.01683f,0.0175932f,0.000695541f,2.45311f,0.00813447f,2.11688f,1.76076f}; // 4.0T-2G
  std::string const flds[] = {"2_0T","3_0T","3_5T","3_8T","4_0T"};
  BCylParam<float> const fpars[]{fpar1,fpar2,fpar3,fpar4,fpar5};

  BCylParam<float> const & findPar(std::string fld) {
    auto f = std::find(flds,flds+5,fld);
    if (f-flds>4)    throw cms::Exception("BadParameters") 
      << "Undefined key - " // abort!\n";
      <<"Defined keys are: \"2_0T\" \"3_0T\" \"3_5T\" \"3_8T\" and \"4_0T\"\n";
    return fpars[f-flds]; 
  }
}

TkBfield::TkBfield(std::string fld) : bcyl(findPar(fld)) {}

void TkBfield::getBrfz(float const  * __restrict__ x, float * __restrict__ Brfz)  const {
  float br; float bz;
  float r2=x[0]*x[0]+x[1]*x[1];
  bcyl(x[0]*x[0]+x[1]*x[1],x[2], br, bz); Brfz[0]=std::sqrt(r2)*br; Brfz[1]=0; Brfz[2]=bz;
}

void TkBfield::getBxyz(float const  * __restrict__ x, float * __restrict__ Bxyz)  const {
  float br; float bz;
  float r2=x[0]*x[0]+x[1]*x[1];
  bcyl(r2, x[2], br,bz);
  Bxyz[0]=br*x[0];
  Bxyz[1]=br*x[1];
  Bxyz[2]=bz;
}
//-------------------------------------------------------------------------------

MagneticField *SignedImpactParameter::paramField_ = 0;

SignedImpactParameter::SignedImpactParameter() {
}

SignedImpactParameter::~SignedImpactParameter() {
}

//Signed 3D IP

Measurement1D 
SignedImpactParameter::signedIP3D(const reco::Track &tk, const reco::Vertex &vtx, const reco::Track::Vector jetdir) const {
    if (paramField_ == 0) paramField_ = new OAEParametrizedMagneticField("3_8T");
    reco::TransientTrack ttk(tk,paramField_);
    return IPTools::signedImpactParameter3D(ttk, GlobalVector(jetdir.X(),jetdir.Y(),jetdir.Z()), vtx).second;
}

Measurement1D 
SignedImpactParameter::signedIP3D(const reco::Track &tk, const reco::VertexCompositePtrCandidate &sv, const reco::Track::Vector jetdir) const {
    reco::Vertex::CovarianceMatrix csv; sv.fillVertexCovariance(csv);
    reco::Vertex svtx(sv.vertex(), csv);
    return signedIP3D(tk, svtx, jetdir);
}

//Signed 2D IP

Measurement1D 
SignedImpactParameter::signedIP2D(const reco::Track &tk, const reco::Vertex &vtx, const reco::Track::Vector jetdir) const {
    if (paramField_ == 0) paramField_ = new OAEParametrizedMagneticField("3_8T");
    reco::TransientTrack ttk(tk,paramField_);
    return IPTools::signedTransverseImpactParameter(ttk, GlobalVector(jetdir.X(),jetdir.Y(),jetdir.Z()), vtx).second;
}

Measurement1D 
SignedImpactParameter::signedIP2D(const reco::Track &tk, const reco::VertexCompositePtrCandidate &sv, const reco::Track::Vector jetdir) const {
    reco::Vertex::CovarianceMatrix csv; sv.fillVertexCovariance(csv);
    reco::Vertex svtx(sv.vertex(), csv);
    return signedIP2D(tk, svtx, jetdir);
}


//3D IP

Measurement1D 
SignedImpactParameter::IP3D(const reco::Track &tk, const reco::Vertex &vtx) const {
    if (paramField_ == 0) paramField_ = new OAEParametrizedMagneticField("3_8T");
    reco::TransientTrack ttk(tk,paramField_);
    return IPTools::absoluteImpactParameter3D(ttk,vtx).second;
}

Measurement1D 
SignedImpactParameter::IP3D(const reco::Track &tk, const reco::VertexCompositePtrCandidate &sv) const {
    reco::Vertex::CovarianceMatrix csv; sv.fillVertexCovariance(csv);
    reco::Vertex svtx(sv.vertex(), csv);
    return IP3D(tk, svtx);
}

//2D IP

Measurement1D 
SignedImpactParameter::IP2D(const reco::Track &tk, const reco::Vertex &vtx) const {
    if (paramField_ == 0) paramField_ = new OAEParametrizedMagneticField("3_8T");
    reco::TransientTrack ttk(tk,paramField_);
    return IPTools::absoluteTransverseImpactParameter(ttk,vtx).second;
}

Measurement1D 
SignedImpactParameter::IP2D(const reco::Track &tk, const reco::VertexCompositePtrCandidate &sv) const {
    reco::Vertex::CovarianceMatrix csv; sv.fillVertexCovariance(csv);
    reco::Vertex svtx(sv.vertex(), csv);
    return IP2D(tk, svtx);
}



std::pair<double,double>
SignedImpactParameter::twoTrackChi2(const reco::Track &tk1, const reco::Track &tk2) const {
    if (paramField_ == 0) paramField_ = new OAEParametrizedMagneticField("3_8T");
    std::vector<reco::TransientTrack> ttks;
    ttks.push_back(reco::TransientTrack(tk1,paramField_));
    ttks.push_back(reco::TransientTrack(tk2,paramField_));
    KalmanVertexFitter vtxFitter;
    TransientVertex myVertex = vtxFitter.vertex(ttks);
    return std::make_pair(myVertex.totalChiSquared(),myVertex.degreesOfFreedom());  
}

//Helping functions
std::vector<reco::TransientTrack> SignedImpactParameter::ttrksf(const reco::Track &trkA, const reco::Track &trkB, const reco::Track &trkC, const reco::Track &trkD, int nlep) const {
  std::vector<reco::TransientTrack> ttrks;
 if (paramField_ == 0) paramField_ = new OAEParametrizedMagneticField("3_8T");
 if(nlep==2){
  ttrks.push_back(reco::TransientTrack(trkA,paramField_));
  ttrks.push_back(reco::TransientTrack(trkB,paramField_));
 }else if(nlep==3){
  ttrks.push_back(reco::TransientTrack(trkA,paramField_));
  ttrks.push_back(reco::TransientTrack(trkB,paramField_));
  ttrks.push_back(reco::TransientTrack(trkC,paramField_));
 }else if(nlep==4){
  ttrks.push_back(reco::TransientTrack(trkA,paramField_));
  ttrks.push_back(reco::TransientTrack(trkB,paramField_));
  ttrks.push_back(reco::TransientTrack(trkC,paramField_));
  ttrks.push_back(reco::TransientTrack(trkD,paramField_));
 }
 return ttrks;
}

std::vector<reco::TransientTrack> SignedImpactParameter::ttrksbuthef(const reco::Track &trkA, const reco::Track &trkB, const reco::Track &trkC, const reco::Track &trkD, int nlep, int iptrk) const {
  std::vector<reco::TransientTrack> ttrks;
 if (paramField_ == 0) paramField_ = new OAEParametrizedMagneticField("3_8T");
 if(nlep==3){
  if(iptrk==0){
   ttrks.push_back(reco::TransientTrack(trkB,paramField_));
   ttrks.push_back(reco::TransientTrack(trkC,paramField_));
  }else if(iptrk==1){
   ttrks.push_back(reco::TransientTrack(trkA,paramField_));
   ttrks.push_back(reco::TransientTrack(trkC,paramField_));
  }else if(iptrk==2){
   ttrks.push_back(reco::TransientTrack(trkA,paramField_));
   ttrks.push_back(reco::TransientTrack(trkB,paramField_));
  } 
 }else if(nlep==4){
  if(iptrk==0){
   ttrks.push_back(reco::TransientTrack(trkB,paramField_));
   ttrks.push_back(reco::TransientTrack(trkC,paramField_));
   ttrks.push_back(reco::TransientTrack(trkD,paramField_));
  }else if(iptrk==1){
   ttrks.push_back(reco::TransientTrack(trkA,paramField_));
   ttrks.push_back(reco::TransientTrack(trkC,paramField_));
   ttrks.push_back(reco::TransientTrack(trkD,paramField_));
  }else if(iptrk==2){
   ttrks.push_back(reco::TransientTrack(trkA,paramField_));
   ttrks.push_back(reco::TransientTrack(trkB,paramField_));
   ttrks.push_back(reco::TransientTrack(trkD,paramField_));
  }else if(iptrk==3){
   ttrks.push_back(reco::TransientTrack(trkA,paramField_));
   ttrks.push_back(reco::TransientTrack(trkB,paramField_));
   ttrks.push_back(reco::TransientTrack(trkC,paramField_));
  }
 }
 return ttrks;
}

reco::TransientTrack SignedImpactParameter::thettrkf(const reco::Track &trkA, const reco::Track &trkB, const reco::Track &trkC, const reco::Track &trkD, int nlep, int iptrk) const {
 reco::TransientTrack thettrk; 
 if (paramField_ == 0) paramField_ = new OAEParametrizedMagneticField("3_8T");
 if(nlep==3){
  if(iptrk==0){
   reco::TransientTrack tmpttrk(trkA,paramField_); thettrk = tmpttrk;   
  }else if(iptrk==1){
   reco::TransientTrack tmpttrk(trkB,paramField_); thettrk = tmpttrk;
  }else if(iptrk==2){
   reco::TransientTrack tmpttrk(trkC,paramField_); thettrk = tmpttrk;
  } 
 }else if(nlep==4){
  if(iptrk==0){
   reco::TransientTrack tmpttrk(trkA,paramField_); thettrk = tmpttrk;
  }else if(iptrk==1){
   reco::TransientTrack tmpttrk(trkB,paramField_); thettrk = tmpttrk;
  }else if(iptrk==2){
   reco::TransientTrack tmpttrk(trkC,paramField_); thettrk = tmpttrk;
  }else if(iptrk==3){
   reco::TransientTrack tmpttrk(trkD,paramField_); thettrk = tmpttrk;
  }
 }
 return thettrk;
}

//Variables related to IP
//Of one lepton w.r.t. the PV of the event
std::pair<double,double> SignedImpactParameter::absIP3D(const reco::Track &trk, const reco::Vertex &pv) const {
 if (paramField_ == 0) paramField_ = new OAEParametrizedMagneticField("3_8T");
 reco::TransientTrack ttrk(trk,paramField_);
 Measurement1D aIP3Dtrk = IPTools::absoluteImpactParameter3D(ttrk, pv).second; 
 return std::make_pair(aIP3Dtrk.value(), aIP3Dtrk.error()); 
}
//Of one lepton w.r.t. the PV of the PV of the other leptons only
std::pair<double,double> SignedImpactParameter::absIP3Dtrkpvtrks(const reco::Track &trkA, const reco::Track &trkB, const reco::Track &trkC, const reco::Track &trkD, int nlep, int iptrk) const {
 //Take the transient tracks
 reco::TransientTrack thettrk = thettrkf(trkA,trkB,trkC,trkD,nlep,iptrk); 
 std::vector<reco::TransientTrack> ttrks = ttrksbuthef(trkA,trkB,trkC,trkD,nlep,iptrk);
 //Build new vertex
 KalmanVertexFitter vtxFitter;
 TransientVertex trkspv = vtxFitter.vertex(ttrks);
 //Measure 3DIP
 Measurement1D aIP3Dtrk = IPTools::absoluteImpactParameter3D(thettrk,reco::Vertex(trkspv)).second;
 return std::make_pair(aIP3Dtrk.value(), aIP3Dtrk.error());
} 

//Variables related to chi2
std::pair<double,double> SignedImpactParameter::chi2pvtrks(const reco::Track &trkA, const reco::Track &trkB, const reco::Track &trkC, const reco::Track &trkD, int nlep) const {
 //Take transient tracks
  std::vector<reco::TransientTrack> ttrks = ttrksf(trkA,trkB,trkC,trkD,nlep);
 //Build new vertex
 KalmanVertexFitter vtxFitter;
 TransientVertex trkspv = vtxFitter.vertex(ttrks);
 //Take interested values
 return std::make_pair(trkspv.totalChiSquared(),trkspv.degreesOfFreedom());  
}

Measurement1D
SignedImpactParameter::vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv) const {
    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}

Measurement1D
SignedImpactParameter::vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv) const {
    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}

float SignedImpactParameter::vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv) const {
    reco::Candidate::Vector p = sv.momentum();
    reco::Candidate::Vector d(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
    return p.Unit().Dot(d.Unit());
}
