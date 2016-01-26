// 
// TFitParticleEt::
// --------------------
//
// Particle with a special parametrization useful in hadron 
// colliders (3 free parameters (Et, Eta, Phi). The parametrization is 
// chosen as follows:
//
// p = (EtCosPhi, EtSinPhi, EtSinhEta)
// E =  EtCoshEta
//

#include <iostream>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEt.h"
#include "TMath.h"
#include <cmath>

//----------------
// Constructor --
//----------------
TFitParticleEt::TFitParticleEt()
  :TAbsFitParticle()  
{
  init(0, 0);
}

TFitParticleEt::TFitParticleEt( const TFitParticleEt& fitParticle )
  :TAbsFitParticle( fitParticle.GetName(), fitParticle.GetTitle() )
{

  _nPar = fitParticle._nPar;
  _u1 = fitParticle._u1;
  _u2 = fitParticle._u2;
  _u3 = fitParticle._u3;
  _covMatrix.ResizeTo(  fitParticle._covMatrix );
  _covMatrix = fitParticle._covMatrix;
  _iniparameters.ResizeTo( fitParticle._iniparameters );
  _iniparameters = fitParticle._iniparameters;
  _parameters.ResizeTo( fitParticle._parameters );
  _parameters = fitParticle._parameters;
  _pini = fitParticle._pini;
  _pcurr = fitParticle._pcurr;

}

TFitParticleEt::TFitParticleEt(TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle()  
{
  init(pini, theCovMatrix);
}

TFitParticleEt::TFitParticleEt(const TString &name, const TString &title, 
			   TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle(name, title)  
{
  init(pini, theCovMatrix);
}

TAbsFitParticle* TFitParticleEt::clone(const TString& newname ) const {
  // Returns a copy of itself
  
  TAbsFitParticle* myclone = new TFitParticleEt( *this );
  if ( newname.Length() > 0 ) myclone->SetName(newname);
  return myclone;

}

//--------------
// Destructor --
//--------------
TFitParticleEt::~TFitParticleEt() {

}

//--------------
// Operations --
//--------------
void TFitParticleEt::init(TLorentzVector* pini, const TMatrixD* theCovMatrix ) {

  _nPar = 1;
  setIni4Vec(pini);
  setCovMatrix(theCovMatrix);

}

TLorentzVector* TFitParticleEt::calc4Vec( const TMatrixD* params ) {
  // Calculates a 4vector corresponding to the given
  // parameter values

  if (params == 0) {
    return 0;
  }

  if ( params->GetNcols() != 1 || params->GetNrows() !=_nPar ) {
    edm::LogError ("WrongMatrixSize")
      << GetName() << "::calc4Vec - Parameter matrix has wrong size.";
    return 0;
  }

  Double_t et = (*params)(0,0);
  
  Double_t X = et*TMath::Cos(_phi);
  Double_t Y = et*TMath::Sin(_phi);
  Double_t Z = et*TMath::SinH(_eta);
  Double_t E = et*TMath::CosH(_eta);
		
  TLorentzVector* vec = new TLorentzVector( X, Y, Z, E );
  return vec;

}

void TFitParticleEt::setIni4Vec(const TLorentzVector* pini) {
  // Set the initial 4vector. Will also set the 
  // inital parameter values 

  if (pini == 0) {

    _u1.SetXYZ(0., 0., 0.);
    _u3.SetXYZ(0., 0., 0.);
    _u2.SetXYZ(0., 0., 0.);
    _pini.SetXYZT(0., 0., 0., 0.);
    _pcurr = _pini;

    _iniparameters.ResizeTo(_nPar,1);
    _iniparameters(0,0) = 0.;
    
    _parameters.ResizeTo(_nPar,1);
    _parameters(0,0) = 0.;
    
    _eta=0;
    _phi=0;
    
  } else {
    
    Double_t et = pini->E()*std::fabs(sin(pini->Theta()));
    Double_t eta = pini->Eta();
    Double_t phi = pini->Phi();
    
    _pini = (*pini);
    _pcurr = _pini;
    
    _iniparameters.ResizeTo(_nPar,1);
    _iniparameters(0,0) = et;
    _eta=eta;
    _phi=phi;
    
    _parameters.ResizeTo(_nPar,1);
    _parameters = _iniparameters;

    _u1.SetXYZ( TMath::Cos(phi), TMath::Sin(phi), 0.); // the base vector of Et
    _u2.SetXYZ( -1.*TMath::Cos(phi)*TMath::TanH(eta), -1.*TMath::Sin(phi)*TMath::TanH(eta), 1./TMath::CosH(eta) );// the base vector of Eta ( same as the base vector for Theta)
    _u3.SetXYZ( -1.*TMath::Sin(phi), TMath::Cos(phi), 0. );// the base vector of Phi

  }

}

TMatrixD* TFitParticleEt::getDerivative() {
  // returns derivative dP/dy with P=(p,E) and y=(et) 
  // the free parameters of the fit. The columns of the matrix contain 
  // (dP/d(et)).

  TMatrixD* DerivativeMatrix = new TMatrixD(4,1);
  (*DerivativeMatrix) *= 0.;

  //1st column: dP/d(et)
  (*DerivativeMatrix)(0,0) = TMath::Cos(_phi);
  (*DerivativeMatrix)(1,0) = TMath::Sin(_phi);
  (*DerivativeMatrix)(2,0) = TMath::SinH(_eta);
  (*DerivativeMatrix)(3,0) = TMath::CosH(_eta);

  return DerivativeMatrix;

}

TMatrixD* TFitParticleEt::transform(const TLorentzVector& vec) {
  // Returns the parameters corresponding to the given 
  // 4vector

  // retrieve parameters
  TMatrixD* tparams = new TMatrixD( _nPar, 1 );
  (*tparams)(0,0) = vec.E()*std::fabs(sin(vec.Theta()));

  return tparams;

}
