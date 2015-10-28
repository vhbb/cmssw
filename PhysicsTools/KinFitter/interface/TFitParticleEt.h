#ifndef TFitParticleEt_hh
#define TFitParticleEt_hh


#include "PhysicsTools/KinFitter/interface/TAbsFitParticle.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"


class TFitParticleEt: public TAbsFitParticle {

public :

  TFitParticleEt();
  TFitParticleEt( const TFitParticleEt& fitParticle );
  TFitParticleEt(TLorentzVector* pini, const TMatrixD* theCovMatrix);
  TFitParticleEt(const TString &name, const TString &title, 
	       TLorentzVector* pini,
	       const TMatrixD* theCovMatrix);
  virtual ~TFitParticleEt();
  virtual TAbsFitParticle* clone(const TString& newname = "" ) const;

  // returns derivative dP/dy with P=(p,E) and y=(et) 
  // the free parameters of the fit. The columns of the matrix contain 
  // (dP/d(et)).
  virtual TMatrixD* getDerivative();
  virtual TMatrixD* transform(const TLorentzVector& vec);
  virtual void setIni4Vec(const TLorentzVector* pini);
  virtual TLorentzVector* calc4Vec( const TMatrixD* params );
  
protected :
  
  void init(TLorentzVector* pini, const TMatrixD* theCovMatrix);
  
  
private:

  double _eta, _phi;
    
};
  
#endif
  
  
