/*
Double _t triggerZnunuCurve(met) {

TF1  f  ("f","1. / (1. +  ( TMath::E()^( 0.059486 * ( 123.27 - x))))", 0., 999999.);
 return f;


}

"1. / (1. +  ( TMath::E()^( 0.059486 * ( 123.27 - x))))
*/

#include "TMath.h"

namespace  TriggerZnunuCurve {

 double trigMet(double *x, double *p){
    return 1. / (1. +  ( TMath::Exp( 0.059486 * ( 123.27 - x[0] ))));
  }


}
