// root -l BTagReshaping.C+
// reshape(1.2, 30, 0.3, 5)

#include "../interface/BTagReshaping.h"

BTagShapeInterface sh("../bin/csvdiscr.root",0,0);

float reshape(float eta, float pt, float csv, int flav)
{
 return  sh.reshape(eta, pt, csv, flav);

}

void BTagReshaping()
{
 
}
