#ifndef MULTITHRESHOLDEFFICIENCY_H
#define MULTITHRESHOLDEFFICIENCY_H
#include <math.h>
#include <iostream>
#include <vector>
using namespace std; 
class MultiThresholdEfficiency 
{
 public:
   MultiThresholdEfficiency(unsigned int nOfThresholds) : thresholds(nOfThresholds) {}
   
   //A vector with, for each object the efficiency of passing threshold N while notpassing N-1 
   template <class Filter> float weight(vector<vector<float> > input);

 private:
  unsigned int thresholds;

};


// NB here sorting is in opposite direction than for btag 0 is tigther , 1 is looser
class Trigger1High2Loose
{
 public:
  static bool filter(std::vector<int> t)
  {
    return t[0] >= 1 && t[1] >= 2;
}
};




template <class Filter> float MultiThresholdEfficiency::weight(vector<vector<float> > objects)
{
 unsigned int nobjects=objects.size();
 std::vector<unsigned int> comb(objects.size());
 for(unsigned int i=0;i < objects.size(); i++) comb[i]=0;
 unsigned int idx=0;
 unsigned int max=thresholds+1; // 
 float p=0,tot=0;
 if(objects.size()==0) return 0.;
 while(comb[objects.size()-1] < max)
 {
//std::cout << std::endl << "New comb" << std::endl;
// for(int i=0;i < objects.size(); i++) {std::cout << comb[i] << " ";}
// std::cout << std::endl;
   std::vector<int> pass;
   for(unsigned int j=0;j<thresholds;j++) pass.push_back(0);

   float thisCombinationProbability=1.;
//   std::cout << "OBJ Probs: ";
   for(size_t j=0;j<nobjects;j++) // loop on objects
    {
       float cumulative=1.;
       for(unsigned int n=0;n<comb[j];n++) cumulative*= (1.-objects[j][n]); // 10 20 70  10/100,20/90   90/100*70/90 
       thisCombinationProbability*=cumulative;
       if(comb[j]< thresholds) {
         thisCombinationProbability*=(objects[j])[comb[j]];     
//         std::cout << cumulative*(objects[j])[comb[j]] << " "; 
       }
/*       else
       { 
              std::cout << cumulative << " " ;
       }*/

 
     for(size_t k=0;k< thresholds; k++ ) // loop on threshold
      {
       bool passed = ( k >= comb[j] ) ;  //k=0, is the tightest, passed only if comb[j] = 0, k=1 pass if comb 0 or 1
       if(passed) pass[k]++;
      }
  }
//   std::cout << endl;
   if(Filter::filter(pass)) p+=thisCombinationProbability;
   tot+=thisCombinationProbability;
//   std::cout << thisCombinationProbability << " " << p << " " << tot<< std::endl;
  while (comb[idx] == max -1  && idx+1 < objects.size()) idx++; // find first object for which we did not already test all configs 
 // next combination
  comb[idx]++;  // test a new config for that jet
  for(unsigned int i=0;i<idx;i++) { comb[i]=0; } // reset the tested configs for all previous objects
  idx=0;
 }
  if( fabs(tot-1.) > 0.01 )
   {
     std::cout << "ERROR, total must be one " << tot << std::endl;
   }
  return p;
}

#endif

