#ifndef MULTITHRESHOLDEFFICIENCY_H
#define MULTITHRESHOLDEFFICIENCY_H
#include <math.h>
#include <iostream>
#include <vector>
using namespace std; 
class MultiThresholdEfficiency 
{
 public:
   MultiThresholdEfficiency(int nOfThresholds) : thresholds(nOfThresholds) {}
   
   //thresholds are from tigther to looser
   virtual bool filter(vector<int> numberOfPassPerThreshold);

   //A vector with, for each object the efficiency of passing threshold N while notpassing N-1 
   float weight(vector<vector<float> > input);

 private:
   int thresholds;

};


bool MultiThresholdEfficiency::filter(std::vector<int> t)
{
 // One high threshold and two low threshold
 return t[0] >= 1 && t[1] >= 2;
}



float MultiThresholdEfficiency::weight(vector<vector<float> > objects)
{
 int nobjects=objects.size();
 std::vector<int> comb(objects.size());
 for(int i=0;i < objects.size(); i++) comb[i]=0;
 int idx=0;
 int max=thresholds+1; // 
 float p=0,tot=0;
 if(objects.size()==0) return 0.;
 while(comb[objects.size()-1] < max)
 {
//std::cout << std::endl << "New comb" << std::endl;
// for(int i=0;i < objects.size(); i++) {std::cout << comb[i] << " ";}
// std::cout << std::endl;
   std::vector<int> pass;
   for(int j=0;j<thresholds;j++) pass.push_back(0);

   float thisCombinationProbability=1.;
//   std::cout << "OBJ Probs: ";
   for(size_t j=0;j<nobjects;j++) // loop on objects
    {
       float cumulative=1.;
       for(int n=0;n<comb[j];n++) cumulative*= (1.-objects[j][n]); // 10 20 70  10/100,20/90   90/100*70/90 
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
   if(filter(pass)) p+=thisCombinationProbability;
   tot+=thisCombinationProbability;
//   std::cout << thisCombinationProbability << " " << p << " " << tot<< std::endl;
  while (comb[idx] == max -1  && idx+1 < objects.size()) idx++; // find first object for which we did not already test all configs 
 // next combination
  comb[idx]++;  // test a new config for that jet
  for(int i=0;i<idx;i++) { comb[i]=0; } // reset the tested configs for all previous objects
  idx=0;
 }
  if( fabs(tot-1.) > 0.01 )
   {
     std::cout << "ERROR, total must be one " << tot << std::endl;
   }
  return p;
}

#endif

