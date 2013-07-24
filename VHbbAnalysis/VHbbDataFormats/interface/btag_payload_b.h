namespace beff{

float ptmin[] = {30, 210, 260, 320, 400, 500};
float ptmax[] = {210,260, 320, 400, 500, 670};

size_t bins=6;

/*
From ttbar for 30-200
Working Points:
CSV_L = 1.01 +/- 0.04
CSV_M = 0.97 +/- 0.04
CSV_T = 0.96 +/- 0.04
*/

// Tagger: CSVL within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
float CSVL_SFb(float x) { if(x < 210) return 1.01 ; return 1.02658*((1.+(0.0195388*x))/(1.+(0.0209145*x))); }
float CSVL_SFb_error[] = {
 0.04,
 0.0257175,
 0.026424,
 0.0264928,
 0.0315127,
 0.030734,
 0.0438259 };

// Tagger: CSVM within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
float CSVM_SFb(float x) {  if(x < 210) return 0.97 ; return  0.6981*((1.+(0.414063*x))/(1.+(0.300155*x)));}
 float CSVM_SFb_error[] = {
 0.04,
 0.0409675,
 0.0420284,
 0.0541299,
 0.0578761,
 0.0655432 };

// Tagger: CSVT within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
float CSVT_SFb(float x) {  if(x < 210) return 0.96 ; return   0.901615*((1.+(0.552628*x))/(1.+(0.547195*x))); }
float CSVT_SFb_error[] = {
 0.04,
 0.0653621,
 0.0712586,
 0.094589,
 0.0777011,
 0.0866563 };
 
}
