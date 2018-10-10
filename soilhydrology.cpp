/*****************************************************************
							Soil water
/***************************************************************** */

#if !defined(HYDROLOFY_H)
  #include "soilhydrology.hpp"
#endif

double HYDM::hydmaet(const int& dm, double& rain, double& snowinf, double& avlh2o1, double& avlh2o2,double& avlh2o3, double& pet, double& awcapmm1, double& awcapmm2,double& awcapmm3, double& daze, double root1, double root2,double root3, const int& vt)
	//   avlh2o1 is pstate[vt][I_AVLW1]

{
double aet, r1, r2, r3,aetrate;
double supply;
double w1, w2,w3, tolrt1, tolrt2, tolrt3, tolaet;


const double edpar =5.0;
double gm, ep, def;
double prob, rbar;
double dsm;

// calculation for the monthly values
   w1=w2=w3 =0.0;

   tolrt1=0.0;
   tolrt2=0.0;
   tolrt3=0.0;

   trans1[vt]=0.0;
   trans2[vt]=0.0;
   trans3[vt]=0.0;

   tolaet=0.0;

 // tolrt1= tolrt1[vt] + root1[vt];
   //tolrt2= tolrt2[vt] + root2[vt];
//  tolrt3= tolrt3[vt] + root3[vt];

  if ((rain+snowinf) >= pet) {
  tolaet=pet;
  }
  else {																										//wsr
  		gm = (1.0 - exp(-edpar * (avlh2o1+avlh2o2)/(awcapmm1+awcapmm2)))/(1.0 - exp(-edpar));
      ep = pet /daze;
      def = ep + (awcapmm1+awcapmm2) - (avlh2o1+avlh2o2);
      prob = 1.0 - exp (-0.005*(rain+snowinf));
      if (prob !=0.0) { rbar = (rain +snowinf)/(daze+prob);}
      else {rbar =0.0;}

      if (rbar !=0.0) {
      dsm = rbar*prob*(gm +((1.0-gm)* exp(-ep/rbar)) -exp (-def/rbar)) -(ep*gm);
      }
      else {
      dsm = -ep*gm;
      }
      dsm *=daze;
      tolaet = rain +snowinf -dsm;
      if (tolaet > pet) {tolaet =pet;}
   }

 // 	w1= avlh2o1 + rain + snowinf - lyperc1[vt][dm];
  // w2= avlh2o2 + lyperc1[vt][dm]- lyperc2[vt][dm];
  // w3 = avlh2o3 + lyperc2[vt][dm];
  


  w1 = avlh2o1;
  w2 = avlh2o2;
 w3 = avlh2o3;

   supply = w1 + w2 +w3; // available water for root for transpitation

   //supply = root1*w1/tolrt1 + root2*w2/tolrt2 + root3*w3/tolrt3;
 // supply = w2 + w3; // available water for root for transpitation

   if (supply <= 0.01) {supply=0.01;}

   if (w1<= 0.0) {w1 = 0.01;}
  if (w2<= 0.0) {w2 = 0.01;}
if (w3<= 0.0) {w3 = 0.01;}

   aetrate= 1;
   aet=tolaet;

   if (aet> supply) {
      aet = supply;
      }

																					//	wsr: for peatland, we only consider the eet for upper 30cm, e.g., root1 and root2. root3 is totally saturated.  
  r1 = root1;
   r2 = root2;
  r3 = root3;

 //  r1 = w1/supply;
 //  r2 = w2/supply;
//   r3 = w3/supply;

   trans1[vt] =  r1*aet;
   trans2[vt] =  r2*aet;
trans3[vt] =  r3*aet;



   return aet;
  };

/*************************************************************
  Percolation data
************************************************************/
void HYDM::soil_perco (double& pctsilt, double& pctclay) {

	double psc;

	psc = (pctsilt + pctclay) * 0.01;
	if (psc < 0.01) {psc = 0.01;}

	if (psc <= 0.20 && psc > 0.0){
		kp=5.0;
		whc =0.11; }
	if (psc <= 0.35 && psc > 0.2){
		kp=4.5;
		whc =0.13; }
	if (psc <= 0.55 && psc > 0.35){
		kp=4.0;
		whc =0.15; }
	if (psc <= 0.65 && psc > 0.55){
		kp=3.5;
		whc =0.135; }
  else{
		kp=3.0;
		whc =0.12; }

	k1= kp;
	k2= 4.0;
   k3= 4.0; // think about this
// kmoss1=5.0;
//   kmoss2=6.0;
};
