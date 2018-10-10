
/*********************************************************
**********************************************************
hydrology.HPP -object describing Hydrology models
**********************************************************
*********************************************************/
class HYDM {

	public:

/******************************************************
			public Functions
******************************************************/
   double hydmaet(const int& dm, double& rain, double& snowinf, double& avlh2o1, double& avlh2o2, double& avlh2o3, double& pet, double& awcapmm1, double& awcapmm2, double& awcapmm3,double& daze, double root1, double root2,double root3, const int& vt);
 	void soil_perco (double& pctsilt, double& pctclay);

/**************************************************************
			Public Variables
**************************************************************/

// three-box soil water variables

		double root1[NVT];
		double root2[NVT];
		double root3[NVT];

      double rtratio1[NVT];
      double rtratio2[NVT];
      double rtratio3[NVT];

      double rtsface;
      double rtindex1[NVT];
      double rtindex2[NVT];
      double rtindex3[NVT];

     	double rootmass[NVT][CYCLE];

      double rteffcnt;
      double kwt1[NVT];
      double kwt2[NVT];
      double kwt3[NVT];

      double prkwt1, prkwt2,prkwt3;

		double trans1[NVT], trans2[NVT],trans3[NVT];

    	double fpc[NVT];
      double fpcmax[NVT];
	  	double sla[NVT];
      double cov[NVT];

// allocation variables
      double alpha;
      double sapallo;
      double alleaf;

		double k1, k2, k3,kmoss1,kmoss2,kp,whc;
		double lyperc[NVT][CYCLE];
      double lyperc1[NVT][CYCLE];
      double lyperc2[NVT][CYCLE];

		double yrlyperc[NVT];
      double yrlyperc1[NVT];
      double yrlyperc2[NVT];

      double wscal;

};

