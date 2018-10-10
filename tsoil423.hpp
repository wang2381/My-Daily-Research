/* **************************************************************
*****************************************************************
TSOIL423.HPP - object describing characteristics of soil used by
	       the Terrestrial Ecosystem Model (TEM)
            - modified by DWK on 20000102
*****************************************************************
************************************************************** */

#if !defined(BIOMS423_H)
  #include "bioms423.hpp"  // Tsoil4 uses Biomass class
#endif

// Tsoil4 also uses the global constants CYCLE, MAXRTIME, MISSING
//   and NUMVEG

#if !defined(TSOIL423_H)
#define TSOIL423_H

class Tsoil4 {

  public:

     Tsoil4(void);

/* **************************************************************
		 Public Functions
************************************************************** */

     void getecd (ofstream& rflog1);
     void getecd (char ecd[80]);
     void getrootz(ofstream& rflog1);
     void getrootz(char ecd[80]);
     void lake(double& tair, double& prec, double& rain, double& snowfall,
	       double& pet, double& eet, int& dm, const int& vt);

	   // added for hydrology model
//     void percol_hydm(double& rain, double& snowinf, double& trans1, double& trans2, double& avlh2o1, double& avlh2o2, double& lyperc1, int& dm);
     void percol_hydm(double& tair, double& rain, double& snowinf, double& trans1, double& trans2,double& trans3, double& avlh2o1, double& avlh2o2,double& avlh2o3, double& lyperc1, double& lyperc2,const int& dm, const int& vt);

     void percol(double& rain, double& snowinf, double& eet, double& avlh2o, const int& dm, const int& vt);
     double rrunoff(const double& rgrndh2o, const double& rperc);
     void showecd(void);
     double snowmelt(const double& elev, const double& tair, const double& prevtair, double& snowpack);
     double srunoff(const double& elev, const double& tair, const double& prevtair, const double& prev2tair,
		    const double& sgrndh2o, const double& sperc);
 // added for hydrology model
     void xtext(int& ez, double& pctsilt, double& pctclay, const int& vt);
     void hydm_xtext(int& ez, double& pctsilt, double& pctclay, const int& vt);
 // added for hydrology model regarding moisture freezing
	 double freeze2(double& pfstate1,double& t1,double& t2,const int& vt,const int& dm);
	 double freeze1(double& pfstate2,double& t1,double& t2,const int& vt,const int& dm);
	 double freeze3(double& pfstate3,double& t1,double& t2,const int& vt,const int& dm);
	 void crit(double& t1,double& t2,const int& vt,const int& dm);
	 void wetlandTAWA(double& VSM1,double& VSM2,double& box1,double& box2, double& porosity1,double& porosity2, const int& vt, const int& dm);
	 void cal_methane(const int& vt, const int& dm, double& watertable,double& npp,double& lowb);
	 double EffectOMD(const double& depthz, const double& rootd, const double& lowb2);
	 double EffectOM(const double& maxfresh, const double& vegNPP);
	 double InterpolatST(const double& dst1, const double& dst2,const double& dst3,
	const double& dst4, const double& dst5, const double& dst6,  double x);
	 double EffectST( double soilt);
	 double EffectPH(const double& soilph);
	 double EffectRX(const double& ehl);
	 double RedoxP(const double& watertable, const double& depthz,const double& lowb2 ,const int& vt);
	 double MethanePR(double fsom, double fcdis, double fmst, double fph, double frx);
	 void cal_diff(double& VSM1,double& VSM2,double& box1,double& box2,const int& vt, const int& dm, double& watertable,const double& tair,double& lowb);
	 double Diffusivity(const double& sand, const double& silt, const double& clay, int satyorn);
	 double BCH4flxandCon(const double& lowb, const double& watertable);
         double EffectMC(const double& ch4con, const double& kc);
         void oxygenC(const double& afp, const double& lowb );
         double  EffectOXY(const double& oxyc, const double& ko);
         double EffectST2(const double& soilt, const double& och4q10, const double& oxiref);
         double EffectMC2(const double& ch4con, const double& kc);
         double EffectSM(const double& mvmin, const double& mvmax, const double& mvopt, const double& soilm);
         double EffectRX2(const double& ehl);
         double MethaneOR(const double& omax, double& fmc, double& foxy, double& fst, double& fsm, double& frx, double& fct );
         double CCH4flxandCon(const double& lowb, const double& watertable, const double& soilt);
         void bulkdensity(double& delta_c, double& sum_depth);
         void getdata_noname(double& sum_deptht);
         void getdata_swanson(double& sum_deptht);
         void getdata_horse(double& sum_deptht);
         void getdata_gasfield(double& sum_deptht);
/* **************************************************************
		 Public Variables
************************************************************** */

     int text;              // soil texture (categorical data)
     int wsoil;             // wetland soil type designation
			    // (categorical data)

     double pctsand;        // %sand
     double pctsilt;        // %silt
     double pctclay;        // %clay
     double psiplusc;       // proportion silt and clay

     double dpwbox1[NVT];        // depth (m) of first-layer soil water box
     double dpwbox2[NVT];        // depth (m) of first-layer soil water box
     double dpwbox3[NVT];        // depth (m) of first-layer soil water box

     double awcapmm1[NVT];        // available water capacity (mm)
	  double awcapmm2[NVT];        // available water capacity (mm)
	  double awcapmm3[NVT];        // available water capacity (mm)


     double awcapmm[NVT];        // available water capacity (mm)
     double fldcap[NVT];         // volume of water at field capacity (mm)
     double wiltpt[NVT];         // volume of water at wilting point  (mm)
     double totpor[NVT];         // volume of total pore space        (mm)

     double fldcap1[NVT];         // volume of water at field capacity (mm)
	  double wiltpt1[NVT];         // volume of water at wilting point  (mm)
	  double totpor1[NVT];         // volume of total pore space        (mm)

	  double fldcap2[NVT];         // volume of water at field capacity (mm)
	  double wiltpt2[NVT];         // volume of water at wilting point  (mm)
	  double totpor2[NVT];

	  double fldcap3[NVT];         // volume of water at field capacity (mm)
	  double wiltpt3[NVT];         // volume of water at wilting point  (mm)
	  double totpor3[NVT];

	  double avlh2o1[NVT][CYCLE];  // available water (mm)
	  double avlh2o2[NVT][CYCLE];  // available water (mm)
     double avlh2o3[NVT][CYCLE];  // available water (mm)

	  double moist1[NVT][CYCLE];   // soil moisture (mm)
	  double moist2[NVT][CYCLE];   // soil moisture (mm)
	  double moist3[NVT][CYCLE];   // soil moisture (mm)
	  double storage2[NVT][CYCLE];
	  double storage3[NVT][CYCLE];


     double avlh2o[NVT][CYCLE];  // available water (mm)
     double moist[NVT][CYCLE];   // soil moisture (mm)
     double pcfc[NVT][CYCLE];    // soil h2o as %field capacity
     double pctp[NVT][CYCLE];    // soil h2o as %total pore space
     double vsm[NVT][CYCLE];     // soil h2o as %rooting depth

     double rperc[NVT][CYCLE];      // rain excess (mm)
     double rgrndh2o[NVT][CYCLE];   // rain runoff storage (mm)
     double rrun[NVT][CYCLE];       // rain runoff (mm)
     double snowpack[NVT][CYCLE];   // snowpack (mm)
     double snowinf[NVT][CYCLE];    // snow melt infiltration (mm)
     double sperc[NVT][CYCLE];      // snow melt excess (mm)
     double sgrndh2o[NVT][CYCLE];   // snowmelt runoff storage (mm)
     double srun[NVT][CYCLE];       // snow runoff (mm)
     double h2oyld[NVT][CYCLE];     // water yield (mm)

	  double pcfc1[NVT][CYCLE];    // soil h2o as %field capacity
	  double pctp1[NVT][CYCLE];    // soil h2o as %total pore space
     double vsm1[NVT][CYCLE];     // soil h2o as %rooting depth

     double pcfc2[NVT][CYCLE];    // soil h2o as %field capacity
	  double pctp2[NVT][CYCLE];    // soil h2o as %total pore space
     double vsm2[NVT][CYCLE];     // soil h2o as %rooting depth

     double pcfc3[NVT][CYCLE];    // soil h2o as %field capacity
	  double pctp3[NVT][CYCLE];    // soil h2o as %total pore space
     double vsm3[NVT][CYCLE];     // soil h2o as %rooting depth

     double rperc1[NVT][CYCLE];
	  double rperc2[NVT][CYCLE];
	  double rperc3[NVT][CYCLE];

	  double sperc1[NVT][CYCLE];
     double sperc2[NVT][CYCLE];
     double sperc3[NVT][CYCLE];
     double init_carbon[1000];																						//wsr: add for soil carbon accumulating
     int callflag;
     double delta_c;
     double delta_c1;
     double sum_depth;
     int rhflag;    // = 1 when in transient simulation for peatland
     int init_carbonflag;
     double carbonthread;
     

     double yravlh2o[NVT];
     double yrrgrndh2o[NVT];
     double yrrperc[NVT];
     double yrrrun[NVT];
     double yrsnowpack[NVT];
     double yrsnowinf[NVT];
     double yrsperc[NVT];
     double yrsgrndh2o[NVT];
     double yrsrun[NVT];
     double yrh2oyld[NVT];

     double yravlh2o1[NVT];
     double yravlh2o2[NVT];
     double yravlh2o3[NVT];

     double yrsmoist1[NVT];
     double yrsmoist2[NVT];
     double yrsmoist3[NVT];

     double yrvsm1[NVT];
	  double yrvsm2[NVT];
	  double yrvsm3[NVT];

	  double yrpctp1[NVT];
	  double yrpctp2[NVT];
	  double yrpctp3[NVT];


     double yrvsm[NVT];
     double yrsmoist[NVT];
     double yrpctp[NVT];

     double meanvsm;             // mean annual volumetric soil moisture (%)
	 double VSM1;
	 double VSM2;

     Biomass org[NVT][CYCLE];         // soil organic matter
     double availn[NVT][CYCLE];       // available nitrogen

     // total nitrogen input to soils (g N / (sq. meter * month))
     double ninput[NVT][CYCLE];

     // total nitrogen input to soils (g N / (sq. meter * month))
     double tninput[MAXRTIME][CYCLE];

     long ninyear[MAXRTIME];

     // total nitrogen lost from soils (g N / (sq. meter * month))
     double nlost[NVT][CYCLE];

     int tndepflag;             // flag for transient N deposition data

     double yrorgc[NVT];
     double yrorgn[NVT];
     double yrc2n[NVT];
     double yravln[NVT];
     double yrnin[NVT];
     double yrnlost[NVT];
	 //added for hydrology model regarding moisture freezing

	 
	 double tolfreeze;
	 double temp;
	 double tair1[CYCLE];
	 double unsatthetaWL[1000];
	 double phi;
	 double watertable;
	 double fcdist;
	 double fsomt;
	 double fmstt;
	 double fpht;
	 double frxt;
	 double ch4pro[1000];
	 double ch4ratesat[1000];
	 double totalch4;
	 double ehlt;
	 double intersoilt[1000];
         double fdfs[1000];
         double ch4flx;
         double ch4flxfinal;
         double fmct;
         double foxyt;
         double fstt;
         double fsmt;
         double frxt2;
         double fctt;
         double ch4oxirate[1000];
         double ch4con_b[1000];
         double oxygenc[1000];
         double rot;
         double ratiot;
         
         
	 //double pfstate[NVT][NUMEQ];
/* **************************************************************
		 Public Parameters
************************************************************** */


     double pctpor[NVT];           //porosity of soil (%soil volume)
     double pctpora;
     double pctporb;
     double pcfldcap[NVT];         //field capacity (%soil volume)
     double fldcapa;
     double fldcapb;
     double pcwiltpt[NVT];         //wilting point (%soil volume)
     double wiltpta;
     double wiltptb;
     double rootz[NVT];            //effective rooting depth (m)
     double rootza[MAXCMNT];
     double rootzb[MAXCMNT];
     double rootzc[MAXCMNT];
     double minrootz[MAXCMNT];
     double rootmx; // added for 2-box


// Calibration Parameters

//total nitrogen input into soil (g N / (square meter * month))
     double nintot[MAXCMNT];

//proportion of available nitrogen lost (g N / (square meter * month))
     double nloss[MAXCMNT];
     double nloss_p[MAXCMNT];
     double nloss_temp[MAXCMNT];

};

#endif






