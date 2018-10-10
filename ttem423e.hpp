/* *************************************************************
****************************************************************
TTEM423.HPP - Terrestrial Ecosystem Model Version 4.2 modified to
           use MPIatms class instead of the Atmosphere class
****************************************************************

Modifications:

19991028 - DWK added bug fix to cpp file
20000102 - DWK added compiler directives

20000715 -QZ added hydrology module for hydrology model
****************************************************************

References:

VERSION 4.1

Tian, H., J.M. Melillo, D.W. Kicklighter, A.D. McGuire and J. Helfrich.  1999.
  The sensitvity of terrestrial carbon storage to historical climate variability
  and atmospheric CO2 in the United States.  Tellus 51B: 414-452.


****************************************************************
************************************************************** */

// Global constants
const int NUMWEQ = 17+29;
const int NUMEEQ = 75+10;
const int NUMEQ = NUMWEQ + NUMEEQ; // added for hydrology model by QZ  July/15/2000, 29 new potential output variables, 10 for STM
const int MAXWSTAT = 7+12; // 12 for 3-box pools  13 new water pools added
const int MAXESTAT = 18;
const int MAXSTATE = MAXWSTAT + MAXESTAT;

const int ACCEPT = 0;
const int REJECT = 1;

// Objects describing basic components of the ecosystem

#if !defined(BIOMS423_H)
  #include "bioms423.hpp"   // TTEM uses Biomass class
#endif

//Objects describing the structure of the ecosystem

#if !defined(HYDROLOGY423_H)
  #include "hydrology423.cpp" // TTEM inherits Odeint4 class
#endif

#if !defined(HYDROLOGY_H)
  #include "soilhydrology.cpp" // TTEM inherits Odeint4 class
#endif

#if !defined(HATMS423_H)
  #include "hatms423.cpp"   // TTEM uses MPIatms class
#endif
#if !defined(TVEG423_H)
  #include "tveg423e.cpp"   // TTEM uses Tveg4 class
#endif
#if !defined(TSOIL423_H)
  #include "tsoil423.cpp"  // TTEM uses Tsoil4 class
#endif
#if !defined(TMCRB423_H)
  #include "tmcrb423.cpp"  // TTEM uses Tmicrobe4 class
#endif

//Objects describing the effects of human activities on the ecosystem

#if !defined(HUMNACT423_H)
  #include "humnact423.cpp" // TTEM uses Humnact class
#endif

#if !defined(QSOILTEMP_H)
   #include "qsoiltemp.cpp"     // added for Soil thermal class
#endif

//Adaptive Runge-Kunge Integrator Module

#if !defined(ODEINT423_H)
  #include "odeint423.cpp" // TTEM inherits Odeint4 class
#endif

#if !defined(TTEM423_H)
#define TTEM423_H

class TTEM : public Odeint4 {

  public:

     TTEM();

     enum temkey { I_VEGC,     I_SOLC,     I_TOTC,     I_VEGN,     I_STRN,
                   I_STON,     I_SOLN,     I_AVLN,     I_AGPRDC,   I_AGPRDN,
                   I_PROD10C,  I_PROD10N,  I_PROD100C, I_PROD100N, I_TOTPRDC,
                   I_TOTPRDN,  I_TOTEC,    I_TOTGC,

                   I_AVLW,     I_RGRW,     I_SNWPCK,   I_SGRW,     I_SM,
                   I_PCTP,     I_VSM,

                   I_AVLW1,    I_AVLW2,    I_AVLW3,    I_SM1,      I_SM2,
                   I_SM3,      I_PCTP1,    I_PCTP2,    I_PCTP3,    I_VSM1,
                   I_VSM2,     I_VSM3,

                   I_INGPP,    I_GPP,      I_INNPP,    I_NPP,      I_GPR,
                   I_RVMNT,    I_RVGRW,    I_LTRC,     I_RH,       I_NEP,

                   I_NINP,     I_INNUP,    I_VNUP,     I_VSUP,     I_VLUP,
                   I_VNMBL,    I_VNRSRB,   I_LTRN,     I_MNUP,     I_NMIN,
                   I_NLST,

                   I_RAIN,     I_RPERC,    I_RRUN,     I_SNWFAL,   I_SNWINF,
                   I_SPERC,    I_SRUN,     I_PET,      I_EET,      I_WYLD,

                   I_RPERC1,   I_RPERC2,   I_RPERC3,   I_SPERC1,   I_SPERC2,
                   I_SPERC3,   I_LYPERC,   I_LYPERC1,  I_LYPERC2,  I_MOIST1,
                   I_MOIST2,   I_MOIST3,

                   I_TRANS,    I_EVAP,     I_SEVAP,    I_SNOWSUB,  I_SUBCAN,

                   I_TSOIL,    I_DST5,     I_DST10,    I_DST20,    I_DST50,
                   I_DST100,   I_DST200,   I_FRONTD,   I_THAWBE,   I_THAWEND,


                   I_UNRMLF,   I_LEAF,     I_FPC,      I_LAI,

                   I_CNVRTC,   I_CNVRTN,   I_SCNVRTC,  I_SCNVRTN,  I_NVRTNT,
                   I_NSRTNT,   I_NRETNT,   I_SLASHC,   I_SLASHN,   I_PRDF10C,
                   I_PRDF10N,  I_PRDF100C, I_PRDF100N,

                   I_AGNPPC,   I_AGNPPN,   I_AGFPRDC,  I_AGFPRDN,  I_AGLTRC,
                   I_AGLTRN,   I_AGFRTN,

                   I_TOTFPRDC, I_TOTFPRDN, I_AGPRDFC,  I_AGPRDFN,  I_PRD10FC,
                   I_PRD10FN,  I_PRD100FC, I_PRD100FN, I_TOTPRDFC, I_TOTPRDFN,

                   I_TOTNPP,   I_CFLX };

                   // added I_TRANS,    I_EVAP, SEVAP,I_SNOWSUB for hydrology model by QZ
/* **************************************************************
			Public Functions
************************************************************** */

     // virtual functions for TEM::adapt()
//     virtual int boundcon(double ptstate[], double err[], double& ptol);
     virtual int boundcon(double ptstate[NVT][NUMEQ], double err[NVT][NUMEQ], double& ptol, const int& vt);
//     virtual void delta(const int& dm, double pstate[], double pdstate[]);
     virtual void delta(const int& dm, double pstate[NVT][NUMEQ], double pdstate[NVT][NUMEQ], const int& vt);

//     void deltaxclm(const int& dcmnt, const double& pcfldcap, const int& dm);
     void deltaxclm(const int& dcmnt, const double& pcfldcap, const int& dm, const int& vt);
//     int ecdqc(const int& dcmnt);
      int ecdqc(const int& dcmnt);
     //     int equilibrium(const int& itype, double& tol);
     int equilibrium(const int& itype, double& tol, const int& vt);
     void getco2(void);
//     void getenviron(const int& dm);
     void getenviron(const int& dm, const int& vt);
     void getsitecd(ofstream& rflog1);
     void getsitecd(const int& numcmnt, ofstream& rflog1);
     void getsitecd_p(const int& numcmnt, ofstream& rflog1);
     void getsitecd(char ecd[80]);
     void getsitecd(const int& dv, char ecd[80]);
     void getsitecd_p(const int& dv, char ecd[80]);
     void transfer(void);
     void recover(void);

//     void setELMNTecd(const int& kdinflg, const int& dcmnt, const double& psiplusc);
      void setELMNTecd(const int& kdinflg, const int& dcmnt, const double& psiplusc, const int& vt);

//     void ECDsetELMNTstate(const int& dcmnt, const double& psiplusc);
     void ECDsetELMNTstate(const int& dcmnt, const double& psiplusc,const int& vt);
// modified for hydrology model by QZ
//     void setELMNTevap(const int& stateflg, const int& dcmnt,
  //                       double pet[CYCLE], double tair[CYCLE]);

 //    void setELMNTevap(const int& stateflg, const int& dcmnt,double pet[CYCLE], double tair[CYCLE],double vap[CYCLE]);
     void setELMNTevap(const int& stateflg, const int& dcmnt,double pet[NVT][CYCLE], double tair[NVT][CYCLE],double vap[NVT][CYCLE],const int& vt);

//     void setELMNTflux(void);
     void setELMNTflux(const int& vt);

     void initrun(ofstream& rflog1, const int& equil);
//     void massbal(double y[NUMEQ], double prevy[NUMEQ]);
      void massbal(double y[NVT][NUMEQ], double prevy[NVT][NUMEQ], const int& vt);

//      void monthxclm(const int& dcmnt, const double& tgppopt, const int& dm);
      void monthxclm(const int& dcmnt, double& tgppopt, const int& dm, const int& vt);

 // void resetODEflux(double y[]);
      void resetODEflux(double y[NVT][NUMEQ], const int& vt);

//     void resetYrFlux(void);
     void resetYrFlux(const int& vt);

//     void setMonth(int& dm, double y[]);
     void setMonth(const int& dm, double y[NVT][NUMEQ], const int& vt);
//     void setPrevState(double prevState[], double currentState[]);
     void setPrevState(double prevState[NVT][NUMEQ], double currentState[NVT][NUMEQ],const int& vt);

// added RTIME to stepyr for soil thermal model
 //    int stepyr(const int& dyr, const int& itype, int& intflag, double& tol);
     int stepyr(const int& dyr, const int& itype, int& intflag, double& tol, const int& vt);

//     int transtepyr(const int& dyr, const int& itype, int& intflag, double& tol, const int& RTIME);
     int transtepyr(const int& dyr, const int& itype, int& intflag, double& tol, const int& RTIME, const int& vt);

//     int transient(const int& dyr, const int& itype, double& tol, const int& RTIME);
     int transient(const int& dyr, const int& itype, double& tol, const int& RTIME, const int& vt,int& init_carbonflag);


/* **************************************************************
			 Public Variables
************************************************************** */
   
     MPIatms atms;
     Tveg42 veg;
     Tsoil4 soil;
     Tmicrobe4 microbe;
     Humanact ag;
     Hydrology hyd; // added for hydrology model by QZ
     HYDM hydm; // added for soil 2 or 3 boxes soil water balance model
     Soilthermal sthermal;   // added for soil thermal model 
     double psstate[NVT][NUMEQ];
     double pddstate[NVT][NUMEQ];
     double elev;              // elevation (m)
     double orgcarbon;         //wsr: initial organic carbon
     double orgcarbon_initave;
     double init_depth;     
//     double nep[CYCLE];      // Net Ecosystem Production (g C / (sq. meter * month))
     double nep[NVT][CYCLE];      // Net Ecosystem Production (g C / (sq. meter * month))
//     double yrnep;           // Annual NEP (g C / (sq. meter * year))
     double yrnep[NVT];           // Annual NEP (g C / (sq. meter * year))
//     double cflux[CYCLE];
     double cflux[NVT][CYCLE];
//     double yrcflux;
     double yrcflux[NVT];
//     double totalc[CYCLE];   // total carbon storage (veg.plant + soil.org) (g C / sq. meter)
     double totalc[NVT][CYCLE];   // total carbon storage (veg.plant + soil.org) (g C / sq. meter)

     char predstr[NUMEQ][9]; // predstr[MAXPRED][9] changed to predstr[NUMEQ][9] by DWK on 20000202

     int ez;       // index for vegetation type (ez = veg.temveg-1)
     static int avlnflag;
     static int nfeed;
     static int rheqflag;
     static int moistlim;
     static int initbase;
     static int baseline;
     static int intflag;
     int nattempt;
     static int maxnrun;
     static int equil;
     static int runsize;
     static int maxyears;
     static int strteq;
     static int endeq;
     static int startyr;
     static int endyr;
     static int diffyr;
     static int wrtyr;
     long totyr;
     double tol;
     static double ctol;
     static double ntol;
     static double wtol;

     double nfert;
     char ecd1[80];
     char ecd2[80];
     int qualcon[MAXRTIME][NUMMSAC];

//     double y[NUMEQ];
     double y[NVT][NUMEQ];
//     double prevy[NUMEQ];
     double prevy[NVT][NUMEQ];

  // added for soil thaw-frozen by Q. Zhuang
     double k_coef;
     double n_frozen;
     double n_thaw;
	 double deltafreeze1;
	 double deltafreeze2;
	 double deltafreeze3;
     // string changed to char[80] by DWK on 20000210
//     char sitename[80];;
	 double soiltemptrans[NVT][CYCLE];
	 double soiltemp[NVT][CYCLE];
/* **************************************************************
			Public Parameters
************************************************************** */

     double vegca[MAXCMNT];
     double vegcb[MAXCMNT];

     double strna[MAXCMNT];
     double strnb[MAXCMNT];

     double solca[MAXCMNT];
     double solcb[MAXCMNT];

     double solna[MAXCMNT];
     double solnb[MAXCMNT];

     double avlna[MAXCMNT];
     double avlnb[MAXCMNT];

     double stona[MAXCMNT];
     double stonb[MAXCMNT];
     
     double vegca_p[MAXCMNT];
     double vegcb_p[MAXCMNT];

     double strna_p[MAXCMNT];
     double strnb_p[MAXCMNT];

     double solca_p[MAXCMNT];
     double solcb_p[MAXCMNT];

     double solna_p[MAXCMNT];
     double solnb_p[MAXCMNT];

     double avlna_p[MAXCMNT];
     double avlnb_p[MAXCMNT];

     double stona_p[MAXCMNT];
     double stonb_p[MAXCMNT];


     double vegca_temp[MAXCMNT];
     double vegcb_temp[MAXCMNT];

     double strna_temp[MAXCMNT];
     double strnb_temp[MAXCMNT];

     double solca_temp[MAXCMNT];
     double solcb_temp[MAXCMNT];

     double solna_temp[MAXCMNT];
     double solnb_temp[MAXCMNT];

     double avlna_temp[MAXCMNT];
     double avlnb_temp[MAXCMNT];

     double stona_temp[MAXCMNT];
     double stonb_temp[MAXCMNT];
     
     

 /* **************************************************************
		     Private Variables
************************************************************** */

  private:

     int initFlag; // = 1 when TEM has been initialized for a grid cell
     int transferflag;  // = 1 when using the peatland community data
    ifstream fecd[MAXCMNT];
//     double gv;
    double gv[NVT];
//     double ksoil;
    double ksoil[NVT];
//     double rhmoist;
    double rhmoist[NVT];
//     double temp;
    double temp[NVT];
//     double respq10;
    double respq10[NVT];
//     double dq10;
     double dq10[NVT];
};

// Initialization of static members

int TTEM::avlnflag = 0;
int TTEM::nfeed = 0;
int TTEM::rheqflag = 0;
int TTEM::moistlim = 0;
int TTEM::initbase = 0;
int TTEM::baseline = 0;
int TTEM::intflag = 0;

int TTEM::maxnrun = 0;
int TTEM::equil = 0;
int TTEM::runsize = 0;
int TTEM::maxyears = 0;
int TTEM::strteq = 0;
int TTEM::endeq = 0;
int TTEM::startyr = 0;
int TTEM::endyr = 0;
int TTEM::diffyr = 0;
int TTEM::wrtyr = 0;

double TTEM::ctol = 1.0;
double TTEM::ntol = 0.02;
double TTEM::wtol = 0.01;

#endif
