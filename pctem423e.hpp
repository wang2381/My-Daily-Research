/* *************************************************************
****************************************************************
PCTEM423.HPP - Class adds DOS-specific interface to the core
                 Terrestrial Ecosystem Model Version 4.2 to allow
     	         calibration of the model on a PC

Modifications:

2000103 - DWK added compiler directives
20000717 - QZ added file name for penman-monteith theory
***************************************************************
************************************************************** */

// Global constants

const int WSY = 5;
const int ESY = 3;

#if !defined(TTEM423_H)
  #include "ttem423e.cpp"      // PCTEM inherits TTEM class
#endif

#if !defined(PCTEM423_H)
#define PCTEM423_H

class PCTTEM : public TTEM {

  public:

    PCTTEM();

    enum seykey { NOEKEY,     GET_GPP,   GET_GPR,    GET_LTRC, GET_RH,
                  GET_LTRN,   GET_NMIN,  GET_NLST,   GET_NINP, GET_NEP,
                  GET_INGPP,  GET_INNPP, GET_INNUP,  GET_VEGC, GET_STRN,
                  GET_AGNPPC, GET_NMBL,  GET_NRSRB,  GET_NPP,  GET_VNUP,
                  GET_VSUP,   GET_VLUP,  GET_POTNPP, GET_L2SN, GET_FPC,
                  GET_LAI,    GET_LEAF,  GET_AGNPPN, GET_AGFRTN, GET_TSOIL,
                  GET_FRONTD }; // addition for Soil thermal model

    enum swykey { NOWKEY,     GET_RAIN,  GET_RPERC, GET_RRUN, GET_SNWFALL,
                  GET_SNWINF, GET_SPERC, GET_SRUN,  GET_PET,  GET_EET,
                  GET_SH2O,   GET_PCTP,  GET_VSM,   GET_WYLD,
                  GET_PCTP1,  GET_PCTP2, GET_PCTP3,
 	               GET_TOP20ST,GET_DST5,  GET_DST10, GET_DST20,
                  GET_DST50,  GET_DST100,GET_DST200,GET_FFRONTD,
                  GET_TTHAWB, GET_TTHAWE };

/* **************************************************************
			Public Functions
************************************************************** */

// modified by Q. Z for STM and HYD

//    inline seykey& next(seykey& s) { return s = (GET_AGFRTN== s) ? GET_GPP : seykey(s+1); }
//    inline seykey& prev(seykey& s) { return s = (GET_GPP == s) ? GET_AGFRTN : seykey(s-1); }

    inline seykey& next(seykey& s) { return s = (GET_FRONTD== s) ? GET_GPP : seykey(s+1); }
    inline seykey& prev(seykey& s) { return s = (GET_GPP == s) ? GET_FRONTD : seykey(s-1); }

//    inline swykey& next(swykey& s) { return s = (GET_WYLD== s) ? GET_RAIN : swykey(s+1); }
//    inline swykey& prev(swykey& s) { return s = (GET_RAIN == s) ? GET_WYLD : swykey(s-1); }

    inline swykey& next(swykey& s) { return s = (GET_TTHAWE== s) ? GET_RAIN : swykey(s+1); }
    inline swykey& prev(swykey& s) { return s = (GET_RAIN == s) ? GET_TTHAWE : swykey(s-1); }

    void displayOptionalEflx(seykey& s);
    void displayOptionalWflx(swykey& s);
    void displayState(const int& dyr, const int& dm, char* Climit, double pstate[NVT][NUMEQ], int& vt);

    double getOptionalEflx(const int& dm, const int& optflx, double pstate[NVT][NUMEQ],int& vt);
    double getOptionalWflx(const int& optflx, double pstate[NVT][NUMEQ],int& vt);

    int pcadapt(const int& numeq, double pstate[NVT][NUMEQ], double& ptol, const int& pdm,int& vt);

    void pcdisplayClm(char* calclm, double clouds[NVT][CYCLE], int& vt);
    void pcdisplayInitState(char* ecd, const double& col, const double& row,int& ez);
    void pcdisplayLimits(int& vt);
    void pcdisplayMonth(const int& dyr, const int& dm, double pstate[NVT][NUMEQ],int& vt);
    void pcdisplayOtherECD(int& ez);
    void pcdisplayPAR(char* calclm, double girr[NVT][CYCLE], int& vt);
    void pcdisplayStext(int& vt);
    void pcdisplayVegECD(int& ez);
    void pcdisplayYrAgTEM(const int& manflag, int& vt);
    void pcdisplayYrTEM(const int& manflag, int& vt);
    void pcdisplayYrWBM(int& vt);

    void pcinitRun(ofstream& rflog1);
//    int pcstepyr(const int& dyr, int& intflag, double& tol);
    int pcstepyr(const int& dyr, int& intflag, double& tol, int& vt);


/* **************************************************************
			 Public Variables
************************************************************** */


    seykey sey[ESY];
    swykey swy[WSY];

    int endgrid;
    int adapttol;
    int intbomb;
    int tolbomb;
    int mintflag;
    int topwind;
    int calwind;
    int firstcal;

// Names of files containing parameter values determined from the
//   literature

    char soilfile[40];
    char rootfile[40];
    char vegfile[40];
    char leaffile[40];
    char mcrvfile[40];
    char agfile[40];
// added for hydrology model by QZ
    char penmanfile[40];

  //added for soil thermal model
     char snowfile[40];
     char soillfile[40];
     char soiltfile[40];

    ofstream flog2;

  private:

/* **************************************************************
			Private Functions
************************************************************** */

    void pcdisplayDT(const double& tottime, const double& deltat);
    void pcdisplayODEerr(const int& test, double pstate[NVT][NUMEQ], int& vt);

    double outflux1;
    double outflux2;
    double outflux3;
    double outflux4;
    double outflux5;
    double outflux6;
    double outflux7;
    double outflux8;

    double eqsolc[NVT];
};

#endif
