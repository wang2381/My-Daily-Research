/* **************************************************************
*****************************************************************
HUMNACT423.HPP - describes human disturbances to natural ecosystems

                20000102 - DWK added compiler directives
*****************************************************************
************************************************************** */

#if !defined(BIOMS423_H)
  #include "bioms423.hpp"  // Humnact uses Biomass class
#endif

// Humnact also uses the global constants CYCLE, NUMMSAC, NUMVEG
//   and MAXRTIME

#if !defined(HUMNACT423_H)
#define HUMNACT423_H

class Humanact {

  public:

  Humanact();

/* **************************************************************
		 Public Functions
************************************************************** */

     void conversion(int& ez, Tveg42& veg, Tsoil4& soil, const int& vt);
     void getecd (ofstream& rflog1);
     void getecd (char ecd[80]);
     void resetPROD(const int& vt);
     void updateyr(const int& dyr, const int& vt);



  /* **************************************************************
		 Public Variables
************************************************************** */

    Biomass PROD1[NVT];

    Biomass initPROD10[NVT][10];
    Biomass PROD10[NVT];

    Biomass initPROD100[NVT][100];
    Biomass PROD100[NVT];

    Biomass TOTPROD[NVT];
    double totgC[NVT][CYCLE];

    Biomass convrtflx[NVT][CYCLE];
    double yrconvrtC[NVT];
    double yrconvrtN[NVT];

    Biomass vconvrtflx[NVT][CYCLE];
    double yrvconvrtC[NVT];
    double yrvconvrtN[NVT];

    Biomass sconvrtflx[NVT][CYCLE];
    double yrsconvrtC[NVT];
    double yrsconvrtN[NVT];



    double totnpp[NVT][CYCLE];
    Biomass formPROD1[NVT];
    Biomass formPROD10[NVT];
    Biomass formPROD100[NVT];
    Biomass formTOTPROD[NVT];

    Biomass PROD1decay[NVT];
    Biomass PROD10decay[NVT];
    Biomass PROD100decay[NVT];
    Biomass TOTPRODdecay[NVT];

    Biomass slash[NVT][CYCLE];
    double yrslashC[NVT];
    double yrslashN[NVT];
    double nretent[NVT][CYCLE];
    double nvretent[NVT][CYCLE];
    double nsretent[NVT][CYCLE];
    double yrnrent[NVT];
    double yrnvrent[NVT];
    double yrnsrent[NVT];
    double potnpp[NVT][CYCLE];
    double natsoil;
    double kd[NVT];
    double tpotnpp[NUMMSAC][MAXRTIME][CYCLE];
    Biomass npp[NVT][CYCLE];
    double yrnppC[NVT];
    double yrnppN[NVT];
    Biomass ltrfal[NVT][CYCLE];
    double yrltrc[NVT];
    double yrltrn[NVT];
    double fertn[NVT][CYCLE];
    double yrfertn[NVT];

    double slashpar[NUMVEG];
    double vconvert[NUMVEG];
    double prod10par[NUMVEG];
    double prod100par[NUMVEG];
    double sconvert[NUMVEG];
    double nvretconv[NUMVEG];
    double nsretconv[NUMVEG];
    double vrespar[NUMVEG];
    double cfall[NUMVEG];
    double nfall[NUMVEG];

    double RAP;
    double tRAP[MAXRTIME];
    double c2n[NVT];

    int state;
    int tstate[MAXRTIME];
    int distflag;
    int tlulcflag;
    int lulcyear[MAXRTIME];

    // Is RAP always going to be zero?
    // (if so, we don't need potential NPP data)
    int RAP0flag;

    int prvstate;
//    char predstr[35][9];
};

#endif

