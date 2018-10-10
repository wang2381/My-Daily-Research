/* **************************************************************
*****************************************************************
CTEM423.CPP - Calibration version of the Terrestrial Ecosystem
			Model Version 4.2

Modifications:

20000717 - QZ changhed for hydrology model
20001012 - Q.Z. bugs fixed besides adding STM and HYD module
20001022 -- Q.Z. to create the version including STM, HYD, and Moss component, and fire disturbance to those modules
*****************************************************************
************************************************************** */

#include<stdio.h>
#include<conio.h>
#include<iostream.h>
#include<fstream.h>
#include <fcntl.h>
#include<iomanip.h>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>
#include<cstring.h>

const int NVT = 2; // number of vegetation types for a grid cell, for instance, here we defile 0 is for tree and 1 for moss, Q. Z. added for moss component
const int CYCLE = 12;
const int MAXRTIME = 600;

//const int MAXPRED = 85;
const int MAXPRED = 85 + 29 +10;
const int MAXCAL = 21;

const double MISSING = -99999.9;

//  added for soil temperature
long int kswitch, fswitch;
int stmflg; // flag for soil thermal model
int subsist[NVT]; // to store plant function types

double ttsoil[12];
double ddst5[12];
double ddst10[12];
double ddst20[12];
double ddst50[12];
double ddst100[12];
double ddst200[12];
double ffrontd[12];
double tthawb[12];
double tthawe[12];

//  end of adding

int treeveg = 0; // 0 for default as tree vegetation

//Modules representing climate and TEM

#include "potclm423e.cpp"      // Potsclm Class
#include "pctem423e.cpp"       // PCTEM Class

Potsclm clm;
PCTTEM tem;


enum ecalkey { C_CFALL, C_C2NB, C_CNEVN, C_CNMN,  C_CMAX,
               C_KRB,   C_KDC,  C_NMAX,  C_NFALL, C_MNUP,
               C_FPCMX, C_SLA,  C_CO2,   C_TOL,   C_LFMX,
               C_KLFC,  C_COV,  C_CNSOL, C_MNFIX, C_NLOS,
               C_MOPT,

               C_STXT,   C_AVLNFLG, C_NFEED, C_BLINE, C_MOISTL,
               C_CLM,    C_SEY0,    C_SEY1,  C_SEY2,  C_RESET,
               C_MANFLG, C_AGSTATE, C_STSTATE, C_TREEVEG, C_CALWIND };

/* *************************************************************
		 Function Declarations
************************************************************* */

inline ecalkey& next(ecalkey& s)
{
  return s = (C_CALWIND == s) ? C_CFALL : ecalkey(s+1);
}

inline ecalkey& prev(ecalkey& s)
{
  return s = (C_CFALL == s) ? C_CALWIND : ecalkey(s-1);
}

void calecdin(void);
void calclmin(int& vt);
void displayECalibPar(ecalkey& dcal, int& vt);
void ecalibrate(int k, int& vt);
void setCalVar(int& vt);
//void showlphn(void);
void showclm(int& vt);
void updateECalibPar(ecalkey& dcal, const double& factor, int& t, int& vt);
void wcalibrate(int k);

// *************************************************************

int reset;
int equilsol;

int dyr;
int newclm = 0;
int newtext = 0;
double ecalvar[21];
ecalkey evar;
int ivar;

int vt; // added for moss component

int kdinflg = 0;
int stateflag = 0;
int manflag = 0;

double col;
double row;
double lon;
double lat;
double elev;

char calclm[80];
ofstream flog1;


/* *************************************************************
**********************START MAIN PROGRAM************************
************************************************************* */

int main()
{

  int i;
  int k;
  int key;

  kswitch =0; // added for soil thermal model
  stmflg = 1; // added for soil thermal model
  fswitch = 0; //

  reset = 1;
  clrscr();


/* Open log file
*  flog1.open("stem.log"); */

  cout << endl << endl << "Initialize conditions for TEM:" << endl << endl;


  tem.pcinitRun(flog1);

  calclmin(vt); // has been modified for moss component by Q.Z.

  calecdin(); // has been modified for moss component by Q.Z.

  for (i = 0; i < CYCLE; i++)
  {
    tem.atms.co2[i] = tem.atms.co2level;
  }

  tem.atms.initco2 = tem.atms.co2level;


 for (vt=0; vt < NVT; vt++) { // added for dealing with moss
   tem.veg.cmnt= subsist[vt];

   showclm(vt);

  //  tem.soil.xtext(tem.veg.cmnt, tem.soil.pctsilt, tem.soil.pctclay);
  // added for 3-box hydrology
  tem.soil.hydm_xtext(tem.veg.cmnt,tem.soil.pctsilt,tem.soil.pctclay, vt); // modified for two boxes hydrology model

//  tem.ECDsetELMNTstate(tem.veg.cmnt, tem.soil.psiplusc);
  tem.ECDsetELMNTstate(tem.veg.cmnt, tem.soil.psiplusc,vt);
//  tem.setELMNTflux();
  tem.setELMNTflux(vt);
//  tem.setPrevState(tem.prevy, tem.y);
  tem.setPrevState(tem.prevy, tem.y, vt);

  cout << "Begin to initialize ecd" << endl;
//  tem.setELMNTecd(kdinflg,tem.veg.cmnt, tem.soil.psiplusc);
  tem.setELMNTecd(kdinflg,tem.veg.cmnt, tem.soil.psiplusc,vt);

  // modified for hydrology model by QZ
//  tem.setELMNTevap(stateflag, tem.veg.cmnt, tem.atms.pet, tem.atms.tair,tem.atms.vap);
  tem.setELMNTevap(stateflag, tem.veg.cmnt, tem.atms.pet, tem.atms.tair,tem.atms.vap,vt);
//  tem.pcdisplayStext();
  tem.pcdisplayStext(vt);
} // end of vt

  tem.mintflag = 0;
  tem.tol = tem.inittol;
  tem.retry = 0;
  tem.endgrid = 0;

  tem.firstcal = 1;
  tem.topwind = 0;

  clrscr();
  cout << "Do you want to start with the WBM screen? ";
  if (toupper(getch()) == 'Y')
  {
    tem.calwind = 1;
    cout << "Y" << endl;
  }
  else
  {
    tem.calwind = 2;
    cout << "N" << endl;
  }

  dyr = 0;
  ivar = 1;
  tem.endeq = 1;

  tem.retry = 0;
  tem.endgrid = 0;
  equilsol = 0;

  tem.ag.state = 0;
  tem.ag.prvstate = 0;

  while (dyr > -10)
  {
   //    tem.pcstepyr(dyr, tem.mintflag, tem.tol);
  for (vt=0; vt < NVT; vt++) {

    tem.pcstepyr(dyr, tem.mintflag, tem.tol, vt);

   // Update annual agricultural product pools and fluxes
  //  tem.ag.updateyr(dyr);
     tem.ag.updateyr(dyr,vt);
  }
    ++dyr;

// Check to see if steady state conditions have been reached.

   for (vt=0; vt < NVT; vt++) {

    if (equilsol == 1) { tem.endgrid = 1; }

    if (tem.equil == 1 && dyr >= tem.strteq && equilsol == 0)
    {
      if (tem.wtol >= fabs(tem.atms.yrrain[vt] + tem.atms.yrsnowfall[vt] \
                      - tem.atms.yreet[vt] - tem.soil.yrrrun[vt] - tem.soil.yrsrun[vt]) \
                      && tem.nfeed == 0 && tem.rheqflag == 0 \
                      && (tem.ctol >= fabs(tem.veg.yrnpp[vt] - tem.veg.yrltrc[vt])) \
                      && (tem.ctol >= fabs(tem.veg.yrnpp[vt] - tem.veg.yrltrc[vt])))
      {
 	equilsol = 1;
      }
      if (tem.wtol >= fabs(tem.atms.yrrain[vt] + tem.atms.yrsnowfall[vt] \
                      - tem.atms.yreet[vt] - tem.soil.yrrrun[vt] - tem.soil.yrsrun[vt]) \
                      && tem.nfeed == 0 && tem.rheqflag == 1 \
                      && (tem.ctol >= fabs(tem.yrnep[vt]))\
                      && (tem.ctol >= fabs(tem.veg.yrnpp[vt] - tem.veg.yrltrc[vt])) \
                      && (tem.ctol >= fabs(tem.veg.yrltrc[vt] - tem.microbe.yrrh[vt])))
      {
	equilsol = 1;
      }
      if (tem.wtol >= fabs(tem.atms.yrrain[vt] + tem.atms.yrsnowfall[vt] \
                      - tem.atms.yreet[vt] - tem.soil.yrrrun[vt] - tem.soil.yrsrun[vt]) \
                      && tem.nfeed == 1 && tem.rheqflag == 1 \
                      && (tem.ntol >= fabs(tem.soil.yrnin[vt] - tem.soil.yrnlost[vt])) \
                      && (tem.ntol >= fabs(tem.veg.yrnup[vt] - tem.veg.yrltrn[vt])) \
                      && (tem.ntol >= fabs(tem.veg.yrnup[vt] - tem.microbe.yrnmin[vt])) \
                      && (tem.ntol >= fabs(tem.veg.yrltrn[vt] - tem.microbe.yrnmin[vt])) \
                      && (tem.ctol >= fabs(tem.yrnep[vt])) \
                      && (tem.ctol >= fabs(tem.veg.yrnpp[vt] - tem.veg.yrltrc[vt])) \
                      && (tem.ctol >= fabs(tem.veg.yrltrc[vt] - tem.microbe.yrrh[vt])))
      {
	equilsol = 1;
      }
    }
   } // end for vt


 if ( treeveg == 0)    // for choosing calibration for tree or moss
    vt = 0;
  else vt =1;

    if (tem.calwind == 1 )
    {
      tem.pcdisplayYrWBM(vt);

      window(46,24,80,24);
      gotoxy(1,1);
      printf(" 4.2.3-15JUN00  YRNETW = %8.6lf  ",
            (tem.atms.yrrain[vt]+tem.atms.yrsnowfall[vt]-tem.atms.yreet[vt]
            -tem.soil.yrrrun[vt]-tem.soil.yrsrun[vt]));
    }
    else
    {
      if (tem.ag.state == 0) { tem.pcdisplayYrTEM(manflag,vt); }
      else { tem.pcdisplayYrAgTEM(manflag,vt); }

      tem.pcdisplayLimits(vt);

      window(46,24,80,24);
      gotoxy(1,1);
      printf("4.2.3-15JUN00 YRNEP = %8.3lf     ",tem.yrnep[vt]);
    }

    if (tem.adapttol == 1 && (dyr+1) == tem.maxyears)
    {
      tem.endgrid = 1;
      tem.tolbomb = 1;
    }

    if (manflag == 1)
    {
      k = 79;
      if (tem.calwind == 1) { wcalibrate(k); }
      else { ecalibrate(k, vt); }
    }

    if (kbhit() != 0)
    {
      key = getch();
      if (key == 0)
      {
	k = getch();
	if (k == 79)
        {
	  if (tem.calwind == 1) {wcalibrate(k);}
	  else { ecalibrate(k, vt); }
	}
      }
    }

    if (tem.endgrid == 1)
    {
      k = 79;
      if (tem.calwind == 1) { wcalibrate(k); }
      else { ecalibrate(k, vt); }
    }
  }

  return 0;

};

/* *************************************************************
*********************END OF MAIN PROGRAM************************
************************************************************* */


/* *************************************************************
************************************************************* */


// modified for hydrology model, adding vapor pressure for calibration by QZ
void calclmin(int& vt)
{

  int i;
  int j;
  ifstream fclm;
  char dummy[12];
  int dumveg;

  tem.soil.getecd(tem.soilfile);

  cout << endl << "Enter file with climatic data: ";
  cin >> calclm;
  fclm.open(calclm,ios::in);

  if (!fclm)
  {
    cerr << "\nCannot open " << calclm << " for data input" << endl;
    exit(-1);
  }

// Read in basic climate data

  fclm >> dummy >> col;
  fclm >> dummy >> row;
  fclm >> dummy >> tem.elev;
  fclm >> dummy >> dumveg;
  fclm >> dummy >> tem.soil.pctsand;
  fclm >> dummy >> tem.soil.pctsilt;
  fclm >> dummy >> tem.soil.pctclay;
  fclm >> dummy >> tem.soil.wsoil;
  fclm >> dummy >> tem.ag.RAP;

  fclm >> dummy >> dummy >> dummy >> dummy >> dummy; // adding dummy for hydrology model

  for (i = 0; i < CYCLE; i++)
  {
//    fclm >> j >> clm.clds[i] >> tem.atms.tair[i] >> tem.atms.prec[i];
    fclm >> j >> clm.clds[vt][i] >> tem.atms.tair[vt][i] >> tem.atms.prec[vt][i] >> tem.atms.vap[vt][i];
  }

 for (i = 0; i < CYCLE; i++) // assign value to vege 2
   {
  clm.clds[1][i] = clm.clds[0][i];
  tem.atms.tair[1][i] = tem.atms.tair[0][i];
  tem.atms.prec[1][i] =  tem.atms.prec[0][i];
  tem.atms.vap[1][i] =   tem.atms.vap[0][i];
  }

  fclm.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void calecdin(void)
{

  char ecd[80];

  tem.soil.getrootz(tem.rootfile);

  tem.veg.getecd(tem.vegfile);
  tem.veg.leafyrs = 10;
  tem.veg.getleafecd(tem.leaffile);

  tem.microbe.getvegecd(tem.mcrvfile);
  tem.ag.getecd(tem.agfile);
// added for hydrology model by QZ
  tem.hyd.getpenmonecd(tem.penmanfile);

    // added for soil thermal model
  tem.sthermal.getsnowecd(tem.snowfile);
  tem.sthermal.getsoillecd(tem.soillfile);
  tem.sthermal.getsoiltecd(tem.soiltfile);

  //  tem.hyd.getecd(tem.penmanfile);

  cout << endl << "Enter name of file with parameter values: ";
  cin >> ecd;
  tem.getsitecd(0, ecd);

  for (vt = 0; vt < NVT; vt++)
  {
   tem.veg.cmnt= subsist[vt];

   tem.veg.adjc2n[vt] = 1.0 + (tem.veg.dc2n[vt] * (tem.atms.co2[11] - tem.atms.initco2));
   tem.veg.cneven[tem.veg.cmnt] = tem.veg.initcneven[tem.veg.cmnt]
                                 * tem.veg.adjc2n[vt];
  tem.soil.hydm_xtext(tem.veg.cmnt,tem.soil.pctsilt,tem.soil.pctclay, vt); // modified for two boxes hydrology model
  }

//  tem.soil.xtext(tem.veg.cmnt, tem.soil.pctsilt, tem.soil.pctclay);
   // added for 3-box hydrology

  for (vt = 0; vt < NVT; vt++)
  {
   tem.veg.cmnt= subsist[vt];
   tem.pcdisplayInitState(ecd,col,row,tem.veg.cmnt );
   tem.pcdisplayVegECD(tem.veg.cmnt);
   tem.pcdisplayOtherECD(tem.veg.cmnt);
  }

};

/***************************************************************
***************************************************************/


/* *************************************************************
************************************************************** */

void displayECalibPar(ecalkey& dcal, int& vt)
{

   if (vt ==0)
   tem.veg.cmnt= subsist[0];
   else    tem.veg.cmnt= subsist[1];

  switch(dcal)
  {
    case C_CFALL:   printf("CFALL = %8.6lf   ",tem.veg.cfall[tem.veg.cmnt]);
                    break;
    case C_C2NB:    printf("PLNTCNB = %8.2lf ",tem.veg.c2nb[tem.veg.cmnt]);
                    break;
    case C_CNEVN:   printf("INCNEV = %8.2lf  ",tem.veg.initcneven[tem.veg.cmnt]);
                    break;
    case C_CNMN:    printf("CNMIN = %8.2lf   ",tem.veg.cnmin[tem.veg.cmnt]);
                    break;
    case C_CMAX:    printf("CMAX = %8.2lf    ",tem.veg.cmax[vt]);
                    break;
    case C_KRB:     printf("KRB = %8.6lf     ",tem.veg.krb[tem.veg.cmnt]);
                    break;
    case C_KDC:     printf("KDC = %8.6lf     ",tem.microbe.kdc[vt]);
                    break;
    case C_NMAX:    printf("NMAX = %8.4lf    ",tem.veg.nmax[vt]);
                    break;
    case C_NFALL:   printf("NFALL = %8.6lf   ",tem.veg.nfall[tem.veg.cmnt]);
                    break;
    case C_MNUP:    printf("MNUP = %8.4lf    ",tem.microbe.nup[vt]);
                    break;
    case C_FPCMX:   printf("FPCMX = %8.6lf   ",tem.veg.fpcmax[tem.veg.cmnt]);
                    break;
    case C_SLA:     printf("SLA = %8.4lf     ",tem.veg.sla[tem.veg.cmnt]);
                    break;
    case C_CO2:     printf("CO2 = %8.1lf     ",tem.atms.co2level);
                    break;
    case C_TOL:     printf("INT TOL = %8.6lf ",tem.tol);
                    break;
    case C_LFMX:    printf("LEAFMXC = %8.2lf ",tem.veg.leafmxc[tem.veg.cmnt]);
                    break;
    case C_KLFC:    printf("KLEAFC = %8.2lf   ",tem.veg.kleafc[tem.veg.cmnt]);
                    break;
    case C_COV:     printf("COV = %8.6lf     ",tem.veg.cov[tem.veg.cmnt]);
                    break;
    case C_CNSOL:   printf("CNSOIL = %8.2lf  ",tem.microbe.cnsoil[tem.veg.cmnt]);
                    break;
    case C_MNFIX:   printf("NFIX   = %8.1lf  ",tem.microbe.nfixpar[tem.veg.cmnt]);
                    break;
    case C_NLOS:    printf("NLOSS  = %8.1lf  ",tem.soil.nloss[tem.veg.cmnt]);
                    break;
    case C_MOPT:    printf("MOPT = %8.4lf    ",tem.microbe.moistopt[tem.veg.cmnt]);
                    break;
    case C_STXT:    printf("NEW TEXTURE? (%4.2lf)",tem.soil.psiplusc);
                    break;
    case C_AVLNFLG: if (tem.avlnflag == 1) { printf("AVAILN = ON"); }
	            if (tem.avlnflag == 0) { printf("AVAILN = OFF"); }
	            break;
    case C_NFEED:   if (tem.nfeed == 1) { printf("NFEED = ON");  }
	            if (tem.nfeed == 0) { printf("NFEED = OFF"); }
	            break;
    case C_BLINE:   if (tem.baseline == 1) { printf ("BLINE = ON"); }
	            if (tem.baseline == 0) { printf ("BLINE = OFF"); }
	            break;
    case C_MOISTL:  if (tem.moistlim == 1) { printf("MOISTLIM = ON");  }
	            if (tem.moistlim == 0) { printf("MOISTLIM = OFF"); }
	            break;
    case C_CLM:     printf("NEW CLIMATE?");
	            break;
    case C_SEY0:    tem.displayOptionalEflx(tem.sey[0]);
	            break;
    case C_SEY1:    tem.displayOptionalEflx(tem.sey[1]);
	            break;
    case C_SEY2:    tem.displayOptionalEflx(tem.sey[2]);
	            break;
    case C_RESET:   if (reset == 1) { cout << "reset = ON" << endl; }
	            else { cout << "reset = OFF" << endl; }
	            break;
    case C_MANFLG:  if (manflag == 1) { cout << "manual = ON" << endl; }
	            else { cout << "manual = OFF" << endl; }
	            break;
    case C_AGSTATE: if (tem.ag.state == 1) { cout << "agstate = ON" << endl; }
	            else { cout << "agstate = OFF" <<endl; }
	            break;
           //added for soil thermal model
    case C_STSTATE: if (stmflg == 1) cout << "ststate = ON" << endl;
             else cout << "ststate = OFF" <<endl;
           break;
    case C_TREEVEG: if (treeveg == 1) cout << "PFT2-cali" << endl;
             else cout << "PFT1-cali" <<endl;   // for choosing calibration for tree or moss
	     break;

  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void ecalibrate(int t, int& vt)
{

  double efactor;
  double invar;

  tem.endgrid = 0;

  window(43,23,80,23);
  printf("AN = %1d, NF = %1d, BL = %1d, MF = %1d, RS = %1d",tem.avlnflag,
        tem.nfeed,tem.baseline,tem.moistlim,reset);

  window(57,1,64,1);
  tem.displayOptionalEflx(tem.sey[0]);

  window(65,1,72,1);
  tem.displayOptionalEflx(tem.sey[1]);

  window(73,1,79,1);
  tem.displayOptionalEflx(tem.sey[2]);

  window(1,23,15,23);
  cout << "CALIBRATION:   " << endl;

  if (tem.firstcal == 1)
  {
    evar = C_CFALL;
    tem.firstcal = 0;
  }

  window(16,23,34,23);
  delline();
  displayECalibPar(evar, vt);

  while (t == 79 || t == 72 || t == 75 || t == 77 || t == 80 || t == 82)
  {

    if (tem.tolbomb == 0 && tem.intbomb == 0)
    {
      t = getch(); t = getch();
    }
    else
    {

      // If simulation reaches RUNSIZE automatically adjust the
      //   tolerance level (i.e. tol)
      if (tem.tolbomb == 1)
      {
        tem.retry = 1;
        tem.tolbomb = 0;
        tem.nattempt += 1;
        evar = C_TOL;
        if (tem.nattempt < tem.maxnrun)
        {
	  tem.tol /= 10.0;
	  t = 0;
        }
        else
        {
	  tem.tol = tem.inittol;
	  tem.nattempt = 0;
        }
      }

      // If simulation "black holes", allow user to adjust tolerance level
      if (tem.intbomb == 1)
      {
        tem.retry = 1;
        tem.intbomb = 0;
        evar = C_TOL;
      }
    }

    if (tem.intbomb == 1) { evar = C_TOL; }

    // Set calvar[] with current TEM parameter values
    setCalVar(vt);

    efactor = 1.0;
    window(35,23,42,23);
    delline();
    //    printf(" ");

    // Select or modify a specific calibration parameter
    //   represented by calvar[]

    switch (t)
    {
      case 79: break;
      case 72: efactor = 1.01; break;
      case 80: efactor = 0.99; break;
      case 77: next(evar); break;
      case 75: prev(evar); break;
      case 82: if (ecalkey(evar) < MAXCAL)
               {
	         cin >> invar;
	         ecalvar[evar] = invar;
	       }
	       break;
      default: t = 0; break;
    }

    window(35,23,42,23);
    delline();

    window(16,23,34,23);
    delline();

//    if (ivar < 1) { ivar = 34; }
//    if (ivar > 34) { ivar = 1; }

    // Update TEM calibration parameters with new calvar[] value

    updateECalibPar(evar, efactor, t, vt);

    if (t == 0) break;
  }


  if (reset == 1 || tem.retry == 1)
  {
    tem.retry = 0;
    dyr = 0;

   for (vt=0; vt < NVT; vt++) { // added for dealing with moss
    tem.veg.cmnt= subsist[vt];

    tem.ECDsetELMNTstate(tem.veg.cmnt, tem.soil.psiplusc,vt);
    tem.setELMNTflux(vt);
    tem.setPrevState(tem.prevy, tem.y,vt);
    tem.prevy[vt][tem.I_UNRMLF] = 0.5;
    tem.prevy[vt][tem.I_LEAF] = 0.5;
     } // end of vt
  }


  if (newtext == 1 || newclm == 1)
  {
    window(1,1,80,24);
    gotoxy(1,1);
    clrscr();
  for (vt=0; vt < NVT; vt++) { // added for dealing with moss
    tem.veg.cmnt= subsist[vt];
    tem.setELMNTecd(kdinflg,tem.veg.cmnt,tem.soil.psiplusc,vt);
    tem.pcdisplayStext(vt);
    tem.topwind = 0;
    newclm = 0;
    newtext = 0;
  } // end of vt
 }


  if (equilsol == 1)
  {
    equilsol = 0;
    dyr = 0;
  }

};

/* *************************************************************
************************************************************** */


/***************************************************************
***************************************************************/

void setCalVar(int& vt)
{

  if (vt ==0) tem.veg.cmnt = subsist[0];
  else tem.veg.cmnt = subsist[1];

  ecalvar[C_CFALL]  =  tem.veg.cfall[tem.veg.cmnt];
  ecalvar[C_C2NB]  =  tem.veg.c2nb[tem.veg.cmnt];
  ecalvar[C_CNEVN]  =  tem.veg.initcneven[tem.veg.cmnt];
  ecalvar[C_CNMN]  =  tem.veg.cnmin[tem.veg.cmnt];
  ecalvar[C_CMAX]  =  tem.veg.cmax[vt];
  ecalvar[C_KRB]  =  tem.veg.krb[tem.veg.cmnt];
  ecalvar[C_KDC]  =  tem.microbe.kdc[vt];
  ecalvar[C_NMAX]  =  tem.veg.nmax[vt];
  ecalvar[C_NFALL]  =  tem.veg.nfall[tem.veg.cmnt];
  ecalvar[C_MNUP]  =  tem.microbe.nup[vt];
  ecalvar[C_FPCMX] =  tem.veg.fpcmax[tem.veg.cmnt];
  ecalvar[C_SLA] =  tem.veg.sla[tem.veg.cmnt];
  ecalvar[C_CO2] =  tem.atms.co2level;
  ecalvar[C_TOL] =  tem.tol;
  ecalvar[C_LFMX] =  tem.veg.leafmxc[tem.veg.cmnt];
  ecalvar[C_KLFC] =  tem.veg.kleafc[tem.veg.cmnt];
  ecalvar[C_COV] =  tem.veg.cov[tem.veg.cmnt];
  ecalvar[C_CNSOL] =  tem.microbe.cnsoil[tem.veg.cmnt];
  ecalvar[C_MNFIX] =  tem.microbe.nfixpar[tem.veg.cmnt];
  ecalvar[C_NLOS] =  tem.soil.nloss[tem.veg.cmnt];
  ecalvar[C_MOPT] =  tem.microbe.moistopt[tem.veg.cmnt];

}

/***************************************************************
***************************************************************/


/***************************************************************
***************************************************************/

void showclm(int& vt)
{

  tem.pcdisplayClm(calclm, clm.clds,vt);

  // Determine monthly gross irradiance (GIRR)

  lat = (double) row;

  clm.yrsumday = 0.0;
  for (int dm = 0; dm < CYCLE; dm++)
  {
    clm.girr[vt][dm] = clm.xgirr(lat,dm,clm.yrsumday);

// Determine net irradiance (NIRR)

    tem.atms.nirr[vt][dm] = clm.xnirr(clm.clds[vt][dm],clm.girr[vt][dm]);

// Determine photosynthetically active radiation (PAR)

    tem.atms.par[vt][dm] = clm.xpar(clm.clds[vt][dm],tem.atms.nirr[vt][dm]);
  }

  tem.pcdisplayPAR(calclm, clm.girr,vt);

 //added for soil thermal model

 // for (vt = 0; vt < NVT; vt++)
//  {
//   tem.veg.cmnt= subsist[vt];
   tem.sthermal.showsnowecd(tem.veg.cmnt); // changed te.ez to tem.veg.cmnt
   tem.sthermal.showsoillecd(tem.veg.cmnt);
   tem.sthermal.showsoiltecd(tem.veg.cmnt);
 // }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void updateECalibPar(ecalkey& dcal, const double& factor, int& t, int& vt)
{
  int i;

   if (vt ==0)
   tem.veg.cmnt= subsist[0];
   else    tem.veg.cmnt= subsist[1];

  switch (dcal)
  {
    case C_CFALL :   tem.veg.cfall[tem.veg.cmnt] = ecalvar[C_CFALL] * factor;
                     printf("CFALL = %8.6lf   ",tem.veg.cfall[tem.veg.cmnt]);
	             break;
    case C_C2NB :    tem.veg.c2nb[tem.veg.cmnt] = ecalvar[C_C2NB] * factor;
       	             printf("PLNTCNB = %8.2lf ",tem.veg.c2nb[tem.veg.cmnt]);
	             break;
    case C_CNEVN :   tem.veg.initcneven[tem.veg.cmnt] = ecalvar[C_CNEVN] * factor;
	             printf("INCNEV = %8.2lf  ",tem.veg.initcneven[tem.veg.cmnt]);
	             break;
    case C_CNMN :    tem.veg.cnmin[tem.veg.cmnt] = ecalvar[C_CNMN] * factor;
	             printf("CNMIN = %8.2lf   ",tem.veg.cnmin[tem.veg.cmnt]);
	             break;
    case C_CMAX :    tem.veg.cmax[vt] = ecalvar[C_CMAX] * factor;
	             printf("CMAX = %8.2lf    ",tem.veg.cmax[vt]);
	             break;
    case C_KRB :     tem.veg.krb[tem.veg.cmnt] = ecalvar[C_KRB] * factor;
	             printf("KRB = %8.6lf     ",tem.veg.krb[tem.veg.cmnt]);
	             break;
    case C_KDC :     tem.microbe.kdc[vt] = ecalvar[C_KDC] * factor;
	             printf("KDC = %8.6lf     ",tem.microbe.kdc[vt]);
	             break;
    case C_NMAX :    tem.veg.nmax[vt] = ecalvar[C_NMAX] * factor;
	             printf("NMAX = %8.4lf    ",tem.veg.nmax[vt]);
	             break;
    case C_NFALL :   tem.veg.nfall[tem.veg.cmnt] = ecalvar[C_NFALL] * factor;
	             printf("NFALL = %8.6lf   ",tem.veg.nfall[tem.veg.cmnt]);
	             break;
    case C_MNUP :    tem.microbe.nup[vt] = ecalvar[C_MNUP] * factor;
	             printf("NUP = %8.4lf     ",tem.microbe.nup[vt]);
	             break;
    case C_FPCMX :   tem.veg.fpcmax[tem.veg.cmnt] = ecalvar[C_FPCMX] * factor;
	             printf("FPCMX = %8.6lf   ",tem.veg.fpcmax[tem.veg.cmnt]);
	             break;
    case C_SLA :     tem.veg.sla[tem.veg.cmnt] = ecalvar[C_SLA] * factor;
	             printf("SLA = %8.4lf    ",tem.veg.sla[tem.veg.cmnt]);
	             break;
    case C_CO2 :     tem.atms.co2level = ecalvar[C_CO2] * factor;
	             printf("CO2 = %8.1lf     ",tem.atms.co2level);
	             for (i = 0; i < CYCLE; i++)
                     {
	               tem.atms.co2[i] = tem.atms.co2level;
	             }
	             break;
    case C_TOL :     tem.tol = ecalvar[C_TOL] * factor;
	             printf("INT TOL = %8.6lf ",tem.tol);
	             break;
    case C_LFMX :    tem.veg.leafmxc[tem.veg.cmnt] = ecalvar[C_LFMX] * factor;
	             printf("LEAFMXC = %8.6lf  ",tem.veg.leafmxc[tem.veg.cmnt]);
	             break;
    case C_KLFC :    tem.veg.kleafc[tem.veg.cmnt] = ecalvar[C_KLFC] * factor;
	             printf("KLEAFC = %8.6lf   ",tem.veg.kleafc[tem.veg.cmnt]);
//             tem.microbe.decay = 0.26299 + (1.14757*tem.microbe.procmntos[tem.veg.cmnt])
//                                 - (0.42956*pow(tem.microbe.procmntos[tem.veg.cmnt],2.0));
	             break;
    case C_COV :     tem.veg.cov[tem.veg.cmnt] = ecalvar[C_COV] * factor;
	             printf("COV = %8.6lf    ",tem.veg.cov[tem.veg.cmnt]);
	             break;
    case C_CNSOL :   tem.microbe.cnsoil[tem.veg.cmnt] = ecalvar[C_CNSOL] * factor;
	             printf("CNSOIL = %8.2lf  ",tem.microbe.cnsoil[tem.veg.cmnt]);
	             break;
    case C_MNFIX :   tem.microbe.nfixpar[tem.veg.cmnt] = ecalvar[C_MNFIX] * factor;
	             printf("NFIX   = %8.1lf  ",tem.microbe.nfixpar[tem.veg.cmnt]);
	             break;
    case C_NLOS :    tem.soil.nloss[tem.veg.cmnt] = ecalvar[C_NLOS] * factor;
	             printf("NLOSS  = %8.1lf  ",tem.soil.nloss[tem.veg.cmnt]);
	             break;
    case C_MOPT :    tem.microbe.moistopt[tem.veg.cmnt] = ecalvar[C_MOPT] * factor;
	             printf("MOPT = %8.4lf    ", tem.microbe.moistopt[tem.veg.cmnt]);
	             break;
    case C_STXT :    if (factor != 1)
                     {
	               window(1,1,80,24);
	               gotoxy(1,1);
	               clrscr();
	               cout << "Enter proportion of silt plus clay: ";
	               cin >> tem.soil.psiplusc;
	               cout << "Enter percent clay: ";
	               cin >> tem.soil.pctclay;
	               tem.soil.pctsilt = (tem.soil.psiplusc * 100.0)
                                          - tem.soil.pctclay;
	               tem.soil.pctsand = (1.0 - tem.soil.psiplusc)* 100.0;
	               tem.soil.xtext(tem.veg.cmnt,
                                      tem.soil.pctsilt, tem.soil.pctclay,vt);
	               newtext = 1;
	               reset = 1;
	               t = 0;
	             }
	             else
                     {
	               printf("NEW TEXTURE? (%4.2lf)",tem.soil.psiplusc);
	             }
	             break;
    case C_AVLNFLG : if (factor != 1) { tem.avlnflag += 1; }
	             if (tem.avlnflag > 1) { tem.avlnflag = 0; }
	             if (tem.avlnflag == 1) { cout << "AVAILN = ON"; }
	             if (tem.avlnflag == 0) { cout << "AVAILN = OFF"; }
	             break;
    case C_NFEED :   if (factor != 1) { tem.nfeed += 1; }
	             if (tem.nfeed > 1) { tem.nfeed = 0; }
	             if (tem.nfeed == 1) { cout << "NFEED = ON"; }
	             if (tem.nfeed == 0) { cout << "NFEED = OFF"; }
	             break;
    case C_BLINE :   if (factor != 1) { tem.baseline += 1; }
	             if (tem.baseline > 1) { tem.baseline = 0; }
	             if (tem.baseline == 1) { cout << "BLINE = ON"; }
	             if (tem.baseline == 0) { cout << "BLINE = OFF"; }
	             break;
    case C_MOISTL :  if (factor != 1) { tem.moistlim += 1; }
	             if(tem.moistlim > 1) { tem.moistlim = 0; }
	             if (tem.moistlim == 1) { cout << "MOISTLIM = ON"; }
	             if (tem.moistlim == 0) { cout << "MOISTLIM = OFF"; }
	             break;
    case C_CLM :     if (factor != 1)
                     {
	               window(1,1,80,24);
	               gotoxy(1,1);
	               clrscr();
	               calclmin(vt);
	               showclm(vt);
	               tem.soil.xtext(tem.veg.cmnt,
                                      tem.soil.pctsilt, tem.soil.pctclay,vt);
	               newclm = 1;
	               reset = 1;
	               t = 0;
	             }
	             else
                     {
	               cout << "NEW CLIMATE?       ";
	             }
	             break;
    case C_SEY0 :    if (factor > 1) { tem.next(tem.sey[0]); }
	             if (factor < 1) { tem.prev(tem.sey[0]); }
	             tem.displayOptionalEflx(tem.sey[0]);
                     break;
    case C_SEY1 :    if (factor > 1) { tem.next(tem.sey[1]); }
      	             if (factor < 1) { tem.prev(tem.sey[1]); }
	             tem.displayOptionalEflx(tem.sey[1]);
	             break;
    case C_SEY2 :    if (factor > 1) { tem.next(tem.sey[2]); }
   	             if (factor < 1) { tem.prev(tem.sey[2]); }
	             tem.displayOptionalEflx(tem.sey[2]);
	             break;
    case C_RESET :   if (factor != 1) { reset += 1; }
                     if (reset > 1) { reset = 0; }
	             if (reset == 1) { cout << "reset = ON"; }
	             else { cout << "reset = OFF"; }
	             break;
    case C_MANFLG :  if (factor != 1) { manflag += 1; }
	             if (manflag > 1) { manflag = 0; }
	             if (manflag == 1) { cout << "manual = ON"; }
	             else { cout << "manual = OFF"; }
	             break;
    case C_AGSTATE : if (factor != 1) { tem.ag.state += 1; }
	             if (tem.ag.state > 1) { tem.ag.state = 0; }
	             if (tem.ag.state == 1) { cout << "agstate = ON"; }
	             else { cout << "agstate = OFF"; }
	             break;
   // adde for soil thermal model
    case C_STSTATE: if (factor != 1) stmflg += 1;
	       if (stmflg > 1) stmflg = 0;
       //        delline();
               if (stmflg == 1) cout << "ststate = ON";
               else cout << "ststate = OFF";
	       break;
    case C_TREEVEG: if (factor != 1) treeveg +=1;
             if (treeveg >1) treeveg =0;
              if (treeveg == 1) cout << "PFT2-cali";   // for choosing calibration for tree or moss
              else cout << "PFT1-cali";
           break;

    case C_CALWIND : if (factor != 1) { tem.calwind += 1; }
	             if (tem.calwind > 2) { tem.calwind = 1; }
	             if (tem.calwind == 1)
                     {
	               tem.firstcal = 1;
	               cout << "screen = WBM" << endl;
	               tem.topwind = 0;
	             }
	             else
                     {
	               tem.firstcal = 0;
	               cout << "screen = TEM" << endl;
	               tem.topwind = 1;
	             }
	             break;
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void wcalibrate(int t)
{

 int i; //
  double invar;
  int calvar[1];
  double factor = 1.0;
  tem.endgrid = 0;

  window(43,23,80,23);
  printf("RS = %1d             WTOL = %8.6lf",reset,tem.wtol);

  window(44,1,50,1);
// fixed bugs by Q. Z. 10/16/2000

//  tem.displayOptionalWflx(tem.swy[ESY+0]);
  tem.displayOptionalWflx(tem.swy[0]);

  window(51,1,57,1);
//  tem.displayOptionalWflx(tem.swy[ESY+1]);
 tem.displayOptionalWflx(tem.swy[1]);

  window(58,1,65,1);
//  tem.displayOptionalWflx(tem.swy[ESY+2]);
  tem.displayOptionalWflx(tem.swy[2]);

  window(66,1,73,1);
//  tem.displayOptionalWflx(tem.swy[ESY+3]);
  tem.displayOptionalWflx(tem.swy[3]);

  window(74,1,80,1);
//  tem.displayOptionalWflx(tem.swy[ESY+4]);
  tem.displayOptionalWflx(tem.swy[4]);

  window(1,23,15,23);
  cout << "CALIBRATION:   " << endl;

  if (tem.firstcal == 1)
  {
    ivar = 1;
    tem.firstcal = 0;
  }

  window(16,23,34,23);
  delline();

  switch(ivar)
  {
    case 1: printf("INT TOL = %8.6lf ",tem.tol); break;
    case 2: printf("NEW TEXTURE? (%4.2lf)",tem.soil.psiplusc); break;
    case 3: printf("NEW CLIMATE?       "); break;
    case 4: tem.displayOptionalWflx(tem.swy[ESY+0]);
	    break;
    case 5: tem.displayOptionalWflx(tem.swy[ESY+1]);
	    break;
    case 6: tem.displayOptionalWflx(tem.swy[ESY+2]);
	    break;
    case 7: tem.displayOptionalWflx(tem.swy[ESY+3]);
	    break;
    case 8: tem.displayOptionalWflx(tem.swy[ESY+4]);
	    break;
    case 9: if (reset == 1) { cout << "reset = ON" << endl; }
	    else { cout << "reset = OFF" << endl; }
            break;
    case 10: if (factor != 1.0) { tem.calwind += 1; }
    	     if (tem.calwind > 2) { tem.calwind = 1; }
	     if (tem.calwind == 2)
             {
	       tem.firstcal = 1;
	       delline();
	       cout << "screen = TEM" << endl;
	       tem.topwind = 0;
	     }
	     else
             {
	       tem.firstcal = 0;
	       delline();
	       cout << "screen = WBM" << endl;
	       tem.topwind = 1;
	     }
	     break;
  }

  while (t == 79 || t == 72 || t == 75 || t == 77 || t == 80 || t == 82)
  {

    if (tem.tolbomb == 0 && tem.intbomb == 0)
    {
      t = getch(); t = getch();
    }
    else
    {
      if (tem.tolbomb == 1)
      {
	tem.retry = 1;
	tem.tolbomb = 0;
	tem.nattempt += 1;
	ivar = 1;
	if (tem.nattempt < tem.maxnrun)
        {
	  tem.tol /= 10.0;
	  t = 0;
	}
	else
        {
	  tem.tol = tem.inittol;
	  tem.nattempt = 0;
	}
      }
      if (tem.intbomb == 1)
      {
	tem.retry = 1;
	tem.intbomb = 0;
	ivar = 1;
      }
    }

    if (tem.intbomb == 1) { ivar = 1; }

    calvar[0] =  tem.tol;

    factor = 1.0;
    window(35,23,42,23);
    printf(" ");
    switch (t)
    {
      case 79: break;
      case 72: factor = 1.01; break;
      case 80: factor = 0.99; break;
      case 77: ++ivar; break;
      case 75: --ivar; break;
      case 82: if (ivar <= 1)
               {
		 cin >> invar;
		 calvar[ivar-1] = invar;
	       }
	       break;
      default: t = 0; break;
    }

    window(35,23,42,23);
    delline();

    window(16,23,34,23);
    delline();

    if (ivar < 1) { ivar = 10; }
    if (ivar > 10) { ivar = 1; }
    switch (ivar)
    {
      case 1:  tem.tol = calvar[ivar - 1] * factor;
               delline();
       	       printf("INT TOL = %8.6lf ",tem.tol);
	       break;
      case 2:  if (factor != 1)
               {
		 window(1,1,80,24);
		 gotoxy(1,1);
		 clrscr();
		 cout << "Enter proportion of silt plus clay: ";
		 cin >> tem.soil.psiplusc;
		 cout << "Enter percent clay: ";
		 cin >> tem.soil.pctclay;
		 tem.soil.pctsilt = (tem.soil.psiplusc * 100.0)
                                    - tem.soil.pctclay;
		 tem.soil.pctsand = (1.0 - tem.soil.psiplusc)* 100.0;
		 tem.soil.xtext(tem.veg.cmnt,
                                tem.soil.pctsilt, tem.soil.pctclay,vt);
		 newtext = 1;
		 reset = 1;
		 t = 0;
	       }
	       else
               {
		 delline();
	         printf("NEW TEXTURE? (%4.2lf)",tem.soil.psiplusc);
	       }
	       break;
      case 3:  if (factor != 1)
               {
		 window(1,1,80,24);
		 gotoxy(1,1);
		 clrscr();
		 calclmin(vt);
		 showclm(vt);
		 tem.soil.xtext(tem.veg.cmnt,
                                tem.soil.pctsilt, tem.soil.pctclay,vt);
		 newclm = 1;
		 reset = 1;
		 t = 0;
	       }
	       else
               {
		 delline();
		 cout << "NEW CLIMATE?       ";
	       }
	       break;
      case 4:  if (factor > 1) { tem.next(tem.swy[ESY+0]); }
	       if (factor < 1) { tem.prev(tem.swy[ESY+0]); }
	       delline();
	       tem.displayOptionalWflx(tem.swy[ESY+0]);
	       break;
      case 5:  if (factor > 1) { tem.next(tem.swy[ESY+1]); }
	       if (factor < 1) { tem.prev(tem.swy[ESY+1]); }
	       delline();
	       tem.displayOptionalWflx(tem.swy[ESY+1]);
	       break;
      case 6:  if (factor > 1) { tem.next(tem.swy[ESY+2]); }
	       if (factor < 1) { tem.prev(tem.swy[ESY+2]); }
	       delline();
	       tem.displayOptionalWflx(tem.swy[ESY+2]);
	       break;
      case 7:  if (factor > 1) { tem.next(tem.swy[ESY+3]); }
	       if (factor < 1) { tem.prev(tem.swy[ESY+3]); }
	       delline();
	       tem.displayOptionalWflx(tem.swy[ESY+3]);
	       break;
      case 8:  if (factor > 1) { tem.next(tem.swy[ESY+4]); }
	       if (factor < 1) { tem.prev(tem.swy[ESY+4]); }
	       delline();
	       tem.displayOptionalWflx(tem.swy[ESY+4]);
	       break;
      case 9:  if (factor != 1) { ++reset; }
	       if (reset > 1) { reset = 0; }
	       delline();
	       if (reset == 1) { cout << "reset = ON"; }
	       else { cout << "reset = OFF"; }
	       break;
      case 10: if (factor != 1) { tem.calwind += 1; }
	       if (tem.calwind > 2) { tem.calwind = 1; }
	       if (tem.calwind == 2)
               {
		 tem.firstcal = 1;
		 delline();
		 cout << "screen = TEM" << endl;
		 tem.topwind = 0;
	       }
	       else
               {
		 tem.firstcal = 0;
		 delline();
		 cout << "screen = WBM" << endl;
		 tem.topwind = 1;
	       }
	       break;
    }

    if (t == 0) break;
  }


  if (reset == 1 || tem.retry == 1)
  {

 // bugs fixed by Q. Z. to correct the error between shifting from TEM to WBM and visa verso, Oct 12 /2000

  for (i = 0; i < CYCLE; i++)
  {
    tem.atms.co2[i] = tem.atms.co2level;
  }
  tem.atms.initco2 = tem.atms.co2level;

 //  tem.soil.xtext(tem.veg.cmnt, tem.soil.pctsilt, tem.soil.pctclay);
 // added for 3-box hydrology
    for (vt=0; vt < NVT; vt++) { // added for dealing with moss
    tem.veg.cmnt= subsist[vt];

  tem.soil.hydm_xtext(tem.veg.cmnt,tem.soil.pctsilt,tem.soil.pctclay,vt); // modified for two boxes hydrology model

  tem.ECDsetELMNTstate(tem.veg.cmnt, tem.soil.psiplusc,vt);

  tem.setELMNTflux(vt);

  tem.setPrevState(tem.prevy, tem.y,vt);

  tem.setELMNTecd(kdinflg,tem.veg.cmnt, tem.soil.psiplusc,vt);
  // modified for hydrology model by QZ
  tem.setELMNTevap(stateflag, tem.veg.cmnt, tem.atms.pet, tem.atms.tair,tem.atms.vap,vt);

  tem.mintflag = 0;
  tem.tol = tem.inittol;
  tem.retry = 0;
  tem.endgrid = 0;

  tem.firstcal = 1;
  tem.topwind = 0;

  dyr = 0;
  ivar = 1;
  tem.endeq = 1;

  tem.retry = 0;
  tem.endgrid = 0;
  equilsol = 0;

  tem.ag.state = 0;
  tem.ag.prvstate = 0;

  tem.prevy[vt][tem.I_UNRMLF] = 0.5;
  tem.prevy[vt][tem.I_LEAF] = 0.5;

 fswitch =0 ;

  } // end of vt
 }

   if (newtext == 1 || newclm == 1)
  {
    window(1,1,80,24);
    gotoxy(1,1);
    clrscr();
    for (vt=0; vt < NVT; vt++) { // added for dealing with moss
    tem.veg.cmnt= subsist[vt];

    tem.setELMNTecd(kdinflg,tem.veg.cmnt,tem.soil.psiplusc,vt);
    tem.pcdisplayStext(vt);
    tem.topwind = 0;
    newclm = 0;
    newtext = 0;
   } // end of vt
  }

  if (equilsol == 1)
  {
    equilsol = 0;
    dyr = 0;
  }

};


