
/* **************************************************************
*****************************************************************
XTEM423.CPP - Terrestrial Ecosystem Model Version 4.2
*****************************************************************

Modifications:

19991028 - DWK added bug fixes found and implemented by Gregg
           Christopher and Jim Long

20000715 -QZ added fvap for hydrology model
*****************************************************************
************************************************************** */
//mao
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<fcntl.h>
#include<iomanip>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>
#include<string.h>
#include<time.h>
//#include<conio.h>
//#include<value.h>


using namespace std;
const int NVT = 2; // number of vegetation types for a grid cell, for instance, here we defile 0 is for tree and 1 for moss, Q. Z. added for moss component


const int CYCLE = 12;
const int MAXRTIME = 600000; // initial 600 by QZ

const int MAXPRED = 92 + 29 +10 ; // changed from 85 to 92 by DWK on 20000202 Change to adding 29 for hydrology model bt QZ; 10 for STM

const int MAXGRID = 1;
// const int MAXNUMROW = 360; // commented out by DWK on 20000202

const double MISSING = -999999.9;
const double ZERO = 0.000000;  // added by DWK on 20000615

long int kswitch; //added for soil thermal model
int stmflg;   // added for soil thermal model
int subsist[NVT]; // to store plant function types
int mossflag; // for fire disturbance, addition by Q.Z. 29/Nov/2000
ifstream fpara;
int years_after_burn =0;
int months=0;
int burned_before=0;
int burned_before1=0;
int months1=0;
int years_after_burn1=0;
int number_of_2=0;

#include "elmnt423.cpp"    //Elmnt Class
#include "telm423e.cpp"     //TEMelmnt Class

ofstream flog1;
Elmnt elmnt;
TEMelmnt telmnt[MAXGRID];

/* *************************************************************
		 Function Declarations
************************************************************* */

int askpred(char pvarname[MAXPRED][9], int spred);
void extrapolate(void);
void startclm(void);
void starttem(void);

// *************************************************************

int equil;
int RTIME;

int spinflag;
int numspin;
int spintime;
int totsptime;
int transtime;

int atmstotyr;
int atmsflag;
int atmsoutfg;
int temflag;
int stateflag;

int cldflag;
int natmspred;
int ntempred;
int totpred;
char predmap[MAXPRED][9];
int kdinflg;
int kdoutflg;

int fatalerr;
int end1;
int icount;
int glob_count;

FILE* flonlat;
FILE* fstxt;
FILE* felev;
FILE* ftveg;
FILE* fclds;
FILE* ftair;
FILE* fprec;
FILE* fsoiloggc;      			//wsr: add soil organic carbon input

//added for hydrology model by QZ
FILE* fvap;

FILE* fnirr;
FILE* fpar;
FILE* flulc;
FILE* fnpp;
FILE* fstate[MAXSTATE];
ofstream fclmpred[NUMATMS];
ofstream ftempred[MAXPRED];

FILE* fkdin;
ofstream fkdout;
KDdata kdparam;

/* *************************************************************
**********************START MAIN PROGRAM************************
************************************************************* */

int main(int argc, char** argv)
{

  int i;
 
  // Open log file
  fpara.open(argv[1]);
  flog1.open(argv[2]);
/* *************************************************************
  Run equilibrium simulation or transient simulation ?
************************************************************* */

  equil = 1;
  cout << endl << "Do you want to run the model only for steady state conditions ? " << endl;
  cout << " Enter 0 for transient simulation" << endl;
  cout << " Enter 1 for steady state simulation" << endl;
  fpara >> equil;

  flog1 << endl << "Do you want to run the model only for steady state conditions ? " << endl;
  flog1 << " Enter 0 for transient simulation" << endl;
  flog1 << " Enter 1 for steady state simulation" << endl;
  flog1 << "equil = " << equil << endl << endl;

  RTIME = 1;
  if (equil == 0)
  {
    // Start transient TEM run from equilibrium conditions (i.e. spinflag == 0)
    //  or with a "spin up" period to remove potential artifacts associated with
    //  changing from equilibrium conditions to transient conditions from model
    //  results

    cout << "Do you want to start the transient with a spin up period? " << endl;
    cout << "Enter 0 for no:" << endl;
    cout << "Enter 1 for yes: ";
    fpara  >> spinflag;

    flog1 << "Do you want to start the transient with a spin up period? " << endl;
    flog1 << "Enter 0 for no:" << endl;
    flog1 << "Enter 1 for yes: " << endl;
    flog1 << "spinflag = " << spinflag << endl << endl;

    totsptime = 0;
    if (spinflag == 1)
    {
      // Specify conditions for initializing TEM with a transient "spin up" period
      //   for a grid cell

      cout << "How many spins do you want in the spin up period? ";
      fpara >> numspin;
      flog1 << "How many spins do you want in the spin up period? " << endl;
      flog1 << "numspin = " << numspin << endl << endl;

      cout << "How many years per spin? ";
      fpara >> spintime;
      flog1 << "How many years per spin? " << endl;
      flog1 << "spintime = " << spintime << endl << endl;

      totsptime = spintime * numspin;
      flog1 << "totsptime = " << totsptime << endl << endl;
      RTIME += totsptime;
    }

    // Specify conditions for the "non-spin up" part of the transient TEM run

    cout << endl << "How many years do you run for transient simulations ? " << endl;
   fpara >> transtime;

    flog1 << endl << "How many years do you run for transient simulations ? " << endl;
    flog1 << "transtime = " << transtime << endl << endl;

    RTIME += transtime;
    flog1 << "RTIME = " << RTIME << endl << endl;
  }

/* *************************************************************
What models do you want to run?  - Open associated files
************************************************************* */

 //  STM (soil thermal model)?

  cout << endl << "Do you want to run the SOIL THERMAL MODEL (STM) model for soil temperatures?" << endl;
  cout << "  Enter 0 for No" << endl;
  cout << "  Enter 1 for Yes" << endl;
  fpara >> stmflg;

  flog1 << endl << "Do you want to run the SOIL THERMAL MODEL (STM) model for soil temperatures?"<< endl;
  flog1 << "  Enter 0 for No" << endl;
  flog1 << "  Enter 1 for Yes" << endl;
  flog1 << " stmflg = " << stmflg << endl << endl;
  if (stmflg == 1) {
// Get snow and soil layer dependent parameters, added for soil thermal model integration
  telmnt[0].tem.sthermal.getsnowecd(flog1);   //qtsp44a2.ecd
  telmnt[0].tem.sthermal.getsoillecd(flog1);  // qtla44a2.ecd
  telmnt[0].tem.sthermal.getsoiltecd(flog1);  // qtst44a2.ecd
  }

// POTSCLM model ?

  cout << endl << "Do you want to run the POTSCLM model for solar radiation variables?" << endl;
  cout << "  Enter 0 for No" << endl;
  cout << "  Enter 1 for Yes" << endl;
 fpara >> atmsflag;

  flog1 << endl << "Do you want to run the POTSCLM model for solar radiation variables?" << endl;
  flog1 << "  Enter 0 for No" << endl;
  flog1 << "  Enter 1 for Yes" << endl;
  flog1 << " atmsflag = " << atmsflag << endl << endl;

  totpred = natmspred = 0;
  if (atmsflag == 1) { startclm(); }

// TEM ?

  cout << endl << "Do you want to run the terrestrial ecosystem model (TEM)?" << endl;
  cout << "  Enter 0 for No" << endl;
  cout << "  Enter 1 for Yes" << endl;
   fpara >> temflag;

  flog1 << endl << "Do you want to run the terrestrial ecosystem model (TEM)?" << endl;
  flog1 << "  Enter 0 for No" << endl;
  flog1 << "  Enter 1 for Yes" << endl;
  flog1 << "temflag = " << temflag << endl << endl;

  ntempred = 0;
  if (temflag == 1) { starttem(); }

  // Use all records in GIS data sets (i.e. desired coverage)?

  elmnt.ask(flog1);

  // Run Potsclm and TTEM modules over desired region (modules could potentially
  // run over several grid cells simultaneously with the use of extrapolate()

  extrapolate();

// Finished processing all elements - close open files

  if (fatalerr != 0)
  {
    if (elmnt.grdcnt != -99 && elmnt.count <= elmnt.grdcnt)
    {
      cout << "FATAL ERROR! Program Terminated" << endl;
    }
    flog1 << "FATAL ERROR! Program Terminated" << endl;
  }
  else
  {
    cout << "Extrapolation successfully completed - Congratulations!" << endl;
    flog1 << "Extrapolation successfully completed - Congratulations!" << endl;
  }

  if (atmsflag == 1)
  {
    fclose(fclds);
    if (telmnt[0].lonlatflag == 0) { fclose(flonlat); }
    if (atmsoutfg == 1)
    {
      for (i = 0; i < natmspred; i++) { fclmpred[i].close(); }
    }
  }

  if (temflag == 1)
  {
    fclose(ftveg);
    fclose(fstxt);
    fclose(felev);
    if (atmsflag == 0)
    {
      fclose(fnirr);
      fclose(fpar);
    }
    fclose(ftair);
    fclose(fprec);
	  //added for hydrology model
    fclose(fvap);
    // end of adding
    fclose(fsoiloggc);

    if(telmnt[0].tem.ag.tlulcflag == 1)
    {
      fclose(flulc);
      if(telmnt[0].tem.ag.RAP0flag == 0) { fclose(fnpp); }
    }
    for (i = 0; i < ntempred; i++) { ftempred[i].close(); }
    if (stateflag == 1)
    {
      for (i = 0; i < MAXSTATE; i++) { fclose(fstate[i]); }
    }
  }
  fpara.close();
  flog1.close();
  return 1;
};

/* *************************************************************
*********************END OF MAIN PROGRAM************************
************************************************************* */


/* **************************************************************
************************************************************** */

int askpred(char pvarname[MAXPRED][9], int spred)
{
  int i;
  int j;
  int k;
  int t;
  int cnt;
  int numpred;
  int length;

  int posspred = MAXPRED;
  if (atmsoutfg == 1 && temflag == 0) { posspred = NUMATMS+1; }

  cout << endl << endl << "           POSSIBLE OUTPUT VARIABLES:" << endl << endl;
  for (i = 0; i < (posspred/5); i++)
  {
    for (t = 0; t < 5; t++)
    {
      cout << pvarname[(5*i)+t] << " ";
    }
    cout << endl;
  }

  flog1 << endl << endl << "           POSSIBLE OUTPUT VARIABLES:" << endl << endl;
  for (i = 0; i < (posspred/5); i++)
  {
     for (t = 0; t < 5; t++)
    {
      flog1 << pvarname[(5*i)+t] << ' ';
    }
    flog1 << endl;
  }

  cout << endl << endl << "How many variables are to be mapped (max " << posspred << ") in output files?  ";
  fpara >> numpred;
  cout << numpred << endl;

  flog1 << endl << endl << "How many variables are to be mapped (max " << posspred << ") in output files?" << numpred << endl << endl;

  cout << "Please enter output variable: " << endl;
  flog1 << "Please enter output variable: " << endl;

  for (i = 0; i < numpred; i++)
  {
    cnt = i + 1;
    k = i+spred;
    cout << cnt << " ";
    fpara >> predmap[k];
    length = strlen(predmap[k]);
    for (j = 0; j < length; j++) { predmap[k][j] = toupper(predmap[k][j]); }
    flog1 << cnt << " " << predmap[k] << endl;
  }

  return numpred;

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void extrapolate(void)
{

  int i;
  int j;
  int k;
  int dm;
  int maxrec;
  int dyr;
  int itype;
//  double totamnt[NUMEQ+2];

  for (i = 1; i < MAXGRID; i++)
  {
    if (atmsflag == 1)
    {
      telmnt[i].clm.tcldsflag = telmnt[0].clm.tcldsflag;
      telmnt[i].lonlatflag = telmnt[0].lonlatflag;
    }

    if (temflag == 1)
    {
      telmnt[i].tem.atms.tco2flag = telmnt[0].tem.atms.tco2flag;
      telmnt[i].lonlatflag = telmnt[0].lonlatflag;
      telmnt[i].tem.atms.ttairflag = telmnt[0].tem.atms.ttairflag;
      telmnt[i].tem.atms.tprecflag = telmnt[0].tem.atms.tprecflag;
// added for hydrology model
      telmnt[i].tem.hyd.vapflag = telmnt[0].tem.hyd.vapflag;

      telmnt[i].tem.ag.tlulcflag = telmnt[0].tem.ag.tlulcflag;

// Soil texture specific TEM parameters

      telmnt[i].tem.soil.pctpora = telmnt[0].tem.soil.pctpora;
      telmnt[i].tem.soil.pctporb = telmnt[0].tem.soil.pctporb;
      telmnt[i].tem.soil.fldcapa = telmnt[0].tem.soil.fldcapa;
      telmnt[i].tem.soil.fldcapb = telmnt[0].tem.soil.fldcapb;
      telmnt[i].tem.soil.wiltpta = telmnt[0].tem.soil.wiltpta;
      telmnt[i].tem.soil.wiltptb = telmnt[0].tem.soil.wiltptb;

// Initialize CO2 for all grid cells

      telmnt[i].tem.atms.initco2 = telmnt[0].tem.atms.initco2;
      telmnt[i].tem.atms.co2level = telmnt[0].tem.atms.co2level;
      if (equil == 0 && telmnt[0].tem.atms.tco2flag == 1)
      {
	for (j = 0; j < RTIME; j++)
        {
	  for (dm = 0; dm < CYCLE; dm++)
          {
 	    telmnt[i].tem.atms.tco2[j][dm] = telmnt[0].tem.atms.tco2[j][dm];
	  }
	}
      }
      for (j = 0; j < NUMVEG; j++)
      {
	telmnt[i].tem.soil.rootza[j] = telmnt[0].tem.soil.rootza[j];
	telmnt[i].tem.soil.rootzb[j] = telmnt[0].tem.soil.rootzb[j];
	telmnt[i].tem.soil.rootzc[j] = telmnt[0].tem.soil.rootzc[j];
	telmnt[i].tem.soil.minrootz[j] = telmnt[0].tem.soil.minrootz[j];
	telmnt[i].tem.veg.numtype[j] = telmnt[0].tem.veg.numtype[j];
	telmnt[i].tem.veg.minleaf[j] = telmnt[0].tem.veg.minleaf[j];
	telmnt[i].tem.veg.aleaf[j] = telmnt[0].tem.veg.aleaf[j];
	telmnt[i].tem.veg.bleaf[j] = telmnt[0].tem.veg.bleaf[j];
	telmnt[i].tem.veg.cleaf[j] = telmnt[0].tem.veg.cleaf[j];

//Site specific TEM parameters

	telmnt[i].tem.vegca[j] = telmnt[0].tem.vegca[j];
	telmnt[i].tem.vegcb[j] = telmnt[0].tem.vegcb[j];
	telmnt[i].tem.strna[j] = telmnt[0].tem.strna[j];
	telmnt[i].tem.strnb[j] = telmnt[0].tem.strnb[j];
	telmnt[i].tem.solca[j] = telmnt[0].tem.solca[j];
	telmnt[i].tem.solcb[j] = telmnt[0].tem.solcb[j];
	telmnt[i].tem.solna[j] = telmnt[0].tem.solna[j];
	telmnt[i].tem.solnb[j] = telmnt[0].tem.solnb[j];
	telmnt[i].tem.avlna[j] = telmnt[0].tem.avlna[j];
	telmnt[i].tem.avlnb[j] = telmnt[0].tem.avlnb[j];
	telmnt[i].tem.stona[j] = telmnt[0].tem.stona[j];
	telmnt[i].tem.stonb[j] = telmnt[0].tem.stonb[j];
	telmnt[i].tem.veg.unleaf12[j] = telmnt[0].tem.veg.unleaf12[j];
	telmnt[i].tem.veg.prvleafmx[j] = telmnt[0].tem.veg.prvleafmx[j];
	telmnt[i].tem.veg.cmaxcut[j] = telmnt[0].tem.veg.cmaxcut[j];
	telmnt[i].tem.veg.cmax1a[j] =  telmnt[0].tem.veg.cmax1a[j];
	telmnt[i].tem.veg.cmax1b[j] = telmnt[0].tem.veg.cmax1b[j];
	telmnt[i].tem.veg.cmax2a[j] = telmnt[0].tem.veg.cmax2a[j];
	telmnt[i].tem.veg.cmax2b[j] = telmnt[0].tem.veg.cmax2b[j];
	telmnt[i].tem.veg.cfall[j] = telmnt[0].tem.veg.cfall[j];
	telmnt[i].tem.veg.kra[j] = telmnt[0].tem.veg.kra[j];
	telmnt[i].tem.veg.krb[j] = telmnt[0].tem.veg.krb[j];
	telmnt[i].tem.microbe.kda[j] = telmnt[0].tem.microbe.kda[j];
	telmnt[i].tem.microbe.kdb[j] = telmnt[0].tem.microbe.kdb[j];
	telmnt[i].tem.microbe.lcclnc[j] = telmnt[0].tem.microbe.lcclnc[j];
	telmnt[i].tem.microbe.propftos[j] = telmnt[0].tem.microbe.propftos[j];
	telmnt[i].tem.veg.vsmmina[j] = telmnt[0].tem.veg.vsmmina[j];
	telmnt[i].tem.veg.vsmminb[j] = telmnt[0].tem.veg.vsmminb[j];
	telmnt[i].tem.veg.nmaxcut[j] = telmnt[0].tem.veg.nmaxcut[j];
	telmnt[i].tem.veg.nmax1a[j] = telmnt[0].tem.veg.nmax1a[j];
	telmnt[i].tem.veg.nmax1b[j] = telmnt[0].tem.veg.nmax1b[j];
	telmnt[i].tem.veg.nmax2a[j] = telmnt[0].tem.veg.nmax2a[j];
	telmnt[i].tem.veg.nmax2b[j] = telmnt[0].tem.veg.nmax2b[j];
	telmnt[i].tem.veg.nfall[j] = telmnt[0].tem.veg.nfall[j];
	telmnt[i].tem.microbe.nupa[j] = telmnt[0].tem.microbe.nupa[j];
	telmnt[i].tem.microbe.nupb[j] = telmnt[0].tem.microbe.nupb[j];
	telmnt[i].tem.soil.nloss[j] = telmnt[0].tem.soil.nloss[j];
	telmnt[i].tem.microbe.nfixpar[j] = telmnt[0].tem.microbe.nfixpar[j];
	telmnt[i].tem.veg.cneven[j] = telmnt[0].tem.veg.cneven[j];
	telmnt[i].tem.veg.cnmin[j] = telmnt[0].tem.veg.cnmin[j];
	telmnt[i].tem.veg.c2na[j] = telmnt[0].tem.veg.c2na[j];
	telmnt[i].tem.veg.c2nb[j] = telmnt[0].tem.veg.c2nb[j];
	telmnt[i].tem.veg.c2nmin[j] = telmnt[0].tem.veg.c2nmin[j];
	telmnt[i].tem.microbe.cnsoil[j] = telmnt[0].tem.microbe.cnsoil[j];

        // Vegetation specific TEM parameters

	telmnt[i].tem.veg.kc[j] = telmnt[0].tem.veg.kc[j];
	telmnt[i].tem.veg.ki[j] = telmnt[0].tem.veg.ki[j];
	telmnt[i].tem.veg.gva[j] = telmnt[0].tem.veg.gva[j];
	telmnt[i].tem.veg.tmin[j] = telmnt[0].tem.veg.tmin[j];
	telmnt[i].tem.veg.toptmin[j] = telmnt[0].tem.veg.toptmin[j];
	telmnt[i].tem.veg.toptmax[j] = telmnt[0].tem.veg.toptmax[j];
	telmnt[i].tem.veg.tmax[j] = telmnt[0].tem.veg.tmax[j];
	telmnt[i].tem.veg.raq10a0[j] = telmnt[0].tem.veg.raq10a0[j];
	telmnt[i].tem.veg.raq10a1[j] = telmnt[0].tem.veg.raq10a1[j];
	telmnt[i].tem.veg.raq10a2[j] = telmnt[0].tem.veg.raq10a2[j];
	telmnt[i].tem.veg.raq10a3[j] = telmnt[0].tem.veg.raq10a3[j];
	telmnt[i].tem.veg.kn1[j] = telmnt[0].tem.veg.kn1[j];
	telmnt[i].tem.veg.labncon[j] = telmnt[0].tem.veg.labncon[j];

        telmnt[i].tem.veg.leafmxc[j] = telmnt[0].tem.veg.leafmxc[j];
	telmnt[i].tem.veg.kleafc[j] = telmnt[0].tem.veg.kleafc[j];
	telmnt[i].tem.veg.sla[j] = telmnt[0].tem.veg.sla[j];
	telmnt[i].tem.veg.cov[j] = telmnt[0].tem.veg.cov[j];
	telmnt[i].tem.veg.fpcmax[j] = telmnt[0].tem.veg.fpcmax[j];

//Vegetation specific TEM parameters for microbes

	telmnt[i].tem.microbe.rhq10[j] = telmnt[0].tem.microbe.rhq10[j];
	telmnt[i].tem.microbe.kn2[j] = telmnt[0].tem.microbe.kn2[j];
	telmnt[i].tem.microbe.moistmin[j] = telmnt[0].tem.microbe.moistmin[j];
	telmnt[i].tem.microbe.moistopt[j] = telmnt[0].tem.microbe.moistopt[j];
	telmnt[i].tem.microbe.moistmax[j] = telmnt[0].tem.microbe.moistmax[j];

//added for hydrology model, penman-monteith parameters

 	telmnt[i].tem.hyd.inter_coef[j] = telmnt[0].tem.hyd.inter_coef[j];
	telmnt[i].tem.hyd.EXT[j] = telmnt[0].tem.hyd.EXT[j];
	telmnt[i].tem.hyd.LWPmin[j] = telmnt[0].tem.hyd.LWPmin[j];
	telmnt[i].tem.hyd.CCmax[j] = telmnt[0].tem.hyd.CCmax[j];
	telmnt[i].tem.hyd.LWPc[j] = telmnt[0].tem.hyd.LWPc[j];
	telmnt[i].tem.hyd.SLOPEcc[j] = telmnt[0].tem.hyd.SLOPEcc[j];

        for (k=0; k < NUMMSAC; k++)
        {
	  telmnt[i].tem.veg.subtype[j][k] = telmnt[0].tem.veg.subtype[j][k];
	  telmnt[i].tem.veg.pcttype[j][k] = telmnt[0].tem.veg.pcttype[j][k];
	}

	telmnt[i].tem.ag.slashpar[j] = telmnt[0].tem.ag.slashpar[j];
	telmnt[i].tem.ag.vconvert[j] = telmnt[0].tem.ag.vconvert[j];
	telmnt[i].tem.ag.vrespar[j] = telmnt[0].tem.ag.vrespar[j];
	telmnt[i].tem.ag.sconvert[j] = telmnt[0].tem.ag.sconvert[j];
	telmnt[i].tem.ag.prod10par[j] = telmnt[0].tem.ag.prod10par[j];
	telmnt[i].tem.ag.prod100par[j] = telmnt[0].tem.ag.prod100par[j];
	telmnt[i].tem.ag.nvretconv[j] = telmnt[0].tem.ag.nvretconv[j];
	telmnt[i].tem.ag.nsretconv[j] = telmnt[0].tem.ag.nsretconv[j];
      }
    }
  }

  end1 = 1;
  fatalerr = 0;

//* *************************************************************

// Skip to the desired record in the GIS data sets

  if (elmnt.strtflag == 0)
  {
    for (i = 0; i < elmnt.numskip; i++)
    {
      if (atmsflag == 1)
      {
	end1 = telmnt[0].atmsgisin(fclds,cldflag,flonlat,telmnt[0].lonlatflag,
                                   numspin,spintime,RTIME,fatalerr,flog1);
	if (end1 == -1) { exit(0); }
      }

      if (temflag == 1)
      {
	end1 = telmnt[0].veggisin(flog1,fatalerr,atmsflag,temflag,ftveg);
	if (end1 == -1) { exit(0); }

        for (itype = 0; itype < telmnt[icount].maxtype; itype++)
        {
	  if (temflag == 1 && fatalerr == 0)
          {

   // mofified for soil thermal model, fvap for hydrological model

	    end1 = telmnt[0].temgisin(flog1,fatalerr,atmsflag,          //wsr: add soiloggc
                                      itype,telmnt[icount].maxtype,
                                      kdinflg,stateflag,
	     	                      fstxt,fsoiloggc,felev,
                                      fnirr,fpar,
                                      ftair,fprec, fvap,
                                      flulc,fnpp,
                                      fkdin,fstate,
                                      numspin,spintime,RTIME);
            if (end1 == -1) { exit(0); }
	  }
	}
      }
    }
  }

//************************************************************** */

  cout << endl;
  flog1 << endl << endl;

  telmnt[0].col = MISSING;
  telmnt[0].row = MISSING;
  telmnt[0].tem.totyr = -99;

  elmnt.show(flog1, telmnt[0].col, telmnt[0].row,
             telmnt[0].tem.totyr, telmnt[0].tem.inittol);

  while (end1 != -1 && fatalerr == 0)
  {

    icount = 0;
    while (icount < MAXGRID && end1 != -1)
    {

// Get spatially-explicit cloudiness data   *********************

      if (atmsflag == 1)
      {
 	end1 = telmnt[icount].atmsgisin(fclds,cldflag,flonlat,
                                        telmnt[icount].lonlatflag,
                                        numspin,spintime,RTIME,fatalerr,flog1);
      }

// Get spatially-explicit vegetation data   *********************

      if (end1 != -1 && temflag == 1)
      {
	end1 = telmnt[icount].veggisin(flog1,fatalerr,atmsflag,temflag,ftveg);
      }

/* *************************************************************
		BEGIN VEGETATION MOSAIC LOOP
************************************************************* */

      if (end1 != -1 && temflag == 1)
      {
        for (itype = 0; itype < telmnt[icount].maxtype; itype++)
        {
          // Get spatially-explicit input data for TEM

	  if (end1 != -1)
          {														//wsr: add soiloggc
	    end1 = telmnt[icount].temgisin(flog1,fatalerr,
                                           atmsflag,itype,
                                           telmnt[icount].maxtype,
                                           kdinflg,stateflag,
					   fstxt,fsoiloggc,felev,
                                           fnirr,fpar,
                                           ftair,fprec,fvap,
                                           flulc,fnpp,
                                           fkdin,fstate,
                                           numspin,spintime,RTIME);
          }
	}
      }

      if (end1 != -1) { ++icount; }
    }

    maxrec = icount;

    for (i = 0; i < maxrec; i++)
    {

      // Run TEM for "i" grid cells
  int vtt=0;
  
//      for (vtt =0 ; vtt< 2; vtt++) { // modified for two plant function types QZ.

      telmnt[i].runtem(flog1,predmap,
                       cldflag,atmsflag,atmsoutfg, natmspred, fclmpred,
                       temflag,kdinflg,ntempred,ftempred,
                       stateflag,equil,totsptime, RTIME, vtt);
  //                     }




    }
	
/* *************************************************************
	  Write georeferenced data to output files
************************************************************* */

    for (i = 0; i < maxrec; i++)
    {
      for (dyr = 0; dyr < telmnt[i].outyr; dyr++)
      {
	     for (itype = 0; itype < telmnt[i].maxtype; itype++)
        {
	       if (temflag == 1)
          {
            // Write out spatially explicit kd parameters

	         if (dyr == 0 && kdoutflg == 1)
            {
	           kdparam.outdel(fkdout, telmnt[i].col, telmnt[i].row,
                             telmnt[i].tem.veg.temveg,
                             telmnt[i].tem.veg.subtype[telmnt[i].mez][itype],
                             telmnt[i].tem.microbe.kdsave[itype],
                             telmnt[i].contnent);
	         }
	       }
	     }
      }

      if (telmnt[i].tem.intflag > 0)
      {
        if (elmnt.count < elmnt.grdcnt)
        {
	       cout << "Integration terminated before attaining tolerance level" << endl;
	     }
	     flog1 << "Integration terminated before attaining tolerance level" << endl;
      }
      elmnt.show(flog1, telmnt[i].col, telmnt[i].row,
                 telmnt[i].tem.totyr, telmnt[i].tem.tol);
      if (elmnt.stopflag == 1) { end1 = -1; }
    }
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void startclm(void)
{

  int i;
  char ifilename[25];
  char clmpredfile[25];

  cldflag = 1;
  cout << "Do you have spatially explicit solar radiation data or cloudiness data?:" << endl;
  cout << "Enter 0 for solar radiation data (W/ sq. m):" << endl;
  cout << "Enter 1 for cloudiness data (percent cloudiness): ";
  fpara >> cldflag;
  flog1 << "Do you have spatially explicit solar radiation data or cloudiness data?:" << endl;
  flog1 << "Enter 0 for solar radiation data (W/ sq. m):" << endl;
  flog1 << "Enter 1 for cloudiness data (percent cloudiness): " << endl;
  flog1 << "cldflag = " << cldflag << endl;

  telmnt[0].clm.tcldsflag = 0;
  if (cldflag == 1)
  {
    cout << "Do you have transient cloudiness data?:" << endl;
    cout << "Enter 0 for no:" << endl;
    cout << "Enter 1 for yes: ";
    fpara >> telmnt[0].clm.tcldsflag;

    flog1 << endl << "Do you have transient cloudiness data?:" << endl;
    flog1 << "Enter 0 for no:" << endl;
    flog1 << "Enter 1 for yes: " << endl;
    flog1 << "telmnt[0].clm.tcldsflag = " << telmnt[0].clm.tcldsflag << endl << endl;

    cout << "Please enter the name of the file containing the mean monthly cloudiness data: " << endl;
    cout << "               (e.g., CLDS.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the mean monthly cloudiness data: " << endl;
    flog1 << "               (e.g., CLDS.GIS) " << endl;
    flog1 << ifilename << endl << endl;
  }
  else
  {
  
      cout << "Do you have transient nirr data?:" << endl;
    cout << "Enter 0 for no:" << endl;
    cout << "Enter 1 for yes: ";
    fpara >> telmnt[0].clm.nirrflag;

    flog1 << endl << "Do you have transient nirr data?:" << endl;
    flog1 << "Enter 0 for no:" << endl;
    flog1 << "Enter 1 for yes: " << endl;
    flog1 << "telmnt[0].clm.nirrflag = " << telmnt[0].clm.nirrflag << endl << endl;
    
    if(telmnt[0].clm.nirrflag==1)
    {
    cout << "Please enter the name of the file containing the mean monthly solar radiation data: " << endl;
    cout << "               (e.g., NIRR.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the mean monthly solar radiation data: " << endl;
    flog1 << "               (e.g., NIRR.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    }
  }
telmnt[0].tem.atms.RGflag = cldflag;

  fclds = fopen(ifilename, "r");

  if (!fclds)
  {
    cerr << "\nCannot open " << ifilename << " for data input" << endl;
    exit(-1);
  }

  cout << "How do you locate your grid cells?" << endl;
  cout << "Enter 0 for column/row:" << endl;
  cout << "Enter 1 for longitude/latitude: ";
  fpara >> telmnt[0].lonlatflag;

  flog1 << "How do you locate your grid cells?" << endl;
  flog1 << "Enter 0 for column/row:" << endl;
  flog1 << "Enter 1 for longitude/latitude: " << telmnt[0].lonlatflag << endl << endl;

  if (telmnt[0].lonlatflag == 0)
  {
    cout << "Please enter the name of the file containing the latitude data: " << endl;
    cout << "               (e.g., LAT.GIS) " << endl;
    fpara >> ifilename;

    flog1 << "Please enter the name of the file containing the latitude data: " << endl;
    flog1 << "               (e.g., LAT.GIS) " << endl;
    flog1 << ifilename << endl << endl;

    flonlat = fopen(ifilename, "r");

    if (!flonlat)
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }
  }

  atmsoutfg = 0;
  cout << endl << "Do you wish output data from the irradiance model? " << endl;
  cout << "  Enter 0 for no" << endl;
  cout << "  Enter 1 for yes: ";
  fpara >> atmsoutfg;

  flog1 << endl << "Do you wish output data from the irradiance model? " << endl;
  flog1 << "  Enter 0 for no" << endl;
  flog1 << "  Enter 1 for yes: " << endl;
  flog1 << "atmsoutfg = " << atmsoutfg << endl << endl;

  if (atmsoutfg == 1)
  {
    totpred = natmspred = askpred(telmnt[0].clm.predstr, totpred);

    for (i = 0; i < natmspred; i++)
    {
      cout << endl << "Enter the name of the OUTPUT file to contain ";
      cout << predmap[i] << ":  ";
      fpara >> clmpredfile;
      flog1 << endl << "Enter the name of the OUTPUT file to contain ";
      flog1 << predmap[i] << ":  " << clmpredfile << endl;
      fclmpred[i].open(clmpredfile, ios::app);
    }
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void starttem(void)
{

  int i;
  int numcmnt;
  char ifilename[40];
  char tempredfile[25];

  telmnt[0].tem.initrun(flog1,equil);
  telmnt[0].tem.ask(flog1);


// Get soil texture dependent parameters

  telmnt[0].tem.soil.getecd(flog1);


// Get vegetation type dependent parameters

  telmnt[0].tem.soil.getrootz(flog1);
  telmnt[0].tem.veg.getecd(flog1);
  telmnt[0].tem.veg.getleafecd(flog1);
  telmnt[0].tem.microbe.getvegecd(flog1);

 // added for penmon equation hydrology model
  telmnt[0].tem.hyd.getecd(flog1);

//Get calibration site specific parameters

  cout << "Please enter the number of community types with calibration data:";
  fpara >> numcmnt;

  flog1 << endl << endl << "Please enter the number of community types with calibration data:";
  flog1 << numcmnt << endl;

  telmnt[0].tem.getsitecd(numcmnt, flog1);

  telmnt[0].tem.getsitecd_p(numcmnt, flog1);

// Get vegetation mosaic information

  telmnt[0].tem.veg.getvtype(flog1);


  cout << "Please enter the name of the file containing the soil texture data:";
  cout << endl;
  cout << "               (e.g., TEXTURE.GIS) " << endl;
  fpara >> ifilename;
  flog1 << "Please enter the name of the file containing the soil texture data:";
  flog1 << endl;
  flog1 << "               (e.g., TEXTURE.GIS) " << endl;
  flog1 << ifilename << endl << endl;
  fstxt = fopen(ifilename, "r");

  if (!fstxt)
  {
    cerr << "\nCannot open " << ifilename << " for data input" << endl;
    exit(-1);
  }

  cout << "Please enter the name of the file containing the vegetation data: " << endl;
  cout << "               (e.g., TEMVEG.GIS) " << endl;
  fpara >> ifilename;
  flog1 << "Please enter the name of the file containing the vegetation data: " << endl;
  flog1 << "               (e.g., TEMVEG.GIS) " << endl;
  flog1 << ifilename << endl << endl;
  ftveg = fopen(ifilename, "r");

  if (!ftveg)
  {
    cerr << "\nCannot open " << ifilename << " for data input" << endl;
    exit(-1);
  }

  cout << "Please enter the name of the file containing the elevation data: " << endl;
  cout << "               (e.g., ELEV.GIS) " << endl;
  fpara >> ifilename;
  flog1 << "Please enter the name of the file containing the elevation data: " << endl;
  flog1 << "               (e.g., ELEV.GIS) " << endl;
  flog1 << ifilename << endl << endl;
  felev = fopen(ifilename, "r");

  if (!felev)
  {
    cerr << "\nCannot open " << ifilename << " for data input" << endl;
    exit(-1);
  }

//****************************************************wsr: enter initial soil carbon data 

  telmnt[0].tem.atms.tpeatflag = 0;
  cout << endl << "Is the vegetation peatland? " << endl;
  cout << "  Enter 0 for no" << endl;
  cout << "  Enter 1 for yes: ";
  fpara >>  telmnt[0].tem.atms.tpeatflag;

  flog1 << endl << "Is the vegetation peatland? " << endl;
  flog1 << "  Enter 0 for no" << endl;
  flog1 << "  Enter 1 for yes: " << endl;
  flog1 << " telmnt[0].tem.atms.tpeatflag = " << telmnt[0].tem.atms.tpeatflag << endl << endl;



if ( telmnt[0].tem.atms.tpeatflag==1)
{

telmnt[0].tem.atms.tsoiloggcflag = 0;

if (equil == 0 )
{
cout << "Do you have soil organic carbon data from previous simulation?" << endl;
cout << "Enter 0 for No:" << endl;
cout << "Enter 1 for Yes:";
fpara >> telmnt[0].tem.atms.tsoiloggcflag;

    flog1 << endl << "Do you have soil organic carbon data from previous simulation?" << endl;
    flog1 << "Enter 0 for No:" << endl;
    flog1 << "Enter 1 for Yes: " << endl;
    flog1 << "telmnt[0].tem.atms.tsoiloggcflag = " << telmnt[0].tem.atms.tsoiloggcflag << endl << endl;
  }
  
  if (telmnt[0].tem.atms.tsoiloggcflag == 1)
{
  
cout << "Please enter the name of the file containing the soil organic carbon data: " << endl;
  cout << "               (e.g., SOILOGGC.GIS) " << endl;
  fpara >> ifilename;
  flog1 << "Please enter the name of the file containing the soil organic carbon data: " << endl;
  flog1 << "               (e.g., SOILOGGC.GIS) " << endl;
  flog1 << ifilename << endl << endl;
  fsoiloggc = fopen(ifilename, "r");

  if (!fsoiloggc)
  {
    cerr << "\nCannot open " << ifilename << " for data input" << endl;
    exit(-1);
  }
}
}//end of peatland flag
//************************************************end of adding


  telmnt[0].tem.atms.ttairflag = 0;
  if (equil == 0)
  {
    cout << "Do you have transient air temperature data?:" << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";
    fpara >> telmnt[0].tem.atms.ttairflag;

    flog1 << endl << "Do you have transient air temperature data?:" << endl;
    flog1 << "Enter 0 for No:" << endl;
    flog1 << "Enter 1 for Yes: " << endl;
    flog1 << "telmnt[0].tem.atms.ttairflag = " << telmnt[0].tem.atms.ttairflag << endl << endl;
  }

  cout << "Please enter the name of the file containing the mean monthly air temperature data: " << endl;
  cout << "               (e.g., TAIR.GIS) " << endl;
  fpara >> ifilename;
  flog1 << "Please enter the name of the file containing the mean monthly air temperature data: " << endl;
  flog1 << "               (e.g., TAIR.GIS) " << endl;
  flog1 << ifilename << endl << endl;
  ftair = fopen(ifilename, "r");

  if (!ftair)
  {
    cerr << "\nCannot open " << ifilename << " for data input" << endl;
    exit(-1);
  }

  telmnt[0].tem.atms.tprecflag = 0;
  if (equil == 0)
  {
    cout << "Do you have transient precipitation data?:" << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";
    fpara >> telmnt[0].tem.atms.tprecflag;

    flog1 << endl << "Do you have transient precipitation data?:" << endl;
    flog1 << "Enter 0 for No:" << endl;
    flog1 << "Enter 1 for Yes: " << endl;
    flog1 << "temnt[0].tem.atms.tprecflag = " << telmnt[0].tem.atms.tprecflag << endl << endl;
  }

  cout << "Please enter the name of the file containing the monthly precipitation data: " << endl;
  cout << "               (e.g., PREC.GIS) " << endl;
  fpara >> ifilename;
  flog1 << "Please enter the name of the file containing the monthly precipitation data: " << endl;
  flog1 << "               (e.g., PREC.GIS) " << endl;
  flog1 << ifilename << endl << endl;
  fprec = fopen(ifilename, "r");

  if (!fprec)
  {
    cerr << "\nCannot open " << ifilename << " for data input" << endl;
    exit(-1);
  }


  // added for hydrology model

  telmnt[0].tem.hyd.vapflag = 0;
  if (equil == 0) {
    cout << "Do you have transient vapor pressure data?:" << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";
    fpara >> telmnt[0].tem.hyd.vapflag;

    flog1 << endl << "Do you have transient vapor pressure data?:" << endl;
    flog1 << "Enter 0 for No:" << endl;
    flog1 << "Enter 1 for Yes: " << endl;
    flog1 << "telmnt[0].tem.hyd.vapflag = " << telmnt[0].tem.hyd.vapflag << endl << endl;
  }

  cout << "Please enter the name of the file containing the mean monthly vapor pressure data: " << endl;
  cout << "               (e.g., tvap.GIS) " << endl;
  fpara >> ifilename;
  flog1 << "Please enter the name of the file containing the mean monthly vapor pressure data: " << endl;
  flog1 << "               (e.g., TVAP.GIS) " << endl;
  flog1 << ifilename << endl << endl;
  fvap = fopen(ifilename, "r");

  if (!fvap) {
    cerr << "\nCannot open " << ifilename << " for data input" << endl;
    exit(-1);
  }

  // end of adding for hydrology model


  if (atmsflag == 0)
  {
    cout << "Please enter the name of the file containing the net irradiance data: " << endl;
    cout << "               (e.g., NIRR.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the net irradiance data: " << endl;
    flog1 << "               (e.g., NIRR.GIS) " << endl;
    flog1 <<  ifilename << endl << endl;
    fnirr = fopen(ifilename, "r");

    if (!fnirr)
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }

    cout << "Please enter the name of the file containing the mean monthly photosynthetically active radiation data: " << endl;
    cout << "               (e.g., PAR.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the mean monthly photosynthetically active radiation data: " << endl;
    flog1 << "               (e.g., PAR.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    fpar = fopen(ifilename, "r");

    if (!fpar)
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }
  }

  cout << endl << endl << "Enter the initial concentration of carbon dioxide in ppmv: ";
  fpara >> telmnt[0].tem.atms.initco2;
  flog1 << endl << endl << "Enter the initial concentration of carbon dioxide in ppmv: " << telmnt[0].tem.atms.initco2 << endl << endl;

  cout << endl << endl << "Enter the final equilibrium concentration of carbon dioxide in ppmv: ";
  fpara >> telmnt[0].tem.atms.co2level;
  flog1 << endl << endl << "Enter the final equilibrium concentration of carbon dioxide in ppmv: " << telmnt[0].tem.atms.co2level << endl << endl;

  telmnt[0].tem.atms.tco2flag = 0;
  if (equil == 0)
  {
    cout << "Do you have transient CO2 data?:" << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";
    fpara >> telmnt[0].tem.atms.tco2flag;

    flog1 << "Do you have transient CO2 data?:" << endl;
    flog1 << "Enter 0 for No:" << endl;
    flog1 << "Enter 1 for Yes: " << endl;
    flog1 << "telmnt[0].tem.atms.tco2flag = " << telmnt[0].tem.atms.tco2flag << endl << endl;

    if (telmnt[0].tem.atms.tco2flag == 1)
    {
      telmnt[0].tem.atms.loadmpitCO2(flog1, totsptime, RTIME);
    }
  }

  cout << endl << endl << "Enter the factor for changing C:N per ppmv of enhanced CO2:" << endl;
  cout << "                     (Enter 0.0 for no change): " << endl;
  fpara >> telmnt[0].tem.veg.dc2n[0];
  telmnt[0].tem.veg.dc2n[1] =telmnt[0].tem.veg.dc2n[0]; // added for plant function 2 which is moss layer

  flog1 << endl << "Enter the factor for changing C:N per ppmv of enhanced CO2:" << endl;
  flog1 << "                     (Enter 0.0 for no change): " << endl;
  flog1 << "telmnt[0].tem.veg.dc2n = " << telmnt[0].tem.veg.dc2n[0] << endl << endl;

  telmnt[0].tem.ag.tlulcflag = 0;
  if (equil == 0)
  {
    cout << "Do you have transient land use data?:" << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";
    fpara >> telmnt[0].tem.ag.tlulcflag;

    flog1 << "Do you have transient land use data?:" << endl;
    flog1 << "Enter 0 for No:" << endl;
    flog1 << "Enter 1 for Yes: ";
    flog1 << "telmnt[0].tem.ag.tlulcflag = " << telmnt[0].tem.ag.tlulcflag << endl << endl;;

    if (telmnt[0].tem.ag.tlulcflag == 1)
    {
      cout << "Please enter the name of the file containing the land use data: " << endl;
      cout << "               (e.g., LULC.GIS) " << endl;
      fpara >> ifilename;
      flog1 << "Please enter the name of the file containing the land use data: " << endl;
      flog1 << "               (e.g., LULC.GIS) " << endl;
      flog1 << ifilename << endl << endl;
      flulc = fopen(ifilename, "r");

      if (!flulc)
      {
	cerr << "\nCannot open " << ifilename << " for data input" << endl;
	exit(-1);
      }

      cout << "Will Relative Agricultural Production (RAP) always be 0.0?" << endl;
      cout << "Enter 0 for No:" << endl;
      cout << "Enter 1 for Yes: ";
      fpara >> telmnt[0].tem.ag.RAP0flag;

      flog1 << "Will Relative Agricultural Production (RAP) always be 0.0?" << endl;
      flog1 << "Enter 0 for No:" << endl;
      flog1 << "Enter 1 for Yes: " << endl;
      flog1 << "telmnt[0].tem.ag.RAP0flag = " << telmnt[0].tem.ag.RAP0flag << endl;

      if(telmnt[0].tem.ag.RAP0flag == 0)
      {
 	cout << "Please enter the name of the file containing the monthly potential NPP data: " << endl;
 	cout << "               (e.g., NPP.GIS) " << endl;
 	fpara >> ifilename;

 	flog1 << "Please enter the name of the file containing the monthly potential NPP data: " << endl;
 	flog1 << "               (e.g., NPP.GIS) " << endl;
	flog1 << ifilename << endl << endl;
	fnpp = fopen(ifilename, "r");

	if (!fnpp)
        {
	  cerr << "\nCannot open " << ifilename << " for data input" << endl;
	  exit(-1);
	}
      }

      telmnt[0].tem.ag.getecd(flog1);
    }
  }

  kdinflg = 0;
  cout << "Do you want to use spatially explicit values of Kd?" << endl;
  cout << "Enter 0 for No:" << endl;
  cout << "Enter 1 for Yes: ";
  fpara >> kdinflg;

  flog1 << "Do you want to use spatially explicit values of Kd?" << endl;
  flog1 << "Enter 0 for No:" << endl;
  flog1 << "Enter 1 for Yes: " << endl;
  flog1 << "kdinflg = " << kdinflg << endl << endl;

  if (kdinflg == 1)
  {
    cout << "Please enter the name of the file containing the Kd values: " << endl;
    cout << "               (e.g., KDIN.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the Kd values: " << endl;
    flog1 << "               (e.g., KDIN.GIS) " << endl;
    flog1 << ifilename << endl << endl;

    fkdin = fopen(ifilename, "r");

    if (!fkdin)
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }
  }

  kdoutflg = 0;
  cout << "Do you want to generate spatially explicit values of Kd?" << endl;
  cout << "Enter 0 for No:" << endl;
  cout << "Enter 1 for Yes: ";
  fpara >> kdoutflg;

  flog1 << endl << "Do you want to generate spatially explicit values of Kd?" << endl;
  flog1 << "Enter 0 for No:" << endl;
  flog1 << "Enter 1 for Yes: " << endl;
  flog1 << "kdoutflg = " << kdoutflg << endl << endl;

  if (kdoutflg == 1)
  {
    cout << "Please enter the name of the file to contain the Kd values: " << endl;
    cout << "               (e.g., KDOUT.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file to contain the Kd values: " << endl;
    flog1 << "               (e.g., KDOUT.GIS) " << endl;
    flog1 << ifilename << endl << endl;

    fkdout.open(ifilename, ios::app);
  }

  ntempred = askpred(telmnt[0].tem.predstr,totpred);
  totpred += ntempred;

  for (i = natmspred; i < totpred; i++)
  {
    cout << endl << "Enter the name of the OUTPUT file to contain " << predmap[i] << ":  ";
    fpara >> tempredfile;
    flog1 << "Enter the name of the OUTPUT file to contain " << predmap[i] << ":  " << tempredfile << endl;
    ftempred[i].open(tempredfile, ios::app);
  }

  stateflag = 0;
  cout << endl << "Do you want to use spatially explicit data for initial conditions? " << endl;
  cout << "  Enter 1 for YES:" << endl;
  cout << "  Enter 0 for NO:" << endl;
  cout << "  NOTE:  If YES, you will need spatially explicit data for VEGC," << endl;
  cout << "         STRUCTN, SOLC, SOLN, AVLN, NSTORE, AVAILH2O, RGRNDH2O," << endl;
  cout << "         SNOWPACK, SGRNDH2O, UNNORMLEAF, PET, EET,(i.e. 13 data files)" << endl;
  fpara >> stateflag;

  flog1 << endl << "Do you want to use spatially explicit data for intial conditions? " << endl;
  flog1 << "  Enter 1 for YES:" << endl;
  flog1 << "  Enter 0 for NO:" << endl;
  flog1 << "  NOTE:  If YES, you will need spatially explicit data for VEGC," << endl;
  flog1 << "         STRUCTN, SOLC, SOLN, AVLN, NSTORE, AVAILH2O, RGRNDH2O," << endl;
  flog1 << "         SNOWPACK, SGRNDH2O, UNNORMLEAF, PET, EET,(i.e. 13 data files)" << endl;
  flog1 << "stateflag = " << stateflag << endl << endl;

  if (stateflag == 1)
  {
    cout << "Please enter the name of the file containing the vegetation carbon data: " << endl;
    cout << "               (e.g., VEGC.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the vegetation carbon data: " << endl;
    flog1 << "               (e.g., VEGC.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    fstate[0] = fopen(ifilename, "r");

    if (!fstate[0])
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }

    cout << "Please enter the name of the file containing the vegetation structural nitrogen data: " << endl;
    cout << "               (e.g., STRN.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the vegetation structural nitrogen data: " << endl;
    flog1 << "               (e.g., STRN.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    fstate[1] = fopen(ifilename, "r");

    if (!fstate[1])
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }

    cout << "Please enter the name of the file containing the soil organic carbon data: " << endl;
    cout << "               (e.g., SOLC.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the soil organic carbon data: " << endl;
    flog1 << "               (e.g., SOLC.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    fstate[2] = fopen(ifilename, "r");

    if (!fstate[2])
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }

    cout << "Please enter the name of the file containing the soil organic nitrogen data: " << endl;
    cout << "               (e.g., SOLN.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the soil organic nitrogen data: " << endl;
    flog1 << "               (e.g., SOLN.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    fstate[3] = fopen(ifilename, "r");

    if (!fstate[3])
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }

    cout << "Please enter the name of the file containing the soil available nitrogen data: " << endl;
    cout << "               (e.g., AVLN.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the soil availablenitrogen data: " << endl;
    flog1 << "               (e.g., AVLN.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    fstate[4] = fopen(ifilename, "r");

    if (!fstate[4])
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }

    cout << "Please enter the name of the file containing the vegetation labile nitrogen data: " << endl;
    cout << "               (e.g., STON.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the vegetation labile nitrogen data: " << endl;
    flog1 << "               (e.g., STON.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    fstate[5] = fopen(ifilename, "r");

    if (!fstate[5])
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }

    cout << "Please enter the name of the file containing the available soil moisture data: " << endl;
    cout << "               (e.g., AVLW.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the available soil moisture data: " << endl;
    flog1 << "               (e.g., AVLW.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    fstate[6] = fopen(ifilename, "r");

    if (!fstate[6])
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }

    cout << "Please enter the name of the file containing the rain-derived groundwater data: " << endl;
    cout << "               (e.g., RGRW.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the rain-derived groundwater data: " << endl;
    flog1 << "               (e.g., RGRW.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    fstate[7] = fopen(ifilename, "r");

    if (!fstate[7])
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }

    cout << "Please enter the name of the file containing the snowpack data: " << endl;
    cout << "               (e.g., SNWP.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the snowpack data: " << endl;
    flog1 << "               (e.g., SNWP.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    fstate[8] = fopen(ifilename, "r");

    if (!fstate[8])
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }

    cout << "Please enter the name of the file containing the snowpack-derived groundwater data: " << endl;
    cout << "               (e.g., SGRW.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the snowpack-derived groundwater data: " << endl;
    flog1 << "               (e.g., SGRW.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    fstate[9] = fopen(ifilename, "r");

    if (!fstate[9])
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }

    cout << "Please enter the name of the file containing the unnormalized leaf phenology data: " << endl;
    cout << "               (e.g., UNLF.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the unnormalized leaf phenology data: " << endl;
    flog1 << "               (e.g., UNLF.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    fstate[10] = fopen(ifilename, "r");

    if (!fstate[10])
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }

    cout << "Please enter the name of the file containing the potential evapotranspiration data: " << endl;
    cout << "               (e.g., PET.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the potential evapotranspiration data: " << endl;
    flog1 << "               (e.g., PET.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    fstate[11] = fopen(ifilename, "r");

    if (!fstate[11])
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }

    cout << "Please enter the name of the file containing the estimated evapotranspiration data: " << endl;
    cout << "               (e.g., EET.GIS) " << endl;
    fpara >> ifilename;
    flog1 << "Please enter the name of the file containing the estimated evapotranspiration data: " << endl;
    flog1 << "               (e.g., EET.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    fstate[12] = fopen(ifilename, "r");

    if (!fstate[12])
    {
      cerr << "\nCannot open " << ifilename << " for data input" << endl;
      exit(-1);
    }
  }

};
