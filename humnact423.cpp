/* **************************************************************
*****************************************************************
HUMANACT42.CPP - describes human disturbances to natural ecosystems

Modifications:

19990821 - DWK changed char dummy[12] to char dummy[30] in void getecd(ofstream& rflog1)
19990821 - DWK changed char dummy[12] to char dummy[30] in void getecd(char ecd[80])
20000102 - DWK added compiler directives

*****************************************************************
************************************************************** */

#if !defined(HUMNACT423_H)
  #include "humnact423.hpp"
#endif

Humanact::Humanact()
{

  int vt = 0; // local use for moss component
  int i;

   c2n[vt] = 54.29;
//  cfall = 0.20;
//  nfall = 0.20;

  state = 0;
  prvstate = 0;

  for (i = 0; i < 10; i++)
  {
    initPROD10[vt][i].carbon = 0.0;
    initPROD10[vt][i].nitrogen = 0.0;
  }
  for (i = 0; i < 100; i++)
  {
    initPROD100[vt][i].carbon = 0.0;
    initPROD100[vt][i].nitrogen = 0.0;
  }

};

/* **************************************************************
************************************************************** */



/* **************************************************************
************************************************************** */

void Humanact::conversion(int& ez, Tveg42& veg, Tsoil4& soil, const int& vt)
{

  int dm;

  formPROD10[vt].carbon  = prod10par[ez] * veg.plant[vt][CYCLE-1].carbon;
  formPROD10[vt].nitrogen  = prod10par[ez] * (veg.strctrl[vt][CYCLE-1].nitrogen + veg.labile[vt][CYCLE-1].nitrogen);
  formPROD100[vt].carbon = prod100par[ez] * veg.plant[vt][CYCLE-1].carbon;
  formPROD100[vt].nitrogen = prod100par[ez] * (veg.strctrl[vt][CYCLE-1].nitrogen + veg.labile[vt][CYCLE-1].nitrogen);

  PROD10[vt].carbon += formPROD10[vt].carbon;
  PROD10[vt].nitrogen += formPROD10[vt].nitrogen;
  PROD100[vt].carbon += formPROD100[vt].carbon;
  PROD100[vt].nitrogen += formPROD100[vt].nitrogen;
  for (dm = 0; dm < CYCLE; dm++)
  {
    slash[vt][dm].carbon = slashpar[ez] * veg.plant[vt][CYCLE-1].carbon / (double) CYCLE;
    slash[vt][dm].nitrogen = slashpar[ez]  * (veg.strctrl[vt][CYCLE-1].nitrogen + veg.labile[vt][CYCLE-1].nitrogen) / (double) CYCLE;
    vconvrtflx[vt][dm].carbon = (vconvert[ez] * veg.plant[vt][CYCLE-1].carbon) / (double) CYCLE;
    sconvrtflx[vt][dm].carbon =  (sconvert[ez] * soil.org[vt][CYCLE-1].carbon)/ (double) CYCLE;
    convrtflx[vt][dm].carbon = vconvrtflx[vt][dm].carbon + sconvrtflx[vt][dm].carbon;
    vconvrtflx[vt][dm].nitrogen = ((1.0 - nvretconv[ez]) * vconvert[ez] * (veg.strctrl[vt][CYCLE-1].nitrogen + veg.labile[vt][CYCLE-1].nitrogen))  / (double) CYCLE;
    sconvrtflx[vt][dm].nitrogen = ((1.0 - nsretconv[ez]) * sconvert[ez] * soil.org[vt][CYCLE-1].nitrogen)  / (double) CYCLE;
    convrtflx[vt][dm].nitrogen = vconvrtflx[vt][dm].nitrogen + sconvrtflx[vt][dm].nitrogen;
    nvretent[vt][dm] = (nvretconv[ez] * vconvert[ez] * (veg.strctrl[vt][CYCLE-1].nitrogen + veg.labile[vt][CYCLE-1].nitrogen)) / (double) CYCLE;
    nsretent[vt][dm] = (nsretconv[ez] * sconvert[ez] * soil.org[vt][CYCLE-1].nitrogen) / (double) CYCLE;
    nretent[vt][dm] = nvretent[vt][dm] + nsretent[vt][dm];

    yrvconvrtC[vt] += vconvrtflx[vt][dm].carbon;
    yrvconvrtN[vt] += vconvrtflx[vt][dm].nitrogen;
    yrnrent[vt] += nretent[vt][dm];
    yrnvrent[vt] += nvretent[vt][dm];

  }

};
 /* **************************************************************
************************************************************** */

void Humanact::getecd(ofstream& rflog1)
{

  char ecd[80];

  cout << "Enter name of the data file (.ECD) with agricultural parameter values: " << endl;
  fpara >> ecd;

  rflog1 << "Enter name of the soil data file (.ECD) with agricultural parameter values: " << ecd << endl;

  getecd(ecd);

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void Humanact::getecd(char ecd[80])
{

  int NUMVAR = 13;
  int i;
  int dcmnt;
  char dummy[30];
  ifstream infile;

  long update;

  infile.open(ecd, ios::in);

  if (!infile)
  {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy; }
  for (dcmnt = 1; dcmnt < MAXCMNT; dcmnt++)
  {
    infile >> dummy >> dummy;
    infile >> slashpar[dcmnt] >> vconvert[dcmnt] >> prod10par[dcmnt] >> prod100par[dcmnt];
    infile >> sconvert[dcmnt] >> nvretconv[dcmnt] >> nsretconv[dcmnt] >> vrespar[dcmnt];
    infile >> cfall[dcmnt] >> nfall[dcmnt];
    infile >> update;
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void Humanact::resetPROD(const int& vt)
{

  PROD1[vt].carbon = 0.0;
  PROD1[vt].nitrogen = 0.0;
  PROD1decay[vt].carbon = 0.0;
  PROD1decay[vt].nitrogen = 0.0;

  formPROD10[vt].carbon  = 0.0;
  formPROD10[vt].nitrogen  = 0.0;
  PROD10decay[vt].carbon  = 0.0;
  PROD10decay[vt].nitrogen  = 0.0;
  PROD10[vt].carbon = 0.0;
  PROD10[vt].nitrogen = 0.0;

  formPROD100[vt].carbon = 0.0;
  formPROD100[vt].nitrogen = 0.0;
  PROD100decay[vt].carbon = 0.0;
  PROD100decay[vt].nitrogen = 0.0;
  PROD100[vt].carbon = 0.0;
  PROD100[vt].nitrogen = 0.0;

  TOTPROD[vt].carbon = 0.0;
  TOTPROD[vt].nitrogen = 0.0;

  formTOTPROD[vt].carbon = 0.0;
  formTOTPROD[vt].nitrogen = 0.0;

  TOTPRODdecay[vt].carbon = 0.0;
  TOTPRODdecay[vt].nitrogen = 0.0;

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void Humanact::updateyr(const int& dyr, const int& vt)
{

  int i;
  int j;
  int k;
  double yrtrash[NVT];

// for (vt =0; vt < NVT; vt++) {
  prvstate = state;
  PROD1[vt].carbon = formPROD1[vt].carbon;
  PROD1[vt].nitrogen = formPROD1[vt].nitrogen;
  PROD1decay[vt].carbon = PROD1[vt].carbon;
  PROD1decay[vt].nitrogen = PROD1[vt].nitrogen;

  j = dyr%10;
  initPROD10[vt][j].carbon = formPROD10[vt].carbon;
  initPROD10[vt][j].nitrogen = formPROD10[vt].nitrogen;
  PROD10decay[vt].carbon  = 0.0;
  PROD10decay[vt].nitrogen  = 0.0;
  for ( i = 0; i < 10; i++)
  {
    yrtrash[vt] = initPROD10[vt][i].carbon * 0.10;
    PROD10decay[vt].carbon += yrtrash[vt];
    yrtrash[vt] = initPROD10[vt][i].nitrogen * 0.10;
    PROD10decay[vt].nitrogen += yrtrash[vt];
  }
  PROD10[vt].carbon -= PROD10decay[vt].carbon;
  PROD10[vt].nitrogen -= PROD10decay[vt].nitrogen;

  k = dyr%100;
  initPROD100[vt][k].carbon = formPROD100[vt].carbon;
  initPROD100[vt][k].nitrogen = formPROD100[vt].nitrogen;
  PROD100decay[vt].carbon = 0.0;
  PROD100decay[vt].nitrogen = 0.0;
  for ( i = 0; i < 100; i++)
  {
    yrtrash[vt] = initPROD100[vt][i].carbon * 0.01;
    PROD100decay[vt].carbon += yrtrash[vt];
    yrtrash[vt] = initPROD100[vt][i].nitrogen * 0.01;
    PROD100decay[vt].nitrogen += yrtrash[vt];
  }

  PROD100[vt].carbon -= PROD100decay[vt].carbon;
  PROD100[vt].nitrogen -= PROD100decay[vt].nitrogen;

  TOTPROD[vt].carbon = PROD1[vt].carbon + PROD10[vt].carbon + PROD100[vt].carbon;
  TOTPROD[vt].nitrogen = PROD1[vt].nitrogen + PROD10[vt].nitrogen + PROD100[vt].nitrogen;

  formTOTPROD[vt].carbon = formPROD1[vt].carbon + formPROD10[vt].carbon + formPROD100[vt].carbon;
  formTOTPROD[vt].nitrogen = formPROD1[vt].nitrogen + formPROD10[vt].nitrogen + formPROD100[vt].nitrogen;

  TOTPRODdecay[vt].carbon = PROD1decay[vt].carbon + PROD10decay[vt].carbon + PROD100decay[vt].carbon;
  TOTPRODdecay[vt].nitrogen = PROD1decay[vt].nitrogen + PROD10decay[vt].nitrogen + PROD100decay[vt].nitrogen;
//} // for vt
};
