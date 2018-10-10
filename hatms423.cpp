/* **************************************************************
*****************************************************************
HATMS423.CPP - uses the physical characteristics of the atmosphere
               as described by TEM, but uses transient CO2, specified
               every 6 months, from an ASCII file
*****************************************************************
************************************************************** */

#if !defined(HATMS423_H)
  #include "hatms423.hpp"
#endif

MPIatms::MPIatms() : Atmosphere() { };

/* *************************************************************
************************************************************** */

void MPIatms::loadfosfuel(ofstream& rflog1, const int& RTIME)
{

  int iyear;
  char ifilename[40];
  int testyear;
  FILE* ffossil;

  cout << endl << endl << "Enter the name of the file containing the fossil fuel data: ";
  fpara >> ifilename;
  rflog1 << endl << endl << "Enter the name of the file containing the fossil fuel data: ";
  rflog1 << ifilename << endl;
  ffossil = fopen(ifilename, "r");

  if (!ffossil)
  {
    cerr << "\nCannot open " << ifilename << " for data input" << endl;
    exit(-1);
  }

  cout << "What year should the data be normalized against?" << endl;
  fpara >> testyear;


// Initialize historical fossil fuel emissions

  for (iyear = 0; iyear < RTIME; iyear++)
  {
    fscanf(ffossil, "%d %lf", ffuelyear+iyear, ffuel+iyear);
    if (ffuelyear[iyear] == testyear) { maxffuel = ffuel[iyear]; }
  }

  fclose(ffossil);

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void MPIatms::loadmpitCO2(ofstream& rflog1, const int& totsptime, const int& RTIME)
{

  int i;
  int j;
  int dm;
  char ifilename[40];
  CO2data co2dat;
  ifstream fco2;
  double mco2[MAXRTIME];         // CO2 concentration in July of year

  cout << endl << endl << "Enter the name of the file containing the CO2 data: ";
  fpara >> ifilename;
  rflog1 << endl << endl << "Enter the name of the file containing the CO2 data: ";
  rflog1 << ifilename << endl;
  fco2.open(ifilename, ios::in);

  if (!fco2)
  {
    cerr << "\nCannot open " << ifilename << " for data input" << endl;
    exit(-1);
  }

// Initialize transient CO2 concentrations

  for (i = (totsptime+1); i < (RTIME+1); i++)
  {
    j = i - (totsptime + 1);
    co2dat.get(fco2);
    co2year[i] = (long) co2dat.year;
    mco2[j] = co2dat.mco2;
  }
  for (i = (totsptime+1); i < RTIME; i++)
  {
    j = i - (totsptime + 1);
    for (dm = 0; dm < CYCLE; dm++)
    {
      if (j == 0) { tco2[i][dm] = mco2[j]; }
      else
      {
        if (dm < 6)
        {
          tco2[i][dm] = mco2[j-1] + ((dm + 6) * (mco2[j] - mco2[j-1]) / (double) CYCLE);
        }
        else
        {
          tco2[i][dm] = mco2[j] + ((dm - 6) * (mco2[j+1] - mco2[j]) / (double) CYCLE);
        }
      }
    }
  }

  fco2.close();

  co2year[0] = co2year[totsptime+1] - totsptime - 1;
  for (dm = 0; dm < CYCLE; dm++) { tco2[0][dm] = co2level; }
  for (i = 1; i < (totsptime+1); i++)
  {
    co2year[i] = co2year[0] + i;
    for (dm = 0; dm < CYCLE; dm++) { tco2[i][dm] = mco2[0]; }
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void MPIatms::loadmpitCO2(char ifilename[80], const int& totsptime, const int& RTIME)
{

  int i;
  int j;
  int dm;
  CO2data co2dat;
  ifstream fco2;
  double mco2[MAXRTIME];         // CO2 concentration in July of year

  fco2.open(ifilename, ios::in);

  if (!fco2)
  {
    cerr << "\nCannot open " << ifilename << " for data input" << endl;
    exit(-1);
  }

// Initialize transient CO2 concentrations

  for (i = (totsptime+1); i < (RTIME+1); i++)
  {
    j = i - (totsptime + 1);
    co2dat.get(fco2);
    co2year[i] = (long) co2dat.year;
    mco2[j] = co2dat.mco2;
  }
  for (i = (totsptime+1); i < RTIME; i++)
  {
    j = i - (totsptime + 1);
    for (dm = 0; dm < CYCLE; dm++)
    {
      if (j == 0) { tco2[i][dm] = mco2[j]; }
      else
      {
        if (dm < 6)
        {
          tco2[i][dm] = mco2[j-1] + ((dm + 6) * (mco2[j] - mco2[j-1]) / (double) CYCLE);
        }
        else
        {
          tco2[i][dm] = mco2[j] + ((dm - 6) * (mco2[j+1] - mco2[j]) / (double) CYCLE);
        }
      }
    }
  }

  fco2.close();

  co2year[0] = co2year[totsptime+1] - totsptime - 1;
  for (dm = 0; dm < CYCLE; dm++) { tco2[0][dm] = co2level; }
  for (i = 1; i < (totsptime+1); i++)
  {
    co2year[i] = co2year[0] + i;
    for (dm = 0; dm < CYCLE; dm++) { tco2[i][dm] = mco2[0]; }
  }

};


double MPIatms::mkclds(const double& girr, const double& nirr)
{

  double clouds;

  if (nirr >= (0.71 * girr)) { return clouds = 0.0; }
  else
  {
    clouds = 1.0 - (((nirr/girr) - 0.23)/0.48);
    clouds *= 100.0;
  }
  if (clouds > 100.0) { clouds = 100.0; }

  return clouds;

};
/*

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

double MPIatms::xgirr(double& lat, int& dm, double sumday)
{

  const double pi = 3.141592654;                // Greek "pi"
  const double sp = 1368.0 * 3600.0 / 41860.0;  // solar constant

  double lambda;
  double sumd;
  double sig;
  double eta;
  double sinbeta;
  double sb;
  double sotd;
  int day;
  int hour;
  double gross;

  ndays[0] = 31;
  ndays[1] = 28;
  ndays[2] = 31;
  ndays[3] = 30;
  ndays[4] = 31;
  ndays[5] = 30;
  ndays[6] = 31;
  ndays[7] = 31;
  ndays[8] = 30;
  ndays[9] = 31;
  ndays[10] = 30;
  ndays[11] = 31;

  lambda = lat * pi / 180.0;

  gross = 0.0;
  for (day = 0; day < ndays[dm]; day++)
  {
    ++sumday;
    sumd = 0;
    sig = -23.4856*cos(2 * pi * (sumday + 10.0)/365.25);
    sig *= pi / 180.0;

    for (hour = 0; hour < 24; hour++)
    {
      eta = (double) ((hour+1) - 12) * pi / 12.0;
      sinbeta = sin(lambda)*sin(sig) + cos(lambda)*cos(sig)*cos(eta);
      sotd = 1 - (0.016729 * cos(0.9856 * (sumday - 4.0) * pi / 180.0));
      sb = sp * sinbeta / pow(sotd,2.0);
      if (sb >= 0.0) { sumd += sb; }
    }

    gross += sumd;
  }

  gross /= (double) ndays[dm];

  return gross;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double MPIatms::xnirr(const double& clds, const double& girr)
{

  double nirr;

  if (clds >= 0.0)
  {
    nirr = girr * (0.251 + (0.509*(1.0 - clds/100.0)));
  }
  else { nirr = MISSING; }
  
  return nirr;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double MPIatms::xpar(const double& clds, const double& nirr)
{

  double par;

  if (clds >= 0.0)
  {
      par = nirr * ((0.2 * clds / 100.0) + 0.45);
  }
  else { par = MISSING; }

  return par;

};


