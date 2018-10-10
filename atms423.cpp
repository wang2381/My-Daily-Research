/* ****************************************************************
ATMS423.CPP - object describes physical characteristics of the
	      atmosphere

           - modified from atms41d.cpp by DWK 20000102
*****************************************************************
************************************************************** */

#if !defined(ATMS423_H)
 #include "atms423.hpp"
#endif

Atmosphere::Atmosphere()
{

  co2level = 0;

// Initialize number of days per month for each month of year

  daze[0] = daze[2] = daze[4] = daze[6] = daze[7] = 31.0;
  daze[9] = daze[11] = 31.0;
  daze[3] = daze[5] = daze[8] = daze[10] = 30.0;
  daze[1] = 28.0;


};

/* **************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

double Atmosphere::petjh(const double& nirr, const double& tair, const int& dm)
{

  double f;
  double rt;
  double pet;

  f = ((9.0/5.0) * tair) + 32.0;
  rt = nirr * 0.016742;
  pet = ((0.014*f) - 0.37) * rt * daze[dm];

  if (pet < 0.0) { pet = 0.0; }

  return pet;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Atmosphere::precsplt(double& prec, double& tair, double& rain, double& snowfall)
{


/* *************************************************************
	Willmott's assumptions on snow/rain split:
************************************************************** */

  if (tair >= 0.0)
  {
    rain = prec;
    snowfall = 0.0;
  }
  else
  {
    rain = 0.0;
    snowfall = prec;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

double Atmosphere::xeet(const double& rain, const double& snowinf, const double& pet, const double& avlh2o, const double& awcapmm, const int& dm)
{

  const double edpar = 5.0;
  double gm;
  double ep;
  double def;
  double prob;
  double rbar;
  double dsm;
  double aet;

  if ((rain+snowinf) >= pet)
  {
    aet = pet;
  }
  else
  {
    gm = (1.0 - exp(-edpar * avlh2o/awcapmm)) / (1.0 - exp(-edpar));
    ep = pet / daze[dm];
    def = ep + awcapmm - avlh2o;
    prob = 1.0 - exp(-0.005*(rain + snowinf));
    if (prob != 0.0) { rbar = (rain + snowinf) / (daze[dm] * prob); }
    else { rbar = 0.0; }

    if (rbar != 0.0)
    {
      dsm = rbar*prob*(gm + ((1.0-gm) * exp(-ep/rbar)) - exp(-def/rbar)) - (ep*gm);
    }
    else {
      dsm = -ep*gm;
    }

    dsm *= daze[dm];

    aet = rain + snowinf - dsm;
    if (aet > pet) { aet = pet; }

  }

  return aet;
};


double Atmosphere::xgirr(double& lat, int& dm, double& sumday)
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
  int ndays[12];
  
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

double Atmosphere::xnirr(const double& clds, const double& girr)
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

double Atmosphere::xpar(const double& clds, const double& nirr)
{

  double par;

  if (clds >= 0.0)
  {
      par = nirr * ((0.2 * clds / 100.0) + 0.45);
  }
  else { par = MISSING; }

  return par;

};
