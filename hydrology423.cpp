
/* ****************************************************************
HYDROLOGY.CPP - object describes physical characteristics of the
	      atmosphere

           - Created by QZ for hydrology model 06162000
*****************************************************************
************************************************************** */

#if !defined(HYDROLOFY423_H)
  #include "hydrology423.hpp"
#endif

Hydrology::Hydrology() {


// Initialize number of days per month for each month of year

  daze[0] = daze[2] = daze[4] = daze[6] = daze[7] = 31.0;
  daze[9] = daze[11] = 31.0;
  daze[3] = daze[5] = daze[8] = daze[10] = 30.0;
  daze[1] = 28.0;


};

void Hydrology::drad(const double& LAI, const double& EXT, const double& nirr, const int& dm, const int& vt)

{

// double nirr[CYCLE];        // Net Irradiance   (cal/(sq. cm * day))
  // drad (MJ m-2 day-1) canopy daily average radiation
  // nirr (MJ m-2 day-1) incoming short-wave readiation
  // LAI (m2 m-2)  leaf area index
  // EXT dimentionless  extinction coefficient of radiation through the canopy
  // 2.2 is the coefficiet changig total LAI to projected LAI (dimentionless)
  double drad;
  double Lsunlit; // the sunlit leaf area index
  double nirrt; // temporary use

   // drad = (nirr * (1-exp(-(LAI/2.2) * EXT))) / ((-1)* EXT * LAI / 2.2);
    drad = (nirr * (1-exp(-(LAI/2.2) * EXT))) / (EXT * LAI / 2.2); // Running and Coughlan 1988.

    nirrt = drad / 0.2388; // J /cm2 /d
    nirrt = nirrt /pow(10,6); //MJ cm-2 d-1
    nirrt = nirrt * pow(10,4); // MJ M-2 D-1
    drad = nirrt;
    hdrad[vt][dm] = drad;

    // solar radiation to the ground (or soil surface layer)
    surfacerad[vt][dm] = drad * exp(-EXT * LAI /2.2); // page 124, Coughlan and Running 1997, landscape ecology, beer's law
  }



double Hydrology::LWP(double& SOILH2O1, double& SOILCAP, const double sand, const double clay, const double silt)

{
   // LWP daily maimum leaf water potential (Mpa)
   // SOILH2O soil water content at time (t) m3
   // SOILCAP soil water capacity, m3
   double LWP;
   double temp_water;
   double soilcapacity;
//   const watercap = 36.0; // m3/m3 Frolking et al, 1996 for NSA-OBS
     																						//wsr
     temp_water = SOILH2O1;
     if (temp_water ==0) temp_water = 0.01; // protected

     LWP = 0.2 /(temp_water / SOILCAP); // Running and Coughlan 1988

     
/* New solution, see:
	Cosby, B.J., G.M. Hornberger, R.B. Clapp, and T.R. Ginn, 1984.  A
	   statistical exploration of the relationships of soil moisture
	   characteristics to the physical properties of soils.  Water Res.
	   Res. 20:682-690.

	Saxton, K.E., W.J. Rawls, J.S. Romberger, and R.I. Papendick, 1986.
		Estimating generalized soil-water characteristics from texture.
		Soil Sci. Soc. Am. J. 50:1031-1036.

   double soil_b = -(3.10 + 0.157*clay - 0.003*sand);
	double vwc_sat = (50.5 - 0.142*sand - 0.037*clay)/100.0;
	double psi_sat = -(exp((1.54 - 0.0095*sand /100 + 0.0063*silt /100)*log(10.0))*9.8e-5);

//	LWP = psi_sat * powl((SOILH2O / vwc_sat), soil_b);
	LWP = psi_sat * powl((SOILH2O /vwc_sat) , soil_b);
   printf(" %5.3f  %5.3f  %5.3f ", soil_b,vwc_sat,psi_sat );

   we may also look at Bachelet et al., 1998, sensitivity of a biogeography model to soil properties about the soil water potential
   */
   return LWP;

   }


double Hydrology::CanopyCond(const double& LWPmin, const double& LWP,const double& CCmax,

const double& LWPc, const double& ABSHD, const double& SLOPEcc)

{
  // CCw canopy H2o conductance (ms-1)
  // CCmax maximum canopy conductance (ms-1)
  // DCCw slope of CCmax vs. LWP (ms-1 Mpa-1) = CCmax / (LWP at stomatal closure - LWP minimum)
  // LWP daily maximum leaf water potential from function of LWP above (MPa)
  // LWPmin minimum leaf water potential (MPa)
  // LWPc LWP at stomatal closure (MPa)
  // CCh canopy conductance with humidity reduction (ms-1)
  // DCCh slope of CC vs. ABSHD (ms-1 g-1m-3) =  CCmax / (SLOPEcc - ABSHD)
  // SLOPEcc slope cc humidity reduction (ms-1/ g-1m-3)
  // ABSHD abolute humidity deficit (gm-3)

  double CCw;
  double DCCw;
  double DCCh;
  double CCh;

  DCCw = CCmax /( LWPc - LWPmin);
  CCw = CCmax - DCCw * (LWP -LWPmin);
  DCCh = CCmax / (SLOPEcc - ABSHD);
  CCh = CCw * (1 - DCCh * ABSHD);

  return CCh;

 }


double Hydrology::slope(const double& airt)
{
 // airt daily or monthly air temperature 0C
 // dewpoint daily or monthly dewpoint or minimum air temperature at night oC

 // calculate temperature offsets for slope estimate
  const double dt = 0.2; // offset to calculate the slope
  double t1, t2;
  double svp1,svp2;
  double slope;

    t1 = airt+dt;
    t2 = airt-dt;
 // calculate saturation vapor pressures at t1 and t2

   svp1 = 6.1078 * exp(( 17.269  * t1) / (237.3 + t1));
   svp2 = 6.1078 * exp(( 17.269  * t2) / (237.3 + t2));
 // calculate slope of pvs vs. T curve, at ta
   slope = (svp1-svp2) / (t1-t2);
 return slope;
}

// to detemine the vapor pressure deficit
void Hydrology::vpd(const double& vaporP, const double& airt, const int& vt)
{

const double Mw = 18; // g/mol
const double R = 8.3143; // J/mol/k

double eadt;
double ea;
double dewpointT; // oC

// vaporP ambient vapor pressure Kpa
// airt oC
// absoluteHD g/m3
// detemine depoint temperature first. see Campbell and Norman, Td = c * ln(ea/ a) /(b - ln(ea/a))
// where c= 240.97 OC, a = 0.611Kpa, b =17.502
dewpointT = (240.97 * logl(vaporP/0.611)) / (17.502 - logl(vaporP/0.611));

eadt = 0.611 * exp((17.502 * airt) / (airt + 240.97));
vpdeficit[vt] =  fabs(eadt - vaporP);
vpdeficit[vt] = (vpdeficit[vt] * 10); // trasfer to mbar
absoluteHD[vt] = fabs((Mw / R) * ( (eadt * 1000 / (273.2 +dewpointT)) - (vaporP * 1000 /  (273.2 + airt))));
}

//calculate density of air (rho) as a function of air temperature
double Hydrology::densityA(const double& airt)

{
// pa   (kg/m3)   density of air
// airt air temperature oC
  double pa;
  pa = 1.292 - (0.00428 * airt);
  return pa;
}


// calculate latent heat of vaporization as a function of ta
double Hydrology::latentV(const double& airt)
{
  // le  (J/kg)        latent heat of vaporization of water
  //airt air temperature oC
   double le;
   le = 2.5023e6 - 2430.54 * airt;
   return le;
}


double Hydrology::penmanmonteith(const double& slope, double drad, const double& pa,
const double& vpd, const double& cch, const double& le)

{

// slope is from slope fuction (dimensionless)
// drad from function drad (MJ m-2 daily)
//CP sepcific heat of air (J Kg -1 oC-1) constant =
 static double CP = 1010.0;   ///(J/kg K) = (J Kg-1 oC-1) specific heat of air, relationship: oC = K - 273.15
// pa density of air (kgm-3)  from function of DensityA
// vpd vapor pressure deficit from canopy to air (mbar)
//ra canopy aerodynamic resistance ( fixed at 5.0 sm-1) typical value 5-10 for forests
 static double ra = 2.0;
// rc canopy resitance to water vapor = 1 / CCh from CanopyCond function (sm-1)
// gamma psychometric constant (mbar oC-1)
 static double gamma = 0.662; // (mbar k-1) = mbar oC-1
// le latent heat of vaporization of water (J kg-1)  from function of LatentV
// lai leaf area index (m2m-2)
//dayl daylength (s day-1)
 static long dayl = 86400; // 24*3600  (second day-1)

 double rc;
 double trans; //  canopy transpiration daily or monthly (m3 day-1)
  rc = 1 / cch;
  drad = drad * 10e6 / (24 * 3600); // convert MJ/m2/day to J/m2/s
  trans = ((((slope * drad) + (CP * pa) * vpd / ra) / (slope + gamma * (1 + rc /ra ))) / (le * 1000)) * dayl;
  if (trans > 0.0)
  {  trans = trans * 1000.0; // transfer the m3/month to mm/month/m2
  }
  else trans = 0.0001;
  return trans;

 }

// evaporation from canopy
double Hydrology::evaporation(const double& water_on_canopy)
{
   // prec_to_canopy; water on canopy will be evaporated at all (mm)
double evaporation;

	if (water_on_canopy > 0.0) evaporation = water_on_canopy;
   else evaporation =0.0;

    return evaporation;
}


// added by Qianlai zhuang for hydrological version  14 June 2000
// after interception, the rest of rainfall goes to soil, and assigned to varible RAIN of WBM, and interception goes to canopy for evaporation modelling
// interception rate is obtained from
//  Helvey (1971) and helvey and Patric (1965)

void Hydrology::intercept(const double& rain, const double& inter_coef, const int& dm, const int& vt)
{

   double max_int; // maximum daily canopy interception mm intercepted/ mm rain in the canopy */
   double tempuse;
//   double  prec_to_canopy;

// 	max_int = inter_coef * rain;
      tempuse = rain /10.0 ; //convert to cm
      max_int = tempuse - ((tempuse*0.77-0.05)+ 0.02 * tempuse);
     max_int = 0.1;									//rate of interception
     max_int = max_int * rain; // covert back to mm, it is the interception potential of leaf

   if (rain  <= max_int)          /* all intercepted */
		{
			prec_to_canopy[vt][dm] = rain;
			rainthrough[vt][dm] = 0.0;
		}
		else                      /* canopy limits interception */
		{
			prec_to_canopy[vt][dm] = max_int;
			rainthrough[vt][dm] = rain - max_int;
		}
//printf("%9.3f\n",0.111111111);
    //  return prec_to_canopy;
};


/* *************************************************************
************************************************************* */

//  Evaporation of the first layer of soil during wet (Maybe moss or bair ground, depeding on the site
double Hydrology::evaporation_layer1(const double& airtemp, const double& pa,  double& lanmuda, double& vapord, double& soilrad)
{

// needed inputs:
// airtemp - air temperature (oc)
// pa - density of air  (kg/m3)
//lamuda - latent heat of vaporization (J/kg)
// vapord - vapor pressure deficit  (Kpa)
// soilr - radiation projected to the soil surface (MJ m-2 day-1)

double evap_soil; // soil evaporation rate mm/s
double eeq; // equlibrium evaporation mm/s
double yipucilong; // the change of latent heat
double gb;  //the boundary layer conductance (which is soil surface) mm/s
double zeta; // as function of serveral parameters
double Tk; // air temperature in degrees Kelvin
double secondpart; // evaporation rate resulted from vapord mm/s

double Gv = 0.462; //  m3 kpa kg-1 K-1, as gas constant for water vapor

soilrad = soilrad * pow(10,6) / (24 * 3600); // transfer MJ m-2 day-1 to J m-2 s-1
soilrad = 0.001 * soilrad; // assume only 0.1% radiation to the soil surface
gb = 10.0; // 10mm/s for soil surface about aerodynamics, see waring et al, 1997, page. 27 	Niger: aero_resis_layer1 = 107 s m-1 (Wallace and Holwill, 1997).
Tk = airtemp + 273.0;  // change airtemp (oc) to kelvin, later on need to replace Tk with soil surface temperature from STM
yipucilong = 0.7185 * exp(0.0544 * airtemp);
zeta = pa * (yipucilong + 1.0) * Gv * Tk;
vapord = vapord * pow(10,2)/pow(10,3);  // bar = 10^5pa; mbar = 10^2pa, trasfer mbar to kpa
vapord = 0.001 * vapord; // assume 0.1% of canopy vapord
/*
vapord = 0.1;
soilrad = 5.0; // w/s-2 = J/s-2/s-1
Tk = 10+273;
gb = 0.1; //mm/s
yipucilong = 0.7185 * exp(0.0544 * 10);
lanmuda =  2.5023e6 - 2430.54 * 10;
zeta = (1.292 - (0.00428 * 10)) * (yipucilong + 1.0) * Gv * Tk;
double pa_temp = 1.292 - (0.00428 * 10);
// need to think the vapord value????
*/

secondpart = vapord * gb /zeta; //mm/s
eeq = soilrad * yipucilong / (pa * lanmuda * (yipucilong + 1.0)) *1000.0; // transfer rate m/s to mm/s
evap_soil = eeq + secondpart; // mm/s total evaporation rate,
evap_soil = evap_soil * 12 * 3600; // daytime evaporation, mm/day
//printf(" Ok = %10.5f", evap_soil);

return evap_soil;
}

/* *************************************************************
************************************************************* */

//Sublimation from ground snowpack based on Coughlan PhD thesis, 1991
double Hydrology::snowsublimation_ground( double snowwatereq, const double& airtemp, double surfacerad)
{
 static double latent_s = 2845.0; // (KJ/kg) latent heat of sublimation
 static double snowabs = 0.6; // radiation absorptivity of snow;
 double incident_rad; //incident radiation to the snowpack KJ/m2/day
 double sublimation; // actual sublimation, mm/day
 double potentialsub; // potential sublimation of snow

 incident_rad = surfacerad * snowabs * pow(10,3); // convert MJ/m2/day to KJ/m2/day
 if ((airtemp < 0.0) && (snowwatereq > 0.0))
 {
   potentialsub = incident_rad / latent_s;
    if (potentialsub > snowwatereq)
     potentialsub = snowwatereq;
     sublimation = potentialsub;
     sublimation *= 1000; // transfer m/day to mm/day
    }
 else sublimation = 0.0;

 return sublimation;
}
/* *************************************************************
************************************************************* */
// sublimation from canopy snow, Coughlan and Running 1997
double Hydrology::sublimation_canopy(double lai, double canopy_rad, const double airt)
{
  double sub_canopy; // mm/day
  double latentv; // latent heat of vaporization (MJmm-1M-2)
  static double snow_int = 0.5; // mm/LAI
  static double k = 3.5e5; // MJmm-1M-2 latent heat fusion

  latentv = 2.5023e6 - 2430.54 * airt; // J/kg
  latentv = (latentv / 10e6) * 10e9 / 10e6; // MJ mm-1M-2;
//  canopy_rad = canopy_rad / 0.2388 / 10e6 * 10e4; // convert cal cm-2 day-1 to MJ m-2 day-1
  //sub_canopy = min((lai * snow_int/2.0), (canopy_rad / (latentv + k))*1000 );

  if ((lai * snow_int/2.0) > ((canopy_rad / (latentv + k))*1000 ))
  sub_canopy = ((canopy_rad / (latentv + k))*1000 );
  else sub_canopy = (lai * snow_int/2.0);

  return sub_canopy;
}
/* *************************************************************
************************************************************* */

//snow melt driven by incident radiation, Coughlan and Running 1997.
double Hydrology::snowmelt_srad(double snow_rad, const double airt)
{
static double albedo = 0.80; // dimensionless Running and Coughlan 1988
static double k = 3.5e5; // MJmm-1M-2 latent heat fusion
double meltmaxday; // melt rate day mm/day
  snow_rad = ((snow_rad /0.2388) / 10e6) * 10e4; // covert Cal cm-1 day to MJm-2day-1
  meltmaxday = albedo * (snow_rad / k);
  return meltmaxday;
}

/* *************************************************************
************************************************************* */


void Hydrology::getecd(char ecd[80])
{
  getpenmonecd(ecd);

}
void Hydrology::getecd(ofstream& rflog1)
{

  char ecd[80];

  cout << "Enter name of the file with Penman-Monteith parameter values (.ECD):" << endl;
  fpara >> ecd;

  rflog1 << "Enter name of the file with Penman-Monteith parameter values (.ECD):";
  rflog1 << ecd << endl << endl;

  getpenmonecd(ecd);

};


/* **************************************************************
************************************************************** */


void Hydrology::getpenmonecd(char ecd[80])
{

  const int NUMVAR = 9;
  char dummy[35][40];
  ifstream infile;
  int i;
  int dcmnt;

  long update[35];

  infile.open(ecd, ios::in);

  if (!infile)
  {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 1; dcmnt < 35; dcmnt++)
  {

  infile >> dummy[dcmnt] >> dummy[dcmnt]>> inter_coef[dcmnt] >> EXT[dcmnt]>> LWPmin[dcmnt] >> CCmax[dcmnt] >> LWPc[dcmnt] >> SLOPEcc[dcmnt]>> update[dcmnt];

   }

  infile.close();

};

/* *************************************************************
************************************************************* */


void Hydrology::getpenmonecd(ofstream& rflog1) {

  const int NUMVAR = 9;
  char dummy[NUMVAR][8];
  ifstream infile;
  int i;
  int dcmnt;
  char ecd[20];

  int update[12];

  infile.open(ecd, ios::in);

  cout << "Enter name of Penman Monteith (.ECD) data file with parameter values:" << endl;
  fpara >> ecd;

  rflog1 << "Enter name of Penman Monteith(.ECD) data file with parameter values:" << ecd << endl << endl;

  infile.open(ecd, ios::in);

  if (!infile) {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 1; dcmnt < 35; dcmnt++)
  {
    infile >> dummy[0] >> dummy[0]>> inter_coef[dcmnt] >> EXT[dcmnt]
    >> LWPmin[dcmnt] >> CCmax[dcmnt] >> LWPc[dcmnt]>>
     SLOPEcc[dcmnt]>> update[dcmnt];
  }


  infile.close();

};


// ******

void Hydrology::showpenmon(int& dcmnt)
{

  cout << endl << "         PARAMETERS FOR THE PENMAN-MONTEITH EQUATION";
  cout << endl << endl;
  printf("     INTER_COEF = %7.5lf     EXT = %7.5lf         LWPmin = %8.5lf\n",
         inter_coef[dcmnt], EXT[dcmnt], LWPmin[dcmnt]);
  printf("   CCmax = %4.2lf       LWPc = %7.4lf   \n",
         CCmax[dcmnt], LWPc[dcmnt]);

  printf("   SLOPEcc = %4.2lf\n",
         SLOPEcc[dcmnt] );
};
