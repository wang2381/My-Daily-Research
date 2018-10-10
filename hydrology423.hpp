/* **************************************************************
*****************************************************************
HYDROLOGY423.HPP - object describes physical characteristics of the
	      atmosphere

           - created by QZ 06162000
*****************************************************************
************************************************************** */

// Class uses global constants CYCLE and MAXRTIME


class Hydrology
{

  public:

/* **************************************************************
		 Public Functions
************************************************************** */
    Hydrology();  // initialization

     void showpenmon(int& dcmnt);
     void getecd(char ecd[80]);
     void getecd(ofstream& rflog1);
     void getpenmonecd(ofstream& rflog1);
     void getpenmonecd(char ecd[80]);

    void drad(const double& LAI, const double& EXT, const double& nirr,const int& dm, const int& vt);
   // double LWP(const double& SOILH2O, const double& SOILCAP);
	double LWP(double& SOILH2O1,double& SOILCAP, const double sand, const double clay, const double silt);

    double CanopyCond(const double& LWPmin, const double& LWP,const double& CCmax,
           const double& LWPc, const double& ABSHD, const double& SLOPEcc);
    double slope(const double& airt);
    double densityA(const double& airt);
    double latentV(const double& airt);
    double penmanmonteith(const double& slope, double drad, const double& pa,
           const double& vpd, const double& cch, const double& le);
    double evaporation(const double& water_on_canopy);
    void intercept(const double& rain, const double& inter_coef, const int& dm, const int& vt);

    void vpd(const double& vaporP, const double& airt, const int& vt);
    double evaporation_layer1(const double& airtemp, const double& pa,  double& lanmuda, double& vapord, double& soilrad);
    double snowsublimation_ground( double snowwatereq, const double& airtemp, double surfacerad);
    double snowmelt_srad(double snow_rad, const double airt);
	double sublimation_canopy(double lai, double canopy_rad, const double airt);


/* **************************************************************
		 Public Variables
************************************************************** */

double yrhevap[NVT];
double yrhtrans[NVT];
double yrsevap[NVT];
double yrsnowsub[NVT];
double yrsubcan[NVT];
 //  pressure Variables

 double vpdeficit[NVT]; // kpa
 double absoluteHD[NVT]; // g/m3
 double vaporpressure[NVT][CYCLE];
 double snowpack_canopy[NVT][CYCLE]; // mm/month
// for soil surface evaporation
 double sub_from_canopy[NVT][CYCLE]; // sublimation from canopy mm/day
 double snowsub[NVT][CYCLE]; // snow sublimation mm/day (water equivalent) from ground
 double soil_evap[NVT][CYCLE]; // actual evaporation of soil surface mm/day
 double ratio;
 double potevap_surface[NVT][CYCLE]; // potential soil surface evaporation mm/day
 int dayssincerain; // days since rain
 int sublimationdays[NVT][CYCLE]; // canopy sublimation days
// end of soil surface

//   double ppet[CYCLE];         // Potential Evapotranspiration (mm)

// Evaporation Variables // added by QZ 14 June 2000
     double hdrad[NVT][CYCLE]; // canopy average daily radiation  (MJ m-2 day-1)
     double hLWP[NVT][CYCLE];  // daily maximum leaf water potential (MPa)
     double hcch[NVT][CYCLE]; // canopy conductance with humidity reduction (ms-1)
     double hslope[NVT][CYCLE]; // slope of the saturation vapor pressure curve
     double hpa[NVT][CYCLE]; // air density
     double hlv[NVT][CYCLE]; // latent heat of vapor
     double htrans[NVT][CYCLE]; // transpiration (mmday-1)
     double hevap[NVT][CYCLE];  // Estimated Actual Evaporation by canopy and bare soil (?)(mm)
     double surfacerad[NVT][CYCLE]; // soil surface radiation MJ m-2 day-1

//   Climatic Variables

 //    double  rain[CYCLE];       // Rainfall (mm)
     double  rainthrough[NVT][CYCLE];       // Rainfall to soil ground (mm)
     double  canopy_water[NVT][CYCLE]; // rain on canopy (mm) // added by Qianlai zhuang 14 June 2000

     double  prec_to_canopy[NVT][CYCLE];

// added for hydrological model
     double inter_coef[35], hvpd[35], EXT[35],dewpoint[NVT][CYCLE];
     double CCmax[35],LWPc[35],SLOPEcc[35];
     double hle,LWPmin[35];

  // added for hydrology model
   int vapflag;

   double  vap[NVT][CYCLE];       // Total vap
   double  tvap[MAXRTIME][CYCLE];
   long    vapyear[MAXRTIME];       // year of vap  data

   double yrvap[NVT];             // Annual sum of total vap
   double yrtvap[MAXRTIME];

   double mxvap;             // Maximum monthly vap
   double mxtvap[MAXRTIME];

   double daze[CYCLE];          // number of days per month

  private:
/* **************************************************************
		      Private Variables
************************************************************** */

};

