/* **************************************************************
*****************************************************************
ATMS423.HPP - object describes physical characteristics of the
	      atmosphere

           - modified from atms41d.hpp by DWK 20000102
           - added vap for hydrology model by QZ 20000717
*****************************************************************
************************************************************** */

// Class uses global constants CYCLE and MAXRTIME

#if !defined(ATMS423_H)
#define ATMS423_H
class Atmosphere
{

  public:

     Atmosphere();

/* **************************************************************
		 Public Functions
************************************************************** */


     void precsplt(double& prec, double& tair, double& rain, double& snowfall);
     double petjh(const double& nirr, const double& tair, const int& dm);
     double xeet(const double& rain, const double& snowinf, const double& pet, const double& avlh2o, const double& awcapmm, const int& dm);
     double xgirr(double& lat, int& dm, double& sumday);
     double xnirr(const double& clds, const double& girr);
     double xpar(const double& clds, const double& nirr);



/* **************************************************************
		 Public Variables
************************************************************** */

// Gaseous components of air

     double initco2;               // initial CO2 concentration (ppmv)
     double co2level;              // constant CO2 concentration (ppmv)
     double co2[CYCLE];            // variable CO2 concentration (ppmv)
     double tco2[MAXRTIME][CYCLE]; // transient CO2 concentration (ppmv)          
     long co2year[MAXRTIME];       // year of CO2 data
     double latt;
    
 //  Light Variables

     double girr[NVT][CYCLE];        // Gross Irradiance (cal/(sq. cm * day))
     double nirr[NVT][CYCLE];        // Net Irradiance   (cal/(sq. cm * day))
     double nirr1[NVT][CYCLE];        // Net Irradiance   (cal/(sq. cm * day))
     double par[NVT][CYCLE];         // Photosynthetically Active Radiation
   //  double par1[NVT][CYCLE];         // Photosynthetically Active Radiation
				      // (cal/(sq.cm * day))
//   Evapotranspiration Variables

     double pet[NVT][CYCLE];         // Potential Evapotranspiration (mm)     
     double eet[NVT][CYCLE];         // Estimated Actual Evapotranspiration (mm)


//   Climatic Variables

     double tair[NVT][CYCLE];        // Surface Air Temperature (degrees C)
     double ttair[MAXRTIME][CYCLE];
     long tairyear[MAXRTIME];       // year of air temperature data

     double frontd[NVT][CYCLE],thawbe[NVT][CYCLE],thawend[NVT][CYCLE]; // front frozen depth
     double tsoil[NVT][CYCLE];       // Top 20cm soil temperature (degrees C)
     double dst5[NVT][CYCLE],dst10[NVT][CYCLE],dst20[NVT][CYCLE],dst50[NVT][CYCLE],dst100[NVT][CYCLE],dst200[NVT][CYCLE]; // for different depth soil temperature

     double clds[NVT][CYCLE];        // Cloudiness (%)
     double tclds[MAXRTIME][CYCLE];
     long cldsyear[MAXRTIME];       // year of cloudiness data

 //    double sun[CYCLE];         // percent sunshine duration

     double  prec[NVT][CYCLE];       // Total Precipitation (mm)
     double  tprec[MAXRTIME][CYCLE];
     long precyear[MAXRTIME];       // year of precipitation data
     
	// double tclds[MAXRTIME][CYCLE];    //wsr: just PAR input

     double  snowfall[NVT][CYCLE];   // Snow (mm)
     double  rain[NVT][CYCLE];       // Rainfall (mm)

     double prevtair;           // Previous Month's Air Temperature
				//   (degrees C)
     double prev2tair;          // Previous 2 Month's Air Temperature
				//   (degrees C)

     double prvpetmx;           // Maximum PET of previous year
     double prveetmx;           // Maximum EET of previous year

     double avetair[NVT];            // Mean annual air temperature (degrees C)
     double mxtair;             // Maximum monthly air temperature (degrees C)
     double mxttair[MAXRTIME];
     double yrprec[NVT];             // Annual sum of total precipitation (mm)
     double yrtprec[MAXRTIME];
     double yrrain[NVT];             // Annual sum of rainfall (mm)
     double yrpet[NVT];              // Annual sum of potential evapotranspiration (mm)
     double yreet[NVT];              // Annual sum of estimated actual evapotranspiration (mm)
     double yrsnowfall[NVT];         // Annual sum of snow

     double yrfrontd[NVT],yrthawbegin[NVT],yrthawend[NVT]; // for soil thermal model
     double yrtsoil[NVT]; // for soil temperature
     double yrdst5[NVT]; // for soil temperature
     double yrdst10[NVT]; // for soil temperature
     double yrdst20[NVT]; // for soil temperature
     double yrdst50[NVT]; // for soil temperature
     double yrdst100[NVT]; // for soil temperature
     double yrdst200[NVT]; // for soil temperature

     int tco2flag;
     int ttairflag;
     int tprecflag;
     int tcldsflag;
     int nirrflag;
	 int tsoiloggcflag;     // wsr: added
	 int tpeatflag;
         int RGflag;
 // added for hydrology model by QZ
     double vap[NVT][CYCLE]; // Vapor pressure (Kpa)

 // private:

/* **************************************************************
		      Private Variables
************************************************************** */

     double daze[CYCLE];          // number of days per month

};

#endif
