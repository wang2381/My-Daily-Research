/* **************************************************************
*****************************************************************
TVEG423.HPP -  Vegetation characteristics used in the Terrestrial
               Ecosystem Model (TEM)
            -  Bug fix added by DWK on 19991028
            -  Additional modifications by DWK on 20000102
20000614 - DWK changes innpp to innpp[CYCLE] and ingpp to ingpp[CYCLE]
*****************************************************************
************************************************************** */

#if !defined(BIOMS423_H)
  #include "bioms423.hpp"   // Tveg42 uses Biomass Class
#endif

#if !defined(TBIOME423_H)
  #include "tbiome423.cpp"   // Tveg42 inherits Biome Class
#endif

// Tveg42 also uses the global constants CYCLE and NUMVEG

#if !defined(TVEG423_H)
#define TVEG423_H

class Tveg42 : public Biome
{

  public:

     Tveg42();

/* **************************************************************
		 Public Functions
************************************************************** */

     double deltaleaf(const int& dcmnt, double& eet, double& prveetmx, double& prvleaf);
     void   getecd(ofstream& rflog1);
     void   getecd(char ecd[80]);
     void   getleafecd(ofstream& rflog1);
     void   getleafecd(char ecd[80]);
//     void   getvtype(ofstream& rflog1);
//     void   getvtype(char ecd[80]);
   // modified for thawing STM, hydrology model, adding thawpercent
     double gppxclm(int& dcmnt, double& co2, double& par, double& temp, double& gv, double& leaf, double& foliage,double& thawpercent, const int& vt);
     void   leafinit(ofstream& rflog1);
     double nupxclm(int& dcmnt, double& soilh2o, double& availn, double& respq10, double& ksoil, double& foliage, const int& vt);
// added for thawing
     double rmxclm(int& dcmnt, double& vegc, double& respq10,double& thawpercent, const int& vt);
     void   showecd(int& dcmnt);
     void   showleaf(int& dcmnt);
     void   updateC2N(const int& dcmnt, double& yreet, double& yrpet, double& currentco2, double& initco2, const int& vt);



/* **************************************************************
		 Public Variables
************************************************************** */

     Biomass plant[NVT][CYCLE];     // whole plant biomass (structural + labile)
     Biomass strctrl[NVT][CYCLE];   // structural plant biomass
     Biomass labile[NVT][CYCLE];    // labile plant biomass
     Biomass ltrfal[NVT][CYCLE];    // litterfall

     double nmobil[NVT][CYCLE];     // N mobilization by plants
     double nresorb[NVT][CYCLE];    // N resorption by plants
     double inuptake[NVT];          // initial N uptake by plants
     double nuptake[NVT][CYCLE];    // N uptake by plants
     double suptake[NVT][CYCLE];    // N uptake by plants for structural N
     double luptake[NVT][CYCLE];    // N uptake by plants for labile N
     double inprodcn[NVT];          // initial C/N of biomass production

     double ingpp[NVT][CYCLE];      // initial gross primary productivity
     double gpp[NVT][CYCLE];        // gross primary productivity (GPP)

     double innpp[NVT][CYCLE];      // initial net primary productivity
     double npp[NVT][CYCLE];        // net primary productivity (NPP)

     double rm[NVT][CYCLE];         // maintenance respiration
     double rg[NVT][CYCLE];         // growth respiration
     double gpr[NVT][CYCLE];        // gross plant respiration (rm + rg)

     double yrcarbon[NVT];          // annual sum of plant.carbon
     double yrnitrogen[NVT];        // annual sum of plant.nitrogen
     double yrstructn[NVT];         // annual sum of strctrl.nitrogen
     double yrstoren[NVT];          // annual sum of labile.nitrogen
     double yrc2n[NVT];             // ratio of yrcarbon to yrnitrogen
     double yrnpp[NVT];             // annual sum of npp
     double yrinnpp[NVT];           // annual sum of innpp
     double yrnup[NVT];             // annual sum of nuptake
     double yrnmobil[NVT];          // annual sum of nmobil
     double yrsup[NVT];             // annual sum of suptake
     double yrlup[NVT];             // annual sum of luptake
     double yrinnup[NVT];           // annual sum of innup
     double yrnrsorb[NVT];          // annual sum of nresorb
     double yrgpp[NVT];             // annual sum of gpp
     double yringpp[NVT];           // annual sum of ingpp
     double yrltrc[NVT];            // annual sum of ltrfal.carbon
     double yrltrn[NVT];            // annual sum of ltrfal.nitrogen
     double yrinpr[NVT];
     double yrprod[NVT];
     double yrunleaf[NVT];          // mean annual unnormalized leaf phenology
     double yrleaf[NVT];             // mean annual normalized leaf phenology

     double unnormleaf[NVT][CYCLE]; // monthly unnormalized leaf phenology
     double leaf[NVT][CYCLE];       // monthly normalized leaf phenology

     double lai[NVT][CYCLE];        // monthly leaf area index
     double fpc[NVT][CYCLE];        // monthly foliar projective cover
     double yrfpc[NVT];
     double alleaf[NVT];
     double foliage[NVT];

     double thawpercent[NVT][CYCLE]; // added for thaw-frown

    // Number of annual iterations for determining monthly phenology

     int leafyrs;

/* **************************************************************
		 Public Parameters
************************************************************** */

     // Biome-specific vegetation C/N parameters

     double cneven[MAXCMNT];
     double cnmin[MAXCMNT];
     double c2n[NVT];
     double c2na[MAXCMNT];
     double c2nb[MAXCMNT];
     double c2nmin[MAXCMNT];
     double initcneven[MAXCMNT];
     double dc2n[NVT];
     double adjc2n[NVT];
     
     double cneven_p[MAXCMNT];
     double cnmin_p[MAXCMNT];
     double c2n_p[NVT];
     double c2na_p[MAXCMNT];
     double c2nb_p[MAXCMNT];
     double c2nmin_p[MAXCMNT];
     double initcneven_p[MAXCMNT];
     double dc2n_p[NVT];
     double adjc2n_p[NVT];
     
     double cneven_temp[MAXCMNT];
     double cnmin_temp[MAXCMNT];
     double c2n_temp[NVT];
     double c2na_temp[MAXCMNT];
     double c2nb_temp[MAXCMNT];
     double c2nmin_temp[MAXCMNT];
     double initcneven_temp[MAXCMNT];
     double dc2n_temp[NVT];
     double adjc2n_temp[NVT];

     // Biome-specific phenology parameters

     double minleaf[MAXCMNT];
     double aleaf[MAXCMNT];
     double bleaf[MAXCMNT];
     double cleaf[MAXCMNT];
     double initleafmx[MAXCMNT];  // Added by DWK on 19991028
     double prvleafmx[MAXCMNT];
     double unleaf12[MAXCMNT];
     
     double minleaf_p[MAXCMNT];
     double aleaf_p[MAXCMNT];
     double bleaf_p[MAXCMNT];
     double cleaf_p[MAXCMNT];
     double initleafmx_p[MAXCMNT];  // Added by DWK on 19991028
     double prvleafmx_p[MAXCMNT];
     double unleaf12_p[MAXCMNT];

     double minleaf_temp[MAXCMNT];
     double aleaf_temp[MAXCMNT];
     double bleaf_temp[MAXCMNT];
     double cleaf_temp[MAXCMNT];
     double initleafmx_temp[MAXCMNT];  // Added by DWK on 19991028
     double prvleafmx_temp[MAXCMNT];
     double unleaf12_temp[MAXCMNT];
     // Biome-specific foliage projection cover parameters

     double fpcmax[MAXCMNT];
     double sla[MAXCMNT];
     double cov[MAXCMNT];

     // Biome-specific allocation parameters

     double kleafc[MAXCMNT];
     double leafmxc[MAXCMNT];
     
     double kleafc_p[MAXCMNT];
     double leafmxc_p[MAXCMNT];
     
     double kleafc_temp[MAXCMNT];
     double leafmxc_temp[MAXCMNT];

     // Biome-specific carbon uptake parameters for function gppxclm

     double cmax[NVT];
     double cmaxcut[MAXCMNT];
     double cmax1a[MAXCMNT];
     double cmax1b[MAXCMNT];
     double cmax2a[MAXCMNT];
     double cmax2b[MAXCMNT];
     
     double cmax_p[NVT];
     double cmaxcut_p[MAXCMNT];
     double cmax1a_p[MAXCMNT];
     double cmax1b_p[MAXCMNT];
     double cmax2a_p[MAXCMNT];
     double cmax2b_p[MAXCMNT];
     
     double cmax_temp[NVT];
     double cmaxcut_temp[MAXCMNT];
     double cmax1a_temp[MAXCMNT];
     double cmax1b_temp[MAXCMNT];
     double cmax2a_temp[MAXCMNT];
     double cmax2b_temp[MAXCMNT];

     // Biome-specific half saturation parameter for function gppxclm
     //   describing the effects of solar atmospheric carbon dioxide
     //   concentration on GPP

     double kc[MAXCMNT];
     double kc_p[MAXCMNT];
     double kc_temp[MAXCMNT];
     // Biome-specific half saturation parameter for function gppxclm
     //   describing the effects of photosybtheically active radiation
     //   on GPP

     double ki[MAXCMNT];

    // Element-specific optimum temperature for GPP

     double topt[NVT];

     double tmin[MAXCMNT];
     double toptmin[MAXCMNT];
     double toptmax[MAXCMNT];
     double tmax[MAXCMNT];

     // Biome-specific parameter to describe the sensitivity of GPP
     //   to evapotranspiration

     double gva[MAXCMNT];

     // Biome-specific respiration parameters for function rmxclm

     double kr[NVT];
     double kra[MAXCMNT];
     double krb[MAXCMNT];
     double kra_p[MAXCMNT];
     double krb_p[MAXCMNT];
     double kra_temp[MAXCMNT];
     double krb_temp[MAXCMNT];
     // Biome-specific parameters for function rq10 to describe the
     //   effect of temperature on plant respiration

     double raq10a0[MAXCMNT];
     double raq10a1[MAXCMNT];
     double raq10a2[MAXCMNT];
     double raq10a3[MAXCMNT];

     // Biome-specific parameters to describe the effect of volumetric
     // soil mositure on GPP

     double vsmmin[NVT];
     double vsmmina[MAXCMNT];
     double vsmminb[MAXCMNT];

     // Biome-specific nitrogen uptake parameters for function nupxclm

     double nmax[NVT];
     double nmaxcut[MAXCMNT];
     double nmax1a[MAXCMNT];
     double nmax1b[MAXCMNT];
     double nmax2a[MAXCMNT];
     double nmax2b[MAXCMNT];
     
     double nmax_p[NVT];
     double nmaxcut_p[MAXCMNT];
     double nmax1a_p[MAXCMNT];
     double nmax1b_p[MAXCMNT];
     double nmax2a_p[MAXCMNT];
     double nmax2b_p[MAXCMNT];
     
     double nmax_temp[NVT];
     double nmaxcut_temp[MAXCMNT];
     double nmax1a_temp[MAXCMNT];
     double nmax1b_temp[MAXCMNT];
     double nmax2a_temp[MAXCMNT];
     double nmax2b_temp[MAXCMNT];

     double kn1[MAXCMNT];
     double kn1_p[MAXCMNT];
     // Biome-specific proportion of vegetation lost as litterfall

     double cfall[MAXCMNT];  // proportion of vegetation carbon
     double nfall[MAXCMNT];  // proportion of vegetation nitrogen
     double cfall_p[MAXCMNT];  // proportion of vegetation carbon
     double nfall_p[MAXCMNT];  // proportion of vegetation nitrogen
     double cfall_temp[MAXCMNT];  // proportion of vegetation carbon
     double nfall_temp[MAXCMNT];  // proportion of vegetation nitrogen


     double labncon[MAXCMNT];

  private:

/* **************************************************************
		 Private Functions
************************************************************** */

     double rq10(int& dcmnt, double& tair);

};

#endif

