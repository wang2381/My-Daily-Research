/* **************************************************************
*****************************************************************
TMCRB423.HPP - object describing characteristics of soil microbes
	       used in the Terrestrial Ecosystem Model (TEM)
*****************************************************************
************************************************************** */

// Tmicrobe4 uses the global constants CYCLE, NUMMSAC and NUMVEG

#if !defined(TMCRB423_H)
#define TMCRB423_H

class Tmicrobe4 {

   public:

/* **************************************************************
		 Public Functions
************************************************************** */

      void   getvegecd(ofstream& rflog1);
      void   getvegecd(char ecd[80]);
      double nminxclm(const int& dcmnt, const int& dm, double& soilh2o,
                      double& soilorgc, double& soilorgn, double& availn,
                      double& decay, double& rh, double& ksoil, const int& vt);
      double rhxclm(double& soilorgc, double& dq10, double& moist, const int& vt);
      void   showecd(const int& dcmnt);
      // ez changed to dcmnt in yrkd() by DWK on 20000130
      double yrkd(int nfeed, double& yrltrc, double& yrltrn, int dcmnt, const int& vt);


/* **************************************************************
		 Public Variables
************************************************************** */

      // monthly heterotrophic respiration (g C / (sq. meter * month))

      double rh[NVT][CYCLE];

      // monthly nitrogen uptake by microbes (g N / (sq. meter * month))

      double nuptake[NVT][CYCLE];

      // monthly net nitrogen mineralization (g N / (sq. meter * month))

      double netnmin[NVT][CYCLE];

      // annual heterotrophic respiration (g C / (sq. meter * year))

      double yrrh[NVT];

      // annual net nitrogen mineralization (g N / (sq. meter * year))

      double yrnmin[NVT];


/* **************************************************************
		 Public Parameters
************************************************************** */

      // Parameter representing the quality of soil organic matter

      double decay[NVT];

      // Biome-specific decomposition parameters for function rhxclm

      double kd[NVT];
      double kda[MAXCMNT];
      double kdb[MAXCMNT];
      double kdc[NVT];
      
      double kda_p[MAXCMNT];
      double kdb_p[MAXCMNT];
      
      double kda_temp[MAXCMNT];
      double kdb_temp[MAXCMNT];

      double kdin[NUMMSAC];     // kd values read in from file
      double kdsave[NUMMSAC];   // kd values saved to a file

      double lcclnc[MAXCMNT];
      double propftos[MAXCMNT];
      double lcclnc_p[MAXCMNT];
      double propftos_p[MAXCMNT];
      double lcclnc_temp[MAXCMNT];
      double propftos_temp[MAXCMNT];

      // Biome-specific microbial nitrogen uptake parameters for function
      //   nminxclm

      double nup[NVT];
      double nupa[MAXCMNT];
      double nupb[MAXCMNT];
      double nupa_p[MAXCMNT];
      double nupb_p[MAXCMNT];
      double nupa_temp[MAXCMNT];
      double nupb_temp[MAXCMNT];

      // Biome-specific N fixation parameter

      double nfixpar[MAXCMNT];
      double nfixpar_p[MAXCMNT];
      double nfixpar_temp[MAXCMNT];
      double cnsoil[MAXCMNT];
      double cnsoil_p[MAXCMNT];
      double cnsoil_temp[MAXCMNT];
      


      double rhq10[MAXCMNT];

      //Biome-specific parameters describing the influence of soil moisture
      // on decomposition (i.e., moist)

      double moistmin[MAXCMNT];
      double moistopt[MAXCMNT];
      double moistmax[MAXCMNT];

      // Biome-specific half saturation parameter for function nminxclm
      //   describing the effect of available nitrogen on microbial
      //   nitrogen uptake

      double kn2[MAXCMNT];

};

#endif

