/* **************************************************************
*****************************************************************
TELM423.HPP - Runs TEM for a single grid cell

Modifications:

20000202 - DWK changed initstat[MAXSTATE] to initstat[MAXSTATE+1];
20000219 - DWK changed ez to cmnt in temgisqc()
20000715 - QZ added or changed for temgisin() for using vapor pressure in hydrology model
*****************************************************************
************************************************************** */

//Modules representing climate and TEM

#if !defined(POTCLM423_H)
  #include "potclm423e.cpp"       // TEMelmnt uses the Potsclm class
#endif

#if !defined(TTEM423_H)
  #include "ttem423e.cpp"      // TEMelmnt uses the TTEM class
#endif



//Modules describing the interface with spatially explicit data sets

#if !defined(TCLMDAT423_H)
  #include "tclmdat423.cpp"  //TEMelmnt uses the Clmdata class
#endif

#if !defined(TVEGDAT423_H)
  #include "tvegdat423.cpp"  //TEMelmnt uses the Vegdata class
#endif

#if !defined(TSOLDAT423_H)
  #include "tsoldat423.cpp"  //TEMelmnt uses the Soildata class
#endif

#if !defined(TELVDAT423_H)
  #include "telvdat423.cpp"  //TEMelmnt uses the Elevdata class
#endif

#if !defined(LULCDAT423_H)
  #include "lulcdat423.cpp"  //TEMelmnt uses the Lulcdata class
#endif

#if !defined(TTEMDAT423_H)
  #include "ttemdat423.cpp"  //TEMelmnt uses the Temdata class
#endif

#if !defined(TKDDAT423_H)
  #include "tkddat423.cpp"   //TEMelmnt uses the KDdata class
#endif

#if !defined(LATDAT423_H)
  #include "latdat423.cpp"   //TEMelmnt uses the Latdata class
#endif

#if !defined(TSOILOGGCDAT_H)
#include "tsoiloggcdat.cpp"   //wsr: add for soil organic carbon input
#endif

#if !defined(TELM423_H)
#define TELM423_H
const int TQCZEROFLAG = 31;
class TEMelmnt {

  public:

     TEMelmnt();

/* *************************************************************
		 Public Function Declarations
************************************************************* */

     int atmsgisin(FILE* fclds, const int& cldflag, FILE* flonlat,
                   const int& lonlatflag, const int& numspin,
                   const int& spintime, const int& RTIME,int& ftlerr, ofstream& flog1);
     void atmswritemiss(ofstream fout[NUMATMS], char predname[MAXPRED][9],
                        const int& dyr, const int& natmspred,
                        const double value);
     void atmswritepred(ofstream fout[NUMATMS], char predname[MAXPRED][9],
                        const int& dyr, const int& natmspred);
     void GISsetELMNTstate(const int& dcmnt, Temdata initstat[NUMMSAC][MAXSTATE+1],const int& vt);

     void runtem(ofstream& rflog1, char predmap[MAXPRED][9],
                 const int& cldflag, const int& atmsflag,
                 const int& atmsoutfg, int& natmspred,
                 ofstream fsradout[NUMATMS],
                 const int& temflag, const int& kdinflg,
                 int& ntempred, ofstream ftemout[MAXPRED],
                 const int& stateflag,
                 const int& equil, const int& totsptime,
                 const int& RTIME, const int& vt);
//added fvap for hydrology model by QZ        //wsr: add fsoiloggc
     int temgisin(ofstream& flog1, int& ftlerr, const int& atmsflag,
		  const int& itype, const int& maxtype, const int& kdinflg,
                  const int& stateflag,
		  FILE* fstxt, FILE*fsoiloggc, FILE*felev, FILE* fnirr, FILE* fpar,
		  FILE* ftair, FILE* fprec, FILE* fvap, FILE* flulc, FILE* fnpp,
                  FILE* fkdin, FILE* fstate[MAXSTATE+1],
		  const int& numspin, const int& spintime, const int& RTIME);
        void temwritemiss(ofstream fout[MAXPRED], char predname[MAXPRED][9],
                            const int& dyr, const int& itype,
                            const int& natmspred, const int& ntempred,
                            const double value);
     void temwritepred(ofstream fout[MAXPRED], char predname[MAXPRED][9],
                       const int& dyr, const int& itype,
                       const int& ntempred, const int& natmspred,const int& vt);
     int veggisin(ofstream& flog1, int& ftlerr, const int& atmsflag,
		  const int& temflag, FILE* ftveg);

     float col;
     float row;
     int lonlatflag;
     double lon;
     double lat;
     int carea;
     char contnent[9];

     int mez;
     int maxtype;

     int outyr;
     int itype;

     long atmstotyr[MAXRTIME];
     long ttotyr[MAXRTIME][NUMMSAC];

     Potsclm clm;
     TTEM tem;


/* *************************************************************
		 Private Function Declarations
************************************************************* */

  private:

     int coregerr(ofstream& rflog1, const char varname1[9], const float& col1,
		  const float& row1, const char varname2[9], const float& col2,
		  const float& row2);
     int loadteclds(FILE* fclds, Clmdata clds, const int& numspin,
                    const int& spintime, const int& RTIME,int& ftlerr, ofstream& flog1);
     int loadtelulc(FILE* flulc, Lulcdata lulc, const int& numspin,
                    const int& spintime, const int& RTIME, int& ftlerr,
		    ofstream& flog1);
     int loadtenpp(FILE* fnpp, Temdata npp, const int& maxtype,
                   const int& RTIME, int& ftlerr, ofstream& flog1);
     int loadteprec(FILE* fprec, Clmdata prec, const int& numspin,
                    const int& spintime, const int& RTIME, int& ftlerr,
		    ofstream& flog1);
     int loadtetair(FILE* ftair, Clmdata tair, const int& numspin,
                    const int& spintime, const int& RTIME, int& ftlerr,
		    ofstream& flog1);

// added for hydrological model

     int loadvap(FILE* fvap, Clmdata vap, const int& numspin,
                  const int& spintime, const int& RTIME, int& ftlerr,
		    ofstream& flog1);


     int temgisqc(const int& stateflag, const double& pctsilt,
                  const double& pctclay, const int& cmnt, const double& elev,
                  double nirr[CYCLE], double par[CYCLE],
		  double tair[CYCLE], double& mxtair,
                  double prec[CYCLE], double& yrprec,
		  Temdata initstat[NUMMSAC][MAXSTATE+1]);
     int transqc(int& maxyears, long& totyr, Biomass plant[CYCLE]);


// *************************************************************

     Temdata initstat[NUMMSAC][MAXSTATE+1];

     int qc;
     int lowtair;
     int noprec;

       //added for hydrological model
     int lowvap;
     int novap;


};

#endif

