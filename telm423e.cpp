/* **************************************************************
*****************************************************************
TELM423.CPP - Runs TEM for a single grid cell

Modifications:

20000214 - DWK changes initstat[itype][0].ave to
           initstat[itype][tem.I_VEGC].ave in runtem()
20000214 - DWK changes tem.ez to tem.veg.cmnt int runtem()
20000214 - DWK added lowtair if condition to temgisin()
20000219 - DWK set mez = 0 when mez < 0 or mez > NUMVEG in veggisin() 
20000219 - DWK changed ez to cmnt in temgisqc()
20000715 -QZ changed for hydrology model
*****************************************************************
************************************************************** */

#if !defined(TELM423_H)
  #include "telm423e.hpp"
#endif

/* ************************************************************ */

TEMelmnt::TEMelmnt()
{

  col = MISSING;
  row = MISSING;
  carea = -999;
  lowtair = 0;
  noprec = 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

int TEMelmnt::atmsgisin(FILE* fclds, const int& cldflag,								//wsr: modify this part in order to input monthly PAR (not clouds) 
                        FILE* flonlat, const int& lonlatflag,
                        const int& numspin, const int& spintime,
                        const int& RTIME,int& ftlerr,ofstream& flog1)
{
int vt;
  int i;
  int gisend=1;
  Clmdata clds;
  Latdata lonlat;

  //if (clm.tcldsflag == 1)															//wsr
  //{
    gisend = loadteclds(fclds, clds, numspin, spintime, RTIME,ftlerr,flog1);
    if (gisend == -1) { return gisend; }
    
//   for (i = 0; i < CYCLE; i++) { for (vt=0; vt< 2; vt++) tem.atms.par[vt][i] = tem.atms.tclds[vt][i]; }
 // }
/*
  else
  {
    gisend = clds.getdel(fclds);
    if (gisend == -1)
    {
      col = MISSING;
      row = MISSING;
      return gisend;
    }
    col = clds.col;
   row = clds.row;
   carea = clds.carea;
    for (i = 0; i < CYCLE; i++)
    {
      if (cldflag == 1) { clm.clds[1][i] = clm.clds[0][i] = clds.mon[i]; }
      else { clm.nirr[1][i] = clm.nirr[0][i] = clds.mon[i]; }
    }
    strcpy(contnent,clds.contnent);
  }
*/
  if (lonlatflag == 1) { lat = (double) row; }
  else
  {
    gisend = lonlat.getdel(flonlat);
    if (gisend != -1) { lat = lonlat.lat; }
    else
    {
      col = MISSING;
      row = MISSING;
      return gisend;
    }
  }

  return gisend;

};
/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TEMelmnt::atmswritemiss(ofstream fout[NUMATMS], char predname[MAXPRED][9],
                             const int& dyr, const int& natmspred,
                             const double value)
{

  int i;
  int dm;
  Clmdata atmspred;

  for (i = 0; i < natmspred; i++)
  {
    for (dm = 0; dm < CYCLE; dm++)
    {
      atmspred.mon[dm] = value;
    }

	 atmspred.outdel(fout[i], col, row, predname[i], carea, atmstotyr[dyr],
                    atmspred.mon, contnent);
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TEMelmnt::atmswritepred(ofstream fout[NUMATMS], char predname[MAXPRED][9],
                             const int& dyr, const int& natmspred)
{

  int i;
  int dm;
  Clmdata atmspred;

  for (i = 0; i < natmspred; i++)
  {
    if (strcmp(predname[i],clm.predstr[clm.I_PAR]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++) { atmspred.mon[dm] = clm.par[0][dm]; }
    }
    else if (strcmp(predname[i],clm.predstr[clm.I_NIRR]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++) { atmspred.mon[dm] = clm.nirr[0][dm]; }
   }
    else if (strcmp(predname[i],clm.predstr[clm.I_GIRR]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++) { atmspred.mon[dm] = clm.girr[0][dm]; }
    }
    else if (strcmp(predname[i],clm.predstr[clm.I_CLDS]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++) { atmspred.mon[dm] = clm.clds[0][dm]; }
    }

	 atmspred.outdel(fout[i], col, row, predname[i], carea, atmstotyr[dyr],
                    atmspred.mon, contnent);
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int TEMelmnt::coregerr(ofstream& rflog1,const char varname1[9],
                       const float& col1, const float& row1,
                       const char varname2[9], const float& col2, const float& row2)
{

  int fatalerr = 0;

  if (col1 != col2 || row1 != row2)
  {
    fatalerr = 1;

    cout << "ERROR:  " << varname1 << " data and ";
    cout << varname2 << "data are not coregistered." << endl;
    cout << "COL = " << col1 << " and ROW = " << row1 << " in ";
    cout << varname1 << " data" << endl;
    cout << "COL = " << col2 << " and ROW = " << row2;
    cout << " in " << varname2 << " data" << endl;

    rflog1 << "ERROR:  " << varname1 << " data and ";
    rflog1 << varname2 << "data are not coregistered." << endl;
    rflog1 << "COL = " << col1 << " and ROW = " << row1 << " in ";
    rflog1 << varname1 << " data" << endl;
    rflog1 << "COL = " << col2 << " and ROW = " << row2 << " in ";
    rflog1 << varname2 << " data" << endl;
  }

  return fatalerr;

};

/* **************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TEMelmnt::GISsetELMNTstate(const int& dcmnt,
                                Temdata initstat[NUMMSAC][MAXSTATE+1],const int& vt)
{

/* **************************************************************
  Function initializes TEM state variables from spatially
explicit data sets
************************************************************** */

  int dm;

  // Initialize ODE carbon, nitrogen, and water state variables
  //   for element



  // Carbon
  tem.y[vt][tem.I_VEGC] = initstat[itype][tem.I_VEGC].mon[CYCLE-1];
  tem.y[vt][tem.I_SOLC] = initstat[itype][tem.I_SOLC].mon[CYCLE-1];

  // Nitrogen
  tem.y[vt][tem.I_STRN] = initstat[itype][tem.I_STRN].mon[CYCLE-1];
  tem.y[vt][tem.I_STON] = initstat[itype][tem.I_STON].mon[CYCLE-1];
  tem.y[vt][tem.I_STON] *= 0.001;
  tem.y[vt][tem.I_SOLN] = initstat[itype][tem.I_SOLN].mon[CYCLE-1];
  tem.y[vt][tem.I_AVLN] = initstat[itype][tem.I_AVLN].mon[CYCLE-1];
  tem.y[vt][tem.I_AVLN] *= 0.001;

  // Water
  tem.y[vt][tem.I_AVLW] = initstat[itype][tem.I_AVLW].mon[CYCLE-1];
  tem.y[vt][tem.I_RGRW] = initstat[itype][tem.I_RGRW].mon[CYCLE-1];
  tem.y[vt][tem.I_SNWPCK] = initstat[itype][tem.I_SNWPCK].mon[CYCLE-1];
  tem.y[vt][tem.I_SGRW] = initstat[itype][tem.I_SGRW].mon[CYCLE-1];
  tem.y[vt][tem.I_SM] = tem.y[vt][tem.I_AVLW] + tem.soil.wiltpt[0];
  tem.y[vt][tem.I_PCTP] = 100.0 * tem.y[vt][tem.I_SM] / tem.soil.totpor[0];
  tem.y[vt][tem.I_VSM] = tem.y[vt][tem.I_SM] / (tem.soil.rootz[0] * 1000.0);
  if (tem.y[vt][tem.I_VSM] <= 0.0) {
    tem.y[vt][tem.I_VSM] = 0.001;
  }

  // for 3-box hydrology
  tem.y[vt][tem.I_AVLW1] = initstat[itype][tem.I_AVLW1].mon[CYCLE-1];
  tem.y[vt][tem.I_AVLW2] = initstat[itype][tem.I_AVLW2].mon[CYCLE-1];
  tem.y[vt][tem.I_AVLW3] = initstat[itype][tem.I_AVLW3].mon[CYCLE-1];

  tem.y[vt][tem.I_SM1] = tem.y[vt][tem.I_AVLW1] + tem.soil.wiltpt1[vt];
  tem.y[vt][tem.I_SM2] = tem.y[vt][tem.I_AVLW2] + tem.soil.wiltpt2[vt];
  tem.y[vt][tem.I_SM3] = tem.y[vt][tem.I_AVLW3] + tem.soil.wiltpt3[vt];

  tem.y[vt][tem.I_PCTP1] =  100.0 * tem.y[vt][tem.I_SM1] / tem.soil.totpor1[vt];
  tem.y[vt][tem.I_PCTP2] =  100.0 * tem.y[vt][tem.I_SM2] / tem.soil.totpor2[vt];
  tem.y[vt][tem.I_PCTP3] =  100.0 * tem.y[vt][tem.I_SM3] / tem.soil.totpor3[vt];

  tem.y[vt][tem.I_VSM1] = tem.y[vt][tem.I_SM1] / (tem.soil.dpwbox1[vt] * 1000.0);
  if (tem.y[vt][tem.I_VSM1] <= 0.0) {
    tem.y[vt][tem.I_VSM1] = 0.001;
  }

  tem.y[vt][tem.I_VSM2] = tem.y[vt][tem.I_SM2] / (tem.soil.dpwbox2[vt]*1000.0);
  if (tem.y[vt][tem.I_VSM2] <= 0.0) {
    tem.y[vt][tem.I_VSM2] = 0.001;
  }

  tem.y[vt][tem.I_VSM3] = tem.y[vt][tem.I_SM3] / (tem.soil.dpwbox3[vt]*1000.0);
  if (tem.y[vt][tem.I_VSM3] <= 0.0) {
    tem.y[vt][tem.I_VSM3] = 0.001;
  }

// end of adding

  // Phenology
  tem.y[vt][tem.I_UNRMLF] = initstat[itype][MAXSTATE].mon[CYCLE-1];
  tem.y[vt][tem.I_UNRMLF] *= 0.01;
  tem.y[vt][tem.I_LEAF] = tem.y[vt][tem.I_UNRMLF] / tem.veg.prvleafmx[dcmnt];


  // Initialize monthly carbon, nitrogen, and water pools for element
  for (dm = 0; dm < CYCLE; dm++)
  {
    // Carbon pools
    tem.veg.plant[vt][dm].carbon = tem.y[vt][tem.I_VEGC];
    tem.soil.org[vt][dm].carbon = tem.y[vt][tem.I_SOLC];
    tem.totalc[vt][dm] = 0.0;

    // Nitrogen pools
    tem.veg.strctrl[vt][dm].nitrogen = tem.y[vt][tem.I_STRN];
    tem.veg.labile[vt][dm].nitrogen = tem.y[vt][tem.I_STON];
    tem.veg.plant[vt][dm].nitrogen = 0.0;
    tem.soil.org[vt][dm].nitrogen = tem.y[vt][tem.I_SOLN];
    tem.soil.availn[vt][dm] = tem.y[vt][tem.I_AVLN];

    // Water pools
    tem.soil.avlh2o[vt][dm] = tem.y[vt][tem.I_AVLW];
    tem.soil.rgrndh2o[vt][dm] = tem.y[vt][tem.I_RGRW];
    tem.soil.snowpack[vt][dm] = tem.y[vt][tem.I_SNWPCK];
    tem.soil.sgrndh2o[vt][dm] = tem.y[vt][tem.I_SGRW];
    tem.soil.moist[vt][dm] = tem.y[vt][tem.I_SM];
    tem.soil.pctp[vt][dm] = tem.y[vt][tem.I_PCTP];
    tem.soil.vsm[vt][dm] = tem.y[vt][tem.I_VSM];

    // for 3-box hydrology
    tem.soil.avlh2o1[vt][dm] = tem.y[vt][tem.I_AVLW1];
    tem.soil.avlh2o2[vt][dm] = tem.y[vt][tem.I_AVLW2];
    tem.soil.avlh2o3[vt][dm] = tem.y[vt][tem.I_AVLW3];

    tem.soil.moist1[vt][dm] = tem.y[vt][tem.I_SM1];
    tem.soil.moist2[vt][dm] = tem.y[vt][tem.I_SM2];
    tem.soil.moist3[vt][dm] = tem.y[vt][tem.I_SM3];

    tem.soil.pctp1[vt][dm] = tem.y[vt][tem.I_PCTP1];
    tem.soil.pctp2[vt][dm] = tem.y[vt][tem.I_PCTP2];
    tem.soil.pctp3[vt][dm] = tem.y[vt][tem.I_PCTP3];

    tem.soil.vsm1[vt][dm] = tem.y[vt][tem.I_VSM1];
    tem.soil.vsm2[vt][dm] = tem.y[vt][tem.I_VSM2];
    tem.soil.vsm3[vt][dm] = tem.y[vt][tem.I_VSM3];
    // end of adding

	


    // adding for soil thermal model

    tem.atms.tsoil[vt][dm] = tem.y[vt][tem.I_TSOIL];
    tem.atms.dst5[vt][dm] = tem.y[vt][tem.I_DST5];
    tem.atms.dst10[vt][dm] = tem.y[vt][tem.I_DST10];
    tem.atms.dst20[vt][dm] = tem.y[vt][tem.I_DST20];
    tem.atms.dst50[vt][dm] = tem.y[vt][tem.I_DST50];
    tem.atms.dst100[vt][dm] = tem.y[vt][tem.I_DST100];
    tem.atms.dst200[vt][dm] = tem.y[vt][tem.I_DST200];
    tem.atms.frontd[vt][dm] = tem.y[vt][tem.I_FRONTD];
    tem.atms.thawbe[vt][dm] = tem.y[vt][tem.I_THAWBE];
    tem.atms.thawend[vt][dm] = tem.y[vt][tem.I_THAWEND];
//end of adding

    // Phenology
    tem.veg.unnormleaf[vt][dm] = tem.y[vt][tem.I_UNRMLF];
    tem.veg.leaf[vt][dm] = tem.y[vt][tem.I_LEAF];

  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

int TEMelmnt::loadteclds(FILE* fclds, Clmdata clds,
                         const int& numspin, const int& spintime,
                         const int& RTIME,int& ftlerr, ofstream& flog1)
{

  int i;
  int j;
  int k;
  int gisend = 1;
  int dm;
  int dyr;
  int totspin;

// Get CLDINESS from transient data set

  totspin = numspin * spintime;
  for (i = totspin; i < RTIME; i++)
  {
    k = i - totspin;
    gisend = clds.getdel(fclds);
    if (gisend == -1)
    {
      col = MISSING;
      row = MISSING;
      return gisend;
   //    ftlerr = coregerr(flog1,"TEMVEG", col, row, "PAR", clds.col, clds.row);
    }

// Determine CLDINESS for equilibrium conditions

    if ( k == 0)
    {
     col = clds.col;
      row = clds.row;
      carea = clds.carea;
      strcpy(contnent,clds.contnent);
     clm.cldsyear[0] = clds.year - totspin;
      for (j = 0; j < CYCLE; j++)
      {
        tem.atms.tclds[0][j] = clds.mon[j];
      }
    }

// Determine CLDINESS for transient conditions

    else
    {
      for (j = 0; j < CYCLE; j++)
      {
        tem.atms.tclds[i][j] = clds.mon[j];
      }
      clm.cldsyear[i] = clds.year;
    }
  }

// Determine CLDINESS during spin up

  for (i = 0; i < numspin; i++)
  {
    dyr = totspin + 1;
    for (j = 0; j < spintime; j++)
    {
      k = (i * spintime) + j + 1;
      clm.cldsyear[k] = clm.cldsyear[0] + k;
      for (dm = 0; dm < CYCLE; dm++)
      {
        tem.atms.tclds[k][dm] = tem.atms.tclds[dyr][dm];
      
      }
      ++dyr;
    }
  }

  return gisend;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

int TEMelmnt::loadtelulc(FILE* flulc, Lulcdata lulc,
                         const int& numspin, const int& spintime,
                         const int& RTIME, int& ftlerr, ofstream& flog1)
{

  int i;
  int j;
  int k;
  int gisend = 1;
  int totspin;

// Get LULC from transient data set

  totspin = numspin*spintime;
  for (i = totspin; i < RTIME; i++)
  {
    k = i - totspin;
    gisend = lulc.getdel(flulc);
    if (gisend == -1) { return gisend; }
    ftlerr = coregerr(flog1,"TEMVEG", col, row, "LULC", lulc.col, lulc.row);

// Determine LULC for equilibrium conditions

    if ( k == 0)
    {
      tem.ag.lulcyear[0] = lulc.year - totspin;
      tem.ag.tstate[0] = lulc.agstate;
      tem.ag.tRAP[0] = lulc.RAP;
    }

// Determine LULC for transient conditions

    else
    {
      tem.ag.lulcyear[i] = lulc.year;
      tem.ag.tstate[i] = lulc.agstate;
      tem.ag.tRAP[i] = lulc.RAP;
    }
  }

// Determine LULC during spin up

  for (i = 0; i < numspin; i++)
  {
    for (j = 0; j < spintime; j++)
    {
      k = (i * spintime) + j + 1;
      tem.ag.lulcyear[k] = tem.ag.lulcyear[0] + k;
      tem.ag.tstate[k] = tem.ag.tstate[0];
      tem.ag.tRAP[k] = tem.ag.tRAP[0];
    }
  }

 return gisend;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

int TEMelmnt::loadtenpp(FILE* fnpp, Temdata npp, const int& maxtype,
                        const int& RTIME, int& ftlerr, ofstream& flog1)
{

// Note: numspin and spintime removed from function call by D. Kicklighter 990724

  int dv;
  int gisend = 1;
  int dm;
  int dyr;
//  int totspin;

// Get POTNPP from transient data set

  for (dyr = 0; dyr < RTIME; dyr++)
  {
    for (dv = 0; dv < maxtype; dv++)
    {
      gisend = npp.getdel(fnpp);
      if (gisend == -1) { return gisend; }
      ftlerr = coregerr(flog1,"TEMVEG", col, row, "POTNPP", npp.col, npp.row);
      for (dm = 0; dm < CYCLE; dm++)
      {
        tem.ag.tpotnpp[dv][dyr][dm] = npp.mon[dm];
      }
    }
  }

  return gisend;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

int TEMelmnt::loadteprec(FILE* fprec, Clmdata prec,
                         const int& numspin, const int& spintime,
                         const int& RTIME, int& ftlerr, ofstream& flog1)
{

  int i;
  int j;
  int k;
  int gisend = 1;
  int dm;
  int dyr;
  int totspin;


// Get PREC from transient data set

  noprec = 0;
  totspin = numspin * spintime;
  for (i = totspin; i < RTIME; i++)
  {
    k = i - totspin;
    gisend = prec.getdel(fprec);
    if (gisend == -1) { return gisend; }
    ftlerr = coregerr(flog1,"TEMVEG", col, row, "PREC", prec.col, prec.row);

// Determine PREC for equilibrium conditions

    if ( k == 0)
    {
      for (j = 0; j < CYCLE; j++) { tem.atms.tprec[0][j] = prec.mon[j]; }
      tem.atms.precyear[0] = prec.year - totspin;
      tem.atms.yrtprec[0] = prec.total;
      if (prec.total <= 0.0) { noprec = 1; }
    }

// Determine PREC for transient conditions

    else
    {
      for (j = 0; j < CYCLE; j++) { tem.atms.tprec[i][j] = prec.mon[j]; }
      tem.atms.precyear[i] = prec.year;
      tem.atms.yrtprec[i] = prec.total;
    }
  }

// Determine PREC during spin up

  for (i = 0; i < numspin; i++)
  {
    dyr = totspin + 1;
    for (j = 0; j < spintime; j++)
    {
      k = (i * spintime) + j + 1;
      tem.atms.precyear[k] = tem.atms.precyear[0] + k;
      tem.atms.yrtprec[k] = tem.atms.yrtprec[dyr];
      for (dm = 0; dm < CYCLE; dm++)
      {
        tem.atms.tprec[k][dm] = tem.atms.tprec[dyr][dm];
      }
      ++dyr;
    }
  }

  return gisend;

};

/* *************************************************************
************************************************************** */

/* *************************************************************
************************************************************** */

int TEMelmnt::loadtetair(FILE* ftair, Clmdata tair,
                         const int& numspin, const int& spintime,
                         const int& RTIME, int& ftlerr, ofstream& flog1)
{

  int i;
  int j;
  int k;
  int gisend = 1;
  int dm;
  int dyr;
  int totspin;

// Get TAIR from transient data set

  lowtair = 0;
  totspin = numspin * spintime;
  for (i = totspin; i < RTIME; i++)
  {
    k = i - totspin;
    gisend = tair.getdel(ftair);
    if (gisend == -1) { return gisend; }
    ftlerr = coregerr(flog1,"TEMVEG", col, row, "TAIR", tair.col, tair.row);

// Determine TAIR for equilibrium conditions

    if ( k == 0)
    {
      tem.atms.tairyear[0] = tair.year - totspin;
      tem.atms.mxttair[0] = tair.max;
      if (tair.max < -1.0) { lowtair = 1; }
      for (j = 0; j < CYCLE; j++)
      {
        tem.atms.ttair[0][j] = tair.mon[j];
      }
    }

// Determine TAIR for transient conditions

    else
    {
      tem.atms.tairyear[i] = tair.year;
      tem.atms.mxttair[i] = tair.max;
      for (j = 0; j < CYCLE; j++) { tem.atms.ttair[i][j] = tair.mon[j]; }
    }
  }
// Determine TAIR during spin up

  for (i = 0; i < numspin; i++)
  {
    dyr = totspin + 1;
    for (j = 0; j < spintime; j++)
    {
      k = (i * spintime) + j + 1;
      tem.atms.tairyear[k] = tem.atms.tairyear[0] + k;
      tem.atms.mxttair[k] = tem.atms.mxttair[dyr];
      for (dm = 0; dm < CYCLE; dm++)
      {
        tem.atms.ttair[k][dm] = tem.atms.ttair[dyr][dm];
      }
      ++dyr;
    }
  }

  return gisend;

};

/* *************************************************************
************************************************************** */

// added for hydrology model

int TEMelmnt::loadvap(FILE* fvap, Clmdata vap,
                           const int& numspin, const int& spintime,
                           const int& RTIME, int& ftlerr, ofstream& flog1)
 {

  int i;
  int j;
  int k;
  int gisend = 1;
  int dm;
  int dyr;
  int totspin;


// Get vap from transient data set

  lowvap = 0;
  totspin = numspin * spintime;
  for (i = totspin; i < RTIME; i++) {
    k = i - totspin;
    gisend = vap.getdel(fvap);
    if (gisend == -1) { return gisend; }
    ftlerr = coregerr(flog1,"TEMVEG", col, row, "VAP", vap.col, vap.row);

// Determine vap for equilibrium conditions

    if ( k == 0) {
      tem.hyd.vapyear[0] = vap.year - totspin;
      tem.hyd.mxtvap[0] = vap.max;
      if (vap.max < -1.0) { lowvap = 1; }
      for (j = 0; j < CYCLE; j++) {
        tem.hyd.tvap[0][j] = vap.mon[j];
      }
    }

// Determine vap for transient conditions

    else {
      tem.hyd.vapyear[i] = vap.year;
      tem.hyd.mxtvap[i] = vap.max;
      for (j = 0; j < CYCLE; j++) { tem.hyd.tvap[i][j] = vap.mon[j]; }
    }
  }

// Determine Vap during spin up

  for (i = 0; i < numspin; i++) {
    dyr = totspin + 1;
    for (j = 0; j < spintime; j++) {
      k = (i * spintime) + j + 1;
      tem.hyd.vapyear[k] = tem.hyd.vapyear[0] + k;
      tem.hyd.mxtvap[k] = tem.hyd.mxtvap[dyr];
      for (dm = 0; dm < CYCLE; dm++) {
        tem.hyd.tvap[k][dm] = tem.hyd.tvap[dyr][dm];
      }
      ++dyr;
    }
  }

  return gisend;

};



/* *************************************************************
************************************************************* */

void TEMelmnt::runtem(ofstream& rflog1, char predmap[MAXPRED][9],
                      const int& cldflag, const int& atmsflag,
                      const int& atmsoutfg, int& natmspred,
                      ofstream fsradout[NUMATMS],
                      const int& temflag, const int& kdinflg,
                      int& ntempred, ofstream ftemout[MAXPRED],
                      const int& stateflag,
                      const int& equil, const int& totsptime,
                      const int& RTIME,const int& vt)
{
int dv;
dv=0;
																						//wsr: re-input the parameters for upland and peatland for each pixel
tem.getsitecd(dv,tem.ecd1);
tem.getsitecd_p(dv,tem.ecd2);
//  int i;
//  int j;
//  int k;
//  int dt;
  int dm;
  int dyr = 0;
  int wrtyr;
//  int maxtime;

  int tqc;
 
 tem.soil.rhflag = 0;                    // = 1 when transient

  qc = ACCEPT;
  tem.totyr = 0;

  if (atmsflag == 1)
  {

// Check cloudiness input for valid data

    for (dm = 0; dm < CYCLE; dm++)
    {
      if (cldflag == 0 && clm.nirr[vt][dm] <= -99.0) { qc = 1; }
      if (cldflag == 1 && tem.atms.tcldsflag == 0 && clm.clds[vt][dm] <= -99.0)
      {
        qc = 2;
      }
      if (cldflag == 1 && tem.atms.tcldsflag == 1 && tem.atms.tclds[vt][dm] <= -99.0)
      {
        qc = 3;
      }
    }

    if (qc == ACCEPT)
    {

    tem.atms.latt = lat;											// wsr: get lattitude for RG to PAR calculation if needed



// Calculate GIRR from POTSCLM model
/*
      clm.yrsumday = 0.0;
      for (dm = 0; dm < CYCLE; dm++)
      {
        clm.girr[vt][dm] = clm.xgirr(lat,dm,clm.yrsumday);
      }
    }
    else { rflog1 << "cldqc = " << qc << endl; }
    */
	}
  }


/* *************************************************************
		BEGIN VEGETATION MOSAIC LOOP
************************************************************* */

  for (itype = 0; itype < maxtype; itype++)
  {
    qc = ACCEPT;
    tem.veg.cmnt = tem.veg.subtype[mez][itype];
/*
// Calculate NIRR, CLDINESS and PAR for Equilibrium Conditions with POTSCLM

    if (atmsflag == 1)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     if (tem.atms.tcldsflag == 1) { clm.clds[vt][dm] = tem.atms.tclds[vt][dm]; }
	     if (cldflag == 1)
        {
	       clm.nirr[vt][dm] = clm.xnirr(clm.clds[vt][dm], clm.girr[vt][dm]);


        }
	     else
        {
	       clm.clds[vt][dm] = clm.mkclds(clm.girr[vt][dm], clm.nirr[vt][dm]);
        }
	//clm.par[vt][dm]  = clm.xpar(clm.clds[vt][dm], clm.nirr[vt][dm]);               //wsr: here we use input PAR instead of cloud
	clm.par[vt][dm] = clm.clds[vt][dm];
      }

// Output POTSCLM steady state results to files

      if (atmsoutfg == 1 && itype == 0)
      {
	     dyr = 0;
	     atmswritepred(fsradout, predmap, dyr, natmspred);
      }
    

*/

// Initialize TEM output variables to missing value

    if (temflag == 1)
    {

      tem.microbe.kdsave[itype] = -9.9999;
      if (kdinflg == 1) { tem.microbe.kdc[vt] = tem.microbe.kdin[itype]; }

// "Hand-off" POTSCLM equilibrium output and equilibrium input data to TEM

      for (dm = 0; dm < CYCLE; dm++)
      {
	     if (atmsflag == 1)
        {
	       tem.atms.nirr[vt][dm] = clm.nirr[vt][dm];
	       tem.atms.par[vt][dm] = tem.atms.tclds[vt][dm];
	     }
	     tem.atms.co2[dm] = tem.atms.co2level;

        // Update TAIR with current year data from transient data files
	     if (tem.atms.ttairflag == 1)
        {
          tem.atms.tair[vt][dm] = tem.atms.ttair[vt][dm];
        }

        // Update PREC with current year data from transient data files
	     if (tem.atms.tprecflag == 1)
        {
          tem.atms.prec[vt][dm] = tem.atms.tprec[vt][dm];
        }

        // added for hydrology model by QZ
       if (tem.hyd.vapflag == 1)
       {
        tem.hyd.vap[vt][dm] = tem.hyd.tvap[vt][dm];
        }


        // Update LULC with current year data from transient data files
   	  if (tem.ag.tlulcflag == 1)
        {
          // if rap is 0 (RAP0flag == 1), then we don't need tpotnpp
	       if(tem.ag.RAP0flag == 0)
          {
  	         tem.ag.potnpp[vt][dm] = tem.ag.tpotnpp[itype][dyr][dm];
	       }
          else { tem.ag.potnpp[vt][dm] = 0.0; }
	     }
      }

      if (tem.atms.ttairflag == 1)
      {
        tem.atms.mxtair = tem.atms.mxttair[vt];
        if (tem.atms.mxtair < -1.0) { lowtair = 1; }
      }
      if (tem.atms.tprecflag == 1)
      {
        tem.atms.yrprec[vt] = tem.atms.yrtprec[vt];
        if (tem.atms.yrprec[vt] <= 0.0) { noprec = 1; }
      }

// Check TEM input for valid data

      qc = temgisqc(stateflag,tem.soil.pctsilt, tem.soil.pctclay,
                    tem.veg.cmnt, tem.elev, tem.atms.nirr[vt], tem.atms.par[vt],
                    tem.atms.tair[vt], tem.atms.mxtair, tem.atms.prec[vt],
                    tem.atms.yrprec[vt], initstat);


// Check TEM parameters for specific vegetation types

      if (qc != ACCEPT) { rflog1 << "temgisqc = " << qc << endl; }
      else { qc = tem.ecdqc(tem.veg.cmnt); }

      if (qc != ACCEPT)
      {

// Skip TEM if input data and/or parameters are missing

        rflog1 << "temecdqc = " << qc << endl;
        if (equil == 0)
        {
	       outyr = ((tem.endyr - tem.startyr + totsptime) / tem.diffyr) + 2;
	       for (dyr = 0; dyr < outyr; dyr++)
          {
	         ttotyr[dyr][itype] = tem.startyr - totsptime - 1
                                 + (dyr * tem.diffyr);
          }
        }
        else
        {
	       outyr = 1;
	       ttotyr[vt][itype] = -999;
        }

        if (mez < 0 || mez >= NUMVEG) { tem.veg.cmnt = 0; }

        // initstat[itype][0].ave changed to initstat[itype][tem.I_VEGC].ave by
        //    DWK on 20000214
        // tem.ez changed to tem.veg.cmnt by DWK on 20000214
        if ((tem.veg.cmnt > 0) && ((noprec == 1) || (lowtair == 1)
           || (stateflag == 1 && initstat[itype][tem.I_VEGC].ave < 0.1
           && initstat[itype][tem.I_VEGC].ave > MISSING)))
        {
	       if (equil == 1) { ttotyr[vt][itype] = 1; }
          for (dyr = 0; dyr < outyr; dyr++)
          {
            temwritemiss(ftemout,predmap,dyr,itype,ntempred,natmspred,ZERO);
 //	       tqc = 31; // to assign zero to all TEM variables during transient
          }
        }
        else
        {
          for (dyr = 0; dyr < outyr; dyr++)
          {
            temwritemiss(ftemout,predmap,dyr,itype,ntempred,natmspred,MISSING);
          }
        }
// bug fixed by Q. Zhuang sept 19 2000, change tem.ez to tem.veg.cmnt
	     if (tem.veg.cmnt < 0)
        {
	       if (clm.nirr[vt][0] >= 0.0 && tem.atms.tair[vt][0] > -99.0)
          {
            for (dyr = 0; dyr < outyr; dyr++)
            {
	           if (equil == 1) { ttotyr[dyr][itype] = tem.totyr; }
	           else
              {
                ttotyr[dyr][itype] = tem.startyr - totsptime - 1
                                     + (dyr * tem.diffyr);
              }
	           temwritepred(ftemout,predmap,dyr,itype,ntempred,natmspred,vt);
 	         }
	       }
	     }
      } // end of "qc = REJECT"

      else
      {

/* *******************************************************************
Start the Terrestrial Ecosystem Model (TEM) for Equilibrium Conditions
******************************************************************* */

 //       tem.soil.xtext(tem.veg.cmnt,tem.soil.pctsilt,tem.soil.pctclay);

        tem.soil.hydm_xtext(tem.veg.cmnt,tem.soil.pctsilt,tem.soil.pctclay,vt); // modified for two boxes hydrology model

// Initialize TEM parameters to grid cell's vegetation type, soil texture,
//   PET and AET

	     tem.setELMNTecd(kdinflg,tem.veg.cmnt, tem.soil.psiplusc,vt);                                                                               //wsr: set parameters for upland in equilibrium run
//        modified for hydrology model by qz adding tem.atms.vap);
//        tem.setELMNTevap(stateflag, tem.veg.cmnt, tem.atms.pet, tem.atms.tair);


      tem.setELMNTevap(stateflag, tem.veg.cmnt, tem.atms.pet, tem.atms.tair,tem.hyd.vap,vt);

// "While" loop to allow adaptive integrator tolerance (i.e. tem.tol) to be
// reduced if chaotic behavior occurs

        tem.ag.state = 0;
        tem.ag.prvstate = 0;

	     for (dyr = 0; dyr < RTIME; dyr++) { tem.qualcon[dyr][itype] = 0; }

	     tem.nattempt = 0;
	     tem.tol = tem.inittol;
	     tem.baseline = tem.initbase;

	    while (tem.nattempt < tem.maxnrun)
       {

// Initialize standing stocks of carbon and nitrogen from spatially-explicit
// data sets (i.e. stateflag = 1) or calibration ("ECD") data

	       if (stateflag == 1) { GISsetELMNTstate(tem.veg.cmnt, initstat,vt); }

		   

	       else {tem.ECDsetELMNTstate(tem.veg.cmnt, tem.soil.psiplusc,vt); }

     	       tem.setELMNTflux(vt);

// Run TEM until steady state conditions occur (equilibrium)

	        tem.nattempt = tem.equilibrium(itype,tem.tol,vt);

	       if (tem.nattempt < tem.maxnrun) { tem.tol /= 10.0; }
        }

// Output TEM steady state results to files

	     outyr = 0;
	     tem.qualcon[outyr][itype] += tem.nattempt;
	     if (tem.nattempt == tem.maxnrun && tem.totyr == tem.maxyears)
        {
	       tem.qualcon[outyr][itype] += 20;
	     }
	     if (equil == 0)
        {
          ttotyr[outyr][itype] = tem.startyr - totsptime - 1;
        }
	     else { ttotyr[outyr][itype] = tem.totyr; }


        temwritepred(ftemout,predmap,dyr,itype,ntempred,natmspred,vt);


 	     if (equil == 0)
        {

/* *******************************************************************
                     Begin Transient Conditions
******************************************************************* */


	       outyr = 1;
	       tem.microbe.kdsave[itype] = tem.microbe.kd[vt];
	       tem.totyr = tem.startyr - totsptime;
	       tqc = transqc(tem.maxyears,tem.totyr, tem.veg.plant[vt]);

	       tqc = ACCEPT;													//wsr:regardless of tqc value (negative carbon)
	       if (tqc == ACCEPT)
          {
	         tem.baseline = 0;
	         tem.totyr = 0;
	         tem.wrtyr = -99;
	         dyr = 1;
	         outyr = 1;
tem.soil.rhflag = 1;
int init_carbonflag;																					//wsr: set carbon flag to zero to input the initial soil carbon for transient step
init_carbonflag = 1;
tem.soil.callflag = 0;
tem.soil.sum_depth = 0;			//initial peat depth = upland depth = 1000 mm

																						       //transfer the peatland parameters into upland
//tem.setELMNTecd(kdinflg,tem.veg.cmnt, tem.soil.psiplusc,vt);                                                                      //set parameters for peatland





	         while (tem.totyr < tem.endyr)
            {
/*
// Calculate NIRR, CLDINESS and PAR for Spin up and Transient Conditions
//   with POTSCLM  -- also, hand off NIRR and PAR to TEM

	           if (atmsflag == 1 && tem.atms.tcldsflag == 1)
              {
                tem.totyr = clm.cldsyear[dyr];
		          for (dm = 0; dm < CYCLE; dm++)
                {
		            clm.clds[vt][dm] = tem.atms.tclds[dyr][dm];
		            if (cldflag == 1)
                  {
		              clm.nirr[vt][dm] = clm.xnirr(clm.clds[vt][dm], clm.girr[vt][dm]);
		            }
		            else
                  {
		              clm.clds[vt][dm] = clm.mkclds(clm.girr[vt][dm], clm.nirr[vt][dm]);
		            }
		           // tem.atms.par[vt][dm]  = clm.xpar(clm.clds[vt][dm], clm.nirr[vt][dm]);										wsr: all shut down
		           tem.atms.par[vt][dm] = clm.clds[vt][dm];
		            tem.atms.nirr[vt][dm] = clm.nirr[vt][dm];
		          }


// Output POTSCLM transient results to files

		          if (atmsoutfg == 1 && itype == 0)
                {
 	               atmswritepred(fsradout, predmap, dyr, natmspred);
		          }
	           }

*/
if(tem.totyr == 1000)     {           tem.transfer();   tem.setELMNTecd(kdinflg,tem.veg.cmnt, tem.soil.psiplusc,vt);    }         //wsr: run 1000 years as the upland ecosystem


      // Run the Terrestrial Ecosystem Model (TEM) under transient conditions
	           wrtyr = tem.transient(dyr,itype, tem.tol, RTIME,vt,init_carbonflag);
	           
	           init_carbonflag = 0;																						//wsr: reset carbonflag to zero 
                   
// Output TEM transient results for specified years to files

	           if ((wrtyr%tem.diffyr) == 0)
              {
		          ttotyr[outyr][itype] = tem.totyr;
//                temwritepred(ftemout,predmap,outyr,itype,ntempred,natmspred);
                temwritepred(ftemout,predmap,outyr,itype,ntempred,natmspred,vt);
		          ++outyr;
	           }
	           ++dyr;
	         } // End of While totyr < endyr
	         
	
	         
	       } // End of tqc = ACCEPT "if" statement
	       else
          {
            rflog1 << "tqc = " << tqc << endl;
	         outyr = ((tem.endyr - tem.startyr + totsptime) / tem.diffyr) + 2;
  	         ttotyr[outyr][itype] = tem.startyr - totsptime - 1
                                   + (dyr * tem.diffyr);
            if (tqc == 31)  // assign zero to all TEM variables
            {
              temwritemiss(ftemout,predmap,outyr,itype,
                           ntempred,natmspred,ZERO);
            }
            else
            {
              temwritemiss(ftemout,predmap,outyr,itype,
                           ntempred,natmspred,MISSING);
            }
	       } // End of tqc = REJECT
	     } // End of Transient else
      } // end of qc = ACCEPT "if" statement
    } // End of the Terrestrial Ecosystem Model (TEM)
  } // End of Vegetation Mosaic Loop

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */
//modified by QZ for hydrology model
int TEMelmnt::temgisin(ofstream& flog1, int& ftlerr, const int& atmsflag,
		       const int& itype, const int& maxtype,
                       const int& kdinflg, const int& stateflag,
                       FILE* fstxt, FILE*fsoiloggc, FILE*felev, FILE* fnirr, FILE* fpar,
		       FILE* ftair, FILE* fprec, FILE* fvap,FILE* flulc, FILE* fnpp,
                       FILE* fkdin, FILE* fstate[MAXESTAT+7],
		       const int& numspin, const int& spintime,
                       const int& RTIME)
{

  int i;
  int gisend = 1;
  int vt;

  Soildata fao;
  Elevdata elv;
  Tsoiloggcdata oggc;                    //wsr: add oggc
  Clmdata nirr;
  Clmdata par;
  Clmdata tair;
  Clmdata prec;
  Lulcdata lulc;
  Temdata npp;
  Temdata leaf;
  KDdata kddat;
 // added for hydrological model  by QZ
  Clmdata vap;

  if (itype == 0)
  {

    noprec = 0;
    lowtair = 0;

    gisend = fao.getdel(fstxt);
    if (gisend == -1) { return gisend; }
    ftlerr = coregerr(flog1,"TEMVEG", col, row, "TEXTURE", fao.col, fao.row);
    carea = fao.carea;
    tem.soil.pctsand = fao.pctsand; // added for hydrology model by QZ
    tem.soil.pctsilt = fao.pctsilt;
    tem.soil.pctclay = fao.pctclay;
    tem.soil.wsoil = fao.wsoil;

    gisend = elv.getdel(felev);
    if (gisend == -1) { return gisend; }
    ftlerr = coregerr(flog1,"TEMVEG", col, row, "ELEV", elv.col, elv.row);
    tem.elev = elv.elev;


if (tem.atms.tsoiloggcflag == 1)                                                                           //wsr: read soil organic carbon data
{
	gisend = oggc.getdel(fsoiloggc);
	if (gisend == -1) {return gisend;  }
	ftlerr = coregerr(flog1, "TEMVEG", col, row, "ORG", oggc.col, oggc.row);
	tem.orgcarbon = oggc.orgcarbon;  //wsr: add for soil organic carbon input

}



    if (atmsflag == 0)
    {
      gisend = nirr.getdel(fnirr);
      if (gisend == -1) { return gisend; }
      ftlerr = coregerr(flog1,"TEMVEG", col, row, "NIRR", nirr.col, nirr.row);
      for (i = 0; i < CYCLE; i++) { for (vt=0; vt< 2; vt++) tem.atms.nirr[vt][i] = nirr.mon[i]; }

    gisend = par.getdel(fpar);
     if (gisend == -1) { return gisend; }
    ftlerr = coregerr(flog1,"TEMVEG", col, row, "PAR", par.col, par.row);
    for (i = 0; i < CYCLE; i++) { for (vt=0; vt< 2; vt++) tem.atms.par[vt][i] = par.mon[i]; }
    }

    if (tem.atms.ttairflag == 1)
    {
      gisend = loadtetair(ftair, tair, numspin, spintime, RTIME, ftlerr, flog1);
      if (gisend == -1) { return gisend; }
      for (i = 0; i < CYCLE; i++) { for (vt=0; vt< 2; vt++) tem.atms.tair[vt][i] = tem.atms.ttair[vt][i]; }
    }
    else
    {
      gisend = tair.getdel(ftair);
      if (gisend == -1) { return gisend; }
      ftlerr = coregerr(flog1,"TEMVEG", col, row, "TAIR", tair.col, tair.row);
      tem.atms.mxtair = tair.max;

      // lowtair condition added by DWK on 20000214
      if (tem.atms.mxtair < -1.0) { lowtair = 1; }
      for (i = 0; i < CYCLE; i++) { for (vt=0; vt< 2; vt++) tem.atms.tair[vt][i] = tair.mon[i]; }
    }

    if (tem.atms.tprecflag == 1)
    {
      gisend = loadteprec(fprec, prec, numspin, spintime, RTIME, ftlerr, flog1);
      if (gisend == -1) { return gisend; }
      for (i = 0; i < CYCLE; i++) { for (vt=0; vt< 2; vt++) tem.atms.prec[vt][i] = tem.atms.tprec[vt][i]; }
    }
    else
    {
      gisend = prec.getdel(fprec);
      if (gisend == -1) { return gisend; }
      ftlerr = coregerr(flog1,"TEMVEG", col, row, "PREC", prec.col, prec.row);
      tem.atms.yrprec[vt] = prec.total;
      if (prec.total <= 0.0) { noprec = 1; }
      for (i = 0; i < CYCLE; i++) { for (vt=0; vt< 2; vt++) tem.atms.prec[vt][i] = prec.mon[i]; }
    }

   // added for hydrological model

    if (tem.hyd.vapflag == 1) {
      gisend = loadvap(fvap, vap, numspin, spintime, RTIME, ftlerr, flog1);
      if (gisend == -1) { return gisend; }
      for (i = 0; i < CYCLE; i++) { for (vt=0; vt< 2; vt++) tem.hyd.vap[vt][i] = tem.hyd.tvap[vt][i]; }
    }
    else {
      gisend = vap.getdel(fvap);
      if (gisend == -1) { return gisend; }
      ftlerr = coregerr(flog1,"TEMVEG", col, row, "VAP", vap.col, vap.row);
      tem.hyd.yrvap[vt] = vap.total;
      if (vap.total <= 0.0) { novap = 1; }
      for (i = 0; i < CYCLE; i++) { for (vt=0; vt< 2; vt++) tem.hyd.vap[vt][i] = vap.mon[i];}
          }

 // end of adding

    if (tem.ag.tlulcflag == 1)
    {
      gisend = loadtelulc(flulc, lulc, numspin, spintime, RTIME, ftlerr, flog1);
      if (gisend == -1) { return gisend; }
      tem.ag.state = tem.ag.tstate[vt];
      tem.ag.RAP = tem.ag.tRAP[vt];

// Note: numspin and spintime removed from function call by D. Kicklighter 990724

      if(tem.ag.RAP0flag == 0)  // only get potnpp if RAP > 0
      {
	gisend = loadtenpp(fnpp, npp, maxtype, RTIME, ftlerr, flog1);
      }
      if (gisend == -1) { return gisend; }
    }
  }

  if (kdinflg == 1)
  {
    kddat.getdel(fkdin);
    ftlerr = coregerr(flog1,"TEMVEG", col, row, "KD", kddat.col, kddat.row);
    tem.microbe.kdin[itype] = kddat.kd;
  }

/* Read in initialization data */

  if (stateflag == 1)
  {
    for (i = 0; i < MAXSTATE; i++)
    {
      gisend = initstat[itype][i].getdel(fstate[i]);
      if (gisend == -1) { return gisend; }
      ftlerr = coregerr(flog1,"TEMVEG", col, row, "Initial", initstat[itype][i].col, initstat[itype][i].row);
	  

      if (ftlerr != 0)
      {
	cout  << "Initial data set " << (i+1) << endl;
	flog1 << "Initial data set " << (i+1) << endl;
      }
    }
   for (vt=0; vt< 2; vt++) {
    tem.veg.prvleafmx[initstat[itype][MAXESTAT+4].subtveg-1] = initstat[itype][MAXESTAT+4].max * 0.01;
    tem.atms.yrpet[vt] = initstat[itype][MAXESTAT+5].total;
    tem.atms.prvpetmx = initstat[itype][MAXESTAT+5].max;
    tem.atms.yreet[vt] = initstat[itype][MAXESTAT+6].total;
    tem.atms.prveetmx = initstat[itype][MAXESTAT+6].max;
    }
    for (i = 0; i < CYCLE; i++)
    {
      for (vt=0; vt< 2; vt++) {
      tem.atms.pet[vt][i] = initstat[itype][MAXESTAT+5].mon[i];
      tem.veg.unnormleaf[vt][i] = initstat[itype][MAXESTAT+4].mon[i] * 0.01;
      }
    }
  }

//  for (i = 0; i < RTIME; i++) {
//    flog1 << tem.atms.tairyear[i];
//    for (int dm = 0; dm < CYCLE; dm++) {
//      flog1 << " " << setprecision(1) << tem.atms.ttair[i][dm];
//    }
//    flog1 << endl;
//}

  return gisend;

};

/* *************************************************************
************************************************************** */


/* * *************************************************************
************************************************************** */

int TEMelmnt::temgisqc(const int& stateflag, const double& pctsilt, const double& pctclay,
		       const int& cmnt, const double& elev, double nirr[CYCLE], double par[CYCLE],
		       double tair[CYCLE], double& mxtair, double prec[CYCLE], double& yrprec,
		       Temdata initstat[NUMMSAC][MAXSTATE+1])
{

  int i;
  int qc;

  qc = ACCEPT;

  if (pctsilt < 0.0) { return qc = 1; }
  if (pctclay < 0.0) { return qc = 2; }
// Following condition modified by DWK on 20000219
//  if (ez < 0 || ez >= NUMVEG) { return qc = 3; }
  if (cmnt < 1 || cmnt > NUMVEG) { return qc = 3; }
  if (elev <= -999.0) { return qc = 4;}
  if (mxtair < -1.0) { return qc = 5; }
  if (yrprec <= 0.0) { return qc = 6; }

  for (i = 0; i < CYCLE; i++)
  {
    if (nirr[i] <= -1.0) { return qc = 7; }
    if (par[i] <= -1.0) { return qc = 8; }
    if (tair[i] <= -99.0) { return qc = 9; }
    if (prec[i] <= -1.0) { return qc = 10; }
  }


  if (stateflag == 1)
  {
    for (i = 0; i < (MAXESTAT+5); i++)
    {
      if (initstat[0][i].mon[11] <= MISSING) { return qc = (11 + i); }
    }
    if (initstat[0][0].ave < 0.1) { return qc = 17 + MAXESTAT; }
  }

  return qc;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TEMelmnt::temwritemiss(ofstream fout[MAXPRED], char predname[MAXPRED][9],
                            const int& dyr, const int& itype,
                            const int& natmspred, const int& ntempred,
                            const double value)
{

  int i;
  int k;
  int dm;
  Temdata tempred;

  for (i = 0; i < ntempred; i++)
  {
    k = i + natmspred;

    for (dm = 0; dm < CYCLE; dm++)
    {
   	tempred.mon[dm] = value;
    }

    // Write output data to files

    tempred.outdel(fout[i], col, row, predname[k],
                   tem.veg.temveg, tem.veg.subtype[mez][itype],
                   (100.0 * tem.soil.psiplusc),
                   tem.qualcon[dyr][itype],
                   carea,
                   ttotyr[dyr][itype],
                   tempred.mon,
                   contnent);

  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TEMelmnt::temwritepred(ofstream fout[MAXPRED], char predname[MAXPRED][9],
                            const int& dyr, const int& itype,
                            const int& ntempred, const int& natmspred, const int& vt)
{

  int i;
  int k;
  int dm;
  Temdata tempred;

  for (i = 0; i < ntempred; i++)
  {
    k = i + natmspred;

    if (strcmp(predname[k],tem.predstr[tem.I_NPP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
    	  tempred.mon[dm] = tem.veg.npp[vt][dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_VEGC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.plant[vt][dm].carbon;  // VEGC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SOLC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.org[vt][dm].carbon;  // SOLC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_STRN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.strctrl[vt][dm].nitrogen;  // VSTRUCTN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_STON]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
  	     tempred.mon[dm] = tem.veg.labile[vt][dm].nitrogen * 1000.0;  // VSTOREN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SOLN]) == 0)
    {
    for (dm = 0; dm < CYCLE; dm++)
      {
     	  tempred.mon[dm] = tem.soil.org[vt][dm].nitrogen;  // SOILORGN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AVLN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.availn[vt][dm] * 1000.0;  // AVAILN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_NMIN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.microbe.netnmin[vt][dm] * 1000.0;  // NETNMIN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_NEP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.nep[vt][dm];  // NEP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_GPP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.gpp[vt][dm];  // GPP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_INGPP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.ingpp[vt][dm];  // INGPP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_INNPP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.innpp[vt][dm];  // INNPP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_RH]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.microbe.rh[vt][dm];  // RH
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_NLST]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.nlost[vt][dm] * 1000.0;  // NLOST
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_NINP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.ninput[vt][dm] * 1000.0;  // NINPUT
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_RVMNT]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.rm[vt][dm];  // RVMAINT
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_RVGRW]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.rg[vt][dm];  // RVGRWTH
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_GPR]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.gpr[vt][dm];  // GPR
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_LTRC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.ltrfal[vt][dm].carbon;  // LTRC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_VNUP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.nuptake[vt][dm] * 1000.0;  // VEGNUP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_LTRN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.ltrfal[vt][dm].nitrogen * 1000.0;  // LTRN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_MNUP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.microbe.nuptake[vt][dm] * 1000.0;  // MICRONUP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_VNMBL]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.nmobil[vt][dm] * 1000.0;  // VNMOBIL
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_VNRSRB]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.nresorb[vt][dm] * 1000.0;  // VNRESORB
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_VSUP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.suptake[vt][dm] * 1000.0;  // VEGSUP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_VLUP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.luptake[vt][dm] * 1000.0;  // VEGLUP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_UNRMLF]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.unnormleaf[vt][dm] * 100.0;  // UNRMLEAF
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_LEAF]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.leaf[vt][dm] * 100.0;      // LEAF
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_CNVRTC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.convrtflx[vt][dm].carbon;   // CONVERTC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_CNVRTN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.convrtflx[vt][dm].nitrogen * 1000.0; // CONVERTN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRDF10C]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.formPROD10[vt].carbon / (double) CYCLE; // PRDF10C
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRDF10N]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.formPROD10[vt].nitrogen / (double) CYCLE) *
                          1000.0; // PRDF10N
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRDF100C]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.formPROD100[vt].carbon / (double) CYCLE; // PRDF100C
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRDF100N]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.formPROD100[vt].nitrogen / (double) CYCLE) *
                          1000.0; // PRDF100N
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PROD10C]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.PROD10[vt].carbon; // PROD10C
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PROD10N]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.PROD10[vt].nitrogen * 1000.0; // PROD10N
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PROD100C]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.PROD100[vt].carbon; // PROD100C
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PROD100N]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.PROD100[vt].nitrogen * 1000.0; // PROD100N
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRD10FC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.PROD10decay[vt].carbon / (double) CYCLE; // PRD10FC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRD10FN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.PROD10decay[vt].nitrogen / (double) CYCLE) *
                          1000.0; // PRD10FN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRD100FC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.PROD100decay[vt].carbon / (double) CYCLE; // PRD100FC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRD100FN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.PROD100decay[vt].nitrogen / (double) CYCLE) *
                          1000.0; // PRD100FN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGNPPC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.npp[vt][dm].carbon; // AGNPPC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGNPPN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.npp[vt][dm].nitrogen * 1000.0; // AGNPPN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGFPRDC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.formPROD1[vt].carbon / (double) CYCLE; // AGFPRODC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGFPRDN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.formPROD1[vt].nitrogen / (double) CYCLE) *
                          1000.0; // AGFPRODN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGFRTN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.fertn[vt][dm] * 1000.0; // AGFERTN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGLTRC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.ltrfal[vt][dm].carbon; // AGLTRFC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGLTRN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.ltrfal[vt][dm].nitrogen * 1000.0; // AGLTRFN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SLASHC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.slash[vt][dm].carbon; // SLASHC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SLASHN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.slash[vt][dm].nitrogen; // SLASHN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGPRDC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.PROD1[vt].carbon; // AGPRODC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGPRDN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.PROD1[vt].nitrogen * 1000.0; // AGPRODN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGPRDFC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.PROD1decay[vt].carbon / (double) CYCLE; // AGPRODFC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGPRDFN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.PROD1decay[vt].nitrogen / (double) CYCLE) *
                          1000.0; // AGPRODFN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTFPRDC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.formTOTPROD[vt].carbon / (double) CYCLE; // TOTFPRDC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTFPRDN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.formTOTPROD[vt].nitrogen / (double) CYCLE) *
                          1000.0; // TOTFPRDN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTPRDC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.TOTPROD[vt].carbon; // TOTPRODC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTPRDN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.TOTPROD[vt].nitrogen * 1000.0; // TOTPRODN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTPRDFC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.TOTPRODdecay[vt].carbon / (double) CYCLE; // TOTPRDFC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTPRDFN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.TOTPRODdecay[vt].nitrogen / (double) CYCLE) *
                          1000.0; // TOTPRDFN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTEC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.totalc[vt][dm]; // TOTEC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTNPP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.totnpp[vt][dm]; // TOTNPP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTGC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.totalc[vt][dm] + tem.ag.TOTPROD[vt].carbon; // TOTGC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_CFLX]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.cflux[vt][dm]; // CFLUX
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SCNVRTC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.sconvrtflx[vt][dm].carbon; // SCONVRTC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SCNVRTN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.sconvrtflx[vt][dm].nitrogen; // SCONVRTN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_NSRTNT]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.nsretent[vt][dm]; // NSRETENT
      }
    }
     else if (strcmp(predname[k],tem.predstr[tem.I_VSM]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.vsm[vt][dm] * 100.0;      // VSM
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PET]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.atms.pet[vt][dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_EET]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.atms.eet[vt][dm];              // EET
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_RAIN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.atms.rain[vt][dm];             // RAIN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SNWINF]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.snowinf[vt][dm];          // SNOWINF
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_WYLD]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.h2oyld[vt][dm];          // H2OYIELD
      }
    }

  else if (strcmp(predname[k],tem.predstr[tem.I_PCTP1]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.pctp1[vt][dm];          // PCTP1
      }
    }

  else if (strcmp(predname[k],tem.predstr[tem.I_PCTP2]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.pctp2[vt][dm] ;          // PCTP2
      }
    }

  else if (strcmp(predname[k],tem.predstr[tem.I_PCTP3]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.pctp3[vt][dm] ;          // PCTP3
      }
    }

  else if (strcmp(predname[k],tem.predstr[tem.I_VSM1]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	 //tempred.mon[dm] = tem.soil.vsm1[vt][dm] * 100.0; 
	    tempred.mon[dm] =( tem.soil.avlh2o1[vt][dm] + tem.soil.wiltpt1[vt] ) * 100 / (tem.soil.dpwbox1[vt]*1000.0) ;      // VSM1
      }
    }
  else if (strcmp(predname[k],tem.predstr[tem.I_VSM2]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	   // tempred.mon[dm] = tem.soil.vsm2[vt][dm] * 100.0; 
	     tempred.mon[dm] =( tem.soil.avlh2o2[vt][dm] + tem.soil.wiltpt2[vt] ) * 100.0 / (tem.soil.dpwbox2[vt]*1000.0);      // VSM2
      }
    }
  else if (strcmp(predname[k],tem.predstr[tem.I_VSM3]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
       tempred.mon[dm] = tem.soil.vsm3[vt][dm] * 100.0; 
	     tempred.mon[dm] = ( tem.soil.avlh2o3[vt][dm] + tem.soil.wiltpt3[vt] ) * 100.0 / (tem.soil.dpwbox3[vt]*1000.0);      // VSM3
      }
    }

  // added for hydrology model
    else if (strcmp(predname[k], tem.predstr[tem.I_EVAP]) == 0) {
      for (dm = 0; dm < CYCLE; dm++) {
	tempred.mon[dm] = tem.hyd.hevap[vt][dm];          // Evaporation
         }
    }
    else if (strcmp(predname[k], tem.predstr[tem.I_TRANS]) == 0) {
      for (dm = 0; dm < CYCLE; dm++) {
	tempred.mon[dm] = tem.hyd.htrans[vt][dm];          // Transpiration
         }
    }
    else if (strcmp(predname[k], tem.predstr[tem.I_SEVAP]) == 0) {
      for (dm = 0; dm < CYCLE; dm++) {
	tempred.mon[dm] = tem.hyd.soil_evap[vt][dm];          // soil surface evaporation
         }
    }
    else if (strcmp(predname[k], tem.predstr[tem.I_SNOWSUB]) == 0) {
      for (dm = 0; dm < CYCLE; dm++) {
	tempred.mon[dm] = tem.hyd.snowsub[vt][dm];          // snow sublimation
         }
    }
    else if (strcmp(predname[k], tem.predstr[tem.I_SUBCAN]) == 0) {
      for (dm = 0; dm < CYCLE; dm++) {
	tempred.mon[dm] = tem.hyd.sub_from_canopy[vt][dm];          // snow sublimation from canopy
         }
    }

 // end of adding

  // for soil temperature
    else if (strcmp(predname[k],tem.predstr[tem.I_TSOIL]) == 0) {
      for (dm = 0; dm < CYCLE; dm++) {
	tempred.mon[dm] = tem.atms.tsoil[vt][dm];          // Tsoil
       }
     }
    else if (strcmp(predname[k],tem.predstr[tem.I_DST5]) == 0) {
      for (dm = 0; dm < CYCLE; dm++) {
	tempred.mon[dm] = tem.atms.dst5[vt][dm];          // Tsoil
       }
     }
  else if (strcmp(predname[k],tem.predstr[tem.I_DST10]) == 0) {
      for (dm = 0; dm < CYCLE; dm++) {
	tempred.mon[dm] = tem.atms.dst10[vt][dm];          // Tsoil
       }
     }
  else if (strcmp(predname[k],tem.predstr[tem.I_DST20]) == 0) {
      for (dm = 0; dm < CYCLE; dm++) {
	tempred.mon[dm] = tem.atms.dst20[vt][dm];          // Tsoil
       }
     }
  else if (strcmp(predname[k],tem.predstr[tem.I_DST50]) == 0) {
      for (dm = 0; dm < CYCLE; dm++) {
	tempred.mon[dm] = tem.atms.dst50[vt][dm];          // Tsoil
       }
     }
  else if (strcmp(predname[k],tem.predstr[tem.I_DST100]) == 0) {
      for (dm = 0; dm < CYCLE; dm++) {
	tempred.mon[dm] = tem.atms.dst100[vt][dm];          // Tsoil
       }
     }
  else if (strcmp(predname[k],tem.predstr[tem.I_DST200]) == 0) {
      for (dm = 0; dm < CYCLE; dm++) {
	tempred.mon[dm] = tem.atms.dst200[vt][dm];          // Tsoil
       }
     }
  else if (strcmp(predname[k],tem.predstr[tem.I_FRONTD]) == 0) {
      for (dm = 0; dm < CYCLE; dm++) {
	tempred.mon[dm] = tem.atms.frontd[vt][dm];          // Tsoil
       }
     }
  else if (strcmp(predname[k],tem.predstr[tem.I_THAWBE]) == 0) {
      for (dm = 0; dm < CYCLE; dm++) {
	tempred.mon[dm] = tem.atms.thawbe[vt][dm];          // Tsoil
       }
     }
  else if (strcmp(predname[k],tem.predstr[tem.I_THAWEND]) == 0) {
      for (dm = 0; dm < CYCLE; dm++) {
	tempred.mon[dm] = tem.atms.thawend[vt][dm];          // Tsoil
       }
     }
      // end of soil temperature

    else if (strcmp(predname[k],tem.predstr[tem.I_PCTP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.pctp[vt][dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SM1]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.avlh2o1[vt][dm] + tem.soil.wiltpt1[vt];                     //SM1
      }
    }
	  else if (strcmp(predname[k],tem.predstr[tem.I_SM2]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.avlh2o2[vt][dm] + tem.soil.wiltpt2[vt];                     //SM2
      }
    }
	  else if (strcmp(predname[k],tem.predstr[tem.I_SM3]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.avlh2o3[vt][dm] + tem.soil.wiltpt3[vt];                      //SM3
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AVLW]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.avlh2o[vt][dm];                              //AVWL
      }
    }
	    else if (strcmp(predname[k],tem.predstr[tem.I_AVLW1]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.avlh2o1[vt][dm];                          //AVWL1
      }
    }
	    else if (strcmp(predname[k],tem.predstr[tem.I_AVLW2]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.avlh2o2[vt][dm];                       //AVLW2
      }
    }
		    else if (strcmp(predname[k],tem.predstr[tem.I_AVLW3]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.avlh2o3[vt][dm];                       //AVLW3
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_RRUN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.rrun[vt][dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SNWFAL]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.atms.snowfall[vt][dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SNWPCK]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.snowpack[vt][dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SRUN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.srun[vt][dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_RGRW]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.rgrndh2o[vt][dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SGRW]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.sgrndh2o[vt][dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_RPERC1]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.rperc1[vt][dm];                          //rperc1
      }
    }
	    else if (strcmp(predname[k],tem.predstr[tem.I_RPERC2]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.rperc2[vt][dm];                          //rperc2
      }
    }
	    else if (strcmp(predname[k],tem.predstr[tem.I_RPERC3]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.rperc3[vt][dm];                          //rperc3
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SPERC1]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.sperc1[vt][dm];                            //sperc1
      }
    }
		    else if (strcmp(predname[k],tem.predstr[tem.I_SPERC2]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.sperc2[vt][dm];                            //sperc2
      }
    }
	    else if (strcmp(predname[k],tem.predstr[tem.I_SPERC3]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.sperc3[vt][dm];                            //sperc3
      }
    }
				    else if (strcmp(predname[k],tem.predstr[tem.I_LYPERC1]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.hydm.lyperc1[vt][dm];                            //lyperc1
      }
    }
			    else if (strcmp(predname[k],tem.predstr[tem.I_LYPERC2]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.hydm.lyperc2[vt][dm];                            //lyperc2
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.totalc[vt][dm];  // TOTALC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_VEGN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
 	     tempred.mon[dm] = tem.veg.plant[vt][dm].nitrogen;  // VEGN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_LAI]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
 	     tempred.mon[dm] = tem.veg.lai[vt][dm];  // LAI
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_FPC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
 	     tempred.mon[dm] = tem.veg.fpc[vt][dm] * 100.0;  // FPC
      }
    }

    // Write output data to files

    if (strcmp(predname[k],tem.predstr[tem.I_VSM]) == 0
        || strcmp(predname[k],tem.predstr[tem.I_PCTP]) == 0
        || strcmp(predname[k],tem.predstr[tem.I_LEAF]) == 0)
    {
// addition for moss and black spruce
     int temp_veg;
     if (vt ==0)  { temp_veg=4;}
     else  { temp_veg = 5;};

      tempred.poutdel(fout[i], col, row, predname[k],
                     temp_veg, temp_veg,
                     (100.0 * tem.soil.psiplusc),
                     tem.qualcon[dyr][itype],
                     carea,
                     ttotyr[dyr][itype],
                     tempred.mon,
                     contnent);
    }
    else
    {
     int temp_veg;
     if (vt ==0)  { temp_veg=4;}
     else  { temp_veg = 5;};

      tempred.outdel(fout[i], col, row, predname[k],
                     temp_veg, temp_veg,
//                   tem.veg.temveg, tem.veg.cmnt,

                     (100.0 * tem.soil.psiplusc),
                     tem.qualcon[dyr][itype],
                     carea,
                     ttotyr[dyr][itype],
                     tempred.mon,
                     contnent);
    }
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

int TEMelmnt::transqc(int& maxyears, long& totyr, Biomass plant[CYCLE])
{

  int i;
  int qc;
  double sumcarbon = 0.0;
  qc = ACCEPT;

  if (totyr < 0 || totyr >= (long) maxyears) { return qc = 30; }
  for (i = 0; i < CYCLE; i++) { sumcarbon += plant[i].carbon; }
  if (sumcarbon <= 0.0) { return qc = 31; }

  return qc;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

int TEMelmnt::veggisin(ofstream& flog1, int& ftlerr, const int& atmsflag,
		       const int& temflag, FILE* ftveg)
{

  int gisend;

  Vegdata tveg;

  gisend = tveg.getdel(ftveg);
  if (gisend == -1) { return gisend; }

  if (atmsflag == 0)
  {
    col = tveg.col;
    row = tveg.row;
  }
  else
  {
    ftlerr = coregerr(flog1,"CLOUDS", col, row, "TEMVEG", tveg.col, tveg.row);
  }

  strcpy(contnent, tveg.contnent);
  mez = tveg.temveg - 1;
  if (temflag == 1)
  {
    maxtype = tem.veg.numtype[mez];
    tem.veg.temveg = tveg.temveg;
  }
  if (mez < 0 || mez >= NUMVEG)
  {
    mez = 0;  // Added by DWK on 20000219
    maxtype = 1;
  }

  return gisend;

};

