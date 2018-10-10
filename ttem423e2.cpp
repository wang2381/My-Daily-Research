/* *************************************************************
****************************************************************
TTEM423.CPP - Terrestrial Ecosystem Model Version 4.2
****************************************************************

Modifications:

19991028 - DWK add bug fixes
20000107 - DWK adds compiler directive
20000107 - DWK renames functions to better describe purpose
20000130 - DWK changes prevleafmx to initleafmx in fourth getsiteecd()
20000130 - DWK changes y[2] to y[I_SOLC] and y[3] to y[I_SOLN]
           at bottom of stepyr()
20000201 - DWK changes pstate[1] to pstate[I_STRN] in delta()
20000207 - DWK added a minimum VEGC and SOLC (0.00001 gC m-2) condition
           to massbal(); Commented out by DWK on 20000210
20000616 - CFD removed Rh and Nmin from boundcon(), added RvMaint to boundcon()
20000616 - CFD restructured massbal()
20000715 -QZ modified tranisentyear () for hydrology model 
****************************************************************
************************************************************** */

#if !defined(TTEM423_H)
  #include "ttem423e.hpp"
#endif

/* *********************************************************** */

TTEM::TTEM() : Odeint4()
{

  nfert = 1.00000;
  tol = inittol;
  totyr = -99;

// Identify potential output variables from TEM

// Ecosystem carbon pools ***************************************

  strcpy(predstr[I_VEGC],"VEGC");       // vegetation carbon
  strcpy(predstr[I_SOLC],"SOILORGC");   // soil organic carbon
  strcpy(predstr[I_TOTC],"TOTALC");     // total carbon


// Ecosystem nitrogen pools *************************************

  // total nitrogen stored in vegetation
  strcpy(predstr[I_VEGN],"VEGN");

  // vegetation structural nitrogen
  strcpy(predstr[I_STRN],"VSTRUCTN");

  // vegetation labile nitrogen
  strcpy(predstr[I_STON],"VSTOREN");

  strcpy(predstr[I_SOLN],"SOILORGN");  // soil organic nitrogen
  strcpy(predstr[I_AVLN],"AVAILN");    // soil available nitrogen


// Carbon and nitrogen pools associated with human products *****

  // carbon in agricultural products
  strcpy(predstr[I_AGPRDC],"AGPRODC");

  // nitrogen in agricultural products
  strcpy(predstr[I_AGPRDN],"AGPRODN");

  // carbon pool of products that decompose in 10 years
  strcpy(predstr[I_PROD10C],"PROD10C");

  // nitrogen pool of products that decompose in 10 years
  strcpy(predstr[I_PROD10N],"PROD10N");

  // carbon pool of products that decompose in 100 years
  strcpy(predstr[I_PROD100C],"PROD100C");

  // nitrogen pool of products that decompose in 100 years
  strcpy(predstr[I_PROD100N],"PROD100N");

  // carbon in all product pools
  strcpy(predstr[I_TOTPRDC],"TOTPRODC");

  // nitrogen in all product pools
  strcpy(predstr[I_TOTPRDN],"TOTPRODN");

  // total carbon pool found in ecosystem excluding products
  strcpy(predstr[I_TOTEC],"TOTEC");

  // total carbon pool found in ecosystem including products
  strcpy(predstr[I_TOTGC],"TOTGC");


// Ecosystem water pools ****************************************

  // available soil moisture
  strcpy(predstr[I_AVLW],"AVAILH2O");

  // groundwater pool resulting from rainfall
  strcpy(predstr[I_RGRW],"RGRNDH2O");

  strcpy(predstr[I_SNWPCK],"SNOWPACK");  // snowpack

  // groundwater pool resulting from snow melt
  strcpy(predstr[I_SGRW],"SGRNDH2O");

  strcpy(predstr[I_SM],"SOILH2O");       // soil moisture

  // soil moisture expressed as percent total porosity
  strcpy(predstr[I_PCTP],"PCTP");

  strcpy(predstr[I_VSM],"VSM");      // volumetric soil moisture

 // *** the following are pools for 3-box hydrology

  strcpy(predstr[I_AVLW1], "AVAIL1");
  strcpy(predstr[I_AVLW2], "AVAIL2");
  strcpy(predstr[I_AVLW3], "AVAIL3");

  strcpy(predstr[I_SM1], "SOILH21");
  strcpy(predstr[I_SM2], "SOILH22");
  strcpy(predstr[I_SM3], "SOILH23");

  strcpy(predstr[I_PCTP1],"PCTP1");
  strcpy(predstr[I_PCTP2],"PCTP2");
  strcpy(predstr[I_PCTP3],"PCTP3");

  strcpy(predstr[I_VSM1],"VSM1");
  strcpy(predstr[I_VSM2],"VSM2");
  strcpy(predstr[I_VSM3],"VSM3");

// end of adding

// Carbon fluxes for natural ecosystems *************************

  // GPP not limited by nutrient availability
  strcpy(predstr[I_INGPP],"VEGINGPP");

  strcpy(predstr[I_GPP],"GPP");      // gross primary production

  // NPP not limited by nutrient availability
  strcpy(predstr[I_INNPP],"VEGINNPP");

  strcpy(predstr[I_NPP],"NPP");      // net primary production
  strcpy(predstr[I_GPR],"GPR");      // gross plant respiration

  // vegetation maintenance respiration
  strcpy(predstr[I_RVMNT],"RVMAINT");

  // vegetation growth respiration
  strcpy(predstr[I_RVGRW],"RVGRWTH");

  strcpy(predstr[I_LTRC],"LTRC");   // litterfall carbon
  strcpy(predstr[I_RH],"RH");       // heterotrophic respiration
  strcpy(predstr[I_NEP],"NEP");     // net ecosystem production


// Nitrogen fluxes for natural ecosystems ***********************

  // total nitrogen inputs into ecosystem
  strcpy(predstr[I_NINP],"NINPUT");

  // VEGNUP not limited by carbon availability
  strcpy(predstr[I_INNUP],"VEGINNUP");

  // nitrogen uptake by vegetation
  strcpy(predstr[I_VNUP],"VEGNUP");

  // vegetation nitrogen uptake for structural components
  strcpy(predstr[I_VSUP],"VEGSUP");

  // vegetation nitrogen uptake for labile components
  strcpy(predstr[I_VLUP],"VEGLUP");

  // nitrogen mobilization by vegetation
  strcpy(predstr[I_VNMBL],"VNMOBIL");

  strcpy(predstr[I_VNRSRB],"VNRESORB"); // nitrogen resorption by vegetation
  strcpy(predstr[I_LTRN],"LTRN");       // litterfall nitrogen
  strcpy(predstr[I_MNUP],"MICRONUP");   // nitrogen uptake by microbes
  strcpy(predstr[I_NMIN],"NETNMIN");    // net nitrogen mineralization
  strcpy(predstr[I_NLST],"NLOST");      // total nitrogen losses from ecosystem


// Water fluxes *************************************************

  strcpy(predstr[I_RAIN],"RAIN");        // rainfall

  // percolation of rainwater through soil profile
  strcpy(predstr[I_RPERC],"RPERC");

  strcpy(predstr[I_RRUN],"RRUN");        // runoff of rainwater
  strcpy(predstr[I_SNWFAL],"SNOWFALL");  // snowfall

  // infiltration into the soil of water from snowmelt
  strcpy(predstr[I_SNWINF],"SNOWINF");

  // percolation of snowmelt through soil profile
  strcpy(predstr[I_SPERC],"SPERC");

  strcpy(predstr[I_SRUN],"SRUN");        // runoff of snowmelt
  strcpy(predstr[I_PET],"PET");          // potential evapotranspiration
  strcpy(predstr[I_EET],"EET");          // estimated evapotranspiration
  strcpy(predstr[I_WYLD],"H2OYIELD");    // water yield

  // added for hydrology model by QZ
  strcpy(predstr[I_EVAP],"EVAP");    // evaporation from canopy
  strcpy(predstr[I_TRANS],"TRANS");    // transpiration of canopy
  strcpy(predstr[I_SEVAP],"SEVAP"); // soil surface evaporation
  strcpy(predstr[I_SNOWSUB],"SNSUB"); // snow sublimation from ground
  strcpy(predstr[I_SUBCAN],"SUBCAN"); // snow sublimation from canopy

// added on 15/NOV/98 for soil temperature
   strcpy(predstr[I_TSOIL],"TSOIL");
   strcpy(predstr[I_DST5],"DST5");
   strcpy(predstr[I_DST10],"DST10");
   strcpy(predstr[I_DST20],"DST20");
   strcpy(predstr[I_DST50],"DST50");
   strcpy(predstr[I_DST100],"DST100");
   strcpy(predstr[I_DST200],"DST200");
   strcpy(predstr[I_FRONTD],"FRONTD");
   strcpy(predstr[I_THAWBE],"THAWBE");
   strcpy(predstr[I_THAWEND],"THAWEND");
// end of ...


//********* added for soil hydrology model

  strcpy(predstr[I_LYPERC], "LYPERC");
  strcpy(predstr[I_LYPERC1],"LYPERC1"); // from moss layer to orgainic
  strcpy(predstr[I_LYPERC2],"LYPERC2"); // from organic to mineral soil layer

  strcpy(predstr[I_RPERC1],"RPERC1");
  strcpy(predstr[I_RPERC2],"RPERC2");
  strcpy(predstr[I_RPERC3],"RPERC3");

  strcpy(predstr[I_SPERC1],"SPERC1");
  strcpy(predstr[I_SPERC2],"SPERC2");
  strcpy(predstr[I_SPERC3],"SPERC3");


// Phenology variables for natural ecosystems *******************

  // un-normalized relative phenology variable
  strcpy(predstr[I_UNRMLF],"UNRMLEAF");

  // normalized relative phenology variable (0 - 1.0)
  strcpy(predstr[I_LEAF],"LEAF");

  strcpy(predstr[I_FPC],"FPC");          // foliar projected cover
  strcpy(predstr[I_LAI],"LAI");          // leaf area index


// Carbon and nitrogen fluxes associated with agricultural
// conversion ***************************************************

  strcpy(predstr[I_CNVRTC],"CONVERTC"); // carbon loss during conversion
  strcpy(predstr[I_CNVRTN],"CONVERTN"); // nitrogen loss during conversion
  strcpy(predstr[I_SCNVRTC],"SCONVRTC");
  strcpy(predstr[I_SCNVRTN],"SCONVRTN");
  strcpy(predstr[I_NVRTNT],"NVRETENT");
  strcpy(predstr[I_NSRTNT],"NSRETENT");
  strcpy(predstr[I_NRETNT],"NRETENT");

  // carbon associated with slash left after conversion
  strcpy(predstr[I_SLASHC],"SLASHC");

  // nitrogen associated with slash left after conversion
  strcpy(predstr[I_SLASHN],"SLASHN");

  // carbon loss to formation of products that decompose in 10 years
  strcpy(predstr[I_PRDF10C],"PRDF10C");

  // nitrogen loss to formation of products that decompose in 10 years
  strcpy(predstr[I_PRDF10N],"PRDF10N");

  // carbon loss to formation of products that decompose in 100 years
  strcpy(predstr[I_PRDF100C],"PRDF100C");

  // nitrogen loss to formation of products that decompose in 100 years
  strcpy(predstr[I_PRDF100N],"PRDF100N");


// Carbon and nitrogen fluxes associated with agriculture *******

  // agricultral net primary production
  strcpy(predstr[I_AGNPPC],"AGNPPC");

  // nitrogen uptake by crops associated with AGNPPC
  strcpy(predstr[I_AGNPPN],"AGNPPN");

  // carbon loss to formation of agricultural products
  strcpy(predstr[I_AGFPRDC],"AGFPRODC");

  // nitrogen loss to formation of agricultural products
  strcpy(predstr[I_AGPRDN],"AGFPRODN");

  strcpy(predstr[I_AGFRTN],"AGFERTN"); // nitrogen fertilization
  strcpy(predstr[I_AGLTRC],"AGLTRFC"); // litterfall carbon from crops
  strcpy(predstr[I_AGLTRN],"AGLTRFN"); // litterfall nitrogen from crops

// Carbon and nitrogen fluxes associated with human products ****

  // carbon loss to the formation of all products

  strcpy(predstr[I_TOTFPRDC],"TOTFPRDC");

  // nitrogen loss to the formation of all products

  strcpy(predstr[I_TOTFPRDN],"TOTFPRDN");

  // carbon loss to resulting from decomposition of agricultural
  //   products

  strcpy(predstr[I_AGPRDFC],"AGPRODFC");

  // nitrogen loss resulting from decomposition of agricultural
  //   products

  strcpy(predstr[I_AGPRDFN],"AGPRODFN");

  // carbon loss resulting from decomposition of PROD10C

  strcpy(predstr[I_PRD10FC],"PRD10FC");

  // nitrogen loss resulting from decomposition of PROD10N

  strcpy(predstr[I_PRD10FN],"PRD10FN");

  // carbon loss resulting from decomposition of PROD100C

  strcpy(predstr[I_PRD100FC],"PRD100FC");

  // nitrogen loss resulting from decomposition of PROD100N

  strcpy(predstr[I_PRD100FN],"PRD100FN");

  // carbon loss resulting from decomposition of all products

  strcpy(predstr[I_TOTPRDFC],"TOTPRDFC");

  // nitrogen loss resulting from decomposition of all products

  strcpy(predstr[I_TOTPRDFN],"TOTPRDFN");


// Integrated carbon fluxes *************************************

  // total net primary production (NPP+AGNPPC)

  strcpy(predstr[I_TOTNPP],"TOTNPP");

  // carbon flux from ecosystem (NEP+CONVERTC)

  strcpy(predstr[I_CFLX],"CFLUX");

};

/* **************************************************************
************************* Public Functions **********************
************************************************************** */


/* **************************************************************
************************************************************** */
//Rh and Nmin removed by CFD 20000616

//int TTEM::boundcon(double ptstate[], double err[], double& ptol)
int TTEM::boundcon(double ptstate[NVT][NUMEQ], double err[NVT][NUMEQ], double& ptol, const int& vt)
{

  int test = ACCEPT;

// Check carbon and nitrogen state variables
  if (err[vt][I_VEGC] > fabs(ptol * ptstate[vt][I_VEGC]))
  {
    return test = temkey(I_VEGC)+1;
  }
  if (nfeed == 1 && err[vt][I_STRN] > fabs(ptol * ptstate[vt][I_STRN]))
  {
    return test = temkey(I_STRN)+1;
  }
  if (err[vt][I_SOLC] > fabs(ptol * ptstate[vt][I_SOLC]))
  {
    return test = temkey(I_SOLC)+1;
  }
  if (nfeed == 1 && err[vt][I_SOLN] > fabs(ptol * ptstate[vt][I_SOLN]))
  {
    return test = temkey(I_SOLN)+1;
  }
  if (nfeed == 1 && err[vt][I_AVLN] > fabs(ptol * ptstate[vt][I_AVLN]))
  {
    return test = temkey(I_AVLN)+1;
  }
  if (err[vt][I_GPP] > fabs(ptol * ptstate[vt][I_GPP]))
  {
    return test = temkey(I_GPP)+1;
  }
  if (err[vt][I_NPP] > fabs(ptol * ptstate[vt][I_NPP]))
  {
    return test = temkey(I_NPP)+1;
  }
  if (nfeed == 1 && err[vt][I_VNUP] > fabs(ptol * ptstate[vt][I_VNUP]))
  {
    return test = temkey(I_VNUP)+1;
  }
  if (nfeed == 1 && err[vt][I_VSUP] > fabs(ptol * ptstate[vt][I_VSUP]))
  {
    return test = temkey(I_VSUP)+1;
  }
  if (nfeed == 1 && err[vt][I_STON] > fabs(ptol * ptstate[vt][I_STON]))
  {
    return test = temkey(I_STON)+1;
  }
  if (nfeed == 1 && err[vt][I_VNMBL] > fabs(ptol * ptstate[vt][I_VNMBL]))
  {
    return test = temkey(I_VNMBL)+1;
  }
  //20000616 - CFD added RvMaint
  if (nfeed == 1 && err[vt][I_RVMNT] > fabs(ptol * ptstate[vt][I_RVMNT]))
  {
    return test = temkey(I_RVMNT)+1;
  }
  //20000616 - CFD removed Rh and Nmin

  // Check water state variables

  if (err[vt][I_AVLW]  > fabs(ptol * ptstate[vt][I_AVLW]))
  {
    return test = temkey(I_AVLW)+1;
  }
  if (err[vt][I_RGRW]  > fabs(ptol * ptstate[vt][I_RGRW]))
  {
    return test = temkey(I_RGRW)+1;
  }
  if (err[vt][I_SNWPCK]  > fabs(ptol * ptstate[vt][I_SNWPCK]))
  {
    return test = temkey(I_SNWPCK)+1;
  }
  if (err[vt][I_SGRW]  > fabs(ptol * ptstate[vt][I_SGRW]))
  {
    return test = temkey(I_SGRW)+1;
  }
  if (err[vt][I_RPERC]  > fabs((ptol) * ptstate[vt][I_RPERC]))
  {
    return test = temkey(I_RPERC)+1;
  }
  if (err[vt][I_EET]  > fabs((ptol) * ptstate[vt][I_EET]))
  {
    return test = temkey(I_EET)+1;
  }

// added for 3-box hydrology model
    if (err[vt][I_AVLW1]  > fabs(ptol * ptstate[vt][I_AVLW1]))
  {
    return test = temkey(I_AVLW1)+1;
  }

    if (err[vt][I_AVLW2]  > fabs(ptol * ptstate[vt][I_AVLW2]))
  {
    return test = temkey(I_AVLW2)+1;
  }
    if (err[vt][I_AVLW3]  > fabs(ptol * ptstate[vt][I_AVLW3]))
  {
    return test = temkey(I_AVLW3)+1;
  }

  return test;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

// modified nminxclm() correspondingly for 3-box hydrology
// to reflect drainage issue

//void TTEM::delta(const int& dm, double pstate[], double pdstate[])

void TTEM::delta(const int& dm, double pstate[NVT][NUMEQ], double pdstate[NVT][NUMEQ], const int& vt)
{
	





  nfert = 1.0000000;

  //added for soil hydrology model
   hydm.soil_perco (soil.pctsilt, soil.pctclay);


 /*  if (atms.tair[vt][dm] < 0.0 && atms.tair[vt][dm-1] >= 0.0)

{
	pstate[vt][I_AVLW2] = pstate[vt][I_AVLW2] * 0.1;
	soil.freeze[vt][dm] = pstate[vt][I_AVLW2];
}
 
   */

/* if  (soil.temp == 2  ||  soil.temp ==3)
 {
	 pstate[vt][I_AVLW2] = 0.1 * pstate[vt][I_AVLW2];  
 
 soil.freeze[vt][dm] = pstate[vt][I_AVLW2];   
 soil.temp = soil.temp +1;}


 */


   //hydm.lyperc[vt][dm] =hyd.daze[dm]*2* pow((pstate[vt][I_AVLW]/soil.awcapmm[vt]), 4);                     wsr
   


//.lyperc1[vt][dm] = hyd.daze[dm]*hydm.k1* pow((pstate[I_AVLW1]/soil.awcapmm1), hydm.k2);
 //moss hydrological conductivity for summer, the maximum value is 0.022mm/s
   
 // printf("%f   %f   %f     %f    \n",atms.tair[vt][dm],pstate[vt][I_AVLW3],pdstate[vt][I_AVLW3],y[vt][I_AVLW3]);



   //add soil moist SM
   
   hydm.lyperc1[vt][dm] = hyd.daze[dm]*6* pow((pstate[vt][I_AVLW1]/soil.awcapmm1[vt]), 4); // k1 will be pretty big 7.0
  // printf("hydm.lyperc1[vt][dm] hydm.2yperc1[vt][dm] \n");
 //  printf("     %f    %f  n",hydm.lyperc1[vt][dm],hydm.lyperc2[vt][dm]);
//   hydm.lyperc1[dm] = pstate[I_AVLW1] * 0.5; // k1 will be pretty big 7.0

   hydm.lyperc2[vt][dm] = hyd.daze[dm]*4* pow((pstate[vt][I_AVLW2]/soil.awcapmm2[vt]), 4);
//   hydm.lyperc2[dm] = hyd.daze[dm]*2.0* pow((pstate[I_AVLW2]/soil.awcapmm2), hydm.k2);

 //added for hydrology model
 
 
 
 
  atms.eet[vt][dm]= hydm.hydmaet(dm, atms.rain[vt][dm], soil.snowinf[vt][dm], pstate[vt][I_AVLW1], pstate[vt][I_AVLW2],pstate[vt][I_AVLW3], atms.pet[vt][dm], soil.awcapmm1[vt], soil.awcapmm2[vt],soil.awcapmm3[vt], hyd.daze[dm], hydm.root1[vt], hydm.root2[vt],hydm.root3[vt], vt);
  
 

  
 // soil.percol(atms.rain[vt][dm],soil.snowinf[vt][dm],atms.eet[vt][dm],pstate[vt][I_AVLW],dm,vt);    wsr: reduce original one box perc
  
  
  
  
  

  soil.percol_hydm(atms.tair[vt][dm], atms.rain[vt][dm], soil.snowinf[vt][dm], hydm.trans1[vt], hydm.trans2[vt], hydm.trans3[vt], pstate[vt][I_AVLW1],pstate[vt][I_AVLW2],pstate[vt][I_AVLW3],hydm.lyperc1[vt][dm],hydm.lyperc2[vt][dm],dm, vt);

// The following is orignal version
  //atms.eet[dm] = atms.xeet(atms.rain[dm], soil.snowinf[dm], atms.pet[dm],
    //                       pstate[I_AVLW], soil.awcapmm, dm);

//  soil.percol(atms.rain[dm], soil.snowinf[dm], atms.eet[dm],
  //            pstate[I_AVLW],dm);

//printf("%f    %f     %f     %f    %f\n",pstate[vt][I_AVLW],pstate[vt][I_AVLW1],pstate[vt][I_AVLW2],pstate[vt][I_AVLW3],atms.eet[vt][dm]);

  if ((pstate[vt][I_AVLW1]+pstate[vt][I_AVLW2]+pstate[vt][I_AVLW3]+ soil.snowinf[vt][dm] + atms.rain[vt][dm] - atms.eet[vt][dm]
    -soil.rperc1[vt][dm] - soil.sperc1[vt][dm] - soil.rperc2[vt][dm] - soil.sperc2[vt][dm] - soil.rperc3[vt][dm] - soil.sperc3[vt][dm]) < 0.0)
  {
   atms.eet[vt][dm] = pstate[vt][I_AVLW1] + pstate[vt][I_AVLW2] + pstate[vt][I_AVLW3] + soil.snowinf[vt][dm] + atms.rain[vt][dm]
                - soil.rperc1[vt][dm] - soil.sperc1[vt][dm] - soil.rperc2[vt][dm] - soil.sperc2[vt][dm] - soil.rperc3[vt][dm] - soil.sperc3[vt][dm];   // wsr: change total AVLW, sperc and rperc to seperate ones
  
  
  }
  
  
  
  
 // printf("%9.3f,%9.3f\n",atms.eet[vt][dm],atms.pet[vt][dm]);
//  printf("%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f\n",atms.pet[vt][dm],hyd.htrans[vt][dm],hyd.hevap[vt][dm],hyd.soil_evap[vt][dm],hyd.snowsub[vt][dm],hyd.sub_from_canopy[vt][dm],atms.eet[vt][dm],atms.rain[vt][dm],soil.snowinf[vt][dm]);
//if((atms.eet[vt][dm] <0 )|| (atms.eet[vt][dm]>250))   {atms.eet[vt][dm]=5.5;}

  //soil.rrun[vt][dm] = soil.rrunoff(pstate[vt][I_RGRW], soil.rperc[vt][dm]);
  //soil.srun[vt][dm] = soil.srunoff(elev, atms.tair[vt][dm], atms.prevtair,
                              // atms.prev2tair, pstate[vt][I_SGRW], soil.sperc[vt][dm]);

 // soil.h2oyld[vt][dm] = soil.rrun[vt][dm] + soil.srun[vt][dm];                                        wsr: reduce all the non-three-box variables


  if (moistlim == 0)
  {
    veg.unnormleaf[vt][dm] = veg.deltaleaf(veg.cmnt, atms.pet[vt][dm],
                                       atms.prvpetmx, prevy[vt][I_UNRMLF]);
  }
  else
  {
    veg.unnormleaf[vt][dm] = veg.deltaleaf(veg.cmnt, atms.eet[vt][dm],
                                       atms.prveetmx, prevy[vt][I_UNRMLF]);
  }
  if (veg.prvleafmx[veg.cmnt] <= 0.0) { veg.leaf[vt][dm] = 0.0; }
  else { veg.leaf[vt][dm] = veg.unnormleaf[vt][dm]/veg.prvleafmx[veg.cmnt]; }
  if (veg.leaf[vt][dm] < veg.minleaf[veg.cmnt])
  {
    veg.leaf[vt][dm] = veg.minleaf[veg.cmnt];
  }
  if (veg.leaf[vt][dm] > 1.0) { veg.leaf[vt][dm] = 1.0; }

  veg.alleaf[vt] = veg.leafmxc[veg.cmnt]/(1.0 + veg.kleafc[veg.cmnt]
               * exp(veg.cov[veg.cmnt]*pstate[vt][I_VEGC]));
  veg.lai[vt][dm] = veg.sla[veg.cmnt] * veg.alleaf[vt];
  veg.foliage[vt] = veg.alleaf[vt] / veg.leafmxc[veg.cmnt];
  veg.fpc[vt][dm] = 1.0;

  deltaxclm(veg.cmnt, soil.pcfldcap[vt], dm,vt);

  if (ag.state == 0)
  {
    // added veg.thawpercent for thaw-frozen mechanism by qianlai 28/08/2000

    veg.ingpp[vt][dm] = veg.gppxclm(veg.cmnt, atms.co2[dm], atms.par[vt][dm],
                            temp[vt], gv[vt], veg.leaf[vt][dm], veg.foliage[vt],veg.thawpercent[vt][dm],vt );
    if (veg.ingpp[vt][dm] < 0.0) { veg.ingpp[vt][dm] = 0.0; }
//printf("%f     %f\n",pstate[vt][I_SM],pstate[vt][I_SM2]);
    veg.inuptake[vt] = veg.nupxclm(veg.cmnt, pstate[vt][I_SM1],pstate[vt][I_SM2],pstate[vt][I_SM3], pstate[vt][I_AVLN],                 //wsr: use VSM instead of soilh2o
                               respq10[vt], ksoil[vt], veg.foliage[vt],vt);
  }
  else
  {
    veg.ingpp[vt][dm] = 0.0;
    veg.inuptake[vt] = 0.0;
  }
  //printf("                                                         %f\n",veg.ingpp[vt][dm]);
   microbe.rh[vt][dm] = microbe.rhxclm(pstate[vt][I_SOLC], dq10[vt], rhmoist[vt],vt);

    if (microbe.rh[vt][dm] < 0.0) { microbe.rh[vt][dm] = 0.0; }

  soil.ninput[vt][dm] = 0.0;

  // modified by QZ Oct 09 2000, using organic layer moisture, replace I_SM with I_SM2

  microbe.netnmin[vt][dm] = microbe.nminxclm(veg.cmnt, dm, pstate[vt][I_SM1],pstate[vt][I_SM2],                        //wsr
                                         pstate[vt][I_SOLC], pstate[vt][I_SOLN],
                                         pstate[vt][I_AVLN], microbe.decay[vt],
                                         microbe.rh[vt][dm], ksoil[vt], vt);

  if (ag.state == 0)
  {
    veg.ltrfal[vt][dm].carbon = veg.cfall[veg.cmnt] * pstate[vt][I_VEGC];
    if (veg.ltrfal[vt][dm].carbon < 0.0) { veg.ltrfal[vt][dm].carbon = 0.0; }
	
    veg.ltrfal[vt][dm].nitrogen = veg.nfall[veg.cmnt] * pstate[vt][I_STRN];
    if (veg.ltrfal[vt][dm].nitrogen < 0.0) { veg.ltrfal[vt][dm].nitrogen = 0.0; }
//printf("%f\n",pstate[vt][I_SOLC]);
   // added for thawpercent
  // general formula is veg.thawpercent = frond[dm] / the maximum depth of active layer
  // for this particular site, we choose the depth is about 2.0
    if ((dm >=5 ) && (dm <=8)) veg.thawpercent[vt][dm] = 0.50;
//    if (dm ==4) veg.thawpercent[dm] = 0.30;

    veg.rm[vt][dm] = veg.rmxclm(veg.cmnt,pstate[vt][I_VEGC],respq10[vt], veg.thawpercent[vt][dm],vt);
    if (veg.rm[vt][dm] < 0.0) { veg.rm[vt][dm] = 0.0; }

    veg.innpp[vt][dm] = veg.ingpp[vt][dm] - veg.rm[vt][dm];

    veg.rg[vt][dm] = 0;
    if (veg.innpp[vt][dm] > 0.0)
    {
      veg.rg[vt][dm]  = 0.2 * veg.innpp[vt][dm];
      veg.innpp[vt][dm] *= 0.8;
    }

    if (veg.inuptake[vt] > pstate[vt][I_AVLN] + (nfert*soil.ninput[vt][dm])
    + microbe.netnmin[vt][dm])
    {
      veg.inuptake[vt] = pstate[vt][I_AVLN] + (nfert*soil.ninput[vt][dm])
                     + microbe.netnmin[vt][dm];
    }

    if (veg.inuptake[vt] < 0.0) { veg.inuptake[vt] = 0.0; }

    veg.nuptake[vt][dm] = veg.inuptake[vt];
    veg.suptake[vt][dm] = veg.nuptake[vt][dm];
    veg.luptake[vt][dm] = 0.0;
    veg.gpp[vt][dm] = veg.ingpp[vt][dm];
    veg.npp[vt][dm] = veg.innpp[vt][dm];
    veg.nmobil[vt][dm] = 0.0;
    veg.nresorb[vt][dm] = 0.0;
	//printf("%f                           %f                  %f\n",veg.gpp[vt][dm],veg.ingpp[vt][dm],veg.npp[vt][dm]);
// Nitrogen feedback of GPP (nfeed == 1)

    if (nfeed == 1)                                                                                													//wsr: close the feedback for debugging 
    {
      if (veg.inuptake[vt] == 0.0) { veg.inuptake[vt] = 0.000001; }
      veg.inprodcn[vt] = veg.innpp[vt][dm] / (veg.inuptake[vt] + pstate[vt][I_STON]);

      if (veg.ltrfal[vt][dm].nitrogen <= veg.ltrfal[vt][dm].carbon/veg.cneven[veg.cmnt])
      {
	veg.nresorb[vt][dm] = veg.ltrfal[vt][dm].carbon/veg.cneven[veg.cmnt]
                              - veg.ltrfal[vt][dm].nitrogen;
      }
      else
      {
	veg.ltrfal[vt][dm].nitrogen = veg.ltrfal[vt][dm].carbon/veg.cneven[veg.cmnt];
	veg.nresorb[vt][dm] = 0.0;
      }
      if (pstate[vt][I_VEGC] > 0.0)
      {
        veg.nresorb[vt][dm] *= (pstate[vt][I_STRN]/pstate[vt][I_VEGC]) * veg.c2n[vt];
      }

      if (veg.inprodcn[vt] > veg.cneven[veg.cmnt])
      {
	veg.npp[vt][dm] = veg.cneven[veg.cmnt] * (veg.nuptake[vt][dm] + pstate[vt][I_STON]);
	if (veg.npp[vt][dm] < 0.0) { veg.npp[vt][dm] = 0.0; }
	veg.rg[vt][dm] = 0.25 * veg.npp[vt][dm];
  	veg.gpp[vt][dm] = veg.npp[vt][dm] + veg.rg[vt][dm] + veg.rm[vt][dm];

	//printf("%f     %f            %f                %f     %f\n",veg.gpp[vt][dm],veg.npp[vt][dm],veg.ingpp[vt][dm],veg.innpp[vt][dm],pstate[vt][I_AVLN]);


	if (veg.gpp[vt][dm] < 0.0) { veg.gpp[vt][dm] = 0.0; }

	veg.nmobil[vt][dm] = pstate[vt][I_STON];
      }

      if (veg.inprodcn[vt] <= veg.cneven[veg.cmnt])
      {
        veg.nuptake[vt][dm] = veg.inuptake[vt] * (veg.inprodcn[vt] - veg.cnmin[veg.cmnt])
                          * (veg.inprodcn[vt] - 2*veg.cneven[veg.cmnt]
                          + veg.cnmin[veg.cmnt]);
	veg.nuptake[vt][dm] /= ((veg.inprodcn[vt] - veg.cnmin[veg.cmnt])
                           * (veg.inprodcn[vt] - 2*veg.cneven[veg.cmnt]
                           + veg.cnmin[veg.cmnt])) - pow(veg.inprodcn[vt]
                           - veg.cneven[veg.cmnt],2.0);
	if (veg.nuptake[vt][dm] < 0.0) { veg.nuptake[vt][dm] = 0.0; }
	if (pstate[vt][I_STON] >= veg.npp[vt][dm]/veg.cneven[veg.cmnt])
        {
	   veg.nmobil[vt][dm] = veg.npp[vt][dm]/veg.cneven[veg.cmnt];
	   if (veg.nmobil[vt][dm] < 0.0 && pstate[vt][I_VEGC] > 0.0)
           {
	     veg.nmobil[vt][dm] *= (pstate[vt][I_STRN]/pstate[vt][I_VEGC]) * veg.c2n[vt];
	   }
	   veg.suptake[vt][dm] = 0.0;
	}
	else
        {
	  veg.nmobil[vt][dm] = pstate[vt][I_STON];
	  veg.suptake[vt][dm] = (veg.npp[vt][dm]/veg.cneven[veg.cmnt])
                             - veg.nmobil[vt][dm];
	  if (veg.suptake[vt][dm] < 0.0) { veg.suptake[vt][dm] = 0.0; }
	  if (veg.suptake[vt][dm] > veg.nuptake[vt][dm])
          {
            veg.suptake[vt][dm] = veg.nuptake[vt][dm];
          }
	}

    // ptate[1] changed to pstate[I_STRN] by DWK on 20000201
	if ((pstate[vt][I_STON] + veg.nuptake[vt][dm] - veg.suptake[vt][dm]
           + veg.nresorb[vt][dm] - veg.nmobil[vt][dm]) < (veg.labncon[veg.cmnt]
           * (pstate[vt][I_STRN] + veg.suptake[vt][dm] - veg.ltrfal[vt][dm].nitrogen
           - veg.nresorb[vt][dm] + veg.nmobil[vt][dm])))
        {
	  veg.luptake[vt][dm] = veg.nuptake[vt][dm] - veg.suptake[vt][dm];
	}
	else
        {
	   veg.luptake[vt][dm] = (veg.labncon[veg.cmnt] * (pstate[vt][I_STRN]
                             + veg.suptake[vt][dm] - veg.ltrfal[vt][dm].nitrogen
                             - veg.nresorb[vt][dm] + veg.nmobil[vt][dm]))
                             - (pstate[vt][I_STON] + veg.nresorb[vt][dm]
                             - veg.nmobil[vt][dm]);
	   if (veg.luptake[vt][dm] < 0.0) { veg.luptake[vt][dm] = 0.0; }
	   veg.nuptake[vt][dm] = veg.suptake[vt][dm] + veg.luptake[vt][dm];
	}
      }
    }

    veg.gpr[vt][dm] = veg.rm[vt][dm] + veg.rg[vt][dm];
    nep[vt][dm] = veg.npp[vt][dm] - microbe.rh[vt][dm];
    cflux[vt][dm] = nep[vt][dm];

    ag.npp[vt][dm].carbon = 0.0;
    ag.npp[vt][dm].nitrogen = 0.0;
    ag.fertn[vt][dm] = 0.0;
    ag.ltrfal[vt][dm].carbon = 0.0;
    ag.ltrfal[vt][dm].nitrogen = 0.0;
    ag.slash[vt][dm].carbon = 0.0;
    ag.slash[vt][dm].nitrogen = 0.0;
    ag.sconvrtflx[vt][dm].carbon = 0.0;
    ag.sconvrtflx[vt][dm].nitrogen = 0.0;
    ag.nsretent[vt][dm] = 0.0;
  }
  else
  {
    veg.ltrfal[vt][dm].carbon = 0.0;
    veg.ltrfal[vt][dm].nitrogen = 0.0;
    veg.gpr[vt][dm] = 0.0;
    veg.rm[vt][dm] = 0.0;
    veg.rg[vt][dm] = 0.0;
    veg.innpp[vt][dm] = 0.0;
    veg.nuptake[vt][dm] = 0.0;
    veg.suptake[vt][dm] = 0.0;
    veg.luptake[vt][dm] = 0.0;
    veg.gpp[vt][dm] = 0.0;
    veg.npp[vt][dm] = 0.0;
    veg.nmobil[vt][dm] = 0.0;
    veg.nresorb[vt][dm] = 0.0;

    ag.npp[vt][dm].carbon = ag.RAP * ag.potnpp[vt][dm];
    if (ag.npp[vt][dm].carbon < 0.0) { ag.npp[vt][dm].carbon = 0.0; }
    ag.npp[vt][dm].nitrogen = ag.npp[vt][dm].carbon / ag.c2n[vt];
    ag.fertn[vt][dm] = 0.0;
    ag.ltrfal[vt][dm].carbon = ag.npp[vt][dm].carbon * ag.cfall[veg.cmnt];
    ag.ltrfal[vt][dm].nitrogen = ag.npp[vt][dm].nitrogen * ag.nfall[veg.cmnt];
    nep[vt][dm] = ag.npp[vt][dm].carbon - microbe.rh[vt][dm];
    cflux[vt][dm] = nep[vt][dm] - ag.convrtflx[vt][dm].carbon;
  }

// Changes in available nitrogen in soil

  if (avlnflag == 1)
  {
    if (ag.state == 0)
    {
      soil.nlost[vt][dm] = pstate[vt][I_AVLN] / ((pstate[vt][I_SM1]*1/3+pstate[vt][I_SM2]*2/3)*10);   //wsr
      soil.nlost[vt][dm] *= ((atms.rain[vt][dm] + soil.snowinf[vt][dm] - atms.eet[vt][dm])                                          
                        + (soil.rootz[vt] * 1000.0)) / (soil.rootz[vt] * 1000.0);
      soil.nlost[vt][dm] *= soil.nloss[veg.cmnt];

      if (soil.nlost[vt][dm] > pstate[vt][I_AVLN] - veg.nuptake[vt][dm]
         + microbe.netnmin[vt][dm] + soil.ninput[vt][dm])
      {
        soil.nlost[vt][dm] = pstate[vt][I_AVLN] - veg.nuptake[vt][dm] + microbe.netnmin[vt][dm]
                         + soil.ninput[vt][dm];
      }
      if (soil.nlost[vt][dm] < 0.0)
      {
        soil.nlost[vt][dm]  = 0.0;
        microbe.netnmin[vt][dm] = soil.nlost[vt][dm] + veg.nuptake[vt][dm] - soil.ninput[vt][dm]
                              - pstate[vt][I_AVLN];
      }
    }
    else
    {
      soil.ninput[vt][dm] = ag.nretent[vt][dm];
      soil.nlost[vt][dm] = 0.0;

      // Condition added to limit addition of agfertn to only
      //   periods of crop growth
      if (ag.npp[vt][dm].nitrogen > 0.0)
      {
        if (pstate[vt][I_AVLN] + soil.ninput[vt][dm] + microbe.netnmin[vt][dm]
            - soil.nlost[vt][dm] < ag.npp[vt][dm].nitrogen)
        {
   	   ag.fertn[vt][dm] = ag.npp[vt][dm].nitrogen + soil.nlost[vt][dm] - pstate[vt][I_AVLN]
           - soil.ninput[vt][dm] - microbe.netnmin[vt][dm];
	   if (ag.fertn[vt][dm] < 0.0)
           {
	      ag.fertn[vt][dm] = 0.0;
	      microbe.netnmin[vt][dm] = soil.nlost[vt][dm] + ag.npp[vt][dm].nitrogen
                                    - soil.ninput[vt][dm] - pstate[vt][I_AVLN];
	   }
	   soil.ninput[vt][dm] += ag.fertn[vt][dm];
        }
      }
    }
  }
  else
  {
    soil.nlost[vt][dm] = soil.ninput[vt][dm] - veg.nuptake[vt][dm] - ag.npp[vt][dm].nitrogen
                     + microbe.netnmin[vt][dm];
  }

// Describe monthly changes to carbon pools and fluxes for ODE state variables
// (i.e., pdstate)

  // Carbon pools in natural ecosystems

  pdstate[vt][I_VEGC] = veg.gpp[vt][dm] - veg.gpr[vt][dm] - veg.ltrfal[vt][dm].carbon;
  pdstate[vt][I_SOLC] = veg.ltrfal[vt][dm].carbon + ag.ltrfal[vt][dm].carbon
                    + ag.slash[vt][dm].carbon - ag.sconvrtflx[vt][dm].carbon
                    - microbe.rh[vt][dm];

  // Nitrogen pools in natural ecosystems

  pdstate[vt][I_STRN] = veg.suptake[vt][dm] - veg.ltrfal[vt][dm].nitrogen - veg.nresorb[vt][dm]
                    + veg.nmobil[vt][dm];
  pdstate[vt][I_SOLN] = veg.ltrfal[vt][dm].nitrogen + ag.ltrfal[vt][dm].nitrogen
                    + ag.slash[vt][dm].nitrogen - ag.sconvrtflx[vt][dm].nitrogen
                    - ag.nsretent[vt][dm] - microbe.netnmin[vt][dm];
  pdstate[vt][I_AVLN] = soil.ninput[vt][dm] - soil.nlost[vt][dm] + microbe.netnmin[vt][dm]
                    - veg.nuptake[vt][dm] - ag.npp[vt][dm].nitrogen;
  pdstate[vt][I_STON] = veg.luptake[vt][dm] + veg.nresorb[vt][dm] - veg.nmobil[vt][dm];

  // Human product pools

  pdstate[vt][I_AGPRDC] = 0.0;
  pdstate[vt][I_AGPRDN] = 0.0;
  pdstate[vt][I_PROD10C] = 0.0;
  pdstate[vt][I_PROD10N] = 0.0;
  pdstate[vt][I_PROD100C] = 0.0;
  pdstate[vt][I_PROD100N] = 0.0;
  pdstate[vt][I_TOTPRDC] = 0.0;
  pdstate[vt][I_TOTPRDN] = 0.0;
  pdstate[vt][I_TOTEC] = 0.0;
  pdstate[vt][I_TOTGC] = 0.0;

  // Water pools
  
 

  //pdstate[vt][I_AVLW] = soil.snowinf[vt][dm] + atms.rain[vt][dm] - atms.eet[vt][dm]
      //              - soil.rperc[vt][dm] - soil.sperc[vt][dm];                                                  wsr
// added for 3-box hydrological model

  pdstate[vt][I_AVLW1] = soil.snowinf[vt][dm] + atms.rain[vt][dm] - hydm.lyperc1[vt][dm] - hydm.trans1[vt] - soil.rperc1[vt][dm] - soil.sperc1[vt][dm]-deltafreeze1;
  if (pstate[vt][I_AVLW1]+pdstate[vt][I_AVLW1] <= 0.0) {
	 pdstate[vt][I_AVLW1] = 0.001 - pstate[vt][I_AVLW1];
  }
 // printf("      %f      \n",soil.snowinf[vt][dm]);
//soil.freeze2(atms.tair[vt][dm],pstate[vt][I_AVLW2],vt,dm);
  
  pdstate[vt][I_AVLW2] = hydm.lyperc1[vt][dm]- hydm.trans2[vt] - soil.rperc2[vt][dm] - soil.sperc2[vt][dm] - hydm.lyperc2[vt][dm]-deltafreeze2;
  if (pstate[vt][I_AVLW2]+pdstate[vt][I_AVLW2] <= 0.0) {
	 pdstate[vt][I_AVLW2] = 0.001 - pstate[vt][I_AVLW2];
  }
 // printf("    %f    %f     %f      %f     %f\n",deltafreeze1,deltafreeze2,deltafreeze3,atms.tair[vt][dm],pdstate[vt][I_AVLW2]);


 //pdstate[vt][I_AVLW3] = hydm.lyperc2[vt][dm]- hydm.trans3[vt] - soil.rperc3[vt][dm] - soil.sperc3[vt][dm];
  pdstate[vt][I_AVLW3] = hydm.lyperc2[vt][dm]-hydm.trans3[vt] - soil.rperc3[vt][dm] - soil.sperc3[vt][dm];
  if (pstate[vt][I_AVLW3]+pdstate[vt][I_AVLW3] <= 0.0) {
	 pdstate[vt][I_AVLW3] = 0.001 - pstate[vt][I_AVLW3];
  }

// end of changed for available water pool

 // pdstate[vt][I_RGRW] = soil.rperc[vt][dm] - soil.rrun[vt][dm];                                    wsr
  pdstate[vt][I_SNWPCK] = atms.snowfall[vt][dm] - soil.snowinf[vt][dm];
  //pdstate[vt][I_SGRW] = soil.sperc[vt][dm] - soil.srun[vt][dm];                                      wsr
 // pdstate[vt][I_SM] = pdstate[vt][I_AVLW];                                                                 wsr


 // influenced by frozen, decreasing soil available water, Sept/19/2000, Q.Z.

  pdstate[vt][I_SM1] = pdstate[vt][I_AVLW1];  // added
  pdstate[vt][I_SM2] = pdstate[vt][I_AVLW2]; // added
  pdstate[vt][I_SM3] = pdstate[vt][I_AVLW3]; // added

 // pdstate[vt][I_PCTP] = 100.0 * pdstate[vt][I_SM]/soil.totpor[vt];                              wsr
  pdstate[vt][I_PCTP1] = 100.0 * pdstate[vt][I_SM1]/soil.totpor1[vt];
 pdstate[vt][I_PCTP2] = 100.0 * pdstate[vt][I_SM2]/soil.totpor2[vt];
  pdstate[vt][I_PCTP3] = 100.0 * pdstate[vt][I_SM3]/soil.totpor3[vt];

 // pdstate[vt][I_VSM] = pdstate[vt][I_SM]/(soil.rootz[vt]*1000.0);                              wsr
 // if (pstate[vt][I_VSM]+pdstate[vt][I_VSM] <= 0.0) {
  //  pdstate[vt][I_VSM] = 0.001 - pstate[vt][I_VSM];
  //}

 // added for 3-box hydrology
pdstate[vt][I_VSM1] = pdstate[vt][I_SM1]/(soil.dpwbox1[vt]*1000.0);
 if (pstate[vt][I_VSM1]+pdstate[vt][I_VSM1] <= 0.0) {
 pdstate[vt][I_VSM1] = 0.001 - pstate[vt][I_VSM1];
 }
 // pdstate[I_VSM2] = pdstate[I_SM2]/(soil.rootmx*1000.0-soil.dpwbox1*1000.0);
  pdstate[vt][I_VSM2] = pdstate[vt][I_SM2]/(soil.dpwbox2[vt]*1000.0);
  if (pstate[vt][I_VSM2]+pdstate[vt][I_VSM2] <= 0.0) {
	 pdstate[vt][I_VSM2] = 0.001 - pstate[vt][I_VSM2];
  }

 // pdstate[I_VSM3] = pdstate[I_SM3]/(soil.rootmx*1000.0-soil.dpwbox1*1000.0); // temporally assume it is equal to layer2
  pdstate[vt][I_VSM3] = pdstate[vt][I_SM3]/(soil.dpwbox3[vt]*1000.0); // temporally assume it is equal to layer2
  if (pstate[vt][I_VSM3]+pdstate[vt][I_VSM3] <= 0.0) {
	 pdstate[vt][I_VSM3] = 0.001 - pstate[vt][I_VSM3];
  }

// end of adding
  //printf("                                     %f\n",veg.gpp[vt][dm]);
  // Carbon fluxes in natural ecosystems

  pdstate[vt][I_INGPP] = veg.ingpp[vt][dm];
  pdstate[vt][I_GPP] = veg.gpp[vt][dm];
  pdstate[vt][I_INNPP] = veg.innpp[vt][dm];
  pdstate[vt][I_NPP] = veg.npp[vt][dm];
  pdstate[vt][I_GPR] = veg.gpr[vt][dm];
  pdstate[vt][I_RVMNT] = veg.rm[vt][dm];
  pdstate[vt][I_RVGRW] = veg.rg[vt][dm];
  pdstate[vt][I_LTRC] = veg.ltrfal[vt][dm].carbon;
  pdstate[vt][I_RH] = microbe.rh[vt][dm];
  pdstate[vt][I_NEP] = nep[vt][dm];

  // Nitrogen fluxes in natural ecosystems

  pdstate[vt][I_NINP] = soil.ninput[vt][dm];
  pdstate[vt][I_INNUP] = veg.inuptake[vt];
  pdstate[vt][I_VNUP] = veg.nuptake[vt][dm];
  pdstate[vt][I_VSUP] = veg.suptake[vt][dm];
  pdstate[vt][I_VLUP] = veg.luptake[vt][dm];
  pdstate[vt][I_VNMBL] = veg.nmobil[vt][dm];
  pdstate[vt][I_VNRSRB] = veg.nresorb[vt][dm];
  pdstate[vt][I_LTRN] = veg.ltrfal[vt][dm].nitrogen;
  pdstate[vt][I_MNUP] = microbe.nuptake[vt][dm];
  pdstate[vt][I_NMIN] = microbe.netnmin[vt][dm];
  pdstate[vt][I_NLST] = soil.nlost[vt][dm];

  // Water fluxes

  pdstate[vt][I_RAIN] = atms.rain[vt][dm];
  //pdstate[vt][I_RPERC] = soil.rperc[vt][dm];                      wsr
  //pdstate[vt][I_RRUN] = soil.rrun[vt][dm];
  pdstate[vt][I_SNWFAL] = atms.snowfall[vt][dm];
  pdstate[vt][I_SNWINF] = soil.snowinf[vt][dm];
  //pdstate[vt][I_SPERC] = soil.sperc[vt][dm];
 // pdstate[vt][I_SRUN] = soil.srun[vt][dm];
  pdstate[vt][I_PET] = atms.pet[vt][dm];
  pdstate[vt][I_EET] = atms.eet[vt][dm];
//  pdstate[vt][I_WYLD] = soil.rrun[vt][dm] + soil.srun[vt][dm];

// added for 2-box soil hydrology model
  //pdstate[vt][I_LYPERC] = hydm.lyperc[vt][dm];
  pdstate[vt][I_LYPERC1] = hydm.lyperc1[vt][dm];
  pdstate[vt][I_LYPERC2] = hydm.lyperc2[vt][dm];

  pdstate[vt][I_RPERC1] = soil.rperc1[vt][dm];
  pdstate[vt][I_RPERC2] = soil.rperc2[vt][dm];
  pdstate[vt][I_RPERC3] = soil.rperc3[vt][dm];

  pdstate[vt][I_SPERC1] = soil.sperc1[vt][dm];
  pdstate[vt][I_SPERC2] = soil.sperc2[vt][dm];
  pdstate[vt][I_SPERC3] = soil.sperc3[vt][dm];
// end of adding


 // added for canopy hydrology model
  pdstate[vt][I_EVAP] = hyd.hevap[vt][dm];
  pdstate[vt][I_TRANS] = hyd.htrans[vt][dm];
  pdstate[vt][I_SEVAP] = hyd.soil_evap[vt][dm];
  pdstate[vt][I_SNOWSUB] = hyd.snowsub[vt][dm];
  pdstate[vt][I_SUBCAN] = hyd.sub_from_canopy[vt][dm];
// end of adding

   // for soil temperature
  pdstate[vt][I_TSOIL] = atms.tsoil[vt][dm];
  pdstate[vt][I_DST5] = atms.dst5[vt][dm];
  pdstate[vt][I_DST10] = atms.dst10[vt][dm];
  pdstate[vt][I_DST20] = atms.dst20[vt][dm];
  pdstate[vt][I_DST50] = atms.dst50[vt][dm];
  pdstate[vt][I_DST100] = atms.dst100[vt][dm];
  pdstate[vt][I_DST200] = atms.dst200[vt][dm];
  pdstate[vt][I_FRONTD] = atms.frontd[vt][dm];
  pdstate[vt][I_THAWBE] = atms.thawbe[vt][dm];
  pdstate[vt][I_THAWEND] = atms.thawend[vt][dm];
// end of ...



  // Phenology

  pdstate[vt][I_UNRMLF] = veg.unnormleaf[vt][dm];
  pdstate[vt][I_LEAF] = veg.leaf[vt][dm];

  // Carbon and nitrogen fluxes from agricultural conversion

  pdstate[vt][I_CNVRTC] = ag.convrtflx[vt][dm].carbon;
  pdstate[vt][I_CNVRTN] = ag.convrtflx[vt][dm].nitrogen;
  pdstate[vt][I_SCNVRTC] = ag.sconvrtflx[vt][dm].carbon;
  pdstate[vt][I_SCNVRTN] = ag.sconvrtflx[vt][dm].nitrogen;
  pdstate[vt][I_NVRTNT] = ag.nvretent[vt][dm];
  pdstate[vt][I_NSRTNT] = ag.nsretent[vt][dm];
  pdstate[vt][I_NRETNT] = ag.nretent[vt][dm];
  pdstate[vt][I_SLASHC] = ag.slash[vt][dm].carbon;
  pdstate[vt][I_SLASHN] = ag.slash[vt][dm].nitrogen;
  pdstate[vt][I_PRDF10C] = ag.formPROD10[vt].carbon / (double) CYCLE;
  pdstate[vt][I_PRDF10N] = ag.formPROD10[vt].nitrogen / (double) CYCLE;
  pdstate[vt][I_PRDF100C] = ag.formPROD100[vt].carbon / (double) CYCLE;
  pdstate[vt][I_PRDF100N] = ag.formPROD100[vt].nitrogen / (double) CYCLE;

  // Carbon and nitrogen fluxes in agricultural ecosystems

  pdstate[vt][I_AGNPPC] = ag.npp[vt][dm].carbon;
  pdstate[vt][I_AGNPPN] = ag.npp[vt][dm].nitrogen;
  pdstate[vt][I_AGFPRDC] = ag.npp[vt][dm].carbon - ag.ltrfal[vt][dm].carbon;
  pdstate[vt][I_AGFPRDN] = ag.npp[vt][dm].nitrogen - ag.ltrfal[vt][dm].nitrogen;
  pdstate[vt][I_AGLTRC] = ag.ltrfal[vt][dm].carbon;
  pdstate[vt][I_AGLTRN] = ag.ltrfal[vt][dm].nitrogen;
  pdstate[vt][I_AGFRTN] = ag.fertn[vt][dm];

  // Carbon and nitrogen fluxes from human product pools

  pdstate[vt][I_TOTFPRDC] = ag.formTOTPROD[vt].carbon / (double) CYCLE;
  pdstate[vt][I_TOTFPRDN] = ag.formTOTPROD[vt].nitrogen / (double) CYCLE;
  pdstate[vt][I_AGPRDFC] = ag.PROD1decay[vt].carbon / (double) CYCLE;
  pdstate[vt][I_AGPRDFN] = ag.PROD1decay[vt].nitrogen / (double) CYCLE;
  pdstate[vt][I_PRD10FC] = ag.PROD10decay[vt].carbon / (double) CYCLE;
  pdstate[vt][I_PRD10FN] = ag.PROD10decay[vt].nitrogen / (double) CYCLE;
  pdstate[vt][I_PRD100FC] = ag.PROD100decay[vt].carbon / (double) CYCLE; 
  pdstate[vt][I_PRD100FN] = ag.PROD100decay[vt].nitrogen / (double) CYCLE;
  pdstate[vt][I_TOTPRDFC] = ag.TOTPRODdecay[vt].carbon / (double) CYCLE;
  pdstate[vt][I_TOTPRDFN] = ag.TOTPRODdecay[vt].nitrogen / (double) CYCLE;

  // Integrated carbon fluxes

  pdstate[vt][I_TOTNPP] = veg.npp[vt][dm] + ag.npp[vt][dm].carbon;
  pdstate[vt][I_CFLX] = cflux[vt][dm];

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

//void TTEM::deltaxclm(const int& dcmnt, const double& pcfldcap, const int& dm)
void TTEM::deltaxclm(const int& dcmnt, const double& pcfldcap, const int& dm, const int& vt)
{

  double vfc;


  vfc = 93 * 0.01;


/* gv: effect of moisture on primary productivity */

  if (moistlim == 0)
  {
    gv[vt] = 1.0;
  }
  else
  {
    if (atms.eet[vt][dm]/atms.pet[vt][dm] <= 0.1)
    {
      gv[vt] = (-10.0 * pow((atms.eet[vt][dm]/atms.pet[vt][dm]),2.0))
           + (2.9 * (atms.eet[vt][dm]/atms.pet[vt][dm]));
      if (gv[vt] < 0.0) { gv[vt] = 0.0; }
    }
    else
    {
      gv[vt] = 0.1 + (0.9 * atms.eet[vt][dm] / atms.pet[vt][dm]);
    }
  }

  //printf("%f     %f   %f      %f\n",atms.tair[vt][dm],atms.eet[vt][dm],atms.pet[vt][dm],gv[vt]);  
/* ksoil: effect of soil moisture on nitrogen uptake by plants
			 and microbes */

  if (moistlim == 0) { ksoil[vt] = pow(vfc,3.0); }                                                                                                  //wsr: change VSM to seperate VSMs
  else { ksoil[vt] = pow((y[vt][I_VSM1]+y[vt][I_VSM2])/2,3.0); }


/* rhmoist: effect of moisture on decomposition */

// modified for 3-box hydrology, replace y[I_VSM] with y[I_VSM2], Oct 11 2000

  if (moistlim == 0)
  {
    rhmoist[vt] = (vfc - microbe.moistmin[dcmnt])
              * (vfc - microbe.moistmax[dcmnt]);
    rhmoist[vt] /= rhmoist[vt] - pow((vfc - microbe.moistopt[dcmnt]),2.0);
  }
  else
  {
    rhmoist[vt] = ((y[vt][I_VSM1]*1/3+y[vt][I_VSM2]*2/3) - microbe.moistmin[dcmnt])                          //wsr: use VSM1 to 2 instead of VSM
              * ((y[vt][I_VSM1]*1/3+y[vt][I_VSM2]*2/3)  - microbe.moistmax[dcmnt]);
//    rhmoist /= rhmoist - pow((y[I_VSM] - microbe.moistopt[dcmnt]),2.0);
    rhmoist[vt] /= rhmoist[vt] - pow(((y[vt][I_VSM1]*1/3+y[vt][I_VSM2]*2/3)- microbe.moistopt[dcmnt]),2.0);
  }
  if (rhmoist[vt] < 0.0) { rhmoist[vt] = 0.0; }

};

/* *************************************************************
************************************************************* */



/* *************************************************************
************************************************************* */

int TTEM::ecdqc(const int& dcmnt)
{

  int qc = ACCEPT;

  if (vegca[dcmnt] <= -999.99) { return qc = REJECT; }
  if (vegcb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (strna[dcmnt] <= -999.99) { return qc = REJECT; }
  if (strnb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (solca[dcmnt] <= -999.99) { return qc = REJECT; }
  if (solcb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (solna[dcmnt] <= -999.99) { return qc = REJECT; }
  if (solnb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (avlna[dcmnt] <= -999.99) { return qc = REJECT; }
  if (avlnb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (stona[dcmnt] <= -999.99) { return qc = REJECT; }
  if (stonb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.unleaf12[dcmnt] <= -9.99) { return qc = REJECT; }
  if (veg.prvleafmx[dcmnt] <= -9.99) { return qc = REJECT; }
  if (veg.cmaxcut[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.cmax1a[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.cmax1b[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.cmax2a[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.cmax2b[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.cfall[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.kra[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.krb[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.kda[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.kdb[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.lcclnc[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.propftos[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.nmaxcut[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.nmax1a[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.nmax1b[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.nmax2a[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.nmax2b[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.nfall[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.nupa[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.nupb[dcmnt] <= -99.99) { return qc = REJECT; }
  if (soil.nloss[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.nfixpar[dcmnt] <= -99.9) { return qc = REJECT; }
  if (veg.cneven[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.cnmin[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.c2na[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.c2nb[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.c2nmin[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.cnsoil[dcmnt] <= -999.99) { return qc = REJECT; }

  return qc;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

//void TTEM::ECDsetELMNTstate(const int& dcmnt, const double& psiplusc)
void TTEM::ECDsetELMNTstate(const int& dcmnt, const double& psiplusc, const int& vt)
{
	
  int dyr;
  int dm;

  for (dm = 0; dm < CYCLE; dm++)
  {
    y[vt][I_AVLW] = soil.avlh2o[vt][dm] = soil.awcapmm[vt];
    y[vt][I_RGRW] = soil.rgrndh2o[vt][dm] = 0.0;
    y[vt][I_SNWPCK] = soil.snowpack[vt][dm] = 0.0;
    y[vt][I_SGRW] = soil.sgrndh2o[vt][dm] = 0.0;
    y[vt][I_SM] = soil.moist[vt][dm] = soil.awcapmm[vt] + soil.wiltpt[vt];
    y[vt][I_PCTP] = soil.pctp[vt][dm] = 100.0 * soil.moist[vt][dm] / soil.totpor[vt];

   
  
  
  soil.vsm[vt][dm] = soil.moist[vt][dm] / (soil.rootz[vt] * 1000.0);
    if (soil.vsm[vt][dm] <= 0.0)
   {
     soil.vsm[vt][dm] = 0.001;
    }

 soil.vsm1[vt][dm] = (soil.moist1[vt][dm]) / (soil.dpwbox1[vt]*1000.0);
 if (soil.vsm1[vt][dm] <= 0.0) {
   soil.vsm1[vt][dm] = 0.001;
   }

//   soil.vsm2[dm] = soil.moist2[dm] / (soil.rootmx*1000.0-soil.dpwbox1*1000.0);

   soil.vsm2[vt][dm] = soil.moist2[vt][dm] / (soil.dpwbox2[vt]*1000.0);
	 if (soil.vsm2[vt][dm] <= 0.0) {
	 soil.vsm2[vt][dm] = 0.001;
  }

//   soil.vsm3[dm] = soil.moist3[dm] / (soil.rootmx*1000.0-soil.dpwbox1*1000.0); // temporally set up as layer 2
   soil.vsm3[vt][dm] = soil.moist3[vt][dm] / (soil.dpwbox3[vt]*1000.0); // temporally set up as layer 2
	 if (soil.vsm3[vt][dm] <= 0.0) {
	 soil.vsm3[vt][dm] = 0.001;
  }
  y[vt][I_VSM] = soil.vsm[vt][dm];
  y[vt][I_VSM1] = soil.vsm1[vt][dm];
  y[vt][I_VSM2] = soil.vsm2[vt][dm];
  y[vt][I_VSM3] = soil.vsm3[vt][dm];

  y[vt][I_AVLW1] = soil.avlh2o1[vt][dm] = soil.awcapmm1[vt];
  y[vt][I_AVLW2] = soil.avlh2o2[vt][dm] = soil.awcapmm2[vt];
  y[vt][I_AVLW3] = soil.avlh2o3[vt][dm] = soil.awcapmm3[vt];


  y[vt][I_SM1] = y[vt][I_MOIST1] = soil.moist1[vt][dm] = (soil.awcapmm1[vt] + soil.wiltpt1[vt]);
  y[vt][I_SM2] = y[vt][I_MOIST2] = soil.moist2[vt][dm] = (soil.awcapmm2[vt] + soil.wiltpt2[vt]);
  y[vt][I_SM3] = y[vt][I_MOIST3] = soil.moist3[vt][dm] = (soil.awcapmm3[vt] + soil.wiltpt3[vt]);

  y[vt][I_PCTP1] = soil.pctp1[vt][dm] = 100.0 * soil.moist1[vt][dm] / soil.totpor1[vt];
  y[vt][I_PCTP2] = soil.pctp2[vt][dm] = 100.0 * soil.moist2[vt][dm] / soil.totpor2[vt];
  y[vt][I_PCTP3] = soil.pctp3[vt][dm] = 100.0 * soil.moist3[vt][dm] / soil.totpor3[vt];




//  hydm.alpha=0.4118/0.5882;
//  hydm.sapallo=0.0;
//  hydm.sla= 0.0072;
//  hydm.cov=0.5;
//  hydm.fpcmax = 0.355;

//  hydm.rtratio1=0.00;  // forest root distributed in second and third layers
//  hydm.rtratio2=0.80;
//  hydm.rtratio3=0.20;

//  hydm.rteffcnt=0.024;    // root surface extraction efficiency
//  hydm.kwt1= 0.25;
//  hydm.kwt2= 0.25;
//  hydm.kwt3= 0.25;
//  hydm.wscal=1.0;
  //hydm.dense=0.5;
//  veg.plant[dm].carbon = vegca[dcmnt] * psiplusc + vegcb[dcmnt];
//  hydm.rootmass[dm] = 0.5*veg.plant[dm].carbon;

//  hydm.rtsface = hydm.rteffcnt*hydm.rootmass[dm];

//  hydm.rtindex1= hydm.rtratio1*hydm.rtsface;
//  hydm.rtindex2= hydm.rtratio2*hydm.rtsface;
//  hydm.rtindex3= hydm.rtratio3*hydm.rtsface;

//  hydm.root1= 1.0-exp(-0.5*hydm.rtindex1+hydm.kwt1);
//  hydm.root2= 1.0-exp(-0.5*hydm.rtindex2+hydm.kwt2);
//  hydm.root3= 1.0-exp(-0.5*hydm.rtindex3+hydm.kwt3);

  hydm.root1[vt]= 0.1;   // assume 1% root at upper layer             based on Haxeltine, 1996
  hydm.root2[vt]= 0.8;   // asume 90% root at second root layer
  hydm.root3[vt]= 0.01;   // asume 90% root at second root layer

  hydm.wscal=1.0; // hydm

 // *************** end of changing

 //   veg.plant[dm].carbon = vegca[dcmnt] * psiplusc + vegcb[dcmnt];

    y[vt][I_VEGC] = veg.plant[vt][dm].carbon;

    veg.strctrl[vt][dm].nitrogen = strna[dcmnt] * psiplusc + strnb[dcmnt];
    y[vt][I_STRN] = veg.strctrl[vt][dm].nitrogen;

    soil.org[vt][dm].carbon = solca[dcmnt] * psiplusc + solcb[dcmnt];
    y[vt][I_SOLC] = soil.org[vt][dm].carbon;

    soil.org[vt][dm].nitrogen = solna[dcmnt] * psiplusc + solnb[dcmnt];
    y[vt][I_SOLN] = soil.org[vt][dm].nitrogen;

    soil.availn[vt][dm] = avlna[dcmnt] * psiplusc + avlnb[dcmnt];
    y[vt][I_AVLN] = soil.availn[vt][dm];

    veg.labile[vt][dm].nitrogen = stona[dcmnt] * psiplusc + stonb[dcmnt];
    y[vt][I_STON] = veg.labile[vt][dm].nitrogen;

    y[vt][I_UNRMLF] = veg.unnormleaf[vt][dm] = veg.unleaf12[dcmnt];

    veg.plant[vt][dm].nitrogen = 0.0;

    veg.alleaf[vt] = veg.leafmxc[dcmnt]/(1.0 + veg.kleafc[dcmnt]
                 * exp(veg.cov[dcmnt]*y[vt][I_VEGC]));
    veg.lai[vt][dm] = veg.sla[dcmnt] * veg.alleaf[vt];
    veg.foliage[vt] = veg.alleaf[vt] / veg.leafmxc[dcmnt];
    veg.fpc[vt][dm] = 1.0;

    y[vt][I_PROD10C] = ag.PROD10[vt].carbon = 0.0;
    y[vt][I_PROD10N] = ag.PROD10[vt].nitrogen = 0.0;
    y[vt][I_PROD100C] = ag.PROD100[vt].carbon = 0.0;
    y[vt][I_PROD100N] = ag.PROD100[vt].nitrogen = 0.0;
    y[vt][I_AGPRDC] = ag.PROD1[vt].carbon = 0.0;
    y[vt][I_AGPRDN] = ag.PROD1[vt].nitrogen = 0.0;
    y[vt][I_TOTPRDC] = ag.TOTPROD[vt].carbon = 0.0;
    y[vt][I_TOTPRDN] = ag.TOTPROD[vt].nitrogen = 0.0;
    y[vt][I_TOTEC] = totalc[vt][dm] = veg.plant[vt][dm].carbon + soil.org[vt][dm].carbon;
    y[vt][I_TOTGC] = totalc[vt][dm];
  }

  for (dyr = 0; dyr < 10; dyr++)
  {
    ag.initPROD10[vt][dyr].carbon = 0.0;
    ag.initPROD10[vt][dyr].nitrogen = 0.0;
  }

  for (dyr = 0; dyr < 100; dyr++)
  {
    ag.initPROD100[vt][dyr].carbon = 0.0;
    ag.initPROD100[vt][dyr].nitrogen = 0.0;
  }
  //soil.tolfreeze[vt][dm] = 0.0;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */
//int TTEM::equilibrium(const int& itype, double& tol)
int TTEM::equilibrium(const int& itype, double& tol, const int& vt)
{

  int dyr = 0;

  setPrevState(prevy,y,vt);

  totyr = 0;
  endeq = 0;
  intflag = 0;
  initFlag = 0;
  while ((dyr < runsize) && (endeq < 2))
  {


    endeq = stepyr(dyr,itype,intflag,tol,vt);

    ++dyr;
    ++totyr;

 // Reset product fluxes and pools to zero

    ag.resetPROD(vt);

// Check to see if steady state conditions have been reached.

    if (dyr >= strteq && endeq == 0)
    {

      if (nfeed == 0 && rheqflag == 0
         
         && ((ctol >= fabs(veg.yrnpp[vt] - veg.yrltrc[vt]))
         || ( 0.001 >= fabs(veg.yrnpp[vt] - veg.yrltrc[vt]))))
      {
	     endeq = 1;
      }
      if (nfeed == 0 && rheqflag == 1 												//wsr: reduce one-box equilibrium also for no nfeedback
         && (ctol >= fabs(yrnep[vt]))
         && (ctol >= fabs(veg.yrnpp[vt] - veg.yrltrc[vt]))
         && (ctol >= fabs(veg.yrltrc[vt] - microbe.yrrh[vt])))
      {
	     endeq = 1;
      }
   /*   if (nfeed == 1 && rheqflag == 1 && (wtol >= fabs(atms.yrrain[vt]
         + atms.yrsnowfall[vt] - atms.yreet[vt] - soil.yrh2oyld[vt]))
         && (ntol >= fabs(soil.yrnin[vt] - soil.yrnlost[vt]))
         && (ntol >= fabs(veg.yrnup[vt] - veg.yrltrn[vt]))
         && (ntol >= fabs(veg.yrnup[vt] - microbe.yrnmin[vt]))
         && (ntol >= fabs(veg.yrltrn[vt] - microbe.yrnmin[vt]))
         && (ctol >= fabs(yrnep[vt])) && (ctol >= fabs(veg.yrnpp[vt] - veg.yrltrc[vt]))
         && (ctol >= fabs(veg.yrltrc[vt] - microbe.yrrh[vt])))                                                                   // wsr: change equilibrium for water to debug
         																						//wsr: reduce one-box water equilibrium
         
         */
         
         if (nfeed == 1 && rheqflag == 1 
         && (ntol >= fabs(soil.yrnin[vt] - soil.yrnlost[vt]))
         && (ntol >= fabs(veg.yrnup[vt] - veg.yrltrn[vt]))
         && (ntol >= fabs(veg.yrnup[vt] - microbe.yrnmin[vt]))
         && (ntol >= fabs(veg.yrltrn[vt] - microbe.yrnmin[vt]))
         && (ctol >= fabs(yrnep[vt])) && (ctol >= fabs(veg.yrnpp[vt] - veg.yrltrc[vt]))
         && (ctol >= fabs(veg.yrltrc[vt] - microbe.yrrh[vt])))
         
         
         
         
         
      {
        endeq = 1;
      }
    }
  }

  if (endeq == 2)
  {
    nattempt = maxnrun;
    initFlag = 1;
  }

  if (dyr >= runsize && endeq < 2) { ++nattempt; }

  return nattempt;

};
/* *************************************************************
************************************************************* */

// modified dfor hydrology model by Qz
/* *************************************************************
************************************************************* */

//void TTEM::getenviron(const int& dm)

void TTEM::getenviron(const int& dm, const int& vt)
{
	//y[vt][I_AVLW2] = soil.freeze(y[vt][I_AVLW2],atms.tair[vt][dm],vt,dm);

	
	//printf("%f   %f   %f     %f       %f     %f    %f    %f   %f    %f   %f\n",soil.snowinf[vt][dm],atms.rain[vt][dm],hyd.snowsub[vt][dm],soil.vsm2[vt][dm],hydm.lyperc1[vt][dm],hydm.lyperc2[vt][dm],soil.storage2[vt][dm],soil.storage3[vt][dm],atms.tair[vt][dm],soil.avlh2o2[vt][dm],soil.avlh2o3[vt][dm]);

	soil.wetlandTAWA(soil.vsm1[vt][dm],soil.vsm2[vt][dm],soil.dpwbox1[vt],soil.dpwbox2[vt],soil.totpor1[vt],soil.totpor2[vt],vt,dm);


	//printf("   %f    \n",soil.vsm2[vt][dm]);

  // Determine monthly potential evapotranspiration
   // Determine if monthly precipitation occurs as rain or snow
  atms.precsplt(atms.prec[vt][dm],atms.tair[vt][dm],atms.rain[vt][dm],atms.snowfall[vt][dm]);

  soil.avlh2o2[vt][dm] =  y[vt][I_AVLW2];    // using the second layer which is organic layer available water to influence the leaf water potential


   //veg.lai[vt][dm] = 2.31; // for calibration version

 //  hyd.hvpd[dm] = 10; // 20gm-3 will result in the close of stomatal need to reconsider

 // Determine monthly potential evapotranspiration
//    hyd.vap[vt][dm] = atms.vap[vt][dm]; // assign static vapor pressure to hyd.vap[dm] instead of using transient data, need to figure out a way to deal with this conflict

  hyd.vaporpressure[vt][dm] = hyd.vap[vt][dm];

  hyd.vpd(hyd.vaporpressure[vt][dm], atms.tair[vt][dm], vt);
  hyd.intercept(atms.rain[vt][dm], hyd.inter_coef[veg.cmnt],dm, vt);
  hyd.canopy_water[vt][dm] =hyd.prec_to_canopy[vt][dm];
  hyd.drad(veg.lai[vt][dm], hyd.EXT[veg.cmnt], atms.nirr[vt][dm], dm, vt);
//  hyd.hLWP[dm] =  hyd.LWP(soil.avlh2o[dm], soil.awcapmm);

  hyd.hLWP[vt][dm] =  hyd.LWP(soil.avlh2o2[vt][dm],soil.fldcap[vt], soil.pctsand, soil.pctclay, soil.pctsilt);    							//wsr

  hyd.hcch[vt][dm] = hyd.CanopyCond(hyd.LWPmin[veg.cmnt], hyd.hLWP[vt][dm],hyd.CCmax[veg.cmnt],hyd.LWPc[veg.cmnt],hyd.absoluteHD[vt],hyd.SLOPEcc[veg.cmnt]);

  hyd.hslope[vt][dm] = hyd.slope(atms.tair[vt][dm]);
  hyd.hpa[vt][dm] = hyd.densityA(atms.tair[vt][dm]);
  hyd.hlv[vt][dm] = hyd.latentV(atms.tair[vt][dm]);

  //rainfall is changed to
  atms.rain[vt][dm] = hyd.rainthrough[vt][dm];

  if (hyd.canopy_water[vt][dm] >0)
   {
     hyd.htrans[vt][dm] = hyd.penmanmonteith(hyd.hslope[vt][dm], hyd.hdrad[vt][dm], hyd.hpa[vt][dm],hyd.vpdeficit[vt],
                    hyd.hcch[vt][dm],hyd.hlv[vt][dm]) * veg.lai[vt][dm];
     hyd.hevap[vt][dm] = hyd.evaporation(hyd.canopy_water[vt][dm]);
     if ((hyd.hevap[vt][dm] / hyd.htrans[vt][dm]) < hyd.daze[dm]) {
              hyd.htrans[vt][dm] = hyd.htrans[vt][dm] * (hyd.daze[dm] - ceil( hyd.hevap[vt][dm]/ hyd.htrans[vt][dm]));
                            }
     else { hyd.htrans[vt][dm] = 0.0; }
    }
   else // there is no canopy water
     {
    hyd.htrans[vt][dm] = hyd.penmanmonteith(hyd.hslope[vt][dm], hyd.hdrad[vt][dm], hyd.hpa[vt][dm],hyd.vpdeficit[vt],
                   hyd.hcch[vt][dm],hyd.hlv[vt][dm]) * veg.lai[vt][dm];
    hyd.hevap[vt][dm]=0.0;
    hyd.htrans[vt][dm] = hyd.htrans[vt][dm] * hyd.daze[dm];
         }

  // Evaporation from soil surface

	// assign days since rain
 	hyd.dayssincerain = 10;

// first calculate potential evaporation, assuming the resistance
 // for vapor transport is equal to the resistance for sensible heat
 //	transport.  That is, no additional resistance for vapor transport to
 //	the soil surface. This represents evaporation from a wet surface with
 //	a specified aerodynamic resistance (= boundary layer resistance).
 //	The aerodynamic resistance is for now set as a constant, and is
 //	taken from observations over bare soil in tiger-bush in south-west
 //	Niger: aero_resis_layer1 = 107 s m-1 (Wallace and Holwill, 1997).


 // calculate potential_evap in mm/m2/day
   hyd.potevap_surface[vt][dm] = hyd.evaporation_layer1(atms.tair[vt][dm], hyd.hpa[vt][dm],hyd.hlv[vt][dm],hyd.vpdeficit[vt], hyd.surfacerad[vt][dm])*hyd.daze[dm];

  // consider only the precipitation flux reaching the soil
 // check for precipitation >= potential evaporation
	if (hyd.rainthrough[vt][dm] >= hyd.potevap_surface[vt][dm])
	{
		// reset days-since-rain parameter
		hyd.dayssincerain = 0.0;

		// soil evaporation proceeds at potential rate
		hyd.soil_evap[vt][dm] = 0.6 * hyd.potevap_surface[vt][dm];
	}
	else
	{
		// increment the days since rain
		hyd.dayssincerain += 1.0;

		// calculate the realized proportion of potential evaporation as a function of the days since rain
		hyd.ratio = 0.3/pow(hyd.dayssincerain,2.0);

		// calculate evaporation for dry days
		hyd.soil_evap[vt][dm] = hyd.ratio * hyd.potevap_surface[vt][dm];

		// for rain events that are smaller than required to reset dayssincerain counter, but larger than dry-day evaporation, all rain is evaporated.
		//In this case, do not advance the drying curve counter.
		//For rain events that are too small to trigger dsr reset, and which
		//are smaller than dry-day evap, there will be more evaporation than
		//rainfall.  In this case the drying curve counter is advanced.
		if (hyd.rainthrough[vt][dm] > hyd.soil_evap[vt][dm])
		{
			hyd.soil_evap[vt][dm] = hyd.rainthrough[vt][dm];
			hyd.dayssincerain -= 1.0;
		}
	}
// available water is changed to
   //soil.avlh2o[vt][dm] =  y[vt][I_AVLW] - hyd.soil_evap[vt][dm];    // ???? need to think

  // end of evaporation from soil surface

    //sublimation from canopy
    //snow intercept
   if (soil.snowpack[vt][dm] > 0.0) {
      hyd.snowpack_canopy[vt][dm] = 0.01 * soil.snowpack[vt][dm]; // 0.5mm/LAI is intercept coefficient,

      hyd.sub_from_canopy[vt][dm] = hyd.sublimation_canopy(veg.lai[vt][dm],hyd.hdrad[vt][dm], atms.tair[vt][dm]);

	  if (hyd.sub_from_canopy[vt][dm] ==0.0)     {  hyd.sub_from_canopy[vt][dm] = 0.001; }// to protect

      hyd.sublimationdays[vt][dm] = hyd.snowpack_canopy[vt][dm] / hyd.sub_from_canopy[vt][dm];
      if (hyd.sublimationdays[vt][dm] > hyd.daze[dm])
	  {  hyd.sub_from_canopy[vt][dm] = hyd.snowpack_canopy[vt][dm]; }
	  else {hyd.sub_from_canopy[vt][dm] = hyd.sub_from_canopy[vt][dm] * hyd.sublimationdays[vt][dm];}
    }

   // snowpack is changed to
    soil.snowpack[vt][dm] = (soil.snowpack[vt][dm] - hyd.sub_from_canopy[vt][dm]);

  // sublimation from ground


  hyd.snowsub[vt][dm] = hyd.snowsublimation_ground(soil.snowpack[vt][dm], atms.tair[vt][dm], hyd.surfacerad[vt][dm]);

  atms.pet[vt][dm] = hyd.htrans[vt][dm] + hyd.hevap[vt][dm] + hyd.soil_evap[vt][dm]+ hyd.snowsub[vt][dm] + hyd.sub_from_canopy[vt][dm];
//  atms.pet[dm] = hyd.htrans[dm] + hyd.hevap[dm] + hyd.soil_evap[dm] + hyd.sub_from_canopy[dm];

printf("%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f\n",atms.pet[vt][dm],hyd.htrans[vt][dm],hyd.hevap[vt][dm],hyd.soil_evap[vt][dm],hyd.snowsub[vt][dm],hyd.sub_from_canopy[vt][dm],atms.eet[vt][dm],hyd.hslope[vt][dm],hyd.hdrad[vt][dm],hyd.hpa[vt][dm],hyd.vpdeficit[vt],hyd.hcch[vt][dm],hyd.hLWP[vt][dm],hyd.absoluteHD[vt]);
   // end of modify
 // assign to atms.pet[dm]

//  atms.pet[dm] = atms.petjh(atms.nirr[dm], atms.tair[dm], dm);

  if (atms.pet[vt][dm] <= 0.0) { atms.pet[vt][dm] = 0.001; }

  // Determine if monthly precipitation occurs as rain or snow

  atms.precsplt(atms.prec[vt][dm],atms.tair[vt][dm],atms.rain[vt][dm],atms.snowfall[vt][dm]);

  // Determine previous two month's air temperatures for following
  //   snowmelt calculations

  switch (dm)
  {
    case 0:  atms.prevtair = atms.tair[vt][CYCLE-1];
             atms.prev2tair = atms.tair[vt][CYCLE-2];
             break;
    case 1:  atms.prevtair = atms.tair[vt][0];
	     atms.prev2tair = atms.tair[vt][CYCLE-1];
	     break;
    default: atms.prevtair = atms.tair[vt][dm-1];
	     atms.prev2tair = atms.tair[vt][dm-2];
	     break;
  }

  // Determine contribution of snowmelt to soil moisture


  soil.snowinf[vt][dm] = soil.snowmelt(elev, atms.tair[vt][dm],
                                   atms.prevtair, y[vt][I_SNWPCK]);

  monthxclm(veg.cmnt, veg.topt[vt], dm, vt);


  	
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::getsitecd(ofstream& rflog1)
{

  char ecd[80];

  cout << "Enter name of the site (.ECD) data file with the parameter values:";
  cout << endl;
  cin >> ecd;

  rflog1 << "Enter name of the site (.ECD) data file with the parameter values:";
  rflog1 << ecd << endl << endl;

  getsitecd(ecd);

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::getsitecd(const int& numcmnt, ofstream& rflog1)
{

  int dv;
  char ecd[80];

  cout << "Enter name of the site (.ECD) data file with the parameter values  for cmnt" << endl;
  rflog1 << "Enter name of the site (.ECD) data file with the parameter values cmnt" << endl;
  for (dv = 0; dv < numcmnt; dv++)
  {
    cout << (dv+1) << ": ";
    cin >> ecd;

    rflog1 << (dv+1) << ": " << ecd << endl;

    getsitecd(dv, ecd);
  }
  rflog1 << endl;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::getsitecd(char ecd[80])
{

  const int NUMVAR = 52;
  char dummy[NUMVAR][10];
  ifstream infile;
  int i;
  int dcmnt;

  int sitetveg[MAXCMNT];
  int sitewsoil[MAXCMNT];
  long update[MAXCMNT];
  char sitevegtype[MAXCMNT][31];
  char sitename[MAXCMNT][17];
  char sitetext[MAXCMNT][14];
  float sitecol[MAXCMNT];
  float siterow[MAXCMNT];

  infile.open(ecd, ios::in);

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 0; dcmnt < MAXCMNT; dcmnt++)
  {
    infile >> sitetveg[dcmnt] >> sitevegtype[dcmnt];
    infile >> sitecol[dcmnt] >> siterow[dcmnt];
    infile >> sitename[dcmnt] >> sitetext[dcmnt] >> sitewsoil[dcmnt];
    infile >> vegca[dcmnt] >> vegcb[dcmnt];
    infile >> strna[dcmnt] >> strnb[dcmnt];
    infile >> solca[dcmnt] >> solcb[dcmnt];
    infile >> solna[dcmnt] >> solnb[dcmnt];
    infile >> avlna[dcmnt] >> avlnb[dcmnt];
    infile >> stona[dcmnt] >> stonb[dcmnt];
    infile >> veg.unleaf12[dcmnt];
    infile >> veg.initleafmx[dcmnt];      // Changed by DWK on 19991028
    infile >> veg.cmaxcut[dcmnt];
    infile >> veg.cmax1a[dcmnt] >> veg.cmax1b[dcmnt];
    infile >> veg.cmax2a[dcmnt] >> veg.cmax2b[dcmnt];
    infile >> veg.cfall[dcmnt];
    infile >> veg.kra[dcmnt] >> veg.krb[dcmnt];
    infile >> microbe.kda[dcmnt] >> microbe.kdb[dcmnt];
    infile >> microbe.lcclnc[dcmnt] >> microbe.propftos[dcmnt];
    infile >> veg.nmaxcut[dcmnt];
    infile >> veg.nmax1a[dcmnt] >> veg.nmax1b[dcmnt];
    infile >> veg.nmax2a[dcmnt] >> veg.nmax2b[dcmnt];
    infile >> veg.nfall[dcmnt];
    infile >> microbe.nupa[dcmnt] >> microbe.nupb[dcmnt];
    infile >> soil.nloss[dcmnt];
    infile >> microbe.nfixpar[dcmnt];
    infile >> veg.initcneven[dcmnt] >> veg.cnmin[dcmnt];
    infile >> veg.c2na[dcmnt] >> veg.c2nb[dcmnt] >> veg.c2nmin[dcmnt];
    infile >> microbe.cnsoil[dcmnt];
    infile >> update[dcmnt];

//    veg.initcneven[i] = veg.cneven[i];
//    veg.adjc2n = 1.0 + (veg.dc2n * (atms.co2[11] - atms.initco2));
//    veg.cneven[i] = veg.initcneven[i] * veg.adjc2n;
  }

  infile.close();

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::getsitecd(const int& dv, char ecd[80])
{
  // Function added by DWK on 20000102

  char dummy[12];
  // string changed to char[80] by DWK on 20000210
  char sitename[80];
  float sitecol;
  float siterow;
  long updated;

  int vt; // local use for moss component

  fecd[dv].open(ecd,ios::in);

  if (!fecd[dv])
  {
    cerr << endl << "Cannot open " << ecd << " for data input" << endl;
    exit(-1);
  }


for (vt=0; vt < NVT; vt++) // read tree first and then moss by Q.Z.
{

// Addition by Q. Z.

  fecd[dv] >> dummy >> veg.cmnt;
  subsist[vt]= veg.cmnt;
  veg.cmnt= subsist[vt];

  fecd[dv] >> dummy >> veg.cmnt_name;
  fecd[dv] >> dummy >> sitename;
  fecd[dv] >> dummy >> sitecol;
  fecd[dv] >> dummy >> siterow;
  fecd[dv] >> dummy >> updated;

  fecd[dv] >> dummy >> vegca[veg.cmnt];
  fecd[dv] >> dummy >> vegcb[veg.cmnt];
  fecd[dv] >> dummy >> strna[veg.cmnt];
  fecd[dv] >> dummy >> strnb[veg.cmnt];
  fecd[dv] >> dummy >> solca[veg.cmnt];
  fecd[dv] >> dummy >> solcb[veg.cmnt];
  fecd[dv] >> dummy >> solna[veg.cmnt];
  fecd[dv] >> dummy >> solnb[veg.cmnt];
  fecd[dv] >> dummy >> avlna[veg.cmnt];
  fecd[dv] >> dummy >> avlnb[veg.cmnt];
  fecd[dv] >> dummy >> stona[veg.cmnt];
  fecd[dv] >> dummy >> stonb[veg.cmnt];

  fecd[dv] >> dummy >> veg.unleaf12[veg.cmnt];
  // veg.prevleafmx changed to veg.initleafmx by DWK on 20000130
  fecd[dv] >> dummy >> veg.initleafmx[veg.cmnt];
  fecd[dv] >> dummy >> veg.cmaxcut[veg.cmnt];
  fecd[dv] >> dummy >> veg.cmax1a[veg.cmnt];
  fecd[dv] >> dummy >> veg.cmax1b[veg.cmnt];
  fecd[dv] >> dummy >> veg.cmax2a[veg.cmnt];
  fecd[dv] >> dummy >> veg.cmax2b[veg.cmnt];
  fecd[dv] >> dummy >> veg.cfall[veg.cmnt];
  fecd[dv] >> dummy >> veg.kra[veg.cmnt];
  fecd[dv] >> dummy >> veg.krb[veg.cmnt];
  fecd[dv] >> dummy >> microbe.kda[veg.cmnt];
  fecd[dv] >> dummy >> microbe.kdb[veg.cmnt];
  fecd[dv] >> dummy >> microbe.lcclnc[veg.cmnt];
  fecd[dv] >> dummy >> microbe.propftos[veg.cmnt];
  fecd[dv] >> dummy >> veg.nmaxcut[veg.cmnt];
  fecd[dv] >> dummy >> veg.nmax1a[veg.cmnt];
  fecd[dv] >> dummy >> veg.nmax1b[veg.cmnt];
  fecd[dv] >> dummy >> veg.nmax2a[veg.cmnt];
  fecd[dv] >> dummy >> veg.nmax2b[veg.cmnt];
  fecd[dv] >> dummy >> veg.nfall[veg.cmnt];
  fecd[dv] >> dummy >> microbe.nupa[veg.cmnt];
  fecd[dv] >> dummy >> microbe.nupb[veg.cmnt];
  fecd[dv] >> dummy >> soil.nloss[veg.cmnt];
  fecd[dv] >> dummy >> microbe.nfixpar[veg.cmnt];
  fecd[dv] >> dummy >> veg.initcneven[veg.cmnt];
  fecd[dv] >> dummy >> veg.cnmin[veg.cmnt];
  fecd[dv] >> dummy >> veg.c2na[veg.cmnt];
  fecd[dv] >> dummy >> veg.c2nb[veg.cmnt];
  fecd[dv] >> dummy >> veg.c2nmin[veg.cmnt];
  fecd[dv] >> dummy >> microbe.cnsoil[veg.cmnt];

 }; // veg.cmnt is changed into last time value (which is moss type)

  fecd[dv].close();

};

/* *************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void TTEM::initrun(ofstream& rflog1, const int& equil)
{

  avlnflag = nfeed = rheqflag = 0;

/* **************************************************************
		  Run Model with Nitrogen Limitation?
************************************************************** */

  cout << endl << "Do you want to allow available N to fluctuate?" << endl;
  cout << "  Enter 0 for No" << endl;
  cout << "  Enter 1 for Yes: ";
  cin >> avlnflag;

  rflog1 << endl << "Do you want to allow available N to fluctuate?" << endl;
  rflog1 << "  Enter 0 for No" << endl;
  rflog1 << "  Enter 1 for Yes: " << endl;
  rflog1 << "avlnflag = " << avlnflag << endl << endl;

  cout << endl << "Do you want nitrogen feedback on GPP?" << endl;
  cout << "  Enter 0 for No" << endl;
  cout << "  Enter 1 for Yes: ";
  cin >> nfeed;

  rflog1 << endl << "Do you want nitrogen feedback on GPP?" << endl;
  rflog1 << "  Enter 0 for No" << endl;
  rflog1 << "  Enter 1 for Yes: " << endl;
  rflog1 << "nfeed = " << nfeed << endl << endl;

  baseline = initbase = 0;
  if (nfeed == 1)
  {
    cout << endl << "Do you want to solve for baseline soil nitrogen?" << endl;
    cout << "  Enter 0 for No" << endl;
    cout << "  Enter 1 for Yes: ";
    cin >> initbase;
    baseline = initbase;

    rflog1 << endl << "Do you want to solve for baseline soil nitrogen?" << endl;
    rflog1 << "  Enter 0 for No" << endl;
    rflog1 << "  Enter 1 for Yes: " << endl;
    rflog1 << "baseline = " << baseline << endl << endl;
  }

/* **************************************************************
			 Run Model with Moisture Limitation?
************************************************************** */

  moistlim = 0;
  cout << endl << "Do you want to run the model with moisture limitation?" << endl;
  cout << "  Enter 0 for No" << endl;
  cout << "  Enter 1 for Yes: ";
  cin >> moistlim;

  rflog1 << endl << "Do you want to run the model with moisture limitation?" << endl;
  rflog1 << "  Enter 0 for No" << endl;
  rflog1 << "  Enter 1 for Yes: " << endl;
  rflog1 << "moistlim = " << moistlim << endl << endl;


/* ***************************************************************
	       Details for Steady State Conditions
************************************************************** */


  maxyears = 0;
  maxnrun = 0;

  cout << endl << "How many years do you want to wait before checking equilibrium conditions? ";
  cin >> strteq;

  rflog1 << endl;
  rflog1 << "How many years do you want to wait before checking equilibrium conditions? ";
  rflog1 << endl;
  rflog1 << "strteq = " << strteq << endl << endl;

  cout << endl << "Enter the maximum number of years for the model to run: ";
  cin >> maxyears;

  rflog1 << endl << "Enter the maximum number of years for the model to run: ";
  rflog1 << endl;
  rflog1 << "maxyears = " << maxyears << endl << endl;

  runsize = maxyears;

  cout << endl << "Enter the maximum number of attempts to reach a solution: ";
  cin >> maxnrun;

  rflog1 << endl;
  rflog1 << "Enter the maximum number of attempts to reach a solution: ";
  rflog1 << endl;
  rflog1 << "maxnrun = " << maxnrun << endl << endl;

  if (nfeed == 0)
  {
    cout << endl << "Do you want decomposition to come into equilibrium? ";
    cout << "  Enter 0 for No" << endl;
    cout << "  Enter 1 for Yes: ";
    cin >> rheqflag;

    rflog1 << endl;
    rflog1 << "Do you want decomposition to come into equilibrium? " << endl;
    rflog1 << "  Enter 0 for No" << endl;
    rflog1 << "  Enter 1 for Yes: " << endl;
    rflog1 << "rheqflag = " << rheqflag << endl << endl;
  }

  wtol = 1000.0;
  cout << endl;
  cout << "What absolute tolerance do you want to use for checking equilibrium";
  cout << endl;
  cout << "of the water cycle? ";
  cin >> wtol;

  rflog1 << endl;
  rflog1 << "What absolute tolerance do you want to use for checking equilibrium";
  rflog1 << endl;
  rflog1 << "of the water cycle? wtol = " << wtol << endl << endl;

  ctol = 1000.0;
  cout << endl;
  cout << "What absolute tolerance do you want to use for checking equilibrium";
  cout << endl;
  cout << "of the carbon cycle? ";
  cin >> ctol;

  rflog1 << endl;
  rflog1 << "What absolute tolerance do you want to use for checking equilibrium";
  rflog1 << endl;
  rflog1 << "of the carbon cycle?" << endl;
  rflog1 << "ctol = " << ctol << endl << endl;

  ntol = 1000.0;
  if (nfeed == 1)
  {
    rheqflag = 1;
    cout << endl;
    cout << "What absolute tolerance do you want to use for checking equilibrium";
    cout << endl;
    cout << "of the nitrogen cycle? ";
    cin >> ntol;

    rflog1 << endl;
    rflog1 << "What absolute tolerance do you want to use for checking equilibrium";
    rflog1 << endl;
    rflog1 << "of the nitrogen cycle?" << endl;
    rflog1 << "ntol = " << ntol << endl << endl;
  }

  if (equil == 0)
  {

    cout << endl << endl;
    cout << "What year do you want to start collecting output data? ";
    cin >> startyr;

    rflog1 << endl << endl;
    rflog1 << "What year do you want to start collecting output data? ";
    rflog1 << "startyr = " << startyr << endl;

    cout << endl << endl;
    cout << "What year do you want to stop collecting output data? ";
    cin >> endyr;

    rflog1 << endl << endl;
    rflog1 << "What year do you want to stop collecting output data? ";
    rflog1 << "endyr = " << endyr << endl;

    cout << "How often (x years) should data be collected after the initial year? ";
    cin >> diffyr;

    rflog1 << "How often (x years) should data be collected after the initial year? ";
    rflog1 << "diffyr = " << diffyr << endl;

  }

};

/* *************************************************************
************************************************************** */

// Modified due to the 2-box hydrology model
/* *************************************************************
************************************************************** */

//void TTEM::massbal(double y[NUMEQ], double prevy[NUMEQ])
void TTEM::massbal(double y[NVT][NUMEQ], double prevy[NVT][NUMEQ], const int& vt)
{


 /* if (y[vt][I_SM] != y[vt][I_AVLW] + soil.wiltpt[vt])
  {
    y[vt][I_SM] = y[vt][I_AVLW] + soil.wiltpt[vt];
  }
// added for hydrology model

   if (y[vt][I_SM1] != y[vt][I_AVLW1] + soil.wiltpt1[vt])
  {
    y[vt][I_SM1] = y[vt][I_AVLW1] + soil.wiltpt1[vt];
  }

   if (y[vt][I_SM2] != y[vt][I_AVLW2] + soil.wiltpt2[vt])
  {
    y[vt][I_SM2] = y[vt][I_AVLW2] + soil.wiltpt2[vt];
  }

    if (y[vt][I_SM3] != y[vt][I_AVLW3] + soil.wiltpt3[vt])
  {
    y[vt][I_SM3] = y[vt][I_AVLW3] + soil.wiltpt3[vt];
      }
	  */
	/*
  if (y[vt][I_PCTP] != 100.0 * y[vt][I_SM] / soil.totpor[vt])
  {
    y[vt][I_PCTP] = 100.0 * y[vt][I_SM] / soil.totpor[vt];
  }
  */
  /*
  if (y[vt][I_PCTP1] != 100.0 * y[vt][I_SM1] / soil.totpor1[vt])
  {
    y[vt][I_PCTP1] = 100.0 * y[vt][I_SM1] / soil.totpor1[vt];
  }

    if (y[vt][I_PCTP2] != 100.0 * y[vt][I_SM2] / soil.totpor2[vt])
  {
    y[vt][I_PCTP2] = 100.0 * y[vt][I_SM2] / soil.totpor2[vt];
  }

   if (y[vt][I_PCTP3] != 100.0 * y[vt][I_SM3] / soil.totpor3[vt])
  {
    y[vt][I_PCTP3] = 100.0 * y[vt][I_SM3] / soil.totpor3[vt];
      }
*/
/*
  if (y[vt][I_VSM] != y[vt][I_SM] / (soil.rootz[vt] * 1000.0))
  {
    y[vt][I_VSM] = y[vt][I_SM] / (soil.rootz[vt] * 1000.0);
    if (y[vt][I_VSM] <= 0.0) { y[vt][I_VSM] = 0.001; }
  }
  */
  if (y[vt][I_VSM1] != y[vt][I_SM1] / (soil.dpwbox1[vt] * 1000.0))
  {
    y[vt][I_VSM1] = y[vt][I_SM1] / (soil.dpwbox1[vt] * 1000.0);
    if (y[vt][I_VSM1] <= 0.0) { y[vt][I_VSM1] = 0.001; }
  }

//   if (y[vt][I_VSM2] != y[vt][I_SM2] / (soil.rootmx[vt] * 1000.0-soil.dpwbox1*1000.0))
   if (y[vt][I_VSM2] != y[vt][I_SM2] / (soil.dpwbox2[vt]*1000.0))
  {
    y[vt][I_VSM2] = y[vt][I_SM2] / (soil.dpwbox2[vt]*1000.0);
    if (y[vt][I_VSM2] <= 0.0) { y[vt][I_VSM2] = 0.001; }
  }

  if (y[vt][I_VSM3] != y[vt][I_SM3] / (soil.dpwbox3[vt]*1000.0))  // temporally set up as layer 2
  {
    y[vt][I_VSM3] = y[vt][I_SM3] / (soil.dpwbox3[vt]*1000.0);
    if (y[vt][I_VSM3] <= 0.0) { y[vt][I_VSM3] = 0.001; }
  }

// end of
  /*
  if ((y[vt][I_SNWPCK] - prevy[vt][I_SNWPCK]) != (y[vt][I_SNWFAL] - y[vt][I_SNWINF]))
  {
    y[vt][I_SNWINF] = y[vt][I_SNWFAL] - y[vt][I_SNWPCK] + prevy[vt][I_SNWPCK];
  }

  if ((y[vt][I_AVLW] - prevy[vt][I_AVLW]) != (y[vt][I_SNWINF] + y[vt][I_RAIN]
       - y[vt][I_RPERC] - y[vt][I_EET] - y[vt][I_SPERC]))
  {
    y[vt][I_SPERC] = y[vt][I_SNWINF] + y[vt][I_RAIN] - y[vt][I_RPERC] - y[vt][I_EET]
                 - y[vt][I_AVLW] + prevy[vt][I_AVLW];
  }

  if ((y[vt][I_RGRW] - prevy[vt][I_RGRW]) != (y[vt][I_RPERC] - y[vt][I_RRUN]))
  {
    y[vt][I_RRUN] = y[vt][I_RPERC] - y[vt][I_RGRW] + prevy[vt][I_RGRW];
  }

  if ((y[vt][I_SGRW] - prevy[vt][I_SGRW]) != (y[vt][I_SPERC] - y[vt][I_SRUN]))
  {
    y[vt][I_SRUN] = y[vt][I_SPERC] - y[vt][I_SGRW] + prevy[vt][I_SGRW];
  }

  if (y[vt][I_WYLD] != y[vt][I_RRUN] + y[vt][I_SRUN])
  {
    y[vt][I_WYLD] = y[vt][I_RRUN] + y[vt][I_SRUN];
  }*/
/************************* Carbon Cycle Balances **************************/
 /* if (y[vt][I_INNPP] < y[vt][I_NPP]) { y[vt][I_INNPP] = y[vt][I_NPP]; }

  if (y[vt][I_INGPP] < y[vt][I_GPP]) { y[vt][I_INGPP] = y[vt][I_GPP]; }

  if (y[vt][I_GPR] != y[vt][I_GPP] - y[vt][I_NPP])
  {
    y[vt][I_GPR] = y[vt][I_GPP] - y[vt][I_NPP];
  }

  if (y[vt][I_GPR] != y[vt][I_RVMNT] + y[vt][I_RVGRW])
  {
    y[vt][I_RVGRW] = y[vt][I_GPR] - y[vt][I_RVMNT];
  }

  if (ag.state == 0)
  {
    // Minimum VEGC added by DWK on 20000207
    //  if (y[I_VEGC] < 0.00001) { y[I_VEGC] = 0.00001; }
    if (y[vt][I_VEGC] - prevy[vt][I_VEGC] != y[vt][I_NPP] - y[vt][I_LTRC])
    {
      y[vt][I_LTRC] = y[vt][I_NPP] - y[vt][I_VEGC] + prevy[vt][I_VEGC];
    }
  }
  else
  {
    if (y[vt][I_AGNPPC] != y[vt][I_AGLTRC] + y[vt][I_AGFPRDC])
    {
      y[vt][I_AGLTRC] = y[vt][I_AGNPPC] - y[vt][I_AGFPRDC];
    }
  }

  // Minimum SOLC added by DWK on 20000207
  //  if (y[I_SOLC] < 0.00001) { y[I_SOLC] = 0.00001; }
  if (y[vt][I_SOLC] - prevy[vt][I_SOLC] != y[vt][I_LTRC] + y[vt][I_AGLTRC]
      + y[vt][I_SLASHC] - y[vt][I_SCNVRTC] - y[vt][I_RH])
  {
    y[vt][I_RH] = y[vt][I_LTRC] + y[vt][I_AGLTRC] + y[vt][I_SLASHC] - y[vt][I_SCNVRTC]
    - y[vt][I_SOLC] + prevy[vt][I_SOLC];
  }

  if (y[vt][I_NEP] != y[vt][I_NPP] + y[vt][I_AGNPPC] - y[vt][I_RH])
  {
    y[vt][I_NEP] = y[vt][I_NPP] + y[vt][I_AGNPPC] - y[vt][I_RH];
  }

  if (y[vt][I_CFLX] != y[vt][I_NEP] - y[vt][I_CNVRTC])
  {
    y[vt][I_CFLX] = y[vt][I_NEP] - y[vt][I_CNVRTC];
  }

  /*********************Nitrogen Cycle Balances**********************/
/*
  if (y[vt][I_VNUP] < 0.0 ) { y[vt][I_VNUP] = 0.0; }

  //if (y[I_VNUP] > y[I_INNUP]) { y[I_VNUP] = y[I_INNUP]; }
  if (y[vt][I_INNUP] < y[vt][I_VNUP]) { y[vt][I_INNUP] = y[vt][I_VNUP]; }

  if (y[vt][I_VSUP] < 0.0) { y[vt][I_VSUP] = 0.0; }

  if (y[vt][I_VSUP] > y[vt][I_VNUP]) { y[vt][I_VSUP] = y[vt][I_VNUP]; }

  if (y[vt][I_VLUP] != y[vt][I_VNUP] - y[vt][I_VSUP])
  {
    y[vt][I_VLUP] = y[vt][I_VNUP] - y[vt][I_VSUP];
  }

  if (ag.state == 0)
  {
    if (y[vt][I_STON] - prevy[vt][I_STON] != y[vt][I_VLUP] + y[vt][I_VNRSRB] - y[vt][I_VNMBL])
    {
      y[vt][I_VNRSRB] = y[vt][I_STON] - prevy[vt][I_STON] + y[vt][I_VNMBL] - y[vt][I_VLUP];
    }

    if (y[vt][I_STRN] - prevy[vt][I_STRN] != y[vt][I_VSUP] - y[vt][I_LTRN]
        - y[vt][I_VNRSRB] + y[vt][I_VNMBL])
    {
      y[vt][I_LTRN] = y[vt][I_VSUP] - y[vt][I_STRN] + prevy[vt][I_STRN]
                  - y[vt][I_VNRSRB] + y[vt][I_VNMBL];
    }
  }
  else
  {
    if (y[vt][I_AGNPPN] != y[vt][I_AGLTRN] + y[vt][I_AGFPRDN])
    {
      y[vt][I_AGLTRN] = y[vt][I_AGNPPN] - y[vt][I_AGFPRDN];
    }
  }

  if (y[vt][I_SOLN] - prevy[vt][I_SOLN] != y[vt][I_LTRN] + y[vt][I_AGLTRN] + y[vt][I_SLASHN]
      - y[vt][I_NMIN] - y[vt][I_SCNVRTN] - y[vt][I_NSRTNT])
  {
    y[vt][I_NMIN] = y[vt][I_LTRN] + y[vt][I_AGLTRN] + y[vt][I_SLASHN] - y[vt][I_SCNVRTN]
                - y[vt][I_NSRTNT] - y[vt][I_SOLN] + prevy[vt][I_SOLN];
  }


  if (y[vt][I_NLST] < 0.0) { y[vt][I_NLST] = 0.0; }

  if (ag.state == 0)
  {
    if (y[vt][I_NINP] < 0.0) { y[vt][I_NINP] = 0.0; }

    if (y[vt][I_AVLN] - prevy[vt][I_AVLN] != y[vt][I_NINP] - y[vt][I_NLST] + y[vt][I_NMIN]
        - y[vt][I_VNUP])
    {
      if (y[vt][I_NINP] + y[vt][I_NMIN] - y[vt][I_VNUP] - y[vt][I_AVLN] + prevy[vt][I_AVLN] > 0.0)
      {
        y[vt][I_NLST] =  y[vt][I_NINP] + y[vt][I_NMIN] - y[vt][I_VNUP] - y[vt][I_AVLN]
                     + prevy[vt][I_AVLN];
      }
      else
      {
        y[vt][I_NINP] = y[vt][I_NLST] - y[vt][I_NMIN] + y[vt][I_VNUP] + y[vt][I_AVLN]
                    - prevy[vt][I_AVLN];
      }
    }
  }

  if(ag.state == 1)
  {

// Following condition added to correct for negative AGFERTN
// D. Kicklighter 19990721

    if (y[vt][I_NRETNT] != y[vt][I_NVRTNT] + y[vt][I_NSRTNT])
    {
      y[vt][I_NRETNT] = y[vt][I_NVRTNT] + y[vt][I_NSRTNT];
    }

    if (y[vt][I_AGFRTN] < 0.0) { y[vt][I_AGFRTN] = 0.0; }

    if (y[vt][I_NINP] != y[vt][I_AGFRTN] + y[vt][I_NRETNT])
    {
      y[vt][I_NINP] = y[vt][I_AGFRTN] + y[vt][I_NRETNT];
    }

    if (y[vt][I_NINP] < 0.0) { y[vt][I_NINP] = 0.0; }

    if (y[vt][I_AVLN] - prevy[vt][I_AVLN] != y[vt][I_NINP] - y[vt][I_NLST] + y[vt][I_NMIN]
       - y[vt][I_AGNPPN])
    {
      y[vt][I_NLST] =  prevy[vt][I_AVLN] - y[vt][I_AVLN] + y[vt][I_NINP] + y[vt][I_NMIN]
                   - y[vt][I_AGNPPN];
    }
  }
*/
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

//void TTEM::monthxclm(const int& dcmnt, const double& tgppopt, const int& dm)
void TTEM::monthxclm(const int& dcmnt, const double& tgppopt, const int& dm, const int& vt)
{

double raq10;


/* temp: effect of temperature on primary productivity */

  if (atms.tair[vt][dm] <= veg.tmin[dcmnt] || atms.tair[vt][dm] >= veg.tmax[dcmnt])
  {
    temp[vt] = 0.0;
  }
  else
  {
    if (atms.tair[vt][dm] >= tgppopt && atms.tair[vt][dm] <= veg.toptmax[dcmnt])
    {
      temp[vt] = 1.0;
    }
    else
    {
      if (atms.tair[vt][dm] > veg.tmin[dcmnt] && atms.tair[vt][dm] < tgppopt)
      {
	temp[vt] = (atms.tair[vt][dm]-veg.tmin[dcmnt])*(atms.tair[vt][dm]-veg.tmax[dcmnt])
               /((atms.tair[vt][dm]-veg.tmin[dcmnt])*(atms.tair[vt][dm]-veg.tmax[dcmnt])
               - pow((atms.tair[vt][dm]-tgppopt),2.0));
      }
      else
      {
	temp[vt] = (atms.tair[vt][dm]-veg.tmin[dcmnt])*(atms.tair[vt][dm]-veg.tmax[dcmnt])
               /((atms.tair[vt][dm]-veg.tmin[dcmnt])*(atms.tair[vt][dm]-veg.tmax[dcmnt])
               - pow((atms.tair[vt][dm]-veg.toptmax[dcmnt]),2.0));
      }
    }
  }


/* respq10: effect of temperature on plant respiration  */

  raq10 = veg.raq10a0[dcmnt] + (veg.raq10a1[dcmnt]*atms.tair[vt][dm])
          + (veg.raq10a2[dcmnt]*pow(atms.tair[vt][dm],2.0))
          + (veg.raq10a3[dcmnt]*pow(atms.tair[vt][dm],3.0));
  respq10[vt] = pow(raq10,atms.tair[vt][dm]/10.0);


/* dq10: effect of temperature on decomposition */

// 19990821-previous year snowpack (soil.snowpack[dm]) changed
//          to previous month snowpack (y[NUMEEQ+2]) by Jim Long

// modified by Q Z. using soil temperature information


  if (stmflg==1) {
     dq10[vt] = pow(microbe.rhq10[dcmnt],atms.tsoil[vt][dm]/10.0);  // Replace air Temperature with Soil Temperature
    }
  else
   {

//  if (y[I_SNWPCK] > 0.0) { dq10 = 1.0; }
//  else
 // {
     dq10[vt] = pow(microbe.rhq10[dcmnt],atms.tair[vt][dm]/10.0);
   }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

//void TTEM::resetODEflux(double y[])
void TTEM::resetODEflux(double y[NVT][NUMEQ], const int& vt)
{
  int i;

  for (i = MAXSTATE; i < NUMEQ; i++) { y[vt][i] = 0.0; }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

//void TTEM::resetYrFlux(void)
void TTEM::resetYrFlux(const int& vt)
{

  int dm;

  // Annual carbon storage
  veg.yrcarbon[vt] = 0.0;
  soil.yrorgc[vt] = 0.0;

  // Annual nitrogen storage
  veg.yrnitrogen[vt] = 0.0;
  veg.yrstructn[vt] = 0.0;
  veg.yrc2n[vt] = 0.0;
  veg.yrstoren[vt] = 0.0;
  soil.yrorgn[vt] = 0.0;
  soil.yrc2n[vt] = 0.0;
  soil.yravln[vt] = 0.0;

  // Annual carbon & nitrogen storage in agricultural ecosystems
  ag.PROD1[vt].carbon = 0.0;
  ag.PROD1[vt].nitrogen = 0.0;

  // Annual water storage
  soil.yravlh2o[vt] = 0.0;
  soil.yrrgrndh2o[vt] = 0.0;
  soil.yrsnowpack[vt] = 0.0;
  soil.yrsgrndh2o[vt] = 0.0;
  soil.yrsmoist[vt] = 0.0;
  soil.yrpctp[vt] = 0.0;
  soil.yrvsm[vt] = 0.0;


  // added for 3-box hydm

  soil.yravlh2o1[vt] = 0.0;
  soil.yravlh2o2[vt] = 0.0;
  soil.yravlh2o3[vt] = 0.0;

  soil.yrsmoist1[vt] = 0.0;
  soil.yrsmoist2[vt] = 0.0;
  soil.yrsmoist3[vt] = 0.0;

  soil.yrpctp1[vt] = 0.0;
  soil.yrpctp2[vt] = 0.0;
  soil.yrpctp3[vt] = 0.0;

  soil.yrvsm1[vt] = 0.0;
  soil.yrvsm2[vt] = 0.0;
  soil.yrvsm3[vt] = 0.0;


  // Annual carbon fluxes
  veg.yringpp[vt] = 0.0;
  veg.yrgpp[vt] = 0.0;
  veg.yrinnpp[vt] = 0.0;
  veg.yrnpp[vt] = 0.0;
  veg.yrltrc[vt] = 0.0;
  microbe.yrrh[vt] = 0.0;
  yrnep[vt] = 0.0;

  // Annual nitrogen fluxes
  soil.yrnin[vt] = 0.0;
  veg.yrinnup[vt] = 0.0;
  veg.yrnup[vt] = 0.0;
  veg.yrsup[vt] = 0.0;
  veg.yrlup[vt] = 0.0;
  veg.yrnmobil[vt] = 0.0;
  veg.yrnrsorb[vt] = 0.0;
  veg.yrltrn[vt] = 0.0;
  microbe.yrnmin[vt] = 0.0;
  soil.yrnlost[vt] = 0.0;

  // Annual water fluxes
  atms.yrrain[vt] = 0.0;
  soil.yrrperc[vt] = 0.0;
  soil.yrrrun[vt] = 0.0;
  atms.yrsnowfall[vt] = 0.0;
  soil.yrsnowinf[vt] = 0.0;
  soil.yrsperc[vt] = 0.0;
  soil.yrsrun[vt] = 0.0;
  atms.yrpet[vt] = 0.0;
  atms.yreet[vt] = 0.0;
  soil.yrh2oyld[vt] = 0.0;

 // added for hydrology model

  hyd.yrhevap[vt] = 0.0;
  hyd.yrhtrans[vt] = 0.0;
  hyd.yrsevap[vt] = 0.0;
  hyd.yrsnowsub[vt] = 0.0;
  hyd.yrsubcan[vt] = 0.0;

  // for soil temperature
  atms.yrfrontd[vt] =0.0;
  atms.yrthawbegin[vt] =0.0;
  atms.yrthawend[vt] =0.0;
  atms.yrtsoil[vt] = 0.0;
  atms.yrdst5[vt] = 0.0;
  atms.yrdst10[vt] = 0.0;
  atms.yrdst20[vt] = 0.0;
  atms.yrdst50[vt] = 0.0;
  atms.yrdst100[vt] = 0.0;
  atms.yrdst200[vt] = 0.0;

  // Phenology
  veg.yrunleaf[vt] = 0.0;
  veg.yrleaf[vt] = 0.0;
  veg.yrfpc[vt] = 0.0;

  // Annual carbon and nitrogen fluxes from agricultural
  // conversion
  ag.yrconvrtC[vt] = 0.0;
  ag.yrvconvrtC[vt] = 0.0;
  ag.yrsconvrtC[vt] = 0.0;
  ag.yrconvrtN[vt] = 0.0;
  ag.yrvconvrtN[vt] = 0.0;
  ag.yrsconvrtN[vt] = 0.0;
  ag.yrnrent[vt] = 0.0;
  ag.yrnvrent[vt] = 0.0;
  ag.yrnsrent[vt] = 0.0;
  ag.yrslashC[vt] = 0.0;
  ag.yrslashN[vt] = 0.0;
  ag.formPROD10[vt].carbon  = 0.0;
  ag.formPROD10[vt].nitrogen  = 0.0;
  ag.formPROD100[vt].carbon = 0.0;
  ag.formPROD100[vt].nitrogen = 0.0;

  // Annual carbon and nitrogen fluxes from agriculture
  ag.yrnppC[vt] = 0.0;
  ag.yrnppN[vt] = 0.0;
  ag.formPROD1[vt].carbon = 0.0;
  ag.formPROD1[vt].nitrogen = 0.0;
  ag.yrltrc[vt] = 0.0;
  ag.yrltrn[vt] = 0.0;
  ag.yrfertn[vt] = 0.0;

  // Annual carbon and nitrogen fluxes from human products
  ag.formTOTPROD[vt].carbon = 0.0;
  ag.formTOTPROD[vt].nitrogen = 0.0;

  // Monthly carbon and nitrogen fluxes
  for (dm = 0; dm < CYCLE; dm++)
  {
    soil.ninput[vt][dm] = 0.0;
    soil.nlost[vt][dm] = 0.0;

    ag.slash[vt][dm].carbon = 0.0;
    ag.slash[vt][dm].nitrogen = 0.0;
    ag.vconvrtflx[vt][dm].carbon = 0.0;
    ag.sconvrtflx[vt][dm].carbon = 0.0;
    ag.convrtflx[vt][dm].carbon = 0.0;
    ag.vconvrtflx[vt][dm].nitrogen = 0.0;
    ag.sconvrtflx[vt][dm].nitrogen = 0.0;
    ag.convrtflx[vt][dm].nitrogen = 0.0;
    ag.nvretent[vt][dm] = 0.0;
    ag.nsretent[vt][dm] = 0.0;
    ag.nretent[vt][dm] = 0.0;
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

//void TTEM::setELMNTecd(const int& kdinflg, const int& dcmnt, const double& psiplusc)
void TTEM::setELMNTecd(const int& kdinflg, const int& dcmnt, const double& psiplusc, const int& vt)
{


// Initialize TEM parameters dependent upon a grid cell's soil texture

  if (psiplusc <= veg.cmaxcut[dcmnt])
  {
    veg.cmax[vt] = (veg.cmax1a[dcmnt] * psiplusc) + veg.cmax1b[dcmnt];
  }
  else
  {
      veg.cmax[vt] = (veg.cmax2a[dcmnt] * psiplusc) + veg.cmax2b[dcmnt];
  }

  if (kdinflg == 0)
  {
    microbe.kdc[vt] = (microbe.kda[dcmnt] / psiplusc) + microbe.kdb[dcmnt];
  }

  if (psiplusc <= veg.nmaxcut[dcmnt])
  {
    veg.nmax[vt] = (veg.nmax1a[dcmnt] * psiplusc) + veg.nmax1b[dcmnt];
  }
  else
  {
    veg.nmax[vt] = (veg.nmax2a[dcmnt] * psiplusc) + veg.nmax2b[dcmnt];
  }

  microbe.nup[vt] = (microbe.nupa[dcmnt] / psiplusc) + microbe.nupb[dcmnt];

// Set initial maximum relative leaf area

  veg.prvleafmx[dcmnt] = veg.initleafmx[dcmnt];  // Changed by DWK on 19991028


// Determine the "decay" parameter

  microbe.decay[vt] = 0.26299 + (1.14757*microbe.propftos[dcmnt])
                  - (0.42956*pow(microbe.propftos[dcmnt],2.0));

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */
// modified for hydrology model by QZ

//void TTEM::setELMNTevap(const int& stateflg, const int& dcmnt, double pet[CYCLE], double tair[CYCLE],double vap[CYCLE])
void TTEM::setELMNTevap(const int& stateflg, const int& dcmnt, double pet[NVT][CYCLE], double tair[NVT][CYCLE],double vap[NVT][CYCLE], const int& vt)
{

  int i;

// Determine initial values for atms.prvpetmx, atms.prveetmx,
//   and veg.topt

  if (stateflg == 1)
  {
    for (i = 0; i < CYCLE; i++)
    {
      if (pet[vt][i] >= atms.prvpetmx)
      {
	     veg.topt[vt] = tair[vt][i];
      }
      if (veg.aleaf[dcmnt] == 0.0 && veg.bleaf[dcmnt] == 0.0
         && veg.cleaf[dcmnt] == 1.0)
      {
	     if (tair[vt][i] > veg.topt[vt]) { veg.topt[vt] = tair[vt][i]; }
      }
      else
      {
	     if (veg.unnormleaf[vt][i] >= veg.prvleafmx[dcmnt])
        {
	       veg.topt[vt] = tair[vt][i];
	     }
      }
    }

    if (veg.topt[vt] > veg.toptmax[dcmnt]) { veg.topt[vt] = veg.toptmax[dcmnt]; }
    if (veg.topt[vt] < veg.toptmin[dcmnt]) { veg.topt[vt] = veg.toptmin[dcmnt]; }
  }

  else
  {
    for (i = 0; i < CYCLE; i++)
    {
    // modified for hydrology model by QZ
   // atms.pet[i] = atms.petjh(atms.nirr[i], tair[i], i);

   soil.avlh2o2[vt][i] =  y[vt][I_AVLW2];    //
 
   //veg.lai[vt][i] = 2.31; // for calibration version

 //  hyd.hvpd[i] = 10; // 20gm-3 will result in the close of stomatal need to reconsider

 // Determine monthly potential evapotranspiration
//    hyd.vap[vt][i] = atms.vap[vt][i]; // assign static vapor pressure to hyd.vap[i] instead of using transient data, need to figure out a way to deal with this conflict

  hyd.vaporpressure[vt][i] = hyd.vap[vt][i];

  hyd.vpd(hyd.vaporpressure[vt][i], atms.tair[vt][i],vt);
  hyd.intercept(atms.rain[vt][i], hyd.inter_coef[veg.cmnt],i,vt);
  hyd.canopy_water[vt][i] =hyd.prec_to_canopy[vt][i];
  hyd.drad(veg.lai[vt][i], hyd.EXT[veg.cmnt], atms.nirr[vt][i], i,vt);
//  hyd.hLWP[i] =  hyd.LWP(soil.avlh2o[i], soil.awcapmm);

  hyd.hLWP[vt][i] =  hyd.LWP(soil.avlh2o2[vt][i], soil.fldcap[vt], soil.pctsand, soil.pctclay, soil.pctsilt);

  hyd.hcch[vt][i] = hyd.CanopyCond(hyd.LWPmin[veg.cmnt], hyd.hLWP[vt][i],hyd.CCmax[veg.cmnt],hyd.LWPc[veg.cmnt],hyd.absoluteHD[vt],hyd.SLOPEcc[veg.cmnt]);

  hyd.hslope[vt][i] = hyd.slope(atms.tair[vt][i]);
  hyd.hpa[vt][i] = hyd.densityA(atms.tair[vt][i]);
  hyd.hlv[vt][i] = hyd.latentV(atms.tair[vt][i]);

  //rainfall is changed to
  atms.rain[vt][i] = hyd.rainthrough[vt][i];

  if (hyd.canopy_water[vt][i] >0)
   {
     hyd.htrans[vt][i] = hyd.penmanmonteith(hyd.hslope[vt][i], hyd.hdrad[vt][i], hyd.hpa[vt][i],hyd.vpdeficit[vt],
                    hyd.hcch[vt][i],hyd.hlv[vt][i]) * veg.lai[vt][i];
     hyd.hevap[vt][i] = hyd.evaporation(hyd.canopy_water[vt][i]);
     if ((hyd.hevap[vt][i] / hyd.htrans[vt][i]) < hyd.daze[i]) {
              hyd.htrans[vt][i] = hyd.htrans[vt][i] * (hyd.daze[i] - ceil( hyd.hevap[vt][i]/ hyd.htrans[vt][i]));
                            }
     else { hyd.htrans[vt][i] = 0.0; }
    }
   else // there is no canopy water
     {
    hyd.htrans[vt][i] = hyd.penmanmonteith(hyd.hslope[vt][i], hyd.hdrad[vt][i], hyd.hpa[vt][i],hyd.vpdeficit[vt],
                   hyd.hcch[vt][i],hyd.hlv[vt][i]) * veg.lai[vt][i];
    hyd.hevap[vt][i]=0.0;
    hyd.htrans[vt][i] = hyd.htrans[vt][i] * hyd.daze[i];
         }

  // Evaporation from soil surface



	// assign days since rain
 	hyd.dayssincerain = 10;

// first calculate potential evaporation, assuming the resistance
 // for vapor transport is equal to the resistance for sensible heat
 //	transport.  That is, no additional resistance for vapor transport to
 //	the soil surface. This represents evaporation from a wet surface with
 //	a specified aerodynamic resistance (= boundary layer resistance).
 //	The aerodynamic resistance is for now set as a constant, and is
 //	taken from observations over bare soil in tiger-bush in south-west
 //	Niger: aero_resis_layer1 = 107 s m-1 (Wallace and Holwill, 1997).


 // calculate potential_evap in mm/m2/day
   hyd.potevap_surface[vt][i] = hyd.evaporation_layer1(atms.tair[vt][i], hyd.hpa[vt][i],hyd.hlv[vt][i],hyd.vpdeficit[vt], hyd.surfacerad[vt][i])*hyd.daze[i];

  // consider only the precipitation flux reaching the soil
 // check for precipitation >= potential evaporation
	if (hyd.rainthrough[vt][i] >= hyd.potevap_surface[vt][i])
	{
		// reset days-since-rain parameter
		hyd.dayssincerain = 0.0;

		// soil evaporation proceeds at potential rate
		hyd.soil_evap[vt][i] = 0.6 * hyd.potevap_surface[vt][i];
	}
	else
	{
		// increment the days since rain
		hyd.dayssincerain += 1.0;

		// calculate the realized proportion of potential evaporation as a function of the days since rain
		hyd.ratio = 0.3/pow(hyd.dayssincerain,2.0);

		// calculate evaporation for dry days
		hyd.soil_evap[vt][i] = hyd.ratio * hyd.potevap_surface[vt][i];

		// for rain events that are smaller than required to reset dayssincerain counter, but larger than dry-day evaporation, all rain is evaporated.
		//In this case, do not advance the drying curve counter.
		//For rain events that are too small to trigger dsr reset, and which
		//are smaller than dry-day evap, there will be more evaporation than
		//rainfall.  In this case the drying curve counter is advanced.
		if (hyd.rainthrough[vt][i] > hyd.soil_evap[vt][i])
		{
			hyd.soil_evap[vt][i] = hyd.rainthrough[vt][i];
			hyd.dayssincerain -= 1.0;
		}
	}
// available water is changed to
   soil.avlh2o[vt][i] =  y[vt][I_AVLW] - hyd.soil_evap[vt][i];    // ???? need to think

  // end of evaporation from soil surface

    //sublimation from canopy
    //snow intercept
   if (soil.snowpack[vt][i] > 0.0) {
      hyd.snowpack_canopy[vt][i] = 0.01 * soil.snowpack[vt][i]; // 0.5mm/LAI is intercept coefficient,

      hyd.sub_from_canopy[vt][i] = hyd.sublimation_canopy(veg.lai[vt][i],hyd.hdrad[vt][i], atms.tair[vt][i]);

	  if (hyd.sub_from_canopy[vt][i] ==0.0)     {  hyd.sub_from_canopy[vt][i] = 0.001;} // to protect

      hyd.sublimationdays[vt][i] = hyd.snowpack_canopy[vt][i] / hyd.sub_from_canopy[vt][i];
      if (hyd.sublimationdays[vt][i] > hyd.daze[i])
	  { hyd.sub_from_canopy[vt][i] = hyd.snowpack_canopy[vt][i];}
	  else {hyd.sub_from_canopy[vt][i] = hyd.sub_from_canopy[vt][i] * hyd.sublimationdays[vt][i];}
    }

   // snowpack is changed to
    soil.snowpack[vt][i] = (soil.snowpack[vt][i] - hyd.sub_from_canopy[vt][i]);

  // sublimation from ground


  hyd.snowsub[vt][i] = hyd.snowsublimation_ground(soil.snowpack[vt][i], atms.tair[vt][i], hyd.surfacerad[vt][i]);

  atms.pet[vt][i] = hyd.htrans[vt][i] + hyd.hevap[vt][i] + hyd.soil_evap[vt][i]+ hyd.snowsub[vt][i] + hyd.sub_from_canopy[vt][i];
 // atms.pet[i] = hyd.htrans[i] + hyd.hevap[i] + hyd.soil_evap[i] + hyd.sub_from_canopy[i];

  // end of modification for hydrology model


      if (i == 0)
      {
	atms.prvpetmx = atms.prveetmx = atms.pet[vt][0];
	veg.topt[vt]  = tair[vt][0];
      }
      else
     {
	if (atms.pet[vt][i] > atms.prvpetmx)
        {
	  atms.prvpetmx = atms.prveetmx = atms.pet[vt][i];
	  veg.topt[vt] = tair[vt][i];
	}
      }
    }

    atms.yrpet[vt] = 1.0;
    atms.yreet[vt] = 1.0;
  }


// Determine vegetation C/N parameter as a function of vegetation type,
// annual PET, and annual EET (annual EET initially equal to yrpet)

  veg.updateC2N(dcmnt,atms.yreet[vt],atms.yrpet[vt],atms.co2[11],atms.initco2,vt);

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

//void TTEM::setELMNTflux(void)
void TTEM::setELMNTflux(const int& vt)
{

  // Initialize carbon, nitrogen and water fluxes including
  // ODE state variables (i.e., y[])
  int dm;

  for (dm = 0; dm < CYCLE; dm++)
  {
    // Initialize carbon fluxes in natural ecosystems to zero

    y[vt][I_INGPP] = veg.ingpp[vt][dm] = 0.0;
    y[vt][I_GPP] = veg.gpp[vt][dm] = 0.0;
    y[vt][I_INNPP] = veg.innpp[vt][dm] = 0.0;
    y[vt][I_NPP] = veg.npp[vt][dm] = 0.0;
    y[vt][I_GPR] = veg.gpr[vt][dm] = 0.0;
    y[vt][I_RVMNT] = veg.rm[vt][dm] = 0.0;
    y[vt][I_RVGRW] = veg.rg[vt][dm] = 0.0;
    y[vt][I_LTRC] = veg.ltrfal[vt][dm].carbon = 0.0;
    y[vt][I_RH] = microbe.rh[vt][dm] = 0.0;
    y[vt][I_NEP] = nep[vt][dm] = 0.0;

    // Initialize nitrogen fluxes in natural ecosystems to zero

    y[vt][I_NINP] = soil.ninput[vt][dm] = 0.0;
    y[vt][I_INNUP] = veg.inuptake[vt] = 0.0;
    y[vt][I_VNUP] = veg.nuptake[vt][dm] = 0.0;
    y[vt][I_VSUP] = veg.suptake[vt][dm] = 0.0;
    y[vt][I_VLUP] = veg.luptake[vt][dm] = 0.0;
    y[vt][I_VNMBL] = veg.nmobil[vt][dm] = 0.0;
    y[vt][I_VNRSRB] = veg.nresorb[vt][dm] = 0.0;
    y[vt][I_LTRN] = veg.ltrfal[vt][dm].nitrogen = 0.0;
    y[vt][I_MNUP] = microbe.nuptake[vt][dm] = 0.0;
    y[vt][I_NMIN] = microbe.netnmin[vt][dm] = 0.0;
    y[vt][I_NLST] = soil.nlost[vt][dm] = 0.0;

    // Initialize water fluxes to zero

    y[vt][I_RAIN] = 0.0;
    y[vt][I_RPERC] = soil.rperc[vt][dm] = 0.0;
    y[vt][I_RRUN] = soil.rrun[vt][dm] = 0.0;
    y[vt][I_SNWFAL] = 0.0;
    y[vt][I_SNWINF] = soil.snowinf[vt][dm] = 0.0;
    y[vt][I_SPERC] = soil.sperc[vt][dm] = 0.0;
    y[vt][I_SRUN] = soil.srun[vt][dm] = 0.0;
    y[vt][I_PET] = 0.0;
    y[vt][I_EET] = atms.eet[vt][dm] = 0.0;
    y[vt][I_WYLD] = soil.h2oyld[vt][dm] = 0.0;

    // added for 3-box hydrology
    y[vt][I_LYPERC] = hydm.lyperc[vt][dm] = 0.0;
    y[vt][I_LYPERC1] = hydm.lyperc1[vt][dm] = 0.0;
    y[vt][I_LYPERC2] = hydm.lyperc2[vt][dm] = 0.0;

    y[vt][I_RPERC1] = soil.rperc1[vt][dm] = 0.0;
    y[vt][I_RPERC2] = soil.rperc2[vt][dm] = 0.0;
    y[vt][I_RPERC3] = soil.rperc3[vt][dm] = 0.0;

    y[vt][I_SPERC1] = soil.sperc1[vt][dm] = 0.0;
    y[vt][I_SPERC2] = soil.sperc2[vt][dm] = 0.0;
    y[vt][I_SPERC3] = soil.sperc3[vt][dm] = 0.0;

  // end of ...

 // added for hydrology model by QZ
    y[vt][I_EVAP] = hyd.hevap[vt][dm] = 0.0;
    y[vt][I_TRANS] = hyd.htrans[vt][dm] = 0.0;
    y[vt][I_SEVAP] = hyd.soil_evap[vt][dm] = 0.0;
    y[vt][I_SNOWSUB] = hyd.snowsub[vt][dm] = 0.0;
    y[vt][I_SUBCAN] = hyd.sub_from_canopy[vt][dm] = 0.0;

  // for soil temperature
  y[vt][I_TSOIL] = atms.tsoil[vt][dm]=0.0;
  y[vt][I_DST5] = atms.dst5[vt][dm]=0.0;
  y[vt][I_DST10] = atms.dst10[vt][dm]=0.0;
  y[vt][I_DST20] = atms.dst20[vt][dm]=0.0;
  y[vt][I_DST50] = atms.dst50[vt][dm]=0.0;
  y[vt][I_DST100] = atms.dst100[vt][dm]=0.0;
  y[vt][I_DST200] = atms.dst200[vt][dm]=0.0;
  y[vt][I_FRONTD] = atms.frontd[vt][dm]=0.0;
  y[vt][I_THAWBE] = atms.thawbe[vt][dm]=0.0;
  y[vt][I_THAWEND] = atms.thawend[vt][dm]=0.0;
// end of ...

  // Initialize carbon and nitrogen fluxes during conversion
  //  to zero

    y[vt][I_CNVRTC] = ag.convrtflx[vt][dm].carbon = 0.0;
    y[vt][I_CNVRTN] = ag.convrtflx[vt][dm].nitrogen = 0.0;
    y[vt][I_SCNVRTC] = ag.sconvrtflx[vt][dm].carbon = 0.0;
    y[vt][I_SCNVRTN] = ag.sconvrtflx[vt][dm].nitrogen = 0.0;
    y[vt][I_NVRTNT] = ag.nvretent[vt][dm] = 0.0;
    y[vt][I_NSRTNT] = ag.nsretent[vt][dm] = 0.0;
    y[vt][I_NRETNT] = ag.nretent[vt][dm] = 0.0;
    y[vt][I_SLASHC] = ag.slash[vt][dm].carbon = 0.0;
    y[vt][I_SLASHN] = ag.slash[vt][dm].nitrogen = 0.0;
    y[vt][I_PRDF10C] = ag.formPROD10[vt].carbon = 0.0;
    y[vt][I_PRDF10N] = ag.formPROD10[vt].nitrogen = 0.0;
    y[vt][I_PRDF100C] = ag.formPROD100[vt].carbon = 0.0;
    y[vt][I_PRDF100N] = ag.formPROD100[vt].nitrogen = 0.0;

    // Initialize carbon and nitrogen in agricultural ecosystems
    //   to zero

    y[vt][I_AGNPPC] = ag.npp[vt][dm].carbon = 0.0;
    y[vt][I_AGNPPN] = ag.npp[vt][dm].nitrogen = 0.0;
    y[vt][I_AGFPRDC] = ag.formPROD1[vt].carbon = 0.0;
    y[vt][I_AGFPRDN] = ag.formPROD1[vt].nitrogen = 0.0;
    y[vt][I_AGLTRC] = ag.ltrfal[vt][dm].carbon = 0.0;
    y[vt][I_AGLTRN] = ag.ltrfal[vt][dm].nitrogen = 0.0;
    y[vt][I_AGFRTN] = ag.fertn[vt][dm] = 0.0;

    // Initialize carbon and nitrogen from human products to zero

    y[vt][I_TOTFPRDC] = ag.formTOTPROD[vt].carbon = 0.0;
    y[vt][I_TOTFPRDN] = ag.formTOTPROD[vt].nitrogen = 0.0;
    y[vt][I_AGPRDFC] = ag.PROD1decay[vt].carbon = 0.0;
    y[vt][I_AGPRDFN] = ag.PROD1decay[vt].nitrogen = 0.0;
    y[vt][I_PRD10FC] = ag.PROD10decay[vt].carbon = 0.0;
    y[vt][I_PRD10FN] = ag.PROD10decay[vt].nitrogen = 0.0;
    y[vt][I_PRD100FC] = ag.PROD100decay[vt].carbon = 0.0;
    y[vt][I_PRD100FN] = ag.PROD100decay[vt].nitrogen = 0.0;
    y[vt][I_TOTPRDFC] = ag.TOTPRODdecay[vt].carbon = 0.0;
    y[vt][I_TOTPRDFN] = ag.TOTPRODdecay[vt].nitrogen = 0.0;

    // Initialize integrated carbon fluxes to zero

    y[vt][I_TOTNPP] = ag.totnpp[vt][dm] = 0.0;
    y[vt][I_CFLX] = cflux[vt][dm] = 0.0;

  }
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

//void TTEM::setMonth(int& dm, double y[])
void TTEM::setMonth(const int& dm, double y[NVT][NUMEQ], const int& vt)
{

//	y[vt][I_AVLW2] = soil.freeze(y[vt][I_AVLW2],atms.tair[vt][dm],vt,dm);
  // Carbon pools
  veg.plant[vt][dm].carbon = y[vt][I_VEGC];
  soil.org[vt][dm].carbon = y[vt][I_SOLC];
  totalc[vt][dm] = veg.plant[vt][dm].carbon + soil.org[vt][dm].carbon;

  // Nitrogen pools
  veg.strctrl[vt][dm].nitrogen = y[vt][I_STRN];
  veg.labile[vt][dm].nitrogen = y[vt][I_STON];
  veg.plant[vt][dm].nitrogen = veg.strctrl[vt][dm].nitrogen + veg.labile[vt][dm].nitrogen;
  soil.org[vt][dm].nitrogen = y[vt][I_SOLN];
  soil.availn[vt][dm] = y[vt][I_AVLN];

  // Water pools
  soil.avlh2o[vt][dm] = y[vt][I_AVLW];
  soil.rgrndh2o[vt][dm] = y[vt][I_RGRW];
  soil.snowpack[vt][dm] = y[vt][I_SNWPCK];
  soil.sgrndh2o[vt][dm] = y[vt][I_SGRW];
  soil.moist[vt][dm] = y[vt][I_SM];
  soil.pctp[vt][dm] = y[vt][I_PCTP];
  soil.vsm[vt][dm] = y[vt][I_VSM];

 /* if (atms.tair[vt][dm] < 0.0 && atms.tair[vt][dm-1] >= 0.0)
  {

	  psstate[vt][dm] = y[vt][I_AVLW2] * 0.1;
  }
  pddstate[vt][dm] = 1.1;

  delta(dm,psstate,pddstate,vt);     */

  soil.avlh2o1[vt][dm] = y[vt][I_AVLW1];
  soil.avlh2o2[vt][dm] = y[vt][I_AVLW2];
  soil.avlh2o3[vt][dm] = y[vt][I_AVLW3];
  //printf("%f    %f\n",atms.tair[vt][dm],soil.avlh2o2[vt][dm]);
  soil.moist1[vt][dm] = y[vt][I_MOIST1];
  soil.moist2[vt][dm] = y[vt][I_MOIST2];
  soil.moist3[vt][dm] = y[vt][I_MOIST3];

  // adding frozen and thaw mechanism to modify the soil moisture *******
  // **** need to modify, we may think to change veg.thawpercent[dm] into other vairable

  if ( ((dm >=0) && (dm <=5))|| ((dm >=10) && (dm <=11)))
   veg.thawpercent[vt][dm] = 0.4276 + 0.0128 * atms.tsoil[vt][dm] + 0.0041 * powl(atms.tsoil[vt][dm],2.0);
  else veg.thawpercent[vt][dm] = 1.0;

  soil.pctp1[vt][dm] = y[vt][I_PCTP1] * veg.thawpercent[vt][dm];
  soil.pctp2[vt][dm] = y[vt][I_PCTP2]* veg.thawpercent[vt][dm];

  veg.thawpercent[vt][dm] = 0.7147 + 0.1063 * atms.dst50[vt][dm] + 0.0167 * powl(atms.dst50[vt][dm],2.0);

  soil.pctp3[vt][dm] = y[vt][I_PCTP3]* veg.thawpercent[vt][dm];

  soil.vsm1[vt][dm] = y[vt][I_VSM1];
  soil.vsm2[vt][dm] = y[vt][I_VSM2];
  soil.vsm3[vt][dm] = y[vt][I_VSM3];


  // Monthly carbon fluxes in natural ecosystems
  veg.ingpp[vt][dm] = y[vt][I_INGPP];
  veg.gpp[vt][dm] = y[vt][I_GPP];                     
  veg.innpp[vt][dm] = y[vt][I_INNPP];
  veg.npp[vt][dm] = y[vt][I_NPP];
  veg.gpr[vt][dm] = y[vt][I_GPR];
  veg.rm[vt][dm] = y[vt][I_RVMNT];
  veg.rg[vt][dm] = y[vt][I_RVGRW];
  veg.ltrfal[vt][dm].carbon = y[vt][I_LTRC];
  microbe.rh[vt][dm] = y[vt][I_RH];
  nep[vt][dm] = y[vt][I_NEP];

  // Monthly nitrogen fluxes in natural ecosystems
  soil.ninput[vt][dm] = y[vt][I_NINP];
  veg.inuptake[vt] = y[vt][I_INNUP];
  veg.nuptake[vt][dm] = y[vt][I_VNUP];
  veg.suptake[vt][dm] = y[vt][I_VSUP];
  veg.luptake[vt][dm] = y[vt][I_VLUP];
  veg.nmobil[vt][dm] = y[vt][I_VNMBL];
  veg.nresorb[vt][dm] = y[vt][I_VNRSRB];
  veg.ltrfal[vt][dm].nitrogen = y[vt][I_LTRN];
  microbe.nuptake[vt][dm] = y[vt][I_MNUP];
  microbe.netnmin[vt][dm] = y[vt][I_NMIN];
  soil.nlost[vt][dm] = y[vt][I_NLST];

  // Monthly water fluxes                                                      wsr
  atms.rain[vt][dm] = y[vt][I_RAIN];
 // soil.rperc[vt][dm] = y[vt][I_RPERC];
 // soil.rrun[vt][dm] = y[vt][I_RRUN];
  atms.snowfall[vt][dm] = y[vt][I_SNWFAL];
  soil.snowinf[vt][dm] = y[vt][I_SNWINF];
  //soil.sperc[vt][dm] = y[vt][I_SPERC];
 // soil.srun[vt][dm] = y[vt][I_SRUN];
  atms.pet[vt][dm] = y[vt][I_PET];
  atms.eet[vt][dm] = y[vt][I_EET];
//  soil.h2oyld[vt][dm] = y[vt][I_WYLD];

 // added for 3-box hydrology
 // hydm.lyperc[vt][dm] = y[vt][I_LYPERC];
  hydm.lyperc1[vt][dm] = y[vt][I_LYPERC1];
  hydm.lyperc2[vt][dm] = y[vt][I_LYPERC2];

  soil.rperc1[vt][dm] = y[vt][I_RPERC1];
  soil.rperc2[vt][dm] = y[vt][I_RPERC2];
  soil.rperc3[vt][dm] = y[vt][I_RPERC3];

  soil.sperc1[vt][dm] = y[vt][I_SPERC1];
  soil.sperc2[vt][dm] = y[vt][I_SPERC2];
  soil.sperc3[vt][dm] = y[vt][I_SPERC3];


 // added for hydrology model by QZ
   hyd.hevap[vt][dm] = y[vt][I_EVAP];
   hyd.htrans[vt][dm] = y[vt][I_TRANS];
   hyd.soil_evap[vt][dm] = y[vt][I_SEVAP];
   hyd.snowsub[vt][dm] = y[vt][I_SNOWSUB];
   hyd.sub_from_canopy[vt][dm] = y[vt][I_SUBCAN];

// added for soil thermal model

  atms.tsoil[vt][dm]=  y[vt][I_TSOIL];
  atms.dst5[vt][dm]=   y[vt][I_DST5];
  atms.dst10[vt][dm]=  y[vt][I_DST10];
  atms.dst20[vt][dm]=  y[vt][I_DST20];
  atms.dst50[vt][dm]=  y[vt][I_DST50];
  atms.dst100[vt][dm]= y[vt][I_DST100];
  atms.dst200[vt][dm]= y[vt][I_DST200];
  atms.frontd[vt][dm]= y[vt][I_FRONTD];
  atms.thawbe[vt][dm]= y[vt][I_THAWBE];
  atms.thawend[vt][dm]= y[vt][I_THAWEND];

// end of adding

  // Monthly phenology
  veg.unnormleaf[vt][dm] = y[vt][I_UNRMLF];
  veg.leaf[vt][dm] = y[vt][I_LEAF];

  // Monthly carbon and nitrogen fluxes associated with
  //  agricultural conversion
  ag.convrtflx[vt][dm].carbon = y[vt][I_CNVRTC];
  ag.convrtflx[vt][dm].nitrogen = y[vt][I_CNVRTN];
  ag.sconvrtflx[vt][dm].carbon = y[vt][I_SCNVRTC];
  ag.sconvrtflx[vt][dm].nitrogen = y[vt][I_SCNVRTN];
  ag.nvretent[vt][dm] = y[vt][I_NVRTNT];
  ag.nsretent[vt][dm] = y[vt][I_NSRTNT];
  ag.nretent[vt][dm] = y[vt][I_NRETNT];

  // Monthly carbon and nitrogen fluxes from agricultural
  //   ecosystems
  ag.npp[vt][dm].carbon = y[vt][I_AGNPPC];
  ag.npp[vt][dm].nitrogen = y[vt][I_AGNPPN];
  ag.ltrfal[vt][dm].carbon = y[vt][I_AGLTRC];
  ag.ltrfal[vt][dm].nitrogen = y[vt][I_AGLTRN];
  ag.fertn[vt][dm] = y[vt][I_AGFRTN];

  // Monthly integrated carbon fluxes
  ag.totnpp[vt][dm] = y[vt][I_TOTNPP];
  cflux[vt][dm] = y[vt][I_CFLX];

  // Update sum of annual carbon storage
  veg.yrcarbon[vt]  += y[vt][I_VEGC];
  soil.yrorgc[vt] += y[vt][I_SOLC];

  // Update sum of annual nitrogen storage
  veg.yrnitrogen[vt]  += y[vt][I_STRN] + y[vt][I_STON];
  veg.yrstructn[vt] += y[vt][I_STRN];
  soil.yrorgn[vt] += y[vt][I_SOLN];
  soil.yravln[vt]  += y[vt][I_AVLN];
  veg.yrstoren[vt] += y[vt][I_STON];

  // Update sum of annual water storage
  soil.yravlh2o[vt] += y[vt][I_AVLW];
  soil.yrrgrndh2o[vt] += y[vt][I_RGRW];
  soil.yrsnowpack[vt] += y[vt][I_SNWPCK];
  soil.yrsgrndh2o[vt] += y[vt][I_SGRW];
  soil.yrsmoist[vt] += y[vt][I_SM];
  soil.yrpctp[vt] += y[vt][I_PCTP];
  soil.yrvsm[vt] += y[vt][I_VSM];

// added for 3-box hydrology water storage
  soil.yravlh2o1[vt] += y[vt][I_AVLW1];
  soil.yravlh2o2[vt] += y[vt][I_AVLW2];
  soil.yravlh2o3[vt] += y[vt][I_AVLW3];

  soil.yrsmoist1[vt] += y[vt][I_SM1];
  soil.yrsmoist2[vt] += y[vt][I_SM2];
  soil.yrsmoist3[vt] += y[vt][I_SM3];

  soil.yrpctp1[vt] += y[vt][I_PCTP1];
  soil.yrpctp2[vt] += y[vt][I_PCTP2];
  soil.yrpctp3[vt] += y[vt][I_PCTP3];

  soil.yrvsm1[vt] += y[vt][I_VSM1];
  soil.yrvsm2[vt] += y[vt][I_VSM2];
  soil.yrvsm3[vt] += y[vt][I_VSM3];

  // Update sum of annual carbon fluxes in natural ecosystems
  veg.yringpp[vt] += y[vt][I_INGPP];
  veg.yrgpp[vt]   += y[vt][I_GPP];
  veg.yrinnpp[vt] += y[vt][I_INNPP];
  veg.yrnpp[vt]   += y[vt][I_NPP];
  veg.yrltrc[vt]  += y[vt][I_LTRC];
  microbe.yrrh[vt]    += y[vt][I_RH];
  yrnep[vt]   += y[vt][I_NEP];

 // Update sum of annual nitrogen fluxes in natural ecosystems
  soil.yrnin[vt]   += y[vt][I_NINP];
  veg.yrinnup[vt] += y[vt][I_INNUP];
  veg.yrnup[vt]   += y[vt][I_VNUP];
  veg.yrsup[vt]    += y[vt][I_VSUP];
  veg.yrlup[vt]    += y[vt][I_VLUP];
  veg.yrnmobil[vt] += y[vt][I_VNMBL];
  veg.yrnrsorb[vt] += y[vt][I_VNRSRB];
  veg.yrltrn[vt]  += y[vt][I_LTRN];
  microbe.yrnmin[vt]  += y[vt][I_NMIN];
  soil.yrnlost[vt] += y[vt][I_NLST];

   // Update sum of annual water fluxes
  atms.yrrain[vt] += y[vt][I_RAIN];
  soil.yrrperc[vt] += y[vt][I_RPERC];
  soil.yrrrun[vt] += y[vt][I_RRUN];
  atms.yrsnowfall[vt] += y[vt][I_SNWFAL];
  soil.yrsnowinf[vt] += y[vt][I_SNWINF];
  soil.yrsperc[vt] += y[vt][I_SPERC];
  soil.yrsrun[vt] += y[vt][I_SRUN];
  atms.yrpet[vt] += y[vt][I_PET];
  atms.yreet[vt] += y[vt][I_EET];
  soil.yrh2oyld[vt] += y[vt][I_WYLD];

  // added for 2- box hydrological model
  hydm.yrlyperc[vt] += y[vt][I_LYPERC];
  hydm.yrlyperc1[vt] += y[vt][I_LYPERC1];
  hydm.yrlyperc2[vt] += y[vt][I_LYPERC2];

/*  soil.rperc1[dm] = y[I_RPERC1];
  soil.rperc2[dm] = y[I_RPERC2];
  soil.rperc3[dm] = y[I_RPERC3];

  soil.sperc1[dm] = y[I_SPERC1];
  soil.sperc2[dm] = y[I_SPERC2];
  soil.sperc3[dm] = y[I_SPERC3];
 */

  // Update sum of annual phenology in natural ecosystems
  veg.yrunleaf[vt] += y[vt][I_UNRMLF];
  veg.yrleaf[vt] += y[vt][I_LEAF];
  veg.yrfpc[vt] += veg.fpc[vt][dm];

 // Update sum of annual carbon and nitrogen fluxes from
 //   agricultural conversion
  ag.yrconvrtC[vt] += y[vt][I_CNVRTC];
  ag.yrconvrtN[vt] += y[vt][I_CNVRTN];
  ag.yrsconvrtC[vt] += y[vt][I_SCNVRTC];
  ag.yrsconvrtN[vt] += y[vt][I_SCNVRTN];
  ag.yrslashC[vt] += y[vt][I_SLASHC];
  ag.yrslashN[vt] += y[vt][I_SLASHN];

 // Update sum of annual carbon and nitrogen fluxes from
 //   agricultural ecosystems
  ag.yrnppC[vt] += y[vt][I_AGNPPC];
  ag.yrnppN[vt] += y[vt][I_AGNPPN];
  ag.formPROD1[vt].carbon += y[vt][I_AGFPRDC];
  ag.formPROD1[vt].nitrogen += y[vt][I_AGFPRDN];
  ag.yrltrc[vt] += y[vt][I_AGLTRC];
  ag.yrltrn[vt] += y[vt][I_AGLTRN];
  ag.yrfertn[vt] += y[vt][I_AGFRTN];

   // Update sum of annual integrated carbon fluxes
  yrcflux[vt] += y[vt][I_CFLX];

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TTEM::setPrevState(double prevState[NVT][NUMEQ],double currentState[NVT][NUMEQ], const int& vt)
{
  for (int i = 0; i < NUMEQ; i++) { prevState[vt][i] = currentState[vt][i]; }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

//int TTEM::stepyr(const int& dyr, const int& itype, int& intflag, double& tol)
int TTEM::stepyr(const int& dyr, const int& itype, int& intflag, double& tol, const int& vt)
{

  int dm;
  int mintflag;
  double soiltemp1;

  if (dyr == 0) { microbe.kd[vt] = microbe.kdc[vt]; }
  else
  {
    if (ag.state == 0 && ag.prvstate == 0)
    {
      microbe.kd[vt] = microbe.yrkd(nfeed, veg.yrltrc[vt], veg.yrltrn[vt], veg.cmnt,vt);
      ag.kd[vt] = microbe.kd[vt];
      ag.natsoil = soil.org[vt][CYCLE-1].carbon;
    }
    else
    {
      if (soil.org[vt][CYCLE-1].carbon < ag.natsoil)
      {
        microbe.kd[vt] = ag.kd[vt] * soil.org[vt][CYCLE-1].carbon / ag.natsoil;
      }
      else { microbe.kd[vt] = ag.kd[vt]; }
    }
  }

  // Reset annual fluxes to zero

  resetYrFlux(vt);

  // Convert natural vegetation to agriculture

  if (ag.state == 1 && ag.prvstate == 0)
  {

    y[vt][I_VEGC] = ag.vrespar[veg.cmnt] * veg.plant[vt][CYCLE-1].carbon;
    y[vt][I_STRN] = ag.vrespar[veg.cmnt] * veg.strctrl[vt][CYCLE-1].nitrogen;
    y[vt][I_STON] = ag.vrespar[veg.cmnt] * veg.labile[vt][CYCLE-1].nitrogen;

    // Create PROD10 and PROD100 from conversion to agriculture

    ag.conversion(veg.cmnt, veg, soil,vt);
    mossflag = 1; // Q.Z. for fire disturbance
  }

  // Revert to natural vegetation after cropland abandonment

  if (ag.state == 0 && ag.prvstate == 1)
  {
    veg.plant[vt][CYCLE-1].carbon = y[vt][I_VEGC];
    veg.strctrl[vt][CYCLE-1].nitrogen = y[vt][I_STRN];
    veg.labile[vt][CYCLE-1].nitrogen = y[vt][I_STON];

    mossflag = 0; // Q.Z. for fire disturbance

  }

  for (dm = 0; dm < CYCLE; dm++)
  {
    if (initFlag == 1)
    {
      if (atms.tco2flag != 0) { atms.co2[dm] = atms.tco2[dyr][dm]; }
      if (atms.ttairflag != 0) { atms.tair[vt][dm] = atms.ttair[dyr][dm]; }
      if (atms.tprecflag != 0) { atms.prec[vt][dm] = atms.tprec[dyr][dm]; }

     //added for hydrology model by QZ
      if (hyd.vapflag != 0) { hyd.vap[vt][dm] = hyd.tvap[dyr][dm]; }

      if (ag.tlulcflag != 0)
      {
        if (ag.RAP0flag == 0)
        {
	  ag.potnpp[vt][dm] = ag.tpotnpp[itype][dyr][dm];
        }
        else { ag.potnpp[vt][dm] = 0.0; }
      }
    }

    // Get environmental conditions for month "dm"
      // calculate snowpack for the specific year

//    atms.pet[dm] = atms.petjh(atms.nirr[dm], atms.tair[dm], dm);
//    if (atms.pet[dm] <= 0.0) { atms.pet[dm] = 0.001; }

	 atms.precsplt(atms.prec[vt][dm], atms.tair[vt][dm], atms.rain[vt][dm], atms.snowfall[vt][dm]);
    switch (dm) {
      case 0:  atms.prevtair = atms.tair[vt][CYCLE-1];
	       atms.prev2tair = atms.tair[vt][CYCLE-2];
          break;
      case 1:  atms.prevtair = atms.tair[vt][0];
	       atms.prev2tair = atms.tair[vt][CYCLE-1];
	       break;
      default: atms.prevtair = atms.tair[vt][dm-1];
	       atms.prev2tair = atms.tair[vt][dm-2];
	       break;
      }
    soil.snowinf[vt][dm] = soil.snowmelt(elev, atms.tair[vt][dm], atms.prevtair, y[vt][I_SNWPCK]);
    soil.snowpack[vt][dm]= atms.snowfall[vt][dm] - soil.snowinf[vt][dm];
    if (soil.snowpack[vt][dm] < 0.0) { soil.snowpack[vt][dm] = 0.0; }

   if (stmflg == 1) {
     switch (dm) {
     case 0: sthermal.airt19= (atms.tair[vt][CYCLE-1]+atms.tair[vt][0])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.tair[vt][0];
             sthermal.airt39= (atms.tair[vt][0]+atms.tair[vt][1])/2.0;
             break;
     case 11:sthermal.airt19= (atms.tair[vt][10]+atms.tair[vt][11])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.tair[vt][11];
             sthermal.airt39= (atms.tair[vt][11]+atms.tair[vt][0])/2.0;
            break;
     default:sthermal.airt19= (atms.tair[vt][dm-1]+atms.tair[vt][dm])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.tair[vt][dm];
             sthermal.airt39= (atms.tair[vt][dm]+atms.tair[vt][dm+1])/2.0;
            break;
         }

    switch (dm) {
     case 0: sthermal.hsnow19=(soil.snowpack[vt][CYCLE-1]+soil.snowpack[vt][0])/2.0/100.0; // satisfy soil thermal model
             sthermal.hsnow29= soil.snowpack[vt][0]/100.0;
             sthermal.hsnow39= (soil.snowpack[vt][1]+soil.snowpack[vt][0])/2.0/100.0;
             break;
     case 11: sthermal.hsnow19=(soil.snowpack[vt][11]+soil.snowpack[vt][10])/2.0/100.0; // satisfy soil thermal model
              sthermal.hsnow29= soil.snowpack[vt][11]/100.0;
              sthermal.hsnow39= (soil.snowpack[vt][11]+soil.snowpack[vt][0])/2.0/100.0;
            break;
     default: sthermal.hsnow19=(soil.snowpack[vt][dm-1]+soil.snowpack[vt][dm])/2.0/100.0; // satisfy soil thermal model
              sthermal.hsnow29= soil.snowpack[vt][dm]/100.0;
              sthermal.hsnow39= (soil.snowpack[vt][dm]+soil.snowpack[vt][dm+1])/2.0/100.0;
            break;
         }
//   cout << airt19 << " " << airt29 << " " << airt39 << " " << hsnow19 << " " << hsnow29 << " " << hsnow39 <<endl;


	sthermal.soiltemp_(&sthermal.airt19, &sthermal.airt29, &sthermal.airt39, &sthermal.hsnow19, &sthermal.hsnow29,
    &sthermal.hsnow39,sthermal.weight9, &sthermal.smass9, &sthermal.is9, &sthermal.ism19, sthermal.x9, sthermal.dx9,
     sthermal.xfb9, sthermal.xfa9,	sthermal.water9, sthermal.t9, &sthermal.tsoil,&sthermal.frontd,&sthermal.thawbegin,
    &sthermal.thawend,sthermal.diffsoilt, veg.cmnt,mossflag);


   //      ++kswitch;

      atms.frontd[vt][dm]=sthermal.frontd;
      atms.thawbe[vt][dm]=sthermal.thawbegin;
      atms.thawend[vt][dm]=sthermal.thawend;

      atms.tsoil[vt][dm]=sthermal.tsoil;
      atms.dst5[vt][dm]=sthermal.diffsoilt[0];
      atms.dst10[vt][dm]=sthermal.diffsoilt[1];
      atms.dst20[vt][dm]=sthermal.diffsoilt[2];
      atms.dst50[vt][dm]=sthermal.diffsoilt[3];
      atms.dst100[vt][dm]=sthermal.diffsoilt[4];
      atms.dst200[vt][dm]=sthermal.diffsoilt[5];

 /*** end of calling soil thermal model ***/
  } // stmflg ===1

 // * calculate the thawing-proportion in early month and later fall

   if (atms.tsoil[vt][dm] >= 0.0) // present month is not frozen
   {
     if ((atms.tsoil[vt][dm-1] <= 0.0))// start to thaw, May issue
     {
       k_coef = abs(atms.tsoil[vt][dm] - atms.tsoil[vt][dm-1]) / atms.daze[dm];
       n_frozen = max(0, abs(atms.tsoil[vt][dm-1] / k_coef));
       veg.thawpercent[vt][dm] = (atms.daze[dm] - n_frozen) / atms.daze[dm];
       veg.thawpercent[vt][dm] = veg.thawpercent[vt][dm] / 2.0; // start from half thawing time to grow up
      }
     else
      {
        if ((atms.tsoil[vt][dm+1] <= 0.0))// start to frozen, October issue
         {
          if (dm ==9) { // for october > 0.0 late fall
                  k_coef = abs(atms.tsoil[vt][dm] - atms.tsoil[vt][dm+1]) / atms.daze[dm];
                  n_thaw = max(0, abs(atms.tsoil[vt][dm] / k_coef));
                  veg.thawpercent[vt][dm] = 0.1+(n_thaw / atms.daze[dm]);
                 }
                 else {
                    k_coef = abs(atms.tsoil[vt][dm] - atms.tsoil[vt][dm+1]) / atms.daze[dm];
                    n_thaw = max(0, abs(atms.tsoil[vt][dm] / k_coef));
                    veg.thawpercent[vt][dm] = n_thaw / atms.daze[dm];
                    }
           }
//        else veg.thawpercent[dm] = -atms.frontd[dm] / 1.5;
        else veg.thawpercent[vt][dm] = 1.0;
    }
   }

  else //if (atms.dst5[dm] < 0.0) // present month is frozen < 0.0
     {
      if ((atms.tsoil[vt][dm-1] >= 0.0))    // start to froze
      {
       if (dm ==9) { // for october < 0.0, only half month could growing
                k_coef = abs(atms.tsoil[vt][dm] - atms.tsoil[vt][dm-1]) / atms.daze[dm];
                n_thaw = max(0, abs(atms.tsoil[vt][dm-1] / k_coef));
                veg.thawpercent[vt][dm] = (n_thaw / atms.daze[dm]) / 2.0;
                }
        else {

               k_coef = abs(atms.tsoil[vt][dm] - atms.tsoil[vt][dm-1]) / atms.daze[dm];
               n_thaw = max(0, abs(atms.tsoil[vt][dm-1] / k_coef));
               veg.thawpercent[vt][dm] = n_thaw / atms.daze[dm];
            }
       }
      else
        {
        if (atms.tsoil[vt][dm+1] >=0.0) { //April
        k_coef = abs(atms.tsoil[vt][dm+1] - atms.tsoil[vt][dm]) / atms.daze[dm];
        n_thaw = max(0, abs(atms.tsoil[vt][dm+1] / k_coef));
//        veg.thawpercent[dm] = (atms.daze[dm] - n_thaw) / atms.daze[dm];
        veg.thawpercent[vt][dm] =  n_thaw / atms.daze[dm];
                                  }
        else  veg.thawpercent[vt][dm] = 0.0;
         }
     }
   // end of adding for thawing-frozen
//} // end of month




 getenviron(dm, vt);
	//printf("                                     %f\n",y[vt][I_INGPP]);
	//soiltemp1 = atms.dst20[vt][dm];
	//soiltemp[vt][dm] = atms.dst20[vt][dm];



	//calculating soil temperature for methane production (1cm)
	int z = 0.0;
	

				for (int i = 0;i < 200;i++)
				{
					z = i+1;

					soil.intersoilt[i] = soil.InterpolatST(atms.dst5[vt][dm],atms.dst10[vt][dm],atms.dst20[vt][dm],
						atms.dst50[vt][dm],atms.dst100[vt][dm],atms.dst200[vt][dm],z);//assign soil temperature at 10cm for the entire array

				}

			



   soil.cal_methane(vt,dm,soil.watertable,veg.npp[vt][dm]);


  soil.cal_diff(soil.vsm1[vt][dm],soil.vsm2[vt][dm],soil.dpwbox1[vt],soil.dpwbox2[vt],vt,dm,soil.watertable,atms.tair[vt][dm]);

//printf("   %f      %f     %f          %f  %f       %f         %f          %f      %f  %f\n",soil.watertable,veg.innpp[vt][dm],veg.ingpp[vt][dm],veg.npp[vt][dm],veg.gpp[vt][dm],veg.rm[vt][dm],y[vt][I_VEGC],y[vt][I_AVLN],veg.nuptake[vt][dm],dq10[vt]);
//printf("%9.3f,%9.3f,%9.3f,%9.3f,%9.3f\n",veg.gpp[vt][dm],veg.npp[vt][dm],veg.ingpp[vt][dm],veg.innpp[vt][dm],y[vt][I_VEGC]);

	deltafreeze1 = soil.freeze1(y[vt][I_AVLW1],atms.tair[vt][dm],atms.prevtair,vt,dm);

	deltafreeze2 = soil.freeze2(y[vt][I_AVLW2],atms.tair[vt][dm],atms.prevtair,vt,dm);

	deltafreeze3 = soil.freeze3(y[vt][I_AVLW3],atms.tair[vt][dm],atms.prevtair,vt,dm);
	
	//printf("                                    %f     \n",deltafreeze);
    mintflag = adapt(NUMEQ,y,tol,dm,vt);
	
	
	//printf("                                                                           %f                      \n",veg.ingpp[vt][dm]);

    if (mintflag == 1) { intflag = 1; }

    if (blackhol != 0) { qualcon[dyr][itype] = 10; }

    massbal(y,prevy,vt);

    setPrevState(prevy,y,vt);

    // Update carbon, nitrogen and water pools and fluxes from
    //  integrator results

    setMonth(dm, y,vt);
    resetODEflux(y,vt);

  } // end of month

  //  Update maximum EET, maximum PET, GPP optimum temperature (veg.topt),
  //  and maximum leaf cover (veg.prvleafmx) for the current year
  for (dm = 0; dm < CYCLE; dm++)
    
  {
    if (dm == 0)
    {
      atms.prveetmx = atms.eet[vt][0];
      atms.prvpetmx = atms.pet[vt][0];
      veg.topt[vt]   = atms.tair[vt][0];
      veg.prvleafmx[veg.cmnt] = veg.unnormleaf[vt][0];
    }
    else
    {
      if (atms.eet[vt][dm] > atms.prveetmx) { atms.prveetmx = atms.eet[vt][dm]; }
      if (atms.pet[vt][dm] > atms.prvpetmx) { atms.prvpetmx = atms.pet[vt][dm]; }
      if (veg.aleaf[veg.cmnt] == 0.0 && veg.bleaf[veg.cmnt] == 0.0
         && veg.cleaf[veg.cmnt] == 1.0)
      {
	if (atms.tair[vt][dm] > veg.topt[vt]) { veg.topt[vt] = atms.tair[vt][dm]; }
      }
      else
      {
	if (veg.unnormleaf[vt][dm] > veg.prvleafmx[veg.cmnt])
        {
	  veg.prvleafmx[veg.cmnt] = veg.unnormleaf[vt][dm];
	  veg.topt[vt] = atms.tair[vt][dm];
	}
      }
    }
  }

  // Update optimum temperature parameters for GPP

  if (veg.topt[vt] > veg.toptmax[veg.cmnt]) { veg.topt[vt] = veg.toptmax[veg.cmnt]; }
  if (veg.topt[vt] < veg.toptmin[veg.cmnt]) { veg.topt[vt] = veg.toptmin[veg.cmnt]; }


// Determine vegetation C/N parameter as a function of vegetation type,
// annual PET, annual EET, CO2 concentration

  veg.updateC2N(veg.cmnt, atms.yreet[vt], atms.yrpet[vt], atms.co2[11], atms.initco2,vt);

  soil.yravlh2o[vt] /= 12.0;
  soil.yrrgrndh2o[vt] /= 12.0;
  soil.yrsnowpack[vt] /= 12.0;
  soil.yrsgrndh2o[vt] /= 12.0;
  soil.yrsmoist[vt] /= 12.0;
  soil.yrpctp[vt] /= 12.0;
  soil.yrvsm[vt] /= 12.0;

 // added for 2-box hydrology

   soil.yravlh2o1[vt] /= 12.0;
   soil.yravlh2o2[vt] /= 12.0;
   soil.yravlh2o3[vt] /= 12.0;

  soil.yrsmoist1[vt] /= 12.0;
  soil.yrsmoist2[vt] /= 12.0;
  soil.yrsmoist3[vt] /= 12.0;

  soil.yrpctp1[vt] /= 12.0;
  soil.yrpctp2[vt] /= 12.0;
  soil.yrpctp3[vt] /= 12.0;

  soil.yrvsm1[vt] /= 12.0;
  soil.yrvsm2[vt] /= 12.0;
  soil.yrvsm3[vt] /= 12.0;

  atms.yrtsoil[vt] /= 12.0; // for soil temperature
  atms.yrfrontd[vt] /= 12.0; // for soil temperature
  atms.yrthawbegin[vt] /= 12.0; // for soil temperature
  atms.yrthawend[vt] /= 12.0; // for soil temperature

  veg.yrcarbon[vt]  /= 12.0;
  veg.yrnitrogen[vt]  /= 12.0;
  veg.yrstructn[vt] /= 12.0;

  if (veg.yrstructn[vt] != 0.0)
  {
    veg.yrc2n[vt]  = veg.yrcarbon[vt] / veg.yrstructn[vt];
  }
//  else { veg.yrc2n = MISSING; }

  soil.yrorgc[vt] /= 12.0;
  soil.yrorgn[vt] /= 12.0;

  if (soil.yrorgn[vt] != 0.0)
  {
    soil.yrc2n[vt] = soil.yrorgc[vt] / soil.yrorgn[vt];
  }
//  else { soil.yrc2n = MISSING; }

  soil.yravln[vt]  /= 12.0;
  veg.yrstoren[vt] /= 12.0;
  veg.yrunleaf[vt] /= 12.0;
  veg.yrleaf[vt] /= 12.0;
  veg.yrfpc[vt] /= 12.0;


  // y[2] changed to y[I_SOLC] by DWK on 20000130 and
  // y[3] changed to y[I_SOLN] by DWK on 20000130
  if (baseline == 1)
  {
    soil.yrnin[vt] = 0.0;
    soil.yrnlost[vt] = 0.0;
    if (y[vt][I_SOLC]/microbe.cnsoil[veg.cmnt] >= y[vt][I_SOLN])
    {
      soil.yrnin[vt] = (y[vt][I_SOLC]/microbe.cnsoil[veg.cmnt]) - y[vt][I_SOLN];
    }
    else
    {
      soil.yrnlost[vt] = y[vt][I_SOLN] - (y[vt][I_SOLC]/microbe.cnsoil[veg.cmnt]);
    }
    y[vt][I_SOLN] = y[vt][I_SOLC]/microbe.cnsoil[veg.cmnt];
  }

  if (endeq > 0) { ++endeq; }

  return endeq;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

//int TTEM::transient(const int& dyr, const int& itype, double& tol, const int& RTIME)
int TTEM::transient(const int& dyr, const int& itype, double& tol, const int& RTIME, const int& vt)
{


  endeq = 0;

  if (atms.tco2flag == 1) { totyr = atms.co2year[dyr]; }
  else if (atms.ttairflag == 1) { totyr = atms.tairyear[dyr]; }
  else if (atms.tprecflag == 1) { totyr = atms.precyear[dyr]; }
//added for hydrology model by QZ
  else if (hyd.vapflag == 1) { totyr = hyd.vapyear[dyr]; }

  if (atms.ttairflag != 0) { atms.mxtair = atms.mxttair[dyr]; }
  if (atms.tprecflag != 0) { atms.yrprec[vt] = atms.yrtprec[dyr]; }
// added for hydrology model by QZ
  if (hyd.vapflag != 0) { hyd.yrvap[vt] = hyd.yrtvap[dyr]; }

  if (ag.tlulcflag == 1)
  {
    ag.state = ag.tstate[dyr];
    ag.RAP = ag.tRAP[dyr];
  }

  transtepyr(dyr,itype, intflag, tol, RTIME,vt);

  // Update annual agricultural product pools and fluxes

  ag.updateyr(dyr,vt);

  if (totyr == startyr) { wrtyr = 0;}
  if (totyr > startyr) {++wrtyr; }

  return wrtyr;

};


// adding RTIME for soil thermal model

//int TTEM::transtepyr(const int& dyr, const int& itype, int& intflag, double& tol, const int& RTIME)
int TTEM::transtepyr(const int& dyr, const int& itype, int& intflag, double& tol, const int& RTIME, const int& vt)
{

  int dm;
  int mintflag;
  double soiltemp;
  if (dyr == 0) { microbe.kd[vt] = microbe.kdc[vt]; }
  else
  {
    if (ag.state == 0 && ag.prvstate == 0)
    {
      microbe.kd[vt] = microbe.yrkd(nfeed, veg.yrltrc[vt], veg.yrltrn[vt], veg.cmnt,vt);
      ag.kd[vt] = microbe.kd[vt];
      ag.natsoil = soil.org[vt][CYCLE-1].carbon;
    }
    else
    {
      if (soil.org[vt][CYCLE-1].carbon < ag.natsoil)
      {
        microbe.kd[vt] = ag.kd[vt] * soil.org[vt][CYCLE-1].carbon / ag.natsoil;
      }
      else { microbe.kd[vt] = ag.kd[vt]; }
    }
  }

  // Reset annual fluxes to zero

  resetYrFlux(vt);

  // Convert natural vegetation to agriculture

 // Q. Z. modify the thickness of moss layer, if no disturbance, we assume there is no change, mossflag = 0; 29/Nov/2000
  if (ag.state == 1 && ag.prvstate == 0)
  {

    y[vt][I_VEGC] = ag.vrespar[veg.cmnt] * veg.plant[vt][CYCLE-1].carbon;
    y[vt][I_STRN] = ag.vrespar[veg.cmnt] * veg.strctrl[vt][CYCLE-1].nitrogen;
    y[vt][I_STON] = ag.vrespar[veg.cmnt] * veg.labile[vt][CYCLE-1].nitrogen;

    // Create PROD10 and PROD100 from conversion to agriculture

    ag.conversion(veg.cmnt, veg, soil,vt);

    mossflag = 1;

  }

  // Revert to natural vegetation after cropland abandonment

  // Q. Z. modify the thickness of moss layer, if no disturbance, we assume there is no change, mossflag = 0; 29/Nov/2000

  if (ag.state == 0 && ag.prvstate == 1)
  {
    veg.plant[vt][CYCLE-1].carbon = y[vt][I_VEGC];
    veg.strctrl[vt][CYCLE-1].nitrogen = y[vt][I_STRN];
    veg.labile[vt][CYCLE-1].nitrogen = y[vt][I_STON];

    mossflag =0; // addition for disturbance

  }

  for (dm = 0; dm < CYCLE; dm++)
  {
    if (initFlag == 1)
    {
      if (atms.tco2flag != 0) { atms.co2[dm] = atms.tco2[dyr][dm]; }
      if (atms.ttairflag != 0) { atms.tair[vt][dm] = atms.ttair[dyr][dm]; }
      if (atms.tprecflag != 0) { atms.prec[vt][dm] = atms.tprec[dyr][dm]; }

     //added for hydrology model by QZ
      if (hyd.vapflag != 0) { hyd.vap[vt][dm] = hyd.tvap[dyr][dm]; }

      if (ag.tlulcflag != 0)
      {
        if (ag.RAP0flag == 0)
        {
	  ag.potnpp[vt][dm] = ag.tpotnpp[itype][dyr][dm];
        }
        else { ag.potnpp[vt][dm] = 0.0; }
      }
    }

    // Get environmental conditions for month "dm"
  // calculate snowpack for the specific year

//  for (dm = 0; dm < CYCLE; dm++) {
 //    atms.pet[dm] = atms.petjh(atms.nirr[dm], atms.tair[dm], dm);
 //   if (atms.pet[dm] <= 0.0) { atms.pet[dm] = 0.001; }

	 atms.precsplt(atms.prec[vt][dm], atms.tair[vt][dm], atms.rain[vt][dm], atms.snowfall[vt][dm]);
     switch (dm) {
      case 0:  atms.prevtair = atms.tair[vt][CYCLE-1];
	       atms.prev2tair = atms.tair[vt][CYCLE-2];
          break;
      case 1:  atms.prevtair = atms.tair[vt][0];
	       atms.prev2tair = atms.tair[vt][CYCLE-1];
	       break;
      default: atms.prevtair = atms.tair[vt][dm-1];
	       atms.prev2tair = atms.tair[vt][dm-2];
	       break;
      }
    soil.snowinf[vt][dm] = soil.snowmelt(elev, atms.tair[vt][dm], atms.prevtair, y[vt][I_SNWPCK]);
    soil.snowpack[vt][dm]= atms.snowfall[vt][dm] - soil.snowinf[vt][dm];
     if (soil.snowpack[vt][dm] < 0.0) { soil.snowpack[vt][dm] = 0.0; }


  if (stmflg==1) {
// call soil thermal subroutine
    switch (dm) {
     case 0: if (dyr==0) {
             sthermal.airt19= atms.ttair[dyr][0];; // satisfy soil thermal model
             sthermal.airt29= (atms.ttair[dyr][1]+atms.ttair[dyr][0])/2.0;
             sthermal.airt39= (atms.ttair[dyr][1]+atms.ttair[dyr][2])/2.0; }
             else {
             sthermal.airt19= (atms.ttair[dyr-1][11]+atms.ttair[dyr][0])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.ttair[dyr][0];
             sthermal.airt39= (atms.ttair[dyr][0]+atms.ttair[dyr][1])/2.0;
              }
             break;
     case 11:if ( dyr<RTIME-3) {
             sthermal.airt19= (atms.ttair[dyr][CYCLE-3]+atms.ttair[dyr][10])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.ttair[dyr][11];
             sthermal.airt39= (atms.ttair[dyr][11]+atms.ttair[dyr+1][0])/2.0; }
             else {
             sthermal.airt19= (atms.ttair[dyr][CYCLE-2]+atms.ttair[dyr][11])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.ttair[dyr][11];
             sthermal.airt39= (atms.ttair[dyr][11]+atms.ttair[dyr][0])/2.0; }
             break;
     default: sthermal.airt19= (atms.ttair[dyr][dm-1]+atms.ttair[dyr][dm])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.ttair[dyr][dm];
             sthermal.airt39= (atms.ttair[dyr][dm]+atms.ttair[dyr][dm+1])/2.0;
             break;
         }

     // using soil.snowpack, which is calculated from the WBM, to drive soil thermal model
  //  printf("snowpack %5.2f ", soil.snowpack[dm]);
  //  getch();

     switch (dm) {
     case 0: sthermal.hsnow19=(soil.snowpack[vt][CYCLE-1]+soil.snowpack[vt][0])/2.0/100.0; // satisfy soil thermal model
             sthermal.hsnow29= soil.snowpack[vt][0]/100.0;
             sthermal.hsnow39= (soil.snowpack[vt][1]+soil.snowpack[vt][0])/2.0/100.0;
             break;
     case 11: sthermal.hsnow19=(soil.snowpack[vt][11]+soil.snowpack[vt][10])/2.0/100.0; // satisfy soil thermal model
              sthermal.hsnow29= soil.snowpack[vt][11]/100.0;
              sthermal.hsnow39= (soil.snowpack[vt][11]+soil.snowpack[vt][0])/2.0/100.0;
            break;
     default: sthermal.hsnow19=(soil.snowpack[vt][dm-1]+soil.snowpack[vt][dm])/2.0/100.0; // satisfy soil thermal model
              sthermal.hsnow29= soil.snowpack[vt][dm]/100.0;
              sthermal.hsnow39= (soil.snowpack[vt][dm]+soil.snowpack[vt][dm+1])/2.0/100.0;
            break;
               }

      sthermal.soiltemp_(&sthermal.airt19, &sthermal.airt29, &sthermal.airt39, &sthermal.hsnow19, &sthermal.hsnow29, &sthermal.hsnow39,sthermal.weight9, &sthermal.smass9, &sthermal.is9, &sthermal.ism19, sthermal.x9, sthermal.dx9, sthermal.xfb9, sthermal.xfa9,
		sthermal.water9, sthermal.t9, &sthermal.tsoil,&sthermal.frontd,&sthermal.thawbegin,&sthermal.thawend,sthermal.diffsoilt,veg.cmnt, mossflag);

    //  ++ integer (kswitch);
      atms.frontd[vt][dm]=sthermal.frontd;
      atms.thawbe[vt][dm]=sthermal.thawbegin;
      atms.thawend[vt][dm]=sthermal.thawend;

      atms.tsoil[vt][dm]=sthermal.tsoil;
      atms.dst5[vt][dm]=sthermal.diffsoilt[0];
      atms.dst10[vt][dm]=sthermal.diffsoilt[1];
      atms.dst20[vt][dm]=sthermal.diffsoilt[2];
      atms.dst50[vt][dm]=sthermal.diffsoilt[3];
      atms.dst100[vt][dm]=sthermal.diffsoilt[4];
      atms.dst200[vt][dm]=sthermal.diffsoilt[5];
     } //stmflg ==1
   /*** end of calling soil thermal model ***/

   // * calculate the thawing-proportion in early month and later fall

   if (atms.tsoil[vt][dm] >= 0.0) // present month is not frozen
   {
     if ((atms.tsoil[vt][dm-1] <= 0.0))// start to thaw, May issue
     {
       k_coef = abs(atms.tsoil[vt][dm] - atms.tsoil[vt][dm-1]) / atms.daze[dm];
       n_frozen = max(0, abs(atms.tsoil[vt][dm-1] / k_coef));
       veg.thawpercent[vt][dm] = (atms.daze[dm] - n_frozen) / atms.daze[dm];
       veg.thawpercent[vt][dm] = veg.thawpercent[vt][dm] / 2.0; // start from half thawing time to grow up
      }
     else
      {
        if ((atms.tsoil[vt][dm+1] <= 0.0))// start to frozen, October issue
         {
          if (dm ==9) { // for october > 0.0 late fall
                  k_coef = abs(atms.tsoil[vt][dm] - atms.tsoil[vt][dm+1]) / atms.daze[dm];
                  n_thaw = max(0, abs(atms.tsoil[vt][dm] / k_coef));
                  veg.thawpercent[vt][dm] = 0.1+(n_thaw / atms.daze[dm]);
                 }
                 else {
                    k_coef = abs(atms.tsoil[vt][dm] - atms.tsoil[vt][dm+1]) / atms.daze[dm];
                    n_thaw = max(0, abs(atms.tsoil[vt][dm] / k_coef));
                    veg.thawpercent[vt][dm] = n_thaw / atms.daze[dm];
                    }
           }
//        else veg.thawpercent[dm] = -atms.frontd[dm] / 1.5;
        else veg.thawpercent[vt][dm] = 1.0;
    }
   }

  else //if (atms.dst5[dm] < 0.0) // present month is frozen < 0.0
     {
      if ((atms.tsoil[vt][dm-1] >= 0.0))    // start to froze
      {
       if (dm ==9) { // for october < 0.0, only half month could growing
                k_coef = abs(atms.tsoil[vt][dm] - atms.tsoil[vt][dm-1]) / atms.daze[dm];
                n_thaw = max(0, abs(atms.tsoil[vt][dm-1] / k_coef));
                veg.thawpercent[vt][dm] = (n_thaw / atms.daze[dm]) / 2.0;
                }
        else {

               k_coef = abs(atms.tsoil[vt][dm] - atms.tsoil[vt][dm-1]) / atms.daze[dm];
               n_thaw = max(0, abs(atms.tsoil[vt][dm-1] / k_coef));
               veg.thawpercent[vt][dm] = n_thaw / atms.daze[dm];
            }
       }
      else
        {
        if (atms.tsoil[vt][dm+1] >=0.0) { //April
        k_coef = abs(atms.tsoil[vt][dm+1] - atms.tsoil[vt][dm]) / atms.daze[dm];
        n_thaw = max(0, abs(atms.tsoil[vt][dm+1] / k_coef));
//        veg.thawpercent[dm] = (atms.daze[dm] - n_thaw) / atms.daze[dm];
        veg.thawpercent[vt][dm] =  n_thaw / atms.daze[dm];
                                  }
        else  veg.thawpercent[vt][dm] = 0.0;
         }
    }
   // end of adding for thawing-frozen
//} // end of monthly

  getenviron(dm,vt);
	
	//printf("%f                      \n",atms.dst20[vt][dm]);
	
	
	//soiltemptrans[vt][dm] = atms.dst20[vt][dm];




		//calculating soil temperature for methane production (1cm)
	int z = 0.0;
	

				for(int i = 0;i < 200;i++)
				{
					z = i+1;

					soil.intersoilt[i] = soil.InterpolatST(atms.dst5[vt][dm],atms.dst10[vt][dm],atms.dst20[vt][dm],
						atms.dst50[vt][dm],atms.dst100[vt][dm],atms.dst200[vt][dm],z);//assign soil temperature at 10cm for the entire array

				}

			



   soil.cal_methane(vt,dm,soil.watertable,veg.npp[vt][dm]);
   
   soil.cal_diff(soil.vsm1[vt][dm],soil.vsm2[vt][dm],soil.dpwbox1[vt],soil.dpwbox2[vt],vt,dm,soil.watertable,atms.tair[vt][dm]);

  //printf("%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f\n",soil.watertable,veg.gpp[vt][dm],veg.npp[vt][dm],veg.ingpp[vt][dm],veg.innpp[vt][dm],y[vt][I_VEGC]);




    deltafreeze1 = soil.freeze1(y[vt][I_AVLW1],atms.tair[vt][dm],atms.prevtair,vt,dm);

deltafreeze2 = soil.freeze2(y[vt][I_AVLW2],atms.tair[vt][dm],atms.prevtair,vt,dm);

	deltafreeze3 = soil.freeze3(y[vt][I_AVLW3],atms.tair[vt][dm],atms.prevtair,vt,dm);

    mintflag = adapt(NUMEQ,y,tol,dm,vt);
	

    if (mintflag == 1) { intflag = 1; }

    if (blackhol != 0) { qualcon[dyr][itype] = 10; }

    massbal(y,prevy,vt);

    setPrevState(prevy,y,vt);

    // Update carbon, nitrogen and water pools and fluxes from
    //  integrator results

    setMonth(dm, y,vt);
    resetODEflux(y,vt);
	
  
  }

  //  Update maximum EET, maximum PET, GPP optimum temperature (veg.topt),
  //  and maximum leaf cover (veg.prvleafmx) for the current year
  for (dm = 0; dm < CYCLE; dm++)
  {
    if (dm == 0)
    {
      atms.prveetmx = atms.eet[vt][0];
      atms.prvpetmx = atms.pet[vt][0];
      veg.topt[vt]   = atms.tair[vt][0];
      veg.prvleafmx[veg.cmnt] = veg.unnormleaf[vt][0];
    }
    else
    {
      if (atms.eet[vt][dm] > atms.prveetmx) { atms.prveetmx = atms.eet[vt][dm]; }
      if (atms.pet[vt][dm] > atms.prvpetmx) { atms.prvpetmx = atms.pet[vt][dm]; }
      if (veg.aleaf[veg.cmnt] == 0.0 && veg.bleaf[veg.cmnt] == 0.0
         && veg.cleaf[veg.cmnt] == 1.0)
      {
	if (atms.tair[vt][dm] > veg.topt[vt]) { veg.topt[vt] = atms.tair[vt][dm]; }
      }
      else
      {
	if (veg.unnormleaf[vt][dm] > veg.prvleafmx[veg.cmnt])
        {
	  veg.prvleafmx[veg.cmnt] = veg.unnormleaf[vt][dm];
	  veg.topt[vt] = atms.tair[vt][dm];
	}
      }
    }
  }

  // Update optimum temperature parameters for GPP

  if (veg.topt[vt] > veg.toptmax[veg.cmnt]) { veg.topt[vt] = veg.toptmax[veg.cmnt]; }
  if (veg.topt[vt] < veg.toptmin[veg.cmnt]) { veg.topt[vt] = veg.toptmin[veg.cmnt]; }


// Determine vegetation C/N parameter as a function of vegetation type,
// annual PET, annual EET, CO2 concentration

  veg.updateC2N(veg.cmnt, atms.yreet[vt], atms.yrpet[vt], atms.co2[11], atms.initco2,vt);

  soil.yravlh2o[vt] /= 12.0;
  soil.yrrgrndh2o[vt] /= 12.0;
  soil.yrsnowpack[vt] /= 12.0;
  soil.yrsgrndh2o[vt] /= 12.0;
  soil.yrsmoist[vt] /= 12.0;
  soil.yrpctp[vt] /= 12.0;
  soil.yrvsm[vt] /= 12.0;

 // added for 2-box hydrology

   soil.yravlh2o1[vt] /= 12.0;
   soil.yravlh2o2[vt] /= 12.0;
   soil.yravlh2o3[vt] /= 12.0;

  soil.yrsmoist1[vt] /= 12.0;
  soil.yrsmoist2[vt] /= 12.0;
  soil.yrsmoist3[vt] /= 12.0;

  soil.yrpctp1[vt] /= 12.0;
  soil.yrpctp2[vt] /= 12.0;
  soil.yrpctp3[vt] /= 12.0;

  soil.yrvsm1[vt] /= 12.0;
  soil.yrvsm2[vt] /= 12.0;
  soil.yrvsm3[vt] /= 12.0;

  atms.yrtsoil[vt] /= 12.0; // for soil temperature
  atms.yrfrontd[vt] /= 12.0; // for soil temperature
  atms.yrthawbegin[vt] /= 12.0; // for soil temperature
  atms.yrthawend[vt] /= 12.0; // for soil temperature

  veg.yrcarbon[vt]  /= 12.0;
  veg.yrnitrogen[vt]  /= 12.0;
  veg.yrstructn[vt] /= 12.0;

  if (veg.yrstructn[vt] != 0.0)
  {
    veg.yrc2n[vt]  = veg.yrcarbon[vt] / veg.yrstructn[vt];
  }
//  else { veg.yrc2n = MISSING; }

  soil.yrorgc[vt] /= 12.0;
  soil.yrorgn[vt] /= 12.0;

  if (soil.yrorgn[vt] != 0.0)
  {
    soil.yrc2n[vt] = soil.yrorgc[vt] / soil.yrorgn[vt];
  }
//  else { soil.yrc2n = MISSING; }

  soil.yravln[vt]  /= 12.0;
  veg.yrstoren[vt] /= 12.0;
  veg.yrunleaf[vt] /= 12.0;
  veg.yrleaf[vt] /= 12.0;
  veg.yrfpc[vt] /= 12.0;


  // y[2] changed to y[I_SOLC] by DWK on 20000130 and
  // y[3] changed to y[I_SOLN] by DWK on 20000130
  if (baseline == 1)
  {
    soil.yrnin[vt] = 0.0;
    soil.yrnlost[vt] = 0.0;
    if (y[vt][I_SOLC]/microbe.cnsoil[veg.cmnt] >= y[vt][I_SOLN])
    {
      soil.yrnin[vt] = (y[vt][I_SOLC]/microbe.cnsoil[veg.cmnt]) - y[vt][I_SOLN];
    }
    else
    {
      soil.yrnlost[vt] = y[vt][I_SOLN] - (y[vt][I_SOLC]/microbe.cnsoil[veg.cmnt]);
    }
    y[I_SOLN][vt] = y[vt][I_SOLC]/microbe.cnsoil[veg.cmnt];
  }

  if (endeq > 0) { ++endeq; }

  return endeq;

};

/* *************************************************************
************************************************************* */

