/* **************************************************************
*****************************************************************
TSOIL423.CPP - object describing general characteristics of soil
            - modified by DWK on 20000102
*****************************************************************
************************************************************** */

#if !defined(TSOIL423_H)
  #include "tsoil423.hpp"
#endif

/* **************************************************************
************************************************************** */

Tsoil4::Tsoil4(void)
{

  text  = -99;
  wsoil = -99;

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tsoil4::getecd(ofstream& rflog1)
{

  char ecd[80];

  cout << "Enter name of the soil (.ECD) data file with parameter values: ";
  cout << endl;
  fpara >> ecd;

  rflog1 << "Enter name of the soil (.ECD) data file with parameter values: ";
  rflog1 << ecd << endl;

  getecd(ecd);

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void Tsoil4::getecd (char ecd[80])
{

  char dummy[12];
  ifstream infile;

  long update;

  infile.open(ecd, ios::in);

  if (!infile)
  {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  infile >> dummy >> dummy >> dummy;
  infile >> dummy >> pctpora >> update;
  infile >> dummy >> pctporb >> update;
  infile >> dummy >> fldcapa >> update;
  infile >> dummy >> fldcapb >> update;
  infile >> dummy >> wiltpta >> update;
  infile >> dummy >> wiltptb >> update;

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil4::getrootz(ofstream& rflog1)
{

  char ecd[80];

  cout << "Enter name of the data file containing the rooting depths:";
  cout << endl;
  cout << "               (e.g., ROOTZVEG.ECD)" << endl;
  fpara >> ecd;

  rflog1 << "Enter name of the data file containing the rooting depths:";
  rflog1 << endl;
  rflog1 << "               (e.g., ROOTZVEG.ECD)" << endl;
  rflog1 << ecd << endl;

  getrootz(ecd);

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil4::getrootz(char ecd[80])
{

  const int NUMVAR = 7;
  char dummy[NUMVAR][10];
  ifstream infile;

  int i;
  int dcmnt;
  int  rootveg[MAXCMNT];
  long update[MAXCMNT];
  char vegname[MAXCMNT][31];

  infile.open(ecd, ios::in);

  if (!infile)
  {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 1; dcmnt < MAXCMNT; dcmnt++)
  {
    infile >> rootveg[dcmnt] >> vegname[dcmnt];
    infile >> rootza[dcmnt] >> rootzb[dcmnt] >> rootzc[dcmnt];
    infile >> minrootz[dcmnt] >> update[dcmnt];

  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil4::lake(double& tair,double& prec,double& rain,double& snowfall,
       		  double& pet, double& eet, int& dm, const int& vt)
{

  rgrndh2o[vt][dm] = 0.0;
  sperc[vt][dm] = 0.0;
  snowpack[vt][dm] = 0.0;
  sgrndh2o[vt][dm] = 0.0;
  moist[vt][dm] = 0.0;

  if (tair >= -1.0)
  {
   rain = prec;
    snowfall = 0.0;
  }
  else
  {
    rain = 0.0;
    snowfall = prec;
  }

  eet = pet;
  h2oyld[vt][dm] = prec - pet;

};

/* *************************************************************
************************************************************* */
/* *************************************************************
************************************************************* */
// added for hydrology model
void Tsoil4::percol_hydm(double& tair, double& rain, double& snowinf, double& trans1, double& trans2,double& trans3,double& avlh2o1, double& avlh2o2,double& avlh2o3, double& lyperc1, double& lyperc2,const int& dm, const int& vt) {

 double extra1, extra2,extra3, extra;
 double recharge;
 
 extra = 0.0;
 extra1= 0.0;
 extra2= 0.0;
 extra3= 0.0;

 rperc[vt][dm]= 0.0;
 sperc[vt][dm]= 0.0;

 rperc1[vt][dm]= 0.0;
 sperc1[vt][dm]= 0.0;

 rperc2[vt][dm]= 0.0;
 sperc2[vt][dm]= 0.0;

 rperc3[vt][dm]= 0.0;
 sperc3[vt][dm]= 0.0;

 


 recharge = rain + snowinf;

 if (recharge <= 0.0) { recharge = 0.001; }
 if ( (avlh2o1 + rain + snowinf -lyperc1 - trans1) > awcapmm1[vt]) {
	  extra1 = rain + snowinf + avlh2o1 - lyperc1 - trans1 - awcapmm1[vt];  // trans1 is an evaporation not transpiration

 if ( avlh2o2+lyperc1- trans2 - lyperc2 > awcapmm2[vt]) {
 extra2=avlh2o2+lyperc1 - trans2 -awcapmm2[vt] - lyperc2;
 }

 if ( avlh2o3+lyperc2- trans3  > awcapmm3[vt]) {
 extra3=avlh2o3+lyperc2 - trans3 -awcapmm3[vt];
 }

 
  
  sperc1[vt][dm]=snowinf*extra1/recharge;
  sperc2[vt][dm]=snowinf*extra2/recharge;
  sperc3[vt][dm]=snowinf*extra3/recharge;

  rperc1[vt][dm]= rain * extra1 / recharge;
  rperc2[vt][dm]= rain * extra2 / recharge;
  rperc3[vt][dm]= rain * extra3 / recharge;

  sperc[vt][dm] = snowinf * extra / recharge;                
  rperc[vt][dm] = rain * extra / recharge;
 }



	//printf("%f   %f      %f       %f             %f\n",tair,recharge,extra3,rperc3[vt][dm],sperc3[vt][dm]);

};

/* *************************************************************
************************************************************* */

/* *************************************************************
************************************************************* */

void Tsoil4::percol(double& rain, double& snowinf, double& eet,
                    double& avlh2o, const int& dm, const int& vt)
{

  double extra;
  double recharge;
  sperc[vt][dm] = 0.0;
  rperc[vt][dm] = 0.0;

  recharge = rain + snowinf;
  if (recharge <= 0.0) { recharge = 0.001; }
  if ((avlh2o + rain + snowinf - eet) > awcapmm[vt])
  {
    extra = rain + snowinf + avlh2o - awcapmm[vt] - eet;
    sperc[vt][dm] = snowinf * extra / recharge;
    rperc[vt][dm] = rain * extra / recharge;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

double Tsoil4::rrunoff(const double& rgrndh2o, const double& rperc)
{

  double rrunof;

  rrunof = 0.5 * (rgrndh2o + rperc);

  return rrunof;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void Tsoil4::showecd(void)
{

  cout << endl << "                   SOIL CHARACTERISTICS OF SITE";
  cout << endl << endl;
  printf("PSAND    = %5.2lf      PSILT = %5.2lf      PCLAY = %5.2lf\n",
         pctsand, pctsilt, pctclay);

  printf("POROSITY = %5.2lf   PCFLDCAP = %5.2lf   PCWILTPT = %5.2lf\n",
         pctpor, pcfldcap, pcwiltpt);

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

double Tsoil4::snowmelt(const double& elev, const double& tair,
                        const double& prevtair, double& snowpack)
{

  double snowflux = 0.0;
/*																						wsr
  if (tair >= 0.0)
  {
    if (elev <= 500.0) { snowflux = snowpack;}
    else
    {
      if (prevtair < 0.0) { snowflux = 0.5 * snowpack; }
      else { snowflux = snowpack; }
    }
  }
*/

if(tair >= 0 && prevtair < 0)
{
if(elev <= 500)  { snowflux = snowpack;   snowpack = 0.0;}

else { snowflux = snowpack * 0.5; }
}
  
  return snowflux;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

double Tsoil4::srunoff(const double& elev, const double& tair,
                       const double& prevtair, const double& prev2tair,
                       const double& sgrndh2o, const double& sperc)
{

  double srunof = 0.0;

  if (tair >= -1.0)
  {
    if (prevtair < -1.0) { srunof = 0.1 * (sgrndh2o + sperc); }
    else
    {
      if (prev2tair < -1)
      {
	if (elev <= 500.0) { srunof = 0.5 * (sgrndh2o + sperc); }
	else { srunof = 0.25 * (sgrndh2o + sperc); }
      }
      else { srunof = 0.5 * (sgrndh2o + sperc); }
    }
  }

  return srunof;

};

/* *************************************************************
************************************************************* */
void Tsoil4::hydm_xtext(int& ez, double& pctsilt, double& pctclay, const int& vt)
 {
 double rootmin;
//double rootmx;

 // totpor[vt] = fldcap[vt] = wiltpt[vt] = MISSING;
 // awcapmm[vt] =  MISSING;

  //rootmin= 99.9;
  //rootmx=-999.9;

  //totpor1[vt] = fldcap1[vt] = wiltpt1[vt] = -999.9;
  //totpor2[vt] = fldcap2[vt] = wiltpt2[vt] = -999.9;
 // totpor3[vt] = fldcap3[vt] = wiltpt3[vt] = -999.9;

 // awcapmm1[vt] =  -999.9;
  //awcapmm2[vt] =  -999.9;
 // awcapmm3[vt] =  -999.9;

  psiplusc = (pctsilt + pctclay) * 0.01;

     
  if (psiplusc < 0.01) { psiplusc = 0.01; }

// unit is meter
//  dpwbox1 =  (-4.721 * pow(psiplusc, 2.0)) + (4.106 * psiplusc) + 0.3003;

  rootz[vt] = (rootza[ez] * pow(psiplusc, 2.0)) + (rootzb[ez] * psiplusc) + rootzc[ez];

  if (rootz[vt] < minrootz[ez]) { rootz[vt] = minrootz[ez]; }
  if (rootz[vt]< rootmin) {rootmin= rootz[vt];}
  if (rootz[vt]> rootmx) { rootmx=rootz[vt];}

  rootz[vt] = 6.1;

//  dpwbox1 = rootmin;
  dpwbox1[vt] = rootz[vt] * 0.1;  // a thin layer (surface laler 0.1m for the first root zone, maybe the moss, Desborough, 1997
  dpwbox2[vt] = rootz[vt] * 0.2; // needed to reconsider it
  dpwbox3[vt] = rootz[vt] * 0.5; // mineral account for 20%

  pctpor[vt] = (pctpora * psiplusc) + pctporb;
  pcfldcap[vt] = (fldcapa * psiplusc) + fldcapb;
  pcwiltpt[vt] = (wiltpta * psiplusc) + wiltptb;
//  totpor  = rootz * pctpor * 10.0;
//  fldcap  = rootz * pcfldcap * 10.0;
//  wiltpt  = rootz * pcwiltpt * 10.0;

  //totpor1[vt]  = 2000;  // assumption for moss layer, we need to define the pore space, field capacity, wilting point
  //fldcap1[vt]  = 600; // assumption
 // wiltpt1[vt]  = 300;  // assumption

  totpor1[vt]  = rootz[vt] * 0.1 * 78 * 10.0 ;
  fldcap1[vt]  = rootz[vt] *  0.1* 72 * 10.0;
  wiltpt1[vt]  = rootz[vt] * 0.1 * 5 * 10.0;

  totpor2[vt]  = rootz[vt] * 0.2 * 78 * 10.0 ;
  fldcap2[vt]  = rootz[vt] * 0.2* 72* 10.0;
  wiltpt2[vt]  = rootz[vt] * 0.2* 6 * 10.0;

  totpor3[vt]  = rootz[vt] * 0.5 * 76 * 10.0 ;
  fldcap3[vt]  = rootz[vt] * 0.5* 72* 10.0;
  wiltpt3[vt]  = rootz[vt] * 0.5* (2) * 10.0;

  totpor[vt]  = totpor1[vt] + totpor2[vt] + totpor3[vt];
  fldcap[vt]  = fldcap1[vt] + fldcap2[vt] + fldcap3[vt];
  wiltpt[vt]  = wiltpt1[vt] + wiltpt2[vt] + wiltpt3[vt];
  // now unit is converted to mm
  /*
// unit is mm
  totpor1  = dpwbox1 * pctpor * 10.0;
  fldcap1  = dpwbox1 * pcfldcap * 10.0;
  wiltpt1  = dpwbox1 * pcwiltpt * 10.0;

  totpor2  = (rootmx-dpwbox1) * pctpor * 10.0;
  fldcap2  = (rootmx-dpwbox1) * pcfldcap * 10.0;
  wiltpt2  = (rootmx-dpwbox1) * pcwiltpt * 10.0;

  totpor3  = (rootmx-dpwbox2) * pctpor * 10.0;
  fldcap3  = (rootmx-dpwbox2) * pcfldcap * 10.0;
  wiltpt3  = (rootmx-dpwbox2) * pcwiltpt * 10.0;
 */
  awcapmm[vt] = fldcap[vt] - wiltpt[vt];

  awcapmm1[vt] = fldcap1[vt] - wiltpt1[vt]; // fldcap1 - wiltpt1;

  awcapmm2[vt] = fldcap2[vt] - wiltpt2[vt];
  awcapmm3[vt] = fldcap3[vt] - wiltpt3[vt];
};

/* *************************************************************
************************************************************** */

void Tsoil4::xtext(int& cmnt, double& pctsilt, double& pctclay, const int& vt)
{

  totpor[vt] = fldcap[vt] = wiltpt[vt] = MISSING;
  awcapmm[vt] =  MISSING;

  psiplusc = (pctsilt + pctclay) * 0.01;
  if (psiplusc < 0.01) { psiplusc = 0.01; }

  rootz[vt] = (rootza[cmnt] * pow(psiplusc, 2.0)) + (rootzb[cmnt] * psiplusc)
          + rootzc[cmnt];
  if (rootz[vt] < minrootz[cmnt]) { rootz[vt] = minrootz[cmnt]; }

  pctpor[vt] = (pctpora * psiplusc) + pctporb;
  pcfldcap[vt] = (fldcapa * psiplusc) + fldcapb;
  pcwiltpt[vt] = (wiltpta * psiplusc) + wiltptb;

  totpor[vt]  = rootz[vt] * pctpor[vt] * 10.0;
  fldcap[vt]  = rootz[vt] * pcfldcap[vt] * 10.0;
  wiltpt[vt]  = rootz[vt] * pcwiltpt[vt] * 10.0;

  awcapmm[vt] = fldcap[vt] - wiltpt[vt];

};

double Tsoil4::freeze1(double& pfstate1, double& t1, double& t2, const int& vt,const int& dm)

{int i;
double pftstate1;	
   
    crit(t1,t2,vt,dm);
	if (temp == 1)
	{
		
		pftstate1 = pfstate1 * 0.8;  //assume the freezing is 0.9 
		
	}
	else if (temp == 0)
	{
	
	
			pftstate1 = -pfstate1/0.3;  //assume the thawing process
		    
	
	}
	else 
	{
		pftstate1 = 1;
	}
	return pftstate1;
};




double Tsoil4::freeze2(double& pfstate2, double& t1,double& t2,const int& vt,const int& dm)

{int i;
double pftstate2;	
   
    crit(t1,t2,vt,dm);
	if (temp == 1)
	{
		
		pftstate2 = pfstate2 * 0.8;  //assume the freezing is 0.9 
		
	}
	else if (temp == 0)
	{
	
	
			pftstate2 = -pfstate2/0.7;  //assume the thawing process
		    
	
	}
	else 
	{
		pftstate2 = 1;
	}
	return pftstate2;
};



double Tsoil4::freeze3(double& pfstate3, double& t1,double& t2, const int& vt,const int& dm)

{int i;
double pftstate3;	
   
    crit(t1,t2,vt,dm);
	if (temp == 1)
	{
		
		pftstate3 = pfstate3 * 0.7;  //assume the freezing is 0.9 
		
	}
	else if (temp == 0)
	{
	
	
			pftstate3 = pfstate3/0.7;  //assume the thawing process
		    
	
	}
	else 
	{
		pftstate3 = 1;
	}
	return pftstate3;
};



void Tsoil4::crit(double& t1,double& t2,const int& vt,const int& dm)

{int i;




   

	if (t1 >= 0.0 && t2 < 0.0)
	{
		temp = 0;

	}
	else if (t1 < 0.0 && t2 >= 0.0)
	{temp = 1;}
	else 
	{
		temp = -1;
	}
};



void Tsoil4::wetlandTAWA(double& VSM1,double& VSM2,double& box1,double& box2, double& porosity1,double& porosity2, const int& vt, const int& dm)
{
	// cf. Granberg1999WRR Fig.2 and Eq 1-5
	//printf("      %f     \n",VSM2);
	// meanVSM: %, mean VSM of the soil profile from vegetation surface to depth zb
	// double phi = 0.9; // soil porosity
	double meanVSM;

	meanVSM = (VSM1 * box1 + VSM2 * 0.2) / (box1 + 0.2);
	double zb = 0.3; // m, maximum water table depth for wetlands
	double Vtot = (meanVSM)*(0.3*1000); // mm, total volume of water in the profile from vegetation surface to depth zb
	double theta_s_min = 0.25; // minimum vegetation surface volumetric water content
	double z_theta_s_min = 0.03; // 0.1; (smaller = shallower water table depth from soil surface, eg., 170mm is shallower than 300mm) // m, theta_s decrease linearly with the decline of water table until z_theta_s_min meter,
				// beyond which (> z_theta_s_min m) theta_s keep constant as theta_s_min (Granberg1999WRR Fig.2)
	double az; // gradient of linear decrease of theta_s
	double zwt1; // water table depth solution 1
	double zwt2; // water table depth solution 2
   	double zwt;  // water table depth
	double theta_s; // vegetation surface volumetric water content
		//phi = (porosity1 + porosity2) / 2;

	phi = 0.78;
	// determine water table depth
	az = (phi - theta_s_min) / z_theta_s_min;
	
        zwt = sqrt(  max(0.000001, 1.8 * (phi * zb - (Vtot/1000.0)) / az)  );       // ? (if z <= z_theta_s_min)  

	zwt2 = 1.5 * (phi * zb - Vtot/1000.0) / (phi - theta_s_min);                // ? (if z > z_theta_s_min)   }

        //zwt = min(zwt1,zwt2);
        
	if (zwt > 0.3) zwt = 0.3;
	//if (zwt < 0.0) zwt = 0.0; // zxd: without considering standing water
	watertable = zwt * 1000.0; // mm
	printf("%f\n",watertable);
	// determine volumetric water content above water table (used for methane calculation)
	theta_s = max(theta_s_min, (phi - az * zwt));
	
	int n = int(zwt * 100.0);
	for (int i = 0; i < n; i++)
	{
		unsatthetaWL[i] = min(phi, theta_s + (phi - theta_s) * (i/100.0/zwt) * (i/100.0/zwt)); // m3/m3
	}	

	//printf("   %f       %f       %f            %f     %f\n",watertable,VSM1,VSM2,snowinf[vt][dm],);
}


void Tsoil4::cal_methane(const int& vt, const int& dm, double& watertable,double& npp,double& lowb)
{
	double lowb1;
	long int i;
	double depthz;
	int satzone;
	int z;
	//double lowb = 1000;                                                            											//wsr: set lowb to 1m equaling to the root zone
	double maxfresh = 450;			//wsr: max npp
                  
	                totalch4 = 0;
lowb1 = lowb+1000;
                        
			satzone = int (ceil(watertable / 10.0));//determine the node where watertable starts
			depthz =0.0;

			// Calculate CH4 production rate
			for ( i = satzone; i < ceil((lowb1)/ 10.0); i++)
			{
			
			        depthz = i * 10.0;
				
				// Soil substrate multiplier
				//  -Calculate SOM distributions
				fcdist = EffectOMD(depthz, rootz[vt], lowb1);
				//  -SOM effect multiplier
				fsomt = EffectOM(maxfresh, npp);

				// Tsoil effect multiplier
				fmstt = EffectST(intersoilt[i]);

				// Soil pH effect multiplier
				fpht = EffectPH(6.8);

				// Calculate redox potential
				ehlt = RedoxP(watertable, depthz,lowb1,vt);
				//  -Redox potential effect multiplier
				frxt = EffectRX(ehlt);

				// Calculate CH4 production rate
				ch4pro[i] = MethanePR(fsomt,fcdist,fmstt,fpht,frxt);// uM h-1
				
				//if (ch4pro[i] < 0.0)
			//	{ ch4pro[i] = 0.0;}
				//integrate methane total production for one month

				totalch4 = totalch4 + ch4pro[i];
				
                          }
				
                               
                                
                                
                                
        totalch4 = totalch4 * pow(10.0,-3.0) / (3.6 * pow(10.0,3.0));
	totalch4 = totalch4 * pow(10.0,4.0); // u mol m-2 s-1;
	totalch4 = totalch4 * pow(10.0,-6.0); // mol m-2 s-1;
	totalch4 = totalch4 * 16; // g m-2 s-1;
	totalch4 = totalch4 * 8.64 * pow(10.0,4.0); // g m-2 day-1;
	totalch4 = totalch4 * pow(10.0,3.0); // mg m-2 day-1;                                          //wsr
			
}




// Effects of organic matter distribution in the soil to methanogenesis
double Tsoil4::EffectOMD(const double& depthz, const double& rootd, const double& lowb2)
{

	double fcdis;	// the index of organic matter distribution as function of the rooting depth and upper and lower boundary
	double Rd;		// root depth, local variable
	// for vegetated area, vegetation type (2-34)
	// hydm_xtext() produce rooting depth, as meter, need to transfer to mm = 1000* rootz
	// lowb -- lower boundary -- parameter say, 1000mm, 1m

	Rd = rootd * 1000.0;  // converse m to mm for rooting depth

	
		if ((depthz <= Rd) && (depthz > 0.0)) {
			fcdis = 1.0;
		}
		if ( (depthz > Rd) && (depthz < lowb2) )
		{
			fcdis = exp( - (depthz - Rd) / 10.0);
		}

	
	return fcdis;
}




// Effects of organic matter to methanogenesis
double Tsoil4::EffectOM(const double& maxfresh, const double& vegNPP)
{


	double fsom;
	double freshorg; // fresh organic matter from the vegetation

	if (vegNPP <= 0.0){
		freshorg = 0.01 * maxfresh;
	}
	else
		freshorg = maxfresh + vegNPP;

	fsom = freshorg / maxfresh; if (fsom > 2) fsom = 2; //zxd: g C m-2 d-1?

	return fsom;

}


	//Interpolating soil temperature to obtain every 1 cm depth soil temperatures with Lagrange method
// Numerical method of computation and FORTRAN lanuage, Xianyi Zheng, 1986
double Tsoil4::InterpolatST(const double& dst1, const double& dst2,const double& dst3,
	const double& dst4, const double& dst5, const double& dst6,  double x)
{
	double f, p;
	int i, j;
	double x0[6], y0[6];

	x0[0] = 5.0;
	x0[1] = 10.0;
	x0[2] = 20.0;
	x0[3] = 50.0;
	x0[4] = 100.0;
	x0[5] = 200.0;

	y0[0] = dst1;
	y0[1] = dst2;
	y0[2] = dst3;
	y0[3] = dst4;
	y0[4] = dst5;
	y0[5] = dst6;

	f = 0.0;
	
	float distance;
	float distance_temp;
	int index1;
	int index2;
	float k; // slope of line determined by first and second nodes

	distance = 10000;
	for (j =0; j<6; j++) // find first nearest node /////////////////////////////////////////////////////
	{
		distance_temp=fabs(x0[j]-x);
		if (distance_temp<distance)
		{
			index1=j;
			distance=distance_temp;
		}
	}

	distance = 10000;
	for (j =0; j<6; j++)  // find second nearest node ////////////////////////////////////////////////////
	{
		if (j!=index1)
		{
			distance_temp=fabs(x0[j]-x);
			if (distance_temp<distance)
			{
				index2=j;
				distance=distance_temp;
			}
		}
	}

	k = (y0[index2] - y0[index1]) / (x0[index2] - x0[index1]);
	f = y0[index1] + k * (x - x0[index1]);

	return f;
	
}



// Effects of soil temperature to methanogenesis
double Tsoil4::EffectST( double soilt)
{
	double fmst;// effect of soil temperature on CH4
	double tt;	// local variable, delta temperature
	// proref: reference temperature for CH4 production
	// pch4q10: Q10 coefficient for CH4 production on soil temperature. 
	double pch4q10 = 6.0;
	double proref = -2.0;


	tt = (soilt - proref) / 10.0;
	fmst = pow(pch4q10, tt);

	return fmst;

}


// Effects of soil PH to methanogenesis, need reading in soil pH, 
// similar to elevation effects of soil PH to methanogenesis
double Tsoil4::EffectPH(const double& soilph)
{
	const double phmin = 3.0;	//min pH = 5.5 in Q.Zhuang et al.(2004)
	const double phmax = 9.0;
	const double phopt = 7.5;

	double fph;
	double v1, v2, v3;

	v1 = soilph - phmin; // readin from external file, similar to elevation
	v2 = soilph - phmax;
	v3 = soilph - phopt;

	fph = v1 * v2 / (v1 * v2 - pow(v3,2.0));// previous coding is wrong
	// fph = v1 * v2 / (v1 * (v2 - pow(v3,2.0)));	
	return fph;

}



// Effects of soil redox potential to methanogenesis
double Tsoil4::EffectRX(const double& ehl)
{
	double frxp;

	if (ehl <= -200.0)	// Ehl is soil redox potential at specific depth at a given time
		frxp = 1.0;
	else if ((ehl > -200) && (ehl <= -100.0))
		frxp = -0.01 * ehl - 1;
	else 
		frxp = 0.0;

	return frxp;
}


// calculate the redox potential
double Tsoil4::RedoxP(const double& watertable, const double& depthz,const double& lowb2 ,const int& vt)
{
	double ehl;

	double wfps;			// 0 - 1.0  as percentage
	const double CR = 200.0;// daily variation of redox potential (mvday-1), consistent with Fig3 in Zhang2002GBC
	double AL;				// plant aerenchyma factor
	double FCA = 0.0013;	// area of the cross section of a typical fine root (cm2), based on McClaugherty et al.(1982) //zxd: Barber and Silberbush, 1984 in Zhang2002GBC ?
	double RLDL = 10; //zxd: cm cm-3, undisturbed boreal forest, see Table 3 in Bauhus1999CJFR // double RLDL = 0.001;	// fine root length density (cm root cm-3 soil)
	double PA;				// between 0 to 1.0, grasses and sedges are good transport =1.0; tree = 0.5; mosses =0.0

	wfps = (fldcap1[vt] + fldcap2[vt])/2;	// wfpsl -- the water filled pore space -- pctp

	
    PA = 0.5;
	

	AL = FCA * PA * RLDL;

	if ((depthz <= lowb2) && (depthz > watertable))
	{
		ehl = CR * (AL - 1.0);// in saturated zone
	}
	else
	{
		ehl = CR * (AL + 1.0 - wfps);// in unsaturated zone //zxd: ehl>0 make full surpress on ch4 production, see Fig3 in Zhang2002GBC
	}

	return ehl;
}


// Methane production rate occurred between the water table and the lower boundary
// the soil column is devided as 5 cm depth interval
double Tsoil4::MethanePR(double fsom, double fcdis, double fmst, double fph, double frx)
{
	double methanepr;// rate of methane production at depth z and time t
	double mgo = 0.22;



	methanepr = mgo  * fmst * fph * frx*fsom*fcdis; 

	return methanepr;
}



//function called to calculate the diffusivity

void Tsoil4::cal_diff(double& VSM1,double& VSM2,double& box1,double& box2,const int& vt, const int& dm, double& watertable, const double& tair,double& lowb)
{

int i;
int sat;
//int lowb= 1000;
int satzone;
double depthz;

double meanVSM = (VSM1 * box1 + VSM2 * 0.2) / (box1 + 0.2);

satzone = int (ceil(watertable / 10.0));//determine the node where watertable starts
			


for (i= 0; i < satzone; i++)
{

sat = 0;

fdfs[i] = Diffusivity(25,30,45,sat);

}


for (i=satzone; i < ceil(lowb/10); i++)

{


sat = 1;
fdfs[i] = Diffusivity(25,30,45,sat);

}

// to redistribute the ch4 in the soil profile

ch4flx = BCH4flxandCon(lowb,watertable);




// oxidation between water  table and soil surface

depthz = 0.0;

for ( i=0; i < satzone; i++)

{

//calculate the CH4 concentration multiplier

fmct = EffectMC2(ch4con_b[i],5.0);      //wsr: set kc = 5.0 for peatland


oxygenC(0.2,lowb);

foxyt = EffectOXY(oxygenc[i],100);



//calculate Tsoil effect multiplier

fstt = EffectST2(intersoilt[i], 1.1,5.5);


//calculate SM effect multiplier

fsmt = EffectSM(0.0,0.7,0.3,meanVSM);                         //wsr ; using mean VSM

//calculate Redox Potential effect multiplier

depthz = (i+1) * 10.0;

ehlt = RedoxP(watertable, depthz,lowb,vt);

frxt2 = EffectRX2(ehlt);

//calculate cultivation intensity multiplier

fctt = 1.0;

//calculate CH4 oxidation rate (hourly)

ch4oxirate[i] = MethaneOR(2.0,fmct,foxyt,fstt,fsmt,frxt,fctt);

}
ch4flxfinal = CCH4flxandCon(lowb, watertable,  tair);

}













double Tsoil4::Diffusivity(const double& sand, const double& silt, const double& clay, int satyorn)
{	// See Walter et al., 2001

	const double d1 = 0.2;					// diffusion coefficient of CH4 in the air and unsaturated soil (cm2 s-1)
	const double d2 = 0.2* pow(10.0,-4.0);	// diffusion coefficient of CH4 in the saturated soil and water (cm2 s-1)
	const double tort = 0.66;	// a tortuosity coefficient, suggesting that the apparent path is about two-thirds the length of the real average path of diffusion in the soil
	const double pvsand = 0.45; // relative volume of coarse pores in sandy soils
	const double pvsilt = 0.20; // silt soils
	const double pvclay = 0.14; // clay soils

	double diffu;	// rate of methane oxidation rate at depth z and time t
	double fcoarse;	//relative volume of the coarse pores in the soils(%)

	fcoarse = (sand * pvsand + silt * pvsilt + clay * pvclay) * 0.01;

	if (satyorn == 0)  diffu = d1 * tort * fcoarse; //diffusivity in unsaturated zone
	else diffu = d2* tort * fcoarse;				//diffusivity in saturated zone

	return diffu;

}




// for anoxic and oxic soils
double Tsoil4::BCH4flxandCon(const double& lowb, const double& watertable)

{

	const int MAXN = 200; //maximum nodes
	int nodes; // number of nodes, which is related to low boundary (lowb)
	const double x = 1.0; // depth step 10mm or 1cm
	// const double aoc = 0.076; // atmospheric ch4 concentration (uM)
	const double aoc = 0.076 * 0.001; // atmospheric ch4 concentration (u mol cm-3)
	// const double aoc = 0.0836 * 0.001; // atmospheric ch4 concentration (u mol cm-3)
	// const double aoc = 0.0684 * 0.001; // atmospheric ch4 concentration (u mol cm-3)

	double a[MAXN],b[MAXN],c[MAXN], d[MAXN], u[MAXN];
	double k[MAXN]; // conductivity
	double co[MAXN]; //concentration
	double z[MAXN]; //depth
	double df[MAXN];//diffusion
        
	int i;
	double ch4flx;
	
                 
	z[0] = 0.0;
	for (i = 0; i < nodes; i++)
	{
		df[i] = fdfs[i];
		z[i+1] = z[i] + x;
	}

	// k[0] = k_zero;
	k[0] = df[0] / 1.0;
	co[0] = aoc;

	for (i =1; i < nodes+1; i++)
	{
		if (i < (watertable /10.0)) u[i]=0.0;
		else u[i] =   ch4ratesat[i] * pow(10.0,-3.0) / (3.6 * pow(10.0,3.0));  // u mol cm-3 s-1

		if (i < nodes+1) k[i] = df[i] / (z[i+1] - z[i]);
		else k[i] =0.0;
		a[i+1] = -k[i];
		b[i] = k[i-1] + k[i];
		c[i] = -k[i];
		d[i] = u[i];
	}

	d[1] = d[1] + k[0] * co[0];
	for (i=1; i < nodes; i++)
	{
		c[i] = c[i] /b[i];
		d[i] = d[i] / b[i];
		b[i+1] = b[i+1] -a[i+1]*c[i];
		d[i+1] = d[i+1] -a[i+1]*d[i];
	}

	co[nodes]= d[nodes]/b[nodes];

	for (i = nodes; i > 0; i--)
	{
		co[i] = d[i] - c[i] * co[i+1];
		ch4con_b[i] = co[i] / 0.001;  // converse to uM
	}

	// surface flux of methane, positive indicates flux to soil

	ch4flx = - df[1] * (co[1] - co[0]) / 1.0;  // u mol cm-2 s-1

	ch4flx = ch4flx * pow(10.0,4.0); // u mol m-2 s-1;
	ch4flx = ch4flx * pow(10.0,-6.0); // mol m-2 s-1;
	ch4flx = ch4flx * 16; // g m-2 s-1;
	ch4flx = ch4flx * 8.64 * pow(10.0,4.0); // g m-2 day-1;
	ch4flx = ch4flx * pow(10.0,3.0); // mg m-2 day-1;

	return ch4flx;
}



// Effects of methane concentration to methanotrophic reaction0
double Tsoil4::EffectMC2(const double& ch4con, const double& kc)
{
	double fmc;

	fmc = ch4con / (kc + ch4con); 
	// Kc is Michaelis-Menten coefficient. (Walter et al., 2000; Bender and Conrad, 1992)
	// ch4con is soil methane concentration at depth z and time t

	return fmc;
}


// calculate the oxygen concentration
void Tsoil4::oxygenC(const double& afp, const double& lowb )	//unit is extremely chaotic
{		
	int nodes; // number of nodes, which is related to low boundary (lowb)
	const int MAXN = 200; //maximum nodes
	const double x = 0.01; // depth step 10mm or 1cm
	//const double x = 0.001;
	const double aoc = 280.0; // atmospheric oxygen concentration (280 gm-3)
	const double k_zero = 0.01; // the atmospheric boundary layer conductance (k(0))
	const double d0 = 0.0000177; // binary diffusion coefficients for O2 (m2s-1) (Campbell, 1977)
	const double cb = 0.90; // constant b, Sallam et al., 1984
	const double cm = 2.36; // constant m, Sallam et al., 1984
	const double azs = 0.0005; // soil surface O2 uptake rate
	double a[MAXN],b[MAXN],c[MAXN], d[MAXN], u[MAXN];
	double k[MAXN]; // conductivity
	double co[MAXN]; //concentration
	double z[MAXN]; //depth
	double df[MAXN];//diffusion
	int i;

	nodes = int (ceil(lowb / 10.0));

	z[0] = 0.0;

	for (i = 0; i < nodes+1; i++)
	{
		z[i+1] = z[i] + x;
	}

	k[0] = k_zero;
	co[0] = aoc;

	for (i =1; i < nodes+1; i++)
	{
		df[i] = d0 * cb * pow(afp,cm);
		u[i] = azs * exp(-z[i] / 0.3) * (z[i+1] -z[i-1])/2.0;

		if (i < nodes+1) k[i] = df[i] / (z[i+1] - z[i]);
		else k[i] =0.0;

		a[i+1] = -k[i];
		b[i] = k[i-1] + k[i];
		c[i] = -k[i];
		d[i] = u[i];
	}

	d[1] = d[1] + k[0] * co[0];
	for (i=1; i < nodes; i++)
	{
		c[i] = c[i] /b[i];
		d[i] = d[i] / b[i];
		b[i+1] = b[i+1] -a[i+1]*c[i];
		d[i+1] = d[i+1] -a[i+1]*d[i];
	}

	co[nodes]= d[nodes]/b[nodes];

	for (i = nodes; i > 0; i--)
	{
		co[i] = d[i] - c[i] * co[i+1];
		oxygenc[i] = co[i];
	}
}



// Effects of soil temperature to methanogenesis
double Tsoil4::EffectST2(const double& soilt, const double& och4q10, const double& oxiref)
{
	double fst;

	fst = pow (och4q10, (soilt - oxiref) / 10.0);
	//och4q10  is a coefficient with a constant, Walter et al., set it as 2.0;
	//oxiref is the reference soil temperature (oC) (e.g. 23oC).

	if (soilt < 0.0) fst =0.0;

	return fst;

}

// Effects of oxygen concentration to oxidation
double  Tsoil4::EffectOXY(const double& oxyc, const double& ko)
{
	double foxy;
	double kk;

	kk = ko * 32 / pow(10.0,3.0);
	foxy = oxyc / (kk + oxyc);
	// Ko is Michaelis-Menten coefficient, typical value ranged from 37 to 200 uM as concentration
	// oxyc denotes the oxygen concentration in the soil

	return foxy;
}








// Effects of soil moisture to methanogenesis
// Similar to TEM for the heterotrophic respiration Tian et al.,
double Tsoil4::EffectSM(const double& mvmin, const double& mvmax, const double& mvopt, const double& soilm)
{

	double fsm;
	double temp1;//temporary variable
	double temp2;
	double temp3;

	temp1 = soilm - mvmin;
	temp2 = soilm - mvmax;
	temp3 = soilm - mvopt;

	fsm = (temp1 * temp2) / (temp1 * temp2 - pow(temp3,2.0));

	return fsm;

}



// Effects of soil redox potential to methanotrophy
double Tsoil4::EffectRX2(const double& ehl)
{
	double frx;

	if ((ehl <= -100.0) && (ehl >= -200))  frx = 0.0075 * ehl + 1.5;
	if ((ehl > -100.0) && (ehl <= 200.0))  frx = 0.00083 * ehl + 0.833;
	if ((ehl > 200) && (ehl <= 600)) frx = 1.0;

	return frx;
}






// Methane oxidation occurred between the soil surface and the water table
// the soil column is devided as 1 cm depth interval
double Tsoil4::MethaneOR(const double& omax, double& fmc, double& foxy, double& fst, double& fsm, double& frx, double& fct )
{
	double methaneor; // methane oxidation rate at depth z and time t

	methaneor = omax * fmc * foxy  * fsm  * fct;

	return methaneor;
}











// for anoxic and oxic soils, eventual flux and concentration profile
double Tsoil4::CCH4flxandCon(const double& lowb, const double& watertable, const double& soilt)

{

	const int MAXN = 200; //maximum nodes
	int nodes; // number of nodes, which is related to low boundary (lowb)
	const double x = 1.0; // depth step 10mm or 1cm
	// const double aoc = 0.076; // atmospheric ch4 concentration (uM)
	const double aoc = 0.076 * 0.001; // atmospheric ch4 concentration (u mol cm-3)
	// const double aoc = 0.0836 * 0.001; // atmospheric ch4 concentration (u mol cm-3)
	// const double aoc = 0.0684 * 0.001; // atmospheric ch4 concentration (u mol cm-3)
        
        double ch4con[MAXN];
	double a[MAXN],b[MAXN],c[MAXN], d[MAXN], u[MAXN];
	double k[MAXN]; // conductivity
	double co[MAXN]; //concentration
	double z[MAXN]; //depth
	double df[MAXN];//diffusion

	int i;
	double ch4flx;

	nodes = int (ceil(lowb / 10.0));

	z[0] = 0.0;
	for (i = 0; i < nodes; i++)
	{
		//   if (soilt < 0.0) { ch4ratesat[i] = 0.0001;}
		df[i] = fdfs[i];
		//   if ( (df[i] < powl(10,-3)) || (df[i] > 1.0)) df[i] = 0.01;

		//  printf(" %5.3f %5.3f", ch4ratesat[i], ch4oxirate1[i]);
		//  if (ch4ratesat[i] > 5.0) ch4ratesat[i] = .010;
		//   ch4oxirate1[i] =0.01;

		//printf(" %5.3f", ch4ratesat[i]);

		z[i+1] = z[i] + x;
	}


	// k[0] = k_zero;
	k[0] = df[0] / 1.0;
	co[0] = aoc;

	for (i =1; i < nodes+1; i++)
	{
		if (i < (watertable /10.0)) u[i]=   -ch4oxirate[i] * pow(10.0,-3.0) / (3.6 * pow(10.0,3.0));  // u mol cm-3 s-1;
		else u[i] =  ch4ratesat[i] * pow(10.0,-3.0) / (3.6 * pow(10.0,3.0));  // u mol cm-3 s-1

		if (i < nodes+1) k[i] = df[i] / (z[i+1] - z[i]);
		else k[i] =0.0;
		a[i+1] = -k[i];
		b[i] = k[i-1] + k[i];
		c[i] = -k[i];
		d[i] = u[i];
	}

	d[1] = d[1] + k[0] * co[0];
	for (i=1; i < nodes; i++)
	{
		c[i] = c[i] /b[i];
		d[i] = d[i] / b[i];
		b[i+1] = b[i+1] -a[i+1]*c[i];
		d[i+1] = d[i+1] -a[i+1]*d[i];
	}

	co[nodes]= d[nodes]/b[nodes];

	for (i = nodes; i > 0; i--)
	{
		co[i] = d[i] - c[i] * co[i+1];
		ch4con[i] = co[i] / 0.001;  // converse to uM
	}

	// surface flux of methane, positive indicates flux to soil

	ch4flx = - df[1] * (co[1] - co[0]) / 1.0;  // u mol cm-2 s-1

	ch4flx = ch4flx * pow(10.0,4.0); // u mol m-2 s-1;
	ch4flx = ch4flx * pow(10.0,-6.0); // mol m-2 s-1;
	ch4flx = ch4flx * 16; // g m-2 s-1;
	ch4flx = ch4flx * 8.64 * pow(10.0,4.0); // g m-2 day-1;
	ch4flx = ch4flx * pow(10.0,3.0); // mg m-2 day-1;

	if (soilt < 0.0) { ch4flx = 0.001;}

	return ch4flx;
}




void Tsoil4::bulkdensity(double& delta_c, double& sum_depth)					//wsr: the bulk density is different between deifferent sites, so specific site date will be input into this function.

{

double ratio;
double mass = 0.0;
double depth = 0.0;												//mm
double ro;															//density distribution dunction   kg/m3
double kk = 1.0;															//accumulation step
double sum_mass = 0.0;                                                   				//total mass for a certain peat layer   g
                                                        
getdata_gasfield(sum_depth);

ratio = ratiot;

ro = rot;                                                                                //kg/m3
								
//			printf("%f,%f\n",ratio,ro);					
while ( fabs(delta_c) /ratio/0.5 - sum_mass > 1)											//g        the percision is 1g

{                  

mass =   ro * 1/10;													//g/m2

sum_mass = sum_mass + mass;														//g/m2

kk = kk + 1;																			

}

depth = kk/10 ;										//mm

if(delta_c < 0.0)

{sum_depth = sum_depth - depth;}
else
{
sum_depth = sum_depth + depth;												//mm
}

}





//         Observed bulk density, percentage of organic matter, and depth of Noname Creek site

void Tsoil4::getdata_noname(double& sum_deptht)

{
int i = 0;
 double depth[] = {25,75,125,175,225,275,325,375,425,475,525,575,625,675,725,775,825,875,925,975,1025,1075,
1125,1175,1225,1275,1325,1375,1425,1475,1525,1575,1625,1675,1725,1775,1825,1875,1925,1975,2025,2075,2125,
2175,2225,2275,2325,2375,2425,2475,2525,2575,2625,2675,2725,2775,2825,2875,2925,2975,3025,3075,3125,3175,3225,3275,3325,3375,3425};

double density[] = {28.743,69.926,84.642,83.889,77.387,81.304,67.632,64.917,66.11,66.605,73.967,124.2,153.69,167.46,256.51,146.46,
172.9,230.37,215.15,176.29,174.52,179.07,181.75,211.4,187.81,189.97,176.12,227.78,204.58,180.48,166.57,176.15,220.68,225.77,184.22,
177.87,123.86,143.01,87.514,93.374,72.406,70.28,60.138,60.72,62.649,61.22,60.39,68.846,83.046,73.317,83.373,74.613,75.805,75.808,70.387,
83.956,66.284,80.942,81.184,79.896,79.236,84.376,87.994,89.043,85.116,81.703,88.329,88.14,100.78};

double percent[] = {0.91932,0.80507,0.82097,0.86217,0.90382,0.88047,0.86199,0.87118,0.90807,0.90884,0.88649,0.58568,0.49432,
0.48197,0.34269,0.59491,0.47064,0.42022,0.41946,0.48421,0.49314,0.46923,0.44922,0.33682,0.41531,0.42021,0.44302,
0.35322,0.39076,0.48154,0.52893,0.48377,0.38331,0.38091,0.47028,0.49951,0.52814,0.52551,0.79055,0.74601,0.84253,
0.88018,0.91477,0.93007,0.92742,0.93024,0.93436,0.91527,0.89112,0.92126,0.91204,0.93988,0.93128,0.94036,0.94392,
0.94684,0.93518,0.93224,0.92671,0.93219,0.93676,0.93132,0.92485,0.88557,0.88463,0.89194,0.89311,0.904,0.8326};


while (sum_deptht > depth[i] )

{

i++;

}


if (i <= 1 ) {ratiot = percent[0];    rot = density[0]; }


ratiot = percent[i]; 

rot = density[i] ;


}

//         Observed bulk density, percentage of organic matter, and depth of Gasfield site

void Tsoil4::getdata_gasfield(double& sum_deptht)

{
int i = 0;
 double depth[] = {39.06,97.136,154.65,215.59,292.88,374.63,454.75,526.54,592.92,669.58,739.36,
811.36,867.36,922.92,978.47,1034,1089.6,1145.1,1200.7,1256.3,1311.8,1367.4,1423.1,1479.4,1535.6,
1591.9,1648.1,1704.4,1760.6,1816.9,1873.1,1929.4,1985.6,2041.9,2098.1,2154.4,2210.6,2266.9,2323.2,
2380.7,2438.4,2496.1,2553.8,2611.4,2669.1,2726.8,2784.5};

double density[] = {132,138.4,145,194.2,116.8,170.8,215.6,208.6,234,175.2,169.6,113.6,128,149.8,262.6,368.6,
366.6,356.8,340,403.4,271.6,234.6,229.8,226.6,236.2,203,231,291,363.2,406.6,331.6,243.2,269.4,300.8,252.8,
215.8,178,184.8,190.4,174.8,171.2,176.2,175.2,193,218,215.6,295.2};

double percent[] = {0.7886,0.7849,0.85691,0.88439,0.91474,0.86218,0.87484,0.7479,0.73422,0.88166,0.88714,
0.88914,0.84711,0.85595,0.63225,0.4614,0.45723,0.44765,0.47638,0.41368,0.64649,0.74049,0.74229,0.7583,
0.70406,0.75254,0.70263,0.58277,0.45469,0.38666,0.48905,0.63586,0.62083,0.57227,0.70011,0.75368,
0.86384,0.87909,0.91183,0.9079,0.90347,0.89412,0.90678,0.88387,0.75451,0.75986,0.58268};



while (sum_deptht > depth[i])

{

i++;
}


if (i <= 1 ) {ratiot = percent[0];    rot = density[0]; }

ratiot = (percent[i-1] + percent[i] ) /2;

rot = (density[i-1] + density[i] ) /2;

}



//         Observed bulk density, percentage of organic matter, and depth of Horse Trail site

void Tsoil4::getdata_horse(double& sum_deptht)

{
int i = 0;
 double depth[] = {50,160,240,320,400,480,560,640,720,800,880,960,1040,1125,
1220,1300,1380,1475,1560,1640,1720,1800,1880,1960,2040,2120,2200,2280,2360,2440,2520,
2600,2680,2760,2840,2920,2995,3040,3080,3120,3160,3200,3240,3280,3320,3360,3400,3440,
3480,3520,3560,3600,3640,3680,3720,3760,3800,3840,3880,3920,3960,4000,4040,4080,4120};

double density[] = {27.175,68.375,58.65,79.725,85.775,86.925,104.7,67.875,93.725,74.55,66.4,
111.03,158.22,210.53,205.55,156.4,102.8,88.2,78.125,104.35,99.525,70.175,65.65,67.425,
81.425,86.075,86.8,97.125,90.9,108.58,111.75,145.18,112.13,139.4,147.65,131.95,94.15,57.475,
63.25,57.7,55.975,45.55,42.225,49.95,57.625,63.747,69.45,67.925,90.15,65.15,
90.75,87.725,87.625,78.65,55.35,84.825,73.125,80.55,82.1,113.75,155.93,175.35,245.3,275.77,228.28};

double percent[] = {0.71094,0.76383,0.71903,0.60398,0.72585,0.72284,0.68059,0.78446,0.77455,0.79782,0.67775,
0.51081,0.38903,0.38275,0.42771,0.52646,0.67496,0.67828,0.67931,0.67869,0.79758,0.84713,0.85085,0.89107,
0.84609,0.81633,0.82534,0.86055,0.83219,0.86819,0.82636,0.85582,0.81764,0.86242,0.85448,0.85224,0.77941,
0.72154,0.71398,0.76657,0.75329,0.81042,0.79105,0.7193,0.71413,0.65535,0.59431,0.60091,0.5388,0.6557,
0.53716,0.53055,0.51162,0.49392,0.57576,0.49274,0.68113,0.32348,0.41556,0.2929,0.31523,0.21958,0.15436,0.14055,0.15711};




while (sum_deptht > depth[i])

{

i++;
}


if (i <= 1 ) {ratiot = percent[0];    rot = density[0]; }

ratiot = (percent[i-1] + percent[i] ) /2;

rot = (density[i-1] + density[i] ) /2;

}






//         Observed bulk density, percentage of organic matter, and depth of Swanson site

void Tsoil4::getdata_swanson(double& sum_deptht)

{
int i = 0;
 double depth[] = {90,130,170,210,250,290,330,370,410,450,490,530,570,610,650,690,730,770,810,
850,890,930,970,1010,1050,1090,1130,1170,1210,1250,1290,1330,1370,1410,1450,1490,
1530,1570,1610,1650,1690,1730,1770,1810,1850,1890,1930,1970,2010,2050,2090,2130,2170,2210,2250,2290,2330,2370,2410,2450};

double density[] = {62.15,247.95,100.85,98.25,159,157.2,157.5,153.85,201.85,102.7,149.6,
210.4,191.35,132.9,162.85,129,85.15,115.75,129.5,126.3,165.95,105.9,107.4,
111.25,103.35,95.8,85.1,99.6,85.15,56.45,102.55,88.75,105.7,79.4,98.7,64,103,96.05,
172.15,146.8,104.3,86.3,118.85,117.2,103.95,76.8,201.3,208.65,336.65,260.5,260.25,
230.65,220.8,179.05,258.65,266,248.75,266.8,227,178};

double percent[] = {0.90587,0.97963,0.81309,0.82901,0.71415,0.63104,0.8546,0.67468,0.30642,
0.7074,0.68917,0.6856,0.6023,0.69827,0.73994,0.79535,0.92543,0.85227,0.91313,0.76326,
0.67671,0.80595,0.88268,0.76944,0.81132,0.83142,0.78613,0.87851,0.85966,0.9194,0.93125,
0.95268,0.89688,0.95025,0.9382,0.94922,0.92718,0.91983,0.90822,0.89271,0.46117,0.89977,
0.87169,0.86433,0.63781,0.5013,0.8157,0.46681,0.42581,0.57255,0.52104,
0.51312,0.46694,0.48981,0.49758,0.5282,0.52844,0.55454,0.52599,0.52247};




while (sum_deptht > depth[i])

{

i++;
}


if (i <= 1 ) {ratiot = percent[0];    rot = density[0]; }

ratiot = (percent[i-1] + percent[i] ) /2;

rot = (density[i-1] + density[i] ) /2;

}










