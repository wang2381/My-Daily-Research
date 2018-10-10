/* **************************************************************
ODEINT423.CPP - object contains functions to integrate ordinary
	       differential equations (ODE)
*****************************************************************

Modifications:

19991102 - DWK deleted code for:
           virtual int boundcon(double ptstate[], double err[], double& ptol)
           virtual void delta(int& pdm, double pstate[], double pdstate[])

20000102 - DWK added compiler directives
************************************************************* */

#if !defined(ODEINT423_H)
  #include "odeint423.hpp"
#endif

/* *************************************************************
************************************************************* */

Odeint4::Odeint4()
{

  syint = 1;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */





int Odeint4::adapt(const int& numeq, double pstate[NVT][NUMEQ], double& ptol, const int& pdm, const int& vt)
{

  int i;
  double ipart;
  double fpart;
  double time = 0.0;
  double dt = 1.0;
  int mflag = 0;
  long nintmon = 0;
  double oldstate[NVT][NUMEQ];
  
  

	rkf(numeq,pstate,dt,pdm,vt);

	//test = boundcon(dum4,error,ptol,vt);
       test = ACCEPT;
	 mflag = 1;

        if (test == ACCEPT)
        {
         for(i = 0; i < numeq;i++) { pstate[vt][i] = dum4[vt][i]; }

       }







  return mflag;

};
  
  
  
/*
int Odeint4::adapt(const int& numeq, double pstate[NVT][NUMEQ], double& ptol, const int& pdm, const int& vt)
{

  int i;
  double ipart;
  double fpart;
  double time = 0.0;
  double dt = 1.0;
  int mflag = 0;
  long nintmon = 0;
  double oldstate[NVT][NUMEQ];


  blackhol = 0;
  while (time != 1.0)
  {
    test = REJECT;
    if (syint == 1)
    {
      while (test != ACCEPT)
      {
	rkf(numeq,pstate,dt,pdm,vt);
	test = boundcon(dum4,error,ptol,vt);
	if (dt <= pow(0.5,maxit))
        {
	  test = ACCEPT;
	  mflag = 1;
          if (nintmon == 0)
          {
            for(i = 0; i < numeq;i++) { oldstate[vt][i] = pstate[vt][i]; }
          }
	  ++nintmon;
	}

        if (test == ACCEPT)
        {
          for(i = 0; i < numeq;i++) { pstate[vt][i] = dum4[vt][i]; }
          time += dt;
          fpart = modf((0.01 + (time/(2.0*dt))),&ipart);
          if ( fpart < 0.1 && dt < 1.0) { dt *= 2.0; }
        }
        else { dt *= 0.500; }

        if (nintmon == maxitmon)
        {
          time = 1.0;
          blackhol = 1;
          for(i = 0; i < numeq;i++) { pstate[vt][i] = oldstate[vt][i]; }
        }
      }
    }    
 }      

	//delta(pdm, pstate, f11,vt);
	//step(numeq, pstate, f11, pstate, dt,vt);

//	rkf(numeq,pstate,dt,pdm,vt);
	//test = boundcon(dum4,error,ptol,vt);
       //test = ACCEPT;
//	 mflag = 1;

    //    if (test == ACCEPT)
       // {
       //  for(i = 0; i < numeq;i++) { pstate[vt][i] = dum4[vt][i]; }

       //}







  return mflag;

};
*/
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Odeint4::ask(ofstream& rflog1)
{


/* **************************************************************
	      Parameters for Adaptive Integrator
************************************************************** */

  cout << endl << "Enter the proportional tolerance for the integrator: ";
  fpara >> inittol;

  rflog1 << endl << "Enter the proportional tolerance for the integrator: " << inittol << endl;

  cout << "Enter the maximum number of iterations in the integrator: ";
  fpara >> maxit;

  rflog1 << "Enter the maximum number of iterations in the integrator: " << maxit << endl;

  cout << "Enter the maximum number of times in a month that the" << endl;
  cout << "integrator can reach the maximum number of iterations: ";
  fpara >> maxitmon;

  rflog1 << "Enter the maximum number of times in a month that the" << endl;
  rflog1 << "integrator can reach the maximum number of iterations: " << maxitmon << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************ */

void Odeint4::rkf(const int& numeq, double pstate[NVT][NUMEQ], double& pdt, const int& pdm, const int& vt)
{

  int i,j;
  double ptdt = 0;

  for (i = 0; i < numeq;i++)
  {
    dum4[vt][i] = dum5[vt][i] = pstate[vt][i];
    yprime[vt][i] = rk45[vt][i] = error[vt][i] = 0.0;
  }

  ptdt = pdt * 0.25 ;
 
	// dum4 = pstate
	// f11 = k1 = f(t_j, y_j)
  delta(pdm,dum4,f11,vt);//printf("%f\n",f11[vt][26]);
  	// yprime = yprime + a1*k1 = 0 + a1*k1 = a1*k1
	// y_j+1 = y_j + a1*k1 = 0 + a1*k1 = a1*k1
	
  //printf("%f\n",dum4[vt][26]);

  step(numeq,yprime,f11,yprime,a1,vt);//printf("%f\n",f11[vt][26]);
  
  	// rk45 = rk45 + b1*k1 = 0 + b1*k1 = b1*k1
	// z_j+1 = z_j + b1*k1 = 0 + b1*k1 = b1*k1
  step(numeq,rk45,f11,rk45,b1,vt);
  
  	// ydum = dum4 + ptdt*f11 = y_j + 1/4*pdt*k1
  step(numeq,dum4,f11,ydum,ptdt,vt);
  
  	// f2 = k2
  delta(pdm,ydum,f2,vt);
  
  	// f13 = a31*k1 + a32*k2
  for (i = 0; i < numeq; i++)
  {
    f13[vt][i] = a31*f11[vt][i] + a32*f2[vt][i];
  }
  
  	// ydum = dum4 + 1*f13 = y_j + a31*k1 + a32*k2
 step(numeq,dum4,f13,ydum,pdt,vt);
 
 	// f3 = k3
  delta(pdm,ydum,f3,vt);
  
  	// yprime = yprime + a3*f3 = a1*k1 + a3*k3
	// y_j+1 = y_j + a1*k1 + a3*k3
  step(numeq,yprime,f3,yprime,a3,vt);
  //for (i = 0; i < numeq; i++)
  //{
  //if (yprime[vt][i] == 0.0 )    { yprime[vt][i] = 0.1;   }
  
  //}
  	// rk45 = rk45 + b3*f3 = b1*k1 + b3*k3
	// z_j+1 = z_j + b3*f3 = b1*k1 + b3*k3
  step(numeq,rk45,f3,rk45,b3,vt);
  
  	// f14 = a41*f11 + a42*f2 + a43*f3 = a41*k1 + a42*k2 + a43*k3
  for (i = 0; i < numeq; i++)
  {
    f14[vt][i] = a41*f11[vt][i] + a42*f2[vt][i] + a43*f3[vt][i];
  }
  
  //	 ydum = dum4 + 1*f14 = y_j + a41*k1 + a42*k2 + a43*k3
  step(numeq,dum4,f14,ydum,pdt,vt);

  	//f4 = k4
  delta(pdm,ydum,f4,vt);

  	// yprime = yprime + a4*f4 = a1*k1 + a3*k3 + a4 * k4
	// y_j+1 = y_j + a4*f4 = a1*k1 + a3*k3 + a4*k4
  step(numeq,yprime,f4,yprime,a4,vt);
    //for (i = 0; i < numeq; i++)
  //{
  //if (yprime[vt][i] == 0.0 || yprime[vt][i] > 1000 || yprime[vt][i] < -99)    { yprime[vt][i] = 0.1;   printf("12%d\n",i);}
  
  //}
  	// rk45 = rk45 + b4*f4 = b1*k1 + b3*k3 + b4*k4
	// z_j+1 = z_j + b4*f4 = b1*k1 + b3*k3 + b4*k4
  step(numeq,rk45,f4,rk45,b4,vt);
  
  	// f15 = a51*f11 + a52*f2 + a53*f3 + a54* f4 = a51*k1 + a52*k2 + a53*k3 + a54*k4
  for (i = 0; i < numeq; i++)
  {
    f15[vt][i] = a51*f11[vt][i] + a52*f2[vt][i] + a53*f3[vt][i] + a54*f4[vt][i];
  }
  
  	// ydum = dum4 + 1*f15 = y_j + a51*k1 + a52*k2 + a53*k3 + a54*k4
  step(numeq,dum4,f15,ydum,pdt,vt);
  
  	// f5 = k5
  delta(pdm,ydum,f5,vt);      
 // printf("%f\n",ydum[vt][26]);
 
 	// yprime = yprime + a5*f5 = a1*k1 + a3*k3 + a4*k4 + a5*k5
	// y_j+1 = y_j + a5*f5 = a1*k1 + a3*k3 + a4*k4 + a5*k5
  step(numeq,yprime,f5,yprime,a5,vt);
  	// rk45 = rk45 + b5*f5 = b1*k1 + b3*k3 + b4*k4 + b5*k5
	// z_j+1 = z_j + b5*f5 = b1*k1 + b3*k3 + b4*k4 + b5*k5
  step(numeq,rk45,f5,rk45,b5,vt);
  	// f16 = b61*f11 + b62*f2 + b63*f3 + b64*f4 + b65*f5 = b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5
  for (i = 0; i < numeq; i++)
  {
    f16[vt][i] = b61*f11[vt][i] + b62*f2[vt][i] + b63*f3[vt][i] + b64*f4[vt][i] + b65*f5[vt][i];
  }
  

	// ydum = dum4 + 1*f16 = y_j + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5
  step(numeq,dum4,f16,ydum,pdt,vt);
  
  	// f6 = k6
  delta(pdm,ydum,f6,vt);
  
  	// rk45 = rk45 + b6*f6 = b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6
	// z_j+1 = z_j + b6*f6 = b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6
  step(numeq,rk45,f6,rk45,b6,vt);
  	// dum4 = dum4 + pdt*yprime = y_j + yprime = y_j + a1*k1 + a3*k3 + a4*k4 + a5*k5
  step(numeq,dum4,yprime,dum4,pdt,vt);
  	// dum5 = dum5 + pdt*rk45 = z_j + rk45 = z_j + b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6
  step(numeq,dum5,rk45,dum5,pdt,vt);
 // printf("%f\n",ydum[vt][26]);
 
 	// Errors between 4-order and 5-order estimates
    for (i = 0; i < numeq; i++)
  {
    error[vt][i] = fabs(dum4[vt][i] - dum5[vt][i]);
  
	}
	
	
//	printf("   %f       %f         %f      %f\n",error[vt][26],ydum[vt][26],dum4[vt][26],dum5[vt][26]);

};







/*

void Odeint4::rkf(const int& numeq, double pstate[NVT][NUMEQ], double& pdt, const int& pdm, const int& vt)
{

  int i,j;
  double ptdt = 0;

  for (i = 0; i < numeq;i++)
  {
    dum4[vt][i] = dum5[vt][i] = pstate[vt][i];
    yprime[vt][i] = rk45[vt][i] = error[vt][i] = 0.0;
  }

  ptdt = pdt * 0.25 ;
 
	// dum4 = pstate
	// f11 = k1 = f(t_j, y_j)
  delta(pdm,dum4,f11,vt);//printf("%f\n",f11[vt][26]);
  	// yprime = yprime + a1*k1 = 0 + a1*k1 = a1*k1
	// y_j+1 = y_j + a1*k1 = 0 + a1*k1 = a1*k1
	
  //printf("%f\n",dum4[vt][26]);

  step(numeq,yprime,f11,yprime,a1,vt);//printf("%f\n",f11[vt][26]);
  
  
  	// ydum = dum4 + ptdt*f11 = y_j + 1/4*pdt*k1
  
  step(numeq,dum4,f11,ydum,ptdt,vt);
  
  	// f2 = k2
  delta(pdm,ydum,f2,vt);
  
  	// f13 = a31*k1 + a32*k2
  for (i = 0; i < numeq; i++)
  {
    f13[vt][i] = a31*f11[vt][i] + a32*f2[vt][i];
 // printf("%d\n",i);  
    
  }
  
  	// ydum = dum4 + 1*f13 = y_j + a31*k1 + a32*k2
 step(numeq,dum4,f13,ydum,pdt,vt);
 
 	// f3 = k3
  delta(pdm,ydum,f3,vt);
  
  	// yprime = yprime + a3*f3 = a1*k1 + a3*k3
	// y_j+1 = y_j + a1*k1 + a3*k3
  step(numeq,yprime,f3,yprime,a3,vt);

  step(numeq,dum4,yprime,dum4,pdt,vt);

	
	


};


*/ 
/***************************************************************
 ***************************************************************/


void Odeint4::step(const int& numeq, double pstate[NVT][NUMEQ], double pdstate[NVT][NUMEQ], double ptstate[NVT][NUMEQ],
			 double& pdt, const int& vt)
{

  for (int i = 0; i < numeq; i++)
  {
    ptstate[vt][i] = pstate[vt][i] + (pdt * pdstate[vt][i]);
  }
     // printf("%f\n",ptstate[vt][26]);
};

