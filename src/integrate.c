/*
Use this code for integrating out the equations outwards with no mathcing. This code also implements Roxana's version of the intergration
routine, which is identical to mine, but has proven to be useful in the debugging process.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NRANSI
#define LEN             1000


typedef struct model
  {
  double        *r;
  double        *Phi;
  double        *X;
  double        *A;
  double        *eta;
  double        *m;
  double        *R;     // Isotropic radius
  double        *f;     // f := R / r
  double        *psi;   // conformal factor, i.e. g_rr = psi^4
  double        *C;     // Compactness defined as m(r)/r
  double        *Psi0;
  double        *phigrav;
  } model;


typedef struct param
  {
  double        rmax;
  int           nint;
  double        omega0;
  double        mphigrav0;
  double        A0;
  double        alpha0;
  double        beta0;
  double        phigrav0;
  double        Phi0;
  int           nzerotarget;
  double        thresh;
  double        thresh_grav;
  double        mpercentage;
  double        rmatchfac;
  double        lambda4;
  double        lambda6;
  double        lambda8;
  double        sigma0;
  char          potential[LEN];
  char          minmax[LEN];
  } param;


// Functions
void    calcIso         ();
void    calcMass        ();
double  calcRadius      ();
int  calcRadius_index      ();
double  findMax         (double*, int, int);
void    initGrid        (int, double);
void    intODE          (double, double, double, int*, int*, double*, double*, int*);
void    intODE_roxana          (double, double, double, int*, int*, double*, double*, int*);
int     mygetline       (char*, FILE*);
void    out1D           (double*, double*, int, int, const char*);
void    printPars       ();
void    readPars        (char*);
void    registerPotential();
void    rescalePhi      ();
void    rhsBSint        (double*, double*, double*, double*, double*, double*, double, double, double, double, double, double, double, double);
double    rhsBSint_X        (double, double, double, double, double, double, double, double);
double    rhsBSint_phigrav        (double, double, double, double, double, double, double, double);
double    rhsBSint_phi       (double, double, double, double, double, double, double, double);
double    rhsBSint_psi0        (double, double, double, double, double, double, double, double);
double    rhsBSint_A        (double, double, double, double, double, double, double, double);
double    rhsBSint_eta        (double, double, double, double, double, double, double, double);
void    rhsIso          (double*, double, double, double, double);
int     iofr            (double);
void    integrate_outwards ();
double  V_series        (double);
double  Vp_series       (double);
double  V_solitonic     (double);
double  Vp_solitonic    (double);
double  F_st               (double);
double  derF_st            (double);
double  W_st               (double);
double  derW_st            (double);
double alpha            (double, double);


// Function pointers
double  (*V)            (double);
double  (*Vp)           (double);
double  (*W)            (double);
double  (*derW)           (double);
double  (*F)            (double);
double  (*derF)           (double);


// Global variables
const double PI = 3.1415926535897932385;
param par;
model star;


/*==========================================================================*/

int main(int argc, char* argv[])
  {
  int i, k;
  int success, converged;
  double omBS, mBS, rBS, CBS;

  if(argc != 2) { printf("Usage:   %s   <parfile>\n\n", argv[0]); exit(0); }

  // printf("I'm reading params file \n");
  // Parameters
  readPars(argv[1]);
  printPars();

  // Register functions
  registerPotential();

  // Allocate memory
  star.r   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.Phi = (double*) malloc((size_t) par.nint * sizeof(double));
  star.X   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.A   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.Psi0   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.phigrav   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.eta = (double*) malloc((size_t) par.nint * sizeof(double));
  star.m   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.R   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.f   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.psi = (double*) malloc((size_t) par.nint * sizeof(double));
  star.C   = (double*) malloc((size_t) par.nint * sizeof(double));  

  printf("Gravitational scalar field guess is set to %22.16g\n", par.phigrav0);

  // Shoot 
  integrate_outwards();

  // IO
  out1D(star.r, star.X, 0, par.nint-1, "X.dat");
  out1D(star.r, star.phigrav, 0, par.nint-1, "phigrav.dat");
  out1D(star.r, star.Psi0, 0, par.nint-1, "Psi0.dat");
  out1D(star.r, star.A, 0, par.nint-1, "A.dat");
  out1D(star.r, star.eta, 0, par.nint-1, "eta.dat");
  out1D(star.r, star.Phi, 0, par.nint-1, "Phi.dat");
  out1D(star.r, star.m, 0, par.nint-1, "m.dat");
  out1D(star.r, star.R, 0, par.nint-1, "RisoofR.dat");
  out1D(star.R, star.f, 0, par.nint-1, "fofriso.dat");
  //out1D(star.R, star.psi, 0, par.nint-1, "psiofriso.dat");
  out1D(star.R, star.A, 0, par.nint-1, "Aofriso.dat");

  printf("\n===================================================\n");
  printf("Physical frequency:   omega = %22.16g\n", par.omega0);
  printf("Total mass:           m     = %22.16g\n", star.m[par.nint-1]);
  printf("===================================================\n");
  }

/*==========================================================================*/

void integrate_outwards()
  {
  int    cnt, success;
  double om, ommin, ommax;
  double phigrav, phigravmin, phigravmax;
  int    nzero, phigravnzero, sigA;
  double rstop, rstopgrav;

  // Grid setup
  initGrid(par.nint, par.rmax);

  intODE(par.A0, par.phigrav0, par.omega0, &nzero, &phigravnzero, &rstop, &rstopgrav, &sigA);
  }

/*==========================================================================*/

void intODE(double A0, double phigrav0, double omega, int* nzero, int* phigravnzero, double* rstop, double* rstopgrav, int* sigAstop)
  {
  int    i, n1, istop, phigravstop;
  double dr, r, X, A, eta, Phi, Psi0, phigrav, om;
  double rhs_X, rhs_A, rhs_eta, rhs_Phi, rhs_Psi0, rhs_phigrav;
  double dX[5], dA[5], dPsi0[5], deta[5], dphigrav[5], dPhi[5];

  n1 = par.nint;
  om = omega;

  *nzero = 0;
  *phigravnzero = 0;
  *rstop = star.r[n1-1];
  *rstopgrav = star.r[n1-1];
  istop  = -1;   // So we have no vacuum region unless we exceed the amplitude
  phigravstop = -1;
  
  // Central values
  star.X[0]   = 1/sqrt(F(phigrav0));
  star.A[0]   = A0;
  star.phigrav[0]   = phigrav0;
  star.Psi0[0]   = 0;
  star.eta[0] = 0;
  star.Phi[0] =  par.Phi0; // Phi has a free constant we later match to Schwarzschild
  //-0.3679003437995622;
  //-2.9658829262999253e-01;  
  //-4.6694715920249441e-01 

  double alpha_value = exp(star.Phi[0])/sqrt(F(star.phigrav[0]));
  printf("\n The value of alpha[0] is: %g\n", alpha_value);

  // RK integration
  for(i = 1; i < n1; i++)
    {
    dr  = star.r[i] - star.r[i-1];

    // 1st RK step
    r   = star.r[i-1];
    X   = star.X[i-1];
    A   = star.A[i-1];
    eta = star.eta[i-1];
    phigrav = star.phigrav[i-1];
    Phi = star.Phi[i-1];
    Psi0 = star.Psi0[i-1];
    rhsBSint(&rhs_X, &rhs_A, &rhs_eta, &rhs_Phi, &rhs_phigrav, &rhs_Psi0, r, X, A, eta, Phi, phigrav, Psi0, om);
    dX[1]   = rhs_X * dr;
    dA[1]   = rhs_A * dr;
    dPsi0[1]   = rhs_Psi0 * dr;
    dphigrav[1] = rhs_phigrav * dr;
    deta[1] = rhs_eta * dr;
    dPhi[1] = rhs_Phi * dr;

    // 2nd RK step
    r   = star.r[i-1] + 0.5 * dr;
    X   = star.X[i-1] + 0.5 * dX[1];
    A   = star.A[i-1] + 0.5 * dA[1];
    eta = star.eta[i-1] + 0.5 * deta[1];
    phigrav = star.phigrav[i-1] +  0.5 * dphigrav[1];
    Phi = star.Phi[i-1] + 0.5 * dPhi[1];
    Psi0 = star.Psi0[i-1] +  0.5 * dPsi0[1];
    rhsBSint(&rhs_X, &rhs_A, &rhs_eta, &rhs_Phi, &rhs_phigrav, &rhs_Psi0, r, X, A, eta, Phi, phigrav, Psi0, om);
    dX[2]   = rhs_X * dr;
    dA[2]   = rhs_A * dr;
    dPsi0[2]   = rhs_Psi0 * dr;
    dphigrav[2] = rhs_phigrav * dr;
    deta[2] = rhs_eta * dr;
    dPhi[2] = rhs_Phi * dr;

    // 3rd RK step
    r   = star.r[i-1] + 0.5 * dr;
    X   = star.X[i-1] + 0.5 * dX[2];
    A   = star.A[i-1] + 0.5 * dA[2];
    eta = star.eta[i-1] + 0.5 * deta[2];
    phigrav = star.phigrav[i-1] +  0.5 * dphigrav[2];
    Phi = star.Phi[i-1] + 0.5 * dPhi[2];
    Psi0 = star.Psi0[i-1] +  0.5 * dPsi0[2];
    rhsBSint(&rhs_X, &rhs_A, &rhs_eta, &rhs_Phi, &rhs_phigrav, &rhs_Psi0, r, X, A, eta, Phi, phigrav, Psi0, om);
    dX[3]   = rhs_X * dr;
    dA[3]   = rhs_A * dr;
    dPsi0[3]   = rhs_Psi0 * dr;
    dphigrav[3] = rhs_phigrav * dr;
    deta[3] = rhs_eta * dr;
    dPhi[3] = rhs_Phi * dr;

    // 4th RK step
    r   = star.r[i];
    X   = star.X[i-1] + dX[3];
    A   = star.A[i-1] + dA[3];
    eta = star.eta[i-1] + deta[3];
    phigrav = star.phigrav[i-1] + dphigrav[3];
    Phi = star.Phi[i-1] + dPhi[3];
    Psi0 = star.Psi0[i-1] + dPsi0[3];
    rhsBSint(&rhs_X, &rhs_A, &rhs_eta, &rhs_Phi, &rhs_phigrav, &rhs_Psi0, r, X, A, eta, Phi, phigrav, Psi0, om);
    dX[4]   = rhs_X * dr;
    dA[4]   = rhs_A * dr;
    dPsi0[4]   = rhs_Psi0 * dr;
    dphigrav[4] = rhs_phigrav * dr;
    deta[4] = rhs_eta * dr;
    dPhi[4] = rhs_Phi * dr;

    // Update variables
    star.X[i]   = star.X[i-1]   + (dX[1]  + 2*dX[2]  + 2*dX[3]  + dX[4] ) / 6.0;
    star.A[i]   = star.A[i-1]   + (dA[1]  + 2*dA[2]  + 2*dA[3]  + dA[4] ) / 6.0;
    star.Psi0[i]   = star.Psi0[i-1]   + (dPsi0[1]  + 2*dPsi0[2]  + 2*dPsi0[3]  + dPsi0[4] ) / 6.0;
    star.phigrav[i]   = star.phigrav[i-1]   + (dphigrav[1]  + 2*dphigrav[2]  + 2*dphigrav[3]  + dphigrav[4] ) / 6.0;
    star.eta[i] = star.eta[i-1] + (deta[1]+ 2*deta[2]+ 2*deta[3]+deta[4]) / 6.0;
    star.Phi[i] = star.Phi[i-1] + (dPhi[1]+ 2*dPhi[2]+ 2*dPhi[3]+dPhi[4]) / 6.0;

    // Analyze: do we cross zero? do we exceed 2*A0?
    if(fabs(star.A[i]) > 2*star.A[0] || star.A[i] != star.A[i])
      {
      *sigAstop = ( (star.A[i-1] > 0) - (star.A[i-1] < 0) );
      istop = i - 1;   // We stop the integration; one point as sanity buffer
      // There is one problem here: We regard the present data point and
      // also the previous one as unrealiable. The previous point,
      // however, might have been counted as a zero crossing. If that
      // happened, we wish to undo this counting.
      if(star.A[i-1] * star.A[i-2] < 0) (*nzero)--;
      break;
      }

    if(star.A[i] * star.A[i-1] < 0) (*nzero)++;   // Zero crossing found
    
    if(fabs(star.phigrav[i]) > 4*star.phigrav[0] || star.phigrav[i] != star.phigrav[i])
      {
        // *phigravsig = ( (star.phigrav[i-1] > 0) - (star.phigrav[i-1] < 0) );
        phigravstop = i - 1;   // We stop the integration; one point as sanity buffer
        // printf("phigravstop: %22.16g \n", star.r[phigravstop]);
        if(star.phigrav[i-1] * star.phigrav[i-2] < 0) 
        {(*phigravnzero)--;}
        // printf("value is %g\n", star.r[i-1]);
        // break;
        // }
      }

      if(star.phigrav[i] * star.phigrav[i-1] < 0) (*phigravnzero)++;   // Zero crossing found
    }

  if(istop < 0)
    {
    *sigAstop = ( (star.A[n1-1] > 0) - (star.A[n1-1] < 0) );
    printf("No need for vacuum \n");
    return;   // Integration went through; no need for vacuum
    }

  // Set to vacuum beyond rstop
  *rstop = star.r[istop];
  *rstopgrav = star.r[phigravstop];

  }

/*==========================================================================*/

void intODE_roxana(double A0, double phigrav0, double omega, int* nzero, int* phigravnzero, double* rstop, double* rstopgrav, int* sigAstop)
  {
  int    i, n1, istop, phigravstop;
  double dr, r, X, A, eta, Phi, Psi0, phigrav, om;

  double  kk1, kk2, kk3, kk4;
  double  ll1, ll2, ll3, ll4;
  double  mm1, mm2, mm3, mm4;
  double  nn1, nn2, nn3, nn4;
  double  ss1, ss2, ss3, ss4;
  double  qq1, qq2, qq3, qq4;

  n1 = par.nint;
  om = omega;

  *nzero = 0;
  *phigravnzero = 0;
  *rstop = star.r[n1-1];
  *rstopgrav = star.r[n1-1];
  istop  = -1;   // So we have no vacuum region unless we exceed the amplitude
  phigravstop = -1;
  
  // Central values
  star.X[0]   = 1/sqrt(F(phigrav0));
  star.A[0]   = A0;
  star.phigrav[0]   = phigrav0;
  star.Psi0[0]   = 0;
  star.eta[0] = 0;
  star.Phi[0] =  par.Phi0; // Phi has a free constant we later match to Schwarzschild
  //-0.3679003437995622;
  //-2.9658829262999253e-01;  
  //-4.6694715920249441e-01 
  
  double alpha_value = exp(star.Phi[0])/sqrt(F(star.phigrav[0]));
  printf("\n The value of alpha[0] is: %g\n", alpha_value);

  // RK integration
  for (i = 1; i < n1; i++)
        {
                dr  = star.r[i] - star.r[i-1];

                kk1=dr*rhsBSint_phi(star.r[i-1], star.X[i-1], star.A[i-1], star.eta[i-1], star.Phi[i-1], star.phigrav[i-1], star.Psi0[i-1], om);
                ll1=dr*rhsBSint_X(star.r[i-1], star.X[i-1], star.A[i-1], star.eta[i-1], star.Phi[i-1], star.phigrav[i-1], star.Psi0[i-1], om);
                nn1=dr*rhsBSint_psi0(star.r[i-1], star.X[i-1], star.A[i-1], star.eta[i-1], star.Phi[i-1], star.phigrav[i-1], star.Psi0[i-1], om);
                ss1=dr*rhsBSint_eta(star.r[i-1], star.X[i-1], star.A[i-1], star.eta[i-1], star.Phi[i-1], star.phigrav[i-1], star.Psi0[i-1], om);
                mm1=dr*rhsBSint_phigrav(star.r[i-1], star.X[i-1], star.A[i-1], star.eta[i-1], star.Phi[i-1], star.phigrav[i-1], star.Psi0[i-1], om);
                qq1=dr*rhsBSint_A(star.r[i-1], star.X[i-1], star.A[i-1], star.eta[i-1], star.Phi[i-1], star.phigrav[i-1], star.Psi0[i-1], om);

                kk2=dr*rhsBSint_phi(star.r[i-1]+dr/2, star.X[i-1]+ll1/2, star.A[i-1]+qq1/2, star.eta[i-1]+ss1/2, star.Phi[i-1]+kk1/2, star.phigrav[i-1]+mm1/2, star.Psi0[i-1]+nn1/2, om);
                ll2=dr*rhsBSint_X(star.r[i-1]+dr/2, star.X[i-1]+ll1/2, star.A[i-1]+qq1/2, star.eta[i-1]+ss1/2, star.Phi[i-1]+kk1/2, star.phigrav[i-1]+mm1/2, star.Psi0[i-1]+nn1/2, om);
                nn2=dr*rhsBSint_psi0(star.r[i-1]+dr/2, star.X[i-1]+ll1/2, star.A[i-1]+qq1/2, star.eta[i-1]+ss1/2, star.Phi[i-1]+kk1/2, star.phigrav[i-1]+mm1/2, star.Psi0[i-1]+nn1/2, om);
                ss2=dr*rhsBSint_eta(star.r[i-1]+dr/2, star.X[i-1]+ll1/2, star.A[i-1]+qq1/2, star.eta[i-1]+ss1/2, star.Phi[i-1]+kk1/2, star.phigrav[i-1]+mm1/2, star.Psi0[i-1]+nn1/2, om);
                mm2=dr*rhsBSint_phigrav(star.r[i-1]+dr/2, star.X[i-1]+ll1/2, star.A[i-1]+qq1/2, star.eta[i-1]+ss1/2, star.Phi[i-1]+kk1/2, star.phigrav[i-1]+mm1/2, star.Psi0[i-1]+nn1/2, om);
                qq2=dr*rhsBSint_A(star.r[i-1]+dr/2, star.X[i-1]+ll1/2, star.A[i-1]+qq1/2, star.eta[i-1]+ss1/2, star.Phi[i-1]+kk1/2, star.phigrav[i-1]+mm1/2, star.Psi0[i-1]+nn1/2, om);

                kk3=dr*rhsBSint_phi(star.r[i-1]+dr/2, star.X[i-1]+ll2/2, star.A[i-1]+qq2/2, star.eta[i-1]+ss2/2, star.Phi[i-1]+kk2/2, star.phigrav[i-1]+mm2/2, star.Psi0[i-1]+nn2/2, om);
                ll3=dr*rhsBSint_X(star.r[i-1]+dr/2, star.X[i-1]+ll2/2, star.A[i-1]+qq2/2, star.eta[i-1]+ss2/2, star.Phi[i-1]+kk2/2, star.phigrav[i-1]+mm2/2, star.Psi0[i-1]+nn2/2, om);
                nn3=dr*rhsBSint_psi0(star.r[i-1]+dr/2, star.X[i-1]+ll2/2, star.A[i-1]+qq2/2, star.eta[i-1]+ss2/2, star.Phi[i-1]+kk2/2, star.phigrav[i-1]+mm2/2, star.Psi0[i-1]+nn2/2, om);
                ss3=dr*rhsBSint_eta(star.r[i-1]+dr/2, star.X[i-1]+ll2/2, star.A[i-1]+qq2/2, star.eta[i-1]+ss2/2, star.Phi[i-1]+kk2/2, star.phigrav[i-1]+mm2/2, star.Psi0[i-1]+nn2/2, om);
                mm3=dr*rhsBSint_phigrav(star.r[i-1]+dr/2, star.X[i-1]+ll2/2, star.A[i-1]+qq2/2, star.eta[i-1]+ss2/2, star.Phi[i-1]+kk2/2, star.phigrav[i-1]+mm2/2, star.Psi0[i-1]+nn2/2, om);
                qq3=dr*rhsBSint_A(star.r[i-1]+dr/2, star.X[i-1]+ll2/2, star.A[i-1]+qq2/2, star.eta[i-1]+ss2/2, star.Phi[i-1]+kk2/2, star.phigrav[i-1]+mm2/2, star.Psi0[i-1]+nn2/2, om);

                kk4=dr*rhsBSint_phi(star.r[i-1]+dr, star.X[i-1]+ll3, star.A[i-1]+qq3, star.eta[i-1]+ss3, star.Phi[i-1]+kk3, star.phigrav[i-1]+mm3, star.Psi0[i-1]+nn3, om);
                ll4=dr*rhsBSint_X(star.r[i-1]+dr, star.X[i-1]+ll3, star.A[i-1]+qq3, star.eta[i-1]+ss3, star.Phi[i-1]+kk3, star.phigrav[i-1]+mm3, star.Psi0[i-1]+nn3, om);
                nn4=dr*rhsBSint_psi0(star.r[i-1]+dr, star.X[i-1]+ll3, star.A[i-1]+qq3, star.eta[i-1]+ss3, star.Phi[i-1]+kk3, star.phigrav[i-1]+mm3, star.Psi0[i-1]+nn3, om);
                ss4=dr*rhsBSint_eta(star.r[i-1]+dr, star.X[i-1]+ll3, star.A[i-1]+qq3, star.eta[i-1]+ss3, star.Phi[i-1]+kk3, star.phigrav[i-1]+mm3, star.Psi0[i-1]+nn3, om);
                mm4=dr*rhsBSint_phigrav(star.r[i-1]+dr, star.X[i-1]+ll3, star.A[i-1]+qq3, star.eta[i-1]+ss3, star.Phi[i-1]+kk3, star.phigrav[i-1]+mm3, star.Psi0[i-1]+nn3, om);
                qq4=dr*rhsBSint_A(star.r[i-1]+dr, star.X[i-1]+ll3, star.A[i-1]+qq3, star.eta[i-1]+ss3, star.Phi[i-1]+kk3, star.phigrav[i-1]+mm3, star.Psi0[i-1]+nn3, om);

                star.Phi[i]   = star.Phi[i-1]  + (kk1+2*kk2+2*kk3+kk4)/6;
                star.X[i]    = star.X[i-1]   + (ll1+2*ll2+2*ll3+ll4)/6;
                star.Psi0[i] = star.Psi0[i-1]+ (nn1+2*nn2+2*nn3+nn4)/6;
                star.eta[i]  = star.eta[i-1] + (ss1+2*ss2+2*ss3+ss4)/6;
                star.phigrav[i]  = star.phigrav[i-1] + (mm1+2*mm2+2*mm3+mm4)/6;
                star.A[i]  = star.A[i-1] + (qq1+2*qq2+2*qq3+qq4)/6;
               
        }
  }
/*==========================================================================*/

void rhsBSint(double* rhs_X, double* rhs_A, double* rhs_eta, double* rhs_Phi, double* rhs_phigrav, double* rhs_Psi0,
              double r, double X, double A, double eta, double Phi, double phigrav, double Psi0, double om)
  {
  if(r < 1.0e-15)
    {
    // My version 
    // *rhs_X   = 0;
    // *rhs_phigrav   = 0;
    // *rhs_Psi0   = (-om*om * A * exp(-2*Phi) + (1/F(phigrav)) * Vp(A) * A) / 3;
    // *rhs_A   = 0;
    // *rhs_eta = (sqrt(F(phigrav)) * derW(phigrav) + 2*PI * derF(phigrav) / (F(phigrav) * sqrt(F(phigrav))) * (om*om * A*A * exp(-2*Phi) * F(phigrav) - 2*V(A))) / 3;
    // *rhs_Phi = 0;

    //Roxana's version
    *rhs_X   = 0;
    *rhs_phigrav   = 0;
    *rhs_Psi0   = (-om*om * A / (alpha(Phi, phigrav)*alpha(Phi, phigrav)) +  Vp(A) * A) / (3*F(phigrav));
    *rhs_A   = 0;
    *rhs_eta = (sqrt(F(phigrav)) * derW(phigrav) + 2*PI * derF(phigrav) / (F(phigrav) * sqrt(F(phigrav))) * (om*om * A*A / (alpha(Phi, phigrav) * alpha(Phi, phigrav)) - 2*V(A))) / 3;
    *rhs_Phi = 0;
    }
  else
    {
    
    //My version
    // *rhs_Phi = 0.5 * (F(phigrav)*X*X - 1) / r - r * F(phigrav) * (X*X) * W(phigrav) + (r/2) * (X*X) * (eta*eta) + 2*PI * r * X*X * (1/F(phigrav)) * (Psi0*Psi0 / (X*X)
    //            + om*om * exp(-2*Phi) * A*A * F(phigrav) - V(A));
    // *rhs_X   = (r/2) * (X*X*X) * (eta*eta) + r * F(phigrav) * (X*X*X) * W(phigrav) - 0.5 * X * (F(phigrav)*X*X - 1) / r - 0.5 * derF(phigrav) * (X*X) * eta + 2*PI * r * X*X*X / F(phigrav)
    //            * (Psi0*Psi0 / (X*X) + om*om * exp(-2*Phi) * A*A * F(phigrav) + V(A));
    // *rhs_phigrav   = X * eta;
    // *rhs_A = Psi0;
    // *rhs_eta = - eta * ((*rhs_Phi) - 0.5 * derF(phigrav) * X * eta) - 2 * eta / r + F(phigrav) * X * derW(phigrav) + 
    //             + 2*PI * X * derF(phigrav) * (1/F(phigrav)) * (om*om * exp(-2*Phi) * A*A * F(phigrav) - Psi0*Psi0 / (X*X) - 2*V(A));
    // *rhs_Psi0 = -2 * Psi0 * (1/r) + Psi0 * ((*rhs_X)/X + 1.5 * derF(phigrav) * X * eta - (*rhs_Phi)) - (X*X) * (om*om) * A * exp(-2*Phi) * F(phigrav) + (X*X) * A * Vp(A);


    //Roxana's version
    *rhs_Phi = (F(phigrav)*pow(X,2)-1)/(2*r)+r*pow(X,2)*pow(eta,2)/2+2*PI*r*pow(X,2)*(pow(Psi0/X,2)+pow(om*A/alpha(Phi, phigrav),2)-V(A))/F(phigrav)-r*F(phigrav)*pow(X,2)*W(phigrav);

    *rhs_X = 2*PI*r*X*X*X*(Psi0*Psi0/(X*X)+(om*om)*(A*A)/(alpha(Phi, phigrav)*alpha(Phi, phigrav))+V(A))/F(phigrav)
            +r*X*X*X*eta*eta/2+r*F(phigrav)*X*X*X*W(phigrav)+(X-F(phigrav)*X*X*X)/(2*r)-derF(phigrav)*X*X*eta/2;
  
    *rhs_phigrav   = X * eta;
    *rhs_A = Psi0;

    *rhs_Psi0 = -om*om*X*X*A/(alpha(Phi, phigrav)*alpha(Phi, phigrav))+eta*X*Psi0*derF(phigrav)+X*X*A*Vp(A)-Psi0*(1+F(phigrav)*X*X)/r
                +4*PI*r*X*X*Psi0*V(A)/F(phigrav)+2*F(phigrav)*r*X*X*Psi0*W(phigrav);

    *rhs_eta = X*F(phigrav)*derW(phigrav)-3*eta/(2*r)-eta*F(phigrav)*X*X/(2*r)-r*eta*eta*eta*X*X/2+eta*eta*X*derF(phigrav)/2
              +F(phigrav)*r*eta*X*X*W(phigrav)+2*PI*X*(derF(phigrav)/F(phigrav))*(-Psi0*Psi0/(X*X)
              +om*om*A*A/(alpha(Phi, phigrav)*alpha(Phi, phigrav))-2*V(A))-2*PI*r*eta*X*X*(Psi0*Psi0/(X*X)
              +om*om*A*A/(alpha(Phi, phigrav)*alpha(Phi, phigrav))-V(A))/F(phigrav);
    }
  }

/*==========================================================================*/

double rhsBSint_X(double r, double X, double A, double eta, double Phi, double phigrav, double Psi0, double om)
  {
  double result;
  if(r < 1.0e-15)
    {
    result  = 0;
    }
  else
    {
    result = 2*PI*r*X*X*X*(Psi0*Psi0/(X*X)+(om*om)*(A*A)/(alpha(Phi, phigrav)*alpha(Phi, phigrav))+V(A))/F(phigrav)
            +r*X*X*X*eta*eta/2+r*F(phigrav)*X*X*X*W(phigrav)+(X-F(phigrav)*X*X*X)/(2*r)-derF(phigrav)*X*X*eta/2;
    }
  return result;
  }

/*==========================================================================*/

double rhsBSint_phigrav(double r, double X, double A, double eta, double Phi, double phigrav, double Psi0, double om)
  {
    double result;
  if(r < 1.0e-15)
    {
    result   = 0;
    }
  else
    {
    result   = X * eta;
    }
  return result;
  }

/*==========================================================================*/

double rhsBSint_phi(double r, double X, double A, double eta, double Phi, double phigrav, double Psi0, double om)
  {
  double result;
  if(r < 1.0e-15)
    {
    result = 0;
    }
  else
    {
    result = (F(phigrav)*pow(X,2)-1)/(2*r)+r*pow(X,2)*pow(eta,2)/2+2*PI*r*pow(X,2)*(pow(Psi0/X,2)+pow(om*A/alpha(Phi, phigrav),2)-V(A))/F(phigrav)-r*F(phigrav)*pow(X,2)*W(phigrav);
    }
  return result;
  }

/*==========================================================================*/

double rhsBSint_psi0(double r, double X, double A, double eta, double Phi, double phigrav, double Psi0, double om)
  {
  double result;
  if(r < 1.0e-15)
    {
    result   = (-om*om * A / (alpha(Phi, phigrav)*alpha(Phi, phigrav)) +  Vp(A) * A) / (3*F(phigrav));
    }
  else
    {
    result = -om*om*X*X*A/(alpha(Phi, phigrav)*alpha(Phi, phigrav))+eta*X*Psi0*derF(phigrav)+X*X*A*Vp(A)-Psi0*(1+F(phigrav)*X*X)/r
                +4*PI*r*X*X*Psi0*V(A)/F(phigrav)+2*F(phigrav)*r*X*X*Psi0*W(phigrav);
    }
  return result;
  }

/*==========================================================================*/

double rhsBSint_A(double r, double X, double A, double eta, double Phi, double phigrav, double Psi0, double om)
  {
  double result;
  if(r < 1.0e-15)
    {
    result   = 0;
    }
  else
    {
    result = Psi0;
    }
  return result;
  }

/*==========================================================================*/

double rhsBSint_eta(double r, double X, double A, double eta, double Phi, double phigrav, double Psi0, double om)
  {
  double result;
  if(r < 1.0e-15)
    {
    result = (sqrt(F(phigrav)) * derW(phigrav) + 2*PI * derF(phigrav) / (F(phigrav) * sqrt(F(phigrav))) * (om*om * A*A / (alpha(Phi, phigrav) * alpha(Phi, phigrav)) - 2*V(A))) / 3;
    }
  else
    {
    result = X*F(phigrav)*derW(phigrav)-3*eta/(2*r)-eta*F(phigrav)*X*X/(2*r)-r*eta*eta*eta*X*X/2+eta*eta*X*derF(phigrav)/2
              +F(phigrav)*r*eta*X*X*W(phigrav)+2*PI*X*(derF(phigrav)/F(phigrav))*(-Psi0*Psi0/(X*X)
              +om*om*A*A/(alpha(Phi, phigrav)*alpha(Phi, phigrav))-2*V(A))-2*PI*r*eta*X*X*(Psi0*Psi0/(X*X)
              +om*om*A*A/(alpha(Phi, phigrav)*alpha(Phi, phigrav))-V(A))/F(phigrav);
    }
  return result;
  }

/*==========================================================================*/
double V_series(double A)
  {
  // Potential function
  return A*A * (1 + par.lambda4 * A*A + par.lambda6 * A*A*A*A
         + par.lambda8 * A*A*A*A*A*A);
  }

/*==========================================================================*/

double Vp_series(double A)
  {
  // Potential derviative dV / d(A^2)
  return 1 + 2*par.lambda4 * A*A + 3*par.lambda6 * A*A*A*A
         + 4*par.lambda8 * A*A*A*A*A*A;
  }

/*==========================================================================*/

double V_solitonic(double A)
  {
  // Solitonic potential function
  return A*A * (1 - 2 * A*A / (par.sigma0*par.sigma0))
             * (1 - 2 * A*A / (par.sigma0*par.sigma0));
  }

/*==========================================================================*/

double Vp_solitonic(double A)
  {
  // Solitonic potential function
  return (1 - 2 * A*A / (par.sigma0*par.sigma0))
       * (1 - 6 * A*A / (par.sigma0*par.sigma0));
  }

/*==========================================================================*/

double F_st(double phigrav)
  {
  // Coupling function of the gravitational scalar field
  return exp(-2 * par.alpha0 * phigrav - par.beta0 * phigrav * phigrav);
  }

/*==========================================================================*/

double derF_st(double phigrav)
  {
  // Function F_{,\varphi} / F
  return -2 * par.alpha0 - 2 * par.beta0 * phigrav;
  }

/*==========================================================================*/

double W_st(double phigrav)
  {
  // Potential for the gravitational scalar field
  return 0.5 * par.mphigrav0 * par.mphigrav0 * phigrav * phigrav;
  }

/*==========================================================================*/

double derW_st(double phigrav)
  {
  // Derivative of the potential for the gravitational scalar field
  return par.mphigrav0 * par.mphigrav0 * phigrav;
  }

/*==========================================================================*/

double alpha(double Phi, double phigrav)
{
  double  result;
		
	result = exp(Phi)/sqrt(F(phigrav));
	return result;
}

/*==========================================================================*/

void initGrid(int n, double rmax)
  {
  int i;


  for(i = 0; i < n; i++)
    star.r[i] = rmax * i / (n - 1.0);
  }

/*==========================================================================*/

void readPars(char* ifil)
  {
  FILE* ifp;
  char  line[LEN];
  int   n, i;


  // First set all parameters to default values
  par.nint      = 401;
  par.rmax      = 1;
  par.A0        = 0.07;
  par.mphigrav0 = 1.0;
  par.alpha0    = 3.0;
  par.beta0     = 0.0;
  par.phigrav0  = 0.0;
  par.Phi0      = -4.3723594315020059e-01;
  par.omega0    = 1.;            // omega0 is always 1, it is not specified
  par.nzerotarget = 0;
  par.thresh    = 2e-16;
  par.thresh_grav    = 2e-6;
  par.mpercentage = 90;
  par.rmatchfac = 1;
  strcpy(par.potential, "series");
  par.lambda4   = 0;
  par.lambda6   = 0;
  par.lambda8   = 0;
  par.sigma0    = 1;
  strcpy(par.minmax, "min");


  // Read parameters from file and overwrite default values
  ifp = fopen(ifil, "r");
  if(ifp == NULL) { printf("Cannot open %s in readPars\n\n", ifil); exit(0); }

  while(mygetline(line, ifp) != 0)
    {
    if(line[0] == '#')
      ;   // Comment line
    else
      {
      n = strlen(line);
      // Replace '=' with white space for syntax flexibility
      for(i = 0; i < n; i++)
        if(line[i] == '=') line[i] = ' ';

      // Analyze parameter values
      if(strstr(line, "nint") != NULL)
        sscanf(line, "nint %d", &(par.nint));
      else if(strstr(line, "rmax") != NULL)
        sscanf(line, "rmax %le", &(par.rmax));
      else if(strstr(line, "A0") != NULL)
        sscanf(line, "A0 %le", &(par.A0));
      else if(strstr(line, "mphigrav0") != NULL)
        sscanf(line, "mphigrav0 %le", &(par.mphigrav0));
      else if(strstr(line, "alpha0") != NULL)
        sscanf(line, "alpha0 %le", &(par.alpha0));
      else if(strstr(line, "beta0") != NULL)
        sscanf(line, "beta0 %le", &(par.beta0));
      else if(strstr(line, "phigrav0") != NULL)
        sscanf(line, "phigrav0 %le", &(par.phigrav0));
      else if(strstr(line, "Phi0") != NULL)
        sscanf(line, "Phi0 %le", &(par.Phi0));
      else if(strstr(line, "omega0") != NULL)
        sscanf(line, "omega0 %le", &(par.omega0));
      else if(strstr(line, "nzerotarget") != NULL)
        sscanf(line, "nzerotarget %d", &(par.nzerotarget));
      else if(strstr(line, "thresh") != NULL)
        sscanf(line, "thresh %le", &(par.thresh));
      else if(strstr(line, "mpercentage") != NULL)
        sscanf(line, "mpercentage %le", &(par.mpercentage));
      else if(strstr(line, "rmatchfac") != NULL)
        sscanf(line, "rmatchfac %le", &(par.rmatchfac));
      else if(strstr(line, "potential") != NULL)
        sscanf(line, "potential %s", par.potential);
      else if(strstr(line, "minmax") != NULL)
        sscanf(line, "minmax %s", par.minmax);
      else if(strstr(line, "lambda4") != NULL)
        sscanf(line, "lambda4 %le", &(par.lambda4));
      else if(strstr(line, "lambda6") != NULL)
        sscanf(line, "lambda6 %le", &(par.lambda6));
      else if(strstr(line, "lambda8") != NULL)
        sscanf(line, "lambda8 %le", &(par.lambda8));
      else if(strstr(line, "sigma0") != NULL)
        sscanf(line, "sigma0 %le", &(par.sigma0));
      }
    }

  fclose(ifp);
  }

/*==========================================================================*/

void printPars()
  {
  printf("=======================================\n");
  printf("nint          = %d\n", par.nint);
  printf("rmax          = %g\n", par.rmax);
  printf("A0            = %g\n", par.A0);
  printf("mphigrav0     = %g\n", par.mphigrav0);
  printf("phigrav0      = %g\n", par.phigrav0);
  printf("omega0      = %g\n", par.omega0);
  printf("alpha0        = %g\n", par.alpha0);
  printf("beta0         = %g\n", par.beta0);
  printf("nzerotarget   = %d\n", par.nzerotarget);
  printf("thresh        = %g\n", par.thresh);
  printf("mpercentage   = %g\n", par.mpercentage);
  printf("rmatchfac     = %g\n", par.rmatchfac);
  printf("potential     = %s\n", par.potential);
  printf("minmax        = %s\n", par.minmax);
  if(strcmp(par.potential, "series") == 0)
    {
    printf("lambda4       = %g\n", par.lambda4);
    printf("lambda6       = %g\n", par.lambda6);
    printf("lambda8       = %g\n", par.lambda8);
    }
  else if(strcmp(par.potential, "solitonic") == 0)
    printf("sigma0        = %g\n", par.sigma0);
  printf("=======================================\n");
  }

/*==========================================================================*/

void registerPotential()
  {
  if(strcmp(par.potential, "series") == 0)
    {
    V  = (double (*)(double))(V_series);
    Vp = (double (*)(double))(Vp_series);
    W = (double (*)(double))(W_st);
    derW = (double (*)(double))(derW_st);
    F = (double (*)(double))(F_st);
    derF = (double (*)(double))(derF_st);
    }
  else if(strcmp(par.potential, "solitonic") == 0)
    {
    V  = (double (*)(double))(V_solitonic);
    Vp = (double (*)(double))(Vp_solitonic);
    W = (double (*)(double))(W_st);
    derW = (double (*)(double))(derW_st);
    F = (double (*)(double))(F_st);
    derF = (double (*)(double))(derF_st);
    }
  else
    { printf("Unknown potential:   %s\n\n", par.potential); exit(0); }
  }

/*==========================================================================*/

int mygetline( char *s, FILE *ifp )
  {
  int c, i;

  for( i = 0; i < LEN && ( (c = getc(ifp)) != EOF ) && (c != '\n'); i++ )
    s[i] = c;

  if( c == '\n' ) s[i++] = c;

  s[i] = '\0';

  return i;
  }

/*==========================================================================*/

void rescalePhi()
  {
  int i, n1;
  double Phitarget, Phiorg;


  // We want the lapse = 1 at infinity. Since we do not go to infinity,
  // our condition on the lapse is F * alpha * X = 1. This fixes Phi
  // at our outer radius to
  //
  //   Phi = -ln(\sqrt{F} X)
  //
  // Note that we also change our time coordinate that way and need
  // to adjust omega such that
  //
  //   omorg / alphaorg = omnew / alphanew
  //
  // So omnew = omorg * exp(Phinew - Phiorg)
  n1 = par.nint;
  Phitarget = -log(sqrt(F(star.phigrav[n1-1])) * star.X[n1-1]);
  Phiorg    = star.Phi[n1-1];
  printf("Phiorg = %g     Phitarget = %g\n", Phiorg, Phitarget);

  for(i = 0; i < n1; i++)
    star.Phi[i] += (Phitarget - Phiorg);

  par.omega0 *= exp(Phitarget - Phiorg);
  }

/*==========================================================================*/

void calcMass()
  {
  int i;
  double r, X, phigrav;


  for(i = 0; i < par.nint; i++)
    {
    r = star.r[i];
    X = star.X[i];
    phigrav = star.phigrav[i];
    star.m[i] = 0.5 * r * (1 - 1 / (F(phigrav)*X*X));
    if(i == 0)
      star.C[i] = 0;
    else
      star.C[i] = star.m[i] / star.r[i];
    }
  }

/*==========================================================================*/

double calcRadius()
  {
  int n1, i;


  n1 = par.nint;
  for(i = n1-2; i >= 0; i--)
    if(star.m[i] < par.mpercentage / 100.0 * star.m[n1-1])
      break;

  printf("total mass = %22.16g\n", star.m[n1-1]);
  printf("star.m[%d] = %22.16g   star.m[%d] = %22.16g\n",
         i, star.m[i], i+1, star.m[i+1]);
  printf("radius = %22.16g\n", star.r[i+1]);

  return star.r[i+1];
  }

/*==========================================================================*/

int calcRadius_index()
  {
  int n1, i;


  n1 = par.nint;
  for(i = n1-2; i >= 0; i--)
    if(star.m[i] < par.mpercentage / 100.0 * star.m[n1-1])
      break;

  return i+1;
  }

/*==========================================================================*/

void calcIso()
  {
  int i, n1;
  double dr, r, R, phigrav, X, f, rhs_f, df[5];
  double m, Rfac;


  n1 = par.nint;

  // Central values;
  star.f[0]   = 1;
  star.R[0]   = 0;
  star.psi[0] = 1;

  // RK integration
  for(i = 1; i < n1; i++)
    {
    dr   = star.r[i] - star.r[i-1];

    // 1st RK step
    r = star.r[i-1];
    X = star.X[i-1];
    f = star.f[i-1];

    rhsIso(&rhs_f, r, phigrav, X, f);
    df[1] = rhs_f * dr;

    // 2nd RK step
    r = 0.5 * (star.r[i-1] + star.r[i]);
    X = 0.5 * (star.X[i-1] + star.X[i]);
    f = star.f[i-1] + 0.5 * df[1];
    rhsIso(&rhs_f, r, phigrav, X, f);
    df[2] = rhs_f * dr;

    // 3rd RK step
    r = 0.5 * (star.r[i-1] + star.r[i]);
    X = 0.5 * (star.X[i-1] + star.X[i]);
    f = star.f[i-1] + 0.5 * df[2];
    rhsIso(&rhs_f, r, phigrav, X, f);
    df[3] = rhs_f * dr;

    // 4th RK step
    r = star.r[i];
    X = star.X[i];
    f = star.f[i-1] + df[3];
    rhsIso(&rhs_f, r, phigrav, X, f);
    df[4] = rhs_f * dr;

    // Update variable and compute new variables
    star.f[i]   = star.f[i-1] + (df[1] + 2*df[2] + 2*df[3] + df[4]) / 6.0;
    star.R[i]   = star.f[i] * star.r[i];
    star.psi[i] = 1 / sqrt(star.f[i] * sqrt(F(star.phigrav[i])));
    }


  // So far our solution is only determined up to an overall constant
  // factor, i.e. c*R is also a solution. We have fixed that by
  // requiring that at r->0, R=r. But ultimately, we would rather
  // obtain identical radii at infinity, i.e. lim_r->\infty (R-r) = 0.
  // This is best realized by recalling that r is the areal radius,
  // i.e. the surface area A(r)=4\pi r^2. From the isotropic line
  // element we likewise find A(R)=4\pi \psi^4 R^2, where \psi is
  // the conformal factor. At large radius, the boson star is
  // virtually vacuum and we can approximate the conformal factor by
  // the Brill-Lindquist value \psi=1-m/(2R). Inserting this, we obtain
  //
  //  Rtarget = (r/sqrt{F}-m) - 1/4 (sqrt{F} m^2)/(r - \sqrt{F}m), which reduced to the usual form in the limit F --> 1 for boson stars.
  //
  // And then rescale over the whole grid
  //
  //   R -> R * Rtarget / R[n-1]

  #if 1
  m = star.m[n1-1];
  r = star.r[n1-1];
  printf("matching radii at   r[%d] = %g\n", n1-1, star.r[n1-1]);
  Rfac = 0.5 * (r / sqrt(F(star.phigrav[n1-1])) - m) * ( 1 + sqrt(1 - F(star.phigrav[n1-1])*m*m / (r - sqrt(F(star.phigrav[n1-1]))*m) / (r - sqrt(F(star.phigrav[n1-1]))*m)) );
  Rfac /= star.R[n1-1];

  for(i = 0; i < n1; i++)
    {
    star.R[i] *= Rfac;
    star.f[i] *= Rfac;
    }
  #endif
  }

/*==========================================================================*/

void rhsIso(double* rhs_f, double r, double phigrav, double X, double f)
  {
  if(r < 1.0e-15)
    {
    // We are at r = 0 and need to use the asymptotic behvaviour.
    // We have f' = f/r * (1 - X) and know X ~ 1 + O(r^2), so that
    // near r = 0, f' = 0.
    *rhs_f = 0;
    }
  else
    {
    *rhs_f = f * (X*sqrt(F(phigrav)) - 1) / r;
    }
  }

/*==========================================================================*/

void out1D(double* x, double* y, int n1, int n2, const char* ofil)
  {
  int i;
  FILE* ofp;


  ofp = fopen(ofil, "w");
  if(ofp == NULL) { printf("Cannot open %s in out1D\n\n", ofil); exit(0); }

  fprintf(ofp, "# %s\n", ofil);

  for(i = n1; i <= n2; i++)
    fprintf(ofp, "%22.16g   %22.16g\n", x[i], y[i]);

  fclose(ofp);
  }

/*==========================================================================*/

double findMax(double* f, int n1, int n2)
  {
  int    i;
  double val;


  val = f[n1];
  for(i = n1+1; i <= n2; i++)
    if(f[i] > val) val = f[i];

  return val;
  }

/*==========================================================================*/
