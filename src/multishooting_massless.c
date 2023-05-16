/*
Complete routine for generating MofR plot for the massless case of BSs in ST theory.
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
  double        *C;     // Compactness defined as m(r)/r
  double        *R;     // Isotropic radius
  double        *f;     // f := R / r
  double        *psi;   // conformal factor, i.e. g_rr = psi^4
  double        *Psi0;
  double        *phigrav;
  } model;

typedef struct param
  {
  double        rmax;
  double        mpercentage;   // mass % that determines r_BS
  double        rmatchfac;     // Rescale matching radius
  int           nint;
  double        omega0;
  double        mphigrav0;
  double        A0;
  double        alpha0;
  double        beta0;
  double        phigrav0;
  double        Amin;
  double        Amax;
  int           nzerotarget;
  int           verbose;
  double        thresh;
  int           nmodels;
  double        lambda4;
  double        lambda6;
  double        lambda8;
  double        sigma0;
  char          potential[LEN];
  } param;


// Functions
void    append2Data     (const char*, double, double);
void    append7Data     (const char*, double, double, double, double, double, double, double);
int    calcExtAsym      ();
void    calcIso         ();
void    calcMass        ();
double  calcRadius      ();
int  calcRadius_index      ();
int     findFreqMinMax  (double*, double*);
void    initGrid        (int, double);
double  findMax         (double*, int, int);
void    intODE          (double, double, double, int*, double*, int*);
int     iofr            (double);
int     mygetline       (char*, FILE*);
int    onemodel        (double*, double*, double*, double*);
void    openFile        (const char*);
void    out1D           (double*, double*, int, int, const char*);
void    printPars       ();
void    readPars        (char*);
void    registerPotential();
void    rescalePhi      ();
void    rhsBSint        (double*, double*, double*, double*, double*, double*, double, double, double, double, double, double, double, double);
void    rhsIso          (double*, double, double, double, double);
int    shoot           ();
double  V_series        (double);
double  Vp_series       (double);
double  V_solitonic     (double);
double  Vp_solitonic    (double);
double  F_st               (double);
double  derF_st            (double);
double  W_st               (double);
double  derW_st            (double);
int check_phigrav      (double);

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
  int i;
  int success;
  double omBS, mBS, rBS, CBS, dA;


  if(argc != 2) { printf("Usage:   %s   <parfile>\n\n", argv[0]); exit(0); }

  // Parameters
  readPars(argv[1]);
  printPars();

  double phigrav_init = par.phigrav0;

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
  star.C   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.R   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.f   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.psi = (double*) malloc((size_t) par.nint * sizeof(double));

  openFile("MofR.dat");
  openFile("MofA.dat");
  openFile("Mofphigrav.dat");
  openFile("AuxiliaryInfo.dat");
  printf("\n\n   A0           phigrav0           omega               mass            radius      max(C)\n");
  printf("==============================================================================================================\n");

  if(par.nmodels == 1) dA = 0;
  else dA = (par.Amax - par.Amin) / (par.nmodels - 1);
  for(i = 0; i < par.nmodels; i++)
    {
    par.phigrav0 = phigrav_init;
    par.omega0 = 1;   // par.omega0 was overwritten with the last solution.
                      // Now we want to reinitialize it for the next model.
    
    par.A0 = par.Amin + dA * i;
    success = onemodel(&omBS, &mBS, &rBS, &CBS);
    if (success==0)
    {
      printf("No solution found for %11.6g\n", par.A0);
    }
    else
    {
    printf("%11.6g  %11.6g   %16.10g   %16.10g   %12.6g   %12.6g\n",
           par.A0, par.phigrav0, omBS, mBS, rBS, CBS);
    append2Data("MofR.dat", rBS, mBS);
    append2Data("MofA.dat", par.A0, mBS);
    append2Data("Mofphigrav.dat", par.phigrav0, mBS);
    double alpha0 = exp(star.Phi[0])/sqrt(F(star.phigrav[0]));
    append7Data("AuxiliaryInfo.dat", omBS, alpha0, par.phigrav0, par.A0, mBS, rBS, CBS);
    }
    }


  printf("===============================================================================\n");
  }

/*==========================================================================*/

int onemodel(double* omBS, double* mBS, double* rBS, double* CBS)
  {
  int k, i;
  int converged, success, check_solution;

  for (k = 0; k < 51; ++k)
  { 
    if (par.verbose)
    {
    printf("I am starting iteration number: %d\n", k);
    printf("Gravitational scalar field guess is set to %22.16g\n", par.phigrav0);
    }

    // Shoot
    success = shoot();
    if (k == 49 && success == 0)
    { 
      printf("Shooting failed for this central boson amplitude, I could not find omin and omax values!\n");
      return 0;
     }
    else
    {
    // Calculate exterior
    converged = calcExtAsym();

    if (converged == 0)
    {par.phigrav0 -= 0.0001;}

    // Rescale the time-time component of the metric such that it
    // corresponds to a Schwarzschild exterior.
    rescalePhi();

    //This is where we define the surface of the boson star to be -- this should be far enough!
    int jj = par.nint - 100;

    // Surface criterion of https://iopscience.iop.org/article/10.1088/0264-9381/33/13/135002/pdf Eq.(3.15).
    double derivative_term = (star.Phi[jj] - star.Phi[jj-1])/(star.r[jj] - star.r[jj-1]);
    double sqrt_expr = sqrt(derivative_term * derivative_term + star.X[jj] * star.X[jj] * star.eta[jj] * star.eta[jj]);
    double varphi_surface = 0.0 - (star.X[jj] * star.eta[jj]) / sqrt_expr * atanh( sqrt_expr / (derivative_term + 1 / star.r[jj]));
    double criterion = star.phigrav[jj] - varphi_surface;
    
    if (par.verbose)
    {
    printf("Surface boundary condition:   %g\n", varphi_surface);
    printf("Gravitational scalar field at the surface boundary condition:   %g\n", star.phigrav[jj]);
    }

    if (converged == 0 && k == 49)
    { 
      printf("I've reached max number of iterations searching for a solution and have not converged\n");
      return 0;
    }

    if (fabs(criterion) < 1e-06 && converged == 1)
      {
        // for(i = 0; i < par.nint; i++)
        // {
        //   check_solution = check_phigrav(star.phigrav[i]);
        //   if (check_solution == 0)
        //   { 
        //     printf("Rejecting solution since it crosses -alpha0/beta0 line!\n");
        //     return 0;
        //   }
        // }
        // Compute diagnostic quantities and transform to isotropic coordinates.
        calcMass();

        *rBS  = calcRadius();
        *mBS  = star.m[par.nint-1];
        *omBS = par.omega0;
        *CBS  = findMax(star.C, 0, par.nint-1); 
        calcIso();

        if (par.verbose)
        {printf("I get the following difference, %22.16g\n", fabs(criterion));}
        break;
      }
      // If surface criterion not satisfied to the above tolerance, use Newton-Raphson to improve the guess for \varphi_c.
      else {
        if (par.verbose)
        {printf("I get the following difference, %22.16g\n", fabs(criterion));}
        par.phigrav0 += 1e-05;

        // Calculate the model all over again
        shoot();
        converged = calcExtAsym();

        // if (converged == 0) 
        // {
        //   printf("Did not find a solution for this central boson amplitude!\n");
        //   break;
        // }

        // else
        // {
          rescalePhi();

          double sqrt_expr_v2 = sqrt(derivative_term * derivative_term + star.X[jj] * star.X[jj] * star.eta[jj] * star.eta[jj]);
          double varphi_surface_v2 = 0.0 - (star.X[jj] * star.eta[jj]) / sqrt_expr_v2
                            * atanh( sqrt_expr_v2 / (derivative_term + 1 / star.r[jj]));
          double criterion_v2 = star.phigrav[jj] - varphi_surface_v2;

          double criterion_derivative = (criterion_v2 - criterion) / (1e-05);
      
          par.phigrav0 -= (1e-05 + criterion/criterion_derivative);
            }
        }
  }

  if(par.verbose)
    {
    printf("\n===================================================\n");
    printf("Physical frequency:   omega = %22.16g\n", par.omega0);
    printf("Total mass:           m     = %22.16g\n", star.m[par.nint-1]);
    printf("===================================================\n");
    }

  return 1;
  }

/*==========================================================================*/

int shoot()
  {
  int    nzero, oldnzero;
  double omcur;
  int    cnt, success;
  double ommax, ommin, om;
  int    sigAmax, sigAmin, sigA;
  double rstopmin, rstopmax, rstop;


  if(par.verbose)
    {
    printf("\n");
    printf("omega               Zero crossings            Rstop        sigA\n");
    printf("===============================================================\n");
    }


  // Grid setup
  initGrid(par.nint, par.rmax);

  // (1) Find an upper and a lower limit for the frequency. This is done
  //     in findFreqMinMax, starting with omega=1 and then doubling or
  //     halving the value; see comments in that function. Should this
  //     search fail, we quit.
  success = findFreqMinMax(&ommin, &ommax);
  if(success == 0)
    { 
      if (par.verbose)
      {printf("Could not find ommin and ommax in shoot\n");}
      return 0;
    }
  else
  {
  if(par.verbose)
    printf("Using   omega in   [%22.16g,   %22.16g]\n", ommin, ommax);


  // (2) Now we finetune the frequency such that it sits right on the
  //     threshold between nzerotarget and nzerotarget+1. This is the
  //     model we seek with nzerotarget zero crossings.
  if(par.verbose)
    {
    printf("\n");
    printf("omega               Zero crossings            Rstop        sigA\n");
    printf("===============================================================\n");
    }

  cnt = 0;
  while( (ommax - ommin) > par.thresh )
    {
    om = 0.5 * (ommin + ommax);
    intODE(par.A0, par.phigrav0, om, &nzero, &rstop, &sigA);
    if(par.verbose)
      printf("%20.15g          %d          %15.7g           %d\n",
             om, nzero, rstop, sigA);
    if( nzero > par.nzerotarget)
      ommax = om;  // the test frequency is still an upper limit
    else
      ommin = om;

    // If we fail to reach the threshold over many iterations, we likely
    // hit round-off error and reduce the accuracy requirement.
    cnt++;
    if(cnt > 100)
      {
      par.thresh *= 10;
      cnt = 0;
      printf("Changed threshold to   %g\n", par.thresh);
      }
    }

  // ommin gives our best estimate for the model, since it has
  // nzero zero crossings whereas ommax has one zero crossing more.
  // We store this frequency in the parameter omega0 for further use.
  par.omega0 = ommin;
  return 1;
  }
  }

/*==========================================================================*/

int  findFreqMinMax(double* ommin, double* ommax)
  {
  double rstop;
  int sigA, nzero, cnt, cntmax;


  // Our strategy is as follows. We use the fact that the number of
  // zero crossings is a non-decreasing function of the frequency;
  // increasing omega0 gives you at least as many zero crossings as
  // before. We have a target, ntarget. We compute the number of zero
  // crossings for omega0 = 1.
  // If the resulting n>ntarget, we half omega0 until n <= ntarget
  // (note that the stable BS is always the limiting case to have
  // one more zero crossing -- the new crossing is at infinity in
  // the limiting case). Once we have n <= ntarget, the corresponding
  // two omega0 values are the brackets.
  // If the resulting n<=ntarget, we double omega0 until we have
  // n > ntarget and the ensuing omega0 values are our brackets.

  intODE(par.A0, par.phigrav0, par.omega0, &nzero, &rstop, &sigA);
  if(par.verbose)
    printf("%15.10g          %d          %15.7g           %d\n",
           par.omega0, nzero, rstop, sigA);

  cnt    = 0;           // cnt is a sanity measure to avoid hanging
  cntmax = 100;         // in the loop forever. It should never be
                        // triggered unless the frequency has values
                        // of 2^{1000} or 2^{-100}...

  if(nzero > par.nzerotarget)
    {
    *ommin = par.omega0;
    while(nzero > par.nzerotarget && cnt < cntmax)
      {
      *ommax = *ommin;
      *ommin /= 2;
      intODE(par.A0, par.phigrav0, *ommin, &nzero, &rstop, &sigA);
      if(par.verbose)
        printf("%15.10g          %d          %15.7g           %d\n",
               *ommin, nzero, rstop, sigA);
      cnt++;
      }
    }
  else
    {
    *ommax = par.omega0;
    while(nzero <= par.nzerotarget && cnt < cntmax)
      {
      *ommin = *ommax;
      *ommax *= 2;
      intODE(par.A0, par.phigrav0, *ommax, &nzero, &rstop, &sigA);
      if(par.verbose)
        printf("%15.10g          %d          %15.7g           %d\n",
               *ommax, nzero, rstop, sigA);
      cnt++;
      }
    }

  if(cnt == cntmax)
    return 0;   // Either upper or lower frequency limit has not been found.
  else
    return 1;
  }

/*==========================================================================*/

void intODE(double A0, double phigrav0, double omega, int* nzero, double* rstop, int* sigAstop)
  {
  int    i, n1, istop;
  double dr, r, X, A, eta, Phi, Psi0, phigrav, om;
  double rhs_X, rhs_A, rhs_eta, rhs_Phi, rhs_Psi0, rhs_phigrav;
  double dX[5], dA[5], dPsi0[5], deta[5], dphigrav[5], dPhi[5];


  n1 = par.nint;
  om = omega;
  *sigAstop = 0;   // Default for no stopping

  *nzero = 0;
  *rstop = star.r[n1-1];
  istop  = -1;   // So we have no vacuum region unless we exceed the amplitude


  // Central values
  star.X[0]   = 1/sqrt(F(phigrav0));
  star.A[0]   = A0;
  star.phigrav[0]   = phigrav0;
  star.Psi0[0]   = 0;
  star.eta[0] = 0;
  star.Phi[0] = 0;   // Phi has a free constant we later match to Schwarzschild

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
      istop = i-1;   // We stop the integration; one point as sanity buffer
      // There is one problem here: We regard the present data point and
      // also the previous one as unrealiable. The previous point,
      // however, might have been counted as a zero crossing. If that
      // happened, we wish to undo this counting.
      if(star.A[i-1] * star.A[i-2] < 0) (*nzero)--;
      break;
      }
    if(star.A[i] * star.A[i-1] < 0) (*nzero)++;
    }

  if(istop < 0)
    {
    *sigAstop = ( (star.A[n1-1] > 0) - (star.A[n1-1] < 0) );
    return;   // Integration went through; no need for vacuum
    }

  // Set to vacuum beyond rstop
  *rstop = star.r[istop];

  // RK integration
  for(i = istop; i < n1; i++)
    {
    dr  = star.r[i] - star.r[i-1];

    // 1st RK step
    r   = star.r[i-1];
    X   = star.X[i-1];
    A   = 0;
    Psi0 = 0;
    phigrav = star.phigrav[i-1];
    eta = star.eta[i-1];
    Phi = star.Phi[i-1];
    rhsBSint(&rhs_X, &rhs_A, &rhs_eta, &rhs_Phi, &rhs_phigrav, &rhs_Psi0, r, X, A, eta, Phi, phigrav, Psi0, om);
    dX[1]   = rhs_X * dr;
    dPhi[1] = rhs_Phi * dr;
    dphigrav[1] = rhs_phigrav * dr;
    deta[1] = rhs_eta * dr;

    // 2nd RK step
    r   = star.r[i-1] + 0.5 * dr;
    X   = star.X[i-1] + 0.5 * dX[1];
    A   = 0;
    Psi0 = 0;
    phigrav = star.phigrav[i-1] + 0.5 * dphigrav[1];
    eta = star.eta[i-1] + 0.5 * deta[1];
    Phi = star.Phi[i-1] + 0.5 * dPhi[1];
    rhsBSint(&rhs_X, &rhs_A, &rhs_eta, &rhs_Phi, &rhs_phigrav, &rhs_Psi0, r, X, A, eta, Phi, phigrav, Psi0, om);
    dX[2]   = rhs_X * dr;
    dPhi[2] = rhs_Phi * dr;
    dphigrav[2] = rhs_phigrav * dr;
    deta[2] = rhs_eta * dr;

    // 3rd RK step
    r   = star.r[i-1] + 0.5 * dr;
    X   = star.X[i-1] + 0.5 * dX[2];
    A   = 0;
    Psi0 = 0;
    phigrav   = star.phigrav[i-1] + 0.5 * dphigrav[2];
    eta = star.eta[i-1] + 0.5 * deta[2];
    Phi = star.Phi[i-1] + 0.5 * dPhi[2];
    rhsBSint(&rhs_X, &rhs_A, &rhs_eta, &rhs_Phi, &rhs_phigrav, &rhs_Psi0, r, X, A, eta, Phi, phigrav, Psi0, om);
    dX[3]   = rhs_X * dr;
    dPhi[3] = rhs_Phi * dr;
    dphigrav[3] = rhs_phigrav * dr;
    deta[3] = rhs_eta * dr;

    // 4th RK step
    r   = star.r[i];
    X   = star.X[i-1] + dX[3];
    A   = 0;
    Psi0 = 0;
    phigrav = star.phigrav[i-1] + dphigrav[3];
    eta = star.eta[i-1] + deta[3];
    Phi = star.Phi[i-1] + dPhi[3];
    rhsBSint(&rhs_X, &rhs_A, &rhs_eta, &rhs_Phi, &rhs_phigrav, &rhs_Psi0, r, X, A, eta, Phi, phigrav, Psi0, om);
    dX[4]   = rhs_X * dr;
    dPhi[4] = rhs_Phi * dr;
    dphigrav[4] = rhs_phigrav * dr;
    deta[4] = rhs_eta * dr;

    // Update variables
    star.X[i]   = star.X[i-1]   + (dX[1]  + 2*dX[2]  + 2*dX[3]  + dX[4] ) / 6.0;
    star.A[i]   = 0;
    star.Psi0[i] = 0;
    star.eta[i]   = star.eta[i-1]   + (deta[1]  + 2*deta[2]  + 2*deta[3]  + deta[4] ) / 6.0;
    star.phigrav[i] = star.phigrav[i-1]   + (dphigrav[1]  + 2*dphigrav[2]  + 2*dphigrav[3]  + dphigrav[4] ) / 6.0;
    star.Phi[i] = star.Phi[i-1] + (dPhi[1]+ 2*dPhi[2]+ 2*dPhi[3]+dPhi[4]) / 6.0;
    }

  }

/*==========================================================================*/

void rhsBSint(double* rhs_X, double* rhs_A, double* rhs_eta, double* rhs_Phi, double* rhs_phigrav, double* rhs_Psi0,
              double r, double X, double A, double eta, double Phi, double phigrav, double Psi0, double om)
  {
  if(r < 1.0e-15)
    {
    // We are at r = 0 and need to use the asymptotic behaviour.
    // For eta we use that near r=0, we have
    //
    //       eta(r) = eta(0) + eta'(0)*r + ... = eta'(0)*r + ...
    //  ==>  2*eta/r = 2*eta'(0) + ...
    //
    //  This gives a factor 1/3 to be applied to the remaining rhs
    //  of the eta equation.
    *rhs_X   = 0;
    *rhs_phigrav   = 0;
    *rhs_Psi0   = (-om*om * A * exp(-2*Phi) + (1/F(phigrav)) * Vp(A) * A) / 3;
    *rhs_A   = 0;
    *rhs_eta = (sqrt(F(phigrav)) * derW(phigrav) + 2*PI * derF(phigrav) / (F(phigrav) * sqrt(F(phigrav))) * (om*om * A*A * exp(-2*Phi) * F(phigrav) - 2*V(A))) / 3;
    *rhs_Phi = 0;
    }
  else
    {
    *rhs_Phi = 0.5 * (F(phigrav)*X*X - 1) / r - r * F(phigrav) * (X*X) * W(phigrav) + (r/2) * (X*X) * (eta*eta) + 2*PI * r * X*X * (1/F(phigrav)) * (Psi0*Psi0 / (X*X)
               + om*om * exp(-2*Phi) * A*A * F(phigrav) - V(A));
    *rhs_X   = (r/2) * (X*X*X) * (eta*eta) + r * F(phigrav) * (X*X*X) * W(phigrav) - 0.5 * X * (F(phigrav)*X*X - 1) / r - 0.5 * derF(phigrav) * (X*X) * eta + 2*PI * r * X*X*X / F(phigrav)
               * (Psi0*Psi0 / (X*X) + om*om * exp(-2*Phi) * A*A * F(phigrav) + V(A));
    *rhs_phigrav   = X * eta;
    *rhs_A = Psi0;
    *rhs_eta = - eta * ((*rhs_Phi) - 0.5 * derF(phigrav) * X * eta) - 2 * eta / r + F(phigrav) * X * derW(phigrav) + 
                + 2*PI * X * derF(phigrav) * (1/F(phigrav)) * (om*om * exp(-2*Phi) * A*A * F(phigrav) - Psi0*Psi0 / (X*X) - 2*V(A));
    *rhs_Psi0 = -2 * Psi0 * (1/r) + Psi0 * ((*rhs_X)/X + 1.5 * derF(phigrav) * X * eta - (*rhs_Phi)) - (X*X) * (om*om) * A * exp(-2*Phi) * F(phigrav) + (X*X) * A * Vp(A);
    }

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

void initGrid(int n, double rmax)
  {
  int i;


  for(i = 0; i < n; i++)
    star.r[i] = rmax * i / (n - 1.0);
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

void printPars()
  {
  printf("=======================================\n");
  printf("nint          = %d\n", par.nint);
  printf("rmax          = %g\n", par.rmax);
  printf("omega0        = %g\n", par.omega0);
  printf("Amin          = %g\n", par.Amin);
  printf("Amax          = %g\n", par.Amax);
  printf("mphigrav0     = %g\n", par.mphigrav0);
  printf("phigrav0      = %g\n", par.phigrav0);
  printf("alpha0        = %g\n", par.alpha0);
  printf("beta0         = %g\n", par.beta0);
  printf("nzerotarget   = %d\n", par.nzerotarget);
  printf("thresh        = %g\n", par.thresh);
  printf("mpercentage   = %g\n", par.mpercentage);
  printf("rmatchfac     = %g\n", par.rmatchfac);
  printf("potential     = %s\n", par.potential);
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

void readPars(char* ifil)
  {
  FILE* ifp;
  char  line[LEN];
  int   n, i;


  // First set all parameters to default values
  par.nint      = 401;
  par.rmax      = 1;
  par.omega0    = 1.0;
  par.Amin      = 0.1;
  par.Amax      = 0.1;
  par.mphigrav0 = 1.0;
  par.alpha0    = 3.0;
  par.beta0     = 0.0;
  par.phigrav0  = 0.0;
  par.nmodels   = 1;
  par.nzerotarget = 0;
  par.verbose   = 0;
  par.thresh    = 2e-16;
  par.mpercentage = 90.0;
  par.rmatchfac = 1;
  strcpy(par.potential, "series");
  par.lambda4   = 0;
  par.lambda6   = 0;
  par.lambda8   = 0;
  par.sigma0    = 1;



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
      else if(strstr(line, "omega0") != NULL)
        sscanf(line, "omega0 %le", &(par.omega0));
      else if(strstr(line, "Amax") != NULL)
        sscanf(line, "Amax %le", &(par.Amax));
      else if(strstr(line, "Amin") != NULL)
        sscanf(line, "Amin %le", &(par.Amin));
        else if(strstr(line, "mphigrav0") != NULL)
        sscanf(line, "mphigrav0 %le", &(par.mphigrav0));
      else if(strstr(line, "alpha0") != NULL)
        sscanf(line, "alpha0 %le", &(par.alpha0));
      else if(strstr(line, "beta0") != NULL)
        sscanf(line, "beta0 %le", &(par.beta0));
      else if(strstr(line, "phigrav0") != NULL)
        sscanf(line, "phigrav0 %le", &(par.phigrav0));
      else if(strstr(line, "nmodels") != NULL)
        sscanf(line, "nmodels %d", &(par.nmodels));
      else if(strstr(line, "nzerotarget") != NULL)
        sscanf(line, "nzerotarget %d", &(par.nzerotarget));
      else if(strstr(line, "verbose") != NULL)
        sscanf(line, "verbose %d", &(par.verbose));
      else if(strstr(line, "thresh") != NULL)
        sscanf(line, "thresh %le", &(par.thresh));
      else if(strstr(line, "mpercentage") != NULL)
        sscanf(line, "mpercentage %le", &(par.mpercentage));
      else if(strstr(line, "rmatchfac") != NULL)
        sscanf(line, "rmatchfac %le", &(par.rmatchfac));
      else if(strstr(line, "potential") != NULL)
        sscanf(line, "potential %s", par.potential);
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

void rhsVac(double* rhs_X, double* rhs_Phi, double r, double X, double Phi)
  {
  if(r < 1.0e-15)
    {
    // At r = 0, the gradients of the metric functions are zero
    *rhs_X   = 0;
    *rhs_Phi = 0;
    }
  else
    {
    *rhs_Phi = 0.5 * (X*X - 1) / r;
    *rhs_X   = -0.5 * X * (X*X - 1) / r;
    }
  }

/*==========================================================================*/

int calcExtAsym()
  {
  int    i, k, nzero, sigA, istop, imatch_phigrav, imatch;
  double rstop, c1, c2, c3, c4, epsilon, delta, mass;
  double r, X, Phi, A, eta, Psi0, phigrav, rhs_X, rhs_Phi, rhs_A, rhs_phigrav, rhs_Psi0, rhs_eta, om, mphigrav, dr;
  double dX[5], deta[5], dphigrav[5], dPhi[5], dA[5], dPsi0[5];



  // At this stage we have the correct frequency from the shooting
  // algorithm. Here we will compute this model and also remove
  // the diverging part in the profile and replace it with a smooth
  // exterior. Note that we are not yet rescaling time and frequency
  // yet, since we will need the complete profile with exterior to do
  // that.
  intODE(par.A0, par.phigrav0, par.omega0, &nzero, &rstop, &sigA);
  if(par.verbose)
    {
    printf("---------------------------------------------------------------\n");
    printf("%20.15g          %d          %15.7g           %d\n",
               par.omega0, nzero, rstop, sigA);
    printf("---------------------------------------------------------------\n");
    }

  // We now have a model that either drops in scalar amplitude all
  // the way to the edge of our grid (unlikely) or will diverge
  // somewhere along our grid towards either positive or negative
  // values. We regularize this divergence and replace it with a
  // smooth exterior solution as follows:
  // (1) Search backwards from the truncation point (rstop -- for
  //     which we need to find the index) until we find a minimum
  //     in abs(A). That is the point from which we construct the
  //     exterior.
  istop = iofr(rstop);
  for(i = istop-1; i > 0; i--)
    if(fabs(star.A[i])<fabs(star.A[i+1]) && fabs(star.A[i])<fabs(star.A[i-1]))
      break;
  if(par.verbose)
    printf("A[%d] = %g   A[%d] = %g   A[%d] = %g\n",
           i-1, star.A[i-1], i, star.A[i], i+1, star.A[i+1]);

  if(i == 1)
    { printf("Search for matching point yielded i = 1\n\n"); exit(0); }

  imatch = iofr(star.r[i]*par.rmatchfac);

  if(par.verbose)
    printf("Matching the bosonic field to exterior at   r[%d] = %g\n\n", imatch, star.r[imatch]);

  // (2) We now match the scalar amplitude to an exterior function
  //
  //       A   = c1 * exp(-r * sqrt(1-om^2 / alpha^2)) / r
  //       eta = c1 * exp(-r * sqrt(1-om^2 / alpha^2)) / r
  //
  //     and determine c1 and c2 through matching to A[imatch] and
  //     eta[imatch]. Note that for asymptotic behaviour of eta,
  //     we ignore here a contribution \propto exp(-...) / r^2.
  //     We could formally integrate eta from A, but that deviates
  //     more strongly from the asymptotic limit \eta -> 0 at infinity.
  mphigrav = par.mphigrav0;
  r   = star.r[imatch];
  Phi = star.Phi[imatch];

  // Need to rescale \omega according to asymptotic condition of the lapse function, \alpha
  om  = par.omega0 * sqrt(F(star.phigrav[imatch])) * exp(-star.Phi[imatch]);

  calcMass();

  mass = star.m[imatch];

  epsilon = mass * (1 - 2*om*om)/sqrt(1-om*om);
  delta = mass * mphigrav;

  //Check that h<2k condition is satisfied in the asymptotics
  double condition = mphigrav - 2 * sqrt(1 - om*om);
  
  if (par.verbose)
  {printf("Condition for asymptotics check = %g \n", condition);}

  if (condition >= 0)
  {printf("Asymptotics not satisfied! %g \n", condition);}

  //Find constants of matching
  c1 = star.A[imatch] * pow(r, 1 + epsilon) * exp(r * sqrt(1 - om*om));
  c2 = star.Psi0[imatch] * pow(r, 1 + epsilon) * exp(r * sqrt(1 - om*om));
  c3 = star.phigrav[imatch] * pow(r, 1 + delta) * exp(r * mphigrav);
  c4 = star.eta[imatch] * pow(r, 1 + delta) * exp(r * mphigrav);

  if (par.verbose)
  {
    printf("\nc1, c2 = %g   %g\n", c1, c2);
    printf("c3, c4 = %g   %g\n", c3, c4);
  }

  for(i = imatch+1; i < par.nint; i++)
    {
    dr  = star.r[i] - star.r[i-1];

    // 1st RK step
    r   = star.r[i-1];
    X   = star.X[i-1];
    Phi = star.Phi[i-1];
    phigrav = star.phigrav[i-1];
    eta = star.eta[i-1];
    A   = c1 * exp(-r * sqrt(1 - om*om)) * pow(r, - 1 - epsilon);
    Psi0 = c2 * exp(-r * sqrt(1 - om*om)) * pow(r, - 1 - epsilon);
    rhsBSint(&rhs_X, &rhs_A, &rhs_eta, &rhs_Phi, &rhs_phigrav, &rhs_Psi0, r, X, A, eta, Phi, phigrav, Psi0, om);
    dX[1]   = rhs_X * dr;
    dPhi[1] = rhs_Phi * dr;
    dphigrav[1] = rhs_phigrav * dr;
    deta[1] = rhs_eta * dr;

    // 2nd RK step
    r   = star.r[i-1]   + 0.5 * dr;
    X   = star.X[i-1]   + 0.5 * dX[1];
    Phi = star.Phi[i-1] + 0.5 * dPhi[1];
    phigrav = star.phigrav[i-1] + 0.5 * dphigrav[1];
    eta = star.eta[i-1] + 0.5 * deta[1];
    A   = c1 * exp(-r * sqrt(1 - om*om)) * pow(r, - 1 - epsilon);
    Psi0 = c2 * exp(-r * sqrt(1 - om*om)) * pow(r, - 1 - epsilon);
    rhsBSint(&rhs_X, &rhs_A, &rhs_eta, &rhs_Phi, &rhs_phigrav, &rhs_Psi0, r, X, A, eta, Phi, phigrav, Psi0, om);
    dX[2]   = rhs_X * dr;
    dPhi[2] = rhs_Phi * dr;
    dphigrav[2] = rhs_phigrav * dr;
    deta[2] = rhs_eta * dr;

    // 3rd RK step
    r   = star.r[i-1]   + 0.5 * dr;
    X   = star.X[i-1]   + 0.5 * dX[2];
    Phi = star.Phi[i-1] + 0.5 * dPhi[2];
    phigrav = star.phigrav[i-1] + 0.5 * dphigrav[2];
    eta = star.eta[i-1] + 0.5 * deta[2];
    A   = c1 * exp(-r * sqrt(1 - om*om)) * pow(r, - 1 - epsilon);
    Psi0 = c2 * exp(-r * sqrt(1 - om*om)) * pow(r, - 1 - epsilon);
    rhsBSint(&rhs_X, &rhs_A, &rhs_eta, &rhs_Phi, &rhs_phigrav, &rhs_Psi0, r, X, A, eta, Phi, phigrav, Psi0, om);
    dX[3]   = rhs_X * dr;
    dPhi[3] = rhs_Phi * dr;
    dphigrav[3] = rhs_phigrav * dr;
    deta[3] = rhs_eta * dr;

    // 4th RK step
    r   = star.r[i];
    X   = star.X[i-1]   + dX[3];
    Phi = star.Phi[i-1] + dPhi[3];
    phigrav = star.phigrav[i-1] + dphigrav[3];
    eta = star.eta[i-1] + deta[3];
    A   = c1 * exp(-r * sqrt(1 - om*om)) * pow(r, - 1 - epsilon);
    Psi0 = c2 * exp(-r * sqrt(1 - om*om)) * pow(r, - 1 - epsilon);
    rhsBSint(&rhs_X, &rhs_A, &rhs_eta, &rhs_Phi, &rhs_phigrav, &rhs_Psi0, r, X, A, eta, Phi, phigrav, Psi0, om);
    dX[4]   = rhs_X * dr;
    dPhi[4] = rhs_Phi * dr;
    dphigrav[4] = rhs_phigrav * dr;
    deta[4] = rhs_eta * dr;

    // Update variables
    star.X[i]   = star.X[i-1]   + (dX[1]  + 2*dX[2]  + 2*dX[3]  + dX[4] ) / 6.0;
    star.Phi[i] = star.Phi[i-1] + (dPhi[1]+ 2*dPhi[2]+ 2*dPhi[3]+dPhi[4]) / 6.0;
    star.phigrav[i] = star.phigrav[i-1] + (dphigrav[1]+ 2*dphigrav[2]+ 2*dphigrav[3]+dphigrav[4]) / 6.0;
    star.eta[i] = star.eta[i-1] + (deta[1]+ 2*deta[2]+ 2*deta[3]+deta[4]) / 6.0;
    r   = star.r[i];
    star.A[i]   = c1 * exp(-r * sqrt(1 - om*om)) * pow(r, - 1 - epsilon);
    star.Psi0[i] = c2 * exp(-r * sqrt(1 - om*om)) * pow(r, - 1 - epsilon);
    }

  //Sanity check for NaNs
  if (isnan(c1)|| isnan(c2) || isnan(c3) || isnan(c4) || isnan(star.Phi[par.nint]) || isnan(star.phigrav[par.nint-1]) || isnan(star.X[par.nint]))
  { if (isnan(star.Phi[par.nint]))
          if (par.verbose)
          {printf("Phi is nan on the grid boundary \n");}
    if (isnan(star.phigrav[par.nint-1]))
          if (par.verbose)
          {printf("Gravitational scalar field is nan on the grid boundary \n");}
    if (isnan(star.X[par.nint]))
          if (par.verbose)
          {printf("X is nan on the grid boundary \n");}
    return 0;}
  else
  { return 1;}

  // Overwrite eta with the derivative of A?
  // This is introduces kinks in the eta profile and we do not use it.
  /*for(i = imatch+1; i < par.nint; i++)
    {
    if(i < par.nint-1)
      star.eta[i] = (star.A[i+1] - star.A[i-1]) / (star.r[i+1] - star.r[i-1])
                    / star.X[i];
    else
      star.eta[i] = (star.A[i] - star.A[i-1]) / (star.r[i] - star.r[i-1])
                    / (0.5 * (star.X[i] + star.X[i-1]));
    }*/
  }

/*==========================================================================*/

int iofr(double rtarget)
  {
  int i;


  // Find largest index i such that r[i] < rtarget.
  i = 0;
  for(i = 1; i < par.nint; i++)
    if(star.r[i-1] <= rtarget && star.r[i+1] > rtarget) break;

  if(par.verbose)
    printf("r[%d] = %g   r[%d] = %g   rtarget = %g\n", i-1, star.r[i-1],
           i, star.r[i], rtarget);

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
  if(par.verbose) printf("Phiorg = %g     Phitarget = %g\n", Phiorg, Phitarget);

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
  //   R = (r-m)/2 * { 1 + sqr[1 - m^2/(r-m)^2] }
  //
  // which becomes for large radii
  //
  //   R \approx r - m - m^2 / [4(r-m)]
  //
  // We can invert this to get approximately
  //
  //   r \approx R + m + m^2 / (4R)
  //
  // which I have tested in a 3D simulation of an isotropic
  // Schwarzschild black hole, calculating the areal radius at R=60M
  // (actually M = 1 in this case) and I get
  //
  //   R = 60
  //   r = 61.00403262009653
  //
  //   R + m + m^2 / (4R) = 61.00416666666666666666
  //
  // In the boson star calculation, we have r and want R. So we can fix
  // at the outer edge of the grid
  //
  //   Rtarget = (r-m)/2 * { 1 + sqr[1 - m^2/(r-m)^2] }
  //
  // And then rescale over the whole grid
  //
  //   R -> R * Rtarget / R[n-1]

  m = star.m[n1-1];
  r = star.r[n1-1];
  Rfac = 0.5 * (r / sqrt(F(star.phigrav[n1-1])) - m) * ( 1 + sqrt(1 - F(star.phigrav[n1-1])*m*m / (r - sqrt(F(star.phigrav[n1-1]))*m) / (r - sqrt(F(star.phigrav[n1-1]))*m)) );
  Rfac /= star.R[n1-1];

  for(i = 0; i < n1; i++)
    star.R[i] *= Rfac;
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

double calcRadius()
  {
  int n1, i;


  n1 = par.nint;
  for(i = n1-2; i >= 0; i--)
    if(star.m[i] < par.mpercentage / 100.0 * star.m[n1-1])
      break;

  if(par.verbose)
    {
    printf("total mass = %22.16g\n", star.m[n1-1]);
    printf("star.m[%d] = %22.16g   star.m[%d] = %22.16g\n",
           i, star.m[i], i+1, star.m[i+1]);
    printf("radius = %22.16g\n", star.r[i+1]);
    }

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

void openFile(const char* ofil)
  {
  FILE* ofp;


  ofp = fopen(ofil, "w");
  if(ofp == NULL)
    { printf("Cannot open file   %s\n", ofil); exit(0); }

  fprintf(ofp, "# %s\n", ofil);

  fclose(ofp);
  }

/*==========================================================================*/

void append2Data(const char* ofil, double x, double y)
  {
  FILE* ofp;


  ofp = fopen(ofil, "a");
  if(ofp == NULL)
    { printf("Cannot open file   %s\n", ofil); exit(0); }

  fprintf(ofp, "%22.16g   %22.16g\n", x, y);

  fclose(ofp);
  }

/*==========================================================================*/

void append7Data(const char* ofil, double x, double y, double z, double q, double w, double e, double r)
  {
  FILE* ofp;


  ofp = fopen(ofil, "a");
  if(ofp == NULL)
    { printf("Cannot open file   %s\n", ofil); exit(0); }

  fprintf(ofp, "%22.16g   %22.16g   %22.16g   %22.16g   %22.16g   %22.16g   %22.16g\n", x, y, z, q, w, e, r);

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

int check_phigrav(double entry)
{
  double ratio = - par.alpha0/par.beta0;
  if (fabs(entry - ratio) < 1e-8)
  {return 0;}
  else {return 1;}
}

/*==========================================================================*/