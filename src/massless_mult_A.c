/*
Complete routine for finding BS families of solutions for the massless scalar-tensor theory. 
Here we assume the control parameter is the BS frequency. 
The code increments the given range for the frequency and for each one of them computes a solution.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NRANSI
#define LEN             1000
#define NR_END 1
#define FREE_ARG char*
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

typedef struct model
  {
  double        *r;
  double        *rJ;
  double        *Phi;
  double        *X;
  double        *A;
  double        *eta;
  double        *m;
  double        *mE;
  double        *q;
  double        *R;     // Isotropic radius
  double        *f;     // f := R / r
  double        *psi;   // conformal factor, i.e. g_rr = psi^4
  double        *C;     // Compactness defined as m(r)/r
  double        *Psi0;
  double        *phigrav;
  double        *absphigrav;
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
  double        omegamin;
  double        omegamax;
  int           nzerotarget;
  double        thresh;
  int           nmodels;
  double        mpercentage;
  double        rmatchfac;
  double        lambda4;
  double        lambda6;
  double        lambda8;
  double        sigma0;
  int           docheck;
  int           verbose;
  char          potential[LEN];
  char          minmax[LEN];
  } param;


// Functions
void    append2Data     (const char*, double, double);
void    appendAllData     (const char*, double, double, double, double, double, double, double, double, double, double, double);
void    calculate_model (double*, double*, double*, double*, double*, double*);
int     calcExtAsym     (double, double);
int     nancheck        (double, double, double, double);
double  calculate_criterion_A (int, double);
double  calculate_criterion_phigrav (int, double);
void    debuginfo_FJ    (double**, double**);
void    calcIso         ();
void    calcMass        ();
void    calcMassEinstein        ();
void    calcRJordan     ();
double  calc99RadiusJordan ();
double  calc99RadiusEinstein ();
int     calcRadius_index      ();
double  findMax         (double*, int, int);
void    initGrid        (int, double);
void    intODE          (double, double, double, int*, double*, int*);
int     mygetline       (char*, FILE*);
int     onemodel        (int*, double*, double*, double*, double*, double*, double*, double*);
void    openFile        (const char*);
void    out1D           (double*, double*, double*, int, int, const char*);
void    out_final_model (double, double, double, double, double, double, double, const char* ofil);
void    out_joint       (double*, double*, double*, double*, double*, double*, double*, double*, double*, int, int, const char*);
void    printPars       ();
void    readPars        (char*);
void    registerPotential();
double    rescalePhi      (double);
void    rhsBSint        (double*, double*, double*, double*, double*, double*, double, double, double, double, double, double, double, double);
void    rhsIso          (double*, double, double, double, double);
int     iofr            (double);
double  V_series        (double);
double  Vp_series       (double);
double  V_solitonic     (double);
double  Vp_solitonic    (double);
double  F_st               (double);
double  derF_st            (double);
double  W_st               (double);
double  derW_st            (double);
void    integrandN        (double*, double, double, double, double, double, double);
void    calculateN      (double);
void    check_phigrav      (double);
void    out_old_solution    (double, const char*);

// Numerical Recipes functions
double  **dmatrix       (long, long, long, long);
void    free_dmatrix    (double**, long, long, long, long);
void    free_ivector    (int*, long, long);
void    gaussj          (double**, int, double**, int);
int     *ivector        (long, long);
void    nrerror         (char[]);
void    test_matrix     ();
 
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
  int i, k, nzero;
  int converged;
  int success;
  double omBS, mBS, mBSE, rBS, rBSE, CBS, phigravmax, omegaini, domega;

  if(argc != 2) { printf("Usage:   %s   <parfile>\n\n", argv[0]); exit(0); }

  // Parameters
  readPars(argv[1]);
  printPars();

  /*printf("V_series = %g\n", V_series(0.1));
  printf("V_solitonic = %g\n", V_solitonic(0.1));
  exit(0);*/

  // Register functions
  registerPotential();

  // Allocate memory
  star.r   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.rJ   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.Phi = (double*) malloc((size_t) par.nint * sizeof(double));
  star.X   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.A   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.Psi0   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.phigrav   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.absphigrav   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.eta = (double*) malloc((size_t) par.nint * sizeof(double));
  star.m   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.mE   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.q   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.R   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.f   = (double*) malloc((size_t) par.nint * sizeof(double));
  star.psi = (double*) malloc((size_t) par.nint * sizeof(double));
  star.C   = (double*) malloc((size_t) par.nint * sizeof(double));  

  // Grid setup
  initGrid(par.nint, par.rmax);

  openFile("MofR.dat");
  openFile("MofA.dat");
  openFile("Mofphigrav.dat");
  openFile("AuxiliaryInfo.dat");
  printf("\n\n   A0           phigrav0           omega               mass            radius      max(C)      nzero      phigravmax\n");
  printf("========================================================================================================================\n");

  if(par.nmodels == 1) domega = 0;
  else domega = (par.omegamax - par.omegamin) / (par.nmodels - 1);
  for(i = 0; i < par.nmodels; i++)
    {
    par.omega0 = par.omegamin + domega * i;
    printf("Setting omega0 to %22.16g\n", par.omega0);
    // printf("Initiliasing A0 to: %g", par.A0);
    success = onemodel(&nzero, &rBS, &rBSE, &omBS, &mBS, &CBS, &phigravmax, &omegaini);
    if (success==0)
    {
      printf("No solution found for %11.6g\n", par.A0);
    }
    else
    {

      printf("%11.6g  %11.6g   %16.10g   %16.10g   %12.6g   %12.6g     %d   %12.6g\n",
           par.A0, par.phigrav0, omBS, mBS, rBS, CBS, nzero, phigravmax);

      append2Data("MofR.dat", rBS, mBS);
      append2Data("MofA.dat", par.A0, mBS);
      append2Data("Aofomega.dat", par.A0, omBS);
      append2Data("Mofphigrav.dat", par.phigrav0, mBS);
      double alpha0 = exp(star.Phi[0])/sqrt(F(star.phigrav[0]));
      double QBS = star.q[par.nint-1];
      appendAllData("AuxiliaryInfo.dat", omBS, alpha0, par.phigrav0, par.A0, mBS, rBS, CBS, phigravmax, QBS, nzero, omegaini);
    }
    }

  }

/*==========================================================================*/

int onemodel(int* nzero, double* rBS, double* rBSE, double* omBS, double* mBS, double* CBS, double* phigravmax, double* omega_init)
  {

  int i, k;
  int converged;
  double criterion_derivative, criterion_derivative_A;

  // Shoot
  // shoot();

  for (k = 0; k < 1001; ++k)
  { 
    if (par.verbose)
    {
    printf("--------------------------------------------------------------------------------------------------\n");
    printf("I am starting iteration number: %d\n", k);
    printf("Gravitational scalar field guess is set to %22.16g\n", par.phigrav0);
    printf("Amplitude of the boson star is set to %22.16g\n", par.A0);
    printf("--------------------------------------------------------------------------------------------------\n");
    }

    if (k==1000)
    {printf("I've reached maximum number of iterations \n"); return 0;}

    // Allocate Newton-Raphson matrices
    double **J, **F;
    J = dmatrix(1,2,1,2);
    F = dmatrix(1,2,1,1);

    int jj = par.nint - 100;

    // Calculate exterior
    converged = calcExtAsym(par.A0, par.phigrav0);

    double criterion_A = calculate_criterion_A(converged, par.omega0);
    double criterion = calculate_criterion_phigrav(jj, par.omega0);

    printf("The matching difference for the gravitational scalar field ...%g\n", criterion);
    printf("The matching difference for the bosonic scalar field ... %g\n", criterion_A);

    if (fabs(criterion_A) < 1e-7 && fabs(criterion) < 1e-04)
    {
      *omega_init = par.omega0;

      if (par.verbose) 
            {
              printf("--------------------------------------------------------------------------------------------------\n");
              printf("I get the following difference for the gravitational scalar field, %22.16g\n", fabs(criterion));
              printf("I get the following difference for the bosonic scalar field, %22.16g\n", fabs(criterion_A));
              printf("--------------------------------------------------------------------------------------------------\n");
            }
          calculate_model(rBS, rBSE, omBS, mBS, CBS, phigravmax);

          *nzero = 0;

          for (i = 0; i < par.nint; i++)
          {
            if(star.A[i] * star.A[i-1] < 0) 
            {
              (*nzero)++;   // Zero crossing found
            }
          }

          for(i = 1; i < par.nint; i++)
           {
            // Check for nan; we terminate then
            if(nancheck(star.A[i], star.phigrav[i], star.Phi[i], par.omega0))
            {
              printf("I have found NaNs in some variables... Oh no!\n");
              printf("A, phigrav, Phi, omega0 = %g %g %g %g\n",star.A[i], star.phigrav[i], star.Phi[i], par.omega0);
              exit(0);
            }
           }
           
          break;
    }

    F[1][1] = criterion_A;
    F[1][2] = criterion;
    
    if (par.verbose)
    {
      out1D(star.r, star.rJ, star.phigrav, 0, par.nint-1, "phigrav_temp.dat");
      out1D(star.r, star.rJ, star.A, 0, par.nint-1, "A_temp.dat");
    }

    double domega = par.omega0 + 1e-14;
    double dA = par.A0 + 1e-14;
    double dphigrav = par.phigrav0 + 1e-14;

    converged = calcExtAsym(dA, par.phigrav0);
    criterion_A = calculate_criterion_A(converged, par.omega0);
    criterion = calculate_criterion_phigrav(jj, par.omega0);

    J[1][1] = (criterion_A - F[1][1]) / 1e-14;
    J[2][1] = (criterion - F[1][2]) / 1e-14;

    converged = calcExtAsym(par.A0, dphigrav);
    criterion_A = calculate_criterion_A(converged, par.omega0);
    criterion = calculate_criterion_phigrav(jj, par.omega0);

    J[1][2] = (criterion_A - F[1][1]) / 1e-14;
    J[2][2] = (criterion   - F[1][2]) / 1e-14;

    debuginfo_FJ(F, J);

    gaussj(J, 2, F, 1);

    debuginfo_FJ(F, J);

    par.A0 -= F[1][1];
    par.phigrav0 -= F[1][2];

    if (par.A0 < 0)
    {par.A0 = fabs(par.A0);}

  }

}

/*==========================================================================*/

void calculate_model(double* rBS, double* rBSE, double* omBS, double* mBS, double* CBS, double* phigravmax)
{
  int i;
  double mBSE;

  *omBS = rescalePhi(par.omega0);

  calculateN(*omBS);
  
  // Convert to Jordan radius
  calcRJordan();
  
  // Compute diagnostic quantities.
  calcMass();
  calcMassEinstein();

  *rBS  = calc99RadiusJordan();
  *rBSE  = calc99RadiusEinstein();
  *mBS  = star.m[par.nint-1];
  mBSE  = star.mE[par.nint-1];
  *CBS  = findMax(star.C, 0, par.nint-1);

  for (i = 0; i < par.nint; ++i)
  {
    star.absphigrav[i] = fabs(star.phigrav[i]);
  }

  *phigravmax = findMax(star.absphigrav, 0, par.nint-1);
  
  // if (par.verbose) {printf("Maximal compactness:   %g\n", CBS);}

  // if (par.verbose)
  // {printf("\n===================================================\n");
  // printf("Physical frequency:   omega = %22.16g\n", omBS);
  // printf("Total mass:           m     = %22.16g\n", star.m[par.nint-1]);
  // printf("===================================================\n");}

}

/*==========================================================================*/

int calcExtAsym(double Amp, double phigrav_guess)
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
  intODE(Amp, phigrav_guess, par.omega0, &nzero, &rstop, &sigA);
  // if (par.verbose)
  // {
  //   printf("-----------------------------------------------------------------------\n");
  //   printf("%22.16g          %d          %15.7g           %d\n",
  //        omega, nzero, rstop, sigA);
  //   printf("-----------------------------------------------------------------------\n");
  // }

  if (rstop == par.rmax)
  {return 1;}

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
  if (par.verbose)
  {
  printf("\nA[%d] = %g   A[%d] = %g   A[%d] = %g\n",
         i-1, star.A[i-1], i, star.A[i], i+1, star.A[i+1]);
  }

  // if(i ==0 || i == 1)
  //   {printf("Search for matching point yielded i = 0 or i = 1\n\n"); exit(0); }

  imatch = iofr(star.r[i]*par.rmatchfac);

  if (par.verbose) {printf("Matching the bosonic field to exterior at   r[%d] = %g\n\n", imatch, star.r[imatch]);}

  // (2) We now match the scalar amplitude to an exterior function

  mphigrav = par.mphigrav0;
  r   = star.r[imatch];
  Phi = star.Phi[imatch];

  // Need to rescale \omega according to asymptotic condition of the lapse function, \alpha
  om  = par.omega0 * sqrt(F(star.phigrav[imatch])) * exp(-star.Phi[imatch]);
  // double Phitarget = -log(sqrt(F(star.phigrav[imatch])) * star.X[imatch]);
  // om = par.omega0 * exp(Phitarget-star.Phi[imatch]);

  calcMassEinstein();

  mass = star.mE[imatch];

  epsilon = mass * (1 - 2*om*om)/sqrt(1-om*om);
  delta = mass * mphigrav;

  //Check that h<2k condition is satisfied in the asymptotics
  double condition = mphigrav - 2 * sqrt(1 - om*om);

  if (par.verbose) {printf("Condition for asymptotics check = %g \n", condition);}

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

  //Integrate outwards once again
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
          {if (par.verbose) {printf("Phi is nan on the grid boundary \n");}}
    if (isnan(star.phigrav[par.nint-1]))
          {if (par.verbose) {printf("Gravitational scalar field is nan on the grid boundary \n");}}
    if (isnan(star.X[par.nint]))
          {if (par.verbose) {printf("X is nan on the grid boundary \n");}}
    return 0;}
  else
  { return imatch;}
  
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

double calculate_criterion_A(int index, double omega)
{
  calcMassEinstein();

  double om = omega * exp(-star.Phi[index]);

  double epsilon = star.mE[index] * (1 - 2*om*om)/sqrt(1-om*om);
  double factor = pow(star.r[index], 2+epsilon) * exp(sqrt(1 - om*om) * star.r[index]);
  double denominator = - 1 - epsilon - star.r[index] * sqrt(1 - om*om);
  double derivative_term = (star.A[index] - star.A[index-1])/(star.r[index] - star.r[index-1]);
  double constantA = factor * derivative_term / denominator;
  double A_surface = constantA * pow(star.r[index], - 1 - epsilon) * exp(-sqrt(1 - om*om) * star.r[index]);
  double criterion_A = star.A[index] - A_surface;

  return criterion_A;
}

/*==========================================================================*/

double calculate_criterion_phigrav(int index, double omega0)
{
  double omega;

  // Surface criterion of https://iopscience.iop.org/article/10.1088/0264-9381/33/13/135002/pdf Eq.(3.15).
  double derivative_term = (star.Phi[index] - star.Phi[index-1])/(star.r[index] - star.r[index-1]);
  double sqrt_expr = sqrt(derivative_term * derivative_term + star.X[index] * star.X[index] * star.eta[index] * star.eta[index]);
  double varphi_surface = 0.0 - (star.X[index] * star.eta[index]) / sqrt_expr * atanh(sqrt_expr / (derivative_term + 1 / star.r[index]));
  double criterion = star.phigrav[index] - varphi_surface;

  return criterion;
}

/*==========================================================================*/

int iofr(double rtarget)
  {
  int i;


  // Find largest index i such that r[i] < rtarget.
  i = 0;
  for(i = 1; i < par.nint; i++)
    if(star.r[i-1] <= rtarget && star.r[i] > rtarget) break;

  //printf("\nMatching to Exterior: r[%d] = %g   r[%d] = %g   rtarget = %g\n\n",
  //       i-1, star.r[i-1], i, star.r[i], rtarget);

  return i;
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

  *nzero = 0;
  *rstop = star.r[n1-1];
  istop  = -1;   // So we have no vacuum region unless we exceed the amplitude

  // Central values
  star.X[0]   = 1/sqrt(F(phigrav0));
  star.A[0]   = A0;
  star.phigrav[0]   = phigrav0;
  star.Psi0[0]   = 0;
  star.eta[0] = 0;
  star.Phi[0] = 0.;   // Phi has a free constant we later match to Schwarzschild

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
    }


  if(istop < 0)
    {
    *sigAstop = ( (star.A[n1-1] > 0) - (star.A[n1-1] < 0) );
    if (par.verbose) {printf("No need for vacuum \n");}
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

void integrandN(double *rhs_q, double r, double A, double omega, double X, double phigrav, double Phi)
  {
    double alpha;

    alpha = exp(Phi)/sqrt(F(phigrav));

    //Integrand in the Noether charge formulae
    *rhs_q = 4 * PI * r * r * (A * A * omega * X) / (F(phigrav) * alpha);
  }

/*==========================================================================*/

void calculateN(double omega)
  {
    int i;
    double dr, q, r, rhs_q, A, X, phigrav, Phi;
    double dq[5];

    star.q[0] = 0;

 	  for (i = 1; i < par.nint; i++)
    	{
        dr  = star.r[i] - star.r[i-1];

        A = star.A[i-1];
        X = star.X[i-1];
        phigrav = star.phigrav[i-1];
        Phi = star.Phi[i-1];

        // 1st RK step
        r = star.r[i-1];
        integrandN(&rhs_q, r, A, omega, X, phigrav, Phi);
		    dq[1] = rhs_q * dr;

        // 2nd RK step
        r   = star.r[i-1] + 0.5 * dr;
				integrandN(&rhs_q, r, A, omega, X, phigrav, Phi);
        dq[2] = rhs_q * dr;

        // 3rd RK step
        r   = star.r[i-1] + 0.5 * dr;
				integrandN(&rhs_q, r, A, omega, X, phigrav, Phi);
        dq[3] = rhs_q * dr;

        // 4th RK step
        r   = star.r[i];
        integrandN(&rhs_q, r, A, omega, X, phigrav, Phi);
        dq[4] = rhs_q * dr;

		    star.q[i] = star.q[i-1] + (dq[1] + 2*dq[2] + 2*dq[3] + dq[4]) / 6;
    	}
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
  par.mphigrav0 = 1.0;
  par.alpha0    = 3.0;
  par.beta0     = 0.0;
  par.phigrav0  = 0.0;
  par.nmodels   = 1;
  par.omega0    = 1.;            // omega0 is always 1, it is not specified
  par.nzerotarget = 0;
  par.thresh    = 2e-16;
  par.mpercentage = 90;
  par.rmatchfac = 1;
  strcpy(par.potential, "series");
  par.lambda4   = 0;
  par.lambda6   = 0;
  par.lambda8   = 0;
  par.sigma0    = 1;
  par.docheck   = 1;
  par.verbose   = 1;
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
      else if(strstr(line, "A0") != NULL)
        sscanf(line, "A0 %le", &(par.A0));
      else if(strstr(line, "omax") != NULL)
        sscanf(line, "omax %le", &(par.omegamax));
      else if(strstr(line, "omin") != NULL)
        sscanf(line, "omin %le", &(par.omegamin));
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
      else if(strstr(line, "docheck") != NULL)
        sscanf(line, "docheck %d", &(par.docheck));
      else if(strstr(line, "verbose") != NULL)
        sscanf(line, "verbose %d", &(par.verbose));
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
  printf("omin          = %g\n", par.omegamin);
  printf("omax          = %g\n", par.omegamax);
  printf("nmodels       = %d\n", par.nmodels);
  printf("mphigrav0     = %g\n", par.mphigrav0);
  printf("phigrav0      = %g\n", par.phigrav0);
  printf("omega0      = %g\n", par.omega0);
  printf("alpha0        = %g\n", par.alpha0);
  printf("beta0         = %g\n", par.beta0);
  printf("nzerotarget   = %d\n", par.nzerotarget);
  printf("docheck        = %d\n", par.docheck);
  printf("thresh        = %g\n", par.thresh);
  printf("mpercentage   = %g\n", par.mpercentage);
  printf("rmatchfac     = %g\n", par.rmatchfac);
  printf("potential     = %s\n", par.potential);
  printf("minmax        = %s\n", par.minmax);
  printf("verbose          = %d\n", par.verbose);
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

double rescalePhi(double omega0)
  {
  int i, n1;
  double Phitarget, Phiorg, omega;


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
  
  // if (par.verbose) {printf("Phiorg = %g     Phitarget = %g\n", Phiorg, Phitarget);}

  for(i = 0; i < n1; i++)
    star.Phi[i] += (Phitarget - Phiorg);

  omega = omega0 * exp(Phitarget - Phiorg);

  return omega;
  }

/*==========================================================================*/

void calcMass()
  {
  int i;
  double r, rJ, X, phigrav;


  for(i = 0; i < par.nint; i++)
    {
    r = star.r[i];
    rJ = star.rJ[i];
    X = star.X[i];
    phigrav = star.phigrav[i];
    star.m[i] = 0.5 * rJ * (1 - 1 / (F(phigrav)*X*X));
    if(i == 0)
      star.C[i] = 0;
    else
      star.C[i] = star.m[i] / star.rJ[i];
    }
  }

/*==========================================================================*/

void calcMassEinstein()
  {
  int i;
  double r, rJ, X, phigrav;


  for(i = 0; i < par.nint; i++)
    {
    r = star.r[i];
    X = star.X[i];
    phigrav = star.phigrav[i];
    star.mE[i] = 0.5 * r * (1 - 1 / (F(phigrav)*X*X));
    }
  }

/*==========================================================================*/

void calcRJordan()
  {
  int   i;

  for(i = 1; i <= par.nint; i++)
	  {
      star.rJ[i] = star.r[i]/sqrt(F(star.phigrav[i]));
	  }

  }

/*==========================================================================*/

double calc99RadiusJordan()
  {
  int n1, i;


  n1 = par.nint;
  for(i = n1-2; i >= 0; i--)
    if(star.m[i] < par.mpercentage / 100.0 * star.m[n1-1])
      break;

  if (par.verbose)
  {
  // printf("total mass = %22.16g\n", star.m[n1-1]);
  // printf("star.m[%d] = %22.16g   star.m[%d] = %22.16g\n",
  //        i, star.m[i], i+1, star.m[i+1]);
     printf("Jordan radius = %22.16g\n", star.rJ[i+1]);
  }

  return star.rJ[i+1];
  }

/*==========================================================================*/

double calc99RadiusEinstein()
  {
  int n1, i;


  n1 = par.nint;
  for(i = n1-2; i >= 0; i--)
    if(star.mE[i] < par.mpercentage / 100.0 * star.mE[n1-1])
      break;

  if (par.verbose)
  {
  // printf("total mass = %22.16g\n", star.m[n1-1]);
  // printf("star.m[%d] = %22.16g   star.m[%d] = %22.16g\n",
  //        i, star.m[i], i+1, star.m[i+1]);
     printf("Einstein radius = %22.16g\n", star.r[i+1]);
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
  // if (par.verbose) {printf("matching radii at   r[%d] = %g\n", n1-1, star.r[n1-1]);}
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

void out1D(double* x, double* xx, double* y, int n1, int n2, const char* ofil)
  {
  int i;
  FILE* ofp;


  ofp = fopen(ofil, "w");
  if(ofp == NULL) { printf("Cannot open %s in out1D\n\n", ofil); exit(0); }

  fprintf(ofp, "# %s\n", ofil);

  for(i = n1; i <= n2; i++)
    fprintf(ofp, "%22.16g   %22.16g   %22.16g\n", x[i], xx[i], y[i]);

  fclose(ofp);
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

void appendAllData(const char* ofil, double x, double y, double z, double q, double w, double e, double r, double u, double j, double i, double o)
  {
  FILE* ofp;


  ofp = fopen(ofil, "a");
  if(ofp == NULL)
    { printf("Cannot open file   %s\n", ofil); exit(0); }

  fprintf(ofp, "%22.16g   %22.16g   %22.16g   %22.16g   %22.16g   %22.16g   %22.16g    %22.16g    %22.16g    %22.16g    %22.16g\n", x, y, z, q, w, e, r, u, j, i, o);

  fclose(ofp);
  }

/*==========================================================================*/

void out_final_model(double x, double y, double z, double w, double a, double s, double r, const char* ofil)
  {
  FILE* ofp;


  ofp = fopen(ofil, "w");
  if(ofp == NULL) { printf("Cannot open %s in out_final_model\n\n", ofil); exit(0); }

  fprintf(ofp, "# %s\n", ofil);
  fprintf(ofp, "#%22.16s   %22.16s   %22.16s   %22.16s   %22.16s   %22.16s   %22.16s\n", "A0", "phigrav0", "omega", "star mass", "compactness", "radius", "phigrav max" );

  fprintf(ofp, "%22.16g   %22.16g   %22.16g   %22.16g   %22.16g   %22.16g   %22.16g\n", x, y, z, w, a, s, r);

  fclose(ofp);
  }

/*==========================================================================*/
void out_joint(double* x, double* xx, double* y, double* q, double* w, double* a, double* s, double* e, double* r, int n1, int n2, const char* ofil)
  {
  int i;
  FILE* ofp;


  ofp = fopen(ofil, "w");
  if(ofp == NULL) { printf("Cannot open %s in out_joint\n\n", ofil); exit(0); }

  fprintf(ofp, "# %s\n", ofil);

  fprintf(ofp, "#%22.16s   %22.16s   %22.16s   %22.16s   %22.16s   %22.16s   %22.16s   %22.16s   %22.16s   %22.16s\n", "radius E", "radius J", "A", "phigrav", "omega", "X", "Phi", "Psi0", "eta", "star mass");
  for(i = n1; i <= n2; i++)
    fprintf(ofp, "%22.16g   %22.16g   %22.16g   %22.16g   %22.16g   %22.16g   %22.16g   %22.16g   %22.16g   %22.16g\n", x[i], xx[i], y[i], q[i], par.omega0, w[i], a[i], s[i], e[i], r[i]);

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

void check_phigrav(double entry)
{
  double ratio = - par.alpha0/par.beta0;
  if (fabs(entry - ratio) < 1e-8)
  {
    printf("Not allowed to cross -a0/b0 line! Exiting... \n");
    exit(0);
  }
  else {return;}
}

/*==========================================================================*/
void out_old_solution(double x, const char* ofil)
{
double phigravmax;
FILE* ofp;

ofp = fopen(ofil, "w");
  if(ofp == NULL) { printf("Cannot open %s in out_old_solution\n\n", ofil); exit(0); }

phigravmax = findMax(star.phigrav, 0, par.nint-1);

if (phigravmax > 0.1)
{
  fprintf(ofp, "%22.16g", x);
  fclose(ofp);
}

}

/*==========================================================================*/

int nancheck(double y1, double y2, double y3, double y4)
  {
  if( (y1 != y1) || (y2 != y2) || (y3 != y3) || (y4 != y4))
    return 1;
  else
    return 0;
  }

/*==========================================================================*/

// Numerical Recipes functions below

/*==========================================================================*/

void nrerror(char error_text[])
  {
  /* Numerical Recipes standard error handler */
  fprintf(stderr, "Numerical Recipes run-time error...\n");
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, "...now exiting to system...\n");
  exit(1);
  }

/*==========================================================================*/

int *ivector(long nl, long nh)
  {
  /* allocate an int vector with subscript range v[nl..nh] */
  int *v;


  v = (int*) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if(!v) nrerror("allocation failure in ivector()");
  return v-nl+NR_END;
  }

/*==========================================================================*/

void free_ivector(int* v, long nl, long nh)
  {
  /* free an int vector allocated with ivector() */
  free((FREE_ARG) (v+nl-NR_END));
  }

/*==========================================================================*/

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
  {
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m = (double**) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
  if(!m) nrerror("allocation failure 1 in dmatrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (double*) malloc((size_t) ((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in dmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i = nrl+1; i <= nrh; i++) m[i] = m[i-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
  }

/*==========================================================================*/

void free_dmatrix(double** m, long nrl, long nrh, long ncl, long nch)
  {
  /* free a double matrix allocated by dmatrix() */
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
  }

/*==========================================================================*/

void test_matrix()
  {
  double **J, **E;


  J = dmatrix(1,3,1,3);
  E = dmatrix(1,3,1,1);

  J[1][1] = -1.0;       J[1][2] = 2.0;          J[1][3] = 3.0;
  J[2][1] = 2.0;        J[2][2] = 3.0;          J[2][3] = 4.0;
  J[3][1] = 3.0;        J[3][2] = 4.0;          J[3][3] = 5.0;

  E[1][1] = 8.0;        E[1][2] = 8.0;          E[1][3] = 10.0;

  // Invert matrix
  gaussj(J, 3, E, 1);

  printf("J[1][...] = %g   %g   %g\n", J[1][1], J[1][2], J[1][3]);
  printf("J[2][...] = %g   %g   %g\n", J[2][1], J[2][2], J[2][3]);
  printf("J[3][...] = %g   %g   %g\n", J[3][1], J[3][2], J[3][3]);
  printf("E[1][...] = %g   %g   %g\n", E[1][1], E[1][2], E[1][3]);
  exit(0);
  }

/*==========================================================================*/

void gaussj(double **a, int n, double **b, int m)
{
        int *indxc,*indxr,*ipiv;
        int i,icol,irow,j,k,l,ll;
        double big,dum,pivinv,temp;

        indxc = (int*) ivector(1,n);
        indxr = (int*) ivector(1,n);
        ipiv  = (int*) ivector(1,n);
        for (j=1;j<=n;j++) ipiv[j]=0;
        for (i=1;i<=n;i++) {
                big=0.0;
                for (j=1;j<=n;j++)
                        if (ipiv[j] != 1)
                                for (k=1;k<=n;k++) {
                                        if (ipiv[k] == 0) {
                                                if (fabs(a[j][k]) >= big) {
                                                        big=fabs(a[j][k]);
                                                        irow=j;
                                                        icol=k;
                                                }
                                        } else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1");
                                }
                ++(ipiv[icol]);
                if (irow != icol) {
                        for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
                        for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
                }
                indxr[i]=irow;
                indxc[i]=icol;
                if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix-2");
                pivinv=1.0/a[icol][icol];
                a[icol][icol]=1.0;
                for (l=1;l<=n;l++) a[icol][l] *= pivinv;
                for (l=1;l<=m;l++) b[icol][l] *= pivinv;
                for (ll=1;ll<=n;ll++)
                        if (ll != icol) {
                                dum=a[ll][icol];
                                a[ll][icol]=0.0;
                                for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
                                for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
                        }
        }
        for (l=n;l>=1;l--) {
                if (indxr[l] != indxc[l])
                        for (k=1;k<=n;k++)
                                SWAP(a[k][indxr[l]],a[k][indxc[l]]);
        }
        free_ivector(ipiv,1,n);
        free_ivector(indxr,1,n);
        free_ivector(indxc,1,n);
}

/*==========================================================================*/

void debuginfo_FJ(double** F, double** J)
  {
  printf("------------------------------------\n");
  printf("F[.][1] = %15.8g   %15.8g\n",
         F[1][1],F[1][2]);
  printf("J[1][.] = %15.8g   %15.8g\n",
         J[1][1],J[1][2]);
  printf("J[2][.] = %15.8g   %15.8g\n",
         J[2][1],J[2][2]);
  printf("------------------------------------\n");
  }

/*==========================================================================*/
