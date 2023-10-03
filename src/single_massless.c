/*
Complete routine for finding a BS solution in the massless scalar-tensor theory with Damour-Esposito-Farese coupling function.
Here we apply a Newton-Raphson to fine-tune the BS frequency and central gravitational scalar field.
Diagnostics such as ADM mass, Noether charge, BS radius, compactness, maximum of the gravitational scalar field and trace of the energy momentum tensor are computed in the process.

Run this programme with single_massless.par file!

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
  double        *r;     //Einstein radius
  double        *rJ;    //Jordan radius
  double        *Phi;   //See Eq.(2.3) of the notes on the infrastructure
  double        *X;     //Metric function
  double        *A;     //Bosonic amplitude
  double        *eta;   //See Eq.(2.1)
  double        *m;     //Mass function in Jordan frame
  double        *mE;    //Mass function in Einstein frame
  double        *q;     //Noether charge
  double        *R;     // Isotropic radius
  double        *f;     // f := R / r
  double        *psi;   // conformal factor, i.e. g_rr = psi^4
  double        *C;     // Compactness defined as m(r)/r
  double        *Psi0;  //See Eq.(2.2)
  double        *phigrav; //Gravitational scalar field
  double        *T;     // trace of energy-momentum tensor
  } model;


typedef struct param
  {
  double        rmax;
  int           nint;
  double        omega0;
  double        mphigrav0;
  double        A0;
  double        Phi0;
  double        alpha0;
  double        beta0;
  double        phigrav0;
  double        thresh;
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
void    calculate_model (double*, double*, double*, double*, double*, double*, double*);
int     calcExtAsym     (double, double, double);
double  calculate_criterion_A (int, double);
int     nancheck        (double, double, double, double);
double  calculate_criterion_phigrav (int, double);
void    debuginfo_FJ    (double**, double**);
void    calcIso         ();
void    calcMass        ();
void    calcMassEinstein        ();
void    calcRJordan     ();
double  calc99RadiusJordan ();
double  calc99RadiusEinstein ();
int     calcRadius_index      ();
void    calculate_trace_t     (double);
double  findMax         (double*, int, int);
void    initGrid        (int, double);
void    intODE          (double, double, double, int*, double*, int*);
int     mygetline       (char*, FILE*);
void    out1D           (double*, double*, double*, int, int, const char*);
void    appendAllData     (const char*, double, double, double, double, double, double, double, double);
void    out_final_model (double, double, double, double, double, double, double, int, double, double, const char* ofil);
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
  int l;
  int i, k, zeros;
  int converged;
  double omBS, mBS, mBSE, rBS, rBSE, CBS, QBS, phigravmax, omegaini;
  double criterion_derivative, criterion_derivative_A, criterion_derivative_omega;

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
  star.T   = (double*) malloc((size_t) par.nint * sizeof(double));  
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

  for (k = 0; k < 1001; ++k)
  { 
    if (par.verbose)
    {
    printf("--------------------------------------------------------------------------------------------------\n");
    printf("I am starting iteration number: %d\n", k);
    printf("Gravitational scalar field guess is set to %22.16g\n", par.phigrav0);
    printf("Amplitude of the boson star is set to %22.16g\n", par.A0);
    printf("Omega of the boson star is set to %22.16g\n", par.omega0);
    printf("--------------------------------------------------------------------------------------------------\n");
    }

    if (k==1000)
    {printf("I've reached maximum number of iterations \n"); exit(0);}

    // Allocate Newton-Raphson matrices
    double **J, **F;
    J = dmatrix(1,2,1,2);
    F = dmatrix(1,2,1,1);

    int jj = par.nint - 100;

    // Calculate exterior
    converged = calcExtAsym(par.omega0, par.A0, par.phigrav0);

    double criterion_A = calculate_criterion_A(converged, par.omega0);
    double criterion = calculate_criterion_phigrav(jj, par.omega0);

    printf("The matching difference for the gravitational scalar field ...%g\n", criterion);
    printf("The matching difference for the bosonic scalar field ... %g\n", criterion_A);

    if (fabs(criterion_A) < 1e-7 && fabs(criterion) < 1e-04)
    {
      if (par.verbose) 
            {
              printf("--------------------------------------------------------------------------------------------------\n");
              printf("I get the following difference for the gravitational scalar field, %22.16g\n", fabs(criterion));
              printf("I get the following difference for the bosonic scalar field, %22.16g\n", fabs(criterion_A));
              printf("--------------------------------------------------------------------------------------------------\n");
            }
          calculate_model(&rBS, &rBSE, &omBS, &mBS, &CBS, &phigravmax, &omegaini);

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

    converged = calcExtAsym(domega, par.A0, par.phigrav0);
    criterion_A = calculate_criterion_A(converged, domega);
    criterion = calculate_criterion_phigrav(jj, domega);

    J[1][1] = (criterion_A - F[1][1]) / 1e-14;
    J[2][1] = (criterion - F[1][2]) / 1e-14;

    converged = calcExtAsym(par.omega0, par.A0, dphigrav);
    criterion_A = calculate_criterion_A(converged, par.omega0);
    criterion = calculate_criterion_phigrav(jj, par.omega0);

    J[1][2] = (criterion_A - F[1][1]) / 1e-14;
    J[2][2] = (criterion   - F[1][2]) / 1e-14;

    //Perform Gaussian elimination
    gaussj(J, 2, F, 1);

    //Update initial guesses
    par.omega0 -= F[1][1];
    par.phigrav0 -= F[1][2];

  }

  if (par.verbose) {printf("Maximal compactness:   %g\n", CBS);}

  // IO
  out1D(star.r, star.rJ, star.X, 0, par.nint-1, "X_massless.dat");
  out1D(star.r, star.rJ, star.phigrav, 0, par.nint-1, "phigrav_massless.dat");
  out1D(star.r, star.rJ, star.T, 0, par.nint-1, "T_massless.dat");
  out1D(star.r, star.rJ, star.Psi0, 0, par.nint-1, "Psi0_massless.dat");
  out1D(star.r, star.rJ, star.A, 0, par.nint-1, "A_massless.dat");
  out1D(star.r, star.rJ, star.eta, 0, par.nint-1, "eta_massless.dat");
  out1D(star.r, star.rJ, star.Phi, 0, par.nint-1, "Phi_massless.dat");
  out1D(star.r, star.rJ,star.m, 0, par.nint-1, "m_massless.dat");
  out1D(star.r, star.rJ, star.q, 0, par.nint-1, "q_massless.dat");
  out1D(star.r, star.rJ,star.mE, 0, par.nint-1, "mE_massless.dat");

  zeros = 0;

  for (l = 0; l < par.nint; l++)
  {
    if(star.A[l] * star.A[l-1] < 0) 
    {
      (zeros)++;   // Zero crossing found
    }
  }

  QBS = star.q[par.nint-1];
  printf("The Noether charge of a BS is %g\n", QBS);

  out_final_model(par.A0, par.phigrav0, omBS, mBS, CBS, rBS, phigravmax, zeros, omegaini, QBS, "final_model_data_massless.dat");  
  out_joint(star.r, star.rJ, star.A, star.phigrav, star.X, star.Phi, star.Psi0, star.eta, star.m, 0, par.nint-1, "joint_profiles.dat");
  double alpha0 = exp(star.Phi[0])/sqrt(F(star.phigrav[0]));
  appendAllData("single_AuxiliaryInfo_massless.dat", omBS, alpha0, par.phigrav0, par.A0, mBS, rBS, CBS, phigravmax);

  if (par.verbose)
  {printf("\n==============================================================\n");
  printf("Physical frequency:            = %22.16g\n", omBS);
  printf("Unrescale omega:               = %22.16g\n", omegaini);
  printf("Total mass:                    = %22.16g\n", star.m[par.nint-1]);
  printf("Number of zero crossings:      =     %d\n", zeros);
  printf("=================================================================\n");}

  }

/*==========================================================================*/
/*
Function for computing a BS model using Newton-Raphson scheme on (\omega, \varphi_c) variables. Until a desired threshold is reached, 
the algorithm keeps on iterating/shooting for the BS frequency and the central gravitational scalar field.
*/

void calculate_model(double* rBS, double* rBSE, double* omBS, double* mBS, double* CBS, double* phigravmax, double* omega_init)
{
  double mBSE;

  *omega_init = par.omega0;

  *omBS = rescalePhi(par.omega0);

  // Calculate the Noether charge

  calculateN(*omBS);

  // Calculate the trace of energy-momentum tensor

  calculate_trace_t(*omBS);

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
  *phigravmax = findMax(star.phigrav, 0, par.nint-1);

}

/*==========================================================================*/
/*
Function for gluing the asymptotics to the outward integrated BS solution. This ensures the BS solution has the appropriate
fall-off and mitigates exponentially growing modes.
*/

int calcExtAsym(double omega, double Amp, double phigrav_guess)
  {
  int    i, k, nzero, sigA, istop, imatch;
  double rstop, c1, c2, c3, c4, epsilon, delta, mass;
  double r, X, Phi, A, eta, Psi0, phigrav, rhs_X, rhs_Phi, rhs_A, rhs_phigrav, rhs_Psi0, rhs_eta, om, mphigrav, dr;
  double dX[5], deta[5], dphigrav[5], dPhi[5], dA[5], dPsi0[5];


  // At this stage we have the correct frequency from the shooting
  // algorithm. Here we will compute this model and also remove
  // the diverging part in the profile and replace it with a smooth
  // exterior. Note that we are not yet rescaling time and frequency
  // yet, since we will need the complete profile with exterior to do
  // that.
  intODE(Amp, phigrav_guess, omega, &nzero, &rstop, &sigA);
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
  om = omega * exp(-star.Phi[imatch]);

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

  // if (par.verbose)
  // {
  // printf("\nc1, c2 = %g   %g\n", c1, c2);
  // printf("c3, c4 = %g   %g\n", c3, c4);
  // }

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
/*
Matching criterion for the bosonic scalar field. This is a simple 'smothness' experssion derived from asymptotic expression 
of the bosonic scalar field and its derivative. 
*/

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
/*
Matching criterion for the gravitational scalar field
*/

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
/*
Find largest index i such that r[i] < rtarget.
*/

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
/*
Infrastructure for outward integration using RK4.
*/

void intODE(double A0, double phigrav0, double omega, int* nzeroshooting, double* rstop, int* sigAstop)
  {
  int    i, n1, istop;
  double dr, r, X, A, eta, Phi, Psi0, phigrav, om;
  double rhs_X, rhs_A, rhs_eta, rhs_Phi, rhs_Psi0, rhs_phigrav;
  double dX[5], dA[5], dPsi0[5], deta[5], dphigrav[5], dPhi[5];

  n1 = par.nint;
  om = omega;

  *nzeroshooting = 0;
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
      istop = i - 1;   // We stop the integration; one point as sanity buffer
      // There is one problem here: We regard the present data point and
      // also the previous one as unrealiable. The previous point,
      // however, might have been counted as a zero crossing. If that
      // happened, we wish to undo this counting.
      if(star.A[i-1] * star.A[i-2] < 0) (*nzeroshooting)--;
      break;
      }
    if(star.A[i] * star.A[i-1] < 0) (*nzeroshooting)++;   // Zero crossing found
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
/*
Partial differential equations to be solved for.
*/

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
/*
General expression for the bosonic potential inclusing self-interaction terms (i.e. V(A)).
*/

double V_series(double A)
  {
  // Potential function
  return A*A * (1 + par.lambda4 * A*A + par.lambda6 * A*A*A*A
         + par.lambda8 * A*A*A*A*A*A);
  }

/*==========================================================================*/
/*
Derivative of the bosonic potential w.r.t A^2, i.e. dV/dA^2.
*/

double Vp_series(double A)
  {
  // Potential derviative dV / d(A^2)
  return 1 + 2*par.lambda4 * A*A + 3*par.lambda6 * A*A*A*A
         + 4*par.lambda8 * A*A*A*A*A*A;
  }

/*==========================================================================*/
/*
Solitonic potential.
*/

double V_solitonic(double A)
  {
  // Solitonic potential function
  return A*A * (1 - 2 * A*A / (par.sigma0*par.sigma0))
             * (1 - 2 * A*A / (par.sigma0*par.sigma0));
  }

/*==========================================================================*/
/*
Derivative of the solitonic potential, dV/dA^2.
*/

double Vp_solitonic(double A)
  {
  // Solitonic potential function
  return (1 - 2 * A*A / (par.sigma0*par.sigma0))
       * (1 - 6 * A*A / (par.sigma0*par.sigma0));
  }

/*==========================================================================*/
/*
Damour-Esposito-Farese coupling function, i.e. F(\varphi).
*/

double F_st(double phigrav)
  {
  // Coupling function of the gravitational scalar field
  return exp(-2 * par.alpha0 * phigrav - par.beta0 * phigrav * phigrav);
  }

/*==========================================================================*/
/*
Derivative of the coupling function, i.e. F_{,\varphi} / F.
*/

double derF_st(double phigrav)
  {
  // Function F_{,\varphi} / F
  return -2 * par.alpha0 - 2 * par.beta0 * phigrav;
  }

/*==========================================================================*/
/*
Gravitational scalar potential (simple quadtratic form).
*/

double W_st(double phigrav)
  {
  // Potential for the gravitational scalar field
  return 0.5 * par.mphigrav0 * par.mphigrav0 * phigrav * phigrav;
  }

/*==========================================================================*/
/*
Derivative of the gravitational scalar potential.
*/

double derW_st(double phigrav)
  {
  // Derivative of the potential for the gravitational scalar field
  return par.mphigrav0 * par.mphigrav0 * phigrav;
  }

/*==========================================================================*/
/*
Integrand of the Noether charge.
*/

void integrandN(double *rhs_q, double r, double A, double omega, double X, double phigrav, double Phi)
  {
    double alpha;

    alpha = exp(Phi)/sqrt(F(phigrav));

    //Integrand in the Noether charge formulae
    *rhs_q = 4 * PI * r * r * (A * A * omega * X) / (F(phigrav) * alpha);
  }

/*==========================================================================*/
/*
Integration of Noether charge using RK4.
*/

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
/*
Initialisation of the grid.
*/

void initGrid(int n, double rmax)
  {
  int i;


  for(i = 0; i < n; i++)
    star.r[i] = rmax * i / (n - 1.0);
  }

/*==========================================================================*/
/*
Function for reading the parameter file.
*/

void readPars(char* ifil)
  {
  FILE* ifp;
  char  line[LEN];
  int   n, i;


  // First set all parameters to default values
  par.nint      = 401;
  par.rmax      = 1;
  par.A0        = 0.07;
  par.Phi0        = -1.1;
  par.mphigrav0 = 1.0;
  par.alpha0    = 3.0;
  par.beta0     = 0.0;
  par.phigrav0  = 0.0;
  par.omega0    = 1.;            // omega0 is always 1, it is not specified
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
      else if(strstr(line, "A0") != NULL)
        sscanf(line, "A0 %le", &(par.A0));
      else if(strstr(line, "Phi0") != NULL)
        sscanf(line, "Phi0 %le", &(par.Phi0));
      else if(strstr(line, "mphigrav0") != NULL)
        sscanf(line, "mphigrav0 %le", &(par.mphigrav0));
      else if(strstr(line, "alpha0") != NULL)
        sscanf(line, "alpha0 %le", &(par.alpha0));
      else if(strstr(line, "beta0") != NULL)
        sscanf(line, "beta0 %le", &(par.beta0));
      else if(strstr(line, "phigrav0") != NULL)
        sscanf(line, "phigrav0 %le", &(par.phigrav0));
      else if(strstr(line, "omega0") != NULL)
        sscanf(line, "omega0 %le", &(par.omega0));
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
/*
Print statements of the parameters read from the file, as a sanity check.
*/

void printPars()
  {
  if (par.verbose)
  {
  printf("=======================================\n");
  printf("nint          = %d\n", par.nint);
  printf("rmax          = %g\n", par.rmax);
  printf("A0            = %g\n", par.A0);
  printf("Phi0            = %g\n", par.Phi0);
  printf("mphigrav0     = %g\n", par.mphigrav0);
  printf("phigrav0      = %g\n", par.phigrav0);
  printf("omega0      = %g\n", par.omega0);
  printf("alpha0        = %g\n", par.alpha0);
  printf("beta0         = %g\n", par.beta0);
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
  }

/*==========================================================================*/
/*
What potential function did the user prefer to use.
*/

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
/*
Useful I/O stuff.
*/

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
/*
Function to rescale the BS frequency to the physical one using the lapse function.
*/

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
/*
Calculate the ADM mass of the BS model in the Jordan frame.
*/

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
/*
Calculate the ADM mass of the BS model in the Einstein frame.
*/

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
/*
Calculate the Jordan radius.
*/

void calcRJordan()
  {
  int   i;

  for(i = 1; i <= par.nint; i++)
	  {
      star.rJ[i] = star.r[i]/sqrt(F(star.phigrav[i]));
	  }

  }

/*==========================================================================*/
/*
Calculate the radius of the BS model in the Jordan frame.
*/

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
/*
Calculate the radius of the BS model in the Einstein frame.
*/

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
/*
Find the index of the radius in an array.
*/

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
/*
Calculate the trace of the energy momentum tensor.
*/

void calculate_trace_t(double omega0)
  {
  int   i;
  double A, phigrav, bracket, factor;

  for(i = 0; i <= par.nint; i++)
	  {
      A = star.A[i];
      phigrav = star.phigrav[i];
      bracket = omega0 * omega0 * A * A * exp(-2*star.Phi[i]) * (F(phigrav) * F(phigrav)) 
                - 2 * V(A)
                - (star.Psi0[i] * star.Psi0[i])/(star.X[i] * star.X[i]);
      factor = 1 / F(phigrav);

      star.T[i] = factor * bracket;  
	  }
  }

/*==========================================================================*/
/*
Convert to isotropic gauge.
*/

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
/*
RHS equation for f in isotropic gauge.
*/

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
/*
Simple file with 3 columns
*/

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
/*
File to write BS data collectively
*/

void appendAllData(const char* ofil, double x, double y, double z, double q, double w, double e, double r, double u)
  {
  FILE* ofp;


  ofp = fopen(ofil, "a");
  if(ofp == NULL)
    { printf("Cannot open file   %s\n", ofil); exit(0); }

  fprintf(ofp, "%22.16g   %22.16g   %22.16g   %22.16g   %22.16g   %22.16g   %22.16g    %22.16g\n", x, y, z, q, w, e, r, u);

  fclose(ofp);
  }

/*==========================================================================*/
/*
File to write the final diagnostics of the BS solution.
*/

void out_final_model(double x, double y, double z, double w, double a, double s, double r, int g, double h, double j, const char* ofil)
  {
  FILE* ofp;


  ofp = fopen(ofil, "w");
  if(ofp == NULL) { printf("Cannot open %s in out_final_model\n\n", ofil); exit(0); }

  fprintf(ofp, "# %s\n", ofil);
  fprintf(ofp, "#%22.16s   %22.16s   %22.16s   %22.16s   %22.16s   %22.16s   %22.16s   %22.16s   %22.16s   %22.16s\n", "A0", "phigrav0", "omega", "star mass", "compactness", "radius", "phigrav max", "nzero", "omegain", "Noether");

  fprintf(ofp, "%22.16g   %22.16g   %22.16g   %22.16g   %22.16g   %22.16g   %22.16g                     %d    %22.16g    %22.16g\n", x, y, z, w, a, s, r, g, h, j);

  fclose(ofp);
  }

/*==========================================================================*/
/*
Simple file to write array data from BS calculation.
*/

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
/*
Function to find the maximum value in the array
*/

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
/*
Check whether gravitational scalar field crosses -\alpha_0/beta_0 line.
*/

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
/*
Do we have any NaNs.
*/

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
  // printf("------------------------------------\n");
  // printf("F[.][1] = %15.8g   %15.8g   %15.8g\n",
  //        F[1][1],F[1][2],F[1][3]);
  // printf("J[1][.] = %15.8g   %15.8g   %15.8g\n",
  //        J[1][1],J[1][2],J[1][3]);
  // printf("J[2][.] = %15.8g   %15.8g   %15.8g\n",
  //        J[2][1],J[2][2],J[2][3]);
  // printf("J[3][.] = %15.8g   %15.8g   %15.8g\n",
  //        J[3][1],J[3][2],J[3][3]);
  // printf("------------------------------------\n");

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
