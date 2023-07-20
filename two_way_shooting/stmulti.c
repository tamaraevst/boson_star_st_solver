#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NRANSI
#define LEN             1000
#define NR_END 1
#define FREE_ARG char*
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}


typedef struct exterior
  {
  long double   *y;
  long double   *Phi;
  long double   *m;
  long double   *mJ;
  long double   *zeta;
  long double   *Pi;
  long double   *rho;
  long double   *Lam;
  long double   *r;
  long double   *rJ;
  long double   *X;
  long double   *A;
  long double   *thet;
  long double   *vph;
  long double   *kap;
  long double   *C;
  long double   *CJ;
  long double   *NQ;
  } exterior;

typedef struct interior
  {
  long double   *r;
  long double   *rJ;
  long double   *Phi;
  long double   *X;
  long double   *A;
  long double   *thet;
  long double   *vph;
  long double   *kap;
  long double   *m;
  long double   *mJ;
  long double   *C;
  long double   *CJ;
  long double   *NQ;
  } interior;


typedef struct param
  {
  long double   accP;
  long double   accm;
  long double   accA;
  long double   accz;
  long double   acco;
  long double   accv;
  long double   accr;
  long double   ymatch;
  long double   ymax;
  int           next;
  int           nint;
  char          potential[LEN];
  long double   Phi0;
  long double   mass0;
  long double   zeta1;
  long double   rho1;
  long double   vph0;
  long double   omega0;
  long double   A0start;
  long double   A0end;
  long double   dA0;
  long double   dA0floor;
  long double   lambda4;
  long double   lambda6;
  long double   lambda8;
  long double   sigma0;
  long double   threshold;
  long double   EPS;
  long double   dfactor;
  int           extrapolate;
  int           dataio;
  long double   massvph;
  long double   alp0;
  long double   beta0;
  int           massless;
  int           shoot;
  int           verbose;
  } param;


// Main functions
void allocateExt        ();
void allocateInt        ();
void calcAuxVarsExt     (long double, long double, long double);
void calcAuxVarsInt     (long double, long double, long double);
void calcModel          (long double, long double, long double, long double,
                         long double, long double, long double, int*,
                         long double*, long double*, long double*, long double*,
                         long double*, long double*, int*);
long double calcNoetherCharge(long double);
long double calcRadius  ();
long double calcRadiusJ ();
int checkDeltas         (long double, long double, long double, long double,
                         long double, long double);
void extODE             (long double, long double, long double, long double,
                         int*, int*);
long double findAbsMax  (long double*, long double*, int, int);
long double findMax     (long double*, long double*, int, int);
void initGrid           ();
void intODE             (long double, long double, long double, long double,
                         int*, int*);
int nancheck            (long double, long double, long double, long double,
                         long double, long double);
void oneIter            (long double*, long double*, long double*, long double*,
                         long double*, long double*, long double*, int*, int*);
void printHeader        (FILE*, FILE*);
void printIniGuess      (long double, long double, long double, long double,
                         long double, long double, long double, int, FILE*);
void printInvHeader     ();
void printOutput        (long double, long double, long double, long double,
                         long double, long double, long double, long double,
                         long double, int, int, FILE*, long double, long double,
                         long double);
void printScreen        (long double, long double, long double, long double,
                         long double, long double, long double, long double,
                         long double, int, int);
void printPars          ();
void readPars           (char*);
void registerPotential  ();
void rhsBSext           (long double*, long double*, long double*, long double*,
                         long double*, long double*, long double, long double,
                         long double, long double, long double, long double,
                         long double, long double, long double, long double,
                         long double);
void rhsBSint           (long double*, long double*, long double*, long double*,
                         long double*, long double*, long double, long double,
                         long double, long double, long double, long double,
                         long double, long double);
void shoot              (long double, long double, long double, long double,
                         long double, long double, long double, long double*,
                         long double*, long double*, long double*, long double*,
                         long double*, int*, int*, int*);
void writeProfiles      ();


// Potential and conformal factor functions
long double V_series    (long double);
long double Vhat_series (long double);
long double Vp2_series  (long double);
long double V_solitonic(long double);
long double Vhat_solitonic(long double);
long double Vp2_solitonic(long double);
long double W           (long double);
long double What        (long double);
long double Wp2         (long double);
long double F_DE        (long double);
long double FpoF_DE     (long double);


// Library functions
void    allout1D        (long double*, long double*, long double*, long double*,
                         int, int, int, int,
                         const char*);
long double  **dmatrix  (long, long, long, long);
void    free_dmatrix    (long double**, long, long, long, long);
void    free_ivector    (int*, long, long);
void    gaussj          (long double**, int, long double**, int);
int     *ivector        (long, long);
void    nrerror         (char[]);
int     mygetline       (char*, FILE*);
void    out1D           (long double*, long double*, int, int, const char*);


// Function pointers
long double  (*V)       (long double);
long double  (*Vhat)    (long double);
long double  (*Vp2)     (long double);


// Global variables
long double  pi3;
param        par;
exterior     ext;
interior     Int;


/*==========================================================================*/

int main(int argc, char* argv[])
  {
  long double A0, dA0, Phi0in, Phi0out, mass0in, mass0out, zeta1in, zeta1out,
              omegain, omegaout, rho1in, rho1out, vph0in, vph0out;
  long double A0_p, Phi0_p, mass0_p, zeta1_p, omega_p, vph0_p, rho1_p;
  long double A0_p_p, Phi0_p_p, mass0_p_p, zeta1_p_p, omega_p_p, vph0_p_p,
              rho1_p_p;
  long double A0_start, A0_end;
  long double Cmax, r99, CJmax, rJ99, vphmax, ncharge;
  int         niter, nzero, success, first;
  FILE        *ofp, *ofpini;


  if(argc != 2) { printf("Usage:   %s   <parfile>\n\n", argv[0]); exit(0); }
  pi3 = 2.0L * acosl(0.0L);

  readPars(argv[1]);
  printPars();

  // Register functions and allocate memory
  registerPotential();
  allocateExt();
  allocateInt();

  // Initialize grids
  initGrid();

  // Initialize shooting parameters
  Phi0in  = par.Phi0;
  A0      = par.A0start;
  mass0in = par.mass0;
  zeta1in = par.zeta1;
  omegain = par.omega0;
  vph0in  = par.vph0;
  rho1in  = par.rho1;

  ofp = fopen("output.dat", "w");
  if(ofp == NULL) { printf("Cannot open output.dat in main.\n\n"); exit(0); }

  ofpini = fopen("iniguess.dat", "w");
  if(ofpini == NULL) {printf("Cannot open iniguess.dat in main.\n\n"); exit(0);}


    /*******************************************************/
   /* Compute the first model and stop if this fails      */
  /*******************************************************/
  printHeader(ofp, ofpini);

  success = 0;
  shoot(A0, Phi0in, mass0in, zeta1in, omegain, vph0in, rho1in,
        &Phi0out, &mass0out, &zeta1out, &omegaout, & vph0out, &rho1out,
        &nzero, &success, &niter);

  printIniGuess(A0, Phi0in, mass0in, zeta1in, omegain, vph0in, rho1in, success,
                ofpini);

  // Diagnostics
  Cmax    = findMax(Int.C, ext.C, par.nint, par.next);
  r99     = calcRadius();
  CJmax   = findMax(Int.CJ, ext.CJ, par.nint, par.next);
  rJ99    = calcRadiusJ();
  vphmax  = findAbsMax(Int.vph, ext.vph, par.nint, par.next);
  ncharge = calcNoetherCharge(omegaout);

  // I/O
  if(par.dataio) writeProfiles();
  printScreen(A0, Phi0out, mass0out, zeta1out, omegaout, vph0out, rho1out,
              Cmax, r99, nzero, niter);
  printOutput(A0, Phi0out, mass0out, zeta1out, omegaout, vph0out, rho1out,
              Cmax, r99, nzero, niter, ofp, vphmax, ncharge, rJ99);

  if(success == -1)
    { printf("Failed to compute first model\n\n"); exit(0); }

  // Store the parameters as old versions
  A0_p_p    = A0_p    = A0;
  Phi0_p_p  = Phi0_p  = Phi0out;
  mass0_p_p = mass0_p = mass0out;
  zeta1_p_p = zeta1_p = zeta1out;
  omega_p_p = omega_p = omegaout;
  vph0_p_p  = vph0_p  = vph0out;
  rho1_p_p  = rho1_p  = rho1out;


    /*******************************************************/
   /* Loop over A0 and compute models                     */
  /*******************************************************/
  first = 1;
  dA0      = par.dA0;
  if(par.dA0 == 0.0L) {printf("dA0 = %Lg not a good choice...\n", dA0);exit(0);}

  for(A0 = par.A0start+dA0; (A0 - par.A0end) / dA0 <= 0.0L; A0 += dA0)
    {
    // Initial guesses for the parameters
    if(first)
      {
      // We only have one preceding model which we use
      Phi0in  = Phi0_p;
      mass0in = mass0_p;
      zeta1in = zeta1_p;
      omegain = omega_p;
      vph0in  = vph0_p;
      rho1in  = rho1_p;
      }
    else
      {
      // Extrapolate linearly from the previous 2 models
      Phi0in  = Phi0_p  + dA0 * (Phi0_p  - Phi0_p_p)  / (A0_p - A0_p_p);
      mass0in = mass0_p + dA0 * (mass0_p - mass0_p_p) / (A0_p - A0_p_p);
      zeta1in = zeta1_p + dA0 * (zeta1_p - zeta1_p_p) / (A0_p - A0_p_p);
      omegain = omega_p + dA0 * (omega_p - omega_p_p) / (A0_p - A0_p_p);
      vph0in  = vph0_p  + dA0 * (vph0_p  - vph0_p_p)  / (A0_p - A0_p_p);
      rho1in  = rho1_p  + dA0 * (rho1_p  - rho1_p_p)  / (A0_p - A0_p_p);
      }

    // Compute this model
    success = 0;

    shoot(A0, Phi0in, mass0in, zeta1in, omegain, vph0in, rho1in,
        &Phi0out, &mass0out, &zeta1out, &omegaout, & vph0out, &rho1out,
        &nzero, &success, &niter);

    printIniGuess(A0, Phi0in, mass0in, zeta1in, omegain, vph0in, rho1in,
                  success, ofpini);
    if(success == -1)
      {
      A0  -= dA0;
      dA0 /= par.dfactor;
      printf("Setting dA0 = %Lg\n", dA0);
      if(fabsl(dA0) < par.dA0floor)
        { printf("par.dA0 = %22.16Lg got too small\n\n", dA0); exit(0); }
      continue;
      }

    // If we get here, we have successfully computed a model.

    // Diagnostics
    Cmax   = findMax(Int.C, ext.C, par.nint, par.next);
    r99    = calcRadius();
    CJmax  = findMax(Int.CJ, ext.CJ, par.nint, par.next);
    rJ99   = calcRadiusJ();
    vphmax = findAbsMax(Int.vph, ext.vph, par.nint, par.next);
    ncharge = calcNoetherCharge(omegaout);

    printScreen(A0, Phi0out, mass0out, zeta1out, omegaout, vph0out, rho1out,
                Cmax, r99, nzero, niter);
    printOutput(A0, Phi0out, mass0out, zeta1out, omegaout, vph0out, rho1out,
                Cmax, r99, nzero, niter, ofp, vphmax, ncharge, rJ99);

    // Rotate the parameter sets; this one becomes the new "_p", the
    // "_p" becomes "_p_p".
    A0_p_p    = A0_p;
    Phi0_p_p  = Phi0_p;
    mass0_p_p = mass0_p;
    zeta1_p_p = zeta1_p;
    omega_p_p = omega_p;
    vph0_p_p  = vph0_p;
    rho1_p_p  = rho1_p;
    A0_p      = A0;
    Phi0_p    = Phi0out;
    mass0_p   = mass0out;
    zeta1_p   = zeta1out;
    omega_p   = omegaout;
    vph0_p    = vph0out;
    rho1_p    = rho1out;

    if(par.extrapolate && first) first = 0;
    }

  printInvHeader(ofp, ofpini);

  fclose(ofpini);
  fclose(ofp);
  }

/*==========================================================================*/

void shoot(long double A0, long double Phi0, long double mass0,
           long double zeta1, long double omega, long double vph0,
           long double rho1, long double* Phi0out, long double* mass0out,
           long double* zeta1out, long double* omegaout, long double* vph0out,
           long double* rho1out, int* nzero, int* success, int* niter)
  {
  *niter = 0;

  // Compute the model for the given initial guesses.
  while(*success == 0)
    {
    if(par.verbose) printf("niter = %d\n", *niter);
    oneIter(&Phi0, &A0, &mass0, &zeta1, &omega, &vph0, &rho1, nzero, success);
    (*niter)++;
    }
  // Return the updated parameters for this model
  *Phi0out  = Phi0;
  *mass0out = mass0;
  *zeta1out = zeta1;
  *omegaout = omega;
  *vph0out  = vph0;
  *rho1out  = rho1;
  if(*success == -1) *Phi0out = NAN;
  }

/*==========================================================================*/

void oneIter(long double* Phi0, long double* A0, long double* mass0,
             long double* zeta1, long double* omega, long double* vph0,
             long double* rho1, int* nzero, int* success)
  {
  long double EPS;
  long double DeltaPhi, Deltam, DeltaA, Deltathet, Deltavph, Deltakap;
  long double dPhi0, dmass0, dzeta1, domega, dvph0, drho1;
  long double **J, **F;


  EPS = par.EPS;

  // Allocate Newton-Raphson matrices
  J = dmatrix(1,6,1,6);
  F = dmatrix(1,6,1,1);

  // Base model -- Recall that A0 is fixed and determines the target model
  calcModel(*Phi0, *A0, *mass0, *zeta1, *omega, *vph0, *rho1, nzero, &DeltaPhi,
            &Deltam, &DeltaA, &Deltathet, &Deltavph, &Deltakap, success);

  // If we require not to shoot, we just compute one model and stop
  if(par.shoot == 0) { *success = 1; return; }

  // Check if the model is good enough to stop iterating
  if(par.verbose) printf("%Lg %Lg %Lg %Lg %Lg %Lg\n",
         DeltaPhi, Deltam, DeltaA, Deltathet, Deltavph, Deltakap);
  *success = checkDeltas(DeltaPhi, Deltam, DeltaA, Deltathet,Deltavph,Deltakap);
  if(*success) return;

  // We search for better parameters
  F[1][1] = DeltaPhi;
  F[1][2] = Deltam;
  F[1][3] = DeltaA;
  F[1][4] = Deltathet;
  F[1][5] = Deltavph;
  F[1][6] = Deltakap;

  // 1. vary Phi0
  if(*Phi0 < EPS) dPhi0 = EPS; else dPhi0 = *Phi0 * EPS;
  calcModel(*Phi0+dPhi0, *A0, *mass0, *zeta1, *omega, *vph0, *rho1, nzero,
            &DeltaPhi, &Deltam, &DeltaA, &Deltathet, &Deltavph, &Deltakap,
            success);
  if(*success == -1) return;
  J[1][1] = (DeltaPhi  - F[1][1]) / dPhi0;
  J[2][1] = (Deltam    - F[1][2]) / dPhi0;
  J[3][1] = (DeltaA    - F[1][3]) / dPhi0;
  J[4][1] = (Deltathet - F[1][4]) / dPhi0;
  J[5][1] = (Deltavph  - F[1][5]) / dPhi0;
  J[6][1] = (Deltakap  - F[1][6]) / dPhi0;

  // 2. vary mass0
  if(*mass0 < EPS) dmass0 = EPS; else dmass0 = *mass0 * EPS;
  calcModel(*Phi0, *A0, *mass0+dmass0, *zeta1, *omega, *vph0, *rho1, nzero,
            &DeltaPhi, &Deltam, &DeltaA, &Deltathet, &Deltavph, &Deltakap,
            success);
  if(*success == -1) return;
  J[1][2] = (DeltaPhi  - F[1][1]) / dmass0;
  J[2][2] = (Deltam    - F[1][2]) / dmass0;
  J[3][2] = (DeltaA    - F[1][3]) / dmass0;
  J[4][2] = (Deltathet - F[1][4]) / dmass0;
  J[5][2] = (Deltavph  - F[1][5]) / dmass0;
  J[6][2] = (Deltakap  - F[1][6]) / dmass0;

  // 3. vary zeta1
  if(*zeta1 < EPS) dzeta1 = EPS; else dzeta1 = *zeta1 * EPS;
  calcModel(*Phi0, *A0, *mass0, *zeta1+dzeta1, *omega, *vph0, *rho1, nzero,
            &DeltaPhi, &Deltam, &DeltaA, &Deltathet, &Deltavph, &Deltakap,
            success);
  if(*success == -1) return;
  J[1][3] = (DeltaPhi  - F[1][1]) / dzeta1;
  J[2][3] = (Deltam    - F[1][2]) / dzeta1;
  J[3][3] = (DeltaA    - F[1][3]) / dzeta1;
  J[4][3] = (Deltathet - F[1][4]) / dzeta1;
  J[5][3] = (Deltavph  - F[1][5]) / dzeta1;
  J[6][3] = (Deltakap  - F[1][6]) / dzeta1;

  // 4. vary omega
  if(*omega < EPS) domega = EPS; else domega = *omega * EPS;
  calcModel(*Phi0, *A0, *mass0, *zeta1, *omega+domega, *vph0, *rho1, nzero,
            &DeltaPhi, &Deltam, &DeltaA, &Deltathet, &Deltavph, &Deltakap,
            success);
  if(*success == -1) return;
  J[1][4] = (DeltaPhi  - F[1][1]) / domega;
  J[2][4] = (Deltam    - F[1][2]) / domega;
  J[3][4] = (DeltaA    - F[1][3]) / domega;
  J[4][4] = (Deltathet - F[1][4]) / domega;
  J[5][4] = (Deltavph  - F[1][5]) / domega;
  J[6][4] = (Deltakap  - F[1][6]) / domega;

  // 5. vary vph0
  if(*vph0 < EPS) dvph0 = EPS; else dvph0 = *vph0 * EPS;
  calcModel(*Phi0, *A0, *mass0, *zeta1, *omega, *vph0+dvph0, *rho1, nzero,
            &DeltaPhi, &Deltam, &DeltaA, &Deltathet, &Deltavph, &Deltakap,
            success);
  if(*success == -1) return;
  J[1][5] = (DeltaPhi  - F[1][1]) / dvph0;
  J[2][5] = (Deltam    - F[1][2]) / dvph0;
  J[3][5] = (DeltaA    - F[1][3]) / dvph0;
  J[4][5] = (Deltathet - F[1][4]) / dvph0;
  J[5][5] = (Deltavph  - F[1][5]) / dvph0;
  J[6][5] = (Deltakap  - F[1][6]) / dvph0;

  // 6. vary rho1
  if(*rho1 < EPS) drho1 = EPS; else drho1 = *rho1 * EPS;
  calcModel(*Phi0, *A0, *mass0, *zeta1, *omega, *vph0, *rho1+drho1, nzero,
            &DeltaPhi, &Deltam, &DeltaA, &Deltathet, &Deltavph, &Deltakap,
            success);
  if(*success == -1) return;
  J[1][6] = (DeltaPhi  - F[1][1]) / drho1;
  J[2][6] = (Deltam    - F[1][2]) / drho1;
  J[3][6] = (DeltaA    - F[1][3]) / drho1;
  J[4][6] = (Deltathet - F[1][4]) / drho1;
  J[5][6] = (Deltavph  - F[1][5]) / drho1;
  J[6][6] = (Deltakap  - F[1][6]) / drho1;

  /*printf("J[:][6] = %Lg %Lg %Lg %Lg %Lg %Lg\n",
         J[1][6], J[2][6], J[3][6], J[4][6], J[5][6], J[6][6]);
  exit(0);*/

  // Calculate correction
  // Note: the function gaussj puts the correction in F[1][:]!!!
  gaussj(J, 6, F, 1);

  *Phi0  -= par.accP * F[1][1];
  *mass0 -= par.accm * F[1][2];
  *zeta1 -= par.accz * F[1][3];
  *omega -= par.acco * F[1][4];
  *vph0  -= par.accv * F[1][5];
  *rho1  -= par.accr * F[1][6];

  if(par.verbose)
    {
    //printf("*Phi0  = %Lg\n", *Phi0);
    //printf("*mass0 = %Lg\n", *mass0);
    //printf("*zeta1 = %Lg\n", *zeta1);
    //printf("*omega = %Lg\n", *omega);
    //printf("*vph0  = %Lg\n", *vph0);
    //printf("*rho1  = %Lg\n", *rho1);
    }

  // For omega we need to sanity check that it does not equal 1 or larger.
  // We also force omega to be positive which should not matter since
  // omega appears exclusively in the form omega^2 in the equations.
  if(*omega < 0.0L)  *omega = 0.01L;
  if(*omega >= 1.0L) *omega = 0.99L;

  free_dmatrix(J,1,6,1,6);
  free_dmatrix(F,1,6,1,1);
  }

/*==========================================================================*/

void calcModel(long double Phi0, long double A0, long double mass0,
               long double zeta1, long double om, long double vph0,
               long double rho1, int* nzero, long double* DeltaPhi,
               long double* Deltam, long double* DeltaA,
               long double* Deltathet, long double* Deltavph,
               long double* Deltakap, int* success)
  {
  int nzeroext, nzeroint;


  // One integration, but we now shoot from the outside and inside
  extODE(mass0, zeta1, om, rho1, &nzeroext, success);

  // Extract remaining variables in the exterior
  calcAuxVarsExt(mass0, om, rho1);
  if(*success == -1) return;

  // One integration from r=0 to the matching radius
  intODE(Phi0, A0, om, vph0, &nzeroint, success);

  // Extract some further variables in the interior
  calcAuxVarsInt(mass0, om, rho1);
  if(*success == -1) return;

  *nzero     = nzeroint + nzeroext;
  *DeltaPhi  = ext.Phi[par.next-1]  - Int.Phi[par.nint-1];
  *Deltam    = ext.m[par.next-1]    - Int.m[par.nint-1];
  *DeltaA    = ext.A[par.next-1]    - Int.A[par.nint-1];
  *Deltathet = ext.thet[par.next-1] - Int.thet[par.nint-1];
  *Deltavph  = ext.vph[par.next-1]  - Int.vph[par.nint-1];
  *Deltakap  = ext.kap[par.next-1]  - Int.kap[par.nint-1];
  }

/*==========================================================================*/

void intODE(long double Phi0, long double A0, long double om, long double vph0,
            int* nzero, int* success)
  {
  int         i;
  long double dr, r, Phi, X, A, thet, vph, kap, F;
  long double rhs_Phi, rhs_X, rhs_A, rhs_thet, rhs_vph, rhs_kap;
  long double DPhi[5], DX[5], DA[5], Dthet[5], Dvph[5], Dkap[5];


  // Central values
  F           = F_DE(vph0);
  Int.Phi[0]  = Phi0;
  Int.X[0]    = 1.0L / sqrt(F);
  Int.A[0]    = A0;
  Int.thet[0] = 0.0L;
  Int.vph[0]  = vph0;
  Int.kap[0]  = 0.0L;

  // RK integration
  for(i = 1; i < par.nint; i++)
    {
    dr   = Int.r[i] - Int.r[i-1];

    // 1st RK step
    r    = Int.r[i-1];
    Phi  = Int.Phi[i-1];
    X    = Int.X[i-1];
    A    = Int.A[i-1];
    thet = Int.thet[i-1];
    vph  = Int.vph[i-1];
    kap  = Int.kap[i-1];
    rhsBSint(&rhs_Phi, &rhs_X, &rhs_A, &rhs_thet, &rhs_vph, &rhs_kap,
             r, Phi, X, A, thet, om, vph, kap);
    DPhi[1]  = rhs_Phi  * dr;
    DX[1]    = rhs_X    * dr;
    DA[1]    = rhs_A    * dr;
    Dthet[1] = rhs_thet * dr;
    Dvph[1]  = rhs_vph  * dr;
    Dkap[1]  = rhs_kap  * dr;

    // 2nd RK step
    r    = Int.r[i-1]    + 0.5L * dr;
    Phi  = Int.Phi[i-1]  + 0.5L * DPhi[1];
    X    = Int.X[i-1]    + 0.5L * DX[1];
    A    = Int.A[i-1]    + 0.5L * DA[1];
    thet = Int.thet[i-1] + 0.5L * Dthet[1];
    vph  = Int.vph[i-1]  + 0.5L * Dvph[1];
    kap  = Int.kap[i-1]  + 0.5L * Dkap[1];
    rhsBSint(&rhs_Phi, &rhs_X, &rhs_A, &rhs_thet, &rhs_vph, &rhs_kap,
             r, Phi, X, A, thet, om, vph, kap);
    DPhi[2]  = rhs_Phi  * dr;
    DX[2]    = rhs_X    * dr;
    DA[2]    = rhs_A    * dr;
    Dthet[2] = rhs_thet * dr;
    Dvph[2]  = rhs_vph  * dr;
    Dkap[2]  = rhs_kap  * dr;

    // 3rd RK step
    r    = Int.r[i-1]    + 0.5L * dr;
    Phi  = Int.Phi[i-1]  + 0.5L * DPhi[2];
    X    = Int.X[i-1]    + 0.5L * DX[2];
    A    = Int.A[i-1]    + 0.5L * DA[2];
    thet = Int.thet[i-1] + 0.5L * Dthet[2];
    vph  = Int.vph[i-1]  + 0.5L * Dvph[2];
    kap  = Int.kap[i-1]  + 0.5L * Dkap[2];
    rhsBSint(&rhs_Phi, &rhs_X, &rhs_A, &rhs_thet, &rhs_vph, &rhs_kap,
             r, Phi, X, A, thet, om, vph, kap);
    DPhi[3]  = rhs_Phi  * dr;
    DX[3]    = rhs_X    * dr;
    DA[3]    = rhs_A    * dr;
    Dthet[3] = rhs_thet * dr;
    Dvph[3]  = rhs_vph  * dr;
    Dkap[3]  = rhs_kap  * dr;

    // 4th RK step
    r    = Int.r[i];
    Phi  = Int.Phi[i-1]  + DPhi[3];
    X    = Int.X[i-1]    + DX[3];
    A    = Int.A[i-1]    + DA[3];
    thet = Int.thet[i-1] + Dthet[3];
    vph  = Int.vph[i-1]  + Dvph[3];
    kap  = Int.kap[i-1]  + Dkap[3];
    rhsBSint(&rhs_Phi, &rhs_X, &rhs_A, &rhs_thet, &rhs_vph, &rhs_kap,
             r, Phi, X, A, thet, om, vph, kap);
    DPhi[4]  = rhs_Phi  * dr;
    DX[4]    = rhs_X    * dr;
    DA[4]    = rhs_A    * dr;
    Dthet[4] = rhs_thet * dr;
    Dvph[4]  = rhs_vph  * dr;
    Dkap[4]  = rhs_kap  * dr;

    // Update variable
    Int.Phi[i]  = Int.Phi[i-1]
                  + (DPhi[1]  + 2.0L*DPhi[2]  + 2.0L*DPhi[3]  + DPhi[4] )/6.0L;
    Int.X[i]    = Int.X[i-1]
                  + (DX[1]    + 2.0L*DX[2]    + 2.0L*DX[3]    + DX[4]   )/6.0L;
    Int.A[i]    = Int.A[i-1]
                  + (DA[1]    + 2.0L*DA[2]    + 2.0L*DA[3]    + DA[4]   )/6.0L;
    Int.thet[i] = Int.thet[i-1]
                  + (Dthet[1] + 2.0L*Dthet[2] + 2.0L*Dthet[3] + Dthet[4])/6.0L;
    Int.vph[i]  = Int.vph[i-1]
                  + (Dvph[1]  + 2.0L*Dvph[2]  + 2.0L*Dvph[3]  + Dvph[4] )/6.0L;
    Int.kap[i]  = Int.kap[i-1]
                  + (Dkap[1]  + 2.0L*Dkap[2]  + 2.0L*Dkap[3]  + Dkap[4] )/6.0L;

    // Check for nan; we then flag failure and return
    if( nancheck(Int.Phi[i], Int.X[i], Int.A[i], Int.thet[i], Int.vph[i],
                 Int.kap[i]) )
      {
      printf("NaN in intODE at r[%d] = %15.8Lg\n", i, Int.r[i]);
      printf("Phi,X,A,thet,vph,kap = %Lg %Lg %Lg %Lg %Lg %Lg\n", Int.Phi[i],
             Int.X[i], Int.A[i], Int.thet[i], Int.vph[i], Int.kap[i]);
      *success = - 1;
      return;
      }

    // Check for zero crossings
    if(Int.A[i-1] * Int.A[i] < 0.0L)
      (*nzero)++;
    }
  }

/*==========================================================================*/

void rhsBSint(long double* rhs_Phi, long double* rhs_X, long double* rhs_A,
              long double* rhs_thet, long double* rhs_vph, long double* rhs_kap,
              long double r, long double Phi, long double X, long double A,
              long double thet, long double om, long double vph,
              long double kap)
  {
  long double F, FpoF, alp;


  if(r < 1.0e-15L)
    {
    // We are at the origin, i.e. r=0, and need the asymptotic behaviour:
    // For thet, kap we use that near r=0, we have
    //
    //       kap(r) = kap(0) + kap'(0)*r + ... = kap'(0)*r + ...
    //  ==>  2*kap/r = 2*kap'(0) + ...
    //
    //  This gives a factor 1/3 to be applied to the remaining rhs
    //  of the kap equation; likewise for thet.
    F    = F_DE(vph);
    FpoF = FpoF_DE(vph);
    alp  = expl(Phi) / sqrtl(F);

    *rhs_Phi  = 0.0L;
    *rhs_X    = 0.0L;
    *rhs_A    = 0.0L;
    *rhs_thet = X / (F * alp) * A * (alp*alp * Vp2(A) - om*om) / 3.0L;
    *rhs_vph  = 0.0L;
    *rhs_kap  = ( 2.0L * alp * X * F * Wp2(vph) * vph + 2.0L*pi3 * X / alp *
                  FpoF / F * (om*om * A*A - 2.0L * alp*alp * V(A)
                  - F*F * thet*thet) ) / 3.0L;

    /*printf("rhs_Phi  = %Lg\n", *rhs_Phi);
    printf("rhs_X    = %Lg\n", *rhs_X);
    printf("rhs_A    = %Lg\n", *rhs_A);
    printf("rhs_thet = %Lg\n", *rhs_thet);
    printf("rhs_vph  = %Lg\n", *rhs_vph);
    printf("rhs_kap  = %Lg\n", *rhs_kap);*/
    }
  else
    {
    F    = F_DE(vph);
    FpoF = FpoF_DE(vph);
    alp  = expl(Phi) / sqrtl(F);

    *rhs_Phi  = 0.5L * (F * X*X - 1.0L) / r - r * F * X*X * W(vph)
                + 0.5L * r * X*X / (alp*alp) * kap*kap
                + 2.0L*pi3 * r * X*X / (F * alp*alp) *
                  ( om*om * A*A - alp*alp * V(A) + F*F * thet*thet );
    *rhs_X    = X * ( -0.5L * (F * X*X - 1.0L) / r + r * F * X*X * W(vph)
                -0.5L * X / alp * FpoF * kap
                + 0.5L * r * X*X / (alp*alp) * kap*kap
                + 2.0L*pi3 * r * X*X / (F * alp*alp) *
                  ( om*om * A*A + alp*alp * V(A) + F*F * thet*thet ) );
    *rhs_A    = X * F * thet / alp;
    *rhs_thet = -2.0L * thet / r + X / (alp * F) * A *
                (alp*alp * Vp2(A) - om*om);
    *rhs_vph  = X * kap / alp;
    *rhs_kap  = -2.0L * kap / r + 2.0L * alp * X * F * Wp2(vph) * vph
                + 2.0L*pi3 * X * FpoF / (alp * F) *
                  ( om*om * A*A - 2.0L * alp*alp * V(A) - F*F * thet*thet );

    /*printf("rhs_Phi  = %Lg\n", *rhs_Phi);
    printf("rhs_X    = %Lg\n", *rhs_X);
    printf("rhs_A    = %Lg\n", *rhs_A);
    printf("rhs_thet = %Lg\n", *rhs_thet);
    printf("rhs_vph  = %Lg\n", *rhs_vph);
    printf("rhs_kap  = %Lg\n", *rhs_kap);
    exit(0);*/
    }
  }

/*==========================================================================*/

void calcAuxVarsInt(long double mass, long double om, long double rho1)
  {
  int   i;
  long double r, X, F, FpoF, alp, Q;


  Int.m[0]  = 0.0L;
  Int.C[0]  = 0.0L;
  Int.rJ[0] = 0.0L;
  Int.mJ[0] = 0.0L;
  Int.CJ[0] = 0.0L;

  for(i = 1; i < par.nint; i++)
    {
    r    = Int.r[i];
    X    = Int.X[i];
    F    = F_DE(Int.vph[i]);
    FpoF = FpoF_DE(Int.vph[i]);
    alp  = expl(Int.Phi[i]) / sqrtl(F);
    Q    = (1.0L - 0.5L*r * FpoF * X / alp * Int.kap[i]) / sqrtl(F);
    Int.m[i]  = 0.5L * r * (1.0L - 1.0L / (F*X*X));
    Int.C[i]  = Int.m[i] / r;
    Int.rJ[i] = r / sqrtl(F);
    Int.mJ[i] = 0.5L * Int.rJ[i] * (1.0L - Q*Q / (X*X));
    Int.CJ[i] = Int.mJ[i] / Int.rJ[i];
    }
  }

/*==========================================================================*/

void extODE(long double mass0, long double zeta1, long double om,
            long double rho1, int* nzero, int* success)
  {
  int         i, n;
  long double dy, y, k, h;
  long double Phi, m, zeta, Pi, rho, Lam;
  long double rhs_Phi, rhs_m, rhs_zeta, rhs_Pi, rhs_rho, rhs_Lam;
  long double DPhi[5], Dm[5], Dzeta[5], DPi[5], Drho[5], DLam[5];


  n      = par.next;
  *nzero = 0;

  // Values at inifinity
  ext.Phi[0]  = 0.0L;
  ext.m[0]    = mass0;
  ext.zeta[0] = 0.0L;
  ext.Pi[0]   = 0.0L;
  ext.rho[0]  = 0.0L;
  ext.Lam[0]  = 0.0L;

  // RK integration
  for(i = 1; i < n; i++)
    {
    dy = ext.y[i] - ext.y[i-1];
    y  = ext.y[i-1];

    // Asymptotic expansion or numerical integration?
    if(ext.y[i] < par.ymatch)
      {
      k  = sqrt(1.0L - om*om);
      h  = par.massvph;
      y  = ext.y[i];

      ext.Phi[i]  = -mass0 * y;
      ext.zeta[i] = zeta1 * y;
      ext.Pi[i]   = k * zeta1 * y;
      ext.rho[i]  = rho1 * y;

      // The variables m, Lambda have different asymptotics in the
      // case of massless or massive grav.scalar.
      if( par.massless )
        {
        ext.m[i]    = mass0 - 0.5L * rho1*rho1 * y;
        ext.Lam[i]  = rho1 * y*y;
        }
      else
        {
        ext.m[i]    = mass0;
        ext.Lam[i]  = h * rho1 * y + ( (1.0L - mass0 * h / 2.0L)
                      + 1.5L * mass0*mass0 * h*h ) * rho1 * y*y;
        }
      }
    else
      {
      // Now the regular integration. 1st RK step
      y    = ext.y[i-1];
      Phi  = ext.Phi[i-1];
      m    = ext.m[i-1];
      zeta = ext.zeta[i-1];
      Pi   = ext.Pi[i-1];
      rho  = ext.rho[i-1];
      Lam  = ext.Lam[i-1];
      rhsBSext(&rhs_Phi, &rhs_m, &rhs_zeta, &rhs_Pi, &rhs_rho, &rhs_Lam,
               y, Phi, m, zeta, Pi, rho, Lam, mass0, zeta1, om, rho1);
      DPhi[1]  = rhs_Phi  * dy;
      Dm[1]    = rhs_m    * dy;
      Dzeta[1] = rhs_zeta * dy;
      DPi[1]   = rhs_Pi   * dy;
      Drho[1]  = rhs_rho  * dy;
      DLam[1]  = rhs_Lam  * dy;

      // 2nd RK step
      y    = ext.y[i-1]    + 0.5L * dy;
      Phi  = ext.Phi[i-1]  + 0.5L * DPhi[1];
      m    = ext.m[i-1]    + 0.5L * Dm[1];
      zeta = ext.zeta[i-1] + 0.5L * Dzeta[1];
      Pi   = ext.Pi[i-1]   + 0.5L * DPi[1];
      rho  = ext.rho[i-1]  + 0.5L * Drho[1];
      Lam  = ext.Lam[i-1]  + 0.5L * DLam[1];
      rhsBSext(&rhs_Phi, &rhs_m, &rhs_zeta, &rhs_Pi, &rhs_rho, &rhs_Lam,
               y, Phi, m, zeta, Pi, rho, Lam, mass0, zeta1, om, rho1);
      DPhi[2]  = rhs_Phi  * dy;
      Dm[2]    = rhs_m    * dy;
      Dzeta[2] = rhs_zeta * dy;
      DPi[2]   = rhs_Pi   * dy;
      Drho[2]  = rhs_rho  * dy;
      DLam[2]  = rhs_Lam  * dy;

      // 3rd RK step
      y    = ext.y[i-1]    + 0.5L * dy;
      Phi  = ext.Phi[i-1]  + 0.5L * DPhi[2];
      m    = ext.m[i-1]    + 0.5L * Dm[2];
      zeta = ext.zeta[i-1] + 0.5L * Dzeta[2];
      Pi   = ext.Pi[i-1]   + 0.5L * DPi[2];
      rho  = ext.rho[i-1]  + 0.5L * Drho[2];
      Lam  = ext.Lam[i-1]  + 0.5L * DLam[2];
      rhsBSext(&rhs_Phi, &rhs_m, &rhs_zeta, &rhs_Pi, &rhs_rho, &rhs_Lam,
               y, Phi, m, zeta, Pi, rho, Lam, mass0, zeta1, om, rho1);
      DPhi[3]  = rhs_Phi  * dy;
      Dm[3]    = rhs_m    * dy;
      Dzeta[3] = rhs_zeta * dy;
      DPi[3]   = rhs_Pi   * dy;
      Drho[3]  = rhs_rho  * dy;
      DLam[3]  = rhs_Lam  * dy;

      // 4th RK step
      y    = ext.y[i];
      Phi  = ext.Phi[i-1]  + DPhi[3];
      m    = ext.m[i-1]    + Dm[3];
      zeta = ext.zeta[i-1] + Dzeta[3];
      Pi   = ext.Pi[i-1]   + DPi[3];
      rho  = ext.rho[i-1]  + Drho[3];
      Lam  = ext.Lam[i-1]  + DLam[3];
      rhsBSext(&rhs_Phi, &rhs_m, &rhs_zeta, &rhs_Pi, &rhs_rho, &rhs_Lam,
               y, Phi, m, zeta, Pi, rho, Lam, mass0, zeta1, om, rho1);
      DPhi[4]  = rhs_Phi  * dy;
      Dm[4]    = rhs_m    * dy;
      Dzeta[4] = rhs_zeta * dy;
      DPi[4]   = rhs_Pi   * dy;
      Drho[4]  = rhs_rho  * dy;
      DLam[4]  = rhs_Lam  * dy;

      // Update variables
      ext.Phi[i]  = ext.Phi[i-1] +
                    (DPhi[1]  + 2.0L*DPhi[2]  + 2.0L*DPhi[3]  + DPhi[4] )/6.0L;
      ext.m[i]    = ext.m[i-1] +
                    (Dm[1]    + 2.0L*Dm[2]    + 2.0L*Dm[3]    + Dm[4]   )/6.0L;
      ext.zeta[i] = ext.zeta[i-1] +
                    (Dzeta[1] + 2.0L*Dzeta[2] + 2.0L*Dzeta[3] + Dzeta[4])/6.0L;
      ext.Pi[i]   = ext.Pi[i-1] +
                    (DPi[1]   + 2.0L*DPi[2]   + 2.0L*DPi[3]   + DPi[4]  )/6.0L;
      ext.rho[i]  = ext.rho[i-1] +
                    (Drho[1]  + 2.0L*Drho[2]  + 2.0L*Drho[3]  + Drho[4] )/6.0L;
      ext.Lam[i]  = ext.Lam[i-1] +
                    (DLam[1]  + 2.0L*DLam[2]  + 2.0L*DLam[3]  + DLam[4] )/6.0L;

      // Check for nan; we then flag failure and return
      if( nancheck(ext.Phi[i], ext.m[i], ext.zeta[i], ext.Pi[i], ext.rho[i],
                   ext.Lam[i]) )
        {
        printf("NaN in extODE at y[%d] = %15.8Lg\n", i, ext.y[i]);
        printf("Phi,m,zeta,Pi,rho,Lam = %Lg %Lg %Lg %Lg %Lg %Lg\n", ext.Phi[i],
               ext.m[i], ext.zeta[i], ext.Pi[i], ext.rho[i], ext.Lam[i]);
        *success = - 1;
        return;
        }

      // Analyze; do we cross zero?
      if(ext.zeta[i] * ext.zeta[i-1] < 0.0L || ext.zeta[i-1] == 0.0L)
        (*nzero)++;
      }
    }
  }

/*==========================================================================*/

void rhsBSext(long double* rhs_Phi, long double* rhs_m, long double* rhs_zeta,
              long double* rhs_Pi, long double* rhs_rho, long double* rhs_Lam,
              long double y, long double Phi, long double m, long double zeta,
              long double Pi, long double rho, long double Lam,
              long double mass0, long double zeta1, long double om,
              long double rho1)
  {
  long double k, h, eps, delta, eky, e2ky, ehy, e2hy, A, vph, F, FpoF, X2, X,
              alp;


  if(y < 1.0e-15L)
    {
    // We are at infinity, i.e. y=0, and need the asymptotic behaviour
    // to integrate out of the singularity. The boundary conditions are
    //
    //   Phi,y  = -M
    //   m,y    = 0 (massive)   or   -rho1^2/2 (massless)
    //   zeta,y = zeta1
    //   Pi,y   = k * zeta1
    //   rho,y  = rho1
    //   Lam,y  = h * rho1 (massive)   or   0 (massless)
    //
    k = sqrt(1.0L - om*om);
    h = par.massvph;

    *rhs_Phi  = -mass0;
    if( par.massless )
      *rhs_m  = -0.5L * rho1*rho1;
    else
      *rhs_m  = 0;
    *rhs_zeta = zeta1;
    *rhs_Pi   = zeta1 * k;
    *rhs_rho  = rho1;
    *rhs_Lam  = h * rho1;
    }
  else
    {
    k     = sqrt(1.0L - om*om);
    h     = par.massvph;
    delta = mass0 * h;
    eky   = expl(k/y);
    e2ky  = eky*eky;
    if(par.massless) eps = ( mass0 * (2.0L*k*k - 1.0L) + par.alp0 * rho1 ) / k;
    else             eps =   mass0 * (2.0L*k*k - 1.0L) / k;
    ehy   = expl(h/y);
    e2hy  = ehy*ehy;
    A     = powl(y,eps) * zeta / eky;
    vph   = powl(y,delta) * rho / ehy;
    F     = F_DE(vph);
    FpoF  = FpoF_DE(vph);
    X2    = 1.0L / ( F * (1.0L - 2.0L*m*y) );
    X     = sqrtl(X2);
    alp   = expl(Phi) / sqrtl(F);

    *rhs_Phi  = -0.5L * (F * X2 - 1.0L) / y
                + powl(y,2.0L*delta-3.0L) / e2hy * X2 *
                  ( F * What(vph) * rho*rho - 0.5L * Lam*Lam / (alp*alp) )
                + 2.0L*pi3 * powl(y,2.0L*eps-3.0L) * X2 / (F * alp*alp) / e2ky *
                  ( (alp*alp * Vhat(A) - om*om) * zeta*zeta - F*F * Pi*Pi );
    *rhs_m    = -powl(y,2.0L*delta-4.0L) / e2hy *
                  ( What(vph) * rho*rho + 0.5L * Lam*Lam / (F * alp*alp) )
                - 2.0L*pi3 * powl(y,2.0L*eps-4.0L) / (F*F * alp*alp * e2ky) *
                  ( (alp*alp * Vhat(A) + om*om) * zeta*zeta + F*F * Pi*Pi );
    *rhs_rho  = (X * Lam / alp - (y * delta + h) * rho) / (y*y);
    *rhs_Lam  = ((2.0L - delta) * y - h) * Lam / (y*y)
                + 2.0L * alp * X * F / (y*y) * Wp2(vph) * rho
                + 2.0L*pi3 / powl(y,2.0L+delta-2.0L*eps) * ehy / e2ky *
                  X / alp * FpoF / F *
                  ( (om*om - 2.0L*alp*alp*Vhat(A)) * zeta*zeta - F*F * Pi*Pi );
    *rhs_zeta = (F * X * Pi / alp - (eps*y + k) * zeta)/ (y*y);
    *rhs_Pi   = ((2.0L - eps) * y - k) * Pi / (y*y)
                + X * zeta / (alp * F * y*y) * (alp*alp * Vp2(A) - om*om);
    }

    // Phi version (instead of alpha); we leave that for now
    /*
    *rhs_Phi  = -0.5L * (F * X2 - 1.0L) / y
                + powl(y,2.0L*delta-3.0L) / e2hy * X2 * F
                  * ( What(vph) * rho*rho - 0.5L * Lam*Lam * expl(-2.0L*Phi) )
                + 2.0L*pi3*X2 * powl(y,2.0L*eps-3.0L) / e2ky
                * ( zeta*zeta * Vhat(A) / F - expl(-2.0L*Phi) * (F*F * Pi*Pi
                    + zeta*zeta * om*om) );
    *rhs_m    = -powl(y,2.0L*delta-4.0L) / e2hy * ( What(vph) * rho*rho
                + 0.5L * Lam*Lam * exp(
                -2.0L*pi3 / (powl(y,4.0L-2.0L*eps) * e2ky) *
                ( Pi*Pi + zeta*zeta * (om*om * expl(-2.0L*Phi) + Vhat(A)) );
    */
  }

/*==========================================================================*/

void calcAuxVarsExt(long double mass0, long double om, long double rho1)
  {
  int         i;
  long double y, m, zeta, Pi, rho, Lam;
  long double k, eky, eps, h, ehy, delta, F, FpoF, vph, X2, alp, Q;


  for(i = 0; i < par.next; i++)
    {
    if(ext.y[i] <= 0.0L)
      {
      ext.r[i]    = 2.0L / ext.y[i+1];   // Use large but finite r for infinity
      ext.X[i]    = 1.0L;
      ext.C[i]    = 0.0L;
      ext.A[i]    = 0.0L;
      ext.thet[i] = 0.0L;
      ext.vph[i]  = 0.0L;
      ext.kap[i]  = 0.0L;
      ext.rJ[i]   = ext.r[i];
      ext.mJ[i]   = ext.m[i];
      ext.CJ[i]   = 0.0L;
      }
    else
      {
      y    = ext.y[i];
      m    = ext.m[i];
      zeta = ext.zeta[i];
      Pi   = ext.Pi[i];
      rho  = ext.rho[i];
      Lam  = ext.Lam[i];

      k     = sqrt(1.0L - om*om);
      eky   = expl(k/y);
      if(par.massless) eps = (mass0 * (2.0L*k*k - 1.0L) + par.alp0 * rho1) / k;
      else             eps =  mass0 * (2.0L*k*k - 1.0L) / k;
      h     = par.massvph;
      ehy   = expl(h/y);
      delta = mass0 * h;

      vph   = powl(y,delta) * rho / ehy;
      F     = F_DE(vph);
      FpoF  = FpoF_DE(vph);
      alp   = expl(ext.Phi[i]) / sqrtl(F);
      X2    = 1.0L / ( F * (1.0L - 2.0L*m*y) );
      Q     = (1.0L + 0.5L*powl(y,delta-1.0L) * FpoF * sqrtl(X2)/alp/ehy * Lam)
              / sqrtl(F);

      ext.r[i]    = 1.0L / ext.y[i];
      ext.C[i]    = m * y;
      ext.X[i]    = sqrtl(X2);
      ext.vph[i]  = vph;
      ext.kap[i]  = -powl(y,delta) * Lam / ehy;
      ext.A[i]    =  powl(y,eps) * zeta / eky;
      ext.thet[i] = -powl(y,eps) * Pi   / eky;

      ext.rJ[i]   = ext.r[i] / sqrtl(F);
      ext.mJ[i]   = 0.5L * ext.rJ[i] * (1.0L - Q*Q / X2);
      ext.CJ[i]   = ext.mJ[i] / ext.rJ[i];
      }
    }
  }

/*==========================================================================*/

long double findAbsMax(long double* y1, long double* y2, int n1, int n2)
  {
  int         i;
  long double max;


  // Find maximum over two input arrays
  max = y1[0];

  for(i = 1; i < n1; i++)
    if(fabsl(y1[i]) > fabsl(max)) max = y1[i];

  for(i = 0; i < n2; i++)
    if(fabsl(y2[i]) > fabsl(max)) max = y2[i];

  return max;
  }

/*==========================================================================*/

long double findMax(long double* y1, long double* y2, int n1, int n2)
  {
  int         i;
  long double max;


  // Find maximum over two input arrays
  max = y1[0];

  for(i = 1; i < n1; i++)
    if(y1[i] > max) max = y1[i];

  for(i = 0; i < n2; i++)
    if(y2[i] > max) max = y2[i];

  return max;
  }

/*==========================================================================*/

long double calcRadius()
  {
  int i;


  for(i = 1; i < par.next; i++)
    if(ext.m[i] < 0.99L * ext.m[0])
      return 2.0L / (ext.y[i] + ext.y[i-1]);

  for(i = par.nint-2; i >= 0; i--)
    if(Int.m[i] < 0.99L * ext.m[0])
      return 0.5L * (Int.r[i] + Int.r[i+1]);

  // Return an unphysical value if we haven't found a radius
  return -1;
  }

/*==========================================================================*/

long double calcRadiusJ()
  {
  int i;


  for(i = 1; i < par.next; i++)
    if(ext.mJ[i] < 0.99L * ext.mJ[0])
      return 0.5L * (ext.rJ[i] + ext.rJ[i-1]);

  for(i = par.nint-2; i >= 0; i--)
    if(Int.mJ[i] < 0.99L * ext.mJ[0])
      return 0.5L * (Int.rJ[i] + Int.rJ[i+1]);

  // Return an unphysical value if we haven't found a radius
  return -1;
  }

/*==========================================================================*/

long double calcNoetherCharge(long double omega)
  {
  int i;
  long double r, dr, Phi, X, A, vph, F, alp;


  Int.NQ[0] = 0.0L;

  // Interior contribution
  for(i = 1; i < par.nint; i++)
    {
    r   = 0.5L * (Int.r[i] + Int.r[i-1]);
    dr  = Int.r[i] - Int.r[i-1];
    Phi = 0.5L * (Int.Phi[i] + Int.Phi[i-1]);
    X   = 0.5L * (Int.X[i]   + Int.X[i-1]  );
    A   = 0.5L * (Int.A[i]   + Int.A[i-1]  );
    vph = 0.5L * (Int.vph[i] + Int.vph[i-1]);
    F   = F_DE(vph);
    alp = expl(Phi) / sqrtl(F);

    Int.NQ[i] = Int.NQ[i-1] + 4.0L*pi3 * X * r*r * omega * A*A / (F*alp) * dr;
    }

  // At the matching radius...
  ext.NQ[par.next-1] = Int.NQ[par.nint-1];

  // Exterior contribution
  for(i = par.next-1; i > 0; i--)
    {
    r   = 0.5L * (ext.r[i-1] + ext.r[i]);
    dr  = ext.r[i-1] - ext.r[i];
    Phi = 0.5L * (ext.Phi[i-1] + ext.Phi[i]);
    X   = 0.5L * (ext.X[i-1]   + ext.X[i]  );
    A   = 0.5L * (ext.A[i-1]   + ext.A[i]  );
    vph = 0.5L * (ext.vph[i-1] + ext.vph[i]);
    F   = F_DE(vph);
    alp = expl(Phi) / sqrtl(F);

    ext.NQ[i-1] = ext.NQ[i] + 4.0L*pi3 * X * r*r * omega * A*A / (F*alp) * dr;
    }

  return ext.NQ[0];
  }

/*==========================================================================*/

void printScreen(long double A0, long double Phi0out, long double mass0out,
                 long double zeta1out, long double omegaout,
                 long double vph0out, long double rho1out, long double Cmax,
                 long double r99, int nzero, int niter)
  {
  printf("%11.9Lg", A0);
  printf(" %13.8Lg", Phi0out);
  printf(" %13.8Lg", mass0out);
  printf(" %13.8Lg", zeta1out);
  printf(" %13.8Lg", omegaout);
  printf(" %13.8Lg", vph0out);
  printf(" %13.8Lg", rho1out);
  printf(" %11.8Lg",  Cmax);
  printf(" %11.8Lg",  r99);
  printf(" %4d",      nzero);
  printf(" %5d\n",    niter);

  fflush(stdout);
  }

/*==========================================================================*/

void printOutput(long double A0, long double Phi0out, long double mass0out,
                 long double zeta1out, long double omegaout,
                 long double vph0out, long double rho1out, long double Cmax,
                 long double r99, int nzero, int niter, FILE* ofp,
                 long double vphmax, long double ncharge, long double rJ99)
  {
  fprintf(ofp, "%12.10Lg", A0);
  fprintf(ofp, " %20.17Lg", Phi0out);
  fprintf(ofp, " %20.17Lg", mass0out);
  fprintf(ofp, " %20.17Lg", zeta1out);
  fprintf(ofp, " %20.17Lg", omegaout);
  fprintf(ofp, " %20.17Lg", vph0out);
  fprintf(ofp, " %20.17Lg", rho1out);
  fprintf(ofp, " %8.6Lg",  Cmax);
  fprintf(ofp, " %7.5Lg",  r99);
  fprintf(ofp, " %3d",      nzero);
  fprintf(ofp, " %5d",      niter);
  fprintf(ofp, " %9.7Lg",  vphmax);
  fprintf(ofp, " %8.6Lg",  ncharge);
  fprintf(ofp, " %7.5Lg\n",  rJ99);

  fflush(ofp);
  }

/*==========================================================================*/

void printIniGuess(long double A0, long double Phi0in, long double mass0in,
                   long double zeta1in, long double omegain, long double vph0in,
                   long double rho1in, int success, FILE* ofpini)
  {
  fprintf(ofpini, "%13.11Lg", A0);
  fprintf(ofpini, " %29.25Lg", Phi0in);
  fprintf(ofpini, " %29.25Lg", mass0in);
  fprintf(ofpini, " %29.25Lg", zeta1in);
  fprintf(ofpini, " %29.25Lg", omegain);
  fprintf(ofpini, " %29.25Lg", vph0in);
  fprintf(ofpini, " %29.25Lg", rho1in);
  fprintf(ofpini, " %4d\n", success);

  fflush(ofpini);
  }

/*==========================================================================*/

void printHeader(FILE* ofp, FILE* ofpini)
  {
  printf("\n");
  printf("     A0       ");
  printf("     Phi0     ");
  printf("     mass0    ");
  printf("     zeta1    ");
  printf("     omega    ");
  printf("     vph0     ");
  printf("     rho1     ");
  printf("     Cmax           r99    nzero  niter\n");
  printf("----------------------------------------");
  printf("----------------------------------------");
  printf("----------------------------------------");
  printf("-----------------\n");

  fprintf(ofp, "#      A0    ");
  fprintf(ofp, "          Phi0       ");
  fprintf(ofp, "         mass0         ");
  fprintf(ofp, "         zeta1        ");
  fprintf(ofp, "        omega         ");
  fprintf(ofp, "         vph0          ");
  fprintf(ofp, "         rho1          ");
  fprintf(ofp, " Cmax      r99  nzero niter  vphmax   ncharge    rJ99\n");
  fprintf(ofp, "#----------------------------------------");
  fprintf(ofp, "----------------------------------------");
  fprintf(ofp, "----------------------------------------");
  fprintf(ofp, "----------------------------------------");
  fprintf(ofp, "------------------------------------------\n");

  fprintf(ofpini, "#      A0    ");
  fprintf(ofpini, "             Phi0              ");
  fprintf(ofpini, "             mass0             ");
  fprintf(ofpini, "             zeta1             ");
  fprintf(ofpini, "             omega             ");
  fprintf(ofpini, "             vph0              ");
  fprintf(ofpini, "             rho1           ");
  fprintf(ofpini, "success\n");
  fprintf(ofpini, "#----------------------------------------");
  fprintf(ofpini, "----------------------------------------");
  fprintf(ofpini, "----------------------------------------");
  fprintf(ofpini, "----------------------------------------");
  fprintf(ofpini, "------------------------------------------\n");
  }

/*==========================================================================*/

void printInvHeader()
  {
  printf("----------------------------------------");
  printf("----------------------------------------");
  printf("----------------------------------------");
  printf("-----------------\n");
  printf("     A0       ");
  printf("     Phi0     ");
  printf("     mass0    ");
  printf("     zeta1    ");
  printf("     omega    ");
  printf("     vph0     ");
  printf("     rho1     ");
  printf("     Cmax           r99    nzero  niter\n");
  }

/*==========================================================================*/

void initGrid()
  {
  int         i;
  long double rmax;


  // Exterior grid
  for(i = 0; i < par.next; i++)
    ext.y[i] = par.ymax * i / (par.next - 1.0L);

  // Interior grid
  if(par.ymax <= 0.0L)
    { printf("ymax = %Lg not allowed\n\n", par.ymax); exit(0); }
  else
    rmax = 1.0L / par.ymax;

  for(i = par.nint-1; i >= 0; i--)
    Int.r[i] = rmax * i /(par.nint - 1.0L);
  }

/*==========================================================================*/

int checkDeltas(long double DeltaPhi, long double Deltam, long double DeltaA,
                long double Deltathet, long double Deltavph,
                long double Deltakap)
  {
  if( fabsl(DeltaPhi)  < par.threshold &&
      fabsl(Deltam)    < par.threshold &&
      fabsl(DeltaA)    < par.threshold &&
      fabsl(Deltathet) < par.threshold &&
      fabsl(Deltavph)  < par.threshold &&
      fabsl(Deltakap)  < par.threshold )
    return 1;
  else
    return 0;
  }

/*==========================================================================*/

void allocateInt()
  {
  Int.r     = (long double*) malloc((size_t) par.nint * sizeof(long double));
  Int.Phi   = (long double*) malloc((size_t) par.nint * sizeof(long double));
  Int.X     = (long double*) malloc((size_t) par.nint * sizeof(long double));
  Int.A     = (long double*) malloc((size_t) par.nint * sizeof(long double));
  Int.thet  = (long double*) malloc((size_t) par.nint * sizeof(long double));
  Int.vph   = (long double*) malloc((size_t) par.nint * sizeof(long double));
  Int.kap   = (long double*) malloc((size_t) par.nint * sizeof(long double));

  Int.rJ    = (long double*) malloc((size_t) par.nint * sizeof(long double));
  Int.m     = (long double*) malloc((size_t) par.nint * sizeof(long double));
  Int.mJ    = (long double*) malloc((size_t) par.nint * sizeof(long double));
  Int.C     = (long double*) malloc((size_t) par.nint * sizeof(long double));
  Int.CJ    = (long double*) malloc((size_t) par.nint * sizeof(long double));
  Int.NQ    = (long double*) malloc((size_t) par.nint * sizeof(long double));
  }

/*==========================================================================*/

void allocateExt()
  {
  ext.y     = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.Phi   = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.m     = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.zeta  = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.Pi    = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.rho   = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.Lam   = (long double*) malloc((size_t) par.next * sizeof(long double));

  ext.r     = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.rJ    = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.mJ    = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.X     = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.A     = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.thet  = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.vph   = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.kap   = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.C     = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.CJ    = (long double*) malloc((size_t) par.next * sizeof(long double));
  ext.NQ    = (long double*) malloc((size_t) par.next * sizeof(long double));
  }
/*==========================================================================*/

void registerPotential()
  {
  if(strcmp(par.potential, "series") == 0)
    {
    V    = (long double (*)(long double))(V_series);
    Vhat = (long double (*)(long double))(Vhat_series);
    Vp2  = (long double (*)(long double))(Vp2_series);
    }
  else if(strcmp(par.potential, "solitonic") == 0)
    {
    V    = (long double (*)(long double))(V_solitonic);
    Vhat = (long double (*)(long double))(Vhat_solitonic);
    Vp2  = (long double (*)(long double))(Vp2_solitonic);
    }
  else
    { printf("Unknown potential:   %s\n\n", par.potential); exit(0); }
  }

/*==========================================================================*/

long double F_DE(long double vph)
  {
  return exp(-2.0L * par.alp0 * vph - par.beta0 * vph*vph);
  }

/*==========================================================================*/

long double FpoF_DE(long double vph)
  {
  return -2.0L * (par.alp0 + par.beta0 * vph);
  }

/*==========================================================================*/

long double V_series(long double A)
  {
  // Potential function for the bosonic scalar
  return A*A * (1.0L + par.lambda4 * A*A + par.lambda6 * A*A*A*A
                + par.lambda8 * A*A*A*A*A*A);
  }

/*==========================================================================*/

long double Vhat_series(long double A)
  {
  // Potential function divided by A^2
  return (1.0L + par.lambda4 * A*A + par.lambda6 * A*A*A*A
          + par.lambda8 * A*A*A*A*A*A);
  }

/*==========================================================================*/

long double Vp2_series(long double A)
  {
  // Potential derivative dV / d(A^2)
  return 1.0L + 2.0L*par.lambda4 * A*A + 3.0L*par.lambda6 * A*A*A*A
         + 4.0L*par.lambda8 * A*A*A*A*A*A;
  }

/*==========================================================================*/

long double V_solitonic(long double A)
  {
  // Solitonic potential function for the bosonic scalar
  return A*A * (1.0L - 2.0L * A*A / (par.sigma0*par.sigma0))
             * (1.0L - 2.0L * A*A / (par.sigma0*par.sigma0));
  }

/*==========================================================================*/

long double Vhat_solitonic(long double A)
  {
  // Solitonic potential function divided by A^2
  return (1.0L - 2.0L * A*A / (par.sigma0*par.sigma0))
       * (1.0L - 2.0L * A*A / (par.sigma0*par.sigma0));
  }

/*==========================================================================*/

long double Vp2_solitonic(long double A)
  {
  // Solitonic potential derivative dV/d(A^2)
  return (1.0L - 2.0L * A*A / (par.sigma0*par.sigma0))
       * (1.0L - 6.0L * A*A / (par.sigma0*par.sigma0));
  }

/*==========================================================================*/

long double W(long double vph)
  {
  // Potential function for the gravitational scalar.
  // Note the difference in convention compared to V:
  // Here we have a factor 1/2.
  return 0.5L * par.massvph*par.massvph * vph*vph;
  }

/*==========================================================================*/

long double What(long double vph)
  {
  // Potential function for the gravitational scalar divided by vph^2
  return 0.5L * par.massvph*par.massvph;
  }

/*==========================================================================*/

long double Wp2(long double vph)
  {
  // Derivative of the potential function for the gravitational scalar.
  // Note, this is d/d(varphi), not d/d(varphi^2).
  return 0.5L * par.massvph*par.massvph;
  }

/*==========================================================================*/

void readPars(char* ifil)
  {
  FILE* ifp;
  char  line[LEN];
  int   i, n;


  // First set all parameters to default values
  par.ymax              = 1.0L;
  par.ymatch            = 0.0L;
  par.next              = 401;
  par.nint              = 401;
  par.Phi0              = -0.3785086583297151L; // Kaup limit
  par.mass0             = 0.6329974872153032L;  // Kaup limit
  par.zeta1             = 0.45064138L;          // Kaup limit
  par.rho1              = 0.0L;
  par.vph0              = 0.0L;
  par.A0start           = 0.07650501672240803L; // Kaup limit
  par.A0end             = 0.07650501672240803L; // Kaup limit
  par.dA0               = 0.1L;
  par.dA0floor          = 0.0001L;
  par.omega0            = 0.8525193102859463L;  // Kaup limit
  par.extrapolate       = 1;
  par.dataio            = 0;
  par.dfactor           = 10.0L;

  strcpy(par.potential, "series");
  par.lambda4   = 0.0L;
  par.lambda6   = 0.0L;
  par.lambda8   = 0.0L;
  par.sigma0    = 1.0L;
  par.massvph   = 0.5L;
  par.alp0      = 0.0L;
  par.beta0     = 0.0L;

  par.threshold = 1.0e-10L;
  par.EPS       = 1.0e-06L;
  par.accP      = 1.0L;
  par.accm      = 1.0L;
  par.accA      = 1.0L;
  par.accz      = 1.0L;
  par.acco      = 1.0L;
  par.accv      = 1.0L;
  par.accr      = 1.0L;
  par.shoot     = 1;
  par.verbose   = 0;


  // Read parameters from file and overwrite default values
  ifp = fopen(ifil, "r");
  if(ifp == NULL) { printf("Cannot open %s in readPars\n\n", ifil); exit(0); }

  while(mygetline(line, ifp) != 0)
    {
    if(line[0] == '#')
      ;         // Comment line
    else
      {
      n = strlen(line);
      // Replace '=' with white space for syntax flexibility
      for(i = 0; i < n; i++)
        if(line[i] == '=') line[i] = ' ';

      // Analyze parameter values
      if(strstr(line, "ymax ") != NULL)
        sscanf(line, "ymax %Lf", &(par.ymax));
      else if(strstr(line, "ymatch ") != NULL)
        sscanf(line, "ymatch %Lf", &(par.ymatch));
      else if(strstr(line, "next ") != NULL)
        sscanf(line, "next %d", &(par.next));
      else if(strstr(line, "nint ") != NULL)
        sscanf(line, "nint %d", &(par.nint));
      else if(strstr(line, "Phi0 ") != NULL)
        sscanf(line, "Phi0 %Lf", &(par.Phi0));
      else if(strstr(line, "mass0 ") != NULL)
        sscanf(line, "mass0 %Lf", &(par.mass0));
      else if(strstr(line, "zeta1 ") != NULL)
        sscanf(line, "zeta1 %Lf", &(par.zeta1));
      else if(strstr(line, "rho1 ") != NULL)
        sscanf(line, "rho1 %Lf", &(par.rho1));
      else if(strstr(line, "vph0 ") != NULL)
        sscanf(line, "vph0 %Lf", &(par.vph0));
      else if(strstr(line, "omega0 ") != NULL)
        sscanf(line, "omega0 %Lf", &(par.omega0));
      else if(strstr(line, "A0start ") != NULL)
        sscanf(line, "A0start %Lf", &(par.A0start));
      else if(strstr(line, "A0end ") != NULL)
        sscanf(line, "A0end %Lf", &(par.A0end));
      else if(strstr(line, "dA0 ") != NULL)
        sscanf(line, "dA0 %Lf", &(par.dA0));
      else if(strstr(line, "dA0floor ") != NULL)
        sscanf(line, "dA0floor %Lf", &(par.dA0floor));
      else if(strstr(line, "extrapolate ") != NULL)
        sscanf(line, "extrapolate %d", &(par.extrapolate));
      else if(strstr(line, "dataio ") != NULL)
        sscanf(line, "dataio %d", &(par.dataio));
      else if(strstr(line, "dfactor ") != NULL)
        sscanf(line, "dfactor %Lf", &(par.dfactor));
      else if(strstr(line, "potential ") != NULL)
        sscanf(line, "potential %s", par.potential);
      else if(strstr(line, "lambda4 ") != NULL)
        sscanf(line, "lambda4 %Lf", &(par.lambda4));
      else if(strstr(line, "lambda6 ") != NULL)
        sscanf(line, "lambda6 %Lf", &(par.lambda6));
      else if(strstr(line, "lambda8 ") != NULL)
        sscanf(line, "lambda8 %Lf", &(par.lambda8));
      else if(strstr(line, "sigma0 ") != NULL)
        sscanf(line, "sigma0 %Lf", &(par.sigma0));
      else if(strstr(line, "massvph ") != NULL)
        sscanf(line, "massvph %Lf", &(par.massvph));
      else if(strstr(line, "alp0 ") != NULL)
        sscanf(line, "alp0 %Lf", &(par.alp0));
      else if(strstr(line, "beta0 ") != NULL)
        sscanf(line, "beta0 %Lf", &(par.beta0));
      else if(strstr(line, "threshold ") != NULL)
        sscanf(line, "threshold %Lf", &(par.threshold));
      else if(strstr(line, "EPS ") != NULL)
        sscanf(line, "EPS %Lf", &(par.EPS));
      else if(strstr(line, "accP ") != NULL)
        sscanf(line, "accP %Lf", &(par.accP));
      else if(strstr(line, "accm ") != NULL)
        sscanf(line, "accm %Lf", &(par.accm));
      else if(strstr(line, "accA ") != NULL)
        sscanf(line, "accA %Lf", &(par.accA));
      else if(strstr(line, "accz ") != NULL)
        sscanf(line, "accz %Lf", &(par.accz));
      else if(strstr(line, "acco ") != NULL)
        sscanf(line, "acco %Lf", &(par.acco));
      else if(strstr(line, "accv ") != NULL)
        sscanf(line, "accv %Lf", &(par.accv));
      else if(strstr(line, "accr ") != NULL)
        sscanf(line, "accr %Lf", &(par.accr));
      else if(strstr(line, "shoot ") != NULL)
        sscanf(line, "shoot %d", &(par.shoot));
      else if(strstr(line, "verbose ") != NULL)
        sscanf(line, "verbose %d", &(par.verbose));
      }
    }

  if( fabsl(par.massvph) < 1.0e-15L )
    {
    printf("par.massvph = %Lg is converted to 0 (massless)!!!\n", par.massvph);
    par.massvph = 0.0L;
    par.massless = 1;
    }
  else
    par.massless = 0;
  }

/*==========================================================================*/

void printPars()
  {
  printf("\nPARAMETERS\n");
  printf("=======================================\n");
  printf("shoot         = %d\n", par.shoot);
  printf("A0start       = %22.16Lg\n", par.A0start);
  printf("A0end         = %22.16Lg\n", par.A0end);
  printf("dA0           = %22.16Lg\n", par.dA0);
  printf("dA0floor      = %22.16Lg\n", par.dA0floor);
  printf("Phi0          = %22.16Lg\n", par.Phi0);
  printf("mass0         = %22.16Lg\n", par.mass0);
  printf("vph0          = %22.16Lg\n", par.vph0);
  printf("zeta1         = %22.16Lg\n", par.zeta1);
  printf("rho1          = %22.16Lg\n", par.rho1);
  printf("omega0        = %22.16Lg\n", par.omega0);
  printf("\n");
  printf("ymax          = %22.16Lg\n", par.ymax);
  printf("ymatch        = %22.16Lg\n", par.ymatch);
  printf("next          = %d\n", par.next);
  printf("nint          = %d\n", par.nint);
  printf("extrapolate   = %d\n", par.extrapolate);
  printf("dataio        = %d\n", par.dataio);
  printf("threshold     = %22.16Lg\n", par.threshold);
  printf("EPS           = %22.16Lg\n", par.EPS);
  printf("dfactor       = %22.16Lg\n", par.dfactor);
  printf("\n");
  printf("potential     = %s\n", par.potential);
  if(strcmp(par.potential, "series") == 0)
    {
    printf("lambda4       = %22.16Lg\n", par.lambda4);
    printf("lambda6       = %22.16Lg\n", par.lambda6);
    printf("lambda8       = %22.16Lg\n", par.lambda8);
    }
  else if(strcmp(par.potential, "solitonic") == 0)
    printf("sigma0        = %22.16Lg\n", par.sigma0);
  printf("massvph       = %22.16Lg\n", par.massvph);
  printf("alp0          = %22.16Lg\n", par.alp0);
  printf("beta0         = %22.16Lg\n", par.beta0);
  printf("=======================================\n");
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

long double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a long double matrix with subscript range m[nrl..nrh][ncl..nch] */
  {
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  long double **m;

  /* allocate pointers to rows */
  m = (long double**) malloc((size_t) ((nrow+NR_END)*sizeof(long double*)));
  if(!m) nrerror("allocation failure 1 in dmatrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (long double*) malloc((size_t)
                                 ((nrow*ncol+NR_END)*sizeof(long double)));
  if (!m[nrl]) nrerror("allocation failure 2 in dmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i = nrl+1; i <= nrh; i++) m[i] = m[i-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
  }

/*==========================================================================*/

void free_dmatrix(long double** m, long nrl, long nrh, long ncl, long nch)
  {
  /* free a long double matrix allocated by dmatrix() */
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
  }

/*==========================================================================*/

void gaussj(long double **a, int n, long double **b, int m)
{
        int *indxc,*indxr,*ipiv;
        int i,icol,irow,j,k,l,ll;
        long double big,dum,pivinv,temp;

        indxc = (int*) ivector(1,n);
        indxr = (int*) ivector(1,n);
        ipiv  = (int*) ivector(1,n);
        for (j=1;j<=n;j++) ipiv[j]=0;
        for (i=1;i<=n;i++) {
                big=0.0L;
                for (j=1;j<=n;j++)
                        if (ipiv[j] != 1)
                                for (k=1;k<=n;k++) {
                                        if (ipiv[k] == 0) {
                                                if (fabsl(a[j][k]) >= big) {
                                                        big=fabsl(a[j][k]);
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
                if (a[icol][icol] == 0.0L) nrerror("gaussj: Singular Matrix-2");
                pivinv=1.0L/a[icol][icol];
                a[icol][icol]=1.0L;
                for (l=1;l<=n;l++) a[icol][l] *= pivinv;
                for (l=1;l<=m;l++) b[icol][l] *= pivinv;
                for (ll=1;ll<=n;ll++)
                        if (ll != icol) {
                                dum=a[ll][icol];
                                a[ll][icol]=0.0L;
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

void out1D(long double* x, long double* y, int n1, int n2, const char* ofil)
  {
  int   i;
  FILE* ofp;


  ofp = fopen(ofil, "w");
  if(ofp == NULL) { printf("Cannot open %s in out1D\n\n", ofil); exit(0); }

  for(i = n1; i <= n2; i++)
    fprintf(ofp, "%22.16Lg   %22.16Lg\n", x[i], y[i]);

  fclose(ofp);
  }

/*==========================================================================*/

int nancheck(long double y1, long double y2, long double y3, long double y4,
             long double y5, long double y6)
  {
  if( (y1!=y1) || (y2!=y2) || (y3!=y3) || (y4!=y4) || (y5!=y5) || (y6!=y6))
    return 1;
  else
    return 0;
  }

/*==========================================================================*/

void writeProfiles()
  {
  int nint, next;


  nint = par.nint;
  next = par.next;

  out1D(ext.y, ext.Phi,  0, next-1, "yPhi.dat");
  out1D(ext.y, ext.m,    0, next-1, "ym.dat");
  out1D(ext.y, ext.zeta, 0, next-1, "yzeta.dat");
  out1D(ext.y, ext.Pi,   0, next-1, "yPi.dat");
  out1D(ext.y, ext.rho,  0, next-1, "yrho.dat");
  out1D(ext.y, ext.Lam,  0, next-1, "yLam.dat");
  out1D(ext.y, ext.X,    0, next-1, "yX.dat");
  out1D(ext.y, ext.A,    0, next-1, "yA.dat");
  out1D(ext.y, ext.thet, 0, next-1, "ythet.dat");
  out1D(ext.y, ext.vph,  0, next-1, "yvph.dat");
  out1D(ext.y, ext.kap,  0, next-1, "ykap.dat");

  out1D(Int.r, Int.Phi,  0, nint-1, "rPhi.dat");
  out1D(Int.r, Int.X,    0, nint-1, "rX.dat");
  out1D(Int.r, Int.A,    0, nint-1, "rA.dat");
  out1D(Int.r, Int.thet, 0, nint-1, "rthet.dat");
  out1D(Int.r, Int.vph,  0, nint-1, "rvph.dat");
  out1D(Int.r, Int.kap,  0, nint-1, "rkap.dat");
  out1D(Int.r, Int.m,    0, nint-1, "rm.dat");

  allout1D(Int.r, ext.r, Int.Phi,  ext.Phi,  0, nint-1, 0, next-1, "Phi.dat");
  allout1D(Int.r, ext.r, Int.m,    ext.m,    0, nint-1, 0, next-1, "m.dat");
  allout1D(Int.r, ext.r, Int.X,    ext.X,    0, nint-1, 0, next-1, "X.dat");
  allout1D(Int.r, ext.r, Int.A,    ext.A,    0, nint-1, 0, next-1, "A.dat");
  allout1D(Int.r, ext.r, Int.thet, ext.thet, 0, nint-1, 0, next-1, "thet.dat");
  allout1D(Int.r, ext.r, Int.vph,  ext.vph,  0, nint-1, 0, next-1, "vph.dat");
  allout1D(Int.r, ext.r, Int.kap,  ext.kap,  0, nint-1, 0, next-1, "kap.dat");
  allout1D(Int.r, ext.r, Int.rJ,   ext.rJ,   0, nint-1, 0, next-1, "rJ.dat");
  allout1D(Int.r, ext.r, Int.mJ,   ext.mJ,   0, nint-1, 0, next-1, "mJ.dat");
  allout1D(Int.r, ext.r, Int.NQ,   ext.NQ,   0, nint-1, 0, next-1, "NQ.dat");
  }

/*==========================================================================*/

void allout1D(long double* xi, long double* xe, long double* yi,
              long double* ye, int nimin, int nimax, int nemin, int nemax,
              const char* ofil)
  {
  int   i;
  FILE* ofp;


  ofp = fopen(ofil, "w");
  if(ofp == NULL) { printf("Cannot open %s in out1D\n\n", ofil); exit(0); }

  fprintf(ofp, "# %s\n", ofil);

  for(i = nimin; i <= nimax; i++)
    fprintf(ofp, "%22.16Lg   %22.16Lg\n", xi[i], yi[i]);

  // Note that the exterior is written in reverse order. Also, the first
  // point on the y grid is the same as the last point on the inner
  // grid, so we don't need to write it again and therefore start with
  // nemax-1.
  for(i = nemax-1; i >= nemin; i--)
    fprintf(ofp, "%22.16Lg   %22.16Lg\n", xe[i], ye[i]);

  fclose(ofp);
  }

/*==========================================================================*/
