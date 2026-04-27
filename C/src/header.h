/*
File: header.h

This file contains libraries, function declarations, and other necessary definitions.
*/

// Basic libraries
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <unistd.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <errno.h>

#ifdef _WIN32
    #include <direct.h>
    #define MAKE_DIR(path) _mkdir(path)
#else
    #include <sys/stat.h>
    #include <sys/types.h>
    #define MAKE_DIR(path) mkdir((path), 0755)
#endif

// Constants
#define Pi (3.1415926535897932384626433832795028)
#define EulerGamma (0.5772156649015328606065121)

// List of options for the type of motion
enum{
    Motion_Circular,
    Motion_EquatorialEccentric,
    Motion_Generic,
    Motion_NOPT
};
static const char* const motion_opt[] = {"circular", "equatorial-eccentric", "generic", "undefined"};

// List of options for the integration scheme
enum{
    ODESolver_RK4,
    ODESolver_RKGL6,
    ODESolver_NOPT
};
static const char* const solver_opt[] = {"RK4", "RK-GL6", "undefined"};

// List of options for the coordinates
enum{
    Coords_Standard,
    Coords_Canonical,
    Coords_NOPT
};
static const char* const coords_opt[] = {"standard", "canonical", "undefined"};

// List of options for the time step
enum{
    Step_Standard,
    Step_Transformed,
    Step_NOPT
};
static const char* const step_opt[] = {"standard", "transformed", "undefined"};

// List of options for the potentials  
enum{
    Pots_NonResummed,
    Pots_Resummed,
    Pots_NOPT
};
static const char* const pots_opt[] = {"non-resummed", "resummed", "undefined"};

// List of options for the Poincare section surface
enum{
    PoincareSurface_Z,
    PoincareSurface_X,
    PoincareSurface_NOPT
};
static const char* const poincare_surface_opt[] = {"z", "x", "undefined"};

// List of options for the Poincare section crossing direction
enum{
    PoincareDirection_Positive,
    PoincareDirection_Negative,
    PoincareDirection_Both,
    PoincareDirection_NOPT
};
static const char* const poincare_direction_opt[] = {"positive", "negative", "both", "undefined"};

// Data type for the parameters
typedef struct Parameters {

    // Binary parameters
    double nu; // symmetric mass ratio: m1*m2 / (m1 + m2)^2
    double q; // mass ratio m1/m2 (with m1 > m2, so q >= 1)
    double X1; // ratio m1/M
    double X2; // ratio m2/M
    double chi1x0, chi1y0, chi1z0; // spin components, with chi1 = S1 / (m1^2)
    double chi2x0, chi2y0, chi2z0; // spin components, with chi1 = S2 / (m2^2)

    // Orbital parameters
    double x0, y0, z0;

    // Type of motion
    int motion;

    // Constants of motion
    double E0, pphi0;
    
    // Options for functions within the Hamiltonian
    int pots;

    // Type of coordinates (standard with non-canonical spins, or canonical)
    int coords;

    // ODE solver settings
    int solver; // Choose the type of solver
    int step; // Standard or rescaled
    double dt;
    double tmax_traj;
    double tmax_poincare;
    double max_iter_RKGL6;
    double tol_RKGL6;
    // double reltol;
    // double abstol;

    // Poincare section settings
    int poincare_on;
    int poincare_surface;
    int poincare_direction;
    double poincare_value;

    // Poincare scan settings
    int scan_on;
    double scan_rmin;
    double scan_rmax;
    int scan_nr;

    // Rotation-number post-processing settings
    int rotation_min_crossings_for_center;
    int rotation_n_center_orbits;

} Parameters;

extern Parameters *pars; 
// *pars variable is of type Parameters; defined in pars.c

// Functions in pars.c
void AllocateParameters(Parameters **pars);
void SetDefaults(Parameters *pars);
void ParseCommandLine(int argc, char *argv[], Parameters *pars);
void ReadParfile(char *filename, Parameters *pars);
void AssignValues(Parameters *pars, char *key, char *value);
void SetParameters(void); // no arguments
void FreeParameters(Parameters *pars);
void WriteMetadataFile(Parameters *pars, char *filepath);

// Functions in dyn.c
void Hamiltonian(double r[], double p[], double nu, double chi1[], double chi2[], double *Heff, double *H, double *dHeff, double *dH);
void get_rhs(double t, double w[], double dw[]);
void get_rhs_canonical(double t, double W[], double dW[]);
void get_rhs_transformed(double s, double w[], double dw[]);
void get_rhs_canonical_transformed(double s, double W[], double dW[]);
double CarterLikeConstant(double r[], double p[], double chi1[], double chi2[]);
double PhotonPotentialCondition(double modr);
void AdiabaticLightRing(double *rLR);

// Functions in metric.c
extern void (*Potentials)(double r[], double p[], double nu, double chi1[], double chi2[], double *A, double *Bp, double *Bnp, double *Benp, double *dA, double *dBp, double *dBnp, double *dBenp);
void NonResummedPotentials(double r[], double p[], double nu, double chi1[], double chi2[], double *A, double *Bp, double *Bnp, double *Benp, double *dA, double *dBp, double *dBnp, double *dBenp);
void ResummedPotentials(double r[], double p[], double nu, double chi1[], double chi2[], double *A, double *Bp, double *Bnp, double *Benp, double *dA, double *dBp, double *dBnp, double *dBenp);
void NLO_SpinSpinContributions(double r[], double nu, double chi1[], double chi2[], double *AchiQ, double *AnchiQ, double *BchiQ, double *BnchiQ, double *dAchiQ, double *dAnchiQ, double *dBchiQ, double *dBnchiQ);
void get_Aeq(double r[], double nu, double chi1[], double chi2[], double *Aeq, double *dAeq); 
void get_AeqB(double r[], double nu, double chi1[], double chi2[], double *AeqB, double *dAeqB);
void get_Aorb(double modr, double nu, double *Aorb, double *dAorb); 
void get_Dorb(double modr, double nu, double *Dorb, double *dDorb); 
void get_Q4(double r[], double p[], double nu, double chi1[], double chi2[], double *Q4, double *dQ4); 
void get_GSGSstar(double r[], double p[], double nu, double chi1[], double chi2[], double *GS, double *GSs, double *dGS, double *dGSs); 

// Functions in ics.c 
void InitialGuess_KerrAngularMomentum(double *l);
double get_dHeffdx(double pphi);
void CircularICs(double *py0);
void KerrEquatorialTurningPoints(double *tp);
double EnergyCondition(double x);
double get_dHeffdx_Equatorial(double x); 
void EquatorialICs(double *x0);
double EnergyCondition_Generic(double pz);
double get_dHeffdpz(double pz);
void GenericICs(double *p0);
int TryGenericICsAtRadius(double r0, double *p0);

// Functions in utils.c
double dot(double a[], double b[], int n);
void cross(double a[], double b[], double *c);
double get_mod(double x[]);
double get_nu(double q);
void get_n(double r[], double *n, double *dn);
void get_pr(double r[], double p[], double *pr, double *dpr);
void get_X1X2(double nu, double *X1, double *X2);
void get_rc(double r[], double nu, double chi1[], double chi2[], double *rc, double *drc);
void get_Delta(double r[], double nu, double chi1[], double chi2[], double *Delta, double *dDelta);
void get_a0(double nu, double chi1[], double chi2[], double *a0);
void get_SSstar(double nu, double chi1[], double chi2[], double *S, double *Ss);
void get_R4(double r[], double nu, double chi1[], double chi2[], double *R4, double *dR4);
void get_l(double r[], double p[], double *l, double *dl);
void get_lSSs(double r[], double p[], double nu, double chi1[], double chi2[], double *lS, double *dlS, double *lSs, double *dlSs);
void trimString(char *str);
void swap(double *xp, double *yp);
void bubbleSort(double arr[], int n);
void CartesianToSpherical(double r[], double p[], double *Q, double *P);

// Functions in solvers.c
extern void (*ODESolver)(double t, double *Y, double h, int N_eq, void (*RHS)(double, double*, double*));
void RK_Fourth(double t, double *Y, double h, int N_eq, void (*RHS)(double, double*, double*));
void RK_GaussLegendre6(double t, double *Y, double h, int N_eq, void (*RHS)(double, double*, double*));
void Secant(double(*fun)(double), double x0, double tol, double *zero);
void Secant_v2(double(*fun)(double), double a, double b, double tol, double *zero);
double Newton(double(*fun)(double), double(*fun_der)(double), double a, double b, double tol);
void QuarticPolyRoots(double coefs[], double roots[]);

// Functions in run.c
void BuildOutputFolderName(char *folder, size_t folder_size);
int RunScanMode(const char *folder);
int RunSingleMode(const char *folder);

// Functions in rotation.c
int PostProcessRotationScan(const char *folder, const char *poincare_scan_path);