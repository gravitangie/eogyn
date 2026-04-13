/*
File: ics.c

This file contains functions for evaluating the initial conditions.
*/

#include "header.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                                     Circular Initial Conditions                                    //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

void InitialGuess_KerrAngularMomentum(double *l)
{
    // All input parameters
    double nu = pars->nu;
    double x0 = pars->x0;
    double y0 = pars->y0;
    double z0 = pars->z0;

    double r0[3];
    r0[0] = x0;
    r0[1] = y0;
    r0[2] = z0;

    double chi1[3], chi2[3];
    chi1[0] = pars->chi1x0;
    chi1[1] = pars->chi1y0;
    chi1[2] = pars->chi1z0;
    chi2[0] = pars->chi2x0;
    chi2[1] = pars->chi2y0;
    chi2[2] = pars->chi2z0;

    double a0[3];
        get_a0(nu, chi1, chi2, a0);

    // |a0| as an effective Kerr parameter
    double a = sqrt(dot(a0, a0, 3));

    // inverse radius
    double u     = 1./get_mod(r0);
    double u2    = u*u;
    double u3by2 = sqrt(u2*u);

    // Kerr angular momentum
    double num = 1. - 2.*a*u3by2 + a*a*u2;
    double den = sqrt(u*(1. - 3.*u + 2.*a*u3by2));
    *l = num/den;
}

// Provides the derivative of Heff to check for circular motion - used in CirculaICs
double get_dHeffdx(double pphi)
{
    double nu = pars->nu;

    double r0[3];
    r0[0] = pars->x0;
    r0[1] = pars->y0;
    r0[2] = pars->z0;

    // FIXME: if (y0 != 0. || z0 != 0.) add error message

    double chi1[3], chi2[3];
    chi1[0] = pars->chi1x0;
    chi1[1] = pars->chi1y0;
    chi1[2] = pars->chi1z0;
    chi2[0] = pars->chi2x0;
    chi2[1] = pars->chi2y0;
    chi2[2] = pars->chi2z0;

    double p[3]; 
    p[0] = 0.;
    p[1] = pphi/r0[0];
    p[2] = 0.;
    
    double Heff, H, dHeff[12];
    Hamiltonian(r0, p, nu, chi1, chi2, &Heff, &H, dHeff, NULL);

    return dHeff[0] - (pphi/(r0[0]*r0[0]))*dHeff[4];

}

void CircularICs(double *pphi0)
{
    double l;
    InitialGuess_KerrAngularMomentum(&l);
    printf("Initial guess for Secant method: %.16f\n", l);

    double zero;
    Secant(get_dHeffdx, l, 1e-10, &zero);
    printf("Secant method found pphi0: %.16f\n", zero);

    *pphi0 = zero;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                                   Eccentric Initial Conditions                                     //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

// This function evaluates the Kerr apoapsis for equatorial motion:
// it exploits the input energy and angular momentum, together with the effective Kerr parameter, within the Kerr radial function;
// then it finds the roots of such function (which is a quartic polynomial) and takes the largest.
void KerrEquatorialTurningPoints(double *tp) {

    double E   = pars->E0;
    double E2  = E*E;
    double Lz  = pars->pphi0;
    double Lz2 = Lz*Lz;

    double chi1[3], chi2[3];
    chi1[0] = pars->chi1x0;
    chi1[1] = pars->chi1y0;
    chi1[2] = pars->chi1z0;
    chi2[0] = pars->chi2x0;
    chi2[1] = pars->chi2y0;
    chi2[2] = pars->chi2z0;

    double nu = pars->nu;
    
    double a0[3];
        get_a0(nu, chi1, chi2, a0);

    // |a0| as an effective Kerr parameter
    double a  = sqrt(dot(a0, a0, 3));
    double a2 = a*a;

    printf("Effective Kerr parameter: %.10f\n", a);

    // the radial function for equatorial motion (Q = 0) is a quartic polynomial:
    // a x^4 + b x^3 + c x^2 + d x + e 
    // with e = 0

    double coefs[5];
    coefs[0] = - 1. + E2;
    coefs[1] = 2.;
    coefs[2] = - a2 + a2*E2 - Lz2;
    coefs[3] = 2.*a2*E2 - 4.*a*E*Lz + 2.*Lz2;
    coefs[4] = 0.;

    double roots[4];
        QuarticPolyRoots(coefs, roots);

    // Sort the roots
    bubbleSort(roots, 4);

    printf("Roots of the Kerr radial function: %.10f, %.10f, %.10f, %.10f\n", roots[0], roots[1], roots[2], roots[3]);

    // The turning points for bound motion are the two largest roots (other ones are unstable)
    tp[0] = roots[2];
    tp[1] = roots[3];

}

// Get Heff(x) - E0
double EnergyCondition(double x)
{
    double nu = pars->nu;

    double r0[3];
    r0[0] = x;
    r0[1] = 0.; // y
    r0[2] = 0.; // z


    double pphi = pars->pphi0;

    double chi1[3], chi2[3];
    chi1[0] = pars->chi1x0;
    chi1[1] = pars->chi1y0;
    chi1[2] = pars->chi1z0;
    chi2[0] = pars->chi2x0;
    chi2[1] = pars->chi2y0;
    chi2[2] = pars->chi2z0;

    double p[3]; 
    p[0] = 0.;     // px
    p[1] = pphi/x; // py
    p[2] = 0.;     // pz
    
    double Heff, H;
    Hamiltonian(r0, p, nu, chi1, chi2, &Heff, &H, NULL, NULL);

    // printf("Heff = %.16f, E0 = %.16f, Heff - E0 = %.16f\n", Heff, pars->E0, Heff - (pars->E0));
    return Heff - (pars->E0);

}

// Get the derivative of Heff with respect to x
double get_dHeffdx_Equatorial(double x)
{

    double nu = pars->nu;

    double r0[3];
    r0[0] = x;
    r0[1] = 0.; // y
    r0[2] = 0.; // z

    double pphi = pars->pphi0;

    double chi1[3], chi2[3];
    chi1[0] = pars->chi1x0;
    chi1[1] = pars->chi1y0;
    chi1[2] = pars->chi1z0;
    chi2[0] = pars->chi2x0;
    chi2[1] = pars->chi2y0;
    chi2[2] = pars->chi2z0;

    double p[3]; 
    p[0] = 0.;     // px
    p[1] = pphi/x; // py
    p[2] = 0.;     // pz

    double Heff, H, dHeff[12];
    Hamiltonian(r0, p, nu, chi1, chi2, &Heff, &H, dHeff, NULL);

    // return nu*dH[0];    
    return dHeff[0] - (pphi/(x*x))*dHeff[4];

}

void EquatorialICs(double *x0)
{
    double Kerr_tp[2];
        KerrEquatorialTurningPoints(Kerr_tp);
    
    // printf("Initial guesses (Kerr periapsis and apoapsis): %.16f, %.16f\n", Kerr_tp[0], Kerr_tp[1]);

    // Check that both turning points also exist within EOB
    double zero1, zero2;
    // Periapsis
    zero1 = Newton(EnergyCondition, get_dHeffdx_Equatorial, Kerr_tp[0] - 1., Kerr_tp[0] + 1., 1e-10);
    // Apoapsis
    zero2 = Newton(EnergyCondition, get_dHeffdx_Equatorial, Kerr_tp[1] - 1., Kerr_tp[1] + 1., 1e-10);

    // printf("EOB periapsis and apoapsis: %.16f, %.16f\n", zero1, zero2);

    if (isnan(zero1) || isnan(zero2) || fabs(zero2 - zero1) < 1e-15) {
        perror("Error: Turning points for bound EOB motion not found! Choose different initial energy and angular momentum.\n");
        exit(EXIT_FAILURE);
    }

    printf("Starting motion at the EOB apoapsis: %.16f\n", zero2);

    // Start the motion at the apoapsis
    *x0 = zero2;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                                     Generic initial conditions                                     //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
static int BracketGenericPz(double pz_min, double pz_max, int nscan, double *a, double *b)
{
    const double tol = 1e-12;

    double pz_prev = pz_min;
    double f_prev  = EnergyCondition_Generic(pz_prev);

    if (fabs(f_prev) < tol) {
        *a = pz_prev;
        *b = pz_prev;
        return 1;
    }

    for (int i = 1; i <= nscan; i++) {
        double frac = (double)i / (double)nscan;
        double pz    = pz_min + frac * (pz_max - pz_min);
        double f     = EnergyCondition_Generic(pz);

        if (fabs(f) < tol) {
            *a = pz;
            *b = pz;
            return 1;
        }

        if (f_prev * f < 0.) {
            *a = pz_prev;
            *b = pz;
            return 1;
        }

        pz_prev = pz;
        f_prev  = f;
    }

    return 0;
}

void GenericICs(double *p0)
{
    // Get py with the initial angular momentum
    double py0 = (pars->pphi0) / (pars->x0);

    // Keep the same branch as in the scan: pz <= 0
    double a, b, pz0;
    int ok = BracketGenericPz(-0.5, 0.0, 200, &a, &b);

    if (!ok) {
        perror("Error: could not bracket a valid pz root for generic initial conditions.\n");
        exit(EXIT_FAILURE);
    }

    if (fabs(b - a) < 1e-14) {
        pz0 = a;
    } else {
        Secant_v2(EnergyCondition_Generic, a, b, 1e-10, &pz0);
    }

    double residual = fabs(EnergyCondition_Generic(pz0));
    if (residual > 1e-8) {
        perror("Error: generic initial conditions do not satisfy the energy constraint accurately enough.\n");
        exit(EXIT_FAILURE);
    }

    printf("Secant method found pz0: %.16f\n", pz0);

    p0[0] = 0.;
    p0[1] = py0;
    p0[2] = pz0;
}

int TryGenericICsAtRadius(double r0, double *p0)
{
    // Save current values
    double x0_old = pars->x0;
    double y0_old = pars->y0;
    double z0_old = pars->z0;

    // Scan is along p_r = 0 on the section z = 0, y = 0
    pars->x0 = r0;
    pars->y0 = 0.;
    pars->z0 = 0.;

    double py0 = pars->pphi0 / pars->x0;

    // We keep the same branch as in GenericICs: pz <= 0
    double a, b, pz0;
    int ok = BracketGenericPz(-0.5, 0.0, 200, &a, &b);

    if (!ok) {
        pars->x0 = x0_old;
        pars->y0 = y0_old;
        pars->z0 = z0_old;
        return 0;
    }

    if (fabs(b - a) < 1e-14) {
        pz0 = a;
    } else {
        Secant_v2(EnergyCondition_Generic, a, b, 1e-10, &pz0);
    }

    // Optional residual check
    double residual = fabs(EnergyCondition_Generic(pz0));
    if (residual > 1e-8) {
        pars->x0 = x0_old;
        pars->y0 = y0_old;
        pars->z0 = z0_old;
        return 0;
    }

    p0[0] = 0.;
    p0[1] = py0;
    p0[2] = pz0;

    // Keep scanned radius in pars for the orbit that will be run next
    return 1;
}

// Get Heff as a function of pz (for generic ICs)
double EnergyCondition_Generic(double pz)
{
    double nu = pars->nu;

    double r0[3];
    r0[0] = pars->x0;
    r0[1] = pars->y0;
    r0[2] = pars->z0;

    double chi1[3], chi2[3];
    chi1[0] = pars->chi1x0;
    chi1[1] = pars->chi1y0;
    chi1[2] = pars->chi1z0;
    chi2[0] = pars->chi2x0;
    chi2[1] = pars->chi2y0;
    chi2[2] = pars->chi2z0;

    double p[3]; 
    p[0] = 0.;
    p[1] = (pars->pphi0)/(pars->x0);
    p[2] = pz;

    double Heff, H;
    Hamiltonian(r0, p, nu, chi1, chi2, &Heff, &H, NULL, NULL);

    // return H*nu - (pars->E0);
    return Heff - (pars->E0);

}

// Get the derivative of the effective Hamiltonian with respect to pz (for generic ICs)
double get_dHeffdpz(double pz)
{

    double nu = pars->nu;

    double r0[3];
    r0[0] = pars->x0;
    r0[1] = pars->y0;
    r0[2] = pars->z0;

    double chi1[3], chi2[3];
    chi1[0] = pars->chi1x0;
    chi1[1] = pars->chi1y0;
    chi1[2] = pars->chi1z0;
    chi2[0] = pars->chi2x0;
    chi2[1] = pars->chi2y0;
    chi2[2] = pars->chi2z0;

    double p[3]; 
    p[0] = 0.;
    p[1] = (pars->pphi0)/(pars->x0);
    p[2] = pz;

    double Heff, H, dHeff[12];
    Hamiltonian(r0, p, nu, chi1, chi2, &Heff, &H, dHeff, NULL);

    // return nu*dH[5];    
    return dHeff[5]; 
}