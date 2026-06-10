/*
File: dyn.c

This file provides the implementation of the effective one-body (EOB) Hamiltonian, as described in the paper:
https://doi.org/10.1103/PhysRevD.92.124022.
It also includes the computation of the right-hand side of the Hamiltonian equations.
*/

#include "header.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                                        HAMILTONIAN FUNCTION                                        //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Computes the EOB Hamiltonian
// [in]  r              Radial vector
// [in]  p              Momentum vector
// [in]  nu             Symmetric mass ratio
// [in]  chi1           Spin vector of the first particle
// [in]  chi2           Spin vector of the second particle
// [out] Heff           Effective EOB Hamiltonian (devided by mu)
// [out] H              Real EOB Hamiltonian (devided by mu)
// [out] dHeff          Derivatives of Heff wrt all spacetime variables
// [out] dH             derivatives of H wrt all spacetime variables
void Hamiltonian(double r[], double p[], double nu, double chi1[], double chi2[], double *Heff, double *H, double *dHeff, double *dH)
{

    // ******************
    // * Initialization *
    // ******************

    // Mass fractions
    double X1, X2;
        get_X1X2(nu, &X1, &X2);
        double X[] = {X1, X2};

    // Unit radial vector and derivatives
    double n[3], dn[3][3];
        get_n(r, n, (double *)dn);

    // Centrifugal radius and derivatives
    double rc, drc[4][3];
        get_rc(r, nu, chi1, chi2, &rc, (double *)drc);

    // Rescaled orbital angular momentum
    double l[3], dl[6][3];
        get_l(r, p, l, (double *)dl);

    // Effective spin vector
    double a0[3];
        get_a0(nu, chi1, chi2, a0);

    // Radial momentum
    double pr, dpr[2][3];
        get_pr(r, p, &pr, (double *)dpr);

    // Symmetric spin combinations
    double S[3], Ss[3];
        get_SSstar(nu, chi1, chi2, S, Ss);

    // Dot product of the rescaled orbital angular momentum with the spin combinations
    double lS, lSs, dlS[4][3], dlSs[4][3];
        get_lSSs(r, p, nu, chi1, chi2, &lS, (double *)dlS, &lSs, (double *)dlSs);

    // Delta function and derivatives
    double Delta, dDelta[4][3];
        get_Delta(r, nu, chi1, chi2, &Delta, (double *)dDelta);

    // R4 function and derivatives
    double R4, dR4[4][3];
        get_R4(r, nu, chi1, chi2, &R4, (double *)dR4);

    // Potentials and derivatives
    double A, Bp, Bnp, Benp, dA[4][3], dBp[4][3], dBnp[4][3], dBenp[4][3];
        Potentials(r, p, nu, chi1, chi2, &A, &Bp, &Bnp, &Benp, (double *)dA, (double *)dBp, (double *)dBnp, (double *)dBenp);

    // Q4
    double Q4, dQ4[4][3];
        get_Q4(r, p, nu, chi1, chi2, &Q4, (double *)dQ4);

    // Gyro-gravitomagnetic functions + derivatives
    double GS, GSs, dGS[4][3], dGSs[4][3];
        get_GSGSstar(r, p, nu, chi1, chi2, &GS, &GSs, (double *)dGS, (double *)dGSs);

    // Auxiliary expressions
    double modr     = get_mod(r);     // Modulus of r
    double u        = 1. / modr;      // Invese of |r|
    double u2       = u * u;          // Inverse of |r|^2
    double na0      = dot(n, a0, 3);  // Dot product of n and a0

    double na02     = na0*na0;        // Squared dot product of n and a0
    double uc       = 1. / rc;        // Inverse of rc
    double uc2      = uc * uc;        // Inverse of rc*rc
    double p2       = dot(p, p, 3);   // Squared modulus of p
    double pr2      = pr * pr;        // Squared pr

    // Cross products
    double n_cross_p[3];
        cross(n, p, n_cross_p);   
    double p_cross_a0[3];
        cross(p, a0, p_cross_a0);
    double a0_cross_n[3];
        cross(a0, n, a0_cross_n);

    // ***************
    // * Hamiltonian *
    // ***************

    // Auxiliary expressions
    double a0_n_cross_p  = dot(a0, n_cross_p, 3);
    double a0_n_cross_p2 = a0_n_cross_p*a0_n_cross_p; 
    double BQ = 1. + Bp*p2 + Bnp*pr2 + Benp*a0_n_cross_p2 + Q4; // term with the B's and Q4

    // Orbital part of the effective Hamiltonian
    double Horb = sqrt( A * BQ );
    // Spin-orbit part of the effective Hamiltonian
    double Hso = lS * GS + lSs * GSs;
    // Effective Hamiltonian (This is actually \hat{Heff} = Heff / mu)
    double Heff_val = Horb + Hso;
    *Heff = Heff_val;
    // EOB Hamiltonian (This is actually \hat{HEOB} = HEOB / mu)
    double H_val = sqrt(1. + 2.*nu*(Heff_val - 1.))/nu;
    *H = H_val;

    // **********************************
    // * Derivatives of the Hamiltonian *
    // **********************************

    if (dHeff != NULL)
    {
        // General auxiliary expressions
        double over_Horb = 1. / (Horb); // Inverse of the orbital part of the Hamiltonian
        double da0_n_cross_p2_dr[3];
        double dBQ[4][3]; // Derivatives of BQ

        for (int i = 0; i < 3; i++)
        {
            // ***********************
            // * Spatial derivatives *
            // ***********************

            // Auxiliary expression 
            da0_n_cross_p2_dr[i] = 2. * a0_n_cross_p * dot(dn[i], p_cross_a0, 3); 

            // Derivative of BQ
            dBQ[0][i] = dBp[0][i]*p2 + dBnp[0][i]*pr2 + Bnp*2.*pr*dpr[0][i] + dBenp[0][i]*a0_n_cross_p2 + Benp*da0_n_cross_p2_dr[i] + dQ4[0][i];

            // Spatial derivative of the orbital part of the Hamiltonian
            double dHorbdr = 0.5*over_Horb * ( dA[0][i]*BQ + A*dBQ[0][i] );
            // Spatial derivative of the spin-orbit part of the Hamiltonian
            double dHsodr = dGS[0][i] * lS + GS * dlS[0][i] + dGSs[0][i] * lSs + GSs * dlSs[0][i];
            // Spatial derivative of the effective Hamiltonian
            dHeff[i] = dHorbdr + dHsodr;

            // ************************
            // * Momentum derivatives *
            // ************************

            // Derivative of BQ
            // Notes:
            // - B's have no dependence on p
            // - d(p_vec)/dpx = (1 0 0) = hat{x}, and: a0 dot (n cross hat{x}) = hat{x} dot (a0 cross n) = (a0 cross n)[i]
            dBQ[1][i] = Bp*2.*p[i] + Bnp*2.*pr*dpr[1][i] + Benp*2.*a0_n_cross_p*a0_cross_n[i] + dQ4[1][i];
            
            // Momentum derivative of the orbital part of the Hamiltonian
            double dHorbdp = 0.5*over_Horb * ( A*dBQ[1][i] ); // dA[1][i] = 0, A does not depend on p
            // Momentum derivative of the spin-orbit part of the Hamiltonian
            double dHsodp = dGS[1][i] * lS + GS * dlS[1][i] + dGSs[1][i] * lSs + GSs * dlSs[1][i];
            // Momentum derivative of the effective Hamiltonian
            dHeff[i + 3] = dHorbdp + dHsodp;

            // ********************
            // * Spin derivatives *
            // ********************

            for (int j = 0; j < 2; j++)
            {
                int k = j + 2; // j = 0 -> chi1, d...[2][i], j = 1 -> chi2, d...[3][i]

                // Derivative of BQ
                dBQ[k][i] = dBp[k][i]*p2 + dBnp[k][i]*pr2 + dBenp[k][i]*a0_n_cross_p2 + Benp*2.*a0_n_cross_p*X[j]*n_cross_p[i] + dQ4[k][i];
                
                // Spin derivative of the orbital part of the Hamiltonian
                double dHorbdchi = 0.5*over_Horb * ( dA[k][i]*BQ + A*dBQ[k][i] );
                // Spin derivative of the spin-orbit part of the Hamiltonian
                double dHsodchi = dGS[k][i] * lS + GS * dlS[k][i] + dGSs[k][i] * lSs + GSs * dlSs[k][i];
                // Spin derivative of the effective Hamiltonian
                dHeff[i + 6 + 3*j] = dHorbdchi + dHsodchi;

            }

        }
    }

    if (dH != NULL)
    {
        double over_E = 1./(nu*(H_val)); // 1/E where E is the real energy
        for (int i = 0; i < 12; i++) dH[i] = dHeff[i]*over_E;
    }

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                             RIGHT-HAND SIDE FUNCTION: Standard                                     //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Computes the RHS of the Hamiltonian equations
// [in]  t              time
// [in]  w              Master variable
// [out] dw             Derivatives of w wrt all spacetime variables
void get_rhs(double t, double w[], double dw[])
{
    double nu = pars->nu;

    // Unwrap
    double r[3], p[3], chi1[3], chi2[3];
    for (int i = 0; i < 3; i++) {
        r[i]     = w[i];
        p[i]     = w[i + 3];
        chi1[i]  = w[i + 6];
        chi2[i]  = w[i + 9];
    }

    // Derivatives of the Hamiltonian
    double Heff, H, dHeff[12], dH[12];
        Hamiltonian(r, p, nu, chi1, chi2, &Heff, &H, dHeff, dH);

    // Extract derivatives wrt to the spins chi1, chi2
    double dHdchi1[3], dHdchi2[3];
    for (int i = 0; i < 3; i++)
    {
        dHdchi1[i] = dH[6 + i];
        dHdchi2[i] = dH[9 + i];
    }

    // Cross products
    double dHdchi1xchi1[3], dHdchi2xchi2[3];
        cross(dHdchi1, chi1, dHdchi1xchi1);
        cross(dHdchi2, chi2, dHdchi2xchi2);

    // Spin EOM: dchi_i/dt = (1/m_i^2)(dH_phys/dchi_i) x chi_i. The code uses Hhat = H_phys/mu
    // and dimensionless time t = T/M, so the prefactor on (dHhat/dchi_i) is mu/m_i^2 = nu/X_i^2
    // (M=1) -- NOT 1/X_i^2 (see the spin-EOM discussion in the paper).
    double X1, X2;
        get_X1X2(nu, &X1, &X2);
        double iX1sq = nu/(X1*X1); // prefactor mu/m_i^2 = nu/X_i^2 (M=1); see Sec. on spin EOMs
        double iX2sq = nu/(X2*X2);

    // *******
    // * rhs *
    // *******

    for (int i = 0; i < 3; i++)
    {
        dw[i]     = dH[i + 3];                 // x, y, z
        dw[i + 3] = - dH[i];                   // px, py, pz
        dw[i + 6] = iX1sq*dHdchi1xchi1[i];     // chi1x, chi1y, chi1z
        dw[i + 9] = iX2sq*dHdchi2xchi2[i];     // chi2x, chi2y, chi2z
    }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                        RIGHT-HAND SIDE FUNCTION: Canonical spins                                   //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Computes the RHS of the Hamiltonian equations
// [in]  t              time
// [in]  W              Master variable
// [out] dW             Derivatives of w wrt all spacetime variables
void get_rhs_canonical(double t, double W[], double dW[])
{
    double nu = pars->nu;

    // Unwrap
    double r[3], p[3], alpha1, xi1, alpha2, xi2;
    r[0]   = W[0];
    r[1]   = W[1];
    r[2]   = W[2];
    p[0]   = W[3];
    p[1]   = W[4];
    p[2]   = W[5];
    alpha1 = W[6];
    xi1    = W[7];
    alpha2 = W[8];
    xi2    = W[9];
    
    // FIXME: this could be changed and updated during the evolution by defining pars->modchi and updating it in the dynamics.
    // However even with the non-canonical evolution, the error on the conservation of the spin magnitude is of the order of 10^-14.
    double modchi1 = sqrt((pars->chi1x0)*(pars->chi1x0) + (pars->chi1y0)*(pars->chi1y0) + (pars->chi1z0)*(pars->chi1z0));
    double modchi2 = sqrt((pars->chi2x0)*(pars->chi2x0) + (pars->chi2y0)*(pars->chi2y0) + (pars->chi2z0)*(pars->chi2z0));

    // Reconstruct the spins in the (rotated) chart frame, then rotate to lab
    // for the Hamiltonian call. (alpha_i, xi_i) live in the chart frame.
    double chi1_chart[3], chi2_chart[3], chi1[3], chi2[3];
    double sqrt1mxi12 = sqrt(1. - xi1*xi1);
    double sqrt1mxi22 = sqrt(1. - xi2*xi2);

    chi1_chart[0] = modchi1 * sqrt1mxi12 * cos(alpha1);
    chi1_chart[1] = modchi1 * sqrt1mxi12 * sin(alpha1);
    chi1_chart[2] = modchi1 * xi1;
    chi2_chart[0] = modchi2 * sqrt1mxi22 * cos(alpha2);
    chi2_chart[1] = modchi2 * sqrt1mxi22 * sin(alpha2);
    chi2_chart[2] = modchi2 * xi2;

    RotateVecT(pars->R1, chi1_chart, chi1);
    RotateVecT(pars->R2, chi2_chart, chi2);

    // Derivatives of the Hamiltonian (in lab frame)
    double Heff, H, dHeff[12], dH[12];
        Hamiltonian(r, p, nu, chi1, chi2, &Heff, &H, dHeff, dH);

    // Spin derivatives: extract from the lab-frame Hamiltonian output and
    // rotate into the chart frame (dH/dchi' = R dH/dchi) before applying the
    // chain rule for (alpha, xi).
    double dHdchi1_lab[3], dHdchi2_lab[3];
    double dHdchi1[3], dHdchi2[3]; // chart frame
    for (int i = 0; i < 3; i++)
    {
        dHdchi1_lab[i] = dH[6 + i];
        dHdchi2_lab[i] = dH[9 + i];
    }
    RotateVec(pars->R1, dHdchi1_lab, dHdchi1);
    RotateVec(pars->R2, dHdchi2_lab, dHdchi2);

    // *******
    // * rhs *
    // *******

    for (int i = 0; i < 3; i++)
    {
        dW[i]     = dH[i + 3];        // x, y, z
        dW[i + 3] = - dH[i];          // px, py, pz
    }

    // Mass fractions: the canonical momentum conjugate to alpha_i is the physical
    // spin projection S_{i,z} = m_i^2 |chi_i| xi_i, so both equations carry 1/m_i^2 = 1/X_i^2 (M = 1)
    double X1, X2;
        get_X1X2(nu, &X1, &X2);
        double iX1sq = nu/(X1*X1); // prefactor mu/m_i^2 = nu/X_i^2 (M=1); see Sec. on spin EOMs
        double iX2sq = nu/(X2*X2);

    const double eps = 1.0e-30;
    dW[6] = iX1sq*(- xi1*cos(alpha1)*dHdchi1[0] / (eps + sqrt1mxi12) - xi1*sin(alpha1)*dHdchi1[1] / (eps + sqrt1mxi12) + dHdchi1[2]); // dalpha1dt = (1/(m1^2|chi1|)) dH/dxi1
    dW[7] = iX1sq*(- sqrt1mxi12 * (dHdchi1[1]*cos(alpha1) - dHdchi1[0]*sin(alpha1))); // dxi1dt = - (1/(m1^2|chi1|)) dH/dalpha1
    dW[8] = iX2sq*(- xi2*cos(alpha2)*dHdchi2[0] / (eps + sqrt1mxi22) - xi2*sin(alpha2)*dHdchi2[1] / (eps + sqrt1mxi22) + dHdchi2[2]); // dalpha2dt = (1/(m2^2|chi2|)) dH/dxi2
    dW[9] = iX2sq*(- sqrt1mxi22 * (dHdchi2[1]*cos(alpha2) - dHdchi2[0]*sin(alpha2))); // dxi2dt = - (1/(m2^2|chi2|)) dH/dalpha2

    // Chart-pole proximity warning: in chart frame |xi| = |chi_z'/|chi||, so
    // it IS the ratio that must stay below 1. The chart is well-conditioned
    // up to about |xi| ~ 0.99; above 0.95 we warn so the user knows the
    // initial rotation may need to be revisited (or chart switching enabled).
    // Throttled to fire only on a 0.01 increase of the highwater mark.
    {
        static double max1 = 0., max2 = 0.;
        double r1 = fabs(xi1), r2 = fabs(xi2);
        if (r1 > 0.95 && r1 > max1 + 0.01) {
            fprintf(stderr, "Warning [get_rhs_canonical]: |xi1|/|chi1| = %.6f (chart pole approached, t=%.6f)\n", r1, t);
            max1 = r1;
        }
        if (r2 > 0.95 && r2 > max2 + 0.01) {
            fprintf(stderr, "Warning [get_rhs_canonical]: |xi2|/|chi2| = %.6f (chart pole approached, t=%.6f)\n", r2, t);
            max2 = r2;
        }
    }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                        RIGHT-HAND SIDE FUNCTION: Time transformation                               //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Computes the RHS of the Hamiltonian equations
// [in]  t              time
// [in]  w              Master variable (here including w[12] = t and w[13] = pt = - H0)
// [out] dw             Derivatives of w wrt all spacetime variables
void get_rhs_transformed(double s, double w[], double dw[])
{
    // Evolving the system described by the Hamiltonian:
    // K = r * (H + pt) with pt = -H0
    // Corresponds to a time transformation dt = r*ds that scales the time step wisely for eccentric orbits.
    
    double nu = pars->nu;

    // Unwrap Coordinates
    double r[3], p[3], chi1[3], chi2[3];
    for (int i = 0; i < 3; i++) {
        r[i]     = w[i];
        p[i]     = w[i + 3];
        chi1[i]  = w[i + 6];
        chi2[i]  = w[i + 9];
    }
    
    double pt = w[13]; // The energy constraint variable

    // Modulus of r and unit radial vector
    double modr = get_mod(r);
    double n[3];
        get_n(r, n, NULL);

    // Hamiltonian and derivatives
    double Heff, H, dHeff[12], dH[12];
    Hamiltonian(r, p, nu, chi1, chi2, &Heff, &H, dHeff, dH);

    // Extract spin derivatives
    double dHdchi1[3], dHdchi2[3];
    for (int i = 0; i < 3; i++) {
        dHdchi1[i] = dH[6 + i];
        dHdchi2[i] = dH[9 + i];
    }

    // Cross products for spin EOM
    double dHdchi1xchi1[3], dHdchi2xchi2[3];
    cross(dHdchi1, chi1, dHdchi1xchi1);
    cross(dHdchi2, chi2, dHdchi2xchi2);

    // Spin EOM: dchi_i/dt = (1/m_i^2)(dH_phys/dchi_i) x chi_i. The code uses Hhat = H_phys/mu
    // and dimensionless time t = T/M, so the prefactor on (dHhat/dchi_i) is mu/m_i^2 = nu/X_i^2
    // (M=1) -- NOT 1/X_i^2 (see the spin-EOM discussion in the paper).
    double X1, X2;
        get_X1X2(nu, &X1, &X2);
        double iX1sq = nu/(X1*X1); // prefactor mu/m_i^2 = nu/X_i^2 (M=1); see Sec. on spin EOMs
        double iX2sq = nu/(X2*X2);

    // Transformation dt = g ds with g = |r|
    double g, dgdr[3];
    // g = modr; 
    // for (int i = 0; i < 3; i++) dgdr[i] = n[i];

    // Transformation dt = g ds with g = |r|^3, better
    g = modr*modr*modr;
    for (int i = 0; i < 3; i++) dgdr[i] = 3.*modr*modr * n[i];

    // *******
    // * rhs *
    // *******

    // dq/ds = g * (dH/dp)
    // dp/ds = - [ g * (dH/dq) + (H + pt) * (dg/dr) ]
    // dchi/ds = g * {chi, H}
    for (int i = 0; i < 3; i++)
    {
        dw[i]     = g * dH[i + 3];
        dw[i + 3] = - (g * dH[i] + (H + pt) * dgdr[i]); 
        dw[i + 6] = g * iX1sq*dHdchi1xchi1[i];
        dw[i + 9] = g * iX2sq*dHdchi2xchi2[i];
    }

    // dt/ds = dK/dpt = g
    dw[12] = g;

    // dpt/ds = -dK/dt = 0
    dw[13] = 0.;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//          RIGHT-HAND SIDE FUNCTION: Orbital block of the Strang split (spins frozen)               //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// This is the "A" block of the operator split (see StrangStepTransformed). It advances the orbital
// extended phase space (r, p, t, pt) under the transformed Hamiltonian K = g*(H + pt), holding the
// spin VECTORS w[6..11] fixed (their derivatives are set to zero). The spin "B" block then rotates
// the spins exactly (Rodrigues). Spins are still read and fed to the Hamiltonian here, because the
// orbital potentials depend on them through a0.
// Layout: w[0..2]=r, w[3..5]=p, w[6..8]=chi1, w[9..11]=chi2, w[12]=t, w[13]=pt(=-H0).
void get_rhs_orbital_transformed(double s, double w[], double dw[])
{
    double nu = pars->nu;

    // Unwrap coordinates (spins are parameters here)
    double r[3], p[3], chi1[3], chi2[3];
    for (int i = 0; i < 3; i++) {
        r[i]    = w[i];
        p[i]    = w[i + 3];
        chi1[i] = w[i + 6];
        chi2[i] = w[i + 9];
    }

    double pt = w[13];

    // Modulus of r and unit radial vector
    double modr = get_mod(r);
    double n[3];
        get_n(r, n, NULL);

    // Hamiltonian and derivatives
    double Heff, H, dHeff[12], dH[12];
        Hamiltonian(r, p, nu, chi1, chi2, &Heff, &H, dHeff, dH);

    // Time transformation dt = g ds with g = |r|^3
    double g = modr*modr*modr;
    double dgdr[3];
    for (int i = 0; i < 3; i++) dgdr[i] = 3.*modr*modr * n[i];

    // *******
    // * rhs * (orbital variables only; spins frozen)
    // *******
    // NOTE: this uses the correct g*dH[i] in the dp/ds term. The older
    // get_rhs_transformed has a latent typo there (it uses modr instead of g);
    // see the message accompanying this change.
    for (int i = 0; i < 3; i++)
    {
        dw[i]     = g * dH[i + 3];
        dw[i + 3] = - (g * dH[i] + (H + pt) * dgdr[i]);
        dw[i + 6] = 0.;   // chi1 frozen in the orbital block
        dw[i + 9] = 0.;   // chi2 frozen in the orbital block
    }

    // dt/ds = dK/dpt = g
    dw[12] = g;
    // dpt/ds = -dK/dt = 0
    dw[13] = 0.;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                    SPIN BLOCK of the Strang split: exact precession (Rodrigues)                   //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Helper: precession axis vectors Omega_i = (1/X_i^2) dH/dchi_i (so dchi_i/dt = Omega_i x chi_i),
// evaluated at the (frozen) orbital state r, p and the given spins.
static void SpinOmega(double r[], double p[], double chi1[], double chi2[],
                      double Om1[], double Om2[])
{
    double nu = pars->nu;
    double Heff, H, dHeff[12], dH[12];
        Hamiltonian(r, p, nu, chi1, chi2, &Heff, &H, dHeff, dH);
    double X1, X2;
        get_X1X2(nu, &X1, &X2);
    double iX1sq = nu/(X1*X1); // prefactor mu/m_i^2 = nu/X_i^2 (M=1); see Sec. on spin EOMs
    double iX2sq = nu/(X2*X2);
    for (int i = 0; i < 3; i++) {
        Om1[i] = iX1sq * dH[6 + i];
        Om2[i] = iX2sq * dH[9 + i];
    }
}

// Helper: rotate v about the axis Omega by the angle |Omega|*dt (exact, Rodrigues).
static void RotateAboutOmega(double v[], double Om[], double dt, double out[])
{
    double nOm = get_mod(Om);
    if (nOm > 0.) {
        double k[3] = { Om[0]/nOm, Om[1]/nOm, Om[2]/nOm };
        RodriguesRotate(v, k, nOm * dt, out);
    } else {
        for (int i = 0; i < 3; i++) out[i] = v[i];
    }
}

// This is the "B" block. The orbital variables (r, p, t, pt) are held FIXED, so each spin obeys a
// constant-axis precession dchi_i/ds = g * Omega_i x chi_i, with the precession axis vector
//   Omega_i = (1/X_i^2) dH/dchi_i   (the same combination as the cross-product EOM in get_rhs).
// With r frozen, g is constant, so each spin's exact flow over the step is a rigid rotation about
// Omega_i by angle |Omega_i| * g * h = |Omega_i| * (physical time of step), done exactly by
// Rodrigues for ANY angle -- so the stiff 1/X_2^2 secondary precession is absorbed with no
// step-size penalty and |chi_i| is conserved to roundoff. No spin chart, hence no pole.
//
// MIDPOINT REFINEMENT (time-reversible): Omega_i depends weakly on the spins themselves (spin-spin
// coupling via a0, S, S*). Instead of freezing Omega_i at the incoming spins, we evaluate it at the
// HALF-ROTATION midpoint and rotate by the full angle about that axis. A fixed-point iteration
// brings Omega_i to self-consistency at the midpoint:
//     chi_i_mid = R(Omega_i, |Omega_i| g h / 2) chi_i_in ,   Omega_i = Omega_i(chi_1_mid, chi_2_mid)
// Because the half-rotation point is the same whether the step is taken forward or backward, the
// resulting map satisfies B(h) o B(-h) = id (time-reversible), which removes the leading spin-spin
// splitting error while keeping the rotation (and |chi_i|) exact. Convergence is fast: the contraction
// factor scales with (step)*(spin-spin coupling), independent of the fast 1/X_2^2 precession, since
// Omega is essentially insensitive to the (fast, tiny) secondary spin direction.
void SpinPrecessionStep(double w[], double h)
{
    double r[3], p[3], chi1[3], chi2[3];
    for (int i = 0; i < 3; i++) {
        r[i]    = w[i];
        p[i]    = w[i + 3];
        chi1[i] = w[i + 6];
        chi2[i] = w[i + 9];
    }

    // Same time-transformation factor g = |r|^3 as the orbital block; dtau = physical time of step.
    double modr = get_mod(r);
    double g    = modr*modr*modr;
    double dtau = g * h;

    // Initial guess: Omega at the incoming spins.
    double Om1[3], Om2[3];
        SpinOmega(r, p, chi1, chi2, Om1, Om2);

    // Fixed-point iteration for the self-consistent midpoint precession axes.
    const int    maxit = (int) pars->max_iter_RKGL6;
    const double tol   = pars->tol_RKGL6;
    for (int it = 0; it < maxit; it++) {
        double chi1m[3], chi2m[3], Om1n[3], Om2n[3];
        RotateAboutOmega(chi1, Om1, 0.5*dtau, chi1m); // half-rotation midpoint spins
        RotateAboutOmega(chi2, Om2, 0.5*dtau, chi2m);
        SpinOmega(r, p, chi1m, chi2m, Om1n, Om2n);    // re-evaluate Omega at the midpoint

        // Relative convergence on the precession axes (handles the |Om1|~O(1) vs |Om2|~1/X_2 scales)
        double d1 = 0., d2 = 0., n1 = 0., n2 = 0.;
        for (int i = 0; i < 3; i++) {
            d1 += (Om1n[i]-Om1[i])*(Om1n[i]-Om1[i]); n1 += Om1n[i]*Om1n[i];
            d2 += (Om2n[i]-Om2[i])*(Om2n[i]-Om2[i]); n2 += Om2n[i]*Om2n[i];
        }
        for (int i = 0; i < 3; i++) { Om1[i] = Om1n[i]; Om2[i] = Om2n[i]; }
        if (sqrt(d1) <= tol*(sqrt(n1)+1.0) && sqrt(d2) <= tol*(sqrt(n2)+1.0)) break;
    }

    // Final full rotation about the converged midpoint axes.
    double out1[3], out2[3];
    RotateAboutOmega(chi1, Om1, dtau, out1);
    RotateAboutOmega(chi2, Om2, dtau, out2);
    for (int i = 0; i < 3; i++) { w[6 + i] = out1[i]; w[9 + i] = out2[i]; }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                    STRANG-SPLIT STEP: A(h/2) o B(h) o A(h/2)  (second order, symmetric)           //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// One full step of the operator-split integrator on the extended phase space w[14] (spins as
// vectors). The orbital half-steps use the symplectic RKGL6 collocation solver on the orbital-only
// RHS (spins frozen); the spin step rotates the spins exactly (Rodrigues). The composition is
// symmetric/time-reversible and preserves |chi_i|; it is the non-canonically symplectic scheme of
// Lubich-Walther-Bruegmann (PRD 81, 104025) adapted to this EOB Hamiltonian.
void StrangStepTransformed(double s, double w[], double h)
{
    RK_GaussLegendre6(s,          w, 0.5*h, 14, get_rhs_orbital_transformed); // A(h/2)
    SpinPrecessionStep(w, h);                                                 // B(h)
    RK_GaussLegendre6(s + 0.5*h,  w, 0.5*h, 14, get_rhs_orbital_transformed); // A(h/2)
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                             POINCARE SECTION OUTPUT (z = 0 upward crossing)                       //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Detects a crossing of the surface of section z = 0 with z increasing (zp <= 0 < zc) between the
// previous (subscript p) and current (subscript c) states, linearly interpolates the state to the
// crossing, and writes a section point. Spins are converted to the canonical pair
//   alpha_i = atan2(chi_i_y, chi_i_x),  S_{i,z} = X_i^2 * chi_i_z   (lab frame)
// at the crossing -- the pole is harmless here (it is a single atan2 evaluation, not an integration
// singularity). Spin vectors are renormalized to their (conserved) magnitude after interpolation.
// Change the section condition below to use a different plane/direction if desired.
void WriteSectionCrossing(FILE *fps,
                          double tp, double rp[], double pp[], double c1p[], double c2p[],
                          double tc, double rc[], double pc[], double c1c[], double c2c[])
{
    if (fps == NULL) return;
    if (!(rp[2] <= 0. && rc[2] > 0.)) return; // z = 0 upward crossing only

    double dz = rc[2] - rp[2];
    double f  = (dz != 0.) ? (-rp[2] / dz) : 0.; // fraction of the step to the crossing, in [0,1]

    // Linear interpolation of the lab-frame state to the crossing
    double t = tp + f*(tc - tp);
    double r[3], p[3], c1[3], c2[3];
    for (int i = 0; i < 3; i++) {
        r[i]  = rp[i] + f*(rc[i] - rp[i]);
        p[i]  = pp[i] + f*(pc[i] - pp[i]);
        c1[i] = c1p[i] + f*(c1c[i] - c1p[i]);
        c2[i] = c2p[i] + f*(c2c[i] - c2p[i]);
    }

    // Renormalize spins to their conserved magnitude (interpolation slightly shrinks the chord)
    double m1 = get_mod(c1p), m2 = get_mod(c2p);
    double n1 = get_mod(c1),  n2 = get_mod(c2);
    if (n1 > 0.) for (int i = 0; i < 3; i++) c1[i] *= m1/n1;
    if (n2 > 0.) for (int i = 0; i < 3; i++) c2[i] *= m2/n2;

    double X1, X2;
        get_X1X2(pars->nu, &X1, &X2);
    double alpha1 = atan2(c1[1], c1[0]);
    double alpha2 = atan2(c2[1], c2[0]);
    double S1z = X1*X1 * c1[2];
    double S2z = X2*X2 * c2[2];

    fprintf(fps, "%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n",
            t, r[0], r[1], p[0], p[1], p[2], alpha1, S1z, alpha2, S2z);
    fflush(fps);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                RIGHT-HAND SIDE FUNCTION: Canonical spins + Time transformation                     //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Computes the RHS of the Hamiltonian equations
// [in]  t              time
// [in]  w              Master variable (here including w[12] = t and w[13] = pt = - H0)
// [out] dw             Derivatives of w wrt all spacetime variables
void get_rhs_canonical_transformed(double s, double W[], double dW[])
{
    double nu = pars->nu;

    // Unwrap
    double r[3], p[3], alpha1, xi1, alpha2, xi2;
    for (int i = 0; i < 3; i++) {
        r[i]     = W[i];
        p[i]     = W[i + 3];
    }
    alpha1 = W[6];
    xi1    = W[7];
    alpha2 = W[8];
    xi2    = W[9];
    
    // FIXME: this could be changed and updated during the evolution by defining pars->modchi and updating it in the dynamics.
    // However even with the non-canonical evolution, the error on the conservation of the spin magnitude is of the order of 10^-14.
    double modchi1 = sqrt((pars->chi1x0)*(pars->chi1x0) + (pars->chi1y0)*(pars->chi1y0) + (pars->chi1z0)*(pars->chi1z0));
    double modchi2 = sqrt((pars->chi2x0)*(pars->chi2x0) + (pars->chi2y0)*(pars->chi2y0) + (pars->chi2z0)*(pars->chi2z0));

    // Reconstruct the spins in the (rotated) chart frame, then rotate to lab
    // for the Hamiltonian call. (alpha_i, xi_i) live in the chart frame.
    double chi1_chart[3], chi2_chart[3], chi1[3], chi2[3];
    double sqrt1mxi12 = sqrt(1. - xi1*xi1);
    double sqrt1mxi22 = sqrt(1. - xi2*xi2);

    chi1_chart[0] = modchi1 * sqrt1mxi12 * cos(alpha1);
    chi1_chart[1] = modchi1 * sqrt1mxi12 * sin(alpha1);
    chi1_chart[2] = modchi1 * xi1;
    chi2_chart[0] = modchi2 * sqrt1mxi22 * cos(alpha2);
    chi2_chart[1] = modchi2 * sqrt1mxi22 * sin(alpha2);
    chi2_chart[2] = modchi2 * xi2;

    RotateVecT(pars->R1, chi1_chart, chi1);
    RotateVecT(pars->R2, chi2_chart, chi2);

    double pt = W[11]; // = -H0

    // Modulus of r and unit radial vector
    double modr = get_mod(r);
    double n[3];
        get_n(r, n, NULL);

    // Hamiltonian and derivatives (lab frame)
    double Heff, H, dHeff[12], dH[12];
        Hamiltonian(r, p, nu, chi1, chi2, &Heff, &H, dHeff, dH);

    // Spin derivatives: lab frame -> chart frame (dH/dchi' = R dH/dchi)
    double dHdchi1_lab[3], dHdchi2_lab[3];
    double dHdchi1[3], dHdchi2[3]; // chart frame
    for (int i = 0; i < 3; i++)
    {
        dHdchi1_lab[i] = dH[6 + i];
        dHdchi2_lab[i] = dH[9 + i];
    }
    RotateVec(pars->R1, dHdchi1_lab, dHdchi1);
    RotateVec(pars->R2, dHdchi2_lab, dHdchi2);

    // Transformation dt = g ds with g = |r|
    double g, dgdr[3];
    // g = modr; 
    // for (int i = 0; i < 3; i++) dgdr[i] = n[i];

    // Transformation dt = g ds with g = |r|^3, better
    g = modr*modr*modr;
    for (int i = 0; i < 3; i++) dgdr[i] = 3.*modr*modr * n[i];

    // *******
    // * rhs *
    // *******

    for (int i = 0; i < 3; i++)
    {
        dW[i]     = g * dH[i + 3]; 
        dW[i + 3] = - (g * dH[i] + (H + pt) * dgdr[i]);
    }

    // Mass fractions: the canonical momentum conjugate to alpha_i is the physical
    // spin projection S_{i,z} = m_i^2 |chi_i| xi_i, so both equations carry 1/m_i^2 = 1/X_i^2 (M = 1)
    double X1, X2;
        get_X1X2(nu, &X1, &X2);
        double iX1sq = nu/(X1*X1); // prefactor mu/m_i^2 = nu/X_i^2 (M=1); see Sec. on spin EOMs
        double iX2sq = nu/(X2*X2);

    const double eps = 1.0e-30;

    double dalpha1dt = iX1sq*(- xi1*cos(alpha1)*dHdchi1[0] / (eps + sqrt1mxi12) - xi1*sin(alpha1)*dHdchi1[1] / (eps + sqrt1mxi12) + dHdchi1[2]); // dalpha1dt = (1/(m1^2|chi1|)) dH/dxi1
    double dxi1dt = iX1sq*(- sqrt1mxi12 * (dHdchi1[1]*cos(alpha1) - dHdchi1[0]*sin(alpha1))); // dxi1dt = - (1/(m1^2|chi1|)) dH/dalpha1
    double dalpha2dt = iX2sq*(- xi2*cos(alpha2)*dHdchi2[0] / (eps + sqrt1mxi22) - xi2*sin(alpha2)*dHdchi2[1] / (eps + sqrt1mxi22) + dHdchi2[2]); // dalpha2dt = (1/(m2^2|chi2|)) dH/dxi2
    double dxi2dt = iX2sq*(- sqrt1mxi22 * (dHdchi2[1]*cos(alpha2) - dHdchi2[0]*sin(alpha2))); // dxi2dt = - (1/(m2^2|chi2|)) dH/dalpha2

    // Spin variables
    dW[6] = g * dalpha1dt;
    dW[7] = g * dxi1dt;
    dW[8] = g * dalpha2dt;
    dW[9] = g * dxi2dt;

    // dt/ds = dK/dpt = r
    dW[10] = g;

    // dpt/ds = -dK/dt = 0
    dW[11] = 0.;

    // Chart-pole proximity warning (see get_rhs_canonical for explanation).
    // Use the physical time stored in W[10] rather than s.
    {
        static double max1 = 0., max2 = 0.;
        double r1 = fabs(xi1), r2 = fabs(xi2);
        if (r1 > 0.95 && r1 > max1 + 0.01) {
            fprintf(stderr, "Warning [get_rhs_canonical_transformed]: |xi1|/|chi1| = %.6f (chart pole approached, t=%.6f)\n", r1, W[10]);
            max1 = r1;
        }
        if (r2 > 0.95 && r2 > max2 + 0.01) {
            fprintf(stderr, "Warning [get_rhs_canonical_transformed]: |xi2|/|chi2| = %.6f (chart pole approached, t=%.6f)\n", r2, W[10]);
            max2 = r2;
        }
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                                      CARTER-LIKE CONSTANT                                          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Computes the Carter-Like constant described in Damour 2001. Note: currently unused.
// [in]  r              Radial vector
// [in]  p              Momentum vector
// [in]  chi1           Spin vector of the first particle
// [in]  chi2           Spin vector of the second particle
double CarterLikeConstant(double r[], double p[], double chi1[], double chi2[]) 
{
    double nu = pars->nu;

    // Unwrap coordinates & momenta
    double x  = r[0];
    double y  = r[1];
    double z  = r[2];
    double px = p[0];
    double py = p[1];
    double pz = p[2];
    
    // Get spherical coordinates & momenta
    double R   = get_mod(r);
    // double phi = atan(y/x); // check if this is okay
    // double phi = atan(fabs(y)/fabs(x));
    // double sinphi = sin(phi);
    // double cosphi = cos(phi);

    double costheta  = z/R; 
    double costheta2 = costheta*costheta; 
    double sintheta2 = 1. - costheta2;
    double sintheta  = sqrt(sintheta2);
    
    // double ptheta  = R*(- pz * sintheta + costheta * (px * cosphi + py * sinphi));
    double ptheta  = - R * pz * sintheta + px * (x*z/sintheta) + py * (y*z/sintheta); 
    double ptheta2 = ptheta*ptheta;

    // double pphi  = R * sintheta * (py * cosphi - px * sinphi);
    double pphi  = x*py - y*px;
    double pphi2 = pphi*pphi;

    double a0[3];
        get_a0(nu, chi1, chi2, a0);
    double a02 = dot(a0, a0, 3);

    double Heff, H;
    Hamiltonian(r, p, nu, chi1, chi2, &Heff, &H, NULL, NULL);
    double Heff2 = Heff*Heff;
    
    return ptheta2 + costheta2 * (pphi2 / sintheta2 + a02 * (1. - Heff2));
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                                      PHOTON POTENTIAL CONDITION                                    //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Returns the derivative of the photon potential u^2*Aorb with respect to u;
// Used to find the maximum of the photon potential, location of the adiabatic light-ring.
// [in]  modr           Modulus of the radial vector
double PhotonPotentialCondition(double modr)
{
    double nu = pars->nu;

    double Aorb, dAorb;
        get_Aorb(modr, nu, &Aorb, &dAorb);

    double dAorbdu = - modr*modr * dAorb;

    return 2.*Aorb + (1./modr) * dAorbdu;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                                      PHOTON POTENTIAL CONDITION                                    //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Evaluates the adiabatic light-ring as the maximum of the photon potential u^2*Aorb.
// [out] rLR    Adiabatic light-ring radius
void AdiabaticLightRing(double *rLR)
{
    double xa = 1.8, xb = 3.1, zero;
    Secant_v2(PhotonPotentialCondition, xa, xb, 1e-10, &zero);

    *rLR = zero;
}