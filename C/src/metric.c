/*
File: metric.c

This file contains functions for computing non-resummed and resummed potentials, next-to-leading-order (NLO)
spin-spin contributions, and gyro-gravitomagnetic functions.
*/

#include "header.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                                       NONRESUMMED POTENTIALS                                       //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// [in]  r              Radial vector
// [in]  p              Momentum vector
// [in]  nu             Symmetric mass ratio
// [in]  chi1           Spin vector of the first particle
// [in]  chi2           Spin vector of the second particle
// [out] A              NLO-spin-spin potential, denoted as A
// [out] Bp             NLO-spin-spin potential, denoted as Bp
// [out] Bnp            NLO-spin-spin potential, denoted as Bnp
// [out] dA             Derivative of potential A wrt all spacetime variables
// [out] dBp            Derivative of potential Bp wrt all spacetime variables
// [out] dBnp           Derivative of potential Bnp wrt all spacetime variables
void NonResummedPotentials(double r[], double p[], double nu, double chi1[], double chi2[],
                           double *A, double *Bp, double *Bnp, double *Benp, double *dA, double *dBp, double *dBnp, double *dBenp)
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

    // Effective spin vector
    double a0[3];
        get_a0(nu, chi1, chi2, a0);

    // Delta function and derivatives
    double Delta, dDelta[4][3];
        get_Delta(r, nu, chi1, chi2, &Delta, (double *)dDelta);

    // NLO-spin-spin contributions and derivatives
    // AA: TODO: change here [12] -> [4][3]
    double AchiQ, AnchiQ, BchiQ, BnchiQ, dAchiQ[12], dAnchiQ[12], dBchiQ[12], dBnchiQ[12];
        NLO_SpinSpinContributions(r, nu, chi1, chi2, &AchiQ, &AnchiQ, &BchiQ, &BnchiQ, dAchiQ, dAnchiQ, dBchiQ, dBnchiQ);

    // Equatorial "bare" A potential
    double AeqB, dAeqB[4][3];
        get_AeqB(r, nu, chi1, chi2, &AeqB, (double *)dAeqB);

    // D potential as a function of rc and derivative wrt to rc
    double Dorb, dDorbdrc;
        get_Dorb(rc, nu, &Dorb, &dDorbdrc);

    // Derivatives of Dorb wrt to all the variables
    double dDorb[4][3];
    for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 3; i++) {
            dDorb[j][i] = drc[j][i]*dDorbdrc;
        }
    }
    
    // Auxiliary expressions
    double na0       = dot(n, a0, 3);       // Dot product of n and a0
    double over_na0;
    if (na0 < 1e-15) {
        over_na0 = 0.;
    }
    else {
        over_na0  = 1. / na0; // Inverse of na0
    }
    
    double na02      = na0*na0;             // Squared dot product of n and a0
    double u         = 1. / get_mod(r);     // Inverse of |r|
    double u2        = u * u;               // Inverse of |r|^2
    double u3        = u2 * u;              // Inverse of |r|^3
    double u4        = u3 * u;              // Inverse of |r|^4
    double uc        = 1. / rc;             // Inverse of rc
    double over_Dorb = 1. / Dorb;           // Inverse of Dorb

    // **************
    // * Potentials *
    // **************

    // Fragments of the potentials
    double A1 = 1. + na02 * u2;
    double A2 = 1. / (1. + Delta * na02 * u2 * uc * uc);
    double B1 = -1. + AeqB * rc * rc * u2 * over_Dorb;

    // Kerr-like EOB potentials
    double AnuK0   = AeqB * A1 * A2;
    double BpnuK0  = 1. / A1;
    double BnpnuK0 = BpnuK0 * B1;

    // NLO-spin-spin potentials
    *A   = AnuK0   + (AchiQ - AnchiQ) * u4;
    *Bp  = BpnuK0  - BnchiQ * u3;
    *Bnp = BnpnuK0 + BchiQ  * u3;
    *Benp = 0.; // FIXME (will be implemented later)

    // *********************************
    // * Derivatives of the potentials *
    // *********************************

    if (dA != NULL & dBp != NULL & dBnp != NULL)
    {
        // General auxiliary expressions
        double modified_A1 = 2.*(A1 - 1.);       // Expression repeated a few times in the code
        double over_A1     = 1. / A1;            // Inverse of A1
        double modified_A2 = 2.*(A2 - 1.);       // Expression repeated a few times in the code
        double modified_B1 = 1. + B1;            // Expression repeated a few times in the code
        double over_Delta  = 1. / (2. * Delta);  // Inverse of 2xDelta
        double over_Aeqb   = 1. / AeqB;          // Inverse of AeqB
        double BpnuK02     = BpnuK0 * BpnuK0;    // Squared BpnuK0

        for (int i = 0; i < 3; i++)
        {
            // Auxiliary expressions for spatial derivatives
            double r_over_r2       = r[i] * u2;                                 // Expression for r[i]/(|r|^2)
            double r_over_r5       = n[i] * u4;                                 // Expression for r[i]/(|r|^5)
            double drc_over_rc     = drc[0][i] * uc;                            // Expression for drc[0][i]/rc
            double dAeqB_over_AeqB = dAeqB[0][i] * over_Aeqb;                   // Expression for dAeqB[0][i]/AeqB
            double da1dr           = dot(dn[i], a0, 3) * over_na0 - r_over_r2;  // Expression repeated a few times in the code

            // Fragments of spatial derivative
            double dA1dr = modified_A1 * da1dr;
            double dA2dr = modified_A2 * (dDelta[0][i] * over_Delta + da1dr - drc_over_rc);
            double dB1dr = modified_B1 * (dAeqB_over_AeqB - dDorb[0][i] * over_Dorb + 2.*(drc_over_rc - r_over_r2));

            // Spatial derivatives
            double dAnuK0dr   = AnuK0 * (dAeqB_over_AeqB + dA1dr * over_A1 + dA2dr);
            double dBpnuK0dr  = -BpnuK02 * dA1dr;
            double dBnpnuK0dr = B1 * dBpnuK0dr + dB1dr * BpnuK0;

            dA[i]    = dAnuK0dr + (dAchiQ[i] - dAnchiQ[i]) * u4 - 4.*u * r_over_r5 * (AchiQ - AnchiQ);
            dBp[i]   = dBpnuK0dr - dBnchiQ[i] * u3 + 3.*r_over_r5 * BnchiQ;
            dBnp[i]  = dBnpnuK0dr + dBchiQ[i] * u3 - 3.*r_over_r5 * BchiQ;
            dBenp[i] = 0.; // FIXME (will be implemented later)

            for(int j = 2; j < 4; j++)
            {
                // Auxiliary expressions for spin derivatives
                drc_over_rc     = drc[j][i] * uc;            // Expression for drc[j][i]/rc
                dAeqB_over_AeqB = dAeqB[j][i] * over_Aeqb;   // Expression for dAeqB[j][i]/AeqB
                double da1dchi  = over_na0 * n[i] * X[j-2];  // Expression repeated a few times in the code
                // FIXME: the line above might be wrong

                // Fragments of spin derivative
                double dA1dchi = modified_A1 * da1dchi;
                double dA2dchi = modified_A2 * (dDelta[j][i] * over_Delta + da1dchi - drc_over_rc);
                double dB1dchi = modified_B1 * (dAeqB_over_AeqB - dDorb[j][i] * over_Dorb + 2.*drc_over_rc);

                // Spin derivatives
                double dAnuK0dchi   = AnuK0 * (dAeqB_over_AeqB + dA1dchi * over_A1 + dA2dchi);
                double dBpnuK0dchi  = -BpnuK02 * dA1dchi;
                double dBnpnuK0dchi = B1 * dBpnuK0dchi + dB1dchi * BpnuK0;

                dA[i+6 + 3*(j-2)]    = dAnuK0dchi + (dAchiQ[i+6 + 3*(j-2)] - dAnchiQ[i+6 + 3*(j-2)]) * u4;
                dBp[i+6 + 3*(j-2)]   = dBpnuK0dchi - dBnchiQ[i+6 + 3*(j-2)] * u3;
                dBnp[i+6 + 3*(j-2)]  = dBnpnuK0dchi + dBchiQ[i+6 + 3*(j-2)] * u3;
                dBenp[i+6 + 3*(j-2)] = 0.; // FIXME (will be implemented later)
            }
        }
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                                         RESUMMED POTENTIALS                                        //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// [in]  r              Radial vector
// [in]  p              Momentum vector
// [in]  nu             Symmetric mass ratio
// [in]  chi1           Spin vector of the first particle
// [in]  chi2           Spin vector of the second particle
// [out] A              NLO-spin-spin potential, denoted as A
// [out] Bp             NLO-spin-spin potential, denoted as Bp
// [out] Bnp            NLO-spin-spin potential, denoted as Bnp
// [out] dA             Derivative of potential A wrt all spacetime variables
// [out] dBp            Derivative of potential Bp wrt all spacetime variables
// [out] dBnp           Derivative of potential Bnp wrt all spacetime variables
void ResummedPotentials(double r[], double p[], double nu, double chi1[], double chi2[], 
                        double *A, double *Bp, double *Bnp, double *Benp, double *dA, double *dBp, double *dBnp, double *dBenp)
{

    double modr = get_mod(r); 
    double r2   = modr*modr;
    double r3   = r2*modr;
    double r4   = r2*r2;
    double u    = 1./modr;
    double u2   = u*u;
    double u3   = u2*u;
    double u4   = u2*u2;

    double n[3], dndr[3][3]; 
        get_n(r, n, (double *)dndr);

    double Delta, dDelta[4][3] = {0.}; 
        get_Delta(r, nu, chi1, chi2, &Delta, (double *)dDelta);

    double a0[3];
        get_a0(nu, chi1, chi2, a0);
    double a02 = dot(a0, a0, 3);

    // Get centrifugal radius & derivatives 
    double rc, drc[4][3] = {0.};
        get_rc(r, nu, chi1, chi2, &rc, (double *)drc);
    double rc2      = rc*rc;
    double uc       = 1./rc;
    double uc2      = uc*uc; 
    double drcdmodr = (modr - a02*u2)*uc;
    double ducdu    = (uc2*r2)*drcdmodr;

    // R4 function and derivatives
    double R4, dR4[4][3] = {0.};
        get_R4(r, nu, chi1, chi2, &R4, (double *)dR4);
    
    // Projection of the effective Kerr parameter along the radial direction
    double a0n    = dot(n, a0, 3); 
    double a0n2   = a0n*a0n;
    double a0n2u2 = a0n2*u2;

    double X1, X2;
        get_X1X2(nu, &X1, &X2);
    
    double AchiQ, AnchiQ, BchiQ, BnchiQ, dAchiQ[4][3], dAnchiQ[4][3], dBchiQ[4][3], dBnchiQ[4][3];
        NLO_SpinSpinContributions(r, nu, chi1, chi2, &AchiQ, &AnchiQ, &BchiQ, &BnchiQ, (double *)dAchiQ, (double *)dAnchiQ, (double *)dBchiQ, (double *)dBnchiQ);

    // Note: Eq. 2.80 does not specify that Dorb is a funciton of rc, but it most likely is. 
    // Compare with 2.20 & 2.26 for the non-resummed version of the potential Bnp.
    double Dorb, dDorbdrc;
        get_Dorb(rc, nu, &Dorb, &dDorbdrc); 

    double Aeq, dAeq[4][3];
        get_Aeq(r, nu, chi1, chi2, &Aeq, (double *)dAeq);

    double AeqB, dAeqB[4][3];
        get_AeqB(r, nu, chi1, chi2, &AeqB, (double *)dAeqB);

    double numA      = (1. + a0n2u2 - AnchiQ * u4);
    double denA      = (1. + (Delta * a0n2u2) * uc2);
    double denBp     = 1. + a0n2u2 + BnchiQ * u3;
    double BnpI      = 1. / (1. + a0n2u2);
    double numBnpII  = AeqB * (rc2 + u*BchiQ);
    double denBnpII  = Dorb*r2;
    double BnpII     = numBnpII/denBnpII - 1.;
    double BenpI     = - 1. / (1. + a0n2u2); // = - BnpI
    double numBenpII = r2 + 2.*modr + a0n2;
    double denBenpII = R4 + Delta*a0n2;
    double BenpII    = numBenpII / denBenpII;
    
    *A    = Aeq * numA/denA; // Defined by eq. 2.74
    *Bp   = 1. / denBp; // Defined by eq. 2.75
    *Bnp  = BnpI * BnpII; // Defined by eq. 2.80
    *Benp = BenpI * BenpII; // Not defined explicitly, but analogous to the Kerr contribution

    if (dA != NULL & dBp != NULL & dBnp != NULL & dBenp != NULL)
    {
        double da0ndr[3], da0n2dr[3], da0n2u2dr[3], da0n2dchi1[3], da0n2dchi2[3], da0n2u2dchi1[3], da0n2u2dchi2[3];
        double dnumAdr[3], ddenAdr[3], dnumAdchi1[3], ddenAdchi1[3], dnumAdchi2[3], ddenAdchi2[3];
        double Bp2        = (*Bp)*(*Bp);
        double BnpI2      = BnpI*BnpI;
        double BenpI2     = BenpI*BenpI;
        double denBnpII2  = denBnpII*denBnpII;
        double denBenpII2 = denBenpII*denBenpII;
        double ddenBpdr[3], ddenBpdchi1[3], ddenBpdchi2[3];
        double dBnpIdr[3], dnumBnpIIdr[3], ddenBnpIIdr[3], dBnpIIdr[3];
        double dBnpIdchi1[3], dnumBnpIIdchi1[3], ddenBnpIIdchi1[3], dBnpIIdchi1[3];
        double dBnpIdchi2[3], dnumBnpIIdchi2[3], ddenBnpIIdchi2[3], dBnpIIdchi2[3];
        double dBenpIdr[3], dnumBenpIIdr[3], ddenBenpIIdr[3], dBenpIIdr[3];
        double dBenpIdchi1[3], dnumBenpIIdchi1[3], ddenBenpIIdchi1[3], dBenpIIdchi1[3];
        double dBenpIdchi2[3], dnumBenpIIdchi2[3], ddenBenpIIdchi2[3], dBenpIIdchi2[3];

        for (int i = 0; i < 3; i++)
        {
            // Derivatives wrt to x, y, z

            da0ndr[i]    = dot(a0, dndr[i], 3);
            da0n2dr[i]   = 2.*a0n*da0ndr[i];
            da0n2u2dr[i] = 2.*u2*a0n*( -a0n*u*n[i] + da0ndr[i] );

            dnumAdr[i] = da0n2u2dr[i] + u4*( 4.*u*n[i]*AnchiQ - dAnchiQ[0][i] );
            ddenAdr[i] = u2*uc*( uc * (dDelta[0][i]*a0n2 + Delta*da0n2dr[i]) + 2.*Delta*a0n2 * (-n[i]*u) * (uc + u*ducdu) );
            dA[i]      = dAeq[0][i] * numA/denA + Aeq * (dnumAdr[i]*denA - ddenAdr[i]*numA)/(denA*denA);

            ddenBpdr[i] = da0n2u2dr[i] + u3*( -3.*u*n[i]*BnchiQ + dBnchiQ[0][i] );
            dBp[i]      = - ddenBpdr[i]*Bp2;

            dBnpIdr[i]     = - BnpI2*da0n2u2dr[i];
            dnumBnpIIdr[i] = dAeqB[0][i] * (rc2 + u*BchiQ) + AeqB * (2.*rc*drc[0][i] - u2*n[i]*BchiQ + u*dBchiQ[0][i]);
            ddenBnpIIdr[i] = dDorbdrc*drc[0][i]*r2 + Dorb*2.*modr*n[i];
            dBnpIIdr[i]    = (dnumBnpIIdr[i]*denBnpII - ddenBnpIIdr[i]*numBnpII)/denBnpII2;
            dBnp[i]        = dBnpIdr[i]*BnpII + BnpI*dBnpIIdr[i];

            dBenpIdr[i]     = BenpI2 * da0n2u2dr[i];
            dnumBenpIIdr[i] = 2.*n[i]*(modr + 1.) + da0n2dr[i];
            ddenBenpIIdr[i] = dR4[0][i] + dDelta[0][i]*a0n2 + Delta*da0n2dr[i];
            dBenpIIdr[i]    = (dnumBenpIIdr[i]*denBenpII - ddenBenpIIdr[i]*numBenpII)/denBenpII2;
            dBenp[i]        = dBenpIdr[i]*BenpII + BenpI*dBenpIIdr[i];

            // Derivatives wrt to px, py, pz

            dA[i + 3]    = 0.;
            dBp[i + 3]   = 0.;
            dBnp[i + 3]  = 0.;
            dBenp[i + 3] = 0.;

            // Derivatives wrt to chi1x, chi1y, chi1z
            // FIXME: this can be optimized by defining a vector for X1, X2!

            da0n2dchi1[i]   = 2.*a0n*X1*n[i];
            da0n2u2dchi1[i] = u2*da0n2dchi1[i];
            dnumAdchi1[i]   = da0n2u2dchi1[i] - u4*dAnchiQ[2][i];

            ddenAdchi1[i] = (dDelta[2][i]*a0n2 + Delta*2.*a0n*X1*n[i])*u2*uc2 - Delta*a0n2*u2*2.*uc2*uc*drc[2][i]; 
            dA[i + 6]     = dAeq[2][i] * numA/denA + Aeq * (dnumAdchi1[i]*denA - ddenAdchi1[i]*numA)/(denA*denA);

            ddenBpdchi1[i] = da0n2u2dchi1[i] + u3*dBnchiQ[2][i];
            dBp[i + 6]     = - ddenBpdchi1[i]*Bp2;

            dBnpIdchi1[i]     = - BnpI2*da0n2u2dchi1[i];
            dnumBnpIIdchi1[i] = dAeqB[2][i] * (rc2 + u*BchiQ) + AeqB * (2.*rc*drc[2][i] + u*dBchiQ[2][i]);
            ddenBnpIIdchi1[i] = dDorbdrc*drc[2][i]*r2;
            dBnpIIdchi1[i]    = (dnumBnpIIdchi1[i]*denBnpII - ddenBnpIIdchi1[i]*numBnpII)/denBnpII2;
            dBnp[i + 6]       = dBnpIdchi1[i]*BnpII + BnpI*dBnpIIdchi1[i];   
            
            dBenpIdchi1[i]     = BenpI2 * da0n2u2dchi1[i];
            dnumBenpIIdchi1[i] = da0n2dchi1[i];
            ddenBenpIIdchi1[i] = dR4[2][i] + dDelta[2][i]*a0n2 + Delta*da0n2dchi1[i];
            dBenpIIdchi1[i]    = (dnumBenpIIdchi1[i]*denBenpII - ddenBenpIIdchi1[i]*numBenpII)/denBenpII2;
            dBenp[i + 6]       = dBenpIdchi1[i]*BenpII + BenpI*dBenpIIdchi1[i];          

            // Derivatives wrt to chi2x, chi2y, chi2z

            da0n2dchi2[i]   = 2.*a0n*X2*n[i];
            da0n2u2dchi2[i] = u2*da0n2dchi2[i];
            dnumAdchi2[i]   = da0n2u2dchi2[i] - u4*dAnchiQ[3][i];

            ddenAdchi2[i] = (dDelta[3][i]*a0n2 + Delta*2.*a0n*X2*n[i])*u2*uc2 - Delta*a0n2*u2*2.*uc2*uc*drc[3][i];
            dA[i + 9]     = dAeq[3][i] * numA/denA + Aeq * (dnumAdchi2[i]*denA - ddenAdchi2[i]*numA)/(denA*denA);

            ddenBpdchi2[i] = da0n2u2dchi2[i] + u3*dBnchiQ[3][i];
            dBp[i + 9]     = - ddenBpdchi2[i]*Bp2;

            dBnpIdchi2[i]     = - BnpI2*da0n2u2dchi2[i];
            dnumBnpIIdchi2[i] = dAeqB[3][i] * (rc2 + u*BchiQ) + AeqB * (2.*rc*drc[3][i] + u*dBchiQ[3][i]);
            ddenBnpIIdchi2[i] = dDorbdrc*drc[3][i]*r2;
            dBnpIIdchi2[i]    = (dnumBnpIIdchi2[i]*denBnpII - ddenBnpIIdchi2[i]*numBnpII)/denBnpII2;
            dBnp[i + 9]       = dBnpIdchi2[i]*BnpII + BnpI*dBnpIIdchi2[i]; 

            dBenpIdchi2[i]     = BenpI2 * da0n2u2dchi2[i];
            dnumBenpIIdchi2[i] = da0n2dchi2[i];
            ddenBenpIIdchi2[i] = dR4[3][i] + dDelta[3][i]*a0n2 + Delta*da0n2dchi2[i];
            dBenpIIdchi2[i]    = (dnumBenpIIdchi2[i]*denBenpII - ddenBenpIIdchi2[i]*numBenpII)/denBenpII2;
            dBenp[i + 9]       = dBenpIdchi2[i]*BenpII + BenpI*dBenpIIdchi2[i];  

        }
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                            NEXT-TO-LEADING ORDER SPIN-SPIN CONTRIBUTIONS                           //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// [in]  r              Radial vector
// [in]  nu             Symmetric mass ratio
// [in]  chi1           Spin vector of the first particle
// [in]  chi2           Spin vector of the second particle
// [out] AchiQ          NLO-spin-spin contributions, denoted as AchiQ
// [out] AnchiQ         NLO-spin-spin contributions, denoted as AnchiQ
// [out] BchiQ          NLO-spin-spin contributions, denoted as BchiQ
// [out] BnchiQ         NLO-spin-spin contributions, denoted as BnchiQ
// [out] dAchiQ         Derivative of spin-spin contribution AchiQ wrt all spacetime variables
// [out] dAnchiQ        Derivative of spin-spin contribution AnchiQ wrt all spacetime variables
// [out] dBchiQ         Derivative of spin-spin contribution BchiQ wrt all spacetime variables
// [out] dBnchiQ        Derivative of spin-spin contribution BnchiQ wrt all spacetime variables
void NLO_SpinSpinContributions(double r[], double nu, double chi1[], double chi2[],
                               double *AchiQ, double *AnchiQ, double *BchiQ, double *BnchiQ,
                               double *dAchiQ, double *dAnchiQ, double *dBchiQ, double *dBnchiQ)
{
    // Mass fractions
    double X1, X2;
        get_X1X2(nu, &X1, &X2);

    // Unit radial vector + derivatives
    double n[3], dndr[3][3];
        get_n(r, n, (double *)dndr);
    
    // Auxiliary expressions
    double chi12    = dot(chi1, chi1, 3);  // Dot product of chi1 with itself
    double chi22    = dot(chi2, chi2, 3);  // Dot product of chi2 with itself
    double chi1chi2 = dot(chi1, chi2, 3);  // Dot product of chi1 and chi2
    double nchi1    = dot(n, chi1, 3);     // Dot product of n and chi1
    double nchi12   = nchi1*nchi1;         // Squared nchi1
    double nchi2    = dot(n, chi2, 3);     // Dot product of n and chi2
    double nchi22   = nchi2*nchi2;         // Squared nchi2

    // Coefficients
    double a[]  = {6.*X1 - nu, 6.*X2 - nu, 2. - nu};                // Coefficients for AchiQ
    double an[] = {4.*X1 + 5.*nu, 4.*X2 + 5.*nu, 3. - 7.*nu};       // Coefficients for AnchiQ
    double b[]  = {18.*X1 - 7.5*nu, 18.*X2 - 7.5*nu, 6. + 4.5*nu};  // Coefficients for BchiQ and BnchiQ

    // *******************************
    // * NLO spin-spin contributions *
    // *******************************

    if (AchiQ != NULL) *AchiQ = nu * (0.5*a[0] * chi12 + 0.5*a[1] * chi22 + a[2] * chi1chi2);

    if (AnchiQ != NULL) *AnchiQ = nu * (0.5*an[0] * nchi12 + 0.5*an[1] * nchi22 + an[2] * nchi1 * nchi2);

    if (BchiQ != NULL) *BchiQ  = nu * (0.5*b[0] * chi12 + 0.5*b[1] * chi22 + b[2] * chi1chi2);
    
    if (BnchiQ != NULL) *BnchiQ = nu * (0.5*b[0] * nchi12 + 0.5*b[1] * nchi22 + b[2] * nchi1 * nchi2);
    
    // **********************************************
    // * Derivatives of NLO spin-spin contributions *
    // **********************************************

    if (dAchiQ != NULL) {
        for(int i = 0; i < 3; i++) {
            dAchiQ[i]   = 0.;
            dAchiQ[i+3] = 0.;
            dAchiQ[i+6] = nu * (a[0] * chi1[i] + a[2] * chi2[i]);
            dAchiQ[i+9] = nu * (a[1] * chi2[i] + a[2] * chi1[i]);
        }
    }

    if (dAnchiQ != NULL) {
        for(int i = 0; i < 3; i++) {
            double dnchi1 = dot(dndr[i], chi1, 3);
            double dnchi2 = dot(dndr[i], chi2, 3);

            dAnchiQ[i]   = nu * (an[0] * nchi1 * dnchi1 + an[1] * nchi2 * dnchi2 + an[2] * (dnchi1 * nchi2 + dnchi2 * nchi1));
            dAnchiQ[i+3] = 0.;
            dAnchiQ[i+6] = nu * n[i] * (an[0] * nchi1 + an[2] * nchi2);
            dAnchiQ[i+9] = nu * n[i] * (an[1] * nchi2 + an[2] * nchi1);
        }
    }

    if (dBchiQ != NULL) {
        for(int i = 0; i < 3; i++) {
            dBchiQ[i]   = 0.;
            dBchiQ[i+3] = 0.;
            dBchiQ[i+6] = nu * (b[0] * chi1[i] + b[2] * chi2[i]);
            dBchiQ[i+9] = nu * (b[1] * chi2[i] + b[2] * chi1[i]);
        }
    }

    if (dBnchiQ != NULL) {
        for(int i = 0; i < 3; i++) {
            double dnchi1 = dot(dndr[i], chi1, 3);
            double dnchi2 = dot(dndr[i], chi2, 3);

            dBnchiQ[i]   = nu * (b[0] * nchi1 * dnchi1 + b[1] * nchi2 * dnchi2 + b[2] * (dnchi1 * nchi2 + dnchi2 * nchi1));
            dBnchiQ[i+3] = 0.;
            dBnchiQ[i+6] = nu * n[i] * (b[0] * nchi1 + b[2] * nchi2);
            dBnchiQ[i+9] = nu * n[i] * (b[1] * nchi2 + b[2] * nchi1);
        }
    }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                                       OTHER METRIC FUNCTIONS                                       //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

void get_Aeq(double r[], double nu, double chi1[], double chi2[], double *Aeq, double *dAeq)
{
    // As defined in Eq. 2.78
    
    double modr  = get_mod(r);
    double r2    = modr*modr;
    
    double X1, X2;
        get_X1X2(nu, &X1, &X2);

    double rc, drc[4][3];
        get_rc(r, nu, chi1, chi2, &rc, (double *)drc);

    double rc2 = rc*rc;
    double rc3 = rc2*rc;
    double rc4 = rc2*rc2;

    double n[3];
        get_n(r, n, NULL);

    double Aorb, dAorbdrc;
        get_Aorb(rc, nu, &Aorb, &dAorbdrc); // call it as a function of rc

    double AchiQ, dAchiQ[4][3];
        NLO_SpinSpinContributions(r, nu, chi1, chi2, &AchiQ, NULL, NULL, NULL, (double *)dAchiQ, NULL, NULL, NULL);

    double num  = (1. + 2. / rc + AchiQ / rc4);
    double den  = (1. + 2. / modr);
    double frac = num / den;
    *Aeq = Aorb * frac;

    double dnumdr[3], ddendr[3], dfracdr[3], dnumdchi1[3], dnumdchi2[3];

    for (int i = 0; i < 3; i++)
    {
        // Auxiliary expressions
        dnumdr[i]    = -2. * drc[0][i] * (1. + 2. * AchiQ / rc3) / rc2;
        ddendr[i]    = -2. * n[i] / r2;
        dfracdr[i]   = (dnumdr[i]*den - ddendr[i]*num)/(den*den);
        dnumdchi1[i] = -2. * drc[2][i] * (1. + 2. * AchiQ / rc3) / rc2 + dAchiQ[2][i] / rc4;
        dnumdchi2[i] = -2. * drc[3][i] * (1. + 2. * AchiQ / rc3) / rc2 + dAchiQ[3][i] / rc4;

        // Derivatives of Aeq
        dAeq[i]   = drc[0][i] * dAorbdrc * frac + Aorb * dfracdr[i]; 
        dAeq[i+3] = 0.;
        dAeq[i+6] = drc[2][i] * dAorbdrc * frac + Aorb * dnumdchi1[i] / den;
        dAeq[i+9] = drc[3][i] * dAorbdrc * frac + Aorb * dnumdchi2[i] / den;
    }

}

void get_AeqB(double r[], double nu, double chi1[], double chi2[], double *AeqB, double *dAeqB)
{
    // As defined in Eq. 2.81

    double modr = get_mod(r);
    double r2   = modr*modr;
    
    double rc, drc[4][3];
        get_rc(r, nu, chi1, chi2, &rc, (double *)drc);

    double rc2 = rc*rc;
    double rc3 = rc2*rc;
    double rc4 = rc2*rc2;

    double n[3];
        get_n(r, n, NULL);

    double Aorb, dAorbdrc;
        get_Aorb(rc, nu, &Aorb, &dAorbdrc); // call it as a function of rc

    double num  = (1. + 2. / rc);
    double den  = (1. + 2. / modr);
    double frac = num / den;

    *AeqB = Aorb * frac;

    double dnumdr[3], ddendr[3], dfracdr[3], dnumdchi1[3], dnumdchi2[3];

    for (int i = 0; i < 3; i++)
    {
        // Auxiliary expressions
        dnumdr[i]    = -2. * drc[0][i] / rc2;
        ddendr[i]    = -2. * n[i] / r2;
        dfracdr[i]   = (dnumdr[i]*den - ddendr[i]*num)/(den*den);
        dnumdchi1[i] = -2. * drc[2][i] / rc2;
        dnumdchi2[i] = -2. * drc[3][i] / rc2;

        // Derivatives of AeqB
        dAeqB[i]   = drc[0][i] * dAorbdrc * frac + Aorb * dfracdr[i];
        dAeqB[i+3] = 0.;
        dAeqB[i+6] = drc[2][i] * dAorbdrc * frac + Aorb * dnumdchi1[i] / den;
        dAeqB[i+9] = drc[3][i] * dAorbdrc * frac + Aorb * dnumdchi2[i] / den;
    }

}

void get_Aorb(double modr, double nu, double *Aorb, double *dAorb)
{
    
    // This is the 5PN Padé-resummed (1,5) A potential as defined in Damour-Nagar 2014 + the latest determination of a6c.
    // Defined as a function of |r|, but usually called as a function of the centrifugal radius rc.
    
    // Notes: (i) Could also think of implementing the 1SF version (~complete linear in nu info).
    //        (ii) Could also define a Padé function in utils.c that takes as input the Taylor coefs.
    //             Done already in TEOBREsumS + derivatives already implemented. But we are fine like this for now.

    // shorthands
    double pi2  = Pi*Pi;
    double pi4  = pi2*pi2;
    double u    = 1./modr;
    double u2   = u*u;
    double u3   = u2*u;
    double u4   = u3*u;
    double u5   = u4*u;
    double logu = log(u);
    double nu2  = nu*nu;
    double nu3  = nu2*nu;

    // PN coefficients 
    double a4     = 94./3. - (41./32.)*pi2;
    double a5c0   = -4237./60. + (128./5.)*EulerGamma + (2275./512.)*pi2 + (256./5.)*log(2.);
    double a5c1   = -221./6. + (41./32.)*pi2;
    double a5c    = a5c0 + a5c1*nu;
    double a5log  = 64./5.;
    double a5     = a5c + a5log*logu;
    double da5du  = a5log/u;
    double a52    = a5*a5;
    double a6log  = (-7004./105.) + (-144./5.)*nu;

    // latest tuning of a6c
    // Note: Could also take analytical value (to do: find ref).
    double n0    =  5.9951;
    double n1    = -34.4844;
    double n2    = -79.2997;
    double n3    =  713.4451;
    double d1    = -3.167;
    double a6c   = n0*(1. + n1*nu + n2*nu2 + n3*nu3)/(1. + d1*nu);
    double a6    = a6c + a6log*logu;
    double da6du = a6log/u;

    // coefficients of the Padé
    double den  = -32. + (24. + 4.*a4 + a5)*nu;
    double iden = 1./den;
    double N1   = (64. - (64. + 12.*a4 + 4.*a5 + a6)*nu + 4.*nu2)*iden;
    double D1   = ((-2.*(8. + 2.*a4 + a5) - a6)*nu + 4*nu2)*iden;
    double D2   = 2.*D1;
    double D3   = (-4.*(4.*a4 + 2.*a5 + a6)*nu - 2.*(16. + 4.*a4 + a5)*nu2)*iden;
    double D4   = (-8.*(2.*a5 + a6)*nu + (-(a4*(32. + 4.*a4 + a5)) + 2.*(-16. + a6))*nu2 - 8.*nu3)*iden;
    double D5   = (-16.*a6*nu + (-4.*a4*a4 - a5*(16. + a5) + 8.*a6 + a4*(-4.*(8. + a5) + a6))*nu2 - 4.*(8. + a4)*nu3)*iden;

    double Num = 1. + N1*u;
    double Den = 1. + D1*u + D2*u2 + D3*u3 + D4*u4 + D5*u5;
    *Aorb          = Num/Den;

    // First derivatives 
    double dden = da5du*nu;
    double fact = 2.*da5du + da6du;
    double dN1  = (-((4.*da5du + da6du)*nu) - dden*N1)*iden;
    double dD1  = (-fact*nu - dden*D1)*iden;
    double dD2  = (-2.*fact*nu - dden*D2)*iden;
    double dD3  = (-4.*fact*nu - 2*da5du*nu2 - dden*D3)*iden;
    double dD4  = (-8.*fact*nu + (-(a4*da5du) + 2.*da6du)*nu2 - dden*D4)*iden;
    double dD5  = (-16.*da6du*nu + nu2*(8.*da6du + a4*(-4.*da5du + da6du) - da5du*a5 - da5du*(16. + a5)) - dden*D5)*iden;
            
    double dNum  = dN1*u + N1;
    double dDen  = D1 + u*(dD1 + 2.*D2) + u2*(dD2 + 3.*D3) + u3*(dD3 + 4.*D4) + u4*(dD4 + 5.*D5) + dD5*u5;
    
    // Derivative wrt to |r|
    double dAdu = (dNum*Den - dDen*Num)/(Den*Den);
    *dAorb      = -u2*dAdu;

}

void get_Dorb(double modr, double nu, double *Dorb, double *dDorb)
{
    // D potential at 3PN, Padé-resummed (0,3).
    // Defined as a function of |r|, but usually called as a function of the centrifugal radius rc.
    
    // shorthands
    double u  = 1./modr;  
    double u2 = u*u;                            
    double u3 = u2*u;

    double den  = 1. + 2.*(26. - 3.*nu)*nu*u3 + 6.*nu*u2;
    *Dorb       = 1./den; 

    // Derivative wrt to |r|
    double ddendu  = 6.*nu*u*(2. - (3.*nu-26.)*u);
    double dDorbdu = -ddendu/(den*den);
    *dDorb         = -u2*dDorbdu;

}

void get_Q4(double r[], double p[], double nu, double chi1[], double chi2[], double *Q4, double *dQ4)
{

    double r2      = dot(r, r, 3); // square modulus
    double modr    = sqrt(r2); // modulus
    double u2      = 1./r2; // square modulus of the inverse radius
    double u3      = u2/modr;

    double n[3], dndr[3][3]; 
        get_n(r, n, (double *)dndr);

    double pr, dpr[12] = {0.}; // Q: could just use 6 components? other ones are 0 anyway... will think about it 
        get_pr(r, p, &pr, dpr);
    double pr3 = pr*pr*pr;
    double pr4 = pr3*pr;

    // version of Q4 from Damour-Nagar 2014 (for spin-aligned binaries, using prstar)
    /*
    double A, dAdr;
        get_Aorb(modr, nu, &A, &dAdr);

    double D, dDdr;
        get_Dorb(modr, nu, &D, &dDdr);

    double prstar2 = (A*A*D)*pr2;
    double prstar4 = prstar2*prstar2;

    *Q4 = 2.*nu*(4. - 3.*nu)*prstar4*u2; 
    */

    // Damour 2001: first spinning binaries paper, generic orientation of the spins, defined with dot(p, n)
    // This definition of Q4 allows to have a two-body analogous of the Carter constant (could be checked numerically, interesting)
    // See eq. 2.37 in the paper and discussion after eq. 3.20.

    double costheta  = r[2]/modr; // z/|r|
    double costheta2 = costheta*costheta; 
    double a0[3];
        get_a0(nu, chi1, chi2, a0);
    double a02 = dot(a0, a0, 3);

    double Cnu  = 2.*nu*(4. - 3.*nu); // constant only depending on nu (only holds at 3PN)
    double den  = r2 + a02*costheta2;
    double iden = 1./den;
    *Q4         = Cnu*pr4*iden; 

    // using one single 12-component array for the derivatives wrt to:
    // x, y, z, px, py, pz, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z

    // useful pieces
    double dcostheta2dr[3]; 
    double z2       = r[2]*r[2];
    dcostheta2dr[0] = -2.*z2*n[0]*u3; // d(costheta^2)/dx
    dcostheta2dr[1] = -2.*z2*n[1]*u3; // d(costheta^2)/dy
    dcostheta2dr[2] = 2.*costheta*(1./modr - r[2]*n[2]*u2); // d(costheta^2)/dz
    double iden2    = iden*iden;
    double X1, X2;
        get_X1X2(nu, &X1, &X2);

    for (int i = 0; i < 3; i++)
    {
        dQ4[i]     = Cnu*(4.*pr3*dpr[i] - pr4*(2.*modr*n[i] + a02*dcostheta2dr[i])*iden)*iden;  // wrt to x, y, z
        dQ4[i + 3] = Cnu*4.*pr3*dpr[i + 3]*iden; // wrt to px, py, pz
        dQ4[i + 6] = - Cnu*pr4*iden2 * 2.*a0[i]*X1*costheta2;  // wrt to chi1x, chi1y, chi1z
        dQ4[i + 9] = - Cnu*pr4*iden2 * 2.*a0[i]*X2*costheta2;  // wrt to chi2x, chi2y, chi2z
    }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                                   GYRO-GRAVITOMAGNETIC FUNCTIONS                                   //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// [in]  r              Radial vector
// [in]  pr             Tortoise radial momentum
// [in]  nu             Symmetric mass ratio
// [in]  chi1           Spin vector of the first particle
// [in]  chi2           Spin vector of the second particle
// [out] GS             Gyro-gravitomagnetic function in the DJS gauge, denoted as GS
// [out] GSs            Gyro-gravitomagnetic function in the DJS gauge, denoted as GSs
// [out] dGS            Derivative of gyro-gravitomagnetic function GS wrt all spacetime variables
// [out] dGSs           Derivative of gyro-gravitomagnetic function GSs wrt all spacetime variables
void get_GSGSstar(double r[], double p[], double nu, double chi1[], double chi2[], 
                  double *GS, double *GSs, double *dGS, double *dGSs)
{
    // Gyro-gravitomagnetic functions in DJS gauge @N3LO, from PRD 102, 124024 (2020) 

    // Note: when calling the function initialise all the derivatives to zero, namely:
    // double dGS = {0.}, dGSs = {0.}; (mostly for dGSs whose derivatives wrt to chi1,2 are actually 0)
    
    // get radial versor, modulus of r and u
    double n[3], dndr[3][3]; 
        get_n(r, n, (double *)dndr);
    double modr = get_mod(r);
    double r2   = modr*modr;
    double u    = 1./modr;

    // get a0 and Delta
    double Delta, dDelta[12] = {0.}; 
        get_Delta(r, nu, chi1, chi2, &Delta, dDelta);
    double a0[3];
        get_a0(nu, chi1, chi2, a0);
    double a02 = dot(a0, a0, 3);
    
    // Projection of the effective Kerr parameter along the radial direction
    double a0n  = dot(n, a0, 3); 
    double a0n2 = a0n*a0n;

    // getting centrifugal radius & derivatives 
    double rc, drc[12] = {0.};
        get_rc(r, nu, chi1, chi2, &rc, drc);
    double uc  = 1./rc;
    double uc2 = uc*uc; 

    // get radial momentum
    double pr, dpr[12] = {0.};
        get_pr(r, p, &pr, dpr);
    
    // shorthands
    double pi2 = Pi*Pi;
    double u2  = u*u;
    double u3  = u2*u;
    double pr2 = pr*pr;
    double pr4 = pr2*pr2;
    double pr6 = pr4*pr2;
    double nu2 = nu*nu;
    double nu3 = nu2*nu;

    // coefficients of gS
    double c10 = (-5.*nu)/16.;
    double c20 = (-51./8. - nu/16.)*nu;
    double c30 = nu*(-80399./2304. + (379.*nu)/64. - (7.*nu2)/256. + (241.*pi2)/384.);
    double c02 = (-27.*nu)/16.;
    double c12 = (-21.*nu)/4. + (23.*nu2)/16.;
    double c22 = (-5283.*nu)/128. + (1557.*nu2)/32. + (69.*nu3)/128.;
    double c04 = (5.*nu)/16. + (35.*nu2)/16.;
    double c14 = (781.*nu)/256. + (831.*nu2)/64. - (771.*nu3)/256.;
    double c06 = (7.*nu)/256. - (63.*nu2)/64. - (665.*nu3)/256.;

    // coefficients of gSs
    double cs10 = -3./4. - nu/2.;
    double cs20 = -9./8. - (13.*nu)/2. - nu2/8.;
    double cs30 = -135./64. - (7627.*nu)/288. + (237.*nu2)/32. - nu3/16. + (41.*nu*pi2)/48.;
    double cs02 = -5./4. - (3.*nu)/2.;
    double cs12 = 23./8. - (3.*nu)/2. + (19.*nu2)/8.;
    double cs22 = -15./32. - (279.*nu)/16. + (787.*nu2)/16. + (9.*nu3)/8.;
    double cs04 = 35./24. + (5.*nu)/3. + (15.*nu2)/8.;
    double cs14 = -1105./192. - (53.*nu)/96. + (117.*nu2)/32. - (81.*nu3)/16.;
    double cs06 = -105./64. - (175.*nu)/96. - (77.*nu2)/32. - (35.*nu3)/16.;

    // gS & gSs
    double gS  = 1.  + u*(c10  + u*(c20 + u*c30))   + pr2*(c02  + u*(c12 + u*c22))   + pr4*(c04  + u*c14)  + pr6*c06;
    double gSs = 1.  + u*(cs10 + u*(cs20 + u*cs30)) + pr2*(cs02 + u*(cs12 + u*cs22)) + pr4*(cs04 + u*cs14) + pr6*cs06;

    // derivatives wrt to u and pr
    double dgSdu   = c10  + 2.*u*c20 + 3.*u2*c30  + pr2*(c12 + 2.*u*c22)   + pr4*c14; 
    double dgSsdu  = cs10 + 2.*u*cs20 + 3.*u2*cs30 + pr2*(cs12 + 2.*u*cs22) + pr4*cs14; 
    double dgSdpr  = 2.*pr*(c02  + u*(c12 + u*c22)   + 2.*pr2*(c04  + u*c14)  + 3.*pr4*c06);
    double dgSsdpr = 2.*pr*(cs02 + u*(cs12 + u*cs22) + 2.*pr2*(cs04 + u*cs14) + 3.*pr4*cs06);

    // derivatives wrt to x, y, z, and px, py, pz
    // using one single 12-component array for the derivatives wrt to:
    // x, y, z, px, py, pz, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z

    double dgS[12] = {0.}, dgSs[12] = {0.};

    for (int i = 0; i < 3; i++) 
    {
        // wrt to x, y, z
        dgS[i]      = (-n[i]*u2)*dgSdu + dpr[i]*dgSdpr; 
        dgSs[i]     = (-n[i]*u2)*dgSsdu + dpr[i]*dgSsdpr; 
        // wrt to px, py, pz
        dgS[i + 3]  = dgSdpr*n[i];
        dgSs[i + 3] = dgSsdpr*n[i];  
    }

    // prefactors 
    double den0  = 1. + Delta*a0n2*u2*uc2;
    double iden0 = 1./den0;
    double GS0   = 2.*u*uc2*iden0;
    double GSs0  = 3.*u3/2.; 

    // derivatives wrt to u
    double drcdmodr    = (modr - a02*u2)*uc;
    double ducdu       = (uc2*r2)*drcdmodr;
    // double dDeltadmodr = 2.*(modr - 1.);
    // double dDeltadu    = -r2*dDeltadmodr;
    // double a0r         = dot(a0, r, 3); 

    double dGSs0du = 9.*u2/2.;

    // derivatives wrt to x, y, z
    double da0ndr[3], dden0dr[3], dnumGS0dr[3], dGS0dr[3], dGSs0dr[3];
    for (int i = 0; i < 3; i++)
    {
        // GS0
        da0ndr[i]    = dot(a0, dndr[i], 3); 
        dden0dr[i]   = u2*uc*( uc*(dDelta[i]*a0n2 + Delta*2.*a0n*da0ndr[i]) + 2.*Delta*a0n2*(-n[i]*u)*(uc + u*ducdu) );
        dnumGS0dr[i] = 2.*(-n[i]*u2)*uc*(uc + 2.*u*ducdu);
        dGS0dr[i]    = (dnumGS0dr[i]*den0 - dden0dr[i]*2.*u*uc2)*(iden0*iden0); 
        // GSs0
        dGSs0dr[i]   = -u2*dGSs0du*n[i];
    }

    // complete functions
    *GS  = GS0*gS;   // this multiplies dot(pphi, S)
    *GSs = GSs0*gSs; // this multiplies dot(pphi, Ss)

    double dden0dchi1[3], dden0dchi2[3], dGS0dchi1[3], dGS0dchi2[3];
    double X1, X2;
        get_X1X2(nu, &X1, &X2);

    // FIXME: This can be optimized with a single X vector for X1, X2

    for (int i = 0; i < 3; i++){
        dden0dchi1[i] = (dDelta[i + 6]*a0n2 + Delta*2.*a0n*X1*n[i])*u2*uc2 - Delta*a0n2*u2*2.*uc2*uc*drc[i + 6];
        dGS0dchi1[i]  = -2.*u*uc2*iden0*(2.*uc*drc[i + 6] + dden0dchi1[i]*iden0);
        dden0dchi2[i] = (dDelta[i + 9]*a0n2 + Delta*2.*a0n*X2*n[i])*u2*uc2 - Delta*a0n2*u2*2.*uc2*uc*drc[i + 9];
        dGS0dchi2[i]  = -2.*u*uc2*iden0*(2.*uc*drc[i + 9] + dden0dchi2[i]*iden0);
    }

    // Derivatives
    for (int i = 0; i < 3; i++)
    {
        // wrt to x, y, z
        dGS[i]  = dGS0dr[i]*gS + GS0*dgS[i];
        dGSs[i] = dGSs0dr[i]*gSs + GSs0*dgSs[i]; 
         // wrt to px, py, pz 
        dGS[i + 3]  = GS0*dgS[i + 3];
        dGSs[i + 3] = GSs0*dgSs[i + 3]; 
        // wrt to chi1x, chi1y, chi1z
        dGS[i + 6]  = dGS0dchi1[i]*gS;
        dGSs[i + 6] = 0.;  
        // wrt to chi2x, chi2y, chi2z
        dGS[i + 9]  = dGS0dchi2[i]*gS;
        dGSs[i + 9] = 0.;  
    }

}