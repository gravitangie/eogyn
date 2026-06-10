/*
File: utils.c

This file contains functions for derived variables and other utility operations used throughout the code.
*/

#include "header.h"

// Dot product
double dot(double a[], double b[], int n)
{
    double ans = 0.;

    for (int i = 0; i < n; i++)
        ans += a[i]*b[i];

    return ans;
}

// Cross product
void cross(double a[], double b[], double *c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

// Modulus of a three-component vector
double get_mod(double x[])
{
    return sqrt(dot(x, x, 3));
}

// Symmetric mass ratio
double get_nu(double q)
{
    return q/((1. + q)*(1. + q));
}

// Unit radial vector n = r/|r|
void get_n(double r[], double *n, double *dn)
{
    // Modulus of r
    double modr = get_mod(r);

    // Unit radial vector
    for (int i = 0; i < 3; i++)
        n[i] = r[i]/modr;

    // Derivatives of unit radial vector
    if (dn != NULL)
    {
        double u3 = 1./(modr*modr*modr);
        double x2 = r[0]*r[0];
        double y2 = r[1]*r[1];
        double z2 = r[2]*r[2];

        dn[0] = (y2 + z2)*u3;      // dnx/dx
        dn[1] = (- r[0]*r[1])*u3;  // dny/dx
        dn[2] = (- r[0]*r[2])*u3;  // dnz/dx
        dn[3] = dn[1];             // dnx/dy
        dn[4] = (x2 + z2)*u3;      // dny/dy
        dn[5] = (- r[1]*r[2])*u3;  // dnz/dy
        dn[6] = dn[2];             // dnx/dz
        dn[7] = dn[5];             // dny/dz
        dn[8] = (x2 + y2)*u3;      // dnz/dz
    }
}

// Tortoise radial momentum
void get_pr(double r[], double p[], double *pr, double *dpr)
{
    // Unit radial vector and derivatives
    double n[3], dndr[3][3];
        get_n(r, n, (double *)dndr);

    // Radial momentum
    *pr = dot(p, n, 3);

    // Derivatives of radial momentum
    if (dpr != NULL)
    {
        for (int i = 0; i < 3; i++)
        {
            dpr[i]     = dot(p, dndr[i], 3);  // Derivatives wrt to x, y, z
            dpr[i + 3] = n[i];                // Derivatives wrt to px, py, pz
        }
    }
}

// Mass fractions X1, X2
void get_X1X2(double nu, double *X1, double *X2)
{
    *X1 = 0.5*(1. + sqrt(1. - 4.*nu));
    *X2 = 1. - *X1;
}

// Centrifugal radius
void get_rc(double r[], double nu, double chi1[], double chi2[], double *rc, double *drc)
{
    // Modulus of r
    double modr = get_mod(r);
    
    // Effective spin vector
    double a0[3];
        get_a0(nu, chi1, chi2, a0);

    // Unit radial vector
    double n[3];
        get_n(r, n, NULL);

    // Mass fractions
    double X1, X2;
        get_X1X2(nu, &X1, &X2);

    // Auxiliary expressions
    double a02   = dot(a0, a0, 3);
    double modr2 = modr*modr;
    double u     = 1./modr;

    // Centrifugal radius
    *rc = sqrt(modr2 + a02*(1. + 2.*u));

    // Derivatives of centrifugal radius
    if (drc != NULL)
    {
        // Auxiliary expressions
        double u2       = u*u;
        double drcdmodr = (modr - a02*u2)/(*rc);

        for (int i = 0; i < 3; i++)
        {
            drc[i]   = drcdmodr*n[i];               // Derivatives wrt to x, y, z
            drc[i+3] = 0.;                          // Derivatives wrt to px, py, pz
            drc[i+6] = (1. + 2.*u)*X1*a0[i]/(*rc);  // Derivatives wrt to chi1
            drc[i+9] = (1. + 2.*u)*X2*a0[i]/(*rc);  // Derivatives wrt to chi2
        }            
    }
}

// Delta (EOB analogous of the same function in the Kerr metric)
void get_Delta(double r[], double nu, double chi1[], double chi2[], double *Delta, double *dDelta)
{
    // Effective spin vector
    double a0[3];
        get_a0(nu, chi1, chi2, a0);

    // Unit radial vector
    double n[3];
        get_n(r, n, NULL);

    // Auxiliary expressions
    double a02  = dot(a0, a0, 3);
    double modr = get_mod(r);

    // Delta
    *Delta = modr*modr - 2.*modr + a02;

    // Derivatives of Delta
    if (dDelta != NULL)
    {
        double dDeltadmodr = 2.*(modr - 1.);
        double X1, X2;
            get_X1X2(nu, &X1, &X2);

        for (int i = 0; i < 3; i++)
        {
            dDelta[i]   = dDeltadmodr*n[i];  // Derivatives wrt to x, y, z
            dDelta[i+6] = 2.*a0[i]*X1;       // Derivatives wrt to chi1
            dDelta[i+9] = 2.*a0[i]*X2;       // Derivatives wrt to chi2
        }
    }
}

// Effective Kerr parameter a0 = chi0
void get_a0(double nu, double chi1[], double chi2[], double *a0)
{
    // Mass fractions
    double X1, X2;
        get_X1X2(nu, &X1, &X2);

    for (int i = 0; i < 3; i++)
        a0[i] = chi1[i]*X1 + chi2[i]*X2;
}

// Effective spins
void get_SSstar(double nu, double chi1[], double chi2[], double *S, double *Ss)
{
    // Mass fractions
    double X1, X2;
        get_X1X2(nu, &X1, &X2);

    for (int i = 0; i < 3; i++)
    {
        S[i]  = chi1[i]*X1*X1 + chi2[i]*X2*X2;
        Ss[i] = nu*(chi1[i] + chi2[i]);
    }
}

// R4 function
void get_R4(double r[], double nu, double chi1[], double chi2[], double *R4, double *dR4)
{
    // Modulus of r
    double modr = get_mod(r);

    // Unit radial vector
    double n[3];
        get_n(r, n, NULL);

    // Effective spin vector
    double a0[3];
        get_a0(nu, chi1, chi2, a0);

    // Auxiliary expressions
    double a02 = dot(a0, a0, 3);
    double r3  = modr*modr*modr;
    double r4  = r3*modr;

    // R4 function
    *R4  = r4 + a02*modr*(modr + 2.);

    // Derivatives of R4 function
    if (dR4 != NULL)
    {
        // Auxiliary expression
        double dR4dmodr = 4.*r3 + 2.*a02*(modr + 1.);

        // Mass fractions
        double X1, X2;
            get_X1X2(nu, &X1, &X2);

        for (int i = 0; i < 3; i++)
        {
            dR4[i]   = dR4dmodr*n[i];                 // Derivatives wrt to x, y, z
            dR4[i+6] = 2.*a0[i]*modr*(modr + 2.)*X1;  // Derivatives wrt to chi1
            dR4[i+9] = 2.*a0[i]*modr*(modr + 2.)*X2;  // Derivatives wrt to chi2
        }            
    }
}

// Orbital angular momentum
void get_l(double r[], double p[], double *l, double *dl)
{

    cross(r, p, l);

    // Derivatives wrt to x, y, z
    dl[0] = 0.;    // dlx/dx
    dl[1] = -p[2]; // dly/dx
    dl[2] = p[1];  // dlz/dx
    dl[3] = p[2];  // dlx/dy
    dl[4] = 0.;    // dly/dy
    dl[5] = -p[0]; // dlz/dy
    dl[6] = -p[1]; // dlx/dz 
    dl[7] = p[0];  // dly/dz
    dl[8] = 0.;    // dlz/dz

    // Derivatives wrt to px, py, pz
    dl[9]  = 0.;    // dlx/dpx
    dl[10] = r[2];  // dly/dpx
    dl[11] = -r[1]; // dlz/dpx
    dl[12] = -r[2]; // dlx/dpy
    dl[13] = 0.;    // dly/dpy
    dl[14] = r[0];  // dlz/dpy
    dl[15] = r[1];  // dlx/dpz
    dl[16] = -r[0]; // dly/dpz
    dl[17] = 0.;    // dlz/dpz

}

// Dot product of l with the spin variables S, Ss and its derivatives
void get_lSSs(double r[], double p[], double nu, double chi1[], double chi2[], double *lS, double *dlS, double *lSs, double *dlSs)
{
    double l[3], dl[6][3]; 
        get_l(r, p, l, (double *)dl);

    double S[3], Ss[3];
        get_SSstar(nu, chi1, chi2, S, Ss);

    double X1, X2;
        get_X1X2(nu, &X1, &X2);
        double X[] = {X1, X2};

    *lS  = dot(l, S, 3);   
    *lSs = dot(l, Ss, 3);

    for (int i = 0; i < 6; i++)
    {
        dlS[i]  = dot(dl[i], S, 3);
        dlSs[i] = dot(dl[i], Ss, 3);
    }

    for (int j = 0; j < 2; j++)
    {
        for (int i = 0; i < 3; i++)
        {
            dlS[i + 6 + 3*j]  = X[j]*X[j]*l[i];
            dlSs[i + 6 + 3*j] = nu*l[i];
        }
    }

}

// Functions to trim strings
void trimString(char *str)
{
    char *end;
    // Trim leading whitespace
    while (isspace((unsigned char)*str)) {
        str++;
    }
    if (*str == '\0') { // All whitespace
        return;
    }
    // Trim trailing whitespace
    end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) {
        end--;
    }
    // Null-terminate the trimmed string
    *(end + 1) = '\0';
}

void swap(double *xp, double *yp) 
{
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

// An optimized version of Bubble Sort
void bubbleSort(double arr[], int n) 
{
    int i, j;
    bool swapped;
    for (i = 0; i < n - 1; i++) {
        swapped = false;
        for (j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                swap(&arr[j], &arr[j + 1]);
                swapped = true;
            }
        }

        // If no two elements were swapped by inner loop,
        // then break
        if (swapped == false)
            break;
    }
}

void CartesianToSpherical(double r[], double p[], double *Q, double *P) 
{
    // Here Q = {R, theta, phi} and P = {p_r, p_theta, p_phi}

    // Unwrap Cartesian coordinates & momenta
    double x  = r[0];
    double y  = r[1];
    double z  = r[2];
    double px = p[0];
    double py = p[1];
    double pz = p[2];

    double R = get_mod(r);

    double phi = atan(y/x); // check if this is okay
    // double phi = atan(fabs(y)/fabs(x));

    double costheta  = z/R; 
    double theta     = acos(costheta);
    double sintheta  = sqrt(1. - costheta*costheta);

    double pr;
        get_pr(r, p, &pr, NULL);
    
    // p_theta  = R*(- pz * sintheta + costheta * (px * cosphi + py * sinphi));
    double ptheta  = - R * pz * sintheta + px * (x*z/sintheta) + py * (y*z/sintheta); 

    // pphi  = R * sintheta * (py * cosphi - px * sinphi);
    double pphi  = x*py - y*px;

    Q[0] = R;
    Q[1] = theta;
    Q[2] = phi;
    P[0] = pr;
    P[1] = ptheta;
    P[2] = pphi;

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//             Rotation helpers for the canonical spin chart (rotated-chart option)                   //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

// Build a rotation matrix R (lab -> chart frame) such that the chart-frame
// z-axis is orthogonal to the initial spin direction. With this choice the
// initial canonical variables are (alpha = 0, xi = 0), as far as possible
// from the chart pole.
// Convention: rows of R are the chart-frame basis vectors (x', y', z')
// expressed in the lab frame. So:
//      v_chart = R v_lab,    v_lab = R^T v_chart.
// We pick x' = chi/|chi| (initial spin along the chart x-axis), z' as a
// unit vector perpendicular to x' (the chart pole), and y' = z' x x'.
// If |chi| < tol the rotation is set to the identity (no chart needed).
void BuildSpinChartRotation(double chi[], double R[3][3])
{
    double modchi = get_mod(chi);
    if (modchi < 1.0e-15) {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                R[i][j] = (i == j) ? 1.0 : 0.0;
        return;
    }

    // x' = chi / |chi|
    double ex[3];
    for (int i = 0; i < 3; i++) ex[i] = chi[i] / modchi;

    // Pick the lab basis vector farthest from ex (smallest |ex[k]|) to
    // build a numerically well-conditioned cross product.
    int k_min = 0;
    if (fabs(ex[1]) < fabs(ex[k_min])) k_min = 1;
    if (fabs(ex[2]) < fabs(ex[k_min])) k_min = 2;
    double v[3] = {0.0, 0.0, 0.0};
    v[k_min] = 1.0;

    // z' = (ex x v) / |ex x v|  (perpendicular to ex by construction)
    double ez[3];
    cross(ex, v, ez);
    double modez = get_mod(ez);
    for (int i = 0; i < 3; i++) ez[i] /= modez;

    // y' = z' x x' (right-handed, automatically unit)
    double ey[3];
    cross(ez, ex, ey);

    // Pack into R (rows = x', y', z')
    for (int i = 0; i < 3; i++) {
        R[0][i] = ex[i];
        R[1][i] = ey[i];
        R[2][i] = ez[i];
    }
}

// v_out = R v_in (lab -> chart)
void RotateVec(double R[3][3], double v_in[], double v_out[])
{
    for (int i = 0; i < 3; i++) {
        v_out[i] = R[i][0]*v_in[0] + R[i][1]*v_in[1] + R[i][2]*v_in[2];
    }
}

// v_out = R^T v_in (chart -> lab)
void RotateVecT(double R[3][3], double v_in[], double v_out[])
{
    for (int i = 0; i < 3; i++) {
        v_out[i] = R[0][i]*v_in[0] + R[1][i]*v_in[1] + R[2][i]*v_in[2];
    }
}

// Rodrigues' rotation formula: rotate vector v by `angle` (right-handed) about
// the UNIT axis k. This is the exact flow of a constant-axis precession
// dv/dt = Omega x v, with k = Omega/|Omega| and angle = |Omega| Delta t.
// Norm-preserving for any angle (rotations are isometries), which is exactly
// why it tolerates arbitrarily fast precession without a step-size limit.
//   out = v cos(a) + (k x v) sin(a) + k (k.v)(1 - cos(a))
void RodriguesRotate(double v[], double k[], double angle, double out[])
{
    double c = cos(angle);
    double s = sin(angle);
    double kxv[3];
        cross(k, v, kxv);
    double kdotv = dot(k, v, 3);
    for (int i = 0; i < 3; i++)
        out[i] = v[i]*c + kxv[i]*s + k[i]*kdotv*(1.0 - c);
}
