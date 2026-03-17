/*
File: solvers.c

This file implements various ODE solvers and root-finding algorithms.
*/

#include "header.h"

// Fourth order Runge-Kutta solver
void RK_Fourth(double t, double *Y, double h, int N_eq,
               void (*RHS)(double, double*, double*))
{
    int i;
    double Y1[N_eq], Y2[N_eq], Y3[N_eq],
           k1[N_eq], k2[N_eq], k3[N_eq], k4[N_eq];
    RHS(t, Y, k1);
    for(i = 0; i < N_eq; i++)
    {
        Y1[i] = Y[i] + 0.5*h*k1[i];
    }
    RHS(t + 0.5*h, Y1, k2);
    for(i = 0; i < N_eq; i++)
    {
        Y2[i] = Y[i] + 0.5*h*k2[i];
    }
    RHS(t + 0.5*h, Y2, k3);
    for(i = 0; i < N_eq; i++)
    {
        Y3[i] = Y[i] + h*k3[i];
    }
    RHS(t + h, Y3, k4);
    for(i = 0; i < N_eq; i++)
    {
        Y[i] += h*(k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i])/6.;
    }
}

// Sixth-order 3-stage Gauss-Legendre collocation solver
void RK_GaussLegendre6(double t, double *Y, double h, int N_eq,
                       void (*RHS)(double, double*, double*))
{
    const double s15 = sqrt(15.0);

    const double c1  = 0.5 - s15/10.0;
    const double c2  = 0.5;
    const double c3  = 0.5 + s15/10.0;

    const double a11 = 5.0/36.0;
    const double a22 = 2.0/9.0;
    const double a33 = a11;

    const double a12 = a22 - s15/15.0;
    const double a13 = a11 - s15/30.0;
    const double a21 = a11 + s15/24.0;
    const double a23 = a11 - s15/24.0;
    const double a31 = a11 + s15/30.0;
    const double a32 = a22 + s15/15.0;

    const double b1 = 2.*a11;
    const double b2 = 2.*a22;
    const double b3 = 2.*a11;

    // Fixed-point iteration parameters
    const int    max_iter = pars->max_iter_RKGL6; // 50; // 10;
    const double tol      = pars->tol_RKGL6; // 1e-15; // 1e-10;

    // Stage derivatives
    double k1[N_eq],   k2[N_eq],  k3[N_eq];
    double nk1[N_eq], nk2[N_eq], nk3[N_eq];

    // Stage solution values
    double Y1[N_eq], Y2[N_eq], Y3[N_eq];

    // Initial guess: use RHS at the same Y, but shifted in time
    RHS(t + c1*h, Y, k1);
    RHS(t + c2*h, Y, k2);
    RHS(t + c3*h, Y, k3);

    for (int it = 0; it < max_iter; it++)
    {
        // Compute stage states from current stage derivatives
        for (int i = 0; i < N_eq; i++)
        {
            Y1[i] = Y[i] + h*(a11*k1[i] + a12*k2[i] + a13*k3[i]);
            Y2[i] = Y[i] + h*(a21*k1[i] + a22*k2[i] + a23*k3[i]);
            Y3[i] = Y[i] + h*(a31*k1[i] + a32*k2[i] + a33*k3[i]);
        }

        // New stage derivatives
        RHS(t + c1*h, Y1, nk1);
        RHS(t + c2*h, Y2, nk2);
        RHS(t + c3*h, Y3, nk3);

        // Measure convergence: max over all components and stages
        double d1 = 0.0, d2 = 0.0, d3 = 0.0;
        double mK = 0.0;

        for (int i = 0; i < N_eq; i++)
        {
            double ad1 = fabs(nk1[i] - k1[i]);
            double ad2 = fabs(nk2[i] - k2[i]);
            double ad3 = fabs(nk3[i] - k3[i]);

            if (ad1 > d1) d1 = ad1;
            if (ad2 > d2) d2 = ad2;
            if (ad3 > d3) d3 = ad3;

            double ak1 = fabs(k1[i]);
            double ak2 = fabs(k2[i]);
            double ak3 = fabs(k3[i]);

            if (ak1 > mK) mK = ak1;
            if (ak2 > mK) mK = ak2;
            if (ak3 > mK) mK = ak3;
        }

        double scale = 1.0 + mK;
        double maxd  = fmax(fmax(d1, d2), d3);

        // Update k's
        for (int i = 0; i < N_eq; i++)
        {
            k1[i] = nk1[i];
            k2[i] = nk2[i];
            k3[i] = nk3[i];
        }

        if (maxd < tol * scale)
            break;
    }

    // Final update (collocation quadrature)
    for (int i = 0; i < N_eq; ++i)
    {
        Y[i] += h * (b1*k1[i] + b2*k2[i] + b3*k3[i]);
    }
}

// Secant method for root finding 
// Note: this is for one root of one function!
void Secant(double(*fun)(double), double x0, double tol, double *zero)
{
    double a  = x0 - 1.;
    double b  = x0 + 1.;

    // printf("a = %.16f, b = %.16f\n", a, b);

    double fa = fun(a);
    double fb = fun(b);
    double dx;
    int imax = 100;

    // printf("f(a) = %.16f, f(b) = %.16f\n", fa, fb);

    *zero = 0.;

    if(fa*fb > 0.)
    {
        fprintf(stderr, "Failed to find root: the function is not changing sign in the provided interval.\n");
    }
    else
    {
        dx = b - a;
        int i = 0;
        while(fabs(dx) > tol)
        {
            i++;
            dx = fb*(b - a)/(fb - fa);
            a = b;
            fa = fb;
            b -= dx;
            fb = fun(b);
            if(i > imax)
            {
                fprintf(stderr, "Failed to find root: maximum number of iterations reached.\n");
                break;
            }
        }
        *zero = b;
        // if(fabs(zero - 0.) < tol) zero = 0.;
    }
}

// Newton's method for root finding
double Newton(double(*fun)(double), double(*fun_der)(double), double a, double b, double tol)
{
    double fa, fb, dx, m, q, xm, xm_prec, fm, fdm;
    int i = 1, imax = 1000;

    fa = fun(a);   
    fb = fun(b);

    // printf("f(a) = %.16f, f(b) = %.16f\n", fa, fb);
    
    if(fa*fb > 0.) {
        printf("Failed to find root: the function is not changing sign in the provided interval.\n");
        exit(EXIT_FAILURE);
    } else {
        dx = b - a;
        xm = 0.5*(a + b); 
        while(fabs(dx) > tol)
        {
            xm_prec = xm;
            fm = fun(xm_prec);
            fdm = fun_der(xm_prec);
            xm = xm_prec - (fm/fdm);
            dx = xm - xm_prec;
            i++;

            if(i > imax)
            {
                printf("Failed to find root: maximum number of iterations reached.\n");
                exit(EXIT_FAILURE);
            }
        }
        // zero = xm;
        // if(fabs(zero - 0.) < tol) zero = 0.;
    }
    return xm;
}

void QuarticPolyRoots(double coefs[], double roots[]) {

    // Unwrap coefficients of the quartic order polynomial:
    // a x^4 + b x^3 + c x^2 + d x + e 
    double a = coefs[0];
    double b = coefs[1];
    double c = coefs[2];
    double d = coefs[3];
    double e = coefs[4];

    // Shortcuts
    double a2 = a*a;
    double a3 = a2*a;
    double b2 = b*b;
    double b3 = b2*b;
    double c2 = c*c;
    double c3 = c2*c;
    double d2 = d*d;

    double p = (8.*a*c - 3.*b2)/(8.*a2);
    double q = (b3 - 4.*a*b*c + 8.*a2*d)/(8.*a3);

    double Delta0 = c2 - 3.*b*d + 12.*a*e;
    double Delta1 = 2.*c3 - 9.*b*c*d + 27.*b2*e + 27.*a*d2 - 72.*a*c*e;

    double sqrtDelta0 = sqrt(Delta0);
    double phi = acos(Delta1 / (2.*sqrtDelta0*sqrtDelta0*sqrtDelta0));
    double S   = 0.5*sqrt(-2.*p/3. + (2./(3.*a)) * sqrtDelta0 * cos(phi/3.));
    double S2  = S*S;

    roots[0] = - b/(4.*a) - S + 0.5 * sqrt(-4.*S2 - 2.*p + q/S);
    roots[1] = - b/(4.*a) - S - 0.5 * sqrt(-4.*S2 - 2.*p + q/S);
    roots[2] = - b/(4.*a) + S + 0.5 * sqrt(-4.*S2 - 2.*p - q/S);
    roots[3] = - b/(4.*a) + S - 0.5 * sqrt(-4.*S2 - 2.*p - q/S);

}
