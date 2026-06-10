/*
File: main.c

This file contains the main function, which serves as the entry point for the program's execution.
*/

#include "header.h"

int main(int argc, char *argv[]){

    // Setting initial parameters
    // pars is a structure of type Parameters (defined in header.h)
    AllocateParameters(&pars);
    SetDefaults(pars); // setting default values
    ParseCommandLine(argc, argv, pars);
    SetParameters(); 

    printf("q set to %.10e\n", pars->q);
    printf("nu set to %.10e\n", pars->nu);

    // ******************
    // Initial conditions
    // ******************

    double pphi0 = 0., x0 = 0., p0[3] = {0.};
    if (pars->motion == Motion_Circular) {
        // Evaluate initial angular momentum with circular condition
        printf("Finding initial conditions for circular motion.\n");
        CircularICs(&pphi0);
    } else if (pars->motion == Motion_EquatorialEccentric) { 
        // Give error if the energy is larger than one (need to have bound motion)
        if (pars->E0 > 1) {
            perror("Error: You chose eccentric motion - the energy needs to be lower than 1!\n");
            return 1;
        }
        // Evaluate initial starting point with given energy and angular momentum 
        // (here x0 is chosen to be the apoapsis)
        printf("Finding initial conditions for eccentric equatorial motion.\n");
        EquatorialICs(&x0);
    } else if (pars->motion == Motion_Generic) { 
        // Evaluate the initial momentum pz corresponding to the initial energy
        printf("Finding initial conditions for generic motion.\n");
        GenericICs(p0);
    }

    // Folder name 
    char folder[100];

    if (pars->motion == Motion_Circular) {
        if (pars->chi1x0 > 1e-15 || pars->chi1y0 > 1e-15 || pars->chi2x0 > 1e-15 || pars->chi2y0 > 1e-15) {
            printf("******************************************************************************************\n");
            printf("Warning: circular ICs + in-plane components of the spins -> motion will be off equatorial!\n");
            printf("******************************************************************************************\n");
            snprintf(folder, sizeof(folder), "circ_q_%.2f_chi1-0_%.2f_%.2f_%.2f_chi2-0_%.2f_%.2f_%.2f", pars->q, pars->chi1x0, pars->chi1y0, pars->chi1z0, pars->chi2x0, pars->chi2y0, pars->chi2z0);
        } else {
            snprintf(folder, sizeof(folder), "circ_q_%.2f_chi1z0_%.2f_chi2z0_%.2f", pars->q, pars->chi1z0, pars->chi2z0);
        }
    } else if (pars->motion == Motion_EquatorialEccentric) {
        if (pars->chi1x0 > 1e-15 || pars->chi1y0 > 1e-15 || pars->chi2x0 > 1e-15 || pars->chi2y0 > 1e-15) {
            printf("******************************************************************************************************\n");
            printf("Warning: equatorial eccentric ICs + in-plane components of the spins -> motion will be off equatorial!\n");
            printf("******************************************************************************************************\n");
            snprintf(folder, sizeof(folder), "eq_q_%.2f_chi1-0_%.2f_%.2f_%.2f_chi2-0_%.2f_%.2f_%.2f", pars->q, pars->chi1x0, pars->chi1y0, pars->chi1z0, pars->chi2x0, pars->chi2y0, pars->chi2z0);
        } else {
            snprintf(folder, sizeof(folder), "eq_q_%.2f_chi1z0_%.2f_chi2z0_%.2f", pars->q, pars->chi1z0, pars->chi2z0);
        }
    } else if (pars->motion == Motion_Generic) {
        if (pars->chi1x0 > 1e-15 || pars->chi1y0 > 1e-15 || pars->chi2x0 > 1e-15 || pars->chi2y0 > 1e-15) {
            snprintf(folder, sizeof(folder), "gen_q_%.2f_chi1-0_%.2f_%.2f_%.2f_chi2-0_%.2f_%.2f_%.2f", pars->q, pars->chi1x0, pars->chi1y0, pars->chi1z0, pars->chi2x0, pars->chi2y0, pars->chi2z0);
        } else {
            snprintf(folder, sizeof(folder), "gen_q_%.2f_chi1z0_%.2f_chi2z0_%.2f", pars->q, pars->chi1z0, pars->chi2z0);
        }       
    }

    // Create the directory
    if (MAKE_DIR(folder) == -1) {
        if (errno != EEXIST) { // Do not give error if dir already exists
            perror("Error creating directory");
            return 1;
        }
    }

    // File path
    char filepath[150];
    snprintf(filepath, sizeof(filepath), "%s/dyn.txt", folder);
    printf("Filepath will be: %s\n", filepath);

    // Open the file for writing
    FILE* fp;
    fp = fopen(filepath, "w+");
    if (fp == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Write the headers 
    fprintf(fp, "0:t\t 1:x\t 2:y\t 3:z\t 4:px\t 5:py\t 6:pz\t 7:chi1x\t 8:chi1y\t 9:chi1z\t 10:chi2x\t 11:chi2y\t 12:chi2z\t 13:E\t 14:lx\t 15:ly\t 16:lz\n");

    // Wrap the initial values into one single variable
    // FIXME: Use dynamical memory allocation (can have w[12] and W[10] without transformed time step)
    double w[14];
    double W[12];

    if (pars->motion == Motion_Circular) {

        // Here the ICs find pphi0, x0 is given
        w[0]  = pars->x0;
        w[1]  = 0.; // pars->y0;
        w[2]  = 0.; // pars->z0;
        w[3]  = 0.; 
        w[4]  = pphi0/(pars->x0); 
        w[5]  = 0.; 
        
    } else if (pars->motion == Motion_EquatorialEccentric) { 

        // Here the ICs find x0 (fixed at the apoapsis), pphi0 set in the parfile
        w[0]  = x0;
        w[1]  = 0.;
        w[2]  = 0.;
        w[3]  = 0.; 
        w[4]  = (pars->pphi0)/x0; 
        w[5]  = 0.; 

    } else if (pars->motion == Motion_Generic) { 

        // Here x0, y0, z0 are given and the ICs find the momenta given initial
        w[0]  = pars->x0;
        w[1]  = pars->y0;
        w[2]  = pars->z0;
        w[3]  = p0[0]; 
        w[4]  = p0[1]; 
        w[5]  = p0[2]; 
    }

    // Spins
    w[6]  = pars->chi1x0;
    w[7]  = pars->chi1y0;
    w[8]  = pars->chi1z0;
    w[9]  = pars->chi2x0;
    w[10] = pars->chi2y0;
    w[11] = pars->chi2z0;

    // Integration details
    double t    = 0.;
    double s    = 0.;
    double ds   = pars->ds; // step in rescaled time s (parfile key "ds"; default 0.0001)
    double dt   = pars->dt; // time step (for RK4; not used for RKGL6)
    double tmax = pars->tmax;

    // Radii to be used for stopping conditions
    double rLR;
        AdiabaticLightRing(&rLR);
    printf("Adiabatic (nonspinning) light ring: %.10f\n", rLR);
    double rmax = 500;

    // Set up the variables
    double r[3], p[3], chi1[3], chi2[3], H, l[3], dummy, dummyd[12], dummyd2[18]; 
    double modr;

    // Unwrap to call the Hamiltonian
    for (int i = 0; i < 3; i++) {
        r[i]    = w[i];
        p[i]    = w[i + 3];
        chi1[i] = w[i + 6];
        chi2[i] = w[i + 9];
    }

    // Evaluate the Hamiltonian and the angular momentum (to be written in the output file)
    Hamiltonian(r, p, pars->nu, chi1, chi2, &dummy, &H, dummyd, dummyd); 
    get_l(r, p, l, dummyd2); // l = r x p = L/(mu M), dimensionless orbital angular momentum
    double H0 = H; // Save initial value

    // Evaluate the Carter-like constant (Unused now)
    // Q = CarterLikeConstant(r, p, chi1, chi2);
    
    // Get spherical coordinates & momenta (Unused now)
    // double Q[3], P[3];
    //    CartesianToSpherical(r, p, Q, P); 

    fprintf(fp, "%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n",
        t, r[0], r[1], r[2], p[0], p[1], p[2], chi1[0], chi1[1], chi1[2], chi2[0], chi2[1], chi2[2], H*(pars->nu), l[0], l[1], l[2]);  

    // Save the values of the spin magnitudes (they are constant)
    // FIXME: They could be added in the pars structure.
    double modchi1 = get_mod(chi1); 
    double modchi2 = get_mod(chi2); 
    printf("Spin magnitude |chi1| = %.16f\n", modchi1);
    printf("Spin magnitude |chi2| = %.16f\n", modchi2);

    // Choice of coordinates: transform the spins into other coordinates that are canonical. Phase space preserves symplectic structure.
    // If evolving extended phase space, add extra variables.
    const double eps = 1.0e-30;
    double alpha1, xi1, alpha2, xi2, sqrt1mxi12, sqrt1mxi22;

    if (pars->step == Step_Split) {
        // Strang split: spins stay as VECTORS (w[6..11]); no canonical chart,
        // no chart pole, no 1/sqrt(1-xi^2). Just set up the extended phase space.
        w[12] = 0.;   // Initial time
        w[13] = -H0;  // Conjugate momentum to the physical time, pt
        printf("Spin evolution: Strang operator split (RKGL6 orbital + exact Rodrigues spin rotation).\n");
    } else if (pars->coords == Coords_Canonical) {
        // Build per-spin rotated charts so that the chart pole is orthogonal
        // to the initial spin direction (Option A: pole fixed at init).
        // After this, the chart-frame initial spin lies along x', i.e.
        // alpha = 0 and xi = 0 -- as far as possible from the chart pole.
        BuildSpinChartRotation(chi1, pars->R1);
        BuildSpinChartRotation(chi2, pars->R2);

        // Move the initial spins into the chart frame to read off (alpha, xi).
        double chi1_chart[3], chi2_chart[3];
        RotateVec(pars->R1, chi1, chi1_chart);
        RotateVec(pars->R2, chi2, chi2_chart);

        // New canonical variables to be evolved
        for (int i = 0; i < 6; i++) W[i] = w[i]; // x, y, z, px, py, pz
        W[6] = atan2(chi1_chart[1], chi1_chart[0]); // alpha1 in chart frame
        W[7] = chi1_chart[2] / (modchi1 + eps);     // xi1   in chart frame
        W[8] = atan2(chi2_chart[1], chi2_chart[0]); // alpha2 in chart frame
        W[9] = chi2_chart[2] / (modchi2 + eps);     // xi2   in chart frame

        printf("Canonical chart 1: chart pole (lab) = (%.6f, %.6f, %.6f), alpha1_0 = %.6e, xi1_0 = %.6e\n",
               pars->R1[2][0], pars->R1[2][1], pars->R1[2][2], W[6], W[7]);
        printf("Canonical chart 2: chart pole (lab) = (%.6f, %.6f, %.6f), alpha2_0 = %.6e, xi2_0 = %.6e\n",
               pars->R2[2][0], pars->R2[2][1], pars->R2[2][2], W[8], W[9]);

        if (pars->step == Step_Transformed) {
            W[10] = 0.; // Initial time
            W[11] = -H0; // Conjugate momentum to the physical time, pt
        }
    } else {
        if (pars->step == Step_Transformed) {
            w[12] = 0.; // Initial time
            w[13] = -H0; // Conjugate momentum to the physical time, pt
        }
    }

    // ****************
    // Solving the ODEs
    // ****************

    // Poincare section output (only for the Strang split). Records z = 0 upward
    // crossings with the spins written as the canonical pair (alpha_i, S_{i,z}).
    FILE *fps = NULL;
    if (pars->step == Step_Split) {
        char sectionpath[200];
        snprintf(sectionpath, sizeof(sectionpath), "%s/section.txt", folder);
        fps = fopen(sectionpath, "w+");
        if (fps != NULL)
            fprintf(fps, "0:t\t 1:x\t 2:y\t 3:px\t 4:py\t 5:pz\t 6:alpha1\t 7:S1z\t 8:alpha2\t 9:S2z\n");
    }

    // Previous lab-frame state, for section-crossing interpolation
    double tprev = t;
    double rprev[3], pprev[3], c1prev[3], c2prev[3];
    for (int i = 0; i < 3; i++) { rprev[i]=r[i]; pprev[i]=p[i]; c1prev[i]=chi1[i]; c2prev[i]=chi2[i]; }

    do {

        if (pars->step == Step_Split) { // Strang split on vector spins

            StrangStepTransformed(s, w, ds);

            // Update variables (spins are already vectors -- no reconstruction)
            for (int i = 0; i < 3; i++) {
                r[i]    = w[i];
                p[i]    = w[i + 3];
                chi1[i] = w[i + 6];
                chi2[i] = w[i + 9];
            }

            t = w[12];
            s = s + ds;

        } else if (pars->coords == Coords_Canonical) {

            if (pars->step == Step_Transformed) { 

                ODESolver(s, W, ds, 12, get_rhs_canonical_transformed);

                // Update variables
                for (int i = 0; i < 3; i++) {
                    r[i]    = W[i];
                    p[i]    = W[i + 3];
                }
                alpha1 = W[6];
                xi1    = W[7];
                alpha2 = W[8];
                xi2    = W[9];

                sqrt1mxi12 = sqrt(1. - xi1*xi1);
                sqrt1mxi22 = sqrt(1. - xi2*xi2);

                // Reconstruct spins in the rotated chart frame, then rotate to lab.
                {
                    double chi1_chart[3], chi2_chart[3];
                    chi1_chart[0] = modchi1 * sqrt1mxi12 * cos(alpha1);
                    chi1_chart[1] = modchi1 * sqrt1mxi12 * sin(alpha1);
                    chi1_chart[2] = modchi1 * xi1;
                    chi2_chart[0] = modchi2 * sqrt1mxi22 * cos(alpha2);
                    chi2_chart[1] = modchi2 * sqrt1mxi22 * sin(alpha2);
                    chi2_chart[2] = modchi2 * xi2;
                    RotateVecT(pars->R1, chi1_chart, chi1);
                    RotateVecT(pars->R2, chi2_chart, chi2);
                }

                t = W[10];
                s = s + ds;
                // E = W[11];

            } else { // Evolve the standard system

                ODESolver(t, W, dt, 10, get_rhs_canonical);

                // Update variables
                for (int i = 0; i < 3; i++) {
                    r[i]    = W[i];
                    p[i]    = W[i + 3];
                }
                alpha1 = W[6];
                xi1    = W[7];
                alpha2 = W[8];
                xi2    = W[9];

                sqrt1mxi12 = sqrt(1. - xi1*xi1);
                sqrt1mxi22 = sqrt(1. - xi2*xi2);

                // Reconstruct spins in the rotated chart frame, then rotate to lab.
                {
                    double chi1_chart[3], chi2_chart[3];
                    chi1_chart[0] = modchi1 * sqrt1mxi12 * cos(alpha1);
                    chi1_chart[1] = modchi1 * sqrt1mxi12 * sin(alpha1);
                    chi1_chart[2] = modchi1 * xi1;
                    chi2_chart[0] = modchi2 * sqrt1mxi22 * cos(alpha2);
                    chi2_chart[1] = modchi2 * sqrt1mxi22 * sin(alpha2);
                    chi2_chart[2] = modchi2 * xi2;
                    RotateVecT(pars->R1, chi1_chart, chi1);
                    RotateVecT(pars->R2, chi2_chart, chi2);
                }

                t = t + dt; // dt is set in the parfile; default is 0.5

            }
        } else { // Standard coordinates with non-canonical spin variables

            if (pars->step == Step_Transformed) { 

                ODESolver(s, w, ds, 14, get_rhs_transformed);

                // Update variables
                for (int i = 0; i < 3; i++) {
                    r[i]    = w[i];
                    p[i]    = w[i + 3];
                    chi1[i] = w[i + 6];
                    chi2[i] = w[i + 9];
                }

                t = w[12];
                s = s + ds;
                // E = w[13];

            } else { // Evolve the standard system

                ODESolver(t, w, dt, 12, get_rhs);

                // Update variables
                for (int i = 0; i < 3; i++) {
                    r[i]    = w[i];
                    p[i]    = w[i + 3];
                    chi1[i] = w[i + 6];
                    chi2[i] = w[i + 9];
                }

                t = t + dt; // dt is set in the parfile; default is 0.5
            }
        }

        // Poincare section: record z = 0 upward crossings (no-op if fps == NULL).
        WriteSectionCrossing(fps, tprev, rprev, pprev, c1prev, c2prev,
                             t, r, p, chi1, chi2);
        tprev = t;
        for (int i = 0; i < 3; i++) { rprev[i]=r[i]; pprev[i]=p[i]; c1prev[i]=chi1[i]; c2prev[i]=chi2[i]; }

    // Evaluate the Hamiltonian and the orbital angular momentum
    Hamiltonian(r, p, pars->nu, chi1, chi2, &dummy, &H, dummyd, dummyd);
    get_l(r, p, l, dummyd2);

    // Q = CarterLikeConstant(r, p, chi1, chi2);
    // CartesianToSpherical(r, p, Q, P); 

    fprintf(fp, "%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n",
            t, r[0], r[1], r[2], p[0], p[1], p[2], chi1[0], chi1[1], chi1[2], chi2[0], chi2[1], chi2[2], H*(pars->nu), l[0], l[1], l[2]);  
    
    // Stop at the light ring or when the radius exceeds the maximum one
    modr = get_mod(r);
    if (modr < rLR) {
        printf("Stop: reached the adiabatic light ring.\n");
        break;
    } else if (modr > rmax) {
        printf("Stop: reached the maximum radius (%.1f).\n", rmax);
        break;        
    }

    // If radius conditions are not met, stop anyway at the maximum integration time
    if (t >= tmax) printf("Stop: maximum time reached.\n");
    
    } while (t < tmax); 

    fclose(fp);
    if (fps != NULL) fclose(fps);
    printf("Data written successfully to %s\n", filepath);

    // Write metadata file (same folder)
    char metadatafile[150];
    snprintf(metadatafile, sizeof(metadatafile), "%s/metadata.txt", folder);
    WriteMetadataFile(pars, metadatafile);

    // Free memory
    FreeParameters(pars);

}

