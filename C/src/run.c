/*
File: run.c

This file contains ...
*/

#include "header.h"

static int HandlePoincareCrossing(const double t_prev,
                                  const double t,
                                  const double r_prev[3],
                                  const double p_prev[3],
                                  const double chi1_prev[3],
                                  const double chi2_prev[3],
                                  const double r[3],
                                  const double p[3],
                                  const double chi1[3],
                                  const double chi2[3],
                                  const int orbit_id,
                                  const double label_r0,
                                  int *cross_id,
                                  FILE *fp_poincare)
{
    if (!pars->poincare_on || fp_poincare == NULL) return 0;

    int coord_idx;
    if (pars->poincare_surface == PoincareSurface_Z) {
        coord_idx = 2;
    } else if (pars->poincare_surface == PoincareSurface_X) {
        coord_idx = 0;
    } else {
        coord_idx = 2;
    }

    double sec_prev = r_prev[coord_idx] - pars->poincare_value;
    double sec_curr = r[coord_idx]      - pars->poincare_value;

    int crossed = ((sec_prev < 0.0 && sec_curr > 0.0) ||
                   (sec_prev > 0.0 && sec_curr < 0.0));

    if (!crossed) return 0;

    double lambda = sec_prev / (sec_prev - sec_curr);

    double t_cross;
    double r_cross[3], p_cross[3], chi1_cross[3], chi2_cross[3];
    double H_cross, dummy_cross, l_cross[3];
    double dHeff_cross[12], dH_cross[12], dummyd2_cross[18];

    t_cross = t_prev + lambda * (t - t_prev);

    for (int i = 0; i < 3; i++) {
        r_cross[i]    = r_prev[i]    + lambda * (r[i]    - r_prev[i]);
        p_cross[i]    = p_prev[i]    + lambda * (p[i]    - p_prev[i]);
        chi1_cross[i] = chi1_prev[i] + lambda * (chi1[i] - chi1_prev[i]);
        chi2_cross[i] = chi2_prev[i] + lambda * (chi2[i] - chi2_prev[i]);
    }

    Hamiltonian(r_cross, p_cross, pars->nu, chi1_cross, chi2_cross,
                &dummy_cross, &H_cross, dHeff_cross, dH_cross);
    get_l(r_cross, p_cross, l_cross, dummyd2_cross);

    double vsec_cross = dH_cross[3 + coord_idx];

    int write_point = 0;
    if (pars->poincare_direction == PoincareDirection_Positive) {
        if (vsec_cross > 0.0) write_point = 1;
    } else if (pars->poincare_direction == PoincareDirection_Negative) {
        if (vsec_cross < 0.0) write_point = 1;
    } else if (pars->poincare_direction == PoincareDirection_Both) {
        write_point = 1;
    }

    if (!write_point) return 0;

    if (pars->scan_on) {
        if (t_cross >= pars->tmax_traj) {
            double rmod_cross = get_mod(r_cross);
            double pr_cross;
            get_pr(r_cross, p_cross, &pr_cross, NULL);

            fprintf(fp_poincare, "%d\t%.16f\t%d\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n",
                    orbit_id, label_r0, *cross_id,
                    t_cross,
                    r_cross[0], r_cross[1], r_cross[2],
                    p_cross[0], p_cross[1], p_cross[2],
                    chi1_cross[0], chi1_cross[1], chi1_cross[2],
                    chi2_cross[0], chi2_cross[1], chi2_cross[2],
                    H_cross*(pars->nu),
                    l_cross[0], l_cross[1], l_cross[2],
                    rmod_cross, pr_cross);

            (*cross_id)++;
            return 1;
        }
    } else {
        fprintf(fp_poincare, "%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n",
                t_cross,
                r_cross[0], r_cross[1], r_cross[2],
                p_cross[0], p_cross[1], p_cross[2],
                chi1_cross[0], chi1_cross[1], chi1_cross[2],
                chi2_cross[0], chi2_cross[1], chi2_cross[2],
                H_cross*(pars->nu),
                l_cross[0], l_cross[1], l_cross[2]);
        return 1;
    }

    return 0;
}

static void InitializeOrbit(const double r0_init[3], const double p0_init[3],
                            double w[14], double W[12],
                            double *t, double *s, double *ds, double *dt,
                            double *tmax_traj, double *tmax_poincare, double *tmax_total,
                            double *modchi1, double *modchi2,
                            double r[3], double p[3], double chi1[3], double chi2[3],
                            double *H, double l[3])
{
    double dummy, dummyd[12], dummyd2[18];

    for (int i = 0; i < 14; i++) w[i] = 0.;
    for (int i = 0; i < 12; i++) W[i] = 0.;

    // Initial phase-space point
    for (int i = 0; i < 3; i++) {
        w[i]     = r0_init[i];
        w[i + 3] = p0_init[i];
    }

    // Initial spins
    w[6]  = pars->chi1x0;
    w[7]  = pars->chi1y0;
    w[8]  = pars->chi1z0;
    w[9]  = pars->chi2x0;
    w[10] = pars->chi2y0;
    w[11] = pars->chi2z0;

    // Integration details
    *t   = 0.;
    *s   = 0.;
    *ds  = 0.0001;
    *dt  = pars->dt;

    *tmax_traj     = pars->tmax_traj;
    *tmax_poincare = pars->tmax_poincare;
    *tmax_total    = (*tmax_traj) + (pars->poincare_on ? (*tmax_poincare) : 0.0);

    // Current physical variables
    for (int i = 0; i < 3; i++) {
        r[i]    = w[i];
        p[i]    = w[i + 3];
        chi1[i] = w[i + 6];
        chi2[i] = w[i + 9];
    }

    Hamiltonian(r, p, pars->nu, chi1, chi2, &dummy, H, dummyd, dummyd);
    get_l(r, p, l, dummyd2);

    *modchi1 = get_mod(chi1);
    *modchi2 = get_mod(chi2);

    // Canonical variables if needed
    const double eps = 1.0e-30;

    if (pars->coords == Coords_Canonical) {
        for (int i = 0; i < 6; i++) W[i] = w[i];
        W[6] = atan2(w[7], w[6]);
        W[7] = w[8] / ((*modchi1) + eps);
        W[8] = atan2(w[10], w[9]);
        W[9] = w[11] / ((*modchi2) + eps);

        if (pars->step == Step_Transformed) {
            W[10] = 0.;
            W[11] = -(*H);
        }
    } else {
        if (pars->step == Step_Transformed) {
            w[12] = 0.;
            w[13] = -(*H);
        }
    }
}

static void AdvanceOrbitOneStep(double w[14], double W[12],
                                double *t, double *s, const double ds, const double dt,
                                const double modchi1, const double modchi2,
                                double r[3], double p[3], double chi1[3], double chi2[3])
{
    double alpha1, xi1, alpha2, xi2, sqrt1mxi12, sqrt1mxi22;

    if (pars->coords == Coords_Canonical) {

        if (pars->step == Step_Transformed) {

            ODESolver(*s, W, ds, 12, get_rhs_canonical_transformed);

            for (int i = 0; i < 3; i++) {
                r[i] = W[i];
                p[i] = W[i + 3];
            }

            alpha1 = W[6];
            xi1    = W[7];
            alpha2 = W[8];
            xi2    = W[9];

            sqrt1mxi12 = sqrt(1. - xi1*xi1);
            sqrt1mxi22 = sqrt(1. - xi2*xi2);

            chi1[0] = modchi1 * sqrt1mxi12 * cos(alpha1);
            chi1[1] = modchi1 * sqrt1mxi12 * sin(alpha1);
            chi1[2] = modchi1 * xi1;

            chi2[0] = modchi2 * sqrt1mxi22 * cos(alpha2);
            chi2[1] = modchi2 * sqrt1mxi22 * sin(alpha2);
            chi2[2] = modchi2 * xi2;

            *t = W[10];
            *s += ds;

        } else {

            ODESolver(*t, W, dt, 10, get_rhs_canonical);

            for (int i = 0; i < 3; i++) {
                r[i] = W[i];
                p[i] = W[i + 3];
            }

            alpha1 = W[6];
            xi1    = W[7];
            alpha2 = W[8];
            xi2    = W[9];

            sqrt1mxi12 = sqrt(1. - xi1*xi1);
            sqrt1mxi22 = sqrt(1. - xi2*xi2);

            chi1[0] = modchi1 * sqrt1mxi12 * cos(alpha1);
            chi1[1] = modchi1 * sqrt1mxi12 * sin(alpha1);
            chi1[2] = modchi1 * xi1;

            chi2[0] = modchi2 * sqrt1mxi22 * cos(alpha2);
            chi2[1] = modchi2 * sqrt1mxi22 * sin(alpha2);
            chi2[2] = modchi2 * xi2;

            *t += dt;
        }

    } else {

        if (pars->step == Step_Transformed) {

            ODESolver(*s, w, ds, 14, get_rhs_transformed);

            for (int i = 0; i < 3; i++) {
                r[i]    = w[i];
                p[i]    = w[i + 3];
                chi1[i] = w[i + 6];
                chi2[i] = w[i + 9];
            }

            *t = w[12];
            *s += ds;

        } else {

            ODESolver(*t, w, dt, 12, get_rhs);

            for (int i = 0; i < 3; i++) {
                r[i]    = w[i];
                p[i]    = w[i + 3];
                chi1[i] = w[i + 6];
                chi2[i] = w[i + 9];
            }

            *t += dt;
        }
    }
}

static void WriteTrajectoryPoint(FILE *fp_dyn,
                                 double t,
                                 const double r[3],
                                 const double p[3],
                                 const double chi1[3],
                                 const double chi2[3],
                                 double H,
                                 const double l[3])
{
    if (fp_dyn == NULL) return;

    fprintf(fp_dyn, "%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n",
            t, r[0], r[1], r[2], p[0], p[1], p[2],
            chi1[0], chi1[1], chi1[2], chi2[0], chi2[1], chi2[2],
            H*(pars->nu), l[0], l[1], l[2]);
}

static int RunOrbit(double r0_init[3], double p0_init[3], int orbit_id, double label_r0,
                    FILE *fp_dyn, FILE *fp_poincare)
{
    // Master variables
    double w[14];
    double W[12];

    // Integration details
    double t, s, ds, dt;
    double tmax_traj, tmax_poincare, tmax_total;

    // Stopping radii
    double rLR;
    AdiabaticLightRing(&rLR);
    double rmax = 500.0;

    // Current physical variables
    double r[3], p[3], chi1[3], chi2[3], H, l[3];
    double modr;

    // Spin magnitudes
    double modchi1, modchi2;

    InitializeOrbit(r0_init, p0_init,
                    w, W,
                    &t, &s, &ds, &dt,
                    &tmax_traj, &tmax_poincare, &tmax_total,
                    &modchi1, &modchi2,
                    r, p, chi1, chi2,
                    &H, l);

    // Write initial dynamical point if requested
    WriteTrajectoryPoint(fp_dyn, t, r, p, chi1, chi2, H, l);

    // Previous state for Poincare crossing detection
    double t_prev = t;
    double r_prev[3], p_prev[3], chi1_prev[3], chi2_prev[3];
    for (int i = 0; i < 3; i++) {
        r_prev[i]    = r[i];
        p_prev[i]    = p[i];
        chi1_prev[i] = chi1[i];
        chi2_prev[i] = chi2[i];
    }

    int cross_id = 0;

    do {
        AdvanceOrbitOneStep(w, W, &t, &s, ds, dt, modchi1, modchi2, r, p, chi1, chi2);

        double dummy, dummyd[12], dummyd2[18];
        Hamiltonian(r, p, pars->nu, chi1, chi2, &dummy, &H, dummyd, dummyd);
        get_l(r, p, l, dummyd2);

        // Write trajectory only during the trajectory phase
        if (t <= tmax_traj) {
            WriteTrajectoryPoint(fp_dyn, t, r, p, chi1, chi2, H, l);
        }

        HandlePoincareCrossing(t_prev, t,
                               r_prev, p_prev, chi1_prev, chi2_prev,
                               r, p, chi1, chi2,
                               orbit_id, label_r0, &cross_id, fp_poincare);

        // Update previous state
        t_prev = t;
        for (int i = 0; i < 3; i++) {
            r_prev[i]    = r[i];
            p_prev[i]    = p[i];
            chi1_prev[i] = chi1[i];
            chi2_prev[i] = chi2[i];
        }

        // Stopping conditions
        modr = get_mod(r);
        if (modr < rLR) {
            printf("Stop: reached the adiabatic light ring.\n");
            break;
        } else if (modr > rmax) {
            printf("Stop: reached the maximum radius (%.1f).\n", rmax);
            break;
        }

        if (t >= tmax_total) printf("Stop: maximum time reached.\n");

    } while (t < tmax_total);

    return cross_id;
}

void BuildOutputFolderName(char *folder, size_t folder_size)
{
    if (pars->motion == Motion_Circular) {
        if (pars->chi1x0 > 1e-15 || pars->chi1y0 > 1e-15 ||
            pars->chi2x0 > 1e-15 || pars->chi2y0 > 1e-15) {

            printf("******************************************************************************************\n");
            printf("Warning: circular ICs + in-plane components of the spins -> motion will be off equatorial!\n");
            printf("******************************************************************************************\n");

            snprintf(folder, folder_size,
                     "circ_q_%.2f_chi1-0_%.2f_%.2f_%.2f_chi2-0_%.2f_%.2f_%.2f",
                     pars->q,
                     pars->chi1x0, pars->chi1y0, pars->chi1z0,
                     pars->chi2x0, pars->chi2y0, pars->chi2z0);
        } else {
            snprintf(folder, folder_size,
                     "circ_q_%.2f_chi1z0_%.2f_chi2z0_%.2f",
                     pars->q, pars->chi1z0, pars->chi2z0);
        }

    } else if (pars->motion == Motion_EquatorialEccentric) {
        if (pars->chi1x0 > 1e-15 || pars->chi1y0 > 1e-15 ||
            pars->chi2x0 > 1e-15 || pars->chi2y0 > 1e-15) {

            printf("******************************************************************************************************\n");
            printf("Warning: equatorial eccentric ICs + in-plane components of the spins -> motion will be off equatorial!\n");
            printf("******************************************************************************************************\n");

            snprintf(folder, folder_size,
                     "eq_q_%.2f_chi1-0_%.2f_%.2f_%.2f_chi2-0_%.2f_%.2f_%.2f",
                     pars->q,
                     pars->chi1x0, pars->chi1y0, pars->chi1z0,
                     pars->chi2x0, pars->chi2y0, pars->chi2z0);
        } else {
            snprintf(folder, folder_size,
                     "eq_q_%.2f_chi1z0_%.2f_chi2z0_%.2f",
                     pars->q, pars->chi1z0, pars->chi2z0);
        }

    } else if (pars->motion == Motion_Generic) {
        if (pars->chi1x0 > 1e-15 || pars->chi1y0 > 1e-15 ||
            pars->chi2x0 > 1e-15 || pars->chi2y0 > 1e-15) {

            snprintf(folder, folder_size,
                     "gen_q_%.2f_chi1-0_%.2f_%.2f_%.2f_chi2-0_%.2f_%.2f_%.2f",
                     pars->q,
                     pars->chi1x0, pars->chi1y0, pars->chi1z0,
                     pars->chi2x0, pars->chi2y0, pars->chi2z0);
        } else {
            snprintf(folder, folder_size,
                     "gen_q_%.2f_chi1z0_%.2f_chi2z0_%.2f",
                     pars->q, pars->chi1z0, pars->chi2z0);
        }

    } else {
        snprintf(folder, folder_size, "output");
    }
}

static FILE *OpenTrajectoryFile(const char *folder, char *filepath, size_t filepath_size)
{
    char local_filepath[150];

    if (filepath != NULL && filepath_size > 0) {
        snprintf(filepath, filepath_size, "%s/dyn.txt", folder);
        printf("Filepath will be: %s\n", filepath);
    } else {
        snprintf(local_filepath, sizeof(local_filepath), "%s/dyn.txt", folder);
        printf("Filepath will be: %s\n", local_filepath);
        filepath = local_filepath;
    }

    FILE *fp = fopen(filepath, "w+");
    if (fp == NULL) {
        perror("Error opening file");
        return NULL;
    }

    fprintf(fp, "0:t\t 1:x\t 2:y\t 3:z\t 4:px\t 5:py\t 6:pz\t 7:chi1x\t 8:chi1y\t 9:chi1z\t 10:chi2x\t 11:chi2y\t 12:chi2z\t 13:E\t 14:lx\t 15:ly\t 16:lz\n");

    return fp;
}

static FILE *OpenPoincareFile(const char *folder, char *filepath, size_t filepath_size)
{
    if (pars->scan_on) {
        snprintf(filepath, filepath_size, "%s/poincare_scan.txt", folder);
    } else {
        snprintf(filepath, filepath_size, "%s/poincare.txt", folder);
    }

    printf("Poincare filepath will be: %s\n", filepath);

    FILE *fp = fopen(filepath, "w+");
    if (fp == NULL) {
        perror("Error opening Poincare file");
        return NULL;
    }

    if (pars->scan_on) {
        fprintf(fp, "0:orbit_id\t1:r0_init\t2:cross_id\t3:t\t4:x\t5:y\t6:z\t7:px\t8:py\t9:pz\t10:chi1x\t11:chi1y\t12:chi1z\t13:chi2x\t14:chi2y\t15:chi2z\t16:E\t17:lx\t18:ly\t19:lz\t20:r\t21:pr\n");
    } else {
        fprintf(fp, "0:t\t 1:x\t 2:y\t 3:z\t 4:px\t 5:py\t 6:pz\t 7:chi1x\t 8:chi1y\t 9:chi1z\t 10:chi2x\t 11:chi2y\t 12:chi2z\t 13:E\t 14:lx\t 15:ly\t 16:lz\n");
    }

    return fp;
}

int RunScanMode(const char *folder)
{
    double x0_input = pars->x0;
    double y0_input = pars->y0;
    double z0_input = pars->z0;

    if (pars->motion != Motion_Generic) {
        fprintf(stderr, "Error: scan mode currently supports only motion = generic.\n");
        return 1;
    }

    if (!pars->poincare_on) {
        fprintf(stderr, "Error: scan mode requires poincare_on = 1.\n");
        return 1;
    }

    if (pars->poincare_surface != PoincareSurface_Z) {
        fprintf(stderr, "Error: scan mode currently supports only poincare_surface = z.\n");
        return 1;
    }

    if (pars->scan_nr < 2) {
        fprintf(stderr, "Error: scan_nr must be at least 2.\n");
        return 1;
    }

    char poincarepath[150];
    FILE *fps = OpenPoincareFile(folder, poincarepath, sizeof(poincarepath));
    if (fps == NULL) {
        return 1;
    }

    double dr = (pars->scan_rmax - pars->scan_rmin) / (double)(pars->scan_nr - 1);

    int orbit_id = 0;
    for (int i = 0; i < pars->scan_nr; i++) {
        double r0_scan = pars->scan_rmin + i * dr;
        double p0_scan[3];

        if (!TryGenericICsAtRadius(r0_scan, p0_scan)) {
            printf("[%d/%d] r0 = %.10f -> invalid ICs\n", i + 1, pars->scan_nr, r0_scan);
            continue;
        }

        double r0_vec[3] = {r0_scan, 0., 0.};
        int n_cross = RunOrbit(r0_vec, p0_scan, orbit_id, r0_scan, NULL, fps);

        printf("[%d/%d] orbit_id = %d, r0 = %.10f, crossings = %d\n",
               i + 1, pars->scan_nr, orbit_id, r0_scan, n_cross);

        orbit_id++;
    }

    fclose(fps);

    if (PostProcessRotationScan(folder, poincarepath) != 0) {
        pars->x0 = x0_input;
        pars->y0 = y0_input;
        pars->z0 = z0_input;
        return 1;
    }

    pars->x0 = x0_input;
    pars->y0 = y0_input;
    pars->z0 = z0_input;

    char metadatafile[150];
    snprintf(metadatafile, sizeof(metadatafile), "%s/metadata.txt", folder);
    WriteMetadataFile(pars, metadatafile);

    printf("Scan finished. Valid orbits: %d\n", orbit_id);

    return 0;
}

int RunSingleMode(const char *folder)
{
    // Open trajectory file only if not scanning
    FILE *fp = OpenTrajectoryFile(folder, NULL, 0);
    if (fp == NULL) {
        return 1;
    }

    // Open the Poincare section file for writing, if requested
    FILE *fps = NULL;
    if (pars->poincare_on) {
        char poincarepath[150];
        fps = OpenPoincareFile(folder, poincarepath, sizeof(poincarepath));
        if (fps == NULL) {
            fclose(fp);
            return 1;
        }
    }

    // Initial conditions for one orbit
    double pphi0 = 0.0, x0 = 0.0, p0[3] = {0.0, 0.0, 0.0};

    if (pars->motion == Motion_Circular) {
        printf("Finding initial conditions for circular motion.\n");
        CircularICs(&pphi0);
    } else if (pars->motion == Motion_EquatorialEccentric) {
        if (pars->E0 > 1) {
            perror("Error: You chose eccentric motion - the energy needs to be lower than 1!\n");
            if (fps != NULL) fclose(fps);
            fclose(fp);
            return 1;
        }
        printf("Finding initial conditions for eccentric equatorial motion.\n");
        EquatorialICs(&x0);
    } else if (pars->motion == Motion_Generic) {
        printf("Finding initial conditions for generic motion.\n");
        GenericICs(p0);
    }

    double r0_vec[3] = {0.0, 0.0, 0.0};

    if (pars->motion == Motion_Circular) {
        r0_vec[0] = pars->x0;
        r0_vec[1] = 0.;
        r0_vec[2] = 0.;

        p0[0] = 0.;
        p0[1] = pphi0 / pars->x0;
        p0[2] = 0.;

    } else if (pars->motion == Motion_EquatorialEccentric) {
        r0_vec[0] = x0;
        r0_vec[1] = 0.;
        r0_vec[2] = 0.;

        p0[0] = 0.;
        p0[1] = pars->pphi0 / x0;
        p0[2] = 0.;

    } else if (pars->motion == Motion_Generic) {
        r0_vec[0] = pars->x0;
        r0_vec[1] = pars->y0;
        r0_vec[2] = pars->z0;
    }

    RunOrbit(r0_vec, p0, 0, r0_vec[0], fp, fps);

    fclose(fp);
    if (fps != NULL) fclose(fps);

    char metadatafile[150];
    snprintf(metadatafile, sizeof(metadatafile), "%s/metadata.txt", folder);
    WriteMetadataFile(pars, metadatafile);

    printf("Data written successfully.\n");

    return 0;
}