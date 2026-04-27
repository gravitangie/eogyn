/*
File: rotation.c

This file contains post-processing for rotation numbers computed from Poincare scan data.
*/

#include "header.h"

typedef struct {
    int orbit_id;
    double r0_init;

    int n_cross;
    int capacity;

    double *rsec;
    double *prsec;

    double r_min;
    double r_max;
    double r_mean;
    double pr_mean;
    double delta_r;

    double rho;
    int valid_for_center;
} RotationOrbit;

typedef struct {
    int idx;
    double delta_r;
} CenterCandidate;

static void *CheckedRealloc(void *ptr, size_t size)
{
    void *out = realloc(ptr, size);
    if (out == NULL) {
        fprintf(stderr, "Error: realloc failed in rotation post-processing.\n");
        exit(EXIT_FAILURE);
    }
    return out;
}

static double wrap_to_2pi(double angle)
{
    angle = fmod(angle, 2.0 * Pi);
    if (angle < 0.0) angle += 2.0 * Pi;
    return angle;
}

static void InitRotationOrbit(RotationOrbit *orb, int orbit_id, double r0_init)
{
    orb->orbit_id = orbit_id;
    orb->r0_init  = r0_init;

    orb->n_cross  = 0;
    orb->capacity = 0;

    orb->rsec  = NULL;
    orb->prsec = NULL;

    orb->r_min = 0.0;
    orb->r_max = 0.0;
    orb->r_mean = 0.0;
    orb->pr_mean = 0.0;
    orb->delta_r = 0.0;

    orb->rho = NAN;
    orb->valid_for_center = 0;
}

static void FreeRotationOrbit(RotationOrbit *orb)
{
    if (orb->rsec  != NULL) free(orb->rsec);
    if (orb->prsec != NULL) free(orb->prsec);

    orb->rsec = NULL;
    orb->prsec = NULL;
    orb->n_cross = 0;
    orb->capacity = 0;
}

static int FindOrbitIndex(const RotationOrbit *orbits, int norbits, int orbit_id)
{
    for (int i = 0; i < norbits; i++) {
        if (orbits[i].orbit_id == orbit_id) return i;
    }
    return -1;
}

static void AppendPointToOrbit(RotationOrbit *orb, double r, double pr)
{
    if (orb->n_cross >= orb->capacity) {
        int new_capacity = (orb->capacity == 0) ? 128 : 2 * orb->capacity;
        orb->rsec  = (double *)CheckedRealloc(orb->rsec,  new_capacity * sizeof(double));
        orb->prsec = (double *)CheckedRealloc(orb->prsec, new_capacity * sizeof(double));
        orb->capacity = new_capacity;
    }

    orb->rsec[orb->n_cross]  = r;
    orb->prsec[orb->n_cross] = pr;
    orb->n_cross++;
}

static void FinalizeOrbitStats(RotationOrbit *orb)
{
    if (orb->n_cross <= 0) {
        orb->r_mean = NAN;
        orb->pr_mean = NAN;
        orb->delta_r = NAN;
        orb->rho = NAN;
        orb->valid_for_center = 0;
        return;
    }

    double rmin = orb->rsec[0];
    double rmax = orb->rsec[0];
    double sumr = 0.0;
    double sumpr = 0.0;

    for (int i = 0; i < orb->n_cross; i++) {
        double r  = orb->rsec[i];
        double pr = orb->prsec[i];

        if (r < rmin) rmin = r;
        if (r > rmax) rmax = r;

        sumr  += r;
        sumpr += pr;
    }

    orb->r_min   = rmin;
    orb->r_max   = rmax;
    orb->r_mean  = sumr  / (double)orb->n_cross;
    orb->pr_mean = sumpr / (double)orb->n_cross;
    orb->delta_r = rmax - rmin;

    orb->valid_for_center = (orb->n_cross >= pars->rotation_min_crossings_for_center);
}

static int CompareCandidates(const void *a, const void *b)
{
    const CenterCandidate *ca = (const CenterCandidate *)a;
    const CenterCandidate *cb = (const CenterCandidate *)b;

    if (ca->delta_r < cb->delta_r) return -1;
    if (ca->delta_r > cb->delta_r) return 1;
    return 0;
}

static int EstimateIslandCenter(const RotationOrbit *orbits, int norbits,
                                double *center_r, double *center_pr, int *n_used)
{
    CenterCandidate *cand = NULL;
    int ncand = 0;
    int cap = 0;

    for (int i = 0; i < norbits; i++) {
        if (!orbits[i].valid_for_center) continue;

        if (ncand >= cap) {
            int new_cap = (cap == 0) ? 16 : 2 * cap;
            cand = (CenterCandidate *)CheckedRealloc(cand, new_cap * sizeof(CenterCandidate));
            cap = new_cap;
        }

        cand[ncand].idx = i;
        cand[ncand].delta_r = orbits[i].delta_r;
        ncand++;
    }

    if (ncand == 0) {
        free(cand);
        return 0;
    }

    qsort(cand, ncand, sizeof(CenterCandidate), CompareCandidates);

    int requested = pars->rotation_n_center_orbits;
    if (requested < 1) requested = 1;

    int use = (ncand < requested) ? ncand : requested;

    double sum_r = 0.0;
    double sum_pr = 0.0;

    for (int i = 0; i < use; i++) {
        int idx = cand[i].idx;
        sum_r  += orbits[idx].r_mean;
        sum_pr += orbits[idx].pr_mean;
    }

    *center_r  = sum_r  / (double)use;
    *center_pr = sum_pr / (double)use;
    *n_used = use;

    free(cand);
    return 1;
}

static double ComputeRotationNumber(const double *rsec,
                                    const double *prsec,
                                    int npts,
                                    double center_r,
                                    double center_pr)
{
    if (npts < 2) return NAN;

    double total_angle = 0.0;

    for (int i = 0; i < npts - 1; i++) {
        double v1x = rsec[i]     - center_r;
        double v1y = prsec[i]    - center_pr;
        double v2x = rsec[i + 1] - center_r;
        double v2y = prsec[i + 1]- center_pr;

        double theta1 = atan2(v1y, v1x);
        double theta2 = atan2(v2y, v2x);

        double gamma = wrap_to_2pi(theta1 - theta2);
        total_angle += gamma;
    }

    return total_angle / (2.0 * Pi * (double)(npts - 1));
}

static int ReadPoincareScan(const char *filepath, RotationOrbit **orbits_out, int *norbits_out)
{
    FILE *fp = fopen(filepath, "r");
    if (fp == NULL) {
        perror("Error opening poincare_scan.txt for rotation post-processing");
        return 1;
    }

    RotationOrbit *orbits = NULL;
    int norbits = 0;
    int orbit_cap = 0;

    char line[4096];

    while (fgets(line, sizeof(line), fp) != NULL) {
        int orbit_id, cross_id;
        double r0_init, t;
        double x, y, z, px, py, pz;
        double chi1x, chi1y, chi1z, chi2x, chi2y, chi2z;
        double E, lx, ly, lz, r, pr;

        int nread = sscanf(line,
                           "%d%lf%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                           &orbit_id, &r0_init, &cross_id, &t,
                           &x, &y, &z, &px, &py, &pz,
                           &chi1x, &chi1y, &chi1z,
                           &chi2x, &chi2y, &chi2z,
                           &E, &lx, &ly, &lz,
                           &r, &pr);

        if (nread != 22) {
            continue; // header or malformed line
        }

        int idx = FindOrbitIndex(orbits, norbits, orbit_id);
        if (idx < 0) {
            if (norbits >= orbit_cap) {
                int new_cap = (orbit_cap == 0) ? 32 : 2 * orbit_cap;
                orbits = (RotationOrbit *)CheckedRealloc(orbits, new_cap * sizeof(RotationOrbit));
                orbit_cap = new_cap;
            }

            InitRotationOrbit(&orbits[norbits], orbit_id, r0_init);
            idx = norbits;
            norbits++;
        }

        AppendPointToOrbit(&orbits[idx], r, pr);
    }

    fclose(fp);

    for (int i = 0; i < norbits; i++) {
        FinalizeOrbitStats(&orbits[i]);
    }

    *orbits_out = orbits;
    *norbits_out = norbits;

    return 0;
}

static int WriteRotationScanFile(const char *folder,
                                 const RotationOrbit *orbits,
                                 int norbits,
                                 double center_r,
                                 double center_pr)
{
    char filepath[256];
    snprintf(filepath, sizeof(filepath), "%s/rotation_scan.txt", folder);

    FILE *fp = fopen(filepath, "w+");
    if (fp == NULL) {
        perror("Error opening rotation_scan.txt");
        return 1;
    }

    fprintf(fp, "# estimated_center_r = %.16f\n", center_r);
    fprintf(fp, "# estimated_center_pr = %.16f\n", center_pr);
    fprintf(fp, "0:orbit_id\t1:r0_init\t2:n_cross\t3:r_min\t4:r_max\t5:delta_r\t6:r_mean\t7:pr_mean\t8:center_r\t9:center_pr\t10:rho\n");

    for (int i = 0; i < norbits; i++) {
        fprintf(fp,
                "%d\t%.16f\t%d\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n",
                orbits[i].orbit_id,
                orbits[i].r0_init,
                orbits[i].n_cross,
                orbits[i].r_min,
                orbits[i].r_max,
                orbits[i].delta_r,
                orbits[i].r_mean,
                orbits[i].pr_mean,
                center_r,
                center_pr,
                orbits[i].rho);
    }

    fclose(fp);
    return 0;
}

int PostProcessRotationScan(const char *folder, const char *poincare_scan_path)
{
        if (pars->rotation_min_crossings_for_center < 1) {
        fprintf(stderr, "Error: rotation_min_crossings_for_center must be >= 1.\n");
        return 1;
    }

    if (pars->rotation_n_center_orbits < 1) {
        fprintf(stderr, "Error: rotation_n_center_orbits must be >= 1.\n");
        return 1;
    }

    RotationOrbit *orbits = NULL;
    int norbits = 0;

    if (ReadPoincareScan(poincare_scan_path, &orbits, &norbits) != 0) {
        return 1;
    }

    if (norbits == 0) {
        fprintf(stderr, "Error: no valid orbit data found in poincare_scan.txt.\n");
        free(orbits);
        return 1;
    }

    double center_r, center_pr;
    int n_used = 0;

    if (!EstimateIslandCenter(orbits, norbits, &center_r, &center_pr, &n_used)) {
        fprintf(stderr, "Error: could not estimate island center from scan data.\n");
        for (int i = 0; i < norbits; i++) FreeRotationOrbit(&orbits[i]);
        free(orbits);
        return 1;
    }

    printf("Estimated island center from %d innermost orbit(s): r = %.16f, pr = %.16f\n",
           n_used, center_r, center_pr);

    for (int i = 0; i < norbits; i++) {
        orbits[i].rho = ComputeRotationNumber(orbits[i].rsec, orbits[i].prsec,
                                              orbits[i].n_cross, center_r, center_pr);
    }

    int status = WriteRotationScanFile(folder, orbits, norbits, center_r, center_pr);

    for (int i = 0; i < norbits; i++) {
        FreeRotationOrbit(&orbits[i]);
    }
    free(orbits);

    return status;
}