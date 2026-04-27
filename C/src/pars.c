/*
File: pars.c

This file contains code responsible for managing all parameters.
*/

#include "header.h"

// Global variable for parameters:
// - firstly declared as extern in the header
// - the structure type Parameters is also defined in the header 
Parameters *pars;

// Global variables for function pointers: 
// - needed for functions that have different options 
// - firstly declared as extern in the header 
void (*Potentials)();
void (*ODESolver)();


// Allocate space for the parameters
void AllocateParameters(Parameters **pars)
{
    *pars = (Parameters *) calloc(1, sizeof(Parameters)); 
    // calloc initialises to 0; if there is not memory left, it returns NULL.
    if (pars == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);  
    }
} 

void FreeParameters(Parameters *pars)
{
    if (pars != NULL) free(pars);
}


// Set the default values for the parameters
void SetDefaults(Parameters *pars)
{
    // Binary parameters
    pars->q      = 1.;
    pars->chi1x0 = 0.;
    pars->chi1y0 = 0.;
    pars->chi1z0 = 0.;
    pars->chi2x0 = 0.;
    pars->chi2y0 = 0.;
    pars->chi2z0 = 0.;

    // Orbital params.
    pars->x0 = 10.;
    pars->y0 = 0.;
    pars->z0 = 0.;

    // Type of motion
    pars->motion = Motion_Generic;

    // Constants of motion 
    // Initialized as NAN to check later if they have been defined in the input parfile or not
    pars->pphi0 = NAN;
    pars->E0  = NAN;

    // Options for functions within the Hamiltonian
    pars->pots  = Pots_Resummed;

    // Type of coordinates (standard with non-canonical spins, or canonical)
    pars->coords = Coords_Canonical; 

    // ODE solver settings
    pars->solver         = ODESolver_RKGL6;
    pars->step           = Step_Transformed; // FIXME: Could think of turning this on only for eccentric
    pars->dt             = 0.5;
    pars->tmax_traj      = 1000.;
    pars->tmax_poincare  = 0.;
    pars->max_iter_RKGL6 = 100.;
    pars->tol_RKGL6      = 1e-15;

    // Poincare section settings
    pars->poincare_on        = 0;
    pars->poincare_surface   = PoincareSurface_Z;
    pars->poincare_direction = PoincareDirection_Positive;
    pars->poincare_value     = 0.;

    // Poincare scan settings
    pars->scan_on   = 0;
    pars->scan_rmin = 5.0;
    pars->scan_rmax = 20.0;
    pars->scan_nr   = 100;

    // Rotation-number post-processing settings
    pars->rotation_min_crossings_for_center = 50;
    pars->rotation_n_center_orbits          = 5;

}

// Parse the command line to get the name of the parfile & read it 
void ParseCommandLine(int argc, char *argv[], Parameters *pars) 
{
    int opt;
    while ((opt = getopt(argc, argv, "p:")) != -1) {
        switch (opt) {
            case 'p':
                ReadParfile(optarg, pars);
                break;
            default:
                fprintf(stderr, "Usage: %s [-p parfile]\n", argv[0]);  //[-param1 value] [-param2 value]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }
}

// Read the parfile & assign the values to the related variables
void ReadParfile(char *filename, Parameters *pars) 
{
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open parameter file.\n");
        return;
    }

    char line[256];
    while (fgets(line, sizeof(line), file)) {
        // Remove comments
        char *comment = strchr(line, '#');
        if (comment) {
            *comment = '\0';  // Terminate the string at the start of the comment
        }

        // Trim whitespace
        char *key = strtok(line, " =\t\n");
        char *value = strtok(NULL, " =\t\n");

        // Assign values if key and value are present
        if (key && value) {
            AssignValues(pars, key, value);
        }
    }
    fclose(file);
}


// Assign the values read from the parfile
void AssignValues(Parameters *pars, char *key, char *value) 
{
    if (strcmp(key, "q") == 0) {
        pars->q = atof(value);
    } else if (strcmp(key, "chi1x0") == 0) {
        pars->chi1x0 = atof(value);
    } else if (strcmp(key, "chi1y0") == 0) {
        pars->chi1y0 = atof(value);
    } else if (strcmp(key, "chi1z0") == 0) {
        pars->chi1z0 = atof(value);
    } else if (strcmp(key, "chi2x0") == 0) {
        pars->chi2x0 = atof(value);
    } else if (strcmp(key, "chi2y0") == 0) {
        pars->chi2y0 = atof(value);
    } else if (strcmp(key, "chi2z0") == 0) {
        pars->chi2z0 = atof(value);
    } else if (strcmp(key, "x0") == 0) {
        pars->x0 = atof(value);
    } else if (strcmp(key, "y0") == 0) {
        pars->y0 = atof(value);
    } else if (strcmp(key, "z0") == 0) {
        pars->z0 = atof(value);
    } else if (strcmp(key, "E0") == 0) {
        pars->E0 = atof(value);
    } else if (strcmp(key, "pphi0") == 0) {
        pars->pphi0 = atof(value);
    } else if (strcmp(key, "dt") == 0) {
        pars->dt = atof(value);
    } else if (strcmp(key, "tmax_traj") == 0) {
        pars->tmax_traj = atof(value);
    } else if (strcmp(key, "tmax_poincare") == 0) {
        pars->tmax_poincare = atof(value);
    } else if (strcmp(key, "max_iter_RKGL6") == 0) {
        pars->max_iter_RKGL6 = atof(value);
    } else if (strcmp(key, "tol_RKGL6") == 0) {
        pars->tol_RKGL6 = atof(value);
    } else if (strcmp(key, "coords") == 0) {
        trimString(value);
        for (pars->coords = 0; pars->coords <= Coords_NOPT; pars->coords++) { 
          if (pars->coords == Coords_NOPT) {
            pars->coords = Coords_Canonical; // default
            printf("No option chosen for the coordinates, thus set to the default: '%s'\n", coords_opt[pars->coords]);
            break;
          }
          if (strcmp(value, coords_opt[pars->coords]) == 0) break;
        }
    } else if (strcmp(key, "solver") == 0) {
        trimString(value);
        for (pars->solver = 0; pars->solver <= ODESolver_NOPT; pars->solver++) { 
          if (pars->solver == ODESolver_NOPT) {
            pars->solver = ODESolver_RKGL6; // default
            printf("No option chosen for the ODE solver, thus set to the default: '%s'\n", solver_opt[pars->solver]);
            break;
          }
          if (strcmp(value, solver_opt[pars->solver]) == 0) break;
        }
    } else if (strcmp(key, "step") == 0) {
        trimString(value);
        for (pars->step = 0; pars->step <= Step_NOPT; pars->step++) { 
          if (pars->step == Step_NOPT) {
            pars->step = Step_Standard; // default
            printf("No option chosen for the time step, thus set to the default: '%s'\n", step_opt[pars->step]);
            break;
          }
          if (strcmp(value, step_opt[pars->step]) == 0) break;
        }
    } else if (strcmp(key, "pots") == 0) {
        trimString(value);
        for (pars->pots = 0; pars->pots <= Pots_NOPT; pars->pots++) { 
          if (pars->pots == Pots_NOPT) {
            pars->pots = Pots_Resummed; // default
            printf("No option chosen for the potentials, thus set to the default: '%s'\n", pots_opt[pars->pots]);
            break;
          }
          if (strcmp(value, pots_opt[pars->pots]) == 0) break;
        }
    } else if (strcmp(key, "motion") == 0) {
        trimString(value);
        for (pars->motion = 0; pars->motion <= Motion_NOPT; pars->motion++) {
          if (pars->motion == Motion_NOPT) {
            pars->motion = Motion_Generic; // default
            printf("No option chosen for the type of motion, thus set to the default: '%s'\n", motion_opt[pars->motion]);
            // Give error if E0 and pphi0 have not been set 
            if (isnan(pars->E0) || isnan(pars->pphi0)) {
                printf("Error: the initial energy and angular momentum must be set for generic motion.\n");
                exit(EXIT_FAILURE);
            } else {
            break;
            }
          }
          if (strcmp(value, motion_opt[pars->motion]) == 0) break;
        }
    } else if (strcmp(key, "poincare_on") == 0) {
        pars->poincare_on = atoi(value);
    } else if (strcmp(key, "poincare_value") == 0) {
        pars->poincare_value = atof(value);
    } else if (strcmp(key, "poincare_surface") == 0) {
        trimString(value);
        for (pars->poincare_surface = 0; pars->poincare_surface <= PoincareSurface_NOPT; pars->poincare_surface++) {
            if (pars->poincare_surface == PoincareSurface_NOPT) {
                pars->poincare_surface = PoincareSurface_Z; // default
                printf("No option chosen for the Poincare surface, thus set to the default: '%s'\n",
                       poincare_surface_opt[pars->poincare_surface]);
                break;
            }
            if (strcmp(value, poincare_surface_opt[pars->poincare_surface]) == 0) break;
        }
    } else if (strcmp(key, "poincare_direction") == 0) {
        trimString(value);
        for (pars->poincare_direction = 0; pars->poincare_direction <= PoincareDirection_NOPT; pars->poincare_direction++) {
            if (pars->poincare_direction == PoincareDirection_NOPT) {
                pars->poincare_direction = PoincareDirection_Positive; // default
                printf("No option chosen for the Poincare direction, thus set to the default: '%s'\n",
                       poincare_direction_opt[pars->poincare_direction]);
                break;
            }
            if (strcmp(value, poincare_direction_opt[pars->poincare_direction]) == 0) break;
        }
    } else if (strcmp(key, "scan_on") == 0) {
        pars->scan_on = atoi(value);
    } else if (strcmp(key, "scan_rmin") == 0) {
        pars->scan_rmin = atof(value);
    } else if (strcmp(key, "scan_rmax") == 0) {
        pars->scan_rmax = atof(value);
    } else if (strcmp(key, "scan_nr") == 0) {
        pars->scan_nr = atoi(value);
    } else if (strcmp(key, "rotation_min_crossings_for_center") == 0) {
        pars->rotation_min_crossings_for_center = atoi(value);
    } else if (strcmp(key, "rotation_n_center_orbits") == 0) {
        pars->rotation_n_center_orbits = atoi(value);
    }
}

// Set additional parameters not determined in parfile
void SetParameters()
{

    // Evaluate parameters that depend on the input ones
    pars->nu = get_nu(pars->q);
    double X1, X2;
      get_X1X2(pars->nu, &X1, &X2);
    pars->X1 = X1;
    pars->X2 = X2;

    // For circular or generic motion, start along the x axis
    if ((pars->motion == Motion_Circular || pars->motion == Motion_Generic) && (pars->y0 > 1e-15)) {
        printf("For circular/generic motion, setting the initial orbital phase to zero and the in-plane component of the radius along the x axis.\n");
        double modr = sqrt((pars->x0)*(pars->x0) + (pars->y0)*(pars->y0));
        pars->x0 = modr;
        pars->y0 = 0.;
        printf("New x0 = %.10f, new y0 = 0.\n", pars->x0);
    }

    // Set the function pointer for the integration scheme
    if (pars->solver == ODESolver_RK4) {
      ODESolver = &RK_Fourth;
    } else if (pars->solver == ODESolver_RKGL6) {
      ODESolver = &RK_GaussLegendre6;
    }

    // Set the function pointer for the potentials
    if (pars->pots == Pots_NonResummed) {
      Potentials = &NonResummedPotentials;
    } else if (pars->pots == Pots_Resummed) {
      Potentials = &ResummedPotentials;
    }
  
    // Choose type of motion depending on the input parameters:
    // if on the equatorial plane & spin components only along the z axis, set motion to circular
    /*
    if (pars->z0 < 1e-15 & pars->chi1x0 < 1e-15 & pars->chi1y0 < 1e-15 & pars->chi2x0 < 1e-15 & pars->chi2y0 < 1e-15) {
        pars->motion = Motion_Circular; 

        if (pars->y0 > 1e-15) {
            printf("For circular motion, setting the initial orbital phase to zero and the radius along the x axis.\n");
            double modr = sqrt((pars->x0)*(pars->x0) + (pars->y0)*(pars->y0));
            pars->x0 = modr;
            pars->y0 = 0.;
            printf("New x0 = %.10f, new y0 = 0.\n", pars->x0);
        }
    } else { // Otherwise, set to generic motion
        pars->motion = Motion_Generic; 
        if (isnan(pars->E0) || isnan(pars->pphi0)) {
            printf("Error: the initial energy and angular momentum must be set for generic motion.\n");
            exit(EXIT_FAILURE);
        }
    }*/

    // Check if E0, pphi0 have been set in the parameter file
    /*
    if (isnan(pars->E0) || isnan(pars->pphi0)) {
        pars->motion = Motion_Circular; 

        // Motion starts with initial orbital phase zero and radius along the x axis
        if (pars->y0 > 1e-15) {
            printf("For circular motion, setting the initial orbital phase to zero and the radius along the x axis.\n");
            double modr = sqrt((pars->x0)*(pars->x0) + (pars->y0)*(pars->y0));
            pars->x0 = modr;
            pars->y0 = 0.;
            printf("New x0 = %.10f, new y0 = 0.\n", pars->x0);
        }

    }
    */

}

void WriteMetadataFile(Parameters *pars, char *filepath) 
{
    // Open the file for writing
    FILE* fp;
    fp = fopen(filepath, "w+");
    if (fp == NULL) {
        perror("Error opening metadata file\n");
        exit(EXIT_FAILURE);
    }

    // Binary parameters
    fprintf(fp, "q = %.16f\n", pars->q);
    fprintf(fp, "nu = %.16f\n", pars->nu); 
    fprintf(fp, "chi1x0 = %.16f\n", pars->chi1x0); 
    fprintf(fp, "chi1y0 = %.16f\n", pars->chi1y0); 
    fprintf(fp, "chi1z0 = %.16f\n", pars->chi1z0); 
    fprintf(fp, "chi2x0 = %.16f\n", pars->chi2x0); 
    fprintf(fp, "chi2y0 = %.16f\n", pars->chi2y0); 
    fprintf(fp, "chi2z0 = %.16f\n", pars->chi2z0); 

    // Type of motion
    fprintf(fp, "motion = %s\n", motion_opt[pars->motion]); 

    // Orbital parameters 
    if (pars->motion == Motion_Circular) {
        fprintf(fp, "x0 = %.16f\n", pars->x0); 
        fprintf(fp, "y0 = %.16f\n", pars->y0);
    } else if (pars->motion == Motion_EquatorialEccentric) {
        fprintf(fp, "E0 = %.16f\n", pars->E0); 
        fprintf(fp, "pphi0 = %.16f\n", pars->pphi0);
    } else if (pars->motion == Motion_Generic) {
        fprintf(fp, "x0 = %.16f\n", pars->x0); 
        fprintf(fp, "y0 = %.16f\n", pars->y0);
        fprintf(fp, "z0 = %.16f\n", pars->z0);
        fprintf(fp, "E0 = %.16f\n", pars->E0); 
        fprintf(fp, "pphi0 = %.16f\n", pars->pphi0);
    } 

    // Options for functions within the Hamiltonian
    fprintf(fp, "pots = %s\n", pots_opt[pars->pots]); 

    // Type of coordinates
    fprintf(fp, "coords = %s\n", coords_opt[pars->coords]); 

    // ODE solver settings
    fprintf(fp, "solver = %s\n", solver_opt[pars->solver]); 
    if (pars->solver == ODESolver_RKGL6) {
        fprintf(fp, "max_iter_RKGL6 = %.16f\n", pars->max_iter_RKGL6);
        fprintf(fp, "tol_RKGL6 = %.16f\n", pars->tol_RKGL6);
    }
    fprintf(fp, "step = %s\n", step_opt[pars->step]); 
    fprintf(fp, "dt = %.16f\n", pars->dt); 
    fprintf(fp, "tmax_traj = %.16f\n", pars->tmax_traj);
    fprintf(fp, "tmax_poincare = %.16f\n", pars->tmax_poincare);

    // Poincare section settings
    fprintf(fp, "poincare_on = %d\n", pars->poincare_on);
    fprintf(fp, "poincare_surface = %s\n", poincare_surface_opt[pars->poincare_surface]);
    fprintf(fp, "poincare_direction = %s\n", poincare_direction_opt[pars->poincare_direction]);
    fprintf(fp, "poincare_value = %.16f\n", pars->poincare_value);

    // Poincare scan settings
    fprintf(fp, "scan_on = %d\n", pars->scan_on);
    fprintf(fp, "scan_rmin = %.16f\n", pars->scan_rmin);
    fprintf(fp, "scan_rmax = %.16f\n", pars->scan_rmax);
    fprintf(fp, "scan_nr = %d\n", pars->scan_nr);

    // Rotation-number post-processing settings
    fprintf(fp, "rotation_min_crossings_for_center = %d\n",
            pars->rotation_min_crossings_for_center);
    fprintf(fp, "rotation_n_center_orbits = %d\n",
            pars->rotation_n_center_orbits);

}