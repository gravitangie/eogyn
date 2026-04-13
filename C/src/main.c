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

    // Folder name
    char folder[100];
    BuildOutputFolderName(folder, sizeof(folder));

    // Create the directory
    if (MAKE_DIR(folder) == -1) {
        if (errno != EEXIST) { // Do not give error if dir already exists
            perror("Error creating directory");
            FreeParameters(pars);
            return 1;
        }
    }

    if (pars->scan_on) {
        int status = RunScanMode(folder);
        FreeParameters(pars);
        return status;
    }

    int status = RunSingleMode(folder);
    FreeParameters(pars);
    return status;
}