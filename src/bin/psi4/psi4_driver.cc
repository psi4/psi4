/*! \file psi4_driver.cc
\defgroup PSI4 The new PSI4 driver.
*/
#include <iostream>
#include <fstream>              // file I/O support
#include <libipv1/ip_lib.h>
#include <libparallel/parallel.h>
#include <psi4-dec.h>
#include <string.h>
#include <algorithm>
#include <ctype.h>
#include "psi4.h"

#define MAX_ARGS (20)

using namespace std;

namespace psi {

extern FILE *infile;
extern void setup_driver(Options &options);
extern int read_options(const std::string &jobName, Options &options);

void psiclean(void);

int
psi4_driver(Options & options, int argc, char *argv[])
{
    // Initialize the list of function pointers for each module
    setup_driver(options);

    if (options.get_str("DERTYPE") == "NONE")
        options.set_str("DERTYPE", "ENERGY");

    // Track down the psi.dat file and set up ipv1 to read it
    // unless the user provided an exec array in the input
    std::string psiDatDirName  = Process::environment("PSIDATADIR");
    std::string psiDatFileName = Process::environment("PSIDATADIR") + "/psi.dat";

    FILE * psidat = fopen(psiDatFileName.c_str(), "r");
    if(psidat == NULL){
        throw PsiException("Psi 4 could not locate psi.dat",
                           __FILE__, __LINE__);
    }
    ip_set_uppercase(1);
    ip_append(psidat, stdout);
    ip_cwk_clear();
    ip_cwk_add(const_cast<char*>(":DEFAULT"));
    ip_cwk_add(const_cast<char*>(":PSI"));
    fclose(psidat);

    // Join the job descriptors into one label
    std::string calcType = options.get_str("WFN");
    calcType += ":";
    calcType += options.get_str("DERTYPE");

    if(Communicator::world->me() == 0) {
      std::string wfn = options.get_str("WFN");
      std::string reference = options.get_str("REFERENCE");
      std::string jobtype = options.get_str("JOBTYPE");
      std::string dertype = options.get_str("DERTYPE");
      fprintf(outfile, "    Calculation type = %s\n",calcType.c_str()); fflush(outfile);
      fprintf(outfile, "    Wavefunction     = %s\n",wfn.c_str()); fflush(outfile);
      fprintf(outfile, "    Reference        = %s\n",reference.c_str()); fflush(outfile);
      fprintf(outfile, "    Job type         = %s\n",jobtype.c_str()); fflush(outfile);
      fprintf(outfile, "    Derivative type  = %s\n", dertype.c_str()); fflush(outfile);
    }

    if(check_only)
      fprintf(outfile, "\n    Sanity check requested. Exiting without execution.\n\n");

    char *argv_new[MAX_ARGS];

    // test SCF optimization
    fprintf(outfile,"Hardwired to test SCF optimizations.\n");

    read_options("INPUT", options);
    dispatch_table["INPUT"](options, 0, argv_new);

for (int n=0; n<40; ++n) {

    read_options("CINTS", options);
    dispatch_table["CINTS"](options, 0, argv_new);

    read_options("CSCF", options);
    dispatch_table["CSCF"](options, 0, argv_new);

    read_options("CINTS", options);
    argv_new[0] = new char [strlen("CINTS")+1];
    strcpy(argv_new[0],"CINTS");
    argv_new[1] = new char [strlen("--deriv1")+1];
    strcpy(argv_new[1],"--deriv1");
    dispatch_table["CINTS"](options, 2, argv_new);
    free(argv_new[0]);
    free(argv_new[1]);

    read_options("OPTKING", options);
    dispatch_table["OPT_STEP"](options, 0, argv_new);
}

    return Success;

    char *jobList = const_cast<char*>(calcType.c_str());
    // This version assumes that the array contains only module names, not
    // macros for other job types, like $SCF
    int numTasks = 0;
    int errcod;
    if(ip_exist(const_cast<char*>("EXEC"), 0)){
        // Override psi.dat with the commands in the exec array in input
        errcod = ip_count(const_cast<char*>("EXEC"), &numTasks, 0);
        jobList = const_cast<char*>("EXEC");
    }else{
        errcod = ip_count(jobList, &numTasks, 0);
        if (!ip_exist(jobList, 0)){
            std::string err("Error: jobtype ");
            err += jobList;
            err += " is not defined in psi.dat";
            throw PsiException(err, __FILE__, __LINE__);
        }
        if (errcod != IPE_OK){
            std::string err("Error: trouble reading ");
            err += jobList;
            err += " array from psi.dat";
            throw PsiException(err, __FILE__, __LINE__);
        }
    }

    if(Communicator::world->me() == 0) fprintf(outfile, "    List of tasks to execute:\n");
    for(int n = 0; n < numTasks; ++n) {
        char *thisJob;
        ip_string(jobList, &thisJob, 1, n);
        if(Communicator::world->me() == 0) fprintf(outfile, "    %s\n", thisJob);
        free(thisJob);
    }

    if(check_only) return Success;

    // remove command-line-like arguments from psi.dat string
    // at this point argv[0] is still "psi4" - which we will remove
    int i, argc_new;
    //char *argv_new[MAX_ARGS];
    for (i=1; i<argc; ++i) {
      argv_new[i-1] = new char [strlen(argv[i])+1];
      strcpy(argv_new[i-1],argv[i]);
    }
    --argc;

    // variables to parse string in psi.dat
    string thisJobWithArguments;
    stringstream ss;
    vector<string> tokens;
    string buf;

    for(int n = 0; n < numTasks; ++n){
        char *thisJob;
        errcod = ip_string(jobList, &thisJob, 1, n);

        // tokenize string in psi.dat
        thisJobWithArguments.assign(thisJob);
        ss.clear();
        ss << thisJobWithArguments;
        while (ss >> buf)
          tokens.push_back(buf);

        free(thisJob);
        thisJob = const_cast<char *>(tokens[0].c_str()); // module name for dispatch table

        argc_new = argc + tokens.size();
        for (i=0; i<tokens.size(); ++i)
          argv_new[argc+i] = const_cast<char *>(tokens[i].c_str());

        // Make sure the job name is all upper case
        int length = strlen(thisJob);
        std::transform(thisJob, thisJob + length, thisJob, ::toupper);
        if(Communicator::world->me() == 0) {
          fprintf(outfile, "\n  Job %d is %s\n", n, thisJob); fflush(outfile);
          fprintf(outfile, "  with command-line argument: ");
          for (i=0; i<argc_new; ++i) fprintf(outfile," %s ", argv_new[i]);
        }
        read_options(thisJob, options);
        // If the function call is LMP2, run in parallel
        if(strcmp(thisJob, "LMP2") == 0 || strcmp(thisJob, "DFMP2") == 0) {
            // Needed a barrier before the functions are called
            Communicator::world->sync();
            if (dispatch_table[thisJob](options, argc_new, argv_new) != Success) {
                // Good chance at this time that an error occurred.
                // Report it to the user.
                fprintf(stderr, "%s did not return a Success code.\n", thisJob);
                throw PsiException("Module failed.", __FILE__, __LINE__);
            }
        }
        // If any other functions are called only process 0 runs the function
        else {
            if (dispatch_table.find(thisJob) != dispatch_table.end()) {
                // Needed a barrier before the functions are called
                Communicator::world->sync();
                if(Communicator::world->me() == 0) {
                    if (dispatch_table[thisJob](options, argc_new, argv_new) != Success) {
                        // Good chance at this time that an error occurred.
                        // Report it to the user.
                        fprintf(stderr, "%s did not return a Success code.\n", thisJob);
                        throw PsiException("Module failed.", __FILE__, __LINE__);
                    }
                }
            }
            else {
                std::transform(thisJob, thisJob + length, thisJob, ::tolower);
                // Close the output file.
                if (!outfile_name.empty())
                    fclose(outfile);
                // Attempt to run the external program
                int ret = ::system(thisJob);
                if (!outfile_name.empty())
                    fopen(outfile_name.c_str(), "a");

                if (ret == -1 || ret == 127) {
                    std::string err = "Module ";
                    err += thisJob;
                    err += " is not known to PSI4.  Please update the driver\n";
                    throw PsiException(err, __FILE__, __LINE__);
                }
            }
        }

        tokens.clear();

        fflush(outfile);
    }

    if (!messy) {
        if(Communicator::world->me() == 0)
            psiclean();
    }


    return Success;
}

} // Namespace
