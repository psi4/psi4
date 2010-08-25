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

  extern PsiReturnType execute_bp(std::string & bp, Options & options);
  extern void setup_driver(Options & options);
  extern int read_options(const std::string &jobName, Options &options, bool call_ipv1 = true,
    bool suppress_printing = false);
  extern void psiclean(void);

  int psi4_driver()
  {
      Options& options = Process::environment.options;

      // Initialize the list of function pointers for each module
      setup_driver(options);

      // Track down the psi.dat file and set up ipv1 to read it
      // unless the user provided an exec array in the input
      std::string psiDatDirName  = Process::environment("PSIDATADIR");
      std::string psiDatFileName = Process::environment("PSIDATADIR") + "/psi.dat";

      FILE * psidat = fopen(psiDatFileName.c_str(), "r");
      if(psidat == NULL){
          throw PsiException("Psi 4 could not locate psi.dat", __FILE__, __LINE__);
      }
      ip_set_uppercase(1);
      ip_append(psidat, stdout);
      ip_cwk_clear();
      ip_cwk_add(const_cast<char*>(":DEFAULT"));
      ip_cwk_add(const_cast<char*>(":PSI"));
      fclose(psidat);

      std::string jobtype = options.get_str("JOBTYPE");
      std::string wfn = options.get_str("WFN");
      std::string reference = options.get_str("REFERENCE");
      std::string dertype = options.get_str("DERTYPE");
      std::string calcType;

      // Join the job descriptors into one label
      calcType = wfn;
      calcType += ":";
      calcType += reference;
      calcType += ":";
      calcType +=  dertype;

      if(Communicator::world->me() == 0) {
          fprintf(outfile, "    Job type         = %s\n", jobtype.c_str());   fflush(outfile);
          fprintf(outfile, "    Wavefunction     = %s\n", wfn.c_str());       fflush(outfile);
          fprintf(outfile, "    Reference        = %s\n", reference.c_str()); fflush(outfile);
          fprintf(outfile, "    Derivative type  = %s\n", dertype.c_str());   fflush(outfile);
          fprintf(outfile, "    Calculation type = %s\n", calcType.c_str());  fflush(outfile);
      }

      if(check_only) fprintf(outfile, "\n    Sanity check requested. Exiting without execution.\n\n");

      if (jobtype == "SP") {
          if (!check_only) {
              read_options("INPUT", options);
              dispatch_table["INPUT"](options);
          }
          else fprintf(outfile,"    Tasks to compute enegy: input\n");

          std::string bp_name = options.get_str("WFN") + ":" + "ENERGY";
          execute_bp(bp_name, options);
      }
      else if (jobtype == "OPT") {
          if (dertype == "FIRST") {
              std::string bp_name = options.get_str("WFN") + ":" + "FIRST";

              if (check_only) {
                  fprintf(outfile,"    Geometry optimization requested.\n");

                  read_options("OPTKING", options, true, true); // don't print
                  int nopt = options.get_int("NOPT");
                  fprintf(outfile,"    Tasks include input, up to %d gradients.\n", nopt);

                  fprintf(outfile,"    To compute gradient:\n");
                  execute_bp(bp_name, options);
              }
              else {
                  read_options("INPUT", options);
                  dispatch_table["INPUT"](options);

                  read_options("OPTKING", options, true, true); // don't print
                  int nopt = options.get_int("NOPT");
                  PsiReturnType rval;

                  for (int i=0; i<nopt; ++i) {
                      execute_bp(bp_name, options);
                      read_options("OPTKING", options);
                      rval = dispatch_table["OPT_STEP"](options);
                      if (rval == Endloop) break;
                  }
              }  // do optimization by gradients
          } // dertype = first
          else {
              fprintf(outfile,"Can only do optimizations by gradients at this time.\n");
          }
      }

      // if (!messy) ???
      if (!check_only && Communicator::world->me() == 0)
          psiclean();

      return Success;
  }

} // Namespace

