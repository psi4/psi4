#ifndef __psi_psi4_psi4_h_
#define __psi_psi4_psi4_h_

#include <stdio.h>
#include <string>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>

#ifdef MAIN
#define EXT
#else
#define EXT extern
#endif

namespace psi {

// Define Calculation class - not used for now

class Calculation
{
  enum Jobtype { sp, opt, freq } jobtype;
  enum Wfn { SCF, CCSD, G3 } wfn;
  enum Dertype { none, first, second, unspecified } dertype;

  // sets values of calculation from cstrings - order must match the
  // enumerated lists above - could make more elegant later 
  void set(char *sj, char *sw, char *sd) {
    if      (!strcmp(sj,"SP"))   jobtype = sp;
    else if (!strcmp(sj,"OPT"))  jobtype = opt;
    else if (!strcmp(sj,"FREQ")) jobtype = freq;
    else throw("Calculation::set(): JOBTYPE unknown\n");

    if      (!strcmp(sw,"SCF"))   wfn = SCF;
    else if (!strcmp(sw,"CCSD"))  wfn = CCSD;
    else if (!strcmp(sw,"G3"))    wfn = G3;
    else throw("Calculation::set() WFN unknown\n");

    if      (!strcmp(sd,"NONE"))        dertype = none;
    else if (!strcmp(sd,"FIRST"))       dertype = first;
    else if (!strcmp(sd,"SECOND"))      dertype = second;
    else if (!strcmp(sd,"UNSPECIFIED")) dertype = unspecified;
    else throw("Calculation::set() DERTYPE unknown\n");
  }
  // this constructor is used if calculation type is known
  Calculation(Jobtype j, Wfn w, Dertype d) {
    jobtype = j; wfn = w; dertype = d;
  }
  // this constructor reads calculation type from input.dat
  Calculation(void) {
    char *sj, *sw, *sd;

    if (!ip_exist(const_cast<char*>("JOBTYPE"),0))
      throw("JOBTYPE keyword not found in input\n");
    ip_string(const_cast<char*>("JOBTYPE"),&sj,0);

    if (!ip_exist(const_cast<char*>("WFN"),0))
      throw("WFN keyword not found in input\n");
    ip_string(const_cast<char*>("WFN"),&sw,0);

    if (ip_exist(const_cast<char*>("DERTYPE"),0))
      ip_string(const_cast<char*>("WFN"),&sd,0);
    else sd = const_cast<char*>("UNSPECIFIED");

    set(sj,sw,sd);
  }
};

enum energy_type { SCF, CISD, CCSD, G3 };

EXT char **atom_basis; // basis set label for each atom

// Useful functions for converting between C and Ruby arrays
//extern VALUE create_array(unsigned int count, int *array);
//extern VALUE create_array(unsigned int count, double *array);
//extern int create_array(VALUE arr, int **array);
//extern int create_array(VALUE arr, double **array);
  
EXT FILE* infile;

/*! Verbosity */
EXT bool verbose;

/*! sanity check boolean */
EXT bool check_only;

/* clean-up */
EXT bool clean_only;

}
  
#endif
