#ifndef __psi_psi4_psi4_h_
#define __psi_psi4_psi4_h_

#include <stdio.h>
#include <string>
//#include <ruby.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
//#include "task.h"

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

/*! The name of the input file, could be useful for something */
EXT std::string g_szInputFile;

/*! Name of the output file */
EXT std::string g_szOutputFilename;

/*! Verbosity */
EXT bool verbose;

/*! sanity check boolean */
EXT bool check_only;

/* clean-up */
EXT bool clean_only;

/*! Global task object */
//EXT Task *g_cTask;
//EXT VALUE g_rbTask;

# ifdef MAIN
  /*! All classes that provide a Ruby interface need this to say which module they belong to */
//  EXT VALUE g_rbPsi = Qnil;
  /*! How much to indent the Ruby puts output by. */
  EXT int g_iPutsIndent = 0;
  EXT bool g_bQuietRuby = false;
  EXT bool g_bIRB = false;
  EXT void *g_rbExecNode = NULL;  // Only used in Ruby 1.9+
# else
//  extern VALUE g_rbPsi;
  extern int g_iPutsIndent;
  extern bool g_bQuietRuby;
  extern bool g_bIRB;
  extern void *g_rbExecNode;
# endif
  
// Ruby Pre-1.9 string access
#if !defined(RSTRING_LEN)
#  define RSTRING_LEN(x) (RSTRING(x)->len)
#endif
#if !defined(RSTRING_PTR)
#  define RSTRING_PTR(x) (RSTRING(x)->ptr)
#endif
  
// Helper macros
/*! Cast a C/++ function to a Ruby acceptable form */
#define RUBYCAST(x) (VALUE (*)(...))(x)

/*! Help the user get the Psi object from Ruby. */
#define RUBYPSIDATA(RObj, OType, OPtr) \
	Data_Get_Object(RObj, OType, OPtr);
}
#endif
