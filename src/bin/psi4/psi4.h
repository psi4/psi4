#ifndef __psi_psi4_psi4_h_
#define __psi_psi4_psi4_h_

#include <stdio.h>
#include <string>
#include <ruby.h>             // THIS IS FROM Ruby
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>

#ifdef MAIN
#define EXT
#else
#define EXT extern
#endif

namespace psi {

// Define Calculation class

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

    if (!ip_exist("JOBTYPE",0))
      throw("JOBTYPE keyword not found in input\n");
    ip_string("JOBTYPE",&sj,0);

    if (!ip_exist("WFN",0))
      throw("WFN keyword not found in input\n");
    ip_string("WFN",&sw,0);

    if (ip_exist("DERTYPE",0))
      ip_string("WFN",&sd,0);
    else sd = "UNSPECIFIED";

    set(sj,sw,sd);
  }
};

enum energy_type { SCF, CISD, CCSD, G3 };

char **atom_basis; // basis set label for each atom



  // Useful functions for converting between C and Ruby arrays
  extern VALUE create_array(unsigned int count, int *array);
  extern VALUE create_array(unsigned int count, double *array);
  extern int create_array(VALUE arr, int **array);
  extern int create_array(VALUE arr, double **array);
  
  /*! Namespace for containing global variables.
      Note: Global does not necessary infer global to all PSI */
  namespace Globals {
    /*! All output should be sent to this variable */
    EXT FILE* g_fOutput;
    
    /*! The name of the input file, could be useful for something */
    EXT std::string g_szInputFile;
    
    /*! Name of the output file */
    EXT std::string g_szOutputFilename;
    
    /*! Verbosity */
    EXT bool g_bVerbose;
    
    #ifdef MAIN
    /*! All classes that provide a Ruby interface need this to say which module they belong to */
    EXT VALUE g_rbPsi = Qnil;
    /*! How much to indent the Ruby puts output by. */
    EXT int g_iPutsIndent = 0;
    EXT bool g_bQuietRuby = false;
    EXT bool g_bIRB = false;
    #else
    EXT VALUE g_rbPsi;
    EXT int g_iPutsIndent;
    EXT bool g_bQuietRuby;
    EXT bool g_bIRB;
    #endif
  };
  
  // Helper macros
  /*! Cast a C/++ function to a Ruby acceptable form */
  #define RUBYCAST(x) (VALUE (*)(...))(x)

  /*! Help the user get the Psi object from Ruby. */
  #define RUBYPSIDATA(RObj, OType, OPtr) \
  	Data_Get_Object(RObj, OType, OPtr);
}
#endif
