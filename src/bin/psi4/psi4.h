#ifndef __psi_psi4_psi4_h_
#define __psi_psi4_psi4_h_

#include <stdio.h>
#include <string>
#include <ruby.h>             // THIS IS FROM Ruby
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>

#ifdef MAIN
#define EXT
#else
#define EXT extern
#endif

namespace psi {

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
