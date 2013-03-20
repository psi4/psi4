#ifndef psi_include_psi4_def_h
#define psi_include_psi4_def_h

// Define the psi namespace global variables here (one time); other
// references use extern's (in psi4-dec.h) 
#include <psi4-dec.h>

namespace psi {
  FILE *outfile;
  char *psi_file_prefix;
  std::string outfile_name;
}

#endif

