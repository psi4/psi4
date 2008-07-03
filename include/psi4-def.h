#ifndef psi_include_psi4_def_h
#define psi_include_psi4_def_h

// this header is included by every stand along program
#include <psi4-dec.h>

extern "C" {
  FILE *infile, *outfile;
  char *psi_file_prefix;
}

namespace psi {
  Module module;
}

#endif

