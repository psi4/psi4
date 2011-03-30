#ifndef __psi_psi4_psi4_h_
#define __psi_psi4_psi4_h_

#include <stdio.h>
#include <string>

#ifdef MAIN
#define EXT
#else
#define EXT extern
#endif

namespace psi {

EXT FILE* infile;
EXT std::string outfile_name;

/*! Directory location of the input file */
EXT std::string infile_directory;

/*! Verbosity */
EXT bool verbose;

/*! sanity check boolean */
EXT bool check_only;

/*! Leave psi temp file */
EXT bool messy;

/*! clean-up */
EXT bool clean_only;

}

#endif
