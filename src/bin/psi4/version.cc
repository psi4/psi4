#include <libparallel/parallel.h>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <psiconfig.h>

namespace psi {

/*! Print PSI version information that was set in configure.ac */
void print_version(FILE *myout)
{
  fprintf(myout, "    -----------------------------------------------------------------------\n");
  fprintf(myout, "          PSI4: An Open-Source Ab Initio Electronic Structure Package\n");
  fprintf(myout, "                              PSI %s Driver\n", PSI_VERSION);

  // Are we using git? If so,what version string
  bool using_git = false;
  char git_version[2048];

  std::stringstream HEAD;
  HEAD << PSI_TOP_SRCDIR;
  HEAD << "/.git/HEAD";
  FILE* fh = fopen(HEAD.str().c_str(), "r");
  if (fh) {
    char line[2048];
    if (fgets(line, 2048, fh)) {
      char* HEAD_REF = line + 1; // First space 
      HEAD_REF = (char*) memchr(HEAD_REF, ' ', strlen(line)); // Second space
      if (HEAD_REF) {
        HEAD_REF++; // Eliminate the space
        HEAD_REF[strlen(HEAD_REF) - 1] = '\0'; // Eliminate the trailing newline
        std::stringstream VERSION;
        VERSION << PSI_TOP_SRCDIR;
        VERSION << "/.git/";
        VERSION << HEAD_REF;
        FILE* fh2 = fopen(VERSION.str().c_str(), "r");
        if (fh2) {
          if (fgets(git_version, 2048, fh2)) {
            using_git = true;
          } 
          fclose(fh2);
        }   
      }    
    }
    fclose(fh);
  }

  if (using_git) {
    fprintf(myout, "            Using Git: Rev %s", git_version);
  } else {
    fprintf(myout, "                                  Not using Git\n");
  }

  fprintf(myout, "\n");
  fprintf(myout, "    J. M. Turney, A. C. Simmonett, R. M. Parrish, E. G. Hohenstein,\n");
  fprintf(myout, "    F. Evangelista, J. T. Fermann, B. J. Mintz, L. A. Burns, J. J. Wilke,\n");
  fprintf(myout, "    M. L. Abrams, N. J. Russ, M. L. Leininger, C. L. Janssen, E. T. Seidl,\n");
  fprintf(myout, "    W. D. Allen, H. F. Schaefer, R. A. King, E. F. Valeev, C. D. Sherrill,\n");
  fprintf(myout, "    and T. D. Crawford, WIREs Comput. Mol. Sci., (2011) (doi: 10.1002/wcms.93)\n");

  fprintf(myout, "\n");
  fprintf(myout, "                         Additional Contributions by\n");
  fprintf(myout, "    A. E. DePrince, M. Saitow, U. Bozkaya, A. Yu. Sokolov\n");
  fprintf(myout, "    -----------------------------------------------------------------------\n\n");
  Communicator::world->print(myout);
}

}
