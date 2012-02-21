#include <libparallel/parallel.h>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <psiconfig.h>

#include "gitversion.h"

namespace psi {

/*! Print PSI version information that was set in configure.ac */
void print_version(FILE *myout)
{
  fprintf(myout, "    -----------------------------------------------------------------------\n");
  fprintf(myout, "          PSI4: An Open-Source Ab Initio Electronic Structure Package\n");
  fprintf(myout, "                              PSI %s Driver\n", PSI_VERSION);

  // Are we using git? If so,what version string
#ifdef GIT_VERSION
  fprintf(myout, "\n               Git: Rev " GIT_VERSION "\n");
#endif

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
