/*! \defgroup DF-MP2 df-mp2: Density Fitted MP2 */

/*! \file
 *  \ingroup DF-MP2
 *  \brief Density-fitted MP2
 *
 *  DF-MP2
 *
 *  Density-Fitted MP2 Program (no symmetry)
 *
 *  Justin Turney, University of Georgia
 *  C. David Sherrill, Georgia Institute of Technology
 *  Andy Simmonett, University of Georgia
 *  Edward Hohenstein, Georgia Institute of Technology
 *
 *  Created 7/15/09
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#define DEBUG 0
#define TIME_DF_MP2 1

//#include <omp.h>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>

#include <iostream>

#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/factory.h>
#include <libmints/symmetry.h>
#include <libmints/wavefunction.h>
// #include <libmints/mints.h>


namespace psi{ namespace dfmp2{
  void title(void);
}}

FILE *infile = NULL, *outfile = NULL;
char *psi_file_prefix = NULL;

extern "C" {
  char *gprgid()
    {
      const char *prgid = "DF-MP2";
      return const_cast<char*>(prgid);
    }
}

std::string to_string(const int val);  // in libmints/matrix.cc

using namespace psi;
using namespace psi::dfmp2;

int main(int argc, char * argv[]) {

}

namespace psi{ namespace dfmp2{

void title(void)
{
  fprintf(outfile, "\t\t\t*************************\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*         DF-MP2        *\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*************************\n");
  fflush(outfile);
}

}} // end namespaces

