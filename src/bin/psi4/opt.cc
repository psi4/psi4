/*! \file
    \ingroup OPT09
    \brief Class for simple internal coordinates */

//#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
//#include <cctype>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
//#include <physconst.h>
#include <libpsio/psio.h>
#include <vector>
//#include <libqt/qt.h>
//#include <psifiles.h>

#include <libchkpt/chkpt.h>

extern "C" {
FILE *infile, *outfile;
char *psi_file_prefix;
char *gprgid() { char *prgid = "OPT09"; return(prgid); }
}

namespace psi { namespace opt09 {
void intro(void);
}}

#include "simple.h"
#include "stretch.h"
#include "bend.h"
#include "simple_set.h"

#include "molecular_system.h"

int main(int argc, char **argv) {

  using namespace psi::opt09;
  using namespace std;
  int i, parsed = 1;

  psi_start(&infile,&outfile,&psi_file_prefix,argc-parsed,argv+parsed,0);
  ip_cwk_add(":PSI");
  ip_cwk_add(":OPT09");
  psio_init();
  psio_ipv1_config();

  intro();

  double conv_factor = 1.0; //assume au
  Molecular_system molecules(conv_factor);  // by default, reads geometry from input.dat
  molecules.print();

  double **geom = molecules.get_geom();

  SIMPLE_SET simple_set;
  simple_set.add_simples_from_input();
  simple_set.compute(geom);
  simple_set.compute_s(geom);
  simple_set.print(1);
  simple_set.print_s();

  double * q = simple_set.get_q();
  double ** B = simple_set.get_B(molecules.get_natoms());

  fprintf(outfile,"\nInternal values vector\n");
  for (i=0; i<simple_set.size(); ++i)
    fprintf(outfile,"%15.10lf",q[i]);
  fprintf(outfile,"\n");

  fprintf(outfile,"\nB matrix\n");
  print_mat(B,simple_set.size(),3*molecules.get_natoms(),outfile);

  return 0;

}

namespace psi { namespace opt09 {
void intro(void) {
  fprintf(outfile, "\n\t----------------------------------\n");
  fprintf(outfile, "\t  OPT09: for we're not sure yet   \n");
  fprintf(outfile, "\t       - R. A. King               \n");
  fprintf(outfile, "\t----------------------------------\n");
}

}}

