/*! \file
    \ingroup LMP2
    \brief print the domains
*/
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
//#include <cmath>
//#include <libipv1/ip_lib.h>
//#include <libciomr/libciomr.h>
//#include <libchkpt/chkpt.h>
//#include <libpsio/psio.h>
//#include <libqt/qt.h>
//#include <psifiles.h>
#define EXTERN
#include <libmints/mints.h>
#include "globals.h"

namespace psi{ namespace lmp2{

void LMP2::print_domains(double *s){

  int i, j, cnt, max;
  double domain_tot, domain_ave;

  max = 0;
  for(i=0; i < nocc; i++)
    if(domain_len[i] > max) max = domain_len[i];

//  if(myid == 0) {
    fprintf(outfile, "   Orbital  Domain");
    for(i=0; i < max-2; i++) fprintf(outfile, "   "); /* formatting junk */
    fprintf(outfile, "  Completeness\n");
    fprintf(outfile, "   -------  ------");
    for(i=0; i < max-2; i++) fprintf(outfile, "---"); /* more formatting junk */
    fprintf(outfile, "  ------------\n");
    for(i=0; i < nocc; i++) {
      fprintf(outfile, "      %2d    ",i);
      for(j=0,cnt=0; j < natom; j++) if(domain[i][j]) { fprintf(outfile, " %2d", j); cnt++; }
      if(cnt < max) for(; cnt < max; cnt++) fprintf(outfile, "   ");
      fprintf(outfile, "     %7.5f\n", s[i]);
    }
    domain_tot = 0;
    for(i=0; i < nocc; i++)
      domain_tot += domain_len[i];
    domain_ave = domain_tot/nocc;
    fprintf(outfile, "\n   The average domain length is %4.2lf\n", domain_ave);
    fflush(outfile);
//  }
}

}} // namespace psi::lmp2

