/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {
                                                                                                    
void G_norm(void) {
  dpdfile2 G1;
  dpdbuf4 G;
  double value, value1, dot_IA, dot_ia, dot_AI, dot_ai;
  int G_irr = 0;

  fprintf(outfile,"Calculating overlaps of CC_OEI\n");
  dpd_file2_init(&G1, CC_OEI, G_irr, 0, 1, "DIA");
  dot_IA = dpd_file2_dot_self(&G1);
  dpd_file2_close(&G1);
  dpd_file2_init(&G1, CC_OEI, G_irr, 0, 1, "Dia");
  dot_ia = dpd_file2_dot_self(&G1);
  dpd_file2_close(&G1);
  dpd_file2_init(&G1, CC_OEI, G_irr, 0, 1, "DAI");
  dot_AI = dpd_file2_dot_self(&G1);
  dpd_file2_close(&G1);
  dpd_file2_init(&G1, CC_OEI, G_irr, 0, 1, "Dai");
  dot_ai = dpd_file2_dot_self(&G1);
  dpd_file2_close(&G1);
  /*
  fprintf(outfile,"<DIA|DIA> = %15.10lf\n", dot_IA);
  fprintf(outfile,"<Dia|Dia> = %15.10lf\n", dot_ia);
  fprintf(outfile,"<DAI|DAI> = %15.10lf\n", dot_AI);
  fprintf(outfile,"<Dai|Dai> = %15.10lf\n", dot_ai);
  */
  fprintf(outfile,"\t<Dpq|Dqp>     = %15.10lf\n", dot_IA+dot_ia+dot_AI+dot_ai);

  fprintf(outfile,"Calculating overlaps of CC_GAMMA\n");

  dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "GIJKL");
  value = dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "Gijkl");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 0, 0, 0, 0, 0, "GIjKl");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  fprintf(outfile,"\t<Gijkl|Gijkl> = %15.10lf\n", value);

  dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 10, 2, 10, 0, "GIJKA");
  value = dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 10, 2, 10, 0, "Gijka");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  fprintf(outfile,"\t<Gijka|Gijka> = %15.10lf\n",value);

  dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  value = dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "Gijab");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  fprintf(outfile,"\t<Gijab|Gijab> = %15.10lf\n", value);


  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIBJA");
  value = dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "Gibja");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbJa");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBjA");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbjA");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBJa");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  fprintf(outfile,"\t<Gibja|Gibja> = %15.10lf\n",value);

  dpd_buf4_init(&G, CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  value = dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "Gciab");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  fprintf(outfile,"\t<Gciab|Gciab> = %15.10lf\n",value);

  dpd_buf4_init(&G, CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
  value = dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "Gabcd");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
  value += dpd_buf4_dot_self(&G);
  dpd_buf4_close(&G);
  fprintf(outfile,"\t<Gabcd|Gabcd> = %15.10lf\n", value);

  return;
}

}} // namespace psi::ccdensity
