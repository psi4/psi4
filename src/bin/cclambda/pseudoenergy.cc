/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

double pseudoenergy(struct L_Params L_params)
{
  double LIJAB_energy, Lijab_energy, LIjAb_energy;
  double LIA_energy=0.0, Lia_energy=0.0, tval;
  dpdbuf4 LIJAB, Lijab, LIjAb, D;
  dpdfile2 Lia, LIA, Fme, FME;
  int L_irr;
  L_irr = L_params.irrep;

  if ( L_params.ground || ((L_params.irrep ==0)&&(fabs(L_params.R0)>1e-10)) ) {
    if(params.ref == 0) { /** RHF **/

      Lia_energy = 0.0;
      dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
      dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
      LIA_energy = dpd_file2_dot(&FME,&LIA);
      dpd_file2_close(&LIA);
      dpd_file2_close(&FME);

      LIJAB_energy = 0.0;
      Lijab_energy = 0.0;
      dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
      dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      LIjAb_energy = dpd_buf4_dot(&D, &LIjAb);
      dpd_buf4_close(&LIjAb);
      dpd_buf4_close(&D);
    }
    else if(params.ref == 1) { /** ROHF **/
      dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
      dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
      dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 0, 1, "Lia");
      dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
  
      LIA_energy = dpd_file2_dot(&FME,&LIA);
      Lia_energy = dpd_file2_dot(&Fme,&Lia);
  
      dpd_file2_close(&Lia);
      dpd_file2_close(&LIA);
      dpd_file2_close(&Fme);
      dpd_file2_close(&FME);
  
      dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
      dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      LIJAB_energy = dpd_buf4_dot(&D, &LIJAB);
      dpd_buf4_close(&LIJAB);
      dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
      Lijab_energy = dpd_buf4_dot(&D, &Lijab);
      dpd_buf4_close(&Lijab);
      dpd_buf4_close(&D);
  
      dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      LIjAb_energy = dpd_buf4_dot(&D, &LIjAb);
      dpd_buf4_close(&LIjAb);
      dpd_buf4_close(&D);
    }
    else if(params.ref == 2) { /** UHF **/
  
      dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
      dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
      dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 2, 3, "Lia");
      dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
  
      LIA_energy = dpd_file2_dot(&FME,&LIA);
      Lia_energy = dpd_file2_dot(&Fme,&Lia);
  
      dpd_file2_close(&Lia);
      dpd_file2_close(&LIA);
      dpd_file2_close(&Fme);
      dpd_file2_close(&FME);
  
      dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
      dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      LIJAB_energy = dpd_buf4_dot(&D, &LIJAB);
      dpd_buf4_close(&LIJAB);
      dpd_buf4_close(&D);
  
      dpd_buf4_init(&D, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
      dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
      Lijab_energy = dpd_buf4_dot(&D, &Lijab);
      dpd_buf4_close(&Lijab);
      dpd_buf4_close(&D);
  
      dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
      dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
      LIjAb_energy = dpd_buf4_dot(&D, &LIjAb);
      dpd_buf4_close(&LIjAb);
      dpd_buf4_close(&D);
    }
    /*
      fprintf(outfile, "One A Energy = %20.14f\n", LIA_energy);
      fprintf(outfile, "One B Energy = %20.14f\n", Lia_energy);
      fprintf(outfile, "Two AA Energy = %20.14f\n", LIJAB_energy);
      fprintf(outfile, "Two BB Energy = %20.14f\n", Lijab_energy);
      fprintf(outfile, "Two AB Energy = %20.14f\n", LIjAb_energy);
    */
    return (LIJAB_energy + Lijab_energy + LIjAb_energy);
  }
  else { /* since pseudoenergy is 0 lets compute norm instead */
    if (params.ref <= 1) { /* RHF or ROHF */
      dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 0, 1, "Lia");
      dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
      LIA_energy = dpd_file2_dot_self(&LIA);
      Lia_energy = dpd_file2_dot_self(&Lia);
      dpd_file2_close(&Lia);
      dpd_file2_close(&LIA);
      dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      LIJAB_energy = dpd_buf4_dot_self(&LIJAB);
      dpd_buf4_close(&LIJAB);
      dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
      Lijab_energy = dpd_buf4_dot_self(&Lijab);
      dpd_buf4_close(&Lijab);
      dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      LIjAb_energy = dpd_buf4_dot_self(&LIjAb);
      dpd_buf4_close(&LIjAb);
      tval = LIA_energy + Lia_energy + LIJAB_energy + Lijab_energy + LIjAb_energy;
      tval = sqrt(tval);
      return tval;
    }
    else if (params.ref == 2) { /* UHF */
      dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
      dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 2, 3, "Lia");
      LIA_energy = dpd_file2_dot_self(&LIA);
      Lia_energy = dpd_file2_dot_self(&Lia);
      dpd_file2_close(&Lia);
      dpd_file2_close(&LIA);
      dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      LIJAB_energy = dpd_buf4_dot_self(&LIJAB);
      dpd_buf4_close(&LIJAB);
      dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
      Lijab_energy = dpd_buf4_dot_self(&Lijab);
      dpd_buf4_close(&Lijab);
      dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
      LIjAb_energy = dpd_buf4_dot_self(&LIjAb);
      dpd_buf4_close(&LIjAb);
      tval = LIA_energy + Lia_energy + LIJAB_energy + Lijab_energy + LIjAb_energy;
      tval = sqrt(tval);
      return tval;
    }
  }
}

}} // namespace psi::cclambda
