/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "Params.h"
#define EXTERN
#include "globals.h"
#include "MOInfo.h"

namespace psi { namespace ccenergy {

double rhf_mp2_energy(void);
double uhf_mp2_energy(void);

double mp2_energy(void)
{
	/* Note that if we reach this point and ref=1 (ROHF), then we aren't using 
	 * semicanonical orbitals and so we can't compute a non-iterative MBPT(2) 
	 * energy */
  if(params.ref == 0) return(rhf_mp2_energy());
  else if(params.ref == 2) return(uhf_mp2_energy());
  else return 0.0;
  
}

double rhf_mp2_energy(void)
{
  double T2_energy, T1_energy;
  dpdfile2 F, T1, D1;
  dpdbuf4 T2, D;
  dpdbuf4 S;
  double os_energy, ss_energy, scs_energy;
  
  /* Initialize MP2 T1 Amps (zero for true HF references) */
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_copy(&F, CC_OEI, "tIA (MP2)");
  dpd_file2_close(&F);
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA (MP2)");
  dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "dIA");
  dpd_file2_dirprd(&D1, &T1);
  dpd_file2_close(&D1);  
  dpd_file2_close(&T1);

  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA (MP2)");
  T1_energy = 2.0 * dpd_file2_dot(&F, &T1);
  dpd_file2_close(&F);
  dpd_file2_close(&T1);

  /* Initialize MP2 T2 Amps */
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_copy(&D, CC_TAMPS, "tIjAb (MP2)");
  dpd_buf4_close(&D);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb (MP2)");
  dpd_buf4_init(&D, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
  dpd_buf4_dirprd(&D, &T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&T2);
  
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb (MP2)");
  T2_energy = dpd_buf4_dot(&D, &T2);

  dpd_buf4_init(&S, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  os_energy = dpd_buf4_dot(&S, &T2);
  dpd_buf4_close(&S);
  ss_energy = (T2_energy - os_energy);

  moinfo.emp2_ss = ss_energy;
  moinfo.emp2_os = os_energy;

  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  
  return (T2_energy+T1_energy);
}

double uhf_mp2_energy(void)
{
  double E2AA, E2BB, E2AB, T1A, T1B;
  dpdbuf4 T2, D;
  dpdfile2 T1, F, D1;
  
  /* Initialize MP2 T1 Amps (zero for true HF references) */
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_copy(&F, CC_OEI, "tIA (MP2)");
  dpd_file2_close(&F);
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA (MP2)");
  dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "dIA");
  dpd_file2_dirprd(&D1, &T1);
  dpd_file2_close(&D1);  
  dpd_file2_close(&T1);
  
  dpd_file2_init(&F, CC_OEI, 0, 2, 3, "fia");
  dpd_file2_copy(&F, CC_OEI, "tia (MP2)");
  dpd_file2_close(&F);
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia (MP2)");
  dpd_file2_init(&D1, CC_OEI, 0, 2, 3, "dia");
  dpd_file2_dirprd(&D1, &T1);
  dpd_file2_close(&D1);  
  dpd_file2_close(&T1);
  
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA (MP2)");
  T1A = dpd_file2_dot(&F, &T1);
  dpd_file2_close(&F);
  dpd_file2_close(&T1);

  dpd_file2_init(&F, CC_OEI, 0, 2, 3, "fia");
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia (MP2)");
  T1B = dpd_file2_dot(&F, &T1);
  dpd_file2_close(&F);
  dpd_file2_close(&T1);

  dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
  dpd_buf4_copy(&D, CC_TAMPS, "tIJAB (MP2)");
  dpd_buf4_close(&D);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB (MP2)");
  dpd_buf4_init(&D, CC_DENOM, 0, 2, 7, 2, 7, 0, "dIJAB");
  dpd_buf4_dirprd(&D, &T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&T2);

  dpd_buf4_init(&D, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
  dpd_buf4_copy(&D, CC_TAMPS, "tijab (MP2)");
  dpd_buf4_close(&D);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab (MP2)");
  dpd_buf4_init(&D, CC_DENOM, 0, 12, 17, 12, 17, 0, "dijab");
  dpd_buf4_dirprd(&D, &T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&T2);

  dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  dpd_buf4_copy(&D, CC_TAMPS, "tIjAb (MP2)");
  dpd_buf4_close(&D);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb (MP2)");
  dpd_buf4_init(&D, CC_DENOM, 0, 22, 28, 22, 28, 0, "dIjAb");
  dpd_buf4_dirprd(&D, &T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&T2);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB (MP2)");
  dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
  E2AA = dpd_buf4_dot(&D, &T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&T2);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab (MP2)");
  dpd_buf4_init(&D, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
  E2BB = dpd_buf4_dot(&D, &T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&T2);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb (MP2)");
  dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  E2AB = dpd_buf4_dot(&D, &T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&T2);

  // We define EMP2_SS as the same-spin pair energy, and EMP2_OS as the 
  // opposite-spin pair energy (singles not included)
  moinfo.emp2_ss = E2AA + E2BB;
  moinfo.emp2_os = E2AB;

  return(T1A + T1B + E2AA + E2BB + E2AB);
}
}} // namespace psi::ccenergy
