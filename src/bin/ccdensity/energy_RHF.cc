/*! \file
    \ingroup CCDENSITY
    \brief Calculates the one- and two-electron CC energies using the
    coresponding one- and two-particle density matrices. 
*/
#include <stdio.h>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

    /* ENERGY_RHF(): Compute the RHF CC energy using the one- and two-particle
    ** density matrices.
    **
    */

    void energy_RHF(struct RHO_Params rho_params)
    {
      dpdfile2 D, F;
      dpdbuf4 G, A, B, C, DInts, E, FInts;
      double one_energy=0.0, two_energy=0.0, total_two_energy = 0.0;
      double this_energy;

      fprintf(outfile, "\n\tEnergies re-computed from CC density:\n");
      fprintf(outfile,   "\t-------------------------------------\n");

      dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
      dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fIJ");
      this_energy = 2.0 * dpd_file2_dot(&D, &F);
      dpd_file2_close(&F);
      dpd_file2_close(&D);
      one_energy += this_energy;

      dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
      dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fAB");
      this_energy = 2.0 * dpd_file2_dot(&D, &F);
      dpd_file2_close(&F);
      dpd_file2_close(&D);
      one_energy += this_energy;

      dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
      dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
      this_energy = 2.0 * dpd_file2_dot(&D, &F);
      dpd_file2_close(&F);
      dpd_file2_close(&D);
      one_energy += this_energy;

      dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
      dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
      this_energy = 2.0 * dpd_file2_dot(&D, &F);
      dpd_file2_close(&F);
      dpd_file2_close(&D);
      one_energy += this_energy;

      fprintf(outfile, "\tOne-electron energy        = %20.15f\n", one_energy);
      fflush(outfile);

      total_two_energy = 0.0;


      /* E_ijkl = (2 Gijkl = Gijlk) <ij|kl> */
      dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
      dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gijkl - Gijlk", 2);
      dpd_buf4_sort_axpy(&G, CC_GAMMA, pqsr, 0, 0, "2 Gijkl - Gijlk", -1);
      dpd_buf4_close(&G);
      dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "2 Gijkl - Gijlk");
      dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
      two_energy = dpd_buf4_dot(&A, &G);
      dpd_buf4_close(&A);
      dpd_buf4_close(&G);
      total_two_energy += two_energy;
      fprintf(outfile, "\tIJKL energy                = %20.15f\n", two_energy);
      fflush(outfile);

      two_energy = 0.0;

      /* E_ijka = 2 [ (2 Gijka - Gjika) <ij|ka> ] */
      /* NB: GIjKa is scaled by 1/2 in Gijka.cc. */
      dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
      dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gijka - Gjika", 2);
      dpd_buf4_sort_axpy(&G, CC_GAMMA, qprs, 0, 10, "2 Gijka - Gjika", -1);
      dpd_buf4_close(&G);

      dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "2 Gijka - Gjika");
      dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
      /* The factor of 4 here is necessary because Gijka is multiplied by 1/2 in Gijka.cc */
      two_energy = 4.0 * dpd_buf4_dot(&E, &G);
      dpd_buf4_close(&E);
      dpd_buf4_close(&G);
      total_two_energy += two_energy;
      fprintf(outfile, "\tIJKA energy                = %20.15f\n", two_energy);
      fflush(outfile);

      two_energy = 0.0;

      /* Generate spin-adapted Gijab jut for energy calculation */
      dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
      dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gijab - Gijba", 2);
      dpd_buf4_sort_axpy(&G, CC_GAMMA, pqsr, 0, 5, "2 Gijab - Gijba", -1);
      dpd_buf4_close(&G);
      dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "2 Gijab - Gijba");
      dpd_buf4_init(&DInts, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      two_energy = 2.0 * dpd_buf4_dot(&G, &DInts);
      dpd_buf4_close(&G);
      dpd_buf4_close(&DInts);

      total_two_energy += two_energy;
      fprintf(outfile, "\tIJAB energy                = %20.15f\n", two_energy);
      fflush(outfile);

      two_energy = 0.0;

      dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
      dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
      two_energy += 2.0 * dpd_buf4_dot(&G, &C);
      dpd_buf4_close(&G);

      dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
      dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
      two_energy += 2.0 * dpd_buf4_dot(&G, &C);
      dpd_buf4_close(&G);

      dpd_buf4_init(&DInts, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
      dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
      two_energy -= 2.0 * dpd_buf4_dot(&G, &DInts);
      dpd_buf4_close(&G);
      dpd_buf4_close(&DInts);

      total_two_energy += two_energy;
      fprintf(outfile, "\tIBJA energy                = %20.15f\n", two_energy);
      fflush(outfile);

      two_energy = 0.0;

      /* Generate spin-adapted Gciab just for energy calculation */
      dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
      dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gciab - Gciba", 2);
      dpd_buf4_sort_axpy(&G, CC_GAMMA, pqsr, 11, 5, "2 Gciab - Gciba", -1);
      dpd_buf4_close(&G);

      dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "2 Gciab - Gciba");
      dpd_buf4_init(&FInts, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
      dpd_buf4_sort(&FInts, CC_FINTS, qpsr, 11, 5, "F <ai|bc>");
      dpd_buf4_close(&FInts);
      dpd_buf4_init(&FInts, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
      /* The factor of 4 here is necessary because Gciab is multiplied by 1/2 in Gciab.cc */
      two_energy = 4*dpd_buf4_dot(&FInts, &G);
      dpd_buf4_close(&FInts);
      dpd_buf4_close(&G);

      total_two_energy += two_energy;
      fprintf(outfile, "\tCIAB energy                = %20.15f\n", two_energy);
      fflush(outfile);

      two_energy = 0.0;
      dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");

      dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gabcd - Gabdc", 2);
      dpd_buf4_sort_axpy(&G, CC_GAMMA, pqsr, 5, 5, "2 Gabcd - Gabdc", -1);
      dpd_buf4_close(&G);

      dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "2 Gabcd - Gabdc");
      dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
      two_energy = dpd_buf4_dot(&B, &G);
      dpd_buf4_close(&B);
      dpd_buf4_close(&G);

      total_two_energy += two_energy;
      fprintf(outfile, "\tABCD energy                = %20.15f\n", two_energy);

      fprintf(outfile, "\tTotal two-electron energy  = %20.15f\n", total_two_energy);
      if (params.ground) {
	fprintf(outfile, "\tCCSD correlation energy    = %20.15f\n",
		one_energy + total_two_energy);
	fprintf(outfile, "\tTotal CCSD energy          = %20.15f\n",
		one_energy + total_two_energy + moinfo.eref);
      }
      else {
	fprintf(outfile, "\tTotal EOM CCSD correlation energy        = %20.15f\n",
		one_energy + total_two_energy);
	fprintf(outfile, "\tCCSD correlation + EOM excitation energy = %20.15f\n",
		moinfo.ecc + params.cceom_energy);
	fprintf(outfile, "\tTotal EOM CCSD energy                    = %20.15f\n",
		one_energy + total_two_energy + moinfo.eref);
      }
    }

  }} // namespace psi::ccdensity
