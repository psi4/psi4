/*! \file 
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <libiwl/iwl.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

    /* DEANTI_RHF(): Convert the RHF two-particle density from an
    ** energy expression using antisymmetrized Dirac integrals to one
    ** using simple Diract integrals. The original, Fock-adjusted
    ** density (see the comments in fock.c) corresponds to a
    ** two-electron energy (or energy derivative) expression of the
    ** form:
    **
    ** E(TWO) = 1/4 sum_pqrs Gpqrs <pq||rs>
    **
    ** However, the derivative two-electron integrals are produced in
    ** Mulliken-ordered, symmetric form rather than Dirac-ordered
    ** antisymmetric form.  This code alters the two-particle density
    ** matrix ordering for the energy expression of the form:
    **
    ** E(TWO) = 1/2 sum_pqrs Gpqrs <pq|rs>
    **
    ** The final conversion to Mulliken ordering is taken care of in
    ** dump.c
    **
    ** The second equation above may be derived via
    **
    ** E(TWO) = 1/4 sum_pqrs Gpqrs (<pq|rs> - <pq|sr>)
    **        = 1/4 sum_pqrs Gpqrs <pq|rs> - 1/4 sum_pqrs Gpqrs <pq|sr>
    **        = 1/4 sum_pqrs Gpqrs <pq|rs> - 1/4 sum_pqrs Gpqsr <pq|rs>
    **        = 1/4 sum_pqrs (Gpqrs - Gpqsr) <pq|rs>
    **        = 1/2 sum_pqrs Gpqrs <pq|rs>
    **
    ** Equations for the individual components are given explicitly in
    ** comments below.
    ** */

    void deanti_RHF(struct RHO_Params rho_params)
    {
      dpdbuf4 G1, G2;
      dpdfile2 D, F;
      double one_energy = 0.0, two_energy = 0.0, total_two_energy = 0.0;
      dpdbuf4 A, B, C, DInts, E, FInts;

      if(!params.aobasis) {
	fprintf(outfile, "\n\tEnergies re-computed from Mulliken density:\n");
	fprintf(outfile,   "\t-------------------------------------------\n");

	dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	dpd_file2_init(&F, CC_OEI, 0, 0, 0, "h(i,j)");
	one_energy += 2.0*dpd_file2_dot(&D, &F);
	dpd_file2_close(&F);
	dpd_file2_close(&D);

	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	dpd_file2_init(&F, CC_OEI, 0, 1, 1, "h(a,b)");
	one_energy += 2.0*dpd_file2_dot(&D, &F);
	dpd_file2_close(&F);
	dpd_file2_close(&D);

	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
	one_energy += 2.0*dpd_file2_dot(&D, &F);
	dpd_file2_close(&F);
	dpd_file2_close(&D);

	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
	one_energy += 2.0*dpd_file2_dot(&D, &F);
	dpd_file2_close(&F);
	dpd_file2_close(&D);

	fprintf(outfile, "\tOne-electron energy        = %20.15f\n", one_energy);
	fflush(outfile);
      }

      /* E_ijkl = (2 Gijkl - Gijlk) <ij|kl> */
      dpd_buf4_init(&G1, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
      dpd_buf4_scmcopy(&G1, CC_GAMMA, "2 Gijkl - Gijlk", 2);
      dpd_buf4_sort_axpy(&G1, CC_GAMMA, pqsr, 0, 0, "2 Gijkl - Gijlk", -1);
      dpd_buf4_close(&G1);

      dpd_buf4_init(&G1, CC_GAMMA, 0, 0, 0, 0, 0, 0, "2 Gijkl - Gijlk");
      dpd_buf4_copy(&G1, CC_GAMMA, "GIjKl");
      if(!params.aobasis) {
	dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	two_energy = dpd_buf4_dot(&A, &G1);
	dpd_buf4_close(&A);
	fprintf(outfile, "\tIJKL energy                = %20.15f\n", two_energy);
	total_two_energy += two_energy;
      }
      dpd_buf4_close(&G1);

      /* E_ijka = 2 [ (2 Gijka - Gjika) <ij|ka> ] */
      /* NB: GIjKa is scaled by 1/2 in Gijka.cc. */
      dpd_buf4_init(&G1, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
      dpd_buf4_scmcopy(&G1, CC_GAMMA, "2 Gijka - Gjika", 2);
      dpd_buf4_sort_axpy(&G1, CC_GAMMA, qprs, 0, 10, "2 Gijka - Gjika", -1);
      dpd_buf4_close(&G1);

      dpd_buf4_init(&G1, CC_GAMMA, 0, 0, 10, 0, 10, 0, "2 Gijka - Gjika");
      /* The factor of 2 here is necessary because Gijka is multiplied by 1/2 in Gijka.cc */
      dpd_buf4_scmcopy(&G1, CC_GAMMA, "GIjKa", 2);
      if(!params.aobasis) {
	dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	two_energy = dpd_buf4_dot(&E, &G1);
	dpd_buf4_close(&E);
	/* The factor of 4 here is necessary because Gijka is multiplied by 1/2 in Gijka.cc */
	fprintf(outfile, "\tIJKA energy                = %20.15f\n", 4*two_energy);
	total_two_energy += 4*two_energy;
      }
      dpd_buf4_close(&G1);

      /* E_ijab = [ (2 Gijab - Gijba) - GIBJA - GIbjA ] <ij|ab> */
      dpd_buf4_init(&G1, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
      dpd_buf4_scmcopy(&G1, CC_GAMMA, "2 Gijab - Gijba", 2);
      dpd_buf4_sort_axpy(&G1, CC_GAMMA, pqsr, 0, 5, "2 Gijab - Gijba", -1);
      dpd_buf4_close(&G1);

      dpd_buf4_init(&G2, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
      dpd_buf4_sort_axpy(&G2, CC_GAMMA, prsq, 0, 5, "2 Gijab - Gijba", -1);
      dpd_buf4_close(&G2);
      dpd_buf4_init(&G2, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
      dpd_buf4_sort_axpy(&G2, CC_GAMMA, prsq, 0, 5, "2 Gijab - Gijba", -1);
      dpd_buf4_close(&G2);  

      dpd_buf4_init(&G1, CC_GAMMA, 0, 0, 5, 0, 5, 0, "2 Gijab - Gijba");
      dpd_buf4_scmcopy(&G1, CC_GAMMA, "GIjAb", 2);
      if(!params.aobasis) {
	dpd_buf4_init(&DInts, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	two_energy = 2.0 * dpd_buf4_dot(&DInts, &G1);
	dpd_buf4_close(&DInts);
	fprintf(outfile, "\tIJAB energy                = %20.15f\n", two_energy);
	total_two_energy += two_energy;
      }

      dpd_buf4_close(&G1);
  
      /* G'_IbJa = 2 GIBJA + 2 GIbJa */
      dpd_buf4_init(&G1, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
      dpd_buf4_init(&G2, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
      dpd_buf4_axpy(&G2, &G1, 1);
      dpd_buf4_close(&G2);
      //  dpd_buf4_scm(&G1, 2);
      if(!params.aobasis) {
	dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	two_energy = 2.0 * dpd_buf4_dot(&C, &G1);
	dpd_buf4_close(&C);
	fprintf(outfile, "\tIBJA energy                = %20.15f\n", two_energy);
	total_two_energy += two_energy;
      }

      dpd_buf4_close(&G1);

      /* E_ciab = 2 [ (2 Gciab - Gciba) <ci|ab> ] */
      /* NB: Gciab is scaled by 1/2 in Gciab.cc. */
      dpd_buf4_init(&G1, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
      dpd_buf4_scmcopy(&G1, CC_GAMMA, "2 Gciab - Gciba", 2);
      dpd_buf4_sort_axpy(&G1, CC_GAMMA, pqsr, 11, 5, "2 Gciab - Gciba", -1);
      dpd_buf4_close(&G1);

      dpd_buf4_init(&G1, CC_GAMMA, 0, 11, 5, 11, 5, 0, "2 Gciab - Gciba");
      /* The factor of 2 here is necessary because Gciab is multiplied by 1/2 in Gciab.cc */
      dpd_buf4_scmcopy(&G1, CC_GAMMA, "GCiAb", 2);
      if(!params.aobasis) {
	dpd_buf4_init(&FInts, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_buf4_sort(&FInts, CC_FINTS, qpsr, 11, 5, "F <ai|bc>");
	dpd_buf4_close(&FInts);
	dpd_buf4_init(&FInts, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
	two_energy = dpd_buf4_dot(&FInts, &G1);
	dpd_buf4_close(&FInts);
	/* The factor of 4 here is necessary because Gciab is multiplied by 1/2 in Gciab.cc */
	fprintf(outfile, "\tCIAB energy                = %20.15f\n", 4*two_energy);
	total_two_energy += 4*two_energy;
      }

      dpd_buf4_close(&G1);

      /* E_abcd = (2 Gabcd - Gabdc) <ab|cd> */
      dpd_buf4_init(&G1, CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");

      dpd_buf4_scmcopy(&G1, CC_GAMMA, "2 Gabcd - Gabdc", 2);
      dpd_buf4_sort_axpy(&G1, CC_GAMMA, pqsr, 5, 5, "2 Gabcd - Gabdc", -1);
      dpd_buf4_close(&G1);

      dpd_buf4_init(&G1, CC_GAMMA, 0, 5, 5, 5, 5, 0, "2 Gabcd - Gabdc");
      dpd_buf4_copy(&G1, CC_GAMMA, "GAbCd");
      if(!params.aobasis) {  
	dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
	two_energy = dpd_buf4_dot(&B, &G1);
	dpd_buf4_close(&B);
	fprintf(outfile, "\tABCD energy                = %20.15f\n", two_energy);
	total_two_energy += two_energy;
      }

      dpd_buf4_close(&G1);

      if(!params.aobasis) {
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
    }

  }} // namespace psi::ccdensity
