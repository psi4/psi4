/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <stdio.h>
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
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
	outfile->Printf( "\n\tEnergies re-computed from Mulliken density:\n");
	outfile->Printf(   "\t-------------------------------------------\n");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "h(i,j)");
	one_energy += 2.0*global_dpd_->file2_dot(&D, &F);
	global_dpd_->file2_close(&F);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "h(a,b)");
	one_energy += 2.0*global_dpd_->file2_dot(&D, &F);
	global_dpd_->file2_close(&F);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");
	one_energy += 2.0*global_dpd_->file2_dot(&D, &F);
	global_dpd_->file2_close(&F);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");
	one_energy += 2.0*global_dpd_->file2_dot(&D, &F);
	global_dpd_->file2_close(&F);
	global_dpd_->file2_close(&D);

	outfile->Printf( "\tOne-electron energy        = %20.15f\n", one_energy);

      }

      /* E_ijkl = (2 Gijkl - Gijlk) <ij|kl> */
      global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
      global_dpd_->buf4_scmcopy(&G1, PSIF_CC_GAMMA, "2 Gijkl - Gijlk", 2);
      global_dpd_->buf4_sort_axpy(&G1, PSIF_CC_GAMMA, pqsr, 0, 0, "2 Gijkl - Gijlk", -1);
      global_dpd_->buf4_close(&G1);

      global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "2 Gijkl - Gijlk");
      global_dpd_->buf4_copy(&G1, PSIF_CC_GAMMA, "GIjKl");
      if(!params.aobasis) {
	global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	two_energy = global_dpd_->buf4_dot(&A, &G1);
	global_dpd_->buf4_close(&A);
	outfile->Printf( "\tIJKL energy                = %20.15f\n", two_energy);
	total_two_energy += two_energy;
      }
      global_dpd_->buf4_close(&G1);

      /* E_ijka = 2 [ (2 Gijka - Gjika) <ij|ka> ] */
      /* NB: GIjKa is scaled by 1/2 in Gijka.cc. */
      global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
      global_dpd_->buf4_scmcopy(&G1, PSIF_CC_GAMMA, "2 Gijka - Gjika", 2);
      global_dpd_->buf4_sort_axpy(&G1, PSIF_CC_GAMMA, qprs, 0, 10, "2 Gijka - Gjika", -1);
      global_dpd_->buf4_close(&G1);

      global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "2 Gijka - Gjika");
      /* The factor of 2 here is necessary because Gijka is multiplied by 1/2 in Gijka.cc */
      global_dpd_->buf4_scmcopy(&G1, PSIF_CC_GAMMA, "GIjKa", 2);
      if(!params.aobasis) {
	global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	two_energy = global_dpd_->buf4_dot(&E, &G1);
	global_dpd_->buf4_close(&E);
	/* The factor of 4 here is necessary because Gijka is multiplied by 1/2 in Gijka.cc */
	outfile->Printf( "\tIJKA energy                = %20.15f\n", 4*two_energy);
	total_two_energy += 4*two_energy;
      }
      global_dpd_->buf4_close(&G1);

      /* E_ijab = [ (2 Gijab - Gijba) - GIBJA - GIbjA ] <ij|ab> */
      global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
      global_dpd_->buf4_scmcopy(&G1, PSIF_CC_GAMMA, "2 Gijab - Gijba", 2);
      global_dpd_->buf4_sort_axpy(&G1, PSIF_CC_GAMMA, pqsr, 0, 5, "2 Gijab - Gijba", -1);
      global_dpd_->buf4_close(&G1);

      global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
      global_dpd_->buf4_sort_axpy(&G2, PSIF_CC_GAMMA, prsq, 0, 5, "2 Gijab - Gijba", -1);
      global_dpd_->buf4_close(&G2);
      global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
      global_dpd_->buf4_sort_axpy(&G2, PSIF_CC_GAMMA, prsq, 0, 5, "2 Gijab - Gijba", -1);
      global_dpd_->buf4_close(&G2);

      global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "2 Gijab - Gijba");
      global_dpd_->buf4_scmcopy(&G1, PSIF_CC_GAMMA, "GIjAb", 2);
      if(!params.aobasis) {
	global_dpd_->buf4_init(&DInts, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	two_energy = 2.0 * global_dpd_->buf4_dot(&DInts, &G1);
	global_dpd_->buf4_close(&DInts);
	outfile->Printf( "\tIJAB energy                = %20.15f\n", two_energy);
	total_two_energy += two_energy;
      }

      global_dpd_->buf4_close(&G1);

      /* G'_IbJa = 2 GIBJA + 2 GIbJa */
      global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
      global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
      global_dpd_->buf4_axpy(&G2, &G1, 1);
      global_dpd_->buf4_close(&G2);
      //  dpd_buf4_scm(&G1, 2);
      if(!params.aobasis) {
	global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	two_energy = 2.0 * global_dpd_->buf4_dot(&C, &G1);
	global_dpd_->buf4_close(&C);
	outfile->Printf( "\tIBJA energy                = %20.15f\n", two_energy);
	total_two_energy += two_energy;
      }

      global_dpd_->buf4_close(&G1);

      /* E_ciab = 2 [ (2 Gciab - Gciba) <ci|ab> ] */
      /* NB: Gciab is scaled by 1/2 in Gciab.cc. */
      global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
      global_dpd_->buf4_scmcopy(&G1, PSIF_CC_GAMMA, "2 Gciab - Gciba", 2);
      global_dpd_->buf4_sort_axpy(&G1, PSIF_CC_GAMMA, pqsr, 11, 5, "2 Gciab - Gciba", -1);
      global_dpd_->buf4_close(&G1);

      global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "2 Gciab - Gciba");
      /* The factor of 2 here is necessary because Gciab is multiplied by 1/2 in Gciab.cc */
      global_dpd_->buf4_scmcopy(&G1, PSIF_CC_GAMMA, "GCiAb", 2);
      if(!params.aobasis) {
	global_dpd_->buf4_init(&FInts, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->buf4_sort(&FInts, PSIF_CC_FINTS, qpsr, 11, 5, "F <ai|bc>");
	global_dpd_->buf4_close(&FInts);
	global_dpd_->buf4_init(&FInts, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
	two_energy = global_dpd_->buf4_dot(&FInts, &G1);
	global_dpd_->buf4_close(&FInts);
	/* The factor of 4 here is necessary because Gciab is multiplied by 1/2 in Gciab.cc */
	outfile->Printf( "\tCIAB energy                = %20.15f\n", 4*two_energy);
	total_two_energy += 4*two_energy;
      }

      global_dpd_->buf4_close(&G1);

      /* E_abcd = (2 Gabcd - Gabdc) <ab|cd> */
      global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");

      global_dpd_->buf4_scmcopy(&G1, PSIF_CC_GAMMA, "2 Gabcd - Gabdc", 2);
      global_dpd_->buf4_sort_axpy(&G1, PSIF_CC_GAMMA, pqsr, 5, 5, "2 Gabcd - Gabdc", -1);
      global_dpd_->buf4_close(&G1);

      global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 5, 5, 5, 5, 0, "2 Gabcd - Gabdc");
      global_dpd_->buf4_copy(&G1, PSIF_CC_GAMMA, "GAbCd");
      if(!params.aobasis) {
	global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
	two_energy = global_dpd_->buf4_dot(&B, &G1);
	global_dpd_->buf4_close(&B);
	outfile->Printf( "\tABCD energy                = %20.15f\n", two_energy);
	total_two_energy += two_energy;
      }

      global_dpd_->buf4_close(&G1);

      if(!params.aobasis) {
	outfile->Printf( "\tTotal two-electron energy  = %20.15f\n", total_two_energy);
	if (params.ground) {
	  //outfile->Printf( "\tCCSD correlation energy    = %20.15f\n",
	  //	  one_energy + total_two_energy);
	  //outfile->Printf( "\tTotal CCSD energy          = %20.15f\n",
	  //	  one_energy + total_two_energy + moinfo.eref);
	  outfile->Printf( "\t%7s correlation energy = %20.15f\n", params.wfn == "CCSD_T" ? "CCSD(T)" : params.wfn.c_str(),
                one_energy + total_two_energy);
          outfile->Printf( "\tTotal %7s energy       = %20.15f\n", params.wfn == "CCSD_T" ? "CCSD(T)" : params.wfn.c_str(),
                one_energy + total_two_energy + moinfo.eref);	  
	}
	else {
	  outfile->Printf( "\tTotal EOM CCSD correlation energy        = %20.15f\n",
		  one_energy + total_two_energy);
	  outfile->Printf( "\tCCSD correlation + EOM excitation energy = %20.15f\n",
		  moinfo.ecc + params.cceom_energy);
	  outfile->Printf( "\tTotal EOM CCSD energy                    = %20.15f\n",
		  one_energy + total_two_energy + moinfo.eref);
	}
      }
    }

  }} // namespace psi::ccdensity
