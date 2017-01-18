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
#include <cstdio>
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* DEANTI_ROHF(): Convert the ROHF two-particle density from Dirac to
** Mulliken ordering.  The original, Fock-adjusted density (see the
** comments in fock.c) corresponds to a two-electron energy (or energy
** derivative) expression of the form:
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

void deanti_ROHF(struct RHO_Params rho_params)
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
    one_energy += global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "h(i,j)");
    one_energy += global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "h(a,b)");
    one_energy += global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "h(a,b)");
    one_energy += global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");
    one_energy += global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");
    one_energy += global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");
    one_energy += global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");
    one_energy += global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    outfile->Printf( "\tOne-electron energy        = %20.15f\n", one_energy);

  }

  /* G(Ij,Kl) <-- 1/2 G(IJ,KL) + 1/2 G(ij,kl) + 1/2 G(Ij,Kl) + 1/2 G(iJ,kL) */
  global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  global_dpd_->buf4_sort(&G1, PSIF_CC_TMP0, qprs, 0, 0, "GjIKl");
  global_dpd_->buf4_init(&G2, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "GjIKl");
  global_dpd_->buf4_sort(&G2, PSIF_CC_TMP0, pqsr, 0, 0, "GiJkL");
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "GiJkL");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 0, 0, 2, 2, 0, "GIJKL");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 0, 0, 2, 2, 0, "Gijkl");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_scm(&G1, 0.5);

  if(!params.aobasis) {
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    two_energy = global_dpd_->buf4_dot(&A, &G1);
    global_dpd_->buf4_close(&A);
    outfile->Printf( "\tIJKL energy                = %20.15f\n", two_energy);
    total_two_energy += two_energy;
  }

  global_dpd_->buf4_close(&G1);

  /* G(IJ,KA) <-- G(IJ,KA) + G(ij,ka) + G(Ij,Ka) + G(iJ,kA) */
  global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 0, 10, 2, 10, 0, "GIJKA");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 0, 10, 2, 10, 0, "Gijka");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);

  if(!params.aobasis) {
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    two_energy = global_dpd_->buf4_dot(&E, &G1);
    global_dpd_->buf4_close(&E);
    outfile->Printf( "\tIJKA energy                = %20.15f\n", 2*two_energy);
    total_two_energy += 2*two_energy;
  }

  global_dpd_->buf4_close(&G1);

  /* G(Ij,Ab) <-- G(Ij,Ab) - G(Ib,jA) */
  global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
  global_dpd_->buf4_sort(&G2, PSIF_CC_TMP0, prsq, 0, 5, "GIbjA (Ij,Ab)");
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "GIbjA (Ij,Ab)");
  global_dpd_->buf4_axpy(&G2, &G1, -1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_close(&G1);

  /* G(IJ,AB) <-- G(IJ,AB) - G(IB,JA) */
  global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 0, 5, 2, 7, 0, "GIJAB");
  global_dpd_->buf4_copy(&G1, PSIF_CC_TMP0, "G(IJ,AB)");
  global_dpd_->buf4_close(&G1);
  global_dpd_->buf4_init(&G1, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "G(IJ,AB)");
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
  global_dpd_->buf4_sort(&G2, PSIF_CC_TMP0, prsq, 0, 5, "GIBJA (IJ,AB)");
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "GIBJA (IJ,AB)");
  global_dpd_->buf4_axpy(&G2, &G1, -1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_close(&G1);

  /* G(ij,ab) <-- G(ij,ab) - G(ib,ja) */
  global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 0, 5, 2, 7, 0, "Gijab");
  global_dpd_->buf4_copy(&G1, PSIF_CC_TMP0, "G(ij,ab)");
  global_dpd_->buf4_close(&G1);
  global_dpd_->buf4_init(&G1, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "G(ij,ab)");
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
  global_dpd_->buf4_sort(&G2, PSIF_CC_TMP0, prsq, 0, 5, "Gibja (ij,ab)");
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Gibja (ij,ab)");
  global_dpd_->buf4_axpy(&G2, &G1, -1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_close(&G1);

  /* G(IJ,AB) <-- G(IJ,AB) + G(ij,ab) + G(Ij,Ab) + G(iJ,aB) */
  global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  global_dpd_->buf4_sort(&G1, PSIF_CC_TMP0, qpsr, 0, 5, "G(jI,bA)");
  global_dpd_->buf4_init(&G2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "G(jI,bA)");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "G(IJ,AB)");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "G(ij,ab)");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);

  if(!params.aobasis) {
    global_dpd_->buf4_init(&DInts, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    two_energy = global_dpd_->buf4_dot(&DInts, &G1);
    global_dpd_->buf4_close(&DInts);
    outfile->Printf( "\tIJAB energy                = %20.15f\n", two_energy);
    total_two_energy += two_energy;
  }

  global_dpd_->buf4_close(&G1);

  /* G(IB,JA) <-- G(IB,JA) + G(ib,ja) + G(Ib,Ja) + G(iB,jA) */
  global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);

  if(!params.aobasis) {
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    two_energy = global_dpd_->buf4_dot(&C, &G1);
    global_dpd_->buf4_close(&C);
    outfile->Printf( "\tIBJA energy                = %20.15f\n", two_energy);
    total_two_energy += two_energy;
  }

  global_dpd_->buf4_close(&G1);

  /* G(CI,AB) <-- G(CI,AB) + G(ci,ab) + G(Ci,Ab) + G(cI,aB) */
  global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 11, 5, 11, 7, 0, "GCIAB");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 11, 5, 11, 7, 0, "Gciab");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);

  if(!params.aobasis) {
    global_dpd_->buf4_init(&FInts, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_sort(&FInts, PSIF_CC_TMP0, qprs, 11, 5, "F(CI,BA)");
    global_dpd_->buf4_close(&FInts);
    global_dpd_->buf4_init(&FInts, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "F(CI,BA)");
    global_dpd_->buf4_sort(&FInts, PSIF_CC_TMP1, pqsr, 11, 5, "F(CI,AB)");
    global_dpd_->buf4_close(&FInts);
    global_dpd_->buf4_init(&FInts, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "F(CI,AB)");
    two_energy = global_dpd_->buf4_dot(&FInts, &G1);
    global_dpd_->buf4_close(&FInts);
    outfile->Printf( "\tCIAB energy                = %20.15f\n", 2*two_energy);
    total_two_energy += 2*two_energy;
  }

  global_dpd_->buf4_close(&G1);

  /* G(Ab,Cd) <-- 1/2 G(AB,CD) + 1/2 G(ab,cd) + G(Ab,Cd) */
  global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
  global_dpd_->buf4_sort(&G1, PSIF_CC_TMP0, qprs, 5, 5, "G(bA,Cd)");
  global_dpd_->buf4_init(&G2, PSIF_CC_TMP0, 0, 5, 5, 5, 5, 0, "G(bA,Cd)");
  global_dpd_->buf4_sort(&G2, PSIF_CC_TMP1, pqsr, 5, 5, "G(aB,cD)");
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_TMP1, 0, 5, 5, 5, 5, 0, "G(aB,cD)");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 5, 5, 7, 7, 0, "GABCD");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 5, 5, 7, 7, 0, "Gabcd");
  global_dpd_->buf4_axpy(&G2, &G1, 1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_scm(&G1, 0.5);

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
      outfile->Printf( "\tCCSD correlation energy    = %20.15f\n",
	      one_energy + total_two_energy);
      outfile->Printf( "\tTotal CCSD energy          = %20.15f\n",
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
