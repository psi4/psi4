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
#include <strings.h>
#include <string.h>
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* DEANTI_UHF(): Convert the UHF two-particle density from Dirac to
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
**
** This code is rather different from its RHF/ROHF counterpart in that
** we must keep the three spin components of Gpqrs separate (i.e., no
** spin-adaptation is allowed for UHF cases).
** */

void deanti_UHF(struct RHO_Params rho_params)
{
  dpdfile2 h, d;
  dpdbuf4 G, G2, A, B, C, D, E, F;
  double one_energy=0.0, two_energy=0.0, total_two_energy=0.0;

  outfile->Printf( "\n\tEnergies re-computed from Mulliken density:\n");
  outfile->Printf(   "\t-------------------------------------------\n");

  global_dpd_->file2_init(&d, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
  global_dpd_->file2_init(&h, PSIF_CC_OEI, 0, 0, 0, "h(I,J)");
  one_energy += global_dpd_->file2_dot(&d, &h);
  global_dpd_->file2_close(&h);
  global_dpd_->file2_close(&d);

  global_dpd_->file2_init(&d, PSIF_CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
  global_dpd_->file2_init(&h, PSIF_CC_OEI, 0, 2, 2, "h(i,j)");
  one_energy += global_dpd_->file2_dot(&d, &h);
  global_dpd_->file2_close(&h);
  global_dpd_->file2_close(&d);

  global_dpd_->file2_init(&d, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
  global_dpd_->file2_init(&h, PSIF_CC_OEI, 0, 1, 1, "h(A,B)");
  one_energy += global_dpd_->file2_dot(&d, &h);
  global_dpd_->file2_close(&h);
  global_dpd_->file2_close(&d);

  global_dpd_->file2_init(&d, PSIF_CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
  global_dpd_->file2_init(&h, PSIF_CC_OEI, 0, 3, 3, "h(a,b)");
  one_energy += global_dpd_->file2_dot(&d, &h);
  global_dpd_->file2_close(&h);
  global_dpd_->file2_close(&d);

  global_dpd_->file2_init(&d, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
  global_dpd_->file2_init(&h, PSIF_CC_OEI, 0, 0, 1, "h(I,A)");
  one_energy += global_dpd_->file2_dot(&d, &h);
  global_dpd_->file2_close(&h);
  global_dpd_->file2_close(&d);

  global_dpd_->file2_init(&d, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
  global_dpd_->file2_init(&h, PSIF_CC_OEI, 0, 2, 3, "h(i,a)");
  one_energy += global_dpd_->file2_dot(&d, &h);
  global_dpd_->file2_close(&h);
  global_dpd_->file2_close(&d);

  global_dpd_->file2_init(&d, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
  global_dpd_->file2_init(&h, PSIF_CC_OEI, 0, 0, 1, "h(I,A)");
  one_energy += global_dpd_->file2_dot(&d, &h);
  global_dpd_->file2_close(&h);
  global_dpd_->file2_close(&d);

  global_dpd_->file2_init(&d, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
  global_dpd_->file2_init(&h, PSIF_CC_OEI, 0, 2, 3, "h(i,a)");
  one_energy += global_dpd_->file2_dot(&d, &h);
  global_dpd_->file2_close(&h);
  global_dpd_->file2_close(&d);

  outfile->Printf( "\tOne-electron energy        = %20.15f\n", one_energy);


  /* G(Ij,Kl) = 1/2 G(Ij,Kl) */
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 22, 22, 22, 0, "GIjKl");
  global_dpd_->buf4_scm(&G, 0.5);
  global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
  two_energy = 2*global_dpd_->buf4_dot(&A, &G);
  global_dpd_->buf4_close(&A);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tIjKl energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  /* G(IJ,KL) = 1/2 G(IJ,KL) */
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 2, 2, 0, "GIJKL");
  global_dpd_->buf4_scm(&G, 0.5);
  global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <IJ|KL>");
  two_energy = global_dpd_->buf4_dot(&A, &G);
  global_dpd_->buf4_close(&A);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tIJKL energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  /* G(ij,kl) = 1/2 G(ij,kl) */
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 12, 12, 0, "Gijkl");
  global_dpd_->buf4_scm(&G, 0.5);
  global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 10, 10, 10, 10, 0, "A <ij|kl>");
  two_energy = global_dpd_->buf4_dot(&A, &G);
  global_dpd_->buf4_close(&A);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tijkl energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  /* No change to G(IJ,KA) components */
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
  two_energy = 2*global_dpd_->buf4_dot(&E, &G);
  global_dpd_->buf4_close(&E);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
  two_energy += 2*global_dpd_->buf4_dot(&E, &G);
  global_dpd_->buf4_close(&E);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tIjKa+iJkA energy           = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 20, 2, 20, 0, "GIJKA");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 20, 0, 20, 0, "E <IJ|KA>");
  two_energy = 2*global_dpd_->buf4_dot(&E, &G);
  global_dpd_->buf4_close(&E);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tIJKA energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 30, 12, 30, 0, "Gijka");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 10, 30, 10, 30, 0, "E <ij|ka>");
  two_energy = 2*global_dpd_->buf4_dot(&E, &G);
  global_dpd_->buf4_close(&E);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tijka energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  /* G(Ij,Ab) <-- G(Ij,Ab) - G(Ib,jA) */
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 24, 27, 24, 27, 0, "GIbjA");
  global_dpd_->buf4_sort(&G2, PSIF_CC_TMP0, prsq, 22, 28, "GIbjA (Ij,Ab)");
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "GIbjA (Ij,Ab)");
  global_dpd_->buf4_axpy(&G2, &G, -1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  two_energy = 2 * global_dpd_->buf4_dot(&D, &G);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tIjAb energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  /* G(IJ,AB) <-- G(IJ,AB) - G(IB,JA) */
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 2, 7, 0, "GIJAB");
  global_dpd_->buf4_copy(&G, PSIF_CC_GAMMA, "G(IJ,AB)");
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "G(IJ,AB)");
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 20, 20, 20, 20, 0, "GIBJA");
  global_dpd_->buf4_sort(&G2, PSIF_CC_TMP0, prsq, 0, 5, "GIBJA (IJ,AB)");
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "GIBJA (IJ,AB)");
  global_dpd_->buf4_axpy(&G2, &G, -1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ|AB>");
  two_energy = global_dpd_->buf4_dot(&D, &G);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tIJAB energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  /* G(ij,ab) <-- G(ij,ab) - G(ib,ja) */
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 15, 12, 17, 0, "Gijab");
  global_dpd_->buf4_copy(&G, PSIF_CC_GAMMA, "G(ij,ab)");
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 15, 10, 15, 0, "G(ij,ab)");
  global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 30, 30, 30, 30, 0, "Gibja");
  global_dpd_->buf4_sort(&G2, PSIF_CC_TMP0, prsq, 10, 15, "Gibja (ij,ab)");
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&G2, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "Gibja (ij,ab)");
  global_dpd_->buf4_axpy(&G2, &G, -1.0);
  global_dpd_->buf4_close(&G2);
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij|ab>");
  two_energy = global_dpd_->buf4_dot(&D, &G);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tijab energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  /* No change to G(IB,JA) components */
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 20, 20, 20, 20, 0, "GIBJA");
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA|JB>");
  two_energy = global_dpd_->buf4_dot(&C, &G);
  global_dpd_->buf4_close(&C);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tIBJA energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 30, 30, 30, 30, 0, "Gibja");
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia|jb>");
  two_energy = global_dpd_->buf4_dot(&C, &G);
  global_dpd_->buf4_close(&C);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tibja energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 24, 24, 24, 24, 0, "GIbJa");
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
  two_energy = global_dpd_->buf4_dot(&C, &G);
  global_dpd_->buf4_close(&C);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 27, 27, 27, 27, 0, "GiBjA");
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
  two_energy += global_dpd_->buf4_dot(&C, &G);
  global_dpd_->buf4_close(&C);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tiBjA+IbJa energy           = %20.15f\n", two_energy);
  total_two_energy += two_energy;


  /* No change to G(CI,AB) components */
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  two_energy = 2*global_dpd_->buf4_dot(&F, &G);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
  two_energy += 2*global_dpd_->buf4_dot(&F, &G);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tcIaB+CiAb energy           = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 21, 5, 21, 7, 0, "GCIAB");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 21, 5, 21, 5, 0, "F <AI|BC>");
  two_energy = 2*global_dpd_->buf4_dot(&F, &G);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tCIAB energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 31, 15, 31, 17, 0, "Gciab");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 31, 15, 31, 15, 0, "F <ai|bc>");
  two_energy = 2*global_dpd_->buf4_dot(&F, &G);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tciab energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  /* G(Ab,Cd) = 1/2 G(Ab,Cd) */
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 28, 28, 28, 28, 0, "GAbCd");
  global_dpd_->buf4_scm(&G, 0.5);
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
  two_energy = 2*global_dpd_->buf4_dot(&B, &G);
  global_dpd_->buf4_close(&B);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tAbCd energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  /* G(AB,CD) = 1/2 G(AB,CD) */
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 5, 5, 7, 7, 0, "GABCD");
  global_dpd_->buf4_scm(&G, 0.5);
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <AB|CD>");
  two_energy = global_dpd_->buf4_dot(&B, &G);
  global_dpd_->buf4_close(&B);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tABCD energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  /* G(ab,cd) = 1/2 G(ab,cd) */
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 15, 15, 17, 17, 0, "Gabcd");
  global_dpd_->buf4_scm(&G, 0.5);
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 15, 15, 15, 15, 0, "B <ab|cd>");
  two_energy = global_dpd_->buf4_dot(&B, &G);
  global_dpd_->buf4_close(&B);
  global_dpd_->buf4_close(&G);
  outfile->Printf( "\tabcd energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  outfile->Printf( "\tTotal two-electron energy  = %20.15f\n", total_two_energy);

  outfile->Printf( "\t%-7s correlation energy = %20.15f\n", params.wfn == "CCSD_T" ? "CCSD(T)" : params.wfn.c_str(),
	  one_energy + total_two_energy);
  outfile->Printf( "\tTotal %-7s energy       = %20.15f\n", params.wfn == "CCSD_T" ? "CCSD(T)" : params.wfn.c_str(),
	  one_energy + total_two_energy + moinfo.eref);
}

}} // namespace psi::ccdensity
