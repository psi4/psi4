/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <strings.h>
#include <string.h>
#include <libiwl/iwl.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
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

  fprintf(outfile, "\n\tEnergies re-computed from Mulliken density:\n");
  fprintf(outfile,   "\t-------------------------------------------\n");

  dpd_file2_init(&d, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
  dpd_file2_init(&h, CC_OEI, 0, 0, 0, "h(I,J)");
  one_energy += dpd_file2_dot(&d, &h);
  dpd_file2_close(&h);
  dpd_file2_close(&d);

  dpd_file2_init(&d, CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
  dpd_file2_init(&h, CC_OEI, 0, 2, 2, "h(i,j)");
  one_energy += dpd_file2_dot(&d, &h);
  dpd_file2_close(&h);
  dpd_file2_close(&d);

  dpd_file2_init(&d, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
  dpd_file2_init(&h, CC_OEI, 0, 1, 1, "h(A,B)");
  one_energy += dpd_file2_dot(&d, &h);
  dpd_file2_close(&h);
  dpd_file2_close(&d);

  dpd_file2_init(&d, CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
  dpd_file2_init(&h, CC_OEI, 0, 3, 3, "h(a,b)");
  one_energy += dpd_file2_dot(&d, &h);
  dpd_file2_close(&h);
  dpd_file2_close(&d);

  dpd_file2_init(&d, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
  dpd_file2_init(&h, CC_OEI, 0, 0, 1, "h(I,A)");
  one_energy += dpd_file2_dot(&d, &h);
  dpd_file2_close(&h);
  dpd_file2_close(&d);

  dpd_file2_init(&d, CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
  dpd_file2_init(&h, CC_OEI, 0, 2, 3, "h(i,a)");
  one_energy += dpd_file2_dot(&d, &h);
  dpd_file2_close(&h);
  dpd_file2_close(&d);

  dpd_file2_init(&d, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
  dpd_file2_init(&h, CC_OEI, 0, 0, 1, "h(I,A)");
  one_energy += dpd_file2_dot(&d, &h);
  dpd_file2_close(&h);
  dpd_file2_close(&d);

  dpd_file2_init(&d, CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
  dpd_file2_init(&h, CC_OEI, 0, 2, 3, "h(i,a)");
  one_energy += dpd_file2_dot(&d, &h);
  dpd_file2_close(&h);
  dpd_file2_close(&d);

  fprintf(outfile, "\tOne-electron energy        = %20.15f\n", one_energy);
  fflush(outfile);

  /* G(Ij,Kl) = 1/2 G(Ij,Kl) */
  dpd_buf4_init(&G, CC_GAMMA, 0, 22, 22, 22, 22, 0, "GIjKl");
  dpd_buf4_scm(&G, 0.5);
  dpd_buf4_init(&A, CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
  two_energy = 2*dpd_buf4_dot(&A, &G);
  dpd_buf4_close(&A);
  dpd_buf4_close(&G);
  fprintf(outfile, "\tIjKl energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  /* G(IJ,KL) = 1/2 G(IJ,KL) */
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 2, 2, 0, "GIJKL");
  dpd_buf4_scm(&G, 0.5);
  dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <IJ|KL>");
  two_energy = dpd_buf4_dot(&A, &G);
  dpd_buf4_close(&A);
  dpd_buf4_close(&G);
  fprintf(outfile, "\tIJKL energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  /* G(ij,kl) = 1/2 G(ij,kl) */
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 12, 12, 0, "Gijkl");
  dpd_buf4_scm(&G, 0.5);
  dpd_buf4_init(&A, CC_AINTS, 0, 10, 10, 10, 10, 0, "A <ij|kl>");
  two_energy = dpd_buf4_dot(&A, &G);
  dpd_buf4_close(&A);
  dpd_buf4_close(&G);
  fprintf(outfile, "\tijkl energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  /* No change to G(IJ,KA) components */
  dpd_buf4_init(&G, CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
  dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
  two_energy = 2*dpd_buf4_dot(&E, &G);
  dpd_buf4_close(&E);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
  dpd_buf4_init(&E, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
  two_energy += 2*dpd_buf4_dot(&E, &G);
  dpd_buf4_close(&E);
  dpd_buf4_close(&G);
  fprintf(outfile, "\tIjKa+iJkA energy           = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 20, 2, 20, 0, "GIJKA");
  dpd_buf4_init(&E, CC_EINTS, 0, 0, 20, 0, 20, 0, "E <IJ|KA>");
  two_energy = 2*dpd_buf4_dot(&E, &G);
  dpd_buf4_close(&E);
  dpd_buf4_close(&G);
  fprintf(outfile, "\tIJKA energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 30, 12, 30, 0, "Gijka");
  dpd_buf4_init(&E, CC_EINTS, 0, 10, 30, 10, 30, 0, "E <ij|ka>");
  two_energy = 2*dpd_buf4_dot(&E, &G);
  dpd_buf4_close(&E);
  dpd_buf4_close(&G);
  fprintf(outfile, "\tijka energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  /* G(Ij,Ab) <-- G(Ij,Ab) - G(Ib,jA) */
  dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
  dpd_buf4_init(&G2, CC_GAMMA, 0, 24, 27, 24, 27, 0, "GIbjA");
  dpd_buf4_sort(&G2, CC_TMP0, prsq, 22, 28, "GIbjA (Ij,Ab)");
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_TMP0, 0, 22, 28, 22, 28, 0, "GIbjA (Ij,Ab)");
  dpd_buf4_axpy(&G2, &G, -1.0);
  dpd_buf4_close(&G2);
  dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  two_energy = 2 * dpd_buf4_dot(&D, &G);
  dpd_buf4_close(&D);
  dpd_buf4_close(&G);
  fprintf(outfile, "\tIjAb energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;
  
  /* G(IJ,AB) <-- G(IJ,AB) - G(IB,JA) */
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 2, 7, 0, "GIJAB");
  dpd_buf4_copy(&G, CC_GAMMA, "G(IJ,AB)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "G(IJ,AB)");
  dpd_buf4_init(&G2, CC_GAMMA, 0, 20, 20, 20, 20, 0, "GIBJA");
  dpd_buf4_sort(&G2, CC_TMP0, prsq, 0, 5, "GIBJA (IJ,AB)");
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_TMP0, 0, 0, 5, 0, 5, 0, "GIBJA (IJ,AB)");
  dpd_buf4_axpy(&G2, &G, -1.0); 
  dpd_buf4_close(&G2);
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ|AB>");
  two_energy = dpd_buf4_dot(&D, &G);
  dpd_buf4_close(&D);
  dpd_buf4_close(&G);
  fprintf(outfile, "\tIJAB energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;
  
  /* G(ij,ab) <-- G(ij,ab) - G(ib,ja) */
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 15, 12, 17, 0, "Gijab");
  dpd_buf4_copy(&G, CC_GAMMA, "G(ij,ab)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 15, 10, 15, 0, "G(ij,ab)");
  dpd_buf4_init(&G2, CC_GAMMA, 0, 30, 30, 30, 30, 0, "Gibja");
  dpd_buf4_sort(&G2, CC_TMP0, prsq, 10, 15, "Gibja (ij,ab)");
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_TMP0, 0, 10, 15, 10, 15, 0, "Gibja (ij,ab)");
  dpd_buf4_axpy(&G2, &G, -1.0);
  dpd_buf4_close(&G2);
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij|ab>");
  two_energy = dpd_buf4_dot(&D, &G);
  dpd_buf4_close(&D);
  dpd_buf4_close(&G);
  fprintf(outfile, "\tijab energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;
  
  /* No change to G(IB,JA) components */
  dpd_buf4_init(&G, CC_GAMMA, 0, 20, 20, 20, 20, 0, "GIBJA");
  dpd_buf4_init(&C, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA|JB>");
  two_energy = dpd_buf4_dot(&C, &G);
  dpd_buf4_close(&C);
  dpd_buf4_close(&G);
  fprintf(outfile, "\tIBJA energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;
  
  dpd_buf4_init(&G, CC_GAMMA, 0, 30, 30, 30, 30, 0, "Gibja");
  dpd_buf4_init(&C, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia|jb>");
  two_energy = dpd_buf4_dot(&C, &G);
  dpd_buf4_close(&C);
  dpd_buf4_close(&G); 
  fprintf(outfile, "\tibja energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;
  
  dpd_buf4_init(&G, CC_GAMMA, 0, 24, 24, 24, 24, 0, "GIbJa");
  dpd_buf4_init(&C, CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
  two_energy = dpd_buf4_dot(&C, &G);
  dpd_buf4_close(&C);
  dpd_buf4_close(&G); 

  dpd_buf4_init(&G, CC_GAMMA, 0, 27, 27, 27, 27, 0, "GiBjA");
  dpd_buf4_init(&C, CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
  two_energy += dpd_buf4_dot(&C, &G);
  dpd_buf4_close(&C);
  dpd_buf4_close(&G);
  fprintf(outfile, "\tiBjA+IbJa energy           = %20.15f\n", two_energy);
  total_two_energy += two_energy;
  

  /* No change to G(CI,AB) components */
  dpd_buf4_init(&G, CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
  dpd_buf4_init(&F, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  two_energy = 2*dpd_buf4_dot(&F, &G);
  dpd_buf4_close(&F);
  dpd_buf4_close(&G);
  
  dpd_buf4_init(&G, CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
  dpd_buf4_init(&F, CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
  two_energy += 2*dpd_buf4_dot(&F, &G);
  dpd_buf4_close(&F);
  dpd_buf4_close(&G); 
  fprintf(outfile, "\tcIaB+CiAb energy           = %20.15f\n", two_energy);
  total_two_energy += two_energy;
  
  dpd_buf4_init(&G, CC_GAMMA, 0, 21, 5, 21, 7, 0, "GCIAB");
  dpd_buf4_init(&F, CC_FINTS, 0, 21, 5, 21, 5, 0, "F <AI|BC>");
  two_energy = 2*dpd_buf4_dot(&F, &G);
  dpd_buf4_close(&F);
  dpd_buf4_close(&G); 
  fprintf(outfile, "\tCIAB energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;
  
  dpd_buf4_init(&G, CC_GAMMA, 0, 31, 15, 31, 17, 0, "Gciab");
  dpd_buf4_init(&F, CC_FINTS, 0, 31, 15, 31, 15, 0, "F <ai|bc>");
  two_energy = 2*dpd_buf4_dot(&F, &G);
  dpd_buf4_close(&F);
  dpd_buf4_close(&G);
  fprintf(outfile, "\tciab energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;
  
  /* G(Ab,Cd) = 1/2 G(Ab,Cd) */
  dpd_buf4_init(&G, CC_GAMMA, 0, 28, 28, 28, 28, 0, "GAbCd");
  dpd_buf4_scm(&G, 0.5);
  dpd_buf4_init(&B, CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
  two_energy = 2*dpd_buf4_dot(&B, &G);
  dpd_buf4_close(&B);
  dpd_buf4_close(&G);
  fprintf(outfile, "\tAbCd energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;
  
  /* G(AB,CD) = 1/2 G(AB,CD) */
  dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 7, 7, 0, "GABCD");
  dpd_buf4_scm(&G, 0.5);
  dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <AB|CD>");
  two_energy = dpd_buf4_dot(&B, &G);
  dpd_buf4_close(&B);
  dpd_buf4_close(&G);
  fprintf(outfile, "\tABCD energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;
  
  /* G(ab,cd) = 1/2 G(ab,cd) */
  dpd_buf4_init(&G, CC_GAMMA, 0, 15, 15, 17, 17, 0, "Gabcd");
  dpd_buf4_scm(&G, 0.5);
  dpd_buf4_init(&B, CC_BINTS, 0, 15, 15, 15, 15, 0, "B <ab|cd>");
  two_energy = dpd_buf4_dot(&B, &G);
  dpd_buf4_close(&B);
  dpd_buf4_close(&G);
  fprintf(outfile, "\tabcd energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;
  
  fprintf(outfile, "\tTotal two-electron energy  = %20.15f\n", total_two_energy);

  fprintf(outfile, "\t%-7s correlation energy = %20.15f\n", !strcmp(params.wfn,"CCSD_T") ? "CCSD(T)" : params.wfn,
	  one_energy + total_two_energy);
  fprintf(outfile, "\tTotal %-7s energy       = %20.15f\n", !strcmp(params.wfn,"CCSD_T") ? "CCSD(T)" : params.wfn,
	  one_energy + total_two_energy + moinfo.eref);
}

}} // namespace psi::ccdensity
