/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libiwl/iwl.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

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

void rhf_sf_deanti(void);
void uhf_sf_deanti(void);

void deanti(void)
{
  if(params.ref == 0) rhf_sf_deanti();
  else if(params.ref == 2) uhf_sf_deanti();
}

void rhf_sf_deanti(void)
{
  dpdbuf4 G1, G2;
  dpdfile2 D, F;
  double one_energy = 0.0;
  double two_energy = 0.0; 
  double total_two_energy = 0.0;
  dpdbuf4 A, B, C, DInts, E, FInts;

  fprintf(outfile, "\n\tEnergies re-computed from Mulliken density:\n");
  fprintf(outfile,   "\t-------------------------------------------\n");

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_file2_init(&F, CC_OEI, 0, 0, 0, "h(i,j)");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_file2_init(&F, CC_OEI, 0, 0, 0, "h(i,j)");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_file2_init(&F, CC_OEI, 0, 1, 1, "h(a,b)");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
  dpd_file2_init(&F, CC_OEI, 0, 1, 1, "h(a,b)");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  fprintf(outfile, "\tOne-electron energy        = %20.15f\n", one_energy);
  fflush(outfile);

  /* G(Ij,Kl) <-- 1/2 G(IJ,KL) + 1/2 G(ij,kl) + 1/2 G(Ij,Kl) + 1/2 G(iJ,kL) */
  dpd_buf4_init(&G1, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  dpd_buf4_sort(&G1, CC_TMP0, qprs, 0, 0, "GjIKl");
  dpd_buf4_init(&G2, CC_TMP0, 0, 0, 0, 0, 0, 0, "GjIKl");
  dpd_buf4_sort(&G2, CC_TMP0, pqsr, 0, 0, "GiJkL");
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_TMP0, 0, 0, 0, 0, 0, 0, "GiJkL");
  dpd_buf4_axpy(&G2, &G1, 1.0);
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_GAMMA, 0, 0, 0, 2, 2, 0, "GIJKL");
  dpd_buf4_axpy(&G2, &G1, 1.0);
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_GAMMA, 0, 0, 0, 2, 2, 0, "Gijkl");
  dpd_buf4_axpy(&G2, &G1, 1.0);
  dpd_buf4_close(&G2);
  dpd_buf4_scm(&G1, 0.5);
  
  dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  two_energy = dpd_buf4_dot(&A, &G1);
  dpd_buf4_close(&A);
  fprintf(outfile, "\tIJKL energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  dpd_buf4_close(&G1);

  /* G(IJ,KA) <-- G(IJ,KA) + G(ij,ka) + G(Ij,Ka) + G(iJ,kA) */
  dpd_buf4_init(&G1, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_init(&G2, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
  dpd_buf4_axpy(&G2, &G1, 1.0);
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_GAMMA, 0, 0, 10, 2, 10, 0, "GIJKA");
  dpd_buf4_axpy(&G2, &G1, 1.0);
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_GAMMA, 0, 0, 10, 2, 10, 0, "Gijka");
  dpd_buf4_axpy(&G2, &G1, 1.0);
  dpd_buf4_close(&G2);

  dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  two_energy = dpd_buf4_dot(&E, &G1);
  dpd_buf4_close(&E);
  fprintf(outfile, "\tIJKA energy                = %20.15f\n", 2*two_energy);
  total_two_energy += 2*two_energy;

  dpd_buf4_close(&G1);

  /* G(Ij,Ab) <-- G(Ij,Ab) - G(Ib,jA) */
  dpd_buf4_init(&G1, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_init(&G2, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
  dpd_buf4_sort(&G2, CC_TMP0, prsq, 0, 5, "GIbjA (Ij,Ab)");
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_TMP0, 0, 0, 5, 0, 5, 0, "GIbjA (Ij,Ab)");
  dpd_buf4_axpy(&G2, &G1, -1.0);
  dpd_buf4_close(&G2);
  dpd_buf4_close(&G1);
  
  /* G(IJ,AB) <-- G(IJ,AB) - G(IB,JA) */
  dpd_buf4_init(&G1, CC_GAMMA, 0, 0, 5, 2, 7, 0, "GIJAB");
  dpd_buf4_copy(&G1, CC_TMP0, "G(IJ,AB)");
  dpd_buf4_close(&G1);
  dpd_buf4_init(&G1, CC_TMP0, 0, 0, 5, 0, 5, 0, "G(IJ,AB)");
  dpd_buf4_init(&G2, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
  dpd_buf4_sort(&G2, CC_TMP0, prsq, 0, 5, "GIBJA (IJ,AB)");
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_TMP0, 0, 0, 5, 0, 5, 0, "GIBJA (IJ,AB)");
  dpd_buf4_axpy(&G2, &G1, -1.0); 
  dpd_buf4_close(&G2);
  dpd_buf4_close(&G1);

  /* G(ij,ab) <-- G(ij,ab) - G(ib,ja) */
  dpd_buf4_init(&G1, CC_GAMMA, 0, 0, 5, 2, 7, 0, "Gijab");
  dpd_buf4_copy(&G1, CC_TMP0, "G(ij,ab)");
  dpd_buf4_close(&G1);
  dpd_buf4_init(&G1, CC_TMP0, 0, 0, 5, 0, 5, 0, "G(ij,ab)");
  dpd_buf4_init(&G2, CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
  dpd_buf4_sort(&G2, CC_TMP0, prsq, 0, 5, "Gibja (ij,ab)");
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_TMP0, 0, 0, 5, 0, 5, 0, "Gibja (ij,ab)");
  dpd_buf4_axpy(&G2, &G1, -1.0);
  dpd_buf4_close(&G2);
  dpd_buf4_close(&G1);

  /* G(IJ,AB) <-- G(IJ,AB) + G(ij,ab) + G(Ij,Ab) + G(iJ,aB) */
  dpd_buf4_init(&G1, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_sort(&G1, CC_TMP0, qpsr, 0, 5, "G(jI,bA)");
  dpd_buf4_init(&G2, CC_TMP0, 0, 0, 5, 0, 5, 0, "G(jI,bA)");
  dpd_buf4_axpy(&G2, &G1, 1.0);
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_TMP0, 0, 0, 5, 0, 5, 0, "G(IJ,AB)");
  dpd_buf4_axpy(&G2, &G1, 1.0);
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_TMP0, 0, 0, 5, 0, 5, 0, "G(ij,ab)");
  dpd_buf4_axpy(&G2, &G1, 1.0);
  dpd_buf4_close(&G2);

  dpd_buf4_init(&DInts, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  two_energy = dpd_buf4_dot(&DInts, &G1);
  dpd_buf4_close(&DInts);
  fprintf(outfile, "\tIJAB energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  dpd_buf4_close(&G1);
  
  /* G(IB,JA) <-- G(IB,JA) + G(ib,ja) + G(Ib,Ja) + G(iB,jA) */
  dpd_buf4_init(&G1, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
  dpd_buf4_init(&G2, CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
  dpd_buf4_axpy(&G2, &G1, 1.0);
  dpd_buf4_close(&G2); 
  dpd_buf4_init(&G2, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
  dpd_buf4_axpy(&G2, &G1, 1.0); 
  dpd_buf4_close(&G2); 
  dpd_buf4_init(&G2, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
  dpd_buf4_axpy(&G2, &G1, 1.0); 
  dpd_buf4_close(&G2);

  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  two_energy = dpd_buf4_dot(&C, &G1);
  dpd_buf4_close(&C);
  fprintf(outfile, "\tIBJA energy                = %20.15f\n", two_energy);
  total_two_energy += two_energy;

  dpd_buf4_close(&G1);

  fprintf(outfile,"\tTotal two-electron energy  = %20.15f\n", total_two_energy);
  fprintf(outfile, "\tMP2 correlation energy    = %20.15f\n",
          one_energy + total_two_energy);
  fprintf(outfile, "\tTotal MP2 energy          = %20.15f\n",
          one_energy + total_two_energy + mo.Escf);
}

void uhf_sf_deanti(void)
{

}

}} /* End namespaces */
