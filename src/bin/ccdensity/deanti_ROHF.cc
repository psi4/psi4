/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libiwl/iwl.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
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
    fprintf(outfile, "\n\tEnergies re-computed from Mulliken density:\n");
    fprintf(outfile,   "\t-------------------------------------------\n");

    dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 0, 0, "h(i,j)");
    one_energy += dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 0, 0, "h(i,j)");
    one_energy += dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "h(a,b)");
    one_energy += dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "h(a,b)");
    one_energy += dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
    one_energy += dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
    one_energy += dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
    one_energy += dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
    one_energy += dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    fprintf(outfile, "\tOne-electron energy        = %20.15f\n", one_energy);
    fflush(outfile);
  }

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
  
  if(!params.aobasis) {
    dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    two_energy = dpd_buf4_dot(&A, &G1);
    dpd_buf4_close(&A);
    fprintf(outfile, "\tIJKL energy                = %20.15f\n", two_energy);
    total_two_energy += two_energy;
  }

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

  if(!params.aobasis) {
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    two_energy = dpd_buf4_dot(&E, &G1);
    dpd_buf4_close(&E);
    fprintf(outfile, "\tIJKA energy                = %20.15f\n", 2*two_energy);
    total_two_energy += 2*two_energy;
  }

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

  if(!params.aobasis) {
    dpd_buf4_init(&DInts, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    two_energy = dpd_buf4_dot(&DInts, &G1);
    dpd_buf4_close(&DInts);
    fprintf(outfile, "\tIJAB energy                = %20.15f\n", two_energy);
    total_two_energy += two_energy;
  }

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

  if(!params.aobasis) {
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    two_energy = dpd_buf4_dot(&C, &G1);
    dpd_buf4_close(&C);
    fprintf(outfile, "\tIBJA energy                = %20.15f\n", two_energy);
    total_two_energy += two_energy;
  }

  dpd_buf4_close(&G1);

  /* G(CI,AB) <-- G(CI,AB) + G(ci,ab) + G(Ci,Ab) + G(cI,aB) */
  dpd_buf4_init(&G1, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_init(&G2, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_axpy(&G2, &G1, 1.0);
  dpd_buf4_close(&G2); 
  dpd_buf4_init(&G2, CC_GAMMA, 0, 11, 5, 11, 7, 0, "GCIAB");
  dpd_buf4_axpy(&G2, &G1, 1.0);  
  dpd_buf4_close(&G2); 
  dpd_buf4_init(&G2, CC_GAMMA, 0, 11, 5, 11, 7, 0, "Gciab");
  dpd_buf4_axpy(&G2, &G1, 1.0);
  dpd_buf4_close(&G2);

  if(!params.aobasis) {
    dpd_buf4_init(&FInts, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_sort(&FInts, CC_TMP0, qprs, 11, 5, "F(CI,BA)");
    dpd_buf4_close(&FInts);
    dpd_buf4_init(&FInts, CC_TMP0, 0, 11, 5, 11, 5, 0, "F(CI,BA)");
    dpd_buf4_sort(&FInts, CC_TMP1, pqsr, 11, 5, "F(CI,AB)");
    dpd_buf4_close(&FInts);
    dpd_buf4_init(&FInts, CC_TMP1, 0, 11, 5, 11, 5, 0, "F(CI,AB)");
    two_energy = dpd_buf4_dot(&FInts, &G1);
    dpd_buf4_close(&FInts);
    fprintf(outfile, "\tCIAB energy                = %20.15f\n", 2*two_energy);
    total_two_energy += 2*two_energy;
  }

  dpd_buf4_close(&G1);

  /* G(Ab,Cd) <-- 1/2 G(AB,CD) + 1/2 G(ab,cd) + G(Ab,Cd) */
  dpd_buf4_init(&G1, CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
  dpd_buf4_sort(&G1, CC_TMP0, qprs, 5, 5, "G(bA,Cd)");
  dpd_buf4_init(&G2, CC_TMP0, 0, 5, 5, 5, 5, 0, "G(bA,Cd)");
  dpd_buf4_sort(&G2, CC_TMP1, pqsr, 5, 5, "G(aB,cD)");
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_TMP1, 0, 5, 5, 5, 5, 0, "G(aB,cD)");
  dpd_buf4_axpy(&G2, &G1, 1.0);
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_GAMMA, 0, 5, 5, 7, 7, 0, "GABCD");
  dpd_buf4_axpy(&G2, &G1, 1.0);
  dpd_buf4_close(&G2);
  dpd_buf4_init(&G2, CC_GAMMA, 0, 5, 5, 7, 7, 0, "Gabcd");
  dpd_buf4_axpy(&G2, &G1, 1.0);
  dpd_buf4_close(&G2);
  dpd_buf4_scm(&G1, 0.5);

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
