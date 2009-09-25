/*! 
  \file
  \ingroup MP2
  \brief Check the SCF energy (only RHF working for now)
*/
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void rhf_check_energy(int);
void uhf_check_energy(int);

/*!
** check_energy(): Recompute the SCF energy from the integrals and various
** choices of the density.
**
** \param chk = 1 for MP2 density, 2 for Fock-adjusted MP2 density, 
**   3 = MP2 Mulliken density
** 
** Returns: none
** \ingroup MP2
*/
void check_energy(int chk)
{
  if(params.ref == 0) return(rhf_check_energy(chk));
  else if(params.ref == 2) return(uhf_check_energy(chk));
}

void rhf_check_energy(int chk) 
{
  int h,i,a;
  double E_opdm = 0.0;
  double E_tpdm = 0.0;
  dpdfile2 D, F;
  dpdbuf4 G, I, I2;

  if(chk == 1) {
    fprintf(outfile, "\n\tEnergies re-computed from MP2 density:\n");
    fprintf(outfile,   "\t-------------------------------------\n");
    dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
    dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fIJ");
    E_opdm += dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);
  
    dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fAB");
    E_opdm += dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
    dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    E_tpdm += 2 * dpd_buf4_dot(&G, &I);
    dpd_buf4_close(&I);
    dpd_buf4_close(&G);
  }
  else if(chk == 2) {
    fprintf(outfile, "\n\tEnergies re-computed from Fock-adjusted MP2 density:\n");
    fprintf(outfile,   "\t----------------------------------------------------\n");
    dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
    dpd_file2_init(&F, CC_OEI, 0, 0, 0, "h(i,j)");
    E_opdm += dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "h(a,b)");
    E_opdm += dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 1, 0, "DAI");
    dpd_file2_mat_init(&D);
    dpd_file2_mat_rd(&D);
    dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
    dpd_file2_mat_init(&F);
    dpd_file2_mat_rd(&F);
    for(h=0; h < mo.nirreps; h++)
      for(a=0; a < mo.virtpi[h]; a++)
        for(i=0; i < mo.doccpi[h]; i++)  
          E_opdm += 2*D.matrix[h][a][i]*F.matrix[h][i][a];
    dpd_file2_mat_close(&F);
    dpd_file2_mat_close(&D);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
    dpd_buf4_init(&I, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    dpd_buf4_scmcopy(&I, CC_AINTS, "A 2<ij|kl> - <ij|lk>", 2);
    dpd_buf4_sort_axpy(&I, CC_AINTS, pqsr, 0, 0, "A 2<ij|kl> - <ij|lk>", -1);
    dpd_buf4_close(&I);
    dpd_buf4_init(&I, CC_AINTS, 0, 0, 0, 0, 0, 0, "A 2<ij|kl> - <ij|lk>");
    E_tpdm += 0.5 * dpd_buf4_dot(&I, &G);
    dpd_buf4_close(&I);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, 0, 11, 0, 11, 0, 0, "GAiJk");
    dpd_buf4_init(&I, CC_EINTS, 0, 11, 0, 11, 0, 0, "E 2<ai|jk> - <ai|kj>");
    E_tpdm += 2*dpd_buf4_dot(&I, &G);
    dpd_buf4_close(&I);
    dpd_buf4_close(&G);
    
    dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
    dpd_buf4_init(&I, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    E_tpdm += dpd_buf4_dot(&I, &G);
    dpd_buf4_close(&I);
    dpd_buf4_init(&I, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    E_tpdm += dpd_buf4_dot(&I, &G);
    dpd_buf4_close(&I);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
    dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    E_tpdm += 2*dpd_buf4_dot(&G, &I);
    dpd_buf4_close(&I);
    dpd_buf4_close(&G);
  }
  else if(chk == 3) {
    fprintf(outfile, "\n\tEnergies re-computed from MP2 Mulliken density:\n");
    fprintf(outfile,   "\t-----------------------------------------------\n");
    dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
    dpd_file2_init(&F, CC_OEI, 0, 0, 0, "h(i,j)");
    E_opdm += dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "h(a,b)");
    E_opdm += dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 1, 0, "DAI");
    dpd_file2_mat_init(&D);
    dpd_file2_mat_rd(&D);
    dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
    dpd_file2_mat_init(&F);
    dpd_file2_mat_rd(&F);
    for(h=0; h < mo.nirreps; h++)
      for(a=0; a < mo.virtpi[h]; a++)
        for(i=0; i < mo.doccpi[h]; i++)  
          E_opdm += 2 * D.matrix[h][a][i] * F.matrix[h][i][a];
    dpd_file2_mat_close(&F);
    dpd_file2_mat_close(&D);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
    dpd_buf4_init(&I, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    E_tpdm += 0.5 * dpd_buf4_dot(&I, &G);
    dpd_buf4_close(&I);
    dpd_buf4_close(&G);
    
    dpd_buf4_init(&G, CC_GAMMA, 0, 11, 0, 11, 0, 0, "GAiJk");
    dpd_buf4_init(&I, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    E_tpdm += 2*dpd_buf4_dot(&I, &G);
    dpd_buf4_close(&I);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
    dpd_buf4_init(&I, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    E_tpdm += dpd_buf4_dot(&I, &G);
    dpd_buf4_close(&I);
    dpd_buf4_close(&G);
    
    dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
    dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    E_tpdm += 2*dpd_buf4_dot(&G, &I);
    dpd_buf4_close(&I);
    dpd_buf4_close(&G);
  }
  else {

  }

  fprintf(outfile,"\n");
  fprintf(outfile,"\tE_OPDM                  = %20.15f\n",E_opdm);
  fprintf(outfile,"\tE_TPDM                  = %20.15f\n",E_tpdm);
  fprintf(outfile,"\tMP2 correlation energy  = %20.15f\n",E_opdm+E_tpdm);
  fflush(outfile);
}

void uhf_check_energy(int chk)
{

}

}} /* End namespace */
