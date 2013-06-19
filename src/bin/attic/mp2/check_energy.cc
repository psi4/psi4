/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

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
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    E_opdm += global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);
  
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fAB");
    E_opdm += global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
    global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    E_tpdm += 2 * global_dpd_->buf4_dot(&G, &I);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&G);
  }
  else if(chk == 2) {
    fprintf(outfile, "\n\tEnergies re-computed from Fock-adjusted MP2 density:\n");
    fprintf(outfile,   "\t----------------------------------------------------\n");
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "h(i,j)");
    E_opdm += global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "h(a,b)");
    E_opdm += global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "DAI");
    global_dpd_->file2_mat_init(&D);
    global_dpd_->file2_mat_rd(&D);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");
    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);
    for(h=0; h < mo.nirreps; h++)
      for(a=0; a < mo.virtpi[h]; a++)
        for(i=0; i < mo.doccpi[h]; i++)  
          E_opdm += 2*D.matrix[h][a][i]*F.matrix[h][i][a];
    global_dpd_->file2_mat_close(&F);
    global_dpd_->file2_mat_close(&D);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
    global_dpd_->buf4_init(&I, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    global_dpd_->buf4_scmcopy(&I, PSIF_CC_AINTS, "A 2<ij|kl> - <ij|lk>", 2);
    global_dpd_->buf4_sort_axpy(&I, PSIF_CC_AINTS, pqsr, 0, 0, "A 2<ij|kl> - <ij|lk>", -1);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A 2<ij|kl> - <ij|lk>");
    E_tpdm += 0.5 * global_dpd_->buf4_dot(&I, &G);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 0, 11, 0, 0, "GAiJk");
    global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E 2<ai|jk> - <ai|kj>");
    E_tpdm += 2*global_dpd_->buf4_dot(&I, &G);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&G);
    
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
    global_dpd_->buf4_init(&I, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    E_tpdm += global_dpd_->buf4_dot(&I, &G);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    E_tpdm += global_dpd_->buf4_dot(&I, &G);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
    global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    E_tpdm += 2*global_dpd_->buf4_dot(&G, &I);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&G);
  }
  else if(chk == 3) {
    fprintf(outfile, "\n\tEnergies re-computed from MP2 Mulliken density:\n");
    fprintf(outfile,   "\t-----------------------------------------------\n");
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "h(i,j)");
    E_opdm += global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "h(a,b)");
    E_opdm += global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "DAI");
    global_dpd_->file2_mat_init(&D);
    global_dpd_->file2_mat_rd(&D);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");
    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);
    for(h=0; h < mo.nirreps; h++)
      for(a=0; a < mo.virtpi[h]; a++)
        for(i=0; i < mo.doccpi[h]; i++)  
          E_opdm += 2 * D.matrix[h][a][i] * F.matrix[h][i][a];
    global_dpd_->file2_mat_close(&F);
    global_dpd_->file2_mat_close(&D);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
    global_dpd_->buf4_init(&I, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    E_tpdm += 0.5 * global_dpd_->buf4_dot(&I, &G);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&G);
    
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 0, 11, 0, 0, "GAiJk");
    global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    E_tpdm += 2*global_dpd_->buf4_dot(&I, &G);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
    global_dpd_->buf4_init(&I, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    E_tpdm += global_dpd_->buf4_dot(&I, &G);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&G);
    
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
    global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    E_tpdm += 2*global_dpd_->buf4_dot(&G, &I);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&G);
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
