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

/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

void scf_check_uhf(void);
void scf_check_rhf(void);

void scf_check(void)
{
  if(params.ref == 2) scf_check_uhf();
  else scf_check_rhf();
}

void scf_check_uhf(void)
{
  int h, nirreps, i, j, I, J, IJ, Gi, Gj;
  int *aoccpi, *boccpi, *aocc_off, *bocc_off;
  double E1A, E1B, E2AA, E2BB, E2AB;
  dpdfile2 hIJ, hij;
  dpdbuf4 A;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi;
  boccpi = moinfo.boccpi;
  aocc_off = moinfo.aocc_off;
  bocc_off = moinfo.bocc_off;

  global_dpd_->file2_init(&hIJ, PSIF_CC_OEI, 0, 0, 0, "h(I,J)");
  global_dpd_->file2_init(&hij, PSIF_CC_OEI, 0, 2, 2, "h(i,j)");
  global_dpd_->file2_mat_init(&hIJ);
  global_dpd_->file2_mat_init(&hij);
  global_dpd_->file2_mat_rd(&hIJ);
  global_dpd_->file2_mat_rd(&hij);

  E1A = E1B = 0.0;
  for(h=0; h < nirreps; h++) {
    for(i=0; i < aoccpi[h]; i++)
      E1A += hIJ.matrix[h][i][i];

    for(i=0; i < boccpi[h]; i++)
      E1B += hij.matrix[h][i][i];
  }

  global_dpd_->file2_mat_close(&hIJ);
  global_dpd_->file2_mat_close(&hij);

  global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 1, "A <IJ|KL>");
  E2AA = 0.0;
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&A, h);
    global_dpd_->buf4_mat_irrep_rd(&A, h);
    for(Gi=0; Gi < nirreps; Gi++) {
      Gj = Gi ^ h;
      for(i=0; i < aoccpi[Gi]; i++) {
	I = aocc_off[Gi] + i;
	for(j=0; j < aoccpi[Gj]; j++) {
	  J = aocc_off[Gj] + j;
	  IJ = A.params->rowidx[I][J];
	  E2AA += 0.5 * A.matrix[h][IJ][IJ];
	}
      }
    }
    global_dpd_->buf4_mat_irrep_close(&A, h);
  }
  global_dpd_->buf4_close(&A);

  global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 10, 10, 10, 10, 1, "A <ij|kl>");
  E2BB = 0.0;
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&A, h);
    global_dpd_->buf4_mat_irrep_rd(&A, h);
    for(Gi=0; Gi < nirreps; Gi++) {
      Gj = Gi ^ h;
      for(i=0; i < boccpi[Gi]; i++) {
	I = bocc_off[Gi] + i;
	for(j=0; j < boccpi[Gj]; j++) {
	  J = bocc_off[Gj] + j;
	  IJ = A.params->rowidx[I][J];
	  E2BB += 0.5 * A.matrix[h][IJ][IJ];
	}
      }
    }
    global_dpd_->buf4_mat_irrep_close(&A, h);
  }
  global_dpd_->buf4_close(&A);

  global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
  E2AB = 0.0;
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&A, h);
    global_dpd_->buf4_mat_irrep_rd(&A, h);
    for(Gi=0; Gi < nirreps; Gi++) {
      Gj = Gi ^ h;
      for(i=0; i < aoccpi[Gi]; i++) {
	I = aocc_off[Gi] + i;
	for(j=0; j < boccpi[Gj]; j++) {
	  J = bocc_off[Gj] + j;
	  IJ = A.params->rowidx[I][J];
	  E2AB += A.matrix[h][IJ][IJ];
	}
      }
    }
    global_dpd_->buf4_mat_irrep_close(&A, h);
  }
  global_dpd_->buf4_close(&A);

  moinfo.eref = E1A+ E1B+ E2AA+ E2BB+ E2AB + moinfo.enuc + moinfo.efzc;

  psi::fprintf(outfile, "\tOne-electron energy          =  %20.14f\n", E1A + E1B);
  psi::fprintf(outfile, "\tTwo-electron (AA) energy     =  %20.14f\n", E2AA);
  psi::fprintf(outfile, "\tTwo-electron (BB) energy     =  %20.14f\n", E2BB);
  psi::fprintf(outfile, "\tTwo-electron (AB) energy     =  %20.14f\n", E2AB);
  psi::fprintf(outfile, "\tTwo-electron energy          =  %20.14f\n", E2AA + E2BB + E2AB);
  psi::fprintf(outfile, "\tFrozen-core energy (transqt) =  %20.14f\n", moinfo.efzc);
  psi::fprintf(outfile, "\tReference energy             =  %20.14f\n", moinfo.eref);
}

void scf_check_rhf(void)
{
  int h, Gi, Gj;
  int i, j, I, J, IJ;
  int nirreps;
  int *occpi, *virtpi;
  int *occ_off, *vir_off;
  int *occ_sym, *vir_sym;
  int *openpi;
  dpdbuf4 AInts_anti, AInts;
  dpdfile2 Hoo;
  double E1A, E1B, E2AA, E2BB, E2AB;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  /* One-electron (frozen-core) contributions */
  global_dpd_->file2_init(&Hoo, PSIF_CC_OEI, 0, 0, 0, "h(i,j)");
  global_dpd_->file2_mat_init(&Hoo);
  global_dpd_->file2_mat_rd(&Hoo);

  E1A = E1B = 0.0;
  for(h=0; h < nirreps; h++) {

      for(i=0; i < occpi[h]; i++) 
              E1A += Hoo.matrix[h][i][i];   

      for(i=0; i < (occpi[h]-openpi[h]); i++) 
              E1B += Hoo.matrix[h][i][i];   
    }

  global_dpd_->file2_mat_close(&Hoo);
  global_dpd_->file2_close(&Hoo);

  /* Two-electron contributions */

  /* Prepare the A integral buffers */
  global_dpd_->buf4_init(&AInts_anti, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
  global_dpd_->buf4_init(&AInts, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");

  E2AA = E2BB = E2AB = 0.0;
  for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&AInts_anti, h);
      global_dpd_->buf4_mat_irrep_rd(&AInts_anti, h);

      /* Loop over irreps of the target */
      for(Gi=0; Gi < nirreps; Gi++) {
          Gj=Gi^h;

          /* Loop over the orbitals of the target */
          for(i=0; i < occpi[Gi]; i++) {
              I = occ_off[Gi] + i;
              for(j=0; j < occpi[Gj]; j++) {
                  J = occ_off[Gj] + j;

                  IJ = AInts_anti.params->rowidx[I][J];

                  E2AA += AInts_anti.matrix[h][IJ][IJ];

                }
            }

          /* Loop over the orbitals of the target */
          for(i=0; i < (occpi[Gi] - openpi[Gi]); i++) {
              I = occ_off[Gi] + i;
              for(j=0; j < (occpi[Gj] - openpi[Gj]); j++) {
                  J = occ_off[Gj] + j;

                  IJ = AInts_anti.params->rowidx[I][J];

                  E2BB += AInts_anti.matrix[h][IJ][IJ];
                }
            }

        }
      
      global_dpd_->buf4_mat_irrep_close(&AInts_anti, h);

      global_dpd_->buf4_mat_irrep_init(&AInts, h);
      global_dpd_->buf4_mat_irrep_rd(&AInts, h);

      /* Loop over irreps of the target */
      for(Gi=0; Gi < nirreps; Gi++) {
          Gj=Gi^h;

          /* Loop over the orbitals of the target */
          for(i=0; i < occpi[Gi]; i++) {
              I = occ_off[Gi] + i;
              for(j=0; j < (occpi[Gj] - openpi[Gj]); j++) {
                  J = occ_off[Gj] + j;

                  IJ = AInts.params->rowidx[I][J];

                  E2AB += AInts.matrix[h][IJ][IJ];
                }
            }

        }
      
      global_dpd_->buf4_mat_irrep_close(&AInts, h);

    }

  /* Close the A Integral buffers */
  global_dpd_->buf4_close(&AInts_anti);
  global_dpd_->buf4_close(&AInts);

  /*
  psi::fprintf(outfile, "\n\tEFZC = %20.15f\n", moinfo.efzc);
  psi::fprintf(outfile, "\n\tE1A = %20.15f\n", E1A);
  psi::fprintf(outfile,   "\tE1B = %20.15f\n", E1B);
  psi::fprintf(outfile,   "\tE2AA = %20.15f\n", E2AA);
  psi::fprintf(outfile,   "\tE2BB = %20.15f\n", E2BB);
  psi::fprintf(outfile,   "\tE2AB = %20.15f\n", E2AB);
  */

  moinfo.eref = E1A+E1B+0.5*(E2AA+E2BB)+E2AB + moinfo.enuc + moinfo.efzc;

  psi::fprintf(outfile, "\tNuclear energy               =  %20.14f\n", moinfo.enuc);
  psi::fprintf(outfile, "\tOne-electron energy          =  %20.14f\n", E1A+E1B);
  psi::fprintf(outfile, "\tTwo-electron (AA) energy     =  %20.14f\n", E2AA);
  psi::fprintf(outfile, "\tTwo-electron (BB) energy     =  %20.14f\n", E2BB);
  psi::fprintf(outfile, "\tTwo-electron (AB) energy     =  %20.14f\n", E2AB);
  psi::fprintf(outfile, "\tTwo-electron energy          =  %20.14f\n", 0.5*(E2AA+E2BB)+E2AB);
  psi::fprintf(outfile, "\tFrozen-core energy (transqt) =  %20.14f\n", moinfo.efzc);
  psi::fprintf(outfile, "\tReference energy             =  %20.14f\n", moinfo.eref);

}

}} // namespace psi::ccsort
