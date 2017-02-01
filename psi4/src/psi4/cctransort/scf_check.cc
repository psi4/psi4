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

#include <vector>
#include "psi4/psifiles.h"

#include "psi4/libdpd/dpd.h"

namespace psi { namespace cctransort {

double scf_check(int reference, Dimension &openpi)
{
  dpdfile2 H;
  dpdbuf4 A;
  double E1A, E1B, E2AA, E2BB, E2AB;

  if(reference == 2) { // UHF/semicanonical
    global_dpd_->file2_init(&H, PSIF_CC_OEI, 0, 0, 0, "h(I,J)");
    E1A = global_dpd_->file2_trace(&H);
    global_dpd_->file2_close(&H);
    global_dpd_->file2_init(&H, PSIF_CC_OEI, 0, 2, 2, "h(i,j)");
    E1B = global_dpd_->file2_trace(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, "IJ", "KL", "IJ", "KL", 1, "A <IJ|KL>");
    E2AA = 0.5 * global_dpd_->buf4_trace(&A);
    global_dpd_->buf4_close(&A);
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, "ij", "kl", "ij", "kl", 1, "A <ij|kl>");
    E2BB = 0.5 * global_dpd_->buf4_trace(&A);
    global_dpd_->buf4_close(&A);
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, "Ij", "Kl", "Ij", "Kl", 0, "A <Ij|Kl>");
    E2AB = global_dpd_->buf4_trace(&A);
    global_dpd_->buf4_close(&A);

    outfile->Printf( "\tOne-electron energy          =  %20.14f\n", E1A + E1B);
    outfile->Printf( "\tTwo-electron (AA) energy     =  %20.14f\n", E2AA);
    outfile->Printf( "\tTwo-electron (BB) energy     =  %20.14f\n", E2BB);
    outfile->Printf( "\tTwo-electron (AB) energy     =  %20.14f\n", E2AB);
    outfile->Printf( "\tTwo-electron energy          =  %20.14f\n", E2AA + E2BB + E2AB);

    return E1A + E1B + E2AA + E2BB + E2AB;
  }
  else if(reference == 1) { // ROHF
    global_dpd_->file2_init(&H, PSIF_CC_OEI, 0, 0, 0, "h(i,j)");
    E1A = global_dpd_->file2_trace(&H);
    E1B = 0.0;
    global_dpd_->file2_mat_init(&H);
    global_dpd_->file2_mat_rd(&H);
    for(int h=0; h < H.params->nirreps; h++)
      for(int i=0; i < (H.params->ppi[h] - openpi[h]); i++)
        E1B += H.matrix[h][i][i];
    global_dpd_->file2_mat_close(&H);
    global_dpd_->file2_close(&H);

    //ROHF is more complicated because of the way we store the integrals
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, "ij", "kl", "ij", "kl", 1, "A <ij|kl>");
    E2AA = 0.5 * global_dpd_->buf4_trace(&A);
    E2BB = 0.0;
    for(int h=0; h < A.params->nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&A, h);
      global_dpd_->buf4_mat_irrep_rd(&A, h);
      for(int Gi=0; Gi < A.params->nirreps; Gi++) {
        int Gj = Gi ^ h;
        for(int i=0; i < (A.params->ppi[Gi] - openpi[Gi]); i++) {
          int I = A.params->poff[Gi] + i;
          for(int j=0; j < (A.params->qpi[Gj] - openpi[Gj]); j++) {
            int J = A.params->qoff[Gj] + j;
            int IJ = A.params->rowidx[I][J];
            E2BB += 0.5 * A.matrix[h][IJ][IJ];
          }
        }
      }
      global_dpd_->buf4_mat_irrep_close(&A, h);
    }
    global_dpd_->buf4_close(&A);
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, "ij", "kl", "ij", "kl", 0, "A <ij|kl>");
    E2AB = 0.0;
    for(int h=0; h < A.params->nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&A, h);
      global_dpd_->buf4_mat_irrep_rd(&A, h);
      for(int Gi=0; Gi < A.params->nirreps; Gi++) {
        int Gj = Gi ^ h;
        for(int i=0; i < A.params->ppi[Gi]; i++) {
          int I = A.params->poff[Gi] + i;
          for(int j=0; j < (A.params->qpi[Gj] - openpi[Gj]); j++) {
            int J = A.params->qoff[Gj] + j;
            int IJ = A.params->rowidx[I][J];
            E2AB += A.matrix[h][IJ][IJ];
          }
        }
      }
      global_dpd_->buf4_mat_irrep_close(&A, h);
    }
    global_dpd_->buf4_close(&A);

    outfile->Printf( "\tOne-electron energy          =  %20.14f\n", E1A+E1B);
    outfile->Printf( "\tTwo-electron (AA) energy     =  %20.14f\n", E2AA);
    outfile->Printf( "\tTwo-electron (BB) energy     =  %20.14f\n", E2BB);
    outfile->Printf( "\tTwo-electron (AB) energy     =  %20.14f\n", E2AB);
    outfile->Printf( "\tTwo-electron energy          =  %20.14f\n", E1A+E1B+E2AA+E2BB+E2AB);

    return E1A + E1B + E2AA + E2BB + E2AB;
  }
  else { // RHF
    global_dpd_->file2_init(&H, PSIF_CC_OEI, 0, 0, 0, "h(i,j)");
    E1A = 2.0 * global_dpd_->file2_trace(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, "ij", "kl", 0, "A 2<ij|kl> - <ij|lk>");
    E2AB = global_dpd_->buf4_trace(&A);
    global_dpd_->buf4_close(&A);

    outfile->Printf( "\tOne-electron energy          =  %20.14f\n", E1A);
    outfile->Printf( "\tTwo-electron energy          =  %20.14f\n", E2AB);

    return E1A + E2AB;
  }
  return 0.0;
}

}} // namespace psi::cctransort
