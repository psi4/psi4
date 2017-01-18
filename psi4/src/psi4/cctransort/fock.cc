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

#include "psi4/libdpd/dpd.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"

namespace psi { namespace cctransort {

void fock_uhf(std::shared_ptr<Wavefunction> ref, Dimension &aoccpi, Dimension  &boccpi,
              Dimension &avirpi, Dimension &bvirpi, Dimension &frdocc, int print)
{
  dpdfile2 fa, fb;

  SharedMatrix Fa = ref->Fa()->clone();
  SharedMatrix Fb = ref->Fb()->clone();
  SharedMatrix Ca = ref->Ca();
  SharedMatrix Cb = ref->Cb();
  Fa->transform(Ca);
  Fb->transform(Cb);

  global_dpd_->file2_init(&fa, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_init(&fb, PSIF_CC_OEI, 0, 2, 2, "fij");
  global_dpd_->file2_mat_init(&fa);
  global_dpd_->file2_mat_init(&fb);
  for(int h=0; h < fa.params->nirreps; h++) {
    for(int i=0; i < aoccpi[h]; i++)
      for(int j=0; j < aoccpi[h]; j++)
          fa.matrix[h][i][j] = Fa->get(h,i + frdocc[h],j + frdocc[h]);
    for(int i=0; i < boccpi[h]; i++)
      for(int j=0; j < boccpi[h]; j++)
        fb.matrix[h][i][j] = Fb->get(h,i + frdocc[h],j + frdocc[h]);
  }
  global_dpd_->file2_mat_wrt(&fa);
  global_dpd_->file2_mat_wrt(&fb);
  global_dpd_->file2_mat_close(&fa);
  global_dpd_->file2_mat_close(&fb);
  if(print > 3) {
    global_dpd_->file2_print(&fa, "outfile");
    global_dpd_->file2_print(&fb, "outfile");
  }
  global_dpd_->file2_close(&fa);
  global_dpd_->file2_close(&fb);

  global_dpd_->file2_init(&fa, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_init(&fb, PSIF_CC_OEI, 0, 3, 3, "fab");
  global_dpd_->file2_mat_init(&fa);
  global_dpd_->file2_mat_init(&fb);
  for(int h=0; h < fa.params->nirreps; h++) {
    for(int a=0; a < avirpi[h]; a++)
      for(int b=0; b < avirpi[h]; b++)
        fa.matrix[h][a][b] = Fa->get(h, a+frdocc[h]+aoccpi[h], b+frdocc[h]+aoccpi[h]);
    for(int a=0; a < bvirpi[h]; a++)
      for(int b=0; b < bvirpi[h]; b++)
        fb.matrix[h][a][b] = Fb->get(h, a+frdocc[h]+boccpi[h], b+frdocc[h]+boccpi[h]);
  }
  global_dpd_->file2_mat_wrt(&fa);
  global_dpd_->file2_mat_wrt(&fb);
  global_dpd_->file2_mat_close(&fa);
  global_dpd_->file2_mat_close(&fb);
  if(print > 3) {
    global_dpd_->file2_print(&fa, "outfile");
    global_dpd_->file2_print(&fb, "outfile");
  }
  global_dpd_->file2_close(&fa);
  global_dpd_->file2_close(&fb);

  global_dpd_->file2_init(&fa, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_init(&fb, PSIF_CC_OEI, 0, 2, 3, "fia");
  global_dpd_->file2_mat_init(&fa);
  global_dpd_->file2_mat_init(&fb);
  for(int h=0; h < fa.params->nirreps; h++) {
    for(int i=0; i < aoccpi[h]; i++)
      for(int a=0; a < avirpi[h]; a++)
        fa.matrix[h][i][a] = Fa->get(h, i+frdocc[h], a+frdocc[h]+aoccpi[h]);
    for(int i=0; i < boccpi[h]; i++)
      for(int a=0; a < bvirpi[h]; a++)
        fb.matrix[h][i][a] = Fb->get(h, i+frdocc[h], a+frdocc[h]+boccpi[h]);
  }
  global_dpd_->file2_mat_wrt(&fa);
  global_dpd_->file2_mat_wrt(&fb);
  global_dpd_->file2_mat_close(&fa);
  global_dpd_->file2_mat_close(&fb);
  if(print > 3) {
    global_dpd_->file2_print(&fa, "outfile");
    global_dpd_->file2_print(&fb, "outfile");
  }
  global_dpd_->file2_close(&fa);
  global_dpd_->file2_close(&fb);
}

void fock_rhf(std::shared_ptr<Wavefunction> ref, Dimension &occpi, Dimension &openpi,
              Dimension &virpi, Dimension &frdocc, int print)
{
  dpdfile2 fa, fb;

  SharedMatrix Fa = ref->Fa()->clone();
  SharedMatrix Fb = ref->Fb()->clone();
  SharedMatrix Ca = ref->Ca();
  SharedMatrix Cb = ref->Cb();
  Fa->transform(Ca);
  Fb->transform(Cb);

  global_dpd_->file2_init(&fa, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_init(&fb, PSIF_CC_OEI, 0, 0, 0, "fij");
  global_dpd_->file2_mat_init(&fa);
  global_dpd_->file2_mat_init(&fb);
  for(int h=0; h < fa.params->nirreps; h++) {
    for(int i=0; i < occpi[h]; i++)
      for(int j=0; j < occpi[h]; j++)
        fa.matrix[h][i][j] = Fa->get(h, i+frdocc[h], j+frdocc[h]);
    for(int i=0; i < (occpi[h]-openpi[h]); i++)
      for(int j=0; j < (occpi[h]-openpi[h]); j++)
        fb.matrix[h][i][j] = Fb->get(h, i+frdocc[h], j+frdocc[h]);
  }
  global_dpd_->file2_mat_wrt(&fa);
  global_dpd_->file2_mat_wrt(&fb);
  global_dpd_->file2_mat_close(&fa);
  global_dpd_->file2_mat_close(&fb);
  if(print > 3) {
    global_dpd_->file2_print(&fa, "outfile");
    global_dpd_->file2_print(&fb, "outfile");
  }
  global_dpd_->file2_close(&fa);
  global_dpd_->file2_close(&fb);

  global_dpd_->file2_init(&fa, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_init(&fb, PSIF_CC_OEI, 0, 1, 1, "fab");
  global_dpd_->file2_mat_init(&fa);
  global_dpd_->file2_mat_init(&fb);
  for(int h=0; h < fa.params->nirreps; h++) {
    for(int a=0; a < (virpi[h]-openpi[h]); a++)
      for(int b=0; b < (virpi[h]-openpi[h]); b++)
        fa.matrix[h][a][b] = Fa->get(h, a+frdocc[h]+occpi[h], b+frdocc[h]+occpi[h]);

    for(int a=0; a < (virpi[h]-openpi[h]); a++)
      for(int b=0; b < (virpi[h]-openpi[h]); b++)
        fb.matrix[h][a][b] = Fb->get(h, a+frdocc[h]+occpi[h], b+frdocc[h]+occpi[h]);

    for(int a=0; a < openpi[h]; a++)
      for(int b=0; b < openpi[h]; b++)
        fb.matrix[h][virpi[h]-openpi[h]+a][virpi[h]-openpi[h]+b] =
          Fb->get(h, a+frdocc[h]+occpi[h]-openpi[h], b+frdocc[h]+occpi[h]-openpi[h]);

    for(int a=0; a < (virpi[h]-openpi[h]); a++)
      for(int b=0; b < openpi[h]; b++)
        fb.matrix[h][a][virpi[h]-openpi[h]+b] =
          Fb->get(h, a+frdocc[h]+occpi[h], b+frdocc[h]+occpi[h]-openpi[h]);

    for(int a=0; a < openpi[h]; a++)
      for(int b=0; b < (virpi[h]-openpi[h]); b++)
        fb.matrix[h][virpi[h]-openpi[h]+a][b] =
          Fb->get(h, a+frdocc[h]+occpi[h]-openpi[h], b+frdocc[h]+occpi[h]);
  }
  global_dpd_->file2_mat_wrt(&fa);
  global_dpd_->file2_mat_wrt(&fb);
  global_dpd_->file2_mat_close(&fa);
  global_dpd_->file2_mat_close(&fb);
  if(print > 3) {
    global_dpd_->file2_print(&fa, "outfile");
    global_dpd_->file2_print(&fb, "outfile");
  }
  global_dpd_->file2_close(&fa);
  global_dpd_->file2_close(&fb);

  global_dpd_->file2_init(&fa, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_init(&fb, PSIF_CC_OEI, 0, 0, 1, "fia");
  global_dpd_->file2_mat_init(&fa);
  global_dpd_->file2_mat_init(&fb);

  for(int h=0; h < fa.params->nirreps; h++) {
    for(int i=0; i < occpi[h]; i++)
      for(int a=0; a < (virpi[h]-openpi[h]); a++)
        fa.matrix[h][i][a] = Fa->get(h, i+frdocc[h], a+frdocc[h]+occpi[h]);

    for(int i=0; i < (occpi[h]-openpi[h]); i++)
      for(int a=0; a < (virpi[h]-openpi[h]); a++)
        fb.matrix[h][i][a] = Fb->get(h, i+frdocc[h], a+frdocc[h]+occpi[h]);

    for(int i=0; i < (occpi[h]-openpi[h]); i++)
      for(int a=0; a < openpi[h]; a++)
        fb.matrix[h][i][virpi[h]-openpi[h]+a] = Fb->get(h, i+frdocc[h], a+frdocc[h]+occpi[h]-openpi[h]);
  }

  global_dpd_->file2_mat_wrt(&fa);
  global_dpd_->file2_mat_wrt(&fb);
  global_dpd_->file2_mat_close(&fa);
  global_dpd_->file2_mat_close(&fb);
  if(print > 3) {
    global_dpd_->file2_print(&fa, "outfile");
    global_dpd_->file2_print(&fb, "outfile");
  }
  global_dpd_->file2_close(&fa);
  global_dpd_->file2_close(&fb);
}

}} // namespace psi::ccsort
