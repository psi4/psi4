/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libmoinfo/libmoinfo.h"

#include "blas.h"
#include "mp2_ccsd.h"
#include "matrix.h"

namespace psi {
namespace psimrcc {

void MP2_CCSD::build_W_intermediates() {
    build_W_mNiJ_intermediates();

    build_W_jbme_intermediates();
    build_W_jBmE_intermediates();
    build_W_jbME_intermediates();
}

void MP2_CCSD::build_W_mNiJ_intermediates() {
    wfn_->blas()->solve("W_mNiJ[oO][oO]{u}  = <[oo]|[oo]>");
    wfn_->blas()->solve("W_mNiJ[oO][oO]{u} += #1234# <[ooo]|[v]> 2@2 t1[O][V]{u}");
    wfn_->blas()->solve("W_mNiJ[oO][oO]{u} += #2143# <[ooo]|[v]> 2@2 t1[o][v]{u}");
    wfn_->blas()->solve("W_mNiJ[oO][oO]{u} += <[oo]|[vv]> 2@2 tau[oO][vV]{u}");

    wfn_->blas()->reduce_spaces("W_mNiJ[oO][aA]{u}", "W_mNiJ[oO][oO]{u}");
    wfn_->blas()->reduce_spaces("W_mNiJ[oO][oA]{u}", "W_mNiJ[oO][oO]{u}");
}

void MP2_CCSD::build_W_jbme_intermediates() {
    wfn_->blas()->solve("W_jbme[aa][ov]{u}  = #3241# <[oa]:[va]>");
    wfn_->blas()->solve("W_jbme[aa][ov]{u} += #3241# <[oav]:[v]> 2@2 t1_ov[a][v]{u}");
    wfn_->blas()->solve("W_jbme[aa][ov]{u} += #2314# - t1_ov[o][a]{u} 1@1 <[o]:[oav]>");
    wfn_->blas()->solve("W_jbme[aa][ov]{u} += - tau3_ovov[aa][ov]{u} 2@2 ([ov]:[ov])");
    wfn_->blas()->solve("W_jbme[aa][ov]{u} += 1/2 t2_ovOV[aa][OV]{u} 2@2 ([ov]|[ov])");

    wfn_->blas()->solve("W_jbme[oa][ov]{u}  = #3241# <[oa]:[vo]>");
    wfn_->blas()->solve("W_jbme[oa][ov]{u} += #3241# <[oav]:[v]> 2@2 t1[o][v]{u}");
    wfn_->blas()->solve("W_jbme[oa][ov]{u} += #2314# - t1_ov[o][a]{u} 1@1 <[o]:[oov]>");
    wfn_->blas()->solve("W_jbme[oa][ov]{u} += - tau3_ovov[oa][ov]{u} 2@2 ([ov]:[ov])");
    wfn_->blas()->solve("W_jbme[oa][ov]{u} += 1/2 t2_ovOV[oa][OV]{u} 2@2 ([ov]|[ov])");

    wfn_->blas()->solve("W_jbme[av][ov]{u}  = #3241# <[ov]:[va]>");
    wfn_->blas()->solve("W_jbme[av][ov]{u} += #3241# <[ovv]:[v]> 2@2 t1_ov[a][v]{u}");
    wfn_->blas()->solve("W_jbme[av][ov]{u} += #2314# - t1[o][v]{u} 1@1 <[o]:[oav]>");
    wfn_->blas()->solve("W_jbme[av][ov]{u} += - tau3_ovov[av][ov]{u} 2@2 ([ov]:[ov])");
    wfn_->blas()->solve("W_jbme[av][ov]{u} += 1/2 t2_ovOV[av][OV]{u} 2@2 ([ov]|[ov])");

    // This term uses an extra integral file
    // I will rewrite it as two terms:
    /*  wfn_->blas()->solve("W_jbme[ov][ov]{u} += #3241#   <[v]|[ovv]> 1@2 t1[o][v]{u}");
      wfn_->blas()->solve("W_jbme[ov][ov]{u} += #2431# - ([vvo]|[v]) 2@2 t1[o][v]{u}");*/
    //
}

void MP2_CCSD::build_W_jBmE_intermediates() {
    wfn_->blas()->solve("W_jBmE[aA][oV]{u}  = #3214# - <[oa]|[av]>");
    wfn_->blas()->solve("W_jBmE[aA][oV]{u} += #2431# - ([avo]|[v]) 2@2 t1_ov[a][v]{u}");
    wfn_->blas()->solve("W_jBmE[aA][oV]{u} += #2341#   t1_OV[O][A]{u} 1@1 <[o]|[ova]>");
    wfn_->blas()->solve("W_jBmE[aA][oV]{u} += tau3_oVvO[aA][vO]{u} 2@2 <[ov]|[vo]>");

    wfn_->blas()->solve("W_jBmE[oA][oV]{u}  = #3214# - <[oa]|[ov]>");
    wfn_->blas()->solve("W_jBmE[oA][oV]{u} += #2431# - ([avo]|[v]) 2@2 t1[o][v]{u}");
    wfn_->blas()->solve("W_jBmE[oA][oV]{u} += #2341#   t1_OV[O][A]{u} 1@1 <[o]|[ovo]>");
    wfn_->blas()->solve("W_jBmE[oA][oV]{u} += tau3_oVvO[oA][vO]{u} 2@2 <[ov]|[vo]>");

    wfn_->blas()->solve("W_jBmE[aV][oV]{u}  = #3214# - <[ov]|[av]>");
    wfn_->blas()->solve("W_jBmE[aV][oV]{u} += #2431# - ([vvo]|[v]) 2@2 t1_ov[a][v]{u}");
    wfn_->blas()->solve("W_jBmE[aV][oV]{u} += #2341#   t1[O][V]{u} 1@1 <[o]|[ova]>");
    wfn_->blas()->solve("W_jBmE[aV][oV]{u} += tau3_oVvO[aV][vO]{u} 2@2 <[ov]|[vo]>");
}

void MP2_CCSD::build_W_jbME_intermediates() {
    wfn_->blas()->solve("W_jbME[aa][OV]{u}  = #3241# <[oa]|[va]>");
    wfn_->blas()->solve("W_jbME[aa][OV]{u} += #3241# <[v]|[oav]> 1@2 t1_ov[a][v]{u}");
    wfn_->blas()->solve("W_jbME[aa][OV]{u} += #2314# - t1_ov[o][a]{u} 1@1 <[o]|[oav]>");
    wfn_->blas()->solve("W_jbME[aa][OV]{u} += - tau3_ovov[aa][ov]{u} 2@2 ([ov]|[ov])");
    wfn_->blas()->solve("W_jbME[aa][OV]{u} += 1/2 t2_ovOV[aa][OV]{u} 2@2 ([ov]:[ov])");

    wfn_->blas()->solve("W_jbME[oa][OV]{u}  = #3241# <[oa]|[vo]>");
    wfn_->blas()->solve("W_jbME[oa][OV]{u} += #3241# <[v]|[oav]> 1@2 t1[o][v]{u}");
    wfn_->blas()->solve("W_jbME[oa][OV]{u} += #2314# - t1_ov[o][a]{u} 1@1 <[o]|[oov]>");
    wfn_->blas()->solve("W_jbME[oa][OV]{u} += - tau3_ovov[oa][ov]{u} 2@2 ([ov]|[ov])");
    wfn_->blas()->solve("W_jbME[oa][OV]{u} += 1/2 t2_ovOV[oa][OV]{u} 2@2 ([ov]:[ov])");

    wfn_->blas()->solve("W_jbME[av][OV]{u}  = #3241# <[ov]|[va]>");
    wfn_->blas()->solve("W_jbME[av][OV]{u} += #3241# <[v]|[ovv]> 1@2 t1_ov[a][v]{u}");
    wfn_->blas()->solve("W_jbME[av][OV]{u} += #2314# - t1[o][v]{u} 1@1 <[o]|[oav]>");
    wfn_->blas()->solve("W_jbME[av][OV]{u} += - tau3_ovov[av][ov]{u} 2@2 ([ov]|[ov])");
    wfn_->blas()->solve("W_jbME[av][OV]{u} += 1/2 t2_ovOV[av][OV]{u} 2@2 ([ov]:[ov])");
}

}  // namespace psimrcc
}  // namespace psi
