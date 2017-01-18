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

#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libmoinfo/libmoinfo.h"

#include "blas.h"
#include "mp2_ccsd.h"
#include "debugging.h"
#include "matrix.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void MP2_CCSD::build_W_intermediates()
{
  build_W_mNiJ_intermediates();

  build_W_jbme_intermediates();
  build_W_jBmE_intermediates();
  build_W_jbME_intermediates();
}

void MP2_CCSD::build_W_mNiJ_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the W_mNiJ Intermediates ...");

  );

  blas->solve("W_mNiJ[oO][oO]{u}  = <[oo]|[oo]>");
  blas->solve("W_mNiJ[oO][oO]{u} += #1234# <[ooo]|[v]> 2@2 t1[O][V]{u}");
  blas->solve("W_mNiJ[oO][oO]{u} += #2143# <[ooo]|[v]> 2@2 t1[o][v]{u}");
  blas->solve("W_mNiJ[oO][oO]{u} += <[oo]|[vv]> 2@2 tau[oO][vV]{u}");

  blas->reduce_spaces("W_mNiJ[oO][aA]{u}","W_mNiJ[oO][oO]{u}");
  blas->reduce_spaces("W_mNiJ[oO][oA]{u}","W_mNiJ[oO][oO]{u}");

  DEBUGGING(3,blas->print("W_mNiJ[oO][aA]{u}"););

  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}

void MP2_CCSD::build_W_jbme_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the W_jbme Intermediates ...");

  );

  blas->solve("W_jbme[aa][ov]{u}  = #3241# <[oa]:[va]>");
  blas->solve("W_jbme[aa][ov]{u} += #3241# <[oav]:[v]> 2@2 t1_ov[a][v]{u}");
  blas->solve("W_jbme[aa][ov]{u} += #2314# - t1_ov[o][a]{u} 1@1 <[o]:[oav]>");
  blas->solve("W_jbme[aa][ov]{u} += - tau3_ovov[aa][ov]{u} 2@2 ([ov]:[ov])");
  blas->solve("W_jbme[aa][ov]{u} += 1/2 t2_ovOV[aa][OV]{u} 2@2 ([ov]|[ov])");

  blas->solve("W_jbme[oa][ov]{u}  = #3241# <[oa]:[vo]>");
  blas->solve("W_jbme[oa][ov]{u} += #3241# <[oav]:[v]> 2@2 t1[o][v]{u}");
  blas->solve("W_jbme[oa][ov]{u} += #2314# - t1_ov[o][a]{u} 1@1 <[o]:[oov]>");
  blas->solve("W_jbme[oa][ov]{u} += - tau3_ovov[oa][ov]{u} 2@2 ([ov]:[ov])");
  blas->solve("W_jbme[oa][ov]{u} += 1/2 t2_ovOV[oa][OV]{u} 2@2 ([ov]|[ov])");

  blas->solve("W_jbme[av][ov]{u}  = #3241# <[ov]:[va]>");
  blas->solve("W_jbme[av][ov]{u} += #3241# <[ovv]:[v]> 2@2 t1_ov[a][v]{u}");
  blas->solve("W_jbme[av][ov]{u} += #2314# - t1[o][v]{u} 1@1 <[o]:[oav]>");
  blas->solve("W_jbme[av][ov]{u} += - tau3_ovov[av][ov]{u} 2@2 ([ov]:[ov])");
  blas->solve("W_jbme[av][ov]{u} += 1/2 t2_ovOV[av][OV]{u} 2@2 ([ov]|[ov])");

  DEBUGGING(3,blas->print("W_jbme[aa][ov]{u}"););

  // This term uses an extra integral file
  // I will rewrite it as two terms:
/*  blas->solve("W_jbme[ov][ov]{u} += #3241#   <[v]|[ovv]> 1@2 t1[o][v]{u}");
  blas->solve("W_jbme[ov][ov]{u} += #2431# - ([vvo]|[v]) 2@2 t1[o][v]{u}");*/
  //

  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}

void MP2_CCSD::build_W_jBmE_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the W_jBmE Intermediates ...");

  );

  blas->solve("W_jBmE[aA][oV]{u}  = #3214# - <[oa]|[av]>");
  blas->solve("W_jBmE[aA][oV]{u} += #2431# - ([avo]|[v]) 2@2 t1_ov[a][v]{u}");
  blas->solve("W_jBmE[aA][oV]{u} += #2341#   t1_OV[O][A]{u} 1@1 <[o]|[ova]>");
  blas->solve("W_jBmE[aA][oV]{u} += tau3_oVvO[aA][vO]{u} 2@2 <[ov]|[vo]>");

  blas->solve("W_jBmE[oA][oV]{u}  = #3214# - <[oa]|[ov]>");
  blas->solve("W_jBmE[oA][oV]{u} += #2431# - ([avo]|[v]) 2@2 t1[o][v]{u}");
  blas->solve("W_jBmE[oA][oV]{u} += #2341#   t1_OV[O][A]{u} 1@1 <[o]|[ovo]>");
  blas->solve("W_jBmE[oA][oV]{u} += tau3_oVvO[oA][vO]{u} 2@2 <[ov]|[vo]>");

  blas->solve("W_jBmE[aV][oV]{u}  = #3214# - <[ov]|[av]>");
  blas->solve("W_jBmE[aV][oV]{u} += #2431# - ([vvo]|[v]) 2@2 t1_ov[a][v]{u}");
  blas->solve("W_jBmE[aV][oV]{u} += #2341#   t1[O][V]{u} 1@1 <[o]|[ova]>");
  blas->solve("W_jBmE[aV][oV]{u} += tau3_oVvO[aV][vO]{u} 2@2 <[ov]|[vo]>");

  DEBUGGING(3,
    blas->print("W_jBmE[oV][oV]{u}");
  );

  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}

void MP2_CCSD::build_W_jbME_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the W_jbME Intermediates ...");

  );

  blas->solve("W_jbME[aa][OV]{u}  = #3241# <[oa]|[va]>");
  blas->solve("W_jbME[aa][OV]{u} += #3241# <[v]|[oav]> 1@2 t1_ov[a][v]{u}");
  blas->solve("W_jbME[aa][OV]{u} += #2314# - t1_ov[o][a]{u} 1@1 <[o]|[oav]>");
  blas->solve("W_jbME[aa][OV]{u} += - tau3_ovov[aa][ov]{u} 2@2 ([ov]|[ov])");
  blas->solve("W_jbME[aa][OV]{u} += 1/2 t2_ovOV[aa][OV]{u} 2@2 ([ov]:[ov])");

  blas->solve("W_jbME[oa][OV]{u}  = #3241# <[oa]|[vo]>");
  blas->solve("W_jbME[oa][OV]{u} += #3241# <[v]|[oav]> 1@2 t1[o][v]{u}");
  blas->solve("W_jbME[oa][OV]{u} += #2314# - t1_ov[o][a]{u} 1@1 <[o]|[oov]>");
  blas->solve("W_jbME[oa][OV]{u} += - tau3_ovov[oa][ov]{u} 2@2 ([ov]|[ov])");
  blas->solve("W_jbME[oa][OV]{u} += 1/2 t2_ovOV[oa][OV]{u} 2@2 ([ov]:[ov])");

  blas->solve("W_jbME[av][OV]{u}  = #3241# <[ov]|[va]>");
  blas->solve("W_jbME[av][OV]{u} += #3241# <[v]|[ovv]> 1@2 t1_ov[a][v]{u}");
  blas->solve("W_jbME[av][OV]{u} += #2314# - t1[o][v]{u} 1@1 <[o]|[oav]>");
  blas->solve("W_jbME[av][OV]{u} += - tau3_ovov[av][ov]{u} 2@2 ([ov]|[ov])");
  blas->solve("W_jbME[av][OV]{u} += 1/2 t2_ovOV[av][OV]{u} 2@2 ([ov]:[ov])");


  DEBUGGING(3,
    blas->print("W_jbME[aa][OV]{u}");
  );
  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}

}} /* End Namespace*/
