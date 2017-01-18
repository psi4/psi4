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

/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/
#include "psi4/libmoinfo/libmoinfo.h"
#include "mrcc.h"
#include "matrix.h"
#include "blas.h"
#include "debugging.h"
#include "psi4/libpsi4util/libpsi4util.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void CCMRCC::build_W_intermediates()
{
  build_W_mnij_intermediates();
  build_W_mNiJ_intermediates();
  build_W_MNIJ_intermediates();

  build_W_jbme_intermediates();
  build_W_JBme_intermediates();
  build_W_jBmE_intermediates();
  build_W_jbME_intermediates();
  build_W_JbMe_intermediates();
  build_W_JBME_intermediates();
}

void CCMRCC::build_W_mnij_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the W_mnij Intermediates ...");

  );

  blas->append("W_mnij[oo][oo]{u}  = <[oo]:[oo]>");
  blas->append("W_mnij[oo][oo]{u} += #1234# <[ooo]:[v]> 2@2 t1[o][v]{u}");
  blas->append("W_mnij[oo][oo]{u} += #1243# - <[ooo]:[v]> 2@2 t1[o][v]{u}");
  blas->append("W_mnij[oo][oo]{u} += 1/2 <[oo]:[vv]> 2@2 tau[oo][vv]{u}");

  DEBUGGING(3,blas->print("W_mnij[oo][oo]{u}"););

  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}

void CCMRCC::build_W_mNiJ_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the W_mNiJ Intermediates ...");

  );

  blas->append("W_mNiJ[oO][oO]{u}  = <[oo]|[oo]>");
  blas->append("W_mNiJ[oO][oO]{u} += #1234# <[ooo]|[v]> 2@2 t1[O][V]{u}");
  blas->append("W_mNiJ[oO][oO]{u} += #2143# <[ooo]|[v]> 2@2 t1[o][v]{u}");
  blas->append("W_mNiJ[oO][oO]{u} += <[oo]|[vv]> 2@2 tau[oO][vV]{u}");

  DEBUGGING(3,blas->print("W_mNiJ[oO][oO]{u}"););

  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}

void CCMRCC::build_W_MNIJ_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the W_MNIJ Intermediates ...");

  );
  blas->append("W_MNIJ[OO][OO]{u}  = <[oo]:[oo]>");
  blas->append("W_MNIJ[OO][OO]{u} += #1234# <[ooo]:[v]> 2@2 t1[O][V]{u}");
  blas->append("W_MNIJ[OO][OO]{u} += #1243# - <[ooo]:[v]> 2@2 t1[O][V]{u}");
  blas->append("W_MNIJ[OO][OO]{u} += 1/2 <[oo]:[vv]> 2@2 tau[OO][VV]{u}");

  DEBUGGING(3,blas->print("W_MNIJ[OO][OO]{u}"););

  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}

void CCMRCC::build_W_jbme_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the W_jbme Intermediates ...");

  );

  blas->append("W_jbme[ov][ov]{u}  = #3241# <[ov]:[vo]>");

  // This term uses an extra integral file
  // blas->append("W_jbme[ov][ov]{u} += #3241# <[ovv]:[v]> 2@2 t1[o][v]{u}");
  // I will rewrite it as two terms:
  blas->append("W_jbme[ov][ov]{u} += #3241#   <[v]|[ovv]> 1@2 t1[o][v]{u}");
  blas->append("W_jbme[ov][ov]{u} += #2431# - ([vvo]|[v]) 2@2 t1[o][v]{u}");
  //

  blas->append("W_jbme[ov][ov]{u} += #2314# - t1[o][v]{u} 1@1 <[o]:[oov]>");
  blas->append("W_jbme[ov][ov]{u} += - tau3[ov][ov]{u} 2@2 ([ov]:[ov])");
  blas->append("W_jbme[ov][ov]{u} += 1/2 t2[ov][OV]{u} 2@2 ([ov]|[ov])");

  DEBUGGING(3,blas->print("W_jbme[ov][ov]{u}"););

  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}


void CCMRCC::build_W_JBme_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the W_JBme Intermediates ...");

  );
  // Open-Shell
  blas->append("W_JBme[OV][ov]{o}  = #3241# <[ov]|[vo]>");
  //blas->append("W_JBme[OV][ov]{o} += #3241# <[ovv]|[v]> 2@2 t1[O][V]{o}");
  blas->append("W_JBme[OV][ov]{o} += #3241# <[v]|[ovv]> 1@2 t1[O][V]{o}");
  blas->append("W_JBme[OV][ov]{o} += #2314# - t1[O][V]{o} 1@1 <[o]|[oov]>");
  blas->append("W_JBme[OV][ov]{o} += - tau3[OV][OV]{o} 2@2 ([ov]|[ov])");
  blas->append("W_JBme[OV][ov]{o} += 1/2 t2[ov][OV]{o} 1@2 ([ov]:[ov])");

  DEBUGGING(3,blas->print("W_JBme[OV][ov]{o}"););

  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}

void CCMRCC::build_W_jBmE_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the W_jBmE Intermediates ...");

  );

  blas->append("W_jBmE[oV][oV]{u}  = #3214# - <[ov]|[ov]>");
  blas->append("W_jBmE[oV][oV]{u} += #2431# - ([vvo]|[v]) 2@2 t1[o][v]{u}");
  blas->append("W_jBmE[oV][oV]{u} += #2341#   t1[O][V]{u} 1@1 <[o]|[ovo]>");
  blas->append("W_jBmE[oV][oV]{u} += tau3[oV][vO]{u} 2@2 <[ov]|[vo]>");

  DEBUGGING(3,
    blas->print("W_jBmE[oV][oV]{u}");
  );

  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}

void CCMRCC::build_W_jbME_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the W_jbME Intermediates ...");

  );

  blas->append("W_jbME[ov][OV]{u}  = #3241# <[ov]|[vo]>");
//  blas->append("W_jbME[ov][OV]{u} += #3241# <[ovv]|[v]> 2@2 t1[o][v]{u}");
  blas->append("W_jbME[ov][OV]{u} += #3241# <[v]|[ovv]> 1@2 t1[o][v]{u}");
  blas->append("W_jbME[ov][OV]{u} += #2314# - t1[o][v]{u} 1@1 <[o]|[oov]>");
  blas->append("W_jbME[ov][OV]{u} += - tau3[ov][ov]{u} 2@2 ([ov]|[ov])");
  blas->append("W_jbME[ov][OV]{u} += 1/2 t2[ov][OV]{u} 2@2 ([ov]:[ov])");

  DEBUGGING(3,
    blas->print("W_jbME[ov][OV]{u}");
  );
  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}


void CCMRCC::build_W_JbMe_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the W_JbMe Intermediates ...");

  );
  // Open-shell
  blas->append("W_JbMe[Ov][Ov]{o}  = #3214# - <[ov]|[ov]>");
  blas->append("W_JbMe[Ov][Ov]{o} += #2431# - ([vvo]|[v]) 2@2 t1[O][V]{o}");
  blas->append("W_JbMe[Ov][Ov]{o} += #2341#   t1[o][v]{o} 1@1 <[o]|[ovo]>");
  blas->append("W_JbMe[Ov][Ov]{o} += tau3[Ov][Vo]{o} 2@2 <[ov]|[vo]>");

  DEBUGGING(3,
    blas->print("W_JbMe[Ov][Ov]{o}");
  );
  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}

void CCMRCC::build_W_JBME_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the W_jbme Intermediates ...");

  );

  blas->append("W_JBME[OV][OV]{o}  = #3241# <[ov]:[vo]>");

  // This term uses an extra integral file:
  // blas->append("W_JBME[OV][OV]{o} += #3241# <[ovv]:[v]> 2@2 t1[O][V]{o}");
  // I will rewrite it as two terms:
  blas->append("W_JBME[OV][OV]{o} += #3241#   <[v]|[ovv]> 1@2 t1[O][V]{o}");
  blas->append("W_JBME[OV][OV]{o} += #2431# - ([vvo]|[v]) 2@2 t1[O][V]{o}");
  //

  blas->append("W_JBME[OV][OV]{o} += #2314# - t1[O][V]{o} 1@1 <[o]:[oov]>");
  blas->append("W_JBME[OV][OV]{o} += - tau3[OV][OV]{o} 2@2 ([ov]:[ov])");
  blas->append("W_JBME[OV][OV]{o} += 1/2 t2[ov][OV]{o} 1@2 ([ov]|[ov])");

  DEBUGGING(3,
    blas->print("W_JBME[OV][OV]{o}");
  );
  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}

}} /* End Namespaces */
