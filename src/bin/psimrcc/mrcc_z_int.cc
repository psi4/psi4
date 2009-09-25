/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/
#include "mrcc.h"
#include "matrix.h"
#include "blas.h"
#include <libutil/libutil.h>

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void CCMRCC::build_Z_intermediates()
{
  // I am rewriting
  //   blas->append("Z_ijam[oov][o]{u} = #1234#  1/2 tau[oo][vv]{u} 2@2 <[vo]:[vv]>");
  // as
  blas->append("Z_ijam[oov][o]{u} = #1234#   tau[oo][vv]{u} 2@2 <[vo]|[vv]>");
  //

  blas->append("Z_iJaM[oOv][O]{u} = #1234#   tau[oO][vV]{u} 2@2 <[vo]|[vv]>");

  // I am rewriting
  //   blas->append("Z_iJAm[oOV][o]{u} = #1243# - tau[oO][vV]{u} 2@2 <[ov]|[vv]>");
  // as
  blas->append("Z_iJAm[oOV][o]{u} = #1234# - tau[oO][Vv]{u} 2@2 <[vo]|[vv]>");

  // I am rewriting
  // blas->append("Z_IJAM[OOV][O]{u} = #1234#  1/2 tau[OO][VV]{u} 2@2 <[vo]:[vv]>");
  // as
  blas->append("Z_IJAM[OOV][O]{u} = #1234#   tau[OO][VV]{u} 2@2 <[vo]|[vv]>");
  //
}

}} /* End Namespaces */
