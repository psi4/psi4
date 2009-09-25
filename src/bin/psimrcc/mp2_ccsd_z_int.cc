#include <libutil/libutil.h>
#include <libmoinfo/libmoinfo.h>

#include "blas.h"
#include "mp2_ccsd.h"
#include "debugging.h"
#include "matrix.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void MP2_CCSD::build_Z_intermediates()
{
  blas->solve("Z_iJaM[aAa][O]{u} = #1234#   tau_oOvV[aA][vV]{u} 2@2 <[ao]|[vv]>");
  blas->solve("Z_iJAm[aAA][o]{u} = #1234# - tau_oOVv[aA][Vv]{u} 2@2 <[ao]|[vv]>");

  blas->solve("Z_iJaM[oAa][O]{u} = #1234#   tau_oOvV[oA][vV]{u} 2@2 <[ao]|[vv]>");
  blas->solve("Z_iJAm[oAA][o]{u} = #1234# - tau_oOVv[oA][Vv]{u} 2@2 <[ao]|[vv]>");

  blas->solve("Z_iJaM[aAv][O]{u} = #1234#   tau_oOvV[aA][vV]{u} 2@2 <[vo]|[vv]>");
}

}} /* End Namespace*/
