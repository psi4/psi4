/*! \file
    \ingroup CINTS
    \brief Initialise global variables.
*/
#include<cstdio>
#include<cstdlib>
#include<libint/libint.h>
#include<psifiles.h>
#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

namespace psi { namespace cints {
  
  //! Initialize global variables.
  void init_globals()
  {
    IOUnits.itap30 = 30;
    IOUnits.itap33 = PSIF_SO_TEI;
    IOUnits.itapDSCF = PSIF_DSCF;
    IOUnits.itapG = PSIF_AO_TPDM;
    IOUnits.itapD1ERI_SO = PSIF_SO_D1ERI;
    IOUnits.itapR12 = PSIF_SO_R12;
    IOUnits.itapT1  = PSIF_SO_R12T1;
    IOUnits.itapERI_MO = PSIF_MO_TEI;
    IOUnits.itapR12_MO = PSIF_MO_R12;
    IOUnits.itapR12T2_MO = PSIF_MO_R12T2;
    IOUnits.itapD = PSIF_AO_OPDM;
    
    IOUnits.itapS = PSIF_OEI;
    IOUnits.itapT = PSIF_OEI;
    IOUnits.itapV = PSIF_OEI;
    IOUnits.itapS_AO = PSIF_OEI;
    IOUnits.itapMX_AO = PSIF_OEI;
    IOUnits.itapMY_AO = PSIF_OEI;
    IOUnits.itapMZ_AO = PSIF_OEI;
    IOUnits.itapQXX_AO = PSIF_OEI;
    IOUnits.itapQXY_AO = PSIF_OEI;
    IOUnits.itapQXZ_AO = PSIF_OEI;
    IOUnits.itapQYY_AO = PSIF_OEI;
    IOUnits.itapQYZ_AO = PSIF_OEI;
    IOUnits.itapQZZ_AO = PSIF_OEI;
    IOUnits.itapNablaX_AO = PSIF_OEI;
    IOUnits.itapNablaY_AO = PSIF_OEI;
    IOUnits.itapNablaZ_AO = PSIF_OEI;
    IOUnits.itapOEInt_Misc = PSIF_OEI;
    
    IOUnits.itapdgdB[0] = PSIF_AO_DGDBX;
    IOUnits.itapdgdB[1] = PSIF_AO_DGDBY;
    IOUnits.itapdgdB[2] = PSIF_AO_DGDBZ;
    
    return;
  }
}}
