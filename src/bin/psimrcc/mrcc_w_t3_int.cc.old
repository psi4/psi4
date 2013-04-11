/***************************************************************************
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include "mrcc.h"
#include "blas.h"
#include "debugging.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void CCMRCC::build_W_T3_intermediates()
{
  if(triples_type>ccsd_t){
    build_W_prime_abic_intermediates();
    build_W_prime_aBIc_intermediates();
    build_W_prime_AbiC_intermediates();
    build_W_prime_ABIC_intermediates();

    build_W_prime_ajki_intermediates();
    build_W_prime_AjKi_intermediates();
    build_W_prime_aJkI_intermediates();
    build_W_prime_AJKI_intermediates();

    build_W_kija_intermediates();
    build_W_kiJA_intermediates();
    build_W_KIja_intermediates();
    build_W_KIJA_intermediates();

    build_W_aibc_intermediates();
    build_W_aIbC_intermediates();
    build_W_AiBc_intermediates();
    build_W_AIBC_intermediates();

  }
}

void CCMRCC::build_W_prime_abic_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W'_abic Intermediates ...");
    fflush(outfile);
  );

   blas->solve("W'_abic[vvo][v]{u}  = #4312# <[v]:[ovv]>");

  DEBUGGING(3,blas->print("W'_abic[vvo][v]{u}"););

  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}

void CCMRCC::build_W_prime_aBIc_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W'_aBIc Intermediates ...");
    fflush(outfile);
  );

   blas->solve("W'_aBIc[vVO][v]{u}  = #4312# <[v]|[ovv]>");

  DEBUGGING(3,
    blas->print("W'_aBIc[vVO][v]{u}"););

  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}

void CCMRCC::build_W_prime_AbiC_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W'_AbiC Intermediates ...");
    fflush(outfile);
  );

   blas->solve("W'_AbiC[Vvo][V]{u}  = #4312# <[v]|[ovv]>");

  DEBUGGING(3,
    blas->print("W'_AbiC[Vvo][V]{u}");
  );
  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}

void CCMRCC::build_W_prime_ABIC_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W'_ABIC Intermediates ...");
    fflush(outfile);
  );

   blas->solve("W'_ABIC[VVO][V]{u}  = #4312# <[v]:[ovv]>");

  DEBUGGING(3,
    blas->print("W'_ABIC[VVO][V]{u}");
  );
  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}

void CCMRCC::build_W_prime_ajki_intermediates() // TODO : change indexing
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W'_ajki Intermediates ...");
    fflush(outfile);
  );

  blas->solve("W'_ajki[voo][o]{u}  = #2341# <[oo]:[ov]>");

  DEBUGGING(3,
    blas->print("W'_ajki[voo][o]{u}");
  );
  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}

void CCMRCC::build_W_prime_AjKi_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W'_AjKi Intermediates ...");
    fflush(outfile);
  );

  blas->solve("W'_AjKi[VoO][o]{u}  = #2341# <[oo]|[ov]>");

  DEBUGGING(3,
    blas->print("W'_AjKi[VoO][o]{u}");
  );
  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}

void CCMRCC::build_W_prime_aJkI_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W'_aJkI Intermediates ...");
    fflush(outfile);
  );

  blas->solve("W'_aJkI[vOo][O]{u}  = #2341# <[oo]|[ov]>");

  DEBUGGING(3,
    blas->print("W'_aJkI[vOo][O]{u}");
  );
  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}

void CCMRCC::build_W_prime_AJKI_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W'_AJKI Intermediates ...");
    fflush(outfile);
  );

  blas->solve("W'_AJKI[VOO][O]{u}  = #2341# <[oo]:[ov]>");

  DEBUGGING(3,
    blas->print("W'_AJKI[VOO][O]{u}");
  );
  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}

/**
 * @brief These correspond to the \f$ W_{ijka} \f$ intermediates
 */
void CCMRCC::build_W_kija_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W_kija Intermediates ...");
    fflush(outfile);
  );

  blas->solve("W_kija[o][oov]{u}  = #2314# <[oo]:[ov]>");

  DEBUGGING(3,
    blas->print("W_kija[o][oov]{u}");
  );
  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}

/**
 * @brief These correspond to the \f$ W_{iJkA} \f$ intermediates
 */
void CCMRCC::build_W_kiJA_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W_kiJA Intermediates ...");
    fflush(outfile);
  );

  blas->solve("W_kiJA[o][oOV]{u}  = #2314# <[oo]|[ov]>");

  DEBUGGING(3,
    blas->print("W_kiJA[o][oOV]{u}");
  );
  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}

/**
 * @brief These correspond to the \f$ W_{IjKa} \f$ intermediates
 */
void CCMRCC::build_W_KIja_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W_KIja Intermediates ...");
    fflush(outfile);
  );

  blas->solve("W_KIja[O][Oov]{u}  = #2314# <[oo]|[ov]>");

  DEBUGGING(3,
    blas->print("W_KIja[O][Oov]{u}");
  );
  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}

/**
 * @brief These correspond to the \f$ W_{IJKA} \f$ intermediates
 */
void CCMRCC::build_W_KIJA_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W_KIJA Intermediates ...");
    fflush(outfile);
  );

  blas->solve("W_KIJA[O][OOV]{u}  = #2314# <[oo]:[ov]>");

  DEBUGGING(3,
    blas->print("W_KIJA[O][OOV]{u}");
  );
  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}

/**
 * @brief These correspond to the \f$ W_{aibc} \f$ intermediates
 */
void CCMRCC::build_W_aibc_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W_aibc Intermediates ...");
    fflush(outfile);
  );

  blas->solve("W_aibc[v][ovv]{u}  = <[v]:[ovv]>");

  DEBUGGING(3,
    blas->print("W_aibc[v][ovv]{u}");
  );
  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}

/**
 * @brief These correspond to the \f$ W_{aIbC} \f$ intermediates
 */
void CCMRCC::build_W_aIbC_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W_aIbC Intermediates ...");
    fflush(outfile);
  );

  blas->solve("W_aIbC[v][OvV]{u}  = <[v]|[ovv]>");

  DEBUGGING(3,
    blas->print("W_aIbC[v][OvV]{u}");
  );
  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}

/**
 * @brief These correspond to the \f$ W_{AiBc} \f$ intermediates
 */
void CCMRCC::build_W_AiBc_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W_AiBc Intermediates ...");
    fflush(outfile);
  );

  blas->solve("W_AiBc[V][oVv]{u}  = <[v]|[ovv]>");
  DEBUGGING(3,
    blas->print("W_AiBc[V][oVv]{u}");
  );
  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}


/**
 * @brief These correspond to the \f$ W_{AIBC} \f$ intermediates
 */
void CCMRCC::build_W_AIBC_intermediates()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tBuilding the W_AIBC Intermediates ...");
    fflush(outfile);
  );

  blas->solve("W_AIBC[V][OVV]{u}  = <[v]:[ovv]>");

  DEBUGGING(3,
    blas->print("W_AIBC[V][OVV]{u}");
  );
  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  );
}


}} /* End Namespaces */
