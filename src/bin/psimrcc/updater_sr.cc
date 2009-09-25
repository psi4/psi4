///**
// *  @file updater_sr.cc
// *  @ingroup (PSIMRCC)
// *  @brief Contains methods for updating the CC equations
//*/
//
//#include <string>
////#include <psifiles.h>
//  //#include <libmoinfo/libmoinfo.h>
//  //#include <liboptions/liboptions.hpp>
//  //#include <libpsio/psio.hpp>
//  //#include <libutil/libutil.h>
//  //
//
//  //#include "debugging.h"
//  //#include "index.h"
//  //#include "mrcc.h"
//  //#include "matrix.h"
//
//#include "blas.h"
//#include "updater.h"
//
//namespace psi{ namespace psimrcc{
//
//void CCMRCC::update_amps()
//{
//  update_t1_amps();
//  update_t2_amps();
//  if(triples_type > ccsd_t)
//    update_t3_amps();
//  synchronize_amps();
//
//  // Compute the T-AMPS difference
//  delta_t1_amps=0.0;
//  delta_t2_amps=0.0;
//  for(int n=0;n<moinfo->get_nunique();n++){
//    int m = moinfo->get_ref_number(n,UniqueRefs);
//    delta_t1_amps+=blas->get_scalar("||Delta_t1||",m);
//    delta_t2_amps+=blas->get_scalar("||Delta_t2||",m);
//  }
//  delta_t1_amps=pow(delta_t1_amps,0.5)/((double)moinfo->get_nunique());
//  delta_t2_amps=pow(delta_t2_amps,0.5)/((double)moinfo->get_nunique());
//}
//
//void CCMRCC::update_t1_amps()
//{
//  blas->solve("t1_delta[o][v]{u}  =   t1_eqns[o][v]{u} / d1[o][v]{u} - t1[o][v]{u}");
//  blas->solve("t1_delta[O][V]{u}  =   t1_eqns[O][V]{u} / d1[O][V]{u} - t1[O][V]{u}");
//
//  blas->solve("t1[o][v]{u} = t1_eqns[o][v]{u} / d1[o][v]{u}");
//  blas->solve("t1[O][V]{u} = t1_eqns[O][V]{u} / d1[O][V]{u}");
//
//  blas->solve("t1_norm{u}  = t1[o][v]{u} . t1[o][v]{u}");
//  blas->solve("t1_norm{u} += t1[O][V]{u} . t1[O][V]{u}");
//
//  blas->solve("||Delta_t1||{u}  = t1_delta[o][v]{u} . t1_delta[o][v]{u}");
//  blas->solve("||Delta_t1||{u} += t1_delta[O][V]{u} . t1_delta[O][V]{u}");
//
//  DEBUGGING(3,
//    blas->print("t1[o][v]{u}");
//    blas->print("t1[O][V]{u}");
//  )
//}
//
//void CCMRCC::update_t2_amps()
//{
//  blas->solve("t2_delta[oo][vv]{u} = t2_eqns[oo][vv]{u} / d2[oo][vv]{u} - t2[oo][vv]{u}");
//  blas->solve("t2_delta[oO][vV]{u} = t2_eqns[oO][vV]{u} / d2[oO][vV]{u} - t2[oO][vV]{u}");
//  blas->solve("t2_delta[OO][VV]{u} = t2_eqns[OO][VV]{u} / d2[OO][VV]{u} - t2[OO][VV]{u}");
//
//  blas->solve("t2[oo][vv]{u} = t2_eqns[oo][vv]{u} / d2[oo][vv]{u}");
//  blas->solve("t2[oO][vV]{u} = t2_eqns[oO][vV]{u} / d2[oO][vV]{u}");
//  blas->solve("t2[OO][VV]{u} = t2_eqns[OO][VV]{u} / d2[OO][VV]{u}");
//
//  blas->solve("||Delta_t2||{u}  = t2_delta[oo][vv]{u} . t2_delta[oo][vv]{u}");
//  blas->solve("||Delta_t2||{u} += t2_delta[oO][vV]{u} . t2_delta[oO][vV]{u}");
//  blas->solve("||Delta_t2||{u} += t2_delta[OO][VV]{u} . t2_delta[OO][VV]{u}");
//
//  DEBUGGING(3,
//    blas->print("t2[oo][vv]{u}");
//    blas->print("t2[oO][vV]{u}");
//    blas->print("t2[OO][VV]{u}");
//  );
//}
//
//}}
