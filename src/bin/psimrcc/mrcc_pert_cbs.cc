#include "blas.h"
#include "mrcc.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void CCMRCC::compute_first_order_amps()
{
  /*
   * MP2 - Doubles contributions
   */
  fprintf(outfile,"\n\n  Computing first-order amplitudes for frozen-virtual MOs");

  blas->solve("t2_1[oo][vf]{u}   = <[oo]:[vf]> / d2[oo][vf]{u}");
  blas->solve("t2_1[oO][vF]{u}   = <[oo]|[vf]> / d2[oO][vF]{u}");
  blas->solve("t2_1[OO][VF]{u}   = <[oo]:[vf]> / d2[OO][VF]{u}");

  blas->solve("t2_1[oo][fv]{u}   = <[oo]:[fv]> / d2[oo][fv]{u}");
  blas->solve("t2_1[oO][fV]{u}   = <[oo]|[fv]> / d2[oO][fV]{u}");
  blas->solve("t2_1[OO][FV]{u}   = <[oo]:[fv]> / d2[OO][FV]{u}");

  blas->solve("t2_1[oo][ff]{u}   = <[oo]:[ff]> / d2[oo][ff]{u}");
  blas->solve("t2_1[oO][fF]{u}   = <[oo]|[ff]> / d2[oO][fF]{u}");
  blas->solve("t2_1[OO][FF]{u}   = <[oo]:[ff]> / d2[OO][FF]{u}");

  blas->solve("t2_1[ov][of]{u} = #1324# t2_1[oo][vf]{u}");
//  blas->add_Matrix("t2_1[ov][of]{u}");
  blas->solve("t2_1[ov][OF]{u} = #1324# t2_1[oO][vF]{u}");
//  blas->add_Matrix("t2_1[ov][OF]{u}");
  blas->solve("t2_1[oV][Of]{u} = #1342# t2_1[oO][fV]{u}");
//  blas->add_Matrix("t2_1[oV][Of]{u}");
  blas->solve("t2_1[oF][Ov]{u} = #1342# t2_1[oO][vF]{u}");
//  blas->add_Matrix("t2_1[oF][Ov]{u}");
  blas->solve("t2_1[OV][OF]{u} = #1324# t2_1[OO][VF]{u}");
//  blas->add_Matrix("t2_1[OV][OF]{u}");
  blas->solve("t2_1[of][OV]{u} = #1324# t2_1[oO][fV]{u}");
//  blas->add_Matrix("t2_1[OV][of]{u}");

  blas->solve("t2_1[o][ovf]{u} = #1234# t2_1[oo][vf]{u}");
//  blas->add_Matrix("t2_1[o][ovf]{u}");
  blas->solve("t2_1[o][OvF]{u} = #1234# t2_1[oO][vF]{u}");
//  blas->add_Matrix("t2_1[o][OvF]{u}");
  blas->solve("t2_1[o][OfV]{u} = #1234# t2_1[oO][fV]{u}");
//  blas->add_Matrix("t2_1[o][OfV]{u}");
  blas->solve("t2_1[o][off]{u} = #1234# t2_1[oo][ff]{u}");
//  blas->add_Matrix("t2_1[o][off]{u}");
  blas->solve("t2_1[o][OfF]{u} = #1234# t2_1[oO][fF]{u}");
//  blas->add_Matrix("t2_1[o][OfF]{u}");
  blas->solve("t2_1[v][foo]{u} = #3412# t2_1[oo][vf]{u}");
//  blas->add_Matrix("t2_1[v][foo]{u}");
  blas->solve("t2_1[v][FoO]{u} = #3412# t2_1[oO][vF]{u}");
//  blas->add_Matrix("t2_1[v][FoO]{u}");
}

void CCMRCC::perturbative_cbs()
{

  fprintf(outfile,"\n\n  Computing perturbative corrections for frozen-virtual MOs");
  /*
   * MP2 - Doubles contributions
   */
  blas->solve("t2_1[oo][vf]{u}   = <[oo]:[vf]> / d2[oo][vf]{u}");
  blas->solve("t2_1[oO][vF]{u}   = <[oo]|[vf]> / d2[oO][vF]{u}");
  blas->solve("t2_1[OO][VF]{u}   = <[oo]:[vf]> / d2[OO][VF]{u}");

  blas->solve("Eaaaa{u} = 1/4 t2_1[oo][vf]{u} . <[oo]:[vf]>");
  blas->solve("Eabab{u} =     t2_1[oO][vF]{u} . <[oo]|[vf]>");
  blas->solve("Ebbbb{u} = 1/4 t2_1[OO][VF]{u} . <[oo]:[vf]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_vf = blas->get_scalar("ECCSD",0);


  blas->solve("t2_1[oo][fv]{u}   = <[oo]:[fv]> / d2[oo][fv]{u}");
  blas->solve("t2_1[oO][fV]{u}   = <[oo]|[fv]> / d2[oO][fV]{u}");
  blas->solve("t2_1[OO][FV]{u}   = <[oo]:[fv]> / d2[OO][FV]{u}");

  blas->solve("Eaaaa{u} = 1/4 t2_1[oo][fv]{u} . <[oo]:[fv]>");
  blas->solve("Eabab{u} =     t2_1[oO][fV]{u} . <[oo]|[fv]>");
  blas->solve("Ebbbb{u} = 1/4 t2_1[OO][FV]{u} . <[oo]:[fv]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_fv = blas->get_scalar("ECCSD",0);


  blas->solve("t2_1[oo][ff]{u}   = <[oo]:[ff]> / d2[oo][ff]{u}");
  blas->solve("t2_1[oO][fF]{u}   = <[oo]|[ff]> / d2[oO][fF]{u}");
  blas->solve("t2_1[OO][FF]{u}   = <[oo]:[ff]> / d2[OO][FF]{u}");

  blas->solve("Eaaaa{u} = 1/4 t2_1[oo][ff]{u} . <[oo]:[ff]>");
  blas->solve("Eabab{u} =     t2_1[oO][fF]{u} . <[oo]|[ff]>");
  blas->solve("Ebbbb{u} = 1/4 t2_1[OO][FF]{u} . <[oo]:[ff]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_ff = blas->get_scalar("ECCSD",0);


  fprintf(outfile,"\n        CBS second-order correction (vf) = %20.12f",E_vf);
  fprintf(outfile,"\n        CBS second-order correction (fv) = %20.12f",E_fv);
  fprintf(outfile,"\n        CBS second-order correction (ff) = %20.12f",E_ff);
  fprintf(outfile,"\n        CBS second-order correction      = %20.12f",E_vf + E_fv + E_ff);


  /*
   * MP3 - Doubles contributions
   */

  // Sort amplitudes
  blas->solve("t2[ov][ov]{u} = #1324# t2[oo][vv]{u}");
  blas->solve("t2[ov][OV]{u} = #1324# t2[oO][vV]{u}");
  blas->solve("t2[oV][Ov]{u} = #1342# t2[oO][vV]{u}");
  blas->solve("t2[OV][OV]{u} = #1324# t2[OO][VV]{u}");

  // OOVF Terms
  // PPPP Ladder
  blas->solve("t2_eqns[oo][vf]{u} += 1/2 t2[oo][vv]{u} 2@2 <[vf]:[vv]>");

  blas->solve("t2_eqns[oO][vF]{u} +=     t2[oO][vV]{u} 2@2 <[vf]|[vv]>");

  blas->solve("t2_eqns[OO][VF]{u} += 1/2 t2[OO][VV]{u} 2@2 <[vf]:[vv]>");

  // HPHP Ring
  blas->solve("t2_eqns[oo][vf]{u} += #1342#   t2[ov][ov]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[oo][vf]{u} += #2341# - t2[ov][ov]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[oo][vf]{u} += #1342#   t2[ov][OV]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[oo][vf]{u} += #2341# - t2[ov][OV]{u} 2@2 ([fo]|[ov])");

  blas->solve("t2_eqns[oO][vF]{u} += #1342#   t2[ov][ov]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[oO][vF]{u} += #1342#   t2[ov][OV]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[oO][vF]{u} += #2314# - t2[oV][Ov]{u} 1@2 <[of]|[ov]>");

  blas->solve("t2_eqns[OO][VF]{u} += #1342#   t2[OV][OV]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[OO][VF]{u} += #2341# - t2[OV][OV]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[OO][VF]{u} += #1342#   t2[ov][OV]{u} 1@2 ([fo]|[ov])");
  blas->solve("t2_eqns[OO][VF]{u} += #2341# - t2[ov][OV]{u} 1@2 ([fo]|[ov])");

  // Solve for second-order amplitudes
  blas->solve("t2_2[oo][vf]{u}   = t2_eqns[oo][vf]{u} / d2[oo][vf]{u}");
  blas->solve("t2_2[oO][vF]{u}   = t2_eqns[oO][vF]{u} / d2[oO][vF]{u}");
  blas->solve("t2_2[OO][VF]{u}   = t2_eqns[OO][VF]{u} / d2[OO][VF]{u}");

  // Compute third-order energy
  blas->solve("Eaaaa{u} = 1/4 t2_2[oo][vf]{u} . <[oo]:[vf]>");
  blas->solve("Eabab{u} =     t2_2[oO][vF]{u} . <[oo]|[vf]>");
  blas->solve("Ebbbb{u} = 1/4 t2_2[OO][VF]{u} . <[oo]:[vf]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_vf_3 = blas->get_scalar("ECCSD",0);

  // T1 Contributions
  blas->solve("t2_eqns[oo][vf]{u} += #1234#   t1[o][v]{u} 2@1 <[v]:[ovf]>");
  blas->solve("t2_eqns[oo][vf]{u} += #2134# - t1[o][v]{u} 2@1 <[v]:[ovf]>");
  blas->solve("t2_eqns[oo][vf]{u} += #3412# - t1[o][v]{u} 1@1 <[o]:[foo]>");

  blas->solve("t2_eqns[oO][vF]{u} += #1234#   t1[o][v]{u} 2@1 <[v]|[ovf]>");
  blas->solve("t2_eqns[oO][vF]{u} += #2143#   t1[O][V]{u} 2@1 <[v]|[ofv]>");
  blas->solve("t2_eqns[oO][vF]{u} += #3412# - t1[o][v]{u} 1@1 <[o]|[foo]>");

  blas->solve("t2_eqns[OO][VF]{u} += #1234#   t1[O][V]{u} 2@1 <[v]:[ovf]>");
  blas->solve("t2_eqns[OO][VF]{u} += #2134# - t1[O][V]{u} 2@1 <[v]:[ovf]>");
  blas->solve("t2_eqns[OO][VF]{u} += #3412# - t1[O][V]{u} 1@1 <[o]:[foo]>");

  // Solve for second-order amplitudes
  blas->solve("t2_2[oo][vf]{u}   = t2_eqns[oo][vf]{u} / d2[oo][vf]{u}");
  blas->solve("t2_2[oO][vF]{u}   = t2_eqns[oO][vF]{u} / d2[oO][vF]{u}");
  blas->solve("t2_2[OO][VF]{u}   = t2_eqns[OO][VF]{u} / d2[OO][VF]{u}");

  // Compute third-order energy
  blas->solve("Eaaaa{u} = 1/4 t2_2[oo][vf]{u} . <[oo]:[vf]>");
  blas->solve("Eabab{u} =     t2_2[oO][vF]{u} . <[oo]|[vf]>");
  blas->solve("Ebbbb{u} = 1/4 t2_2[OO][VF]{u} . <[oo]:[vf]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_vf_3_s = blas->get_scalar("ECCSD",0);


  // OOFF Terms
  // PPPP Ladder
  blas->solve("t2_eqns[oo][ff]{u} += 1/2 t2[oo][vv]{u} 2@2 <[ff]:[vv]>");

  blas->solve("t2_eqns[oO][fF]{u} +=     t2[oO][vV]{u} 2@2 <[ff]|[vv]>");

  blas->solve("t2_eqns[OO][FF]{u} += 1/2 t2[OO][VV]{u} 2@2 <[ff]:[vv]>");

  // Solve for second-order amplitudes
  blas->solve("t2_2[oo][ff]{u}   = t2_eqns[oo][ff]{u} / d2[oo][ff]{u}");
  blas->solve("t2_2[oO][fF]{u}   = t2_eqns[oO][fF]{u} / d2[oO][fF]{u}");
  blas->solve("t2_2[OO][FF]{u}   = t2_eqns[OO][FF]{u} / d2[OO][FF]{u}");

  // Compute third-order energy
  blas->solve("Eaaaa{u} = 1/4 t2_2[oo][ff]{u} . <[oo]:[ff]>");
  blas->solve("Eabab{u} =     t2_2[oO][fF]{u} . <[oo]|[ff]>");
  blas->solve("Ebbbb{u} = 1/4 t2_2[OO][FF]{u} . <[oo]:[ff]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_ff_3 = blas->get_scalar("ECCSD",0);

  // T1 Contributions
  blas->solve("t2_eqns[oo][ff]{u} += #1234#   t1[o][v]{u} 2@1 <[v]:[off]>");
  blas->solve("t2_eqns[oo][ff]{u} += #2134# - t1[o][v]{u} 2@1 <[v]:[off]>");

  blas->solve("t2_eqns[oO][fF]{u} += #1234#   t1[o][v]{u} 2@1 <[v]|[off]>");
  blas->solve("t2_eqns[oO][fF]{u} += #2143#   t1[O][V]{u} 2@1 <[v]|[off]>");

  blas->solve("t2_eqns[OO][FF]{u} += #1234#   t1[O][V]{u} 2@1 <[v]:[off]>");
  blas->solve("t2_eqns[OO][FF]{u} += #2134# - t1[O][V]{u} 2@1 <[v]:[off]>");

  // Solve for second-order amplitudes
  blas->solve("t2_2[oo][ff]{u}   = t2_eqns[oo][ff]{u} / d2[oo][ff]{u}");
  blas->solve("t2_2[oO][fF]{u}   = t2_eqns[oO][fF]{u} / d2[oO][fF]{u}");
  blas->solve("t2_2[OO][FF]{u}   = t2_eqns[OO][FF]{u} / d2[OO][FF]{u}");

  // Compute third-order energy
  blas->solve("Eaaaa{u} = 1/4 t2_2[oo][ff]{u} . <[oo]:[ff]>");
  blas->solve("Eabab{u} =     t2_2[oO][fF]{u} . <[oo]|[ff]>");
  blas->solve("Ebbbb{u} = 1/4 t2_2[OO][FF]{u} . <[oo]:[ff]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_ff_3_s = blas->get_scalar("ECCSD",0);

  fprintf(outfile,"\n\n        CBS third-order  correction (vf) = %20.12f (no singles)",E_vf_3);
  fprintf(outfile,"\n        CBS third-order  correction (fv) = %20.12f (no singles)",E_vf_3);
  fprintf(outfile,"\n        CBS third-order  correction (ff) = %20.12f (no singles)",E_ff_3);
  fprintf(outfile,"\n        CBS third-order  correction      = %20.12f (no singles)",E_vf_3 * 2.0 + E_ff_3);
  fprintf(outfile,"\n\n        CBS third-order  correction (vf) = %20.12f",E_vf_3_s);
  fprintf(outfile,"\n        CBS third-order  correction (fv) = %20.12f",E_vf_3_s);
  fprintf(outfile,"\n        CBS third-order  correction (ff) = %20.12f",E_ff_3_s);
  fprintf(outfile,"\n        CBS third-order  correction      = %20.12f",E_vf_3_s * 2.0 + E_ff_3_s);

  fprintf(outfile,"\n\n      * CBS corrected energy (2)         = %20.12f",current_energy + E_ff + E_vf + E_fv);
  fprintf(outfile,"\n      * CBS corrected energy (2+3)       = %20.12f (no singles)",current_energy + E_ff + E_vf + E_fv + E_vf_3 * 2.0 + E_ff_3);
  fprintf(outfile,"\n      * CBS corrected energy (2+3)       = %20.12f",current_energy + E_ff + E_vf + E_fv + E_vf_3_s * 2.0 + E_ff_3_s);

  fflush(outfile);
}

}}

/* End Namespaces */
