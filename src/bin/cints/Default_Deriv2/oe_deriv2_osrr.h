#ifndef _psi_src_bin_cints_Default_Deriv2_oe_deriv2_osrr_h
#define _psi_src_bin_cints_Default_Deriv2_oe_deriv2_osrr_h

/*! \file oe_deriv2_osrr.h
    \ingroup CINTS
*/
namespace psi { 
  namespace CINTS {

void AI_Deriv2_OSrecurs(double ***AI0, double ***AIX, double ***AIY, double ***AIZ,
			double ***AIXX, double ***AIXY, double ***AIXZ,
			double ***AIYY, double ***AIYZ, double ***AIZZ,
			struct coordinates PA, struct coordinates PB,
			struct coordinates PC, double gamma, int iang, int jang);

}}
#endif
