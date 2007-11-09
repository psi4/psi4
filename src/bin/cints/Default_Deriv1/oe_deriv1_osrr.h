#ifndef _psi_src_bin_cints_Default_Deriv1_oe_deriv1_osrr_h
#define _psi_src_bin_cints_Default_Deriv1_oe_deriv1_osrr_h

/*! \file oe_deriv1_osrr.h
    \ingroup (CINTS)
*/namespace psi { namespace CINTS {

void AI_Deriv1_OSrecurs(double ***AI0, double ***AIX, double ***AIY, double ***AIZ, struct coordinates PA, struct coordinates PB,
		 struct coordinates PC, double gamma, int iang, int jang);

}}
#endif
