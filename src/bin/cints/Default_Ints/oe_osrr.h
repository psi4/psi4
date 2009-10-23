#ifndef _psi_src_bin_cints_Default_Ints_oe_osrr_h
#define _psi_src_bin_cints_Default_Ints_oe_osrr_h

/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/namespace psi { namespace cints {

void OI_OSrecurs(double **OIX, double **OIY, double **OIZ, struct coordinates PA, struct coordinates PB,
		 double gamma, int lmaxi, int lmaxj);
void AI_OSrecurs(double ***, struct coordinates, struct coordinates,
                 struct coordinates, double, int, int);
		 }}
		 #endif
