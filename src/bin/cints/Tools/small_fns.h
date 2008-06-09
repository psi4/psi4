#ifndef _psi_src_bin_cints_Tools_small_fns_h
#define _psi_src_bin_cints_Tools_small_fns_h

/*! \file small_fns.h
    \ingroup CINTS
*/

namespace psi { namespace CINTS {

void setup();
void start_io(int argc, char *argv[]);
void stop_io();
void punt(char* message);
double distance_calc(struct coordinates g1, struct coordinates g2);
double ***init_box(int, int, int);
void free_box(double ***, int, int);
};}
#endif
