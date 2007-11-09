#ifndef _psi_src_bin_cints_Tools_read_gen_opdm_h
#define _psi_src_bin_cints_Tools_read_gen_opdm_h

/*! \file read_gen_opdm.h
    \ingroup (CINTS)
*/namespace psi { namespace CINTS {

  void read_gen_opdm();

void read_density(FILE *fpo, double ***Dens, double ***WDens, int num_ao);

}}
#endif
