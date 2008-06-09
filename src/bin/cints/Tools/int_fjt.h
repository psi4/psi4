#ifndef _psi_src_bin_cints_Tools_int_fjt_h
#define _psi_src_bin_cints_Tools_int_fjt_h

/*! \file int_fjt.h
    \ingroup CINTS
*/

namespace psi { namespace CINTS {
void init_fjt(int max);
void free_fjt();
void init_fjt_table(double_array_t *table);
void free_fjt_table(double_array_t *table);
void int_fjt(double_array_t *table, int J, double wval);
};}
#endif
