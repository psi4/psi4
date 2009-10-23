#ifndef _psi_src_bin_cints_Tools_taylor_fm_eval_h
#define _psi_src_bin_cints_Tools_taylor_fm_eval_h

/*! \file taylor_fm_eval.h
    \ingroup CINTS
*/
namespace psi { namespace cints {
void init_Taylor_Fm_Eval(unsigned int mmax, double epsilon);
void taylor_compute_fm(double *F, double T, unsigned int l);
void free_Taylor_Fm_Eval();
}}
#endif
