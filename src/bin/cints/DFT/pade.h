#ifndef _psi_src_bin_cints_DFT_pade_h
#define _psi_src_bin_cints_DFT_pade_h

/*! \file pade.h
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
namespace psi { namespace CINTS {
double Pade_int(double p, double x0, double b, double c, double A, double Q);
double d_Pade_int(double p, double x0, double b, double c, double A, double Q);
}}
#endif
