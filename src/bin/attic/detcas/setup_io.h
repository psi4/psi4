/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_bin_detcas_setup_io_h
#define _psi_src_bin_detcas_setup_io_h

namespace psi { namespace detcas {

extern void init_io(int argc, char *argv[]);
extern void close_io(void);
extern void check(int a, char *errmsg);

}} // end namespace psi::detcas

#endif // header guard
