/*! \file
    \ingroup DETCASMAN
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_bin_detcasman_setup_io_h_
#define _psi_src_bin_detcasman_setup_io_h_

namespace psi { namespace detcasman {

extern void init_io(int argc, char *argv[]);
extern void close_io(void);
extern void check(int a, const char *errmsg);

}} // end namespace psi::detcasman

#endif // header guard
