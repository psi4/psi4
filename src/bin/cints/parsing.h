#ifndef _psi_src_bin_cints_parsing_h
#define _psi_src_bin_cints_parsing_h

/*! \file
    \ingroup CINTS
    \brief parsing header file.
*/
#include <liboptions/liboptions.h>

namespace psi { namespace cints {
void parsing(Options &);
void parsing_cmdline(Options &);
}}
#endif
