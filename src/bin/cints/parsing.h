#ifndef _psi_src_bin_cints_parsing_h
#define _psi_src_bin_cints_parsing_h

/*! \file
    \ingroup CINTS
    \brief parsing header file.
*/
#include <liboptions/liboptions.hpp>

namespace psi { namespace CINTS {
void parsing(Options &);
void parsing_cmdline(int argc, char *argv[]);
}}
#endif
