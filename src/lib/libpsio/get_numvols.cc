/*!
 \file get_numvols.cc
 \ingroup (PSIO)
 */

#include <stdio.h>
#include <stdlib.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

using namespace psi;

extern "C" {
  extern char *gprgid();
}

unsigned int PSIO::get_numvols(unsigned int unit) {
  std::string charnum;
  
  charnum = filecfg_kwd(gprgid(), "NVOLUME", unit);
  if (!charnum.empty())
    return ((unsigned int)atoi(charnum.c_str()));
  charnum = filecfg_kwd(gprgid(), "NVOLUME", -1);
  if (!charnum.empty())
    return ((unsigned int)atoi(charnum.c_str()));
  charnum = filecfg_kwd("PSI", "NVOLUME", unit);
  if (!charnum.empty())
    return ((unsigned int)atoi(charnum.c_str()));
  charnum = filecfg_kwd("PSI", "NVOLUME", -1);
  if (!charnum.empty())
    return ((unsigned int)atoi(charnum.c_str()));
  charnum = filecfg_kwd("DEFAULT", "NVOLUME", unit);
  if (!charnum.empty())
    return ((unsigned int)atoi(charnum.c_str()));
  charnum = filecfg_kwd("DEFAULT", "NVOLUME", -1);
  if (!charnum.empty())
    return ((unsigned int)atoi(charnum.c_str()));
  
  // assume that the default has been provided already
  abort();
}

extern "C" {
  /*!
   ** PSIO_GET_NUMVOLS_DEFAULT(): Get the number of volumes that file 
   ** number 'unit' is split across.
   **
   ** \ingroup (PSIO)
   */
  unsigned int psio_get_numvols_default(void) {
    std::string charnum;
    
    charnum = _default_psio_lib_->filecfg_kwd("PSI", "NVOLUME", -1);
    if (!charnum.empty())
      return ((unsigned int)atoi(charnum.c_str()));
    charnum = _default_psio_lib_->filecfg_kwd("DEFAULT", "NVOLUME", -1);
    if (!charnum.empty())
      return ((unsigned int)atoi(charnum.c_str()));
    
    // assume that the default has been provided already
    abort();
  }
}

