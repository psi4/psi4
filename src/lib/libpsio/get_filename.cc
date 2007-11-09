/*!
 \file get_filename.cc
 \ingroup (PSIO)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

using namespace psi;

extern "C" {
  extern char *gprgid();
}

void PSIO::get_filename(unsigned int unit, char **name) {
  std::string kval;
  kval = filecfg_kwd(gprgid(), "NAME", unit);
  if (!kval.empty()) {
    *name = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd(gprgid(), "NAME", -1);
  if (!kval.empty()) {
    *name = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd("PSI", "NAME", unit);
  if (!kval.empty()) {
    *name = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd("PSI", "NAME", -1);
  if (!kval.empty()) {
    *name = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd("DEFAULT", "NAME", unit);
  if (!kval.empty()) {
    *name = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd("DEFAULT", "NAME", -1);
  if (!kval.empty()) {
    *name = strdup(kval.c_str());
    return;
  }
  
  // assume that the default has been provided already
  abort();
}

extern "C" {
  /*!
   ** PSIO_GET_FILENAME_DEFAULT(): Get the default filename
   */
  int psio_get_filename_default(char **name) {
    std::string kval;
    kval = _default_psio_lib_->filecfg_kwd("PSI", "NAME", -1);
    if (!kval.empty()) {
      *name = strdup(kval.c_str());
      return (1);
    }
    kval = _default_psio_lib_->filecfg_kwd("DEFAULT", "NAME", -1);
    if (!kval.empty()) {
      *name = strdup(kval.c_str());
      return (1);
    }
    
    // assume that the default has been provided already
    abort();
  }
}

