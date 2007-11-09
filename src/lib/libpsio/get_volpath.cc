/*!
 \file get_volpath.cc
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

void PSIO::get_volpath(unsigned int unit, unsigned int volume, char **path) {
  std::string kval;
  char volumeX[20];
  sprintf(volumeX, "VOLUME%u", volume+1);
  
  kval = filecfg_kwd(gprgid(), volumeX, unit);
  if (!kval.empty()) {
    *path = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd(gprgid(), volumeX, -1);
  if (!kval.empty()) {
    *path = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd("PSI", volumeX, unit);
  if (!kval.empty()) {
    *path = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd("PSI", volumeX, -1);
  if (!kval.empty()) {
    *path = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd("DEFAULT", volumeX, unit);
  if (!kval.empty()) {
    *path = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd("DEFAULT", volumeX, -1);
  if (!kval.empty()) {
    *path = strdup(kval.c_str());
    return;
  }
  
  // assume default has been provided
  abort();
}

extern "C" {
  /*
   ** PSIO_GET_VOLPATH_DEFAULT(): Get the default path for the nth volume
   ** of any file.
   **
   ** \ingroup (PSIO)
   */
  int psio_get_volpath_default(unsigned int volume, char **path) {
    std::string kval;
    char volumeX[20];
    sprintf(volumeX, "VOLUME%u", volume+1);
    
    kval = _default_psio_lib_->filecfg_kwd("PSI", volumeX, -1);
    if (!kval.empty()) {
      *path = strdup(kval.c_str());
      return (1);
    }
    kval = _default_psio_lib_->filecfg_kwd("DEFAULT", volumeX, -1);
    if (!kval.empty()) {
      *path = strdup(kval.c_str());
      return (1);
    }
    
    // assume default has been provided
    abort();
  }
}

