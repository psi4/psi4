/*!
 \file
 \ingroup PSIO
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <psi4-dec.h>

namespace psi {

void PSIO::get_volpath(unsigned int unit, unsigned int volume, char **path) {
  std::string kval;
  char volumeX[20];
  sprintf(volumeX, "VOLUME%u", volume+1);
  
  std::string module_name = module.gprgid();
  kval = filecfg_kwd(module_name.c_str(), volumeX, unit);
  if (!kval.empty()) {
    *path = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd(module_name.c_str(), volumeX, -1);
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

