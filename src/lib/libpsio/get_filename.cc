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

std::string PSIO::current_namespace_;

void PSIO::get_filename(unsigned int unit, char **name) {
  std::string kval;
  std::string module_name = module.gprgid();
  std::string dot("."); 
  std::string ns = (current_namespace_ == "") ? "" : dot + current_namespace_;
  kval = filecfg_kwd(module_name.c_str(), "NAME", unit);
  //printf("File namespace is %s\n",(current_namespace_).c_str());
  if (!kval.empty()) {
    kval = kval + ns;
    *name = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd(module_name.c_str(), "NAME", -1);
  if (!kval.empty()) {
    kval = kval + ns;
    *name = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd("PSI", "NAME", unit);
  if (!kval.empty()) {
    kval = kval + ns;
    *name = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd("PSI", "NAME", -1);
  if (!kval.empty()) {
    kval = kval + ns;
    *name = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd("DEFAULT", "NAME", unit);
  if (!kval.empty()) {
    kval = kval + ns;
    *name = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd("DEFAULT", "NAME", -1);
  if (!kval.empty()) {
    kval = kval + ns;
    *name = strdup(kval.c_str());
    return;
  }
  
  // assume that the default has been provided already
  abort();
}

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

