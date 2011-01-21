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

void PSIO::get_filename(unsigned int unit, char **name, bool remove_namespace) {
  std::string kval;
  std::string dot(".");

  std::string ns = (default_namespace_ == "" || remove_namespace) ? "" : dot + default_namespace_;
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

