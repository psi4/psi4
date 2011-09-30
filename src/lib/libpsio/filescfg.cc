#include <boost/shared_ptr.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

#include <string>
#include <map>
#include <sstream>
#include <algorithm>
#include <cctype>

#include<iostream>
#include<ostream>

namespace psi {

std::string fullkwd(const char* kwdgrp, const char* kwd, int unit) {
  std::string unitname;
  if (unit < 0)
    unitname = "DEFAULT";
  else {
    std::ostringstream oss;
    oss << "FILE"<< unit;
    unitname = oss.str();
  }
  const std::string sep(":");

  std::string fkwd = sep + kwdgrp + sep + "FILES"+ sep + unitname + sep + kwd;
  // convert to upper case
  std::transform(fkwd.begin(), fkwd.end(), fkwd.begin(), static_cast<int(*)(int)>(toupper));
  return fkwd;
}

void PSIO::filecfg_kwd(const char* kwdgrp, const char* kwd, int unit,
                       const char* kwdval) {
  std::string fkwd = fullkwd(kwdgrp, kwd, unit);
  files_keywords_[fkwd] = kwdval;
}

const std::string&PSIO::filecfg_kwd(const char* kwdgrp, const char* kwd,
                                    int unit) {
  static std::string nullstr;

  const std::string fkwd = fullkwd(kwdgrp, kwd, unit);
  KWDMap::const_iterator kwd_loc = files_keywords_.find(fkwd);
  if (kwd_loc != files_keywords_.end())
    return kwd_loc->second;
  else
    return nullstr;
}

int psio_set_filescfg_kwd(const char* kwdgrp, const char* kwd, int unit,
                          const char* kwdval) {
  _default_psio_lib_->filecfg_kwd(kwdgrp, kwd, unit, kwdval);
  return 1;
}

const char* psio_get_filescfg_kwd(const char* kwdgrp, const char* kwd,
                                  int unit) {
  return _default_psio_lib_->filecfg_kwd(kwdgrp,kwd,unit).c_str();
}

}

