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

void PSIO::change_file_namespace(unsigned int unit, const std::string & ns1, const std::string & ns2) {
    char *old_name, *new_name, *path, *old_fullpath, *new_fullpath;
    _default_psio_lib_->get_filename(unit, &old_name);
    _default_psio_lib_->get_filename(unit, &new_name);
    _default_psio_lib_->get_volpath(unit, 0, &path);  

    old_fullpath = (char*) malloc( (strlen(path)+strlen(old_name)+80)*sizeof(char));
    new_fullpath = (char*) malloc( (strlen(path)+strlen(new_name)+80)*sizeof(char));
    sprintf(old_fullpath, "%s%s.%s.%u", path, old_name, ns1.c_str(), unit);
    sprintf(new_fullpath, "%s%s.%s.%u", path, new_name, ns2.c_str(), unit);

    //printf("%s\n",old_fullpath);
    //printf("%s\n",new_fullpath);

    ::rename(old_fullpath,new_fullpath);
}
/**
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
**/
}

