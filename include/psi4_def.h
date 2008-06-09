#ifndef psi_include_psi4_def_h
#define psi_include_psi4_def_h

// definitions for modules

extern "C" {
FILE *infile, *outfile;
char *psi_file_prefix;
}

#include <string>

namespace psi {

class Module {
  std::string prgid;

  public:
  Module(std::string s = "PSI4") { prgid = s; }
  void set_prgid (std::string s) { prgid = s; }
  // const std::string & gprgid(void) const { return prgid; }
  char *gprgid(void) const {
    char *name;
    name = new char [prgid.size()+1];
    for (int i=0; i<prgid.size(); ++i)
      name[i] = prgid[i];
    name[prgid.size()] = '\0';
    return name;
  }
};

Module module;

}

#endif
