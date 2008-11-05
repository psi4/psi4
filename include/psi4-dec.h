#ifndef psi_include_psi4_dec_h
#define psi_include_psi4_dec_h

#include <string>

namespace psi {

class Module {
  std::string prgid;

  public:
  Module(std::string s = "PSI4") { prgid = s; }
  void set_prgid (std::string s) { prgid = s; }
  // const std::string & gprgid(void) const { return prgid; }
  std::string gprgid(void) const {
    return prgid;
  }
};

extern Module module;
extern FILE *outfile;
extern char *psi_file_prefix;
}

#endif

