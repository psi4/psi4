#ifndef psi_include_psi4_h
#define psi_include_psi4_h

// declarations for modules
#ifdef __cplusplus

extern "C" {
extern FILE *infile, *outfile;
extern char *psi_file_prefix;
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

extern Module module;

}

#else
extern FILE *infile, *outfile;
extern char *psi_file_prefix;
#endif

#endif

