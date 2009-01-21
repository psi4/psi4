#ifndef psi_include_psi4_dec_h
#define psi_include_psi4_dec_h

#include <string>
// TODO: Change header file to liboptions.h when appropriate
#include <libpsio/psio.hpp>
#include <libutil/libutil.h>
#include <liboptions/liboptions.hpp>

namespace psi {

  class Module {
    std::string prgid;

  public:
    Module(std::string s = "PSI4") { prgid = s; }
    void set_prgid (std::string s) { prgid = s; }
    std::string gprgid(void) const { return prgid; }
  };

  typedef struct ModuleInformation {
    std::string rubyFunctionName;
    int (*entryPoint)();
    void (*registerOptions)();
  };
  
  extern Module module;
  extern FILE *outfile;
  extern Options options;
  extern PSIO *psio;
  extern char *psi_file_prefix;
  extern bool g_bVerbose;
  
  #ifndef NDEBUG
    #ifdef __GNUC__
      #define WHEREAMI() if (g_bVerbose) { fprintf(stderr, "@:%s:%-6d\t-> %s\n", __FILE__, __LINE__, __PRETTY_FUNCTION__ ); }
    #else
      #define WHEREAMI() if (g_bVerbose) { fprintf(stderr, "@:%s:%-6d\n", __FILE__, __LINE__); }
    #endif
  #else
    #define WHEREAMI()
  #endif
}

#endif
