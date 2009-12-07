#ifndef psi_include_psi4_dec_h
#define psi_include_psi4_dec_h

#include <string>
#include <libutil/libutil.h>
#include <liboptions/liboptions.h>
#include <exception.h>
#include <libutil/ref.h>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <libpsio/psio.hpp>

using namespace boost;

namespace psi {

  class Module {
    std::string prgid;
    long int memory;

  public:
    Module(std::string s = "PSI4") { prgid = s; }
    void set_prgid (std::string s) { prgid = s; }
    std::string gprgid(void) const { return prgid; }
    long int get_memory(void) const { return memory; }
    void set_memory(long int newmem) { memory = newmem; }
  };

  struct ModuleInformation {
    std::string rubyFunctionName;
    int (*entryPoint)();
    void (*registerOptions)();
  };

  enum PsiReturnType {Success, Failure, Balk};
  
  extern Module module;
  extern FILE *outfile;
  extern PSIO *psio;
  extern char *psi_file_prefix;
  extern bool g_bVerbose;
  extern std::map<std::string, PsiReturnType(*)(Options &, int argc, char *argv[])> dispatch_table;
  
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
