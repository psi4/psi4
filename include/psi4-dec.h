#ifndef psi_include_psi4_dec_h
#define psi_include_psi4_dec_h

#include <string>
#include <libutil/libutil.h>
#include <liboptions/liboptions.h>
#include <exception.h>
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

  enum PsiReturnType {Success, Failure, Balk, Endloop};

  extern Module module;
  extern FILE *outfile;
  extern PSIO *psio;
  extern char *psi_file_prefix;
  extern bool verbose;

  #ifndef NDEBUG
    #ifdef __GNUC__
      #define WHEREAMI() if (verbose) { fprintf(stderr, "@:%s:%-6d\t-> %s\n", __FILE__, __LINE__, __PRETTY_FUNCTION__ ); }
    #else
      #define WHEREAMI() if (verbose) { fprintf(stderr, "@:%s:%-6d\n", __FILE__, __LINE__); }
    #endif
  #else
    #define WHEREAMI()
  #endif

  class Molecule;

  class Process
  {
  public:
      class Environment
      {
          std::map<std::string, std::string> environment_;

          boost::shared_ptr<Molecule> molecule_;
      public:
          void init(char **envp);

          const std::string& operator()(const std::string& key) const;
          std::string operator()(const std::string& key);
          const std::string& set(const std::string& key, const std::string& value);

          // Set active molecule
          void set_molecule(boost::shared_ptr<Molecule> molecule);
          // Return active molecule
          boost::shared_ptr<Molecule> molecule() const;

          /// Map containing current energies
          std::map<std::string, double> globals;
      };

      class Arguments
      {
          std::vector<std::string> arguments_;

      public:
          void init(int argc, char **argv);

          int argc() const { return arguments_.size(); }

          const std::string& operator()(int argc) const;
          std::string operator()(int argc);
      };

      static Environment environment;
      static Arguments arguments;
  };
}

#endif
