#ifndef psi_include_psi4_dec_h
#define psi_include_psi4_dec_h

#include <boost/shared_ptr.hpp>
#include <boost/current_function.hpp>

#include <string>
#include <list>
#include <map>
#include <liboptions/liboptions.h>
#include <exception.h>
#include <libmints/typedefs.h>

namespace psi {

enum PsiReturnType {Success, Failure, Balk, EndLoop};

extern FILE *outfile;
//  extern PSIO *psio;
extern char *psi_file_prefix;
extern bool verbose;

// Very useful regex for matching floating point numbers
#define NUMBER "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

class Molecule;
class Wavefunction;
class ExternalPotential;

class Process
{
public:
    class Environment
    {
        std::map<std::string, std::string> environment_;
        unsigned long int memory_;
        int nthread_;

        boost::shared_ptr<Molecule> molecule_;
        SharedMatrix gradient_;
        boost::shared_ptr<Wavefunction> reference_wavefunction_;
    public:
        void init(char **envp);

        const std::string& operator()(const std::string& key) const;
        std::string operator()(const std::string& key);
        const std::string& set(const std::string& key, const std::string& value);

        /// Set active molecule
        void set_molecule(const boost::shared_ptr<Molecule>& molecule);
        /// Return active molecule
        boost::shared_ptr<Molecule> molecule() const;

        /// Set reference wavefunction
        void set_reference_wavefunction(const boost::shared_ptr<Wavefunction>& reference_wavefunction);
        /// Get reference wavefunction
        boost::shared_ptr<Wavefunction> reference_wavefunction() const;

        /// Set gradient manually
        void set_gradient(const SharedMatrix gradient) { gradient_ = gradient; }
        /// Get gradient manually
        SharedMatrix gradient() const { return gradient_; }

        /// Map containing current energies
        std::map<std::string, double> globals;

        /** User specified basis files.
           *  These are specific files:  ~/basis/dz.gbs, ~/basis/tz.gbs
           */
        std::list<std::string> user_basis_files;

        /// Number of threads per process
        int get_n_threads() const;
        void set_n_threads(int nthread);

        /// Memory in bytes
        unsigned long int get_memory() const;
        void set_memory(unsigned long int m);

        /// "Global" liboptions object.
        Options options;
    };

    class Arguments
    {
        std::vector<std::string> arguments_;

    public:
        void init(int argc, char **argv);

        int argc() const;

        const std::string& operator()(int argc) const;
        std::string operator()(int argc);
    };

    static Environment environment;
    static Arguments arguments;

    static Environment get_environment();
};
}


#endif
