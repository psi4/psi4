#ifndef _psi_src_lib_libmints_basisset_parser_h_
#define _psi_src_lib_libmints_basisset_parser_h_

#include <vector>
#include <string>
#include <psi4-dec.h>
#include "gshell.h"

// Forward declare boost boost::shared_ptr
namespace boost {
template<class T> class shared_ptr;

}

namespace psi {

class BasisSet;

class BasisSetFileNotFound : public PsiException {
public:
    /**
     * Constructor
     * @param message The message that will be printed by exception
     * @param file The file that threw the exception (use __FILE__ macro)
     * @param line The line number that threw the exception (use __LINE__ macro)
     */
    BasisSetFileNotFound(
        std::string message,
        const char* file,
        int line
        ) throw();

    virtual ~BasisSetFileNotFound() throw();
};

class BasisSetNotFound : public PsiException {
public:
    /**
     * Constructor
     * @param message The message that will be printed by exception
     * @param file The file that threw the exception (use __FILE__ macro)
     * @param line The line number that threw the exception (use __LINE__ macro)
     */
    BasisSetNotFound(
        std::string message,
        const char* file,
        int line
        ) throw();

    virtual ~BasisSetNotFound() throw();
};

/*! @ingroup MINTS
    @class BasisSetParser
    @brief Abstract class for parsing basis sets from a text file.

    Provides an interface for parsing basis sets from a text file.
*/
class BasisSetParser
{
public:
    BasisSetParser();
    virtual ~BasisSetParser();

    /** Load and return the file to be used by parse.
     *  @param basisname If specified only return only lines that pertain to that basis name. (for multi-basisset files)
     *                   Otherwise return the entire file is basisname="".
     */
    std::vector<std::string> load_file(const std::string& filename, const std::string& basisname="");

    //! Take a multiline string and convert it to a vector of strings.
    std::vector<std::string> string_to_vector(const std::string& data);

    /**
     * Given a string, parse for the basis set needed for atom.
     * @param basisset object to add to
     * @param atom atom index to look for in basisset->molecule()
     * @param dataset data set to look through
     */
    virtual std::vector<GaussianShell> parse(const std::string& symbol, const std::string& dataset) {
        return parse(symbol, string_to_vector(dataset));
    }

    /**
     * Given a string, parse for the basis set needed for atom.
     * @param basisset object to add to
     * @param atom atom index to look for in basisset->molecule()
     * @param dataset data set to look through
     */
    virtual std::vector<GaussianShell> parse(const std::string& symbol, const std::vector<std::string>& dataset) = 0;
};

/*! \class Gaussian94BasisSetParser
    \brief Class for reading in basis sets formatted for Gaussian.
*/
class Gaussian94BasisSetParser : public BasisSetParser
{
public:
    Gaussian94BasisSetParser() {}

    virtual std::vector<GaussianShell> parse(const std::string& symbol, const std::vector<std::string>& dataset);
};

} /* end psi namespace */

#endif
