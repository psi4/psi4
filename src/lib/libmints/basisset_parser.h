/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

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
protected:
    std::string filename_;
public:
    //! If the parser needs to force spherical or cartesian (e.g., loading old guess)
    bool force_puream_or_cartesian_;
    //! Is the forced value to use puream?  (Otherwise force Cartesian).
    bool forced_is_puream_;

    BasisSetParser();
    BasisSetParser(bool forced_puream);
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
    Gaussian94BasisSetParser(): BasisSetParser() {}
    Gaussian94BasisSetParser(bool forced_puream): BasisSetParser(forced_puream) {}
    virtual std::vector<GaussianShell> parse(const std::string& symbol, const std::vector<std::string>& dataset);
};

} /* end psi namespace */

#endif
