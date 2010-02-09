#ifndef _psi_src_lib_libmints_basisset_parser_h_
#define _psi_src_lib_libmints_basisset_parser_h_

#include <cstdio>
// Need libint for maximum angular momentum
#include "basisset.h"

namespace psi {

/*! \ingroup MINTS
    \class BasisSetParser
    \brief Abstract class for parsing basis sets from a text file.

    Provides an interface for parsing basis sets from a text file.
*/
class BasisSetParser
{
    std::string searchpath_;
public:
    BasisSetParser(const std::string& searchpath = "");
    virtual ~BasisSetParser();

    const std::string searchpath() const { return searchpath_; }
    /**
     * Takes a basis set and the molecule associated with it and
     * reads in the basis set information. Just basic shell information
     * is to be constructed. The BasisSet object itself is responsible
     * for symmetry adaptation.
     */
    virtual void parse(boost::shared_ptr<BasisSet> &basisset, const std::vector<std::string> &basisnames) = 0;
};

/*! \class Gaussian94BasisSetParser
    \brief Class for reading in basis sets formatted for Gaussian.
*/
class Gaussian94BasisSetParser : public BasisSetParser
{
public:
    Gaussian94BasisSetParser(const std::string& searchpath) : BasisSetParser(searchpath) {}

    void parse(boost::shared_ptr<BasisSet> &basisSet, const std::vector<std::string> &basisnames);
};

} /* end psi namespace */

#endif
