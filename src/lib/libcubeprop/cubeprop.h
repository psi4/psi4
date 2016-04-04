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

#ifndef _psi_src_lib_libcubeprop_cubeprop_h_
#define _psi_src_lib_libcubeprop_cubeprop_h_

#include <map>

#include <libmints/typedefs.h>
#include <libmints/wavefunction.h>

namespace psi {

class CubicScalarGrid;

class CubeProperties {

protected:

    // => Task specification <= //

    /// Global options object
    Options& options_;

    // => Key Member Data <= //

    /// Orbital Basis Set
    boost::shared_ptr<BasisSet> basisset_;
    /// AO-basis (C1) OPDM for alpha electrons
    boost::shared_ptr<Matrix> Da_;
    /// AO-basis (C1) OPDM for beta electrons
    boost::shared_ptr<Matrix> Db_;
    /// AO-basis (C1) SCF orbital coefficients for alpha electrons
    boost::shared_ptr<Matrix> Ca_;
    /// AO-basis (C1) SCF orbital coefficients for beta electrons
    boost::shared_ptr<Matrix> Cb_;
    /// Info for alpha electrons (epsilon,rel. index,irrep)
    std::vector<boost::tuple<double, int, int> > info_a_;
    /// Info for beta electrons (epsilon,rel. index,irrep)
    std::vector<boost::tuple<double, int, int> > info_b_;

    // => Computers <= //

    /// Grid-based property computer
    boost::shared_ptr<CubicScalarGrid> grid_;

    // => Helper Functions <= //

    /// Common setup across all constructors
    void common_init();

public:
    // => Constructors <= //

    /// Construct a CubeProperties object from a Wavefunction (possibly with symmetry in wfn)
    CubeProperties(SharedWavefunction wfn);

    /// Common Destructor
    virtual ~CubeProperties();

    // => High-Level Property Computers <= //

    /// Compute all relevant properties from options object specifications
    void compute_properties();

    // => Low-Level Property Computers (Do not use unless you are an advanced client code) <= //

    /// Obligatory title info
    void print_header();
    /// Compute a density grid task (key.cube)
    void compute_density(boost::shared_ptr<Matrix> D, const std::string& key);
    /// Compute an ESP grid task (Dt.cube and ESP.cube)
    void compute_esp(boost::shared_ptr<Matrix> Dt, const std::vector<double>& nuc_weights = std::vector<double>());
    /// Compute an orbital task (key_N.cube, for 0-based indices of C)
    void compute_orbitals(boost::shared_ptr<Matrix> C, const std::vector<int>& indices, const std::vector<std::string>& labels, const std::string& key);
    /// Compute a basis function task (key_N.cube, for 0-based indices of basisset_)
    void compute_basis_functions(const std::vector<int>& indices, const std::string& key);
    /// Compute a LOL grid task (key.cube)
    void compute_LOL(boost::shared_ptr<Matrix> D, const std::string& key);
    /// Compute an ELF grid task (key.cube)
    void compute_ELF(boost::shared_ptr<Matrix> D, const std::string& key);
};

}

#endif
