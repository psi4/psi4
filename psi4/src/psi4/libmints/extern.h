/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#ifndef _psi_src_lib_libmints_extern_potential_h_
#define _psi_src_lib_libmints_extern_potential_h_

#include <vector>
#include <utility>
#include <string>
#include <tuple>
#include "typedefs.h"

namespace psi {

class Matrix;
class Molecule;
class BasisSet;

/*! \ingroup MINTS
 *  \class ExternalPotential
 *  Stores external potential field, computes external potential matrix
 *  Like standard potential integrals, this is negative definite (electrons are the test charge)
 */
class ExternalPotential {
protected:

    /// Debug flag
    int debug_;
    /// Print flag
    int print_;

    /// Name of potential
    std::string name_;
    /// <Z,x,y,z> array of charges
    std::vector<std::tuple<double,double,double,double> > charges_;
    /// Auxiliary basis sets (with accompanying molecules and coefs) of diffuse charges
    std::vector<std::pair<std::shared_ptr<BasisSet>, SharedVector> > bases_;

public:
    /// Constructur, does nothing
    ExternalPotential();
    /// Destructor, does nothing
    ~ExternalPotential();

    /// Set name
    void setName(const std::string & name) { name_ = name; }

    /// Add a charge Z at (x,y,z)
    void addCharge(double Z,double x, double y, double z);
    /// Add a basis of S auxiliary functions with DF coefficients
    void addBasis(std::shared_ptr<BasisSet> basis, SharedVector coefs);

    /// Reset the field to zero (eliminates all entries)
    void clear();

    /// Compute the external potential matrix in the given basis set
    SharedMatrix computePotentialMatrix(std::shared_ptr<BasisSet> basis);
    /// Compute the gradients due to the external potential
    SharedMatrix computePotentialGradients(std::shared_ptr<BasisSet> basis, std::shared_ptr<Matrix> Dt);
    /// Compute the contribution to the nuclear repulsion energy for the given molecule
    double computeNuclearEnergy(std::shared_ptr<Molecule> mol);

    /// Print a trace of the external potential
    void print(std::string OutFileRMR = "outfile") const;

    /// Python print helper
    void py_print() const { print("outfile"); }

    /// Print flag
    void set_print(int print) { print_ = print; }
    /// Debug flag
    void set_debug(int debug) { debug_ = debug; }

};

}

#endif
