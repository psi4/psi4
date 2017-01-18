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

#ifndef MULTIPOLESYMMETRY_H
#define MULTIPOLESYMMETRY_H

#include <vector>
#include <string>
#include <map>

namespace psi{

class Molecule;
class IntegralFactory;
class MatrixFactory;

class OperatorSymmetry
{
    // The order of the multipole (dipole=1, quadrupole=2, etc...)
    int order_;

    // Variables we need from the user
    std::shared_ptr<Molecule> molecule_;
    std::shared_ptr<IntegralFactory> integral_;
    std::shared_ptr<MatrixFactory> matrix_;

    /**
     * The symmetry of each component of the multipole.
     * Length = INT_NCART(order_)
     */
    std::vector<int> component_symmetry_;

    /**
     * Forms a string that describes what the name of the
     * componet is.
     * (x=1, y=0, z=0) => "x"
     * (x=1, y=0, z=1) => "xz"
     * (x=0, y=2, z=0) => "y2"
     */
    std::string form_suffix(int x, int y, int z);

    void common_init();

public:
    enum Operator {
        Dipole = 1,
        Quadrupole = 2,

        L = -1,
        AngularMomentum = -1,
        P = Dipole,
        Nabla = Dipole
    };

    /** Constructor
     * Constructs an object that determines the symmetry of the different
     * components of the order-th multipole. Orders with negative values
     * have a special meaning (see OperatorSymmetry enum)
     *
     * @param order Order of the multipole (1 = dipole, 2 = quadrupole, etc.)
     * @param mol Molecule the the multipole will be computed for. Needed to obtain
     *            point group object.
     * @param ints Integral factory. Needed for creation of ShellRotation objects.
     * @param mats Matrix factory. Used by create_matrices to create matrices of the
     *             proper size and symmetry.
     */
    OperatorSymmetry(int order,
                     std::shared_ptr<Molecule> mol,
                     std::shared_ptr<IntegralFactory> ints,
                     std::shared_ptr<MatrixFactory> mats);
    //OperatorSymmetry(int order,
    //                 std::shared_ptr<Molecule> mol,
    //                 std::shared_ptr<IntegralFactory> ints);
    virtual ~OperatorSymmetry();

    std::string name_of_component(int i);
    int component_symmetry(int i) const { return component_symmetry_[i]; }

    std::vector<SharedMatrix > create_matrices(const std::string& basename);
};

class MultipoleSymmetry
{
    // The order of the multipole (dipole=1, quadrupole=2, etc...)
    int order_;

    // Variables we need from the user
    std::shared_ptr<Molecule> molecule_;
    std::shared_ptr<IntegralFactory> integral_;
    std::shared_ptr<MatrixFactory> matrix_;

    /**
     * The symmetry of each component of the multipole.
     * Length = INT_NCART(order_)
     */
    std::vector<int> component_symmetry_;

    /**
     * A 3D map to hold the addresses of each {lx, ly, lz} combination
     */
    std::map< int, std::map< int, std::map< int, int > > > addresses_;

    void common_init();

public:

    /** Constructor
     * Constructs an object that determines the symmetry of the different
     * components of all multipoles up to (and including) L=order.  For all componenets of a given order
     * use the OperatorSymmetry class instead
     *
     * @param order Order of the highest multipole (1 = dipole, 2 = quadrupole, etc.)
     * @param mol Molecule the the multipole will be computed for. Needed to obtain
     *            point group object.
     * @param ints Integral factory. Needed for creation of ShellRotation objects.
     * @param mats Matrix factory. Used by create_matrices to create matrices of the
     *             proper size and symmetry.
     */
    MultipoleSymmetry(int order,
                     std::shared_ptr<Molecule> mol,
                     std::shared_ptr<IntegralFactory> ints,
                     std::shared_ptr<MatrixFactory> mats);
    //MultipoleSymmetry(int order,
    //                 std::shared_ptr<Molecule> mol,
    //                 std::shared_ptr<IntegralFactory> ints);
    virtual ~MultipoleSymmetry();

    /**
    * Returns the address in the array of the {lx, ly, lz} moment.
    */
    int address_of_component(int lx, int ly, int lz);
    int component_symmetry(int i) const { return component_symmetry_[i]; }

    std::vector<SharedMatrix > create_matrices(const std::string& basename, bool ignore_symmetry=false);
};
}

#endif // MULTIPOLESYMMETRY_H
