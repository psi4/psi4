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

#ifndef libfock_points_H
#define libfock_points_H

#include <cstdio>
#include <map>
#include <tuple>

#include "psi4/libmints/typedefs.h"

namespace psi {

class PSIO;
class BasisSet;
class Matrix;
class Vector;
class IntVector;
class Vector3;
class BlockOPoints;


class BasisFunctions {

protected:
    /// Basis set for this BasisFunctions
    std::shared_ptr<BasisSet> primary_;
    /// Pure AM or not
    bool puream_;
    /// Maximum number of points in a block
    int max_points_;
    /// Maximum number of functions in a block
    int max_functions_;
    /// Maximum derivative to compute
    int deriv_;
    /// Map of value names to Matrices containing values
    std::map<std::string, SharedMatrix > basis_values_;
    /// Map of temp names to Matrices containing temps
    std::map<std::string, SharedMatrix > basis_temps_;
    /// [L]: pure_index, cart_index, coef
    std::vector<std::vector<std::tuple<int,int,double> > > spherical_transforms_;

    /// Setup spherical_transforms_
    void build_spherical();
    /// Allocate registers
    virtual void allocate();

public:
    // => Constructors <= //

    BasisFunctions(std::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    virtual ~BasisFunctions();

    // => Computers <= //

    void compute_functions(std::shared_ptr<BlockOPoints> block);

    // => Accessors <= //

    SharedMatrix basis_value(const std::string& key);
    std::map<std::string, SharedMatrix>& basis_values() { return basis_values_; }

    int max_functions() const { return max_functions_; }
    int max_points() const { return max_points_; }
    int deriv() const { return deriv_; }

    virtual void print(std::string OutFileRMR = "outfile", int print = 2) const;

    // => Setters <= //

    void set_deriv(int deriv) { deriv_ = deriv; allocate(); }
    void set_max_functions(int max_functions) { max_functions_ = max_functions; allocate(); }
    void set_max_points(int max_points) { max_points_ = max_points; allocate(); }
};

class PointFunctions : public BasisFunctions {

protected:
    /// Ansatz (0 - LSDA, 1 - GGA, 2 - Meta-GGA)
    int ansatz_;
    /// Map of value names to Vectors containing values
    std::map<std::string, std::shared_ptr<Vector> > point_values_;

    // => Orbital Collocation <= //

    /// Map of value names to Matrices containing values
    std::map<std::string, std::shared_ptr<Matrix> > orbital_values_;

public:
    // => Constructors <= //

    PointFunctions(std::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    virtual ~PointFunctions();

    // => Computers <= //

    virtual void compute_points(std::shared_ptr<BlockOPoints> block) = 0;

    // => Accessors <= //

    std::shared_ptr<Vector> point_value(const std::string& key);
    std::map<std::string, SharedVector>& point_values() { return point_values_; }

    virtual std::vector<SharedMatrix> scratch() = 0;
    virtual std::vector<SharedMatrix> D_scratch() = 0;

    int ansatz() const { return ansatz_; }

    // => Setters <= //

    void set_ansatz(int ansatz) { ansatz_ = ansatz; deriv_ = ansatz; allocate(); }
    virtual void set_pointers(SharedMatrix Da_occ_AO) = 0;
    virtual void set_pointers(SharedMatrix Da_occ_AO, SharedMatrix Db_occ_AO) = 0;

    // => Orbital Collocation <= //

    std::shared_ptr<Matrix> orbital_value(const std::string& key);
    std::map<std::string, SharedMatrix>& orbital_values() { return orbital_values_; }

    virtual void compute_orbitals(std::shared_ptr<BlockOPoints> block) = 0;
    virtual void set_Cs(SharedMatrix Cocc) = 0;
    virtual void set_Cs(SharedMatrix Caocc, SharedMatrix Cbocc) = 0;
};

class RKSFunctions : public PointFunctions {

protected:
    // => Pointers <= //

    /// Density matrix, AO
    SharedMatrix D_AO_;

    // => Temps <= //

    /// Buffer for half-transform
    SharedMatrix temp_;
    /// Local D matrix
    SharedMatrix D_local_;

    /// Build temporary work arrays
    void build_temps();
    /// Allocate registers
    void allocate();

    // => Orbital Collocation <= //

    /// Orbital coefficients, AO
    SharedMatrix C_AO_;
    /// Orbital coefficeints, local AO
    SharedMatrix C_local_;

public:
    RKSFunctions(std::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    virtual ~RKSFunctions();

    void set_pointers(SharedMatrix Da_occ_AO);
    void set_pointers(SharedMatrix Da_occ_AO, SharedMatrix Db_occ_AO);

    void compute_points(std::shared_ptr<BlockOPoints> block);

    std::vector<SharedMatrix> scratch();
    std::vector<SharedMatrix> D_scratch();

    void print(std::string OutFileRMR = "outfile", int print = 2) const;

    void compute_orbitals(std::shared_ptr<BlockOPoints> block);
    void set_Cs(SharedMatrix Cocc);
    void set_Cs(SharedMatrix Caocc, SharedMatrix Cbocc);
};

class UKSFunctions : public PointFunctions {

protected:
    // => Pointers <= //

    /// Density matrix, AO
    SharedMatrix Da_AO_;
    /// Density matrix, AO
    SharedMatrix Db_AO_;

    // => Temps <= //

    /// Buffer for half-transform
    SharedMatrix tempa_;
    /// Buffer for half-transform
    SharedMatrix tempb_;
    /// Local D matrix
    SharedMatrix Da_local_;
    /// Local D matrix
    SharedMatrix Db_local_;

    /// Build temporary work arrays
    void build_temps();
    /// Allocate registers
    void allocate();

    // => Orbital Collocation <= //

    /// Orbital coefficients, AO
    SharedMatrix Ca_AO_;
    /// Orbital coefficients, AO
    SharedMatrix Cb_AO_;
    /// Orbital coefficeints, local AO
    SharedMatrix Ca_local_;
    /// Orbital coefficeints, local AO
    SharedMatrix Cb_local_;

public:
    UKSFunctions(std::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    virtual ~UKSFunctions();

    void set_pointers(SharedMatrix Da_occ_AO);
    void set_pointers(SharedMatrix Da_occ_AO, SharedMatrix Db_occ_AO);

    void compute_points(std::shared_ptr<BlockOPoints> block);

    std::vector<SharedMatrix> scratch();
    std::vector<SharedMatrix> D_scratch();

    void print(std::string OutFileRMR = "outfile", int print = 2) const;

    void compute_orbitals(std::shared_ptr<BlockOPoints> block);
    void set_Cs(SharedMatrix Cocc);
    void set_Cs(SharedMatrix Caocc, SharedMatrix Cbocc);
};


}
#endif
