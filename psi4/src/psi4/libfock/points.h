/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef libfock_points_H
#define libfock_points_H

#include "psi4/libmints/typedefs.h"
#include "psi4/pragma.h"

#include <cstdio>
#include <map>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <string>

namespace psi {

class BasisSet;
class Vector;
class BlockOPoints;

class PSI_API BasisFunctions {
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
    /// PHI : values of AOs at points in space
    std::map<std::string, SharedMatrix> basis_values_;
    /// Map of temp names to Matrices containing temps
    std::map<std::string, SharedMatrix> basis_temps_;
    /// Allocate registers
    virtual void allocate();

   public:
    // => Constructors <= //

    BasisFunctions(std::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    virtual ~BasisFunctions();

    // => Computers <= //
    /// Compute value of basis functions (and derivatives) at points in block.
    /// Store in basis_values_. Needed derivatives specified by deriv_.
    void compute_functions(std::shared_ptr<BlockOPoints> block);

    // => Accessors <= //

    SharedMatrix basis_value(const std::string& key) { return basis_values_[key]; }
    std::map<std::string, SharedMatrix>& basis_values() { return basis_values_; }

    int max_functions() const { return max_functions_; }
    int max_points() const { return max_points_; }
    int deriv() const { return deriv_; }

    virtual void print(std::string out_fname = "outfile", int print = 2) const;

    // => Setters <= //

    void set_deriv(int deriv) {
        deriv_ = deriv;
        allocate();
    }
    void set_max_functions(int max_functions) {
        max_functions_ = max_functions;
        allocate();
    }
    void set_max_points(int max_points) {
        max_points_ = max_points;
        allocate();
    }
};

class PointFunctions : public BasisFunctions {
   protected:
    // => Indices <= //

    /// The index of the current referenced block.
    size_t block_index_;

    // Contains a map to the cache the global basis_values
    std::unordered_map<size_t, std::map<std::string, SharedMatrix>>* cache_map_ = nullptr;

    // Contains a pointer to the current map to use for basis_values
    std::map<std::string, SharedMatrix>* current_basis_map_ = nullptr;

    /// Ansatz (0 - LSDA, 1 - GGA, 2 - Meta-GGA)
    int ansatz_;
    /// Map of value names to Vectors containing values
    std::map<std::string, std::shared_ptr<Vector>> point_values_;

    // => Orbital Collocation <= //

    /// Map of value names to Matrices containing values
    std::map<std::string, std::shared_ptr<Matrix>> orbital_values_;

   public:
    // => Constructors <= //

    PointFunctions(std::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    ~PointFunctions() override;

    // => Setters <= //
    void set_cache_map(std::unordered_map<size_t, std::map<std::string, SharedMatrix>>* cache_map) {
        cache_map_ = cache_map;
    }

    // => Computers <= //
    /// Compute needed DFT intermediates, e.g. rho, gamma, at the points in block.
    /// ansatz_ determines which intermediates are needed.
    /// force_compute forces basis function values at points to be re-computed.
    virtual void compute_points(std::shared_ptr<BlockOPoints> block, bool force_compute = true) = 0;

    // => Accessors <= //

    std::shared_ptr<Vector> point_value(const std::string& key);
    std::map<std::string, SharedVector>& point_values() { return point_values_; }

    SharedMatrix basis_value(const std::string& key) { return (*current_basis_map_)[key]; }
    std::map<std::string, SharedMatrix>& basis_values() { return (*current_basis_map_); }

    virtual std::vector<SharedMatrix> scratch() = 0;
    virtual std::vector<SharedMatrix> D_scratch() = 0;

    int ansatz() const { return ansatz_; }

    // => Setters <= //

    void set_ansatz(int ansatz) {
        ansatz_ = ansatz;
        deriv_ = ansatz;
        allocate();
    }
    virtual void set_pointers(SharedMatrix Da_occ_AO) = 0;
    virtual void set_pointers(SharedMatrix Da_occ_AO, SharedMatrix Db_occ_AO) = 0;

    // => Orbital Collocation <= //

    std::shared_ptr<Matrix> orbital_value(const std::string& key);
    std::map<std::string, SharedMatrix>& orbital_values() { return orbital_values_; }

    virtual void compute_orbitals(std::shared_ptr<BlockOPoints> block, bool force_compute = true) = 0;
    virtual void set_Cs(SharedMatrix Cocc) = 0;
    virtual void set_Cs(SharedMatrix Caocc, SharedMatrix Cbocc) = 0;
};

class SAPFunctions : public PointFunctions {
   protected:
    // => Temps <= //

    /// Buffer for half-transform
    SharedMatrix temp_;
    /// Build temporary work arrays
    void build_temps();
    /// Allocate registers
    void allocate() override;

   public:
    SAPFunctions(std::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    ~SAPFunctions() override;
    void compute_points(std::shared_ptr<BlockOPoints> block, bool force_compute = true) override;
    std::vector<SharedMatrix> scratch() override;
    void print(std::string out_fname = "outfile", int print = 2) const override;
    size_t block_index() { return block_index_; }
    // Dummy functions
    std::vector<SharedMatrix> D_scratch() override;
    void compute_orbitals(std::shared_ptr<BlockOPoints> block, bool force_compute = true) override;
    void set_Cs(SharedMatrix Cocc) override;
    void set_Cs(SharedMatrix Caocc, SharedMatrix Cbocc) override;
    void set_pointers(SharedMatrix Da_occ_AO) override;
    void set_pointers(SharedMatrix Da_occ_AO, SharedMatrix Db_occ_AO) override;
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
    void allocate() override;

    // => Orbital Collocation <= //

    /// Orbital coefficients, AO
    SharedMatrix C_AO_;
    /// Orbital coefficeints, local AO
    SharedMatrix C_local_;

   public:
    RKSFunctions(std::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    ~RKSFunctions() override;

    void set_pointers(SharedMatrix Da_occ_AO) override;
    void set_pointers(SharedMatrix Da_occ_AO, SharedMatrix Db_occ_AO) override;

    /// Compute the needed DFT intermediates at the points in the block.
    /// "Which DFT intermediates are needed?" is determined from ansatz_.
    /// RHO_A : Spin-summed density
    /// RHO_AX | RHO_AY | RHO_AZ : Spin-summed density gradient components
    /// GAMMA_AA : Spin-summed density gradient dotted with itself
    void compute_points(std::shared_ptr<BlockOPoints> block, bool force_compute = true) override;

    std::vector<SharedMatrix> scratch() override;
    std::vector<SharedMatrix> D_scratch() override;

    void print(std::string out_fname = "outfile", int print = 2) const override;

    void compute_orbitals(std::shared_ptr<BlockOPoints> block, bool force_compute = true) override;
    void set_Cs(SharedMatrix Cocc) override;
    void set_Cs(SharedMatrix Caocc, SharedMatrix Cbocc) override;
    size_t block_index() { return block_index_; }
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
    void allocate() override;

    // => Orbital Collocation <= //

    /// Orbital coefficients, AO
    SharedMatrix Ca_AO_;
    /// Orbital coefficients, AO
    SharedMatrix Cb_AO_;
    /// Orbital coefficients, local AO
    SharedMatrix Ca_local_;
    /// Orbital coefficients, local AO
    SharedMatrix Cb_local_;

   public:
    UKSFunctions(std::shared_ptr<BasisSet> primary, int max_points, int max_functions);
    ~UKSFunctions() override;

    void set_pointers(SharedMatrix Da_occ_AO) override;
    void set_pointers(SharedMatrix Da_occ_AO, SharedMatrix Db_occ_AO) override;
    void set_cache_map(std::unordered_map<size_t, std::map<std::string, SharedMatrix>>* cache_map) {
        cache_map_ = cache_map;
    }

    /// Compute the needed DFT intermediates at the points in the block.
    /// "Which DFT intermediates are needed?" is determined from ansatz_.
    /// RHO_A | RHO_B : Spindensity
    /// RHO_AX | RHO_AY | RHO_AZ | RHO_BX | RHO_BY | RHO_BZ : Spindensity gradient components
    /// GAMMA_AA | GAMMA_AB | GAMMA_BB : Spindensity gradients dotted together
    void compute_points(std::shared_ptr<BlockOPoints> block, bool force_compute = true) override;

    std::vector<SharedMatrix> scratch() override;
    std::vector<SharedMatrix> D_scratch() override;

    void print(std::string out_fname = "outfile", int print = 2) const override;

    void compute_orbitals(std::shared_ptr<BlockOPoints> block, bool force_compute = true) override;
    void set_Cs(SharedMatrix Cocc) override;
    void set_Cs(SharedMatrix Caocc, SharedMatrix Cbocc) override;
    size_t block_index() { return block_index_; }
};
}  // namespace psi
#endif
