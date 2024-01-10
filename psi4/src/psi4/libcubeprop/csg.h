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

#ifndef _psi_src_lib_libcubeprop_csg_h_
#define _psi_src_lib_libcubeprop_csg_h_

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace psi {

class BasisExtents;
class BasisSet;
class BlockOPoints;
class Matrix;
class Molecule;
class Options;
class RKSFunctions;

class CubicScalarGrid {
   protected:
    // => Input Specification <= //

    /// Options object for overages and voxel spacing
    Options& options_;
    /// Molecule this grid is built around
    std::shared_ptr<Molecule> mol_;
    /// Basis set this grid is built around
    std::shared_ptr<BasisSet> primary_;
    /// The auxiliary basisset for ESP contractions
    std::shared_ptr<BasisSet> auxiliary_;
    /// File path for grid storage
    std::string filepath_;

    // => Physical Grid <= //

    /// Voxel quanta in x, y, z [(N_x+1) x (N_y + 1) x (N_z + 1) points]
    std::array<int, 3> N_;
    /// Voxel spacing in x, y, z
    std::array<double, 3> D_;
    /// Voxel origin in x, y, z
    std::array<double, 3> O_;

    /// number of points of grid
    size_t npoints_;
    /// Sparsity blocking in all cardinal directions
    size_t nxyz_;

    /// x coordinates of grid
    double* x_;
    /// y coordinates of grid
    double* y_;
    /// z coordinates of grid
    double* z_;
    /// w quadrature weights of grid (rectangular)
    double* w_;

    // => Grid Computers <= //

    /// Vector of blocks
    std::vector<std::shared_ptr<BlockOPoints> > blocks_;
    /// Points to basis extents, built internally
    std::shared_ptr<BasisExtents> extents_;
    /// RKS points object
    std::shared_ptr<RKSFunctions> points_;

    // => Helper Routines <= //

    /// Setup grid from info in N_, D_, O_
    void populate_grid();

   public:
    // => Constructors <= //

    CubicScalarGrid(std::shared_ptr<BasisSet> primary, Options& options);
    virtual ~CubicScalarGrid();

    // => High-Level Setup Routines <= //

    /// Build grid with options overages
    void build_grid();
    /// Build grid with specified geometry (e.g., from another grid)
    void build_grid(const std::string filepath, int* N, double* D, double* O);
    /// Build grid from the dimensions in another grid
    void build_grid(std::shared_ptr<CubicScalarGrid> other);
    /// Header info
    void print_header();

    /// Set the auxiliary for ESP if desired
    void set_auxiliary_basis(std::shared_ptr<BasisSet> aux) { auxiliary_ = aux; }

    // => High-Level Set Routines <= //

    /// Set the directory for cube file storage (defaults to "./")
    void set_filepath(const std::string filepath) { filepath_ = filepath; }

    // => High-Level Accessor Routines <= //

    /// Number of voxels in [x,y,z]. Number of points along each dimensions in N_i + 1.
    const std::array<int, 3>& N() const { return N_; }
    /// Voxel width in [x,y,z], in bohr
    const std::array<double, 3>& D() const { return D_; }
    /// Lower-left corner of th grid, in bohr
    const std::array<double, 3>& O() const { return O_; }
    /// Filepath where grid output will be stored
    std::string filepath() const { return filepath_; }

    // => Low-Level Accessor Routines (Use only if you know what you are doing) <= //

    /// Total number of points in grid
    size_t npoints() const { return npoints_; }
    /// Number of points
    size_t nxyz() const { return nxyz_; }

    /// x points in fast ordering
    double* x() const { return x_; }
    /// y points in fast ordering
    double* y() const { return y_; }
    /// z points in fast ordering
    double* z() const { return z_; }
    /// w weights (rectangular) in fast ordering
    double* w() const { return w_; }

    // => Low-Level Write Routines (Use only if you know what you are doing) <= //

    /// Write a general file of the scalar field v (in fast ordering) to filepath/name.ext
    void write_gen_file(double* v, const std::string& name, const std::string& type, const std::string& comment = "");
    /// Write a Gaussian cube file of the scalar field v (in fast ordering) to filepath/name.cube
    void write_cube_file(double* v, const std::string& name, const std::string& comment = "");

    // => Low-Level Scalar Field Computation (Use only if you know what you are doing) <= //

    /// Add a density-type property to the scalar field
    void add_density(double* v, std::shared_ptr<Matrix> D);
    /// Add an ESP-type property to the scalar field (total density matrix, must set DF_BASIS_SCF option)
    void add_esp(double* v, std::shared_ptr<Matrix> D, const std::vector<double>& nuc_weights = std::vector<double>());
    /// Add a basis function property for desired indices to the scalar fields in v (rows are basis functions)
    void add_basis_functions(double** v, const std::vector<int>& indices);
    /// Add orbital property for desired indices to the scalar fields in v (rows are orbitals)
    void add_orbitals(double** v, std::shared_ptr<Matrix> C);
    /// Add a LOL-type property to the scalar field
    void add_LOL(double* v, std::shared_ptr<Matrix> D);
    /// Add an ELF-type property to the scalar field
    void add_ELF(double* v, std::shared_ptr<Matrix> D);

    // => High-Level Scalar Field Computation <= //

    /// Compute a density-type property and drop a file corresponding to name and type
    void compute_density(std::shared_ptr<Matrix> D, const std::string& name, const std::string& type = "CUBE");
    /// Compute an ESP-type property and drop a file corresponding to name and type
    void compute_esp(std::shared_ptr<Matrix> D, const std::vector<double>& nuc_weights, const std::string& name,
                     const std::string& type = "CUBE");
    /// Compute a set of basis function-type properties and drop files corresponding to name, index, and type
    void compute_basis_functions(const std::vector<int>& indices, const std::string& name,
                                 const std::string& type = "CUBE");
    /// Compute a set of orbital-type properties and drop files corresponding to name, index, symmetry label, and type
    void compute_orbitals(std::shared_ptr<Matrix> C, const std::vector<int>& indices,
                          const std::vector<std::string>& labels, const std::string& name,
                          const std::string& type = "CUBE");
    /// Compute a set of orbital-type properties and drop files corresponding to name, index, symmetry label, and type
    void compute_difference(std::shared_ptr<Matrix> C, const std::vector<int>& indices,
                          const std::string& label, bool square = false, const std::string& type = "CUBE");

    /// Compute a LOL-type property and drop a file corresponding to name and type
    void compute_LOL(std::shared_ptr<Matrix> D, const std::string& name, const std::string& type = "CUBE");
    /// Compute an ELF-type property and drop a file corresponding to name and type (TODO: this seems very unstable)
    void compute_ELF(std::shared_ptr<Matrix> D, const std::string& name, const std::string& type = "CUBE");

    /// Compute the isocountour range that capture a given fraction of a property. Exponent is used
    /// to properly compute the density. E.g. for orbitals exponent = 2, for densities exponent = 1
    std::pair<double, double> compute_isocontour_range(double* v2, double exponent);

    /// A helper function to construct ECP comment in the cubefile header.
    std::string ecp_header();
};

}  // namespace psi

#endif
