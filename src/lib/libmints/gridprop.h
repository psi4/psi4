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

#ifndef GRIDPROP_H
#define GRIDPROP_H

#include <libmints/typedefs.h>
#include <libmints/wavefunction.h>
#include <map>
#include <set>

namespace psi {

class Options;
class BasisExtents;
class RKSFunctions;
class BlockOPoints;

class CubicScalarGrid {

protected:

    // => Input Specification <= //
    
    /// Options object for overages and voxel spacing
    Options& options_;
    /// Molecule this grid is built around
    boost::shared_ptr<Molecule> mol_;
    /// Basis set this grid is built around
    boost::shared_ptr<BasisSet> primary_;
    /// File path for grid storage
    std::string filepath_;

    // => Physical Grid <= //

    /// number of points of grid
    size_t npoints_;
    /// Sparsity blocking in all cardinal directions
    size_t nxyz_;

    // Voxel quanta in x, y, z [(N_x+1) x (N_y + 1) x (N_z + 1)]
    int* N_;
    // Voxel spacing in x, y, z
    double* D_;
    // Voxel origin in x, y, z
    double* O_;

    /// x coordinates of grid
    double* x_;
    /// y coordinates of grid
    double* y_;
    /// z coordinates of grid
    double* z_;
    /// Scalar field values on grid
    double* v_;

    // => Grid Computers <= //
    
    /// Vector of blocks 
    std::vector<boost::shared_ptr<BlockOPoints> > blocks_;
    /// Points to basis extents, built internally
    boost::shared_ptr<BasisExtents> extents_;
    /// RKS points object
    boost::shared_ptr<RKSFunctions> points_;

    // => Helper Routines <= //

    /// Setup grid from info in N_, D_, O_
    void populate_grid();

public:
    CubicScalarGrid(boost::shared_ptr<BasisSet> primary);

    virtual ~CubicScalarGrid(); 

    // => Setup Routines <= //

    /// Build grid with options overages
    void build_grid();
    /// Build grid with specified geometry (e.g., from another grid)
    void build_grid(const std::string filepath, int* N, double* D, double* O);
    /// Set the directory for cube file storage
    void set_filepath(const std::string filepath) { filepath_ = filepath; }  
    /// Header info
    void print_header();

    // => Accessor Routines <= //

    int* N() const { return N_; }
    double* D() const { return D_; }
    double* O() const { return O_; }

    double* x() const { return x_; } 
    double* y() const { return y_; } 
    double* z() const { return z_; } 
    double* v() const { return v_; } 

    size_t npoints() const { return npoints_; }
    size_t nxyz() const { return nxyz_; }
    
    // => Low-Level Scalar Field Computation <= //

    /// Zero the scalar field
    void zero();

    /// Add a density-type property to the scalar field
    void compute_density(boost::shared_ptr<Matrix> D, double scale = 1.0);
    /// Add an orbital-type property to the scalar field
    void compute_orbital(boost::shared_ptr<Vector> C, double scale = 1.0);
    /// Add a basis-function-type property to the scalar field
    void compute_basis(int function, double scale = 1.0);
    
    /// Write a cube file with the current scalar field to filepath/name.cube
    void write_cube_file(const std::string& name, double* v = NULL);

    // => High-Level Scalar Field Computation <= //

    /// Compute a Gaussian Cube file corresponding to the density matrix D
    void compute_density_cube(boost::shared_ptr<Matrix> D, const std::string name);
    /// Compute a series of Gaussian Cube files corresponding to the orbitals in indices (or all by default)
    void compute_orbital_cubes(boost::shared_ptr<Matrix> C, const std::string name, const std::vector<int>& indices = std::vector<int>());
    /// Compute a series of Gaussian Cube files corresponding to the basis functions in indices (or all by default)
    void compute_basis_cubes(const std::string name, const std::vector<int>& indices = std::vector<int>());

};

} // End namespace

#endif

