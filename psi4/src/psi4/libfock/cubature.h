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

#ifndef libmints_cubature_H
#define libmints_cubature_H

#include "psi4/psi4-dec.h"
#include "psi4/pragma.h"

#include "psi4/libmints/vector3.h"
#include "psi4/libmints/typedefs.h"
#include "psi4/libpsi4util/exception.h"

#include <map>
#include <vector>

namespace psi {

class BasisSet;
class Matrix;
class Vector;
class BasisExtents;
class BlockOPoints;
class RadialGrid;
class SphericalGrid;
class Options;

// This is an auxiliary structure used internally by the grid-builder class.
// Apparently, for performance reasons, it is not good for the final molecular grid
// to consist of structures like this one.
//
// RMP: You got that right. Calling nuclear weights one point at a time is another
// great way to get a 1000x slowdown. What an incredible smell you've discovered!
//
struct MassPoint {
    double x, y, z, w;
};

class MolecularGrid {
   protected:
    int debug_;

    /// The molecule this grid is built on
    std::shared_ptr<Molecule> molecule_;

    // ==> Fast Grid Specification <== //

    /// Total points for this molecule
    int npoints_;
    /// Maximum number of points in a block
    int max_points_;
    /// Maximum number of functions in a block
    int max_functions_;
    // The total collocation size
    size_t collocation_size_;
    /// Full x points.
    double* x_;
    /// Full y points.
    double* y_;
    /// Full z points.
    double* z_;
    /// Full weights
    double* w_;

    // ==> Clean Grid Specification <== //

    /// Orientation matrix
    std::shared_ptr<Matrix> orientation_;
    /// Radial grids, per atom
    std::vector<std::shared_ptr<RadialGrid>> radial_grids_;
    /// Spherical grids, per atom and radial point
    std::vector<std::vector<std::shared_ptr<SphericalGrid>>> spherical_grids_;
    /// Grid points, per atom. Available for any grid blocking scheme unlike atomic_blocks_.
    std::vector<std::vector<MassPoint>> atomic_grids_;

    /// Vector of blocks
    std::vector<std::shared_ptr<BlockOPoints>> blocks_;
    /// Vector of blocks by atoms. Only available if the correct blocking scheme is chosen.
    std::vector<std::vector<std::shared_ptr<BlockOPoints>>> atomic_blocks_;

    /// Points to basis extents, built internally
    std::shared_ptr<BasisExtents> extents_;
    /// BasisSet from extents_
    std::shared_ptr<BasisSet> primary_;

    /// Sieve and block
    void postProcess(std::shared_ptr<BasisExtents> extents, int max_points, int min_points, double max_radius);
    void remove_distant_points(double Rcut);
    void block(int max_points, int min_points, double max_radius);

   public:
    struct MolecularGridOptions {
        double bs_radius_alpha;
        double pruning_alpha;
        short radscheme;  // Effectively an enumeration
        short prunefunction;
        short nucscheme;
        short namedGrid;  // -1 = None, 0 = SG-0, 1 = SG-1
        bool remove_distant_points;
        int nradpts;
        int nangpts;
        int print;
        int debug;
        int bench;
        double weights_cutoff;
        std::string prunescheme;
        std::string prunetype;
        std::string blockscheme;
    };

   protected:
    /// A copy of the options used, for printing purposes.
    MolecularGridOptions options_;

   public:
    MolecularGrid(std::shared_ptr<Molecule> molecule);
    virtual ~MolecularGrid();

    /// Returns the molecule this grid is built on
    std::shared_ptr<Molecule> molecule() const { return molecule_; }

    /// Build the grid
    void buildGridFromOptions(MolecularGridOptions const& opt);
    /// Build the grid
    void buildGridFromOptions(MolecularGridOptions const& opt,
                              const std::vector<std::vector<double>>& rs,  // Radial nodes,     per atom
                              const std::vector<std::vector<double>>& ws,  // Radial weights,   per atom
                              const std::vector<std::vector<int>>& Ls);    // Spherical orders, per atom

    /// Print information about the grid
    void print(std::string out_fname = "outfile", int print = 2) const;
    void print_details(std::string out_fname = "outfile", int print = 2) const;

    /// Orientation matrix
    std::shared_ptr<Matrix> orientation() const { return orientation_; }
    /// Radial grids, per atom
    const std::vector<std::shared_ptr<RadialGrid>>& radial_grids() const { return radial_grids_; }
    /// Spherical grids, per atom and radial point
    const std::vector<std::vector<std::shared_ptr<SphericalGrid>>>& spherical_grids() const { return spherical_grids_; }

    /// Number of grid points
    int npoints() const { return npoints_; }
    /// Maximum number of grid points in a block
    int max_points() const { return max_points_; }
    /// Maximum number of funtions in a block
    int max_functions() const { return max_functions_; }
    /// Total collocation size of all blocks
    size_t collocation_size() { return collocation_size_; }

    /// The x points. You do not own this
    double* x() const { return x_; }
    /// The y points. You do not own this
    double* y() const { return y_; }
    /// The z points. You do not own this
    double* z() const { return z_; }
    /// The weights, normalized to 1 on R3. You do not own this
    double* w() const { return w_; }

    /// Pointer to basis extents
    std::shared_ptr<BasisExtents> extents() const { return extents_; }
    /// Set of spatially sieved blocks of points, generated by sieve() internally
    const std::vector<std::shared_ptr<BlockOPoints>>& blocks() const { return blocks_; }
    /// Set of spatially sieved blocks of points of a given atom, generated by sieve() internally. Only for blockscheme=atomic!
    const std::vector<std::vector<std::shared_ptr<BlockOPoints>>>& atomic_blocks() const {
        if (atomic_blocks_.empty()) {
            throw PSIEXCEPTION("MolecularGrid: No atomic blocks available. Wrong blockscheme?");
        } else {
            return atomic_blocks_;
        }
    }

    void set_debug(int debug) { debug_ = debug; }
};

class PseudospectralGrid : public MolecularGrid {
   protected:
    /// The primary basis
    std::shared_ptr<BasisSet> primary_;
    /// The filename used to optionally build the grid
    std::string filename_;

    /// The Options object
    Options& options_;

    /// Master builder methods
    void buildGridFromOptions();

   public:
    /// Constructor to use for autogeneration
    PseudospectralGrid(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> primary, Options& options);
    ~PseudospectralGrid() override;
};

class DFTGrid : public MolecularGrid {
   protected:
    /// The primary basis
    std::shared_ptr<BasisSet> primary_;
    /// Master builder methods
    void buildGridFromOptions(std::map<std::string, int> int_opts_map,
                              std::map<std::string, std::string> str_opts_map,
                              std::map<std::string, double> float_opts_map);
    /// The Options object
    Options& options_;

   public:
    std::shared_ptr<BasisSet> primary() const { return primary_; }

    DFTGrid(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> primary, Options& options);
    DFTGrid(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> primary,
            std::map<std::string, int> int_opts_map, std::map<std::string, std::string> str_opts_map, Options& options);
    DFTGrid(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> primary,
            std::map<std::string, int> int_opts_map, std::map<std::string, std::string> str_opts_map, 
            std::map<std::string, double> float_opts_map, Options& options);
    ~DFTGrid() override;
};

class RadialGrid {
   protected:
    /// Scheme
    std::string scheme_;
    /// Number of points in radial grid
    int npoints_;
    /// Alpha scale (for user reference)
    double alpha_;
    /// Nodes (including alpha)
    double* r_;
    /// Weights (including alpha and r^2)
    double* w_;

    // ==> Standard Radial Grids <== //

    /// Build the Becke 1988 radial grid
    static std::shared_ptr<RadialGrid> build_becke(int npoints, double alpha, int Z);
    /// Build the Treutler-Ahlrichs 1995 radial grid (scale power = 0.6)
    static std::shared_ptr<RadialGrid> build_treutler(int npoints, double alpha, int Z);
    // TODO: Add more grids

    /// Protected constructor
    RadialGrid();

   public:
    // ==> Initializers <== //

    /// Destructor
    virtual ~RadialGrid();

    /// Master build routine
    static std::shared_ptr<RadialGrid> build(const std::string& scheme, int npoints, double alpha, int Z);
    /// Hack build routine (TODO: Remove ASAP)
    static std::shared_ptr<RadialGrid> build(const std::string& scheme, int npoints, double* r, double* wr,
                                             double alpha, int Z);

    // ==> Accessors <== //

    /// Scheme of this radial grid
    std::string scheme() const { return scheme_; }
    /// Number of points in radial grid
    int npoints() const { return npoints_; }
    /// Alpha scale (for user reference)
    double alpha() const { return alpha_; }
    /// Radial nodes (including alpha scale). You do not own this.
    double* r() const { return r_; }
    /// Radial weights (including alpha scale and r^2). You do not own this.
    double* w() const { return w_; }

    /// Reflection
    void print(std::string out_fname = "outfile", int level = 1) const;
};

class SphericalGrid {
   protected:
    /// Scheme
    std::string scheme_;
    /// Number of points in radial grid
    int npoints_;
    /// Order of spherical harmonics in spherical grid (integrates products up to L_tot = 2 * order_ + 1)
    int order_;
    /// Spherical nodes, on the unit sphere
    double* x_;
    /// Spherical nodes, on the unit sphere
    double* y_;
    /// Spherical nodes, on the unit sphere
    double* z_;
    /// Spherical weights, normalized to 4pi
    double* w_;

    /// Spherical nodes, in spherical coordinates (azimuth)
    double* phi_;
    /// Spherical nodes, in spherical coordinates (inclination)
    double* theta_;

    // ==> Unique Lebedev Grids (statically stored) <== //

    /// Grid npoints to order map
    static std::map<int, int> lebedev_mapping_;

    /// Print valid Lebedev grids and error out (throws)
    static void lebedev_error();

    // ==> Utility Routines <== //

    /// Build the spherical angles from <x,y,z>, for reference
    void build_angles();

    /// Protected constructor
    SphericalGrid();

   public:
    // ==> Initializers <== //

    /// Destructor
    virtual ~SphericalGrid();

    /// Master build routines
    static std::shared_ptr<SphericalGrid> build(const std::string& scheme, int npoints, const MassPoint* points);

    // ==> Accessors <== //

    /// Scheme of this radial grid
    std::string scheme() const { return scheme_; }
    /// Number of points in radial grid
    int npoints() const { return npoints_; }
    /// Order of spherical harmonics in spherical grid (integrates products up to L_tot = 2 * order_ + 1)
    int order() const { return order_; }
    /// Spherical nodes, on the unit sphere
    double* x() const { return x_; }
    /// Spherical nodes, on the unit sphere
    double* y() const { return y_; }
    /// Spherical nodes, on the unit sphere
    double* z() const { return z_; }
    /// Spherical weights, normalized to 4pi
    double* w() const { return w_; }

    /// Spherical nodes, in spherical coordinates (azimuth)
    double* phi() const { return phi_; }
    /// Spherical nodes, in spherical coordinates (inclination)
    double* theta() const { return theta_; }

    /// Reflection
    void print(std::string out_fname = "outfile", int level = 1) const;
};

class BlockOPoints {
   protected:
    /// number of points in this block
    size_t index_;
    size_t npoints_;
    size_t local_nbf_;

    /// Data holders if requested
    SharedVector xvec_;
    SharedVector yvec_;
    SharedVector zvec_;
    SharedVector wvec_;

    /// Pointer to x (does not own)
    double* x_;
    /// Pointer to y (does not own)
    double* y_;
    /// Pointer to z (does not own)
    double* z_;
    /// Pointer to w (does not own)
    double* w_;
    /// Relevant shells, local -> global
    std::vector<int> shells_local_to_global_;
    /// Relevant functions, local -> global
    std::vector<int> functions_local_to_global_;
    /// Reference to the extents object
    std::shared_ptr<BasisExtents> extents_;

    /// Center of this BlockOPoints
    Vector3 xc_;
    /// Bounding radius of the BlockOPoints
    double R_;

    /// Parent atom of this BlockOPoints. -1 indicated none has been set.
    size_t parent_atom_ = -1;

    /// Populate significant functions given information in extents
    void populate();
    /// Compute bounding sphere
    void bound();

   public:
    BlockOPoints(SharedVector x, SharedVector y, SharedVector z, SharedVector w, std::shared_ptr<BasisExtents> extents);
    BlockOPoints(size_t index, size_t npoints, double* x, double* y, double* z, double* w,
                 std::shared_ptr<BasisExtents> extents);
    virtual ~BlockOPoints();

    /// Refresh populations (if extents_->delta() changes)
    void refresh() { populate(); }

    /// Number of grid points
    size_t npoints() const { return npoints_; }
    /// Number of basis functions in the block
    size_t local_nbf() const { return local_nbf_; }
    /// Index of the currently owned block
    size_t index() const { return index_; }
    /// Print a trace of this BlockOPoints
    void print(std::string out_fname = "outfile", int print = 2);
    /// Set parent atom for the current block
    void set_parent_atom(size_t atom) { parent_atom_ = atom; };
    /// Parent atom if the current block
    size_t parent_atom() const {
        if (parent_atom_ >= 0) {
            return parent_atom_;
        } else {
            throw PSIEXCEPTION("BlockOPoints: no parent atom set! Wrong blockscheme?");
        }
    };

    /// The x points. You do not own this
    double* x() const { return x_; }
    /// The y points. You do not own this
    double* y() const { return y_; }
    /// The z points. You do not own this
    double* z() const { return z_; }
    /// The weights. You do not own this
    double* w() const { return w_; }
    /// The center of the block
    Vector3 center() const { return xc_; }

    /// Relevant shells, local -> global
    const std::vector<int>& shells_local_to_global() const { return shells_local_to_global_; }
    /// Relevant functions, local -> global
    const std::vector<int>& functions_local_to_global() const { return functions_local_to_global_; }
};

class BasisExtents {
   protected:
    /// Basis this corresponds to
    std::shared_ptr<BasisSet> primary_;
    /// Cutoff value for basis values
    double delta_;
    /// Significant extent of shells
    std::shared_ptr<Vector> shell_extents_;
    /// Maximum extent
    double maxR_;

    /// Recompute and shell_extents_
    void computeExtents();

   public:
    BasisExtents(std::shared_ptr<BasisSet> primary, double delta);
    virtual ~BasisExtents();

    /// Print a trace of these extents
    void print(std::string out_fname = "outfile");
    /// Reset delta and recompute extents
    void set_delta(double delta) {
        delta_ = delta;
        computeExtents();
    }

    /// The cutoff value
    double delta() const { return delta_; }
    /// The basis set this BasisExtents is built on
    std::shared_ptr<BasisSet> basis() const { return primary_; }
    /// WCS significant extent of each shell
    std::shared_ptr<Vector> shell_extents() const { return shell_extents_; }
    /// Maximum spatial extent over all atoms
    double maxR() const { return maxR_; }
};
}  // namespace psi
#endif
