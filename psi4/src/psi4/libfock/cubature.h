/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
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

    static std::tuple<MolecularGrid::MolecularGridOptions,
                      std::map<std::string, int>,
                      std::map<std::string, std::string>,
                      std::map<std::string, double>> populateOptions(const Options& options,
                                const std::map<std::string, int>& int_opts_map = {},
                                const std::map<std::string, std::string>& str_opts_map = {},
                                const std::map<std::string, double>& float_opts_map = {});
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

// StandardGridMgr is used to build the SG-0 and SG-1 grids.
class StandardGridMgr {
   public:
    struct PruneGroup {
        short npts;
        short nreps;
    };
    struct PruneSpec {
        const PruneGroup *group;
        short nrad;     // Number of points for the radial grid
        short npts;     // Total number of points at all
        double rparam;  // Oftentimes the Bragg-Slater radius of the atom
    };

    // See S. Chien and P. Gill,  J. Comput. Chem. 27 (2006) 730-739.
    // This is Table 1 from that paper, with {0,0} terminators at the end of each row.
    //
    // We need an 18-point rule; this is not a Lebedev rule, but we
    // borrow it from Abramowitz & Stegun, p. 894 (the 5th-degree rule).
    static constexpr PruneGroup H__grp[] = {{6, 6},   {18, 3}, {26, 1}, {38, 1}, {74, 1}, {110, 1},
                                        {146, 6}, {86, 1}, {50, 1}, {38, 1}, {18, 1}, {0, 0}};
    static constexpr PruneGroup Li_grp[] = {{6, 6},   {18, 3}, {26, 1}, {38, 1}, {74, 1}, {110, 1},
                                        {146, 6}, {86, 1}, {50, 1}, {38, 1}, {18, 1}, {0, 0}};
    static constexpr PruneGroup Be_grp[] = {{6, 4},   {18, 2}, {26, 1}, {38, 2}, {74, 1}, {86, 1}, {110, 2},
                                        {146, 5}, {50, 1}, {38, 1}, {18, 1}, {6, 2},  {0, 0}};
    static constexpr PruneGroup B__grp[] = {{6, 4}, {26, 4}, {38, 3}, {86, 3}, {146, 6}, {38, 1}, {6, 2}, {0, 0}};
    static constexpr PruneGroup C__grp[] = {{6, 6},   {18, 2},  {26, 1},  {38, 2}, {50, 2}, {86, 1}, {110, 1},
                                        {146, 1}, {170, 2}, {146, 2}, {86, 1}, {38, 1}, {18, 1}, {0, 0}};
    static constexpr PruneGroup N__grp[] = {{6, 6},   {18, 3},  {26, 1}, {38, 2}, {74, 2}, {110, 1},
                                        {170, 2}, {146, 3}, {86, 1}, {50, 2}, {0, 0}};
    static constexpr PruneGroup O__grp[] = {{6, 5},   {18, 1}, {26, 2}, {38, 1}, {50, 4}, {86, 1},
                                        {110, 5}, {86, 1}, {50, 1}, {38, 1}, {6, 1},  {0, 0}};
    static constexpr PruneGroup F__grp[] = {{6, 4},   {38, 2}, {50, 4}, {74, 2}, {110, 2}, {146, 2},
                                        {110, 2}, {86, 3}, {50, 1}, {6, 1},  {0, 0}};
    static constexpr PruneGroup Na_grp[] = {{6, 6}, {18, 2}, {26, 3}, {38, 1}, {50, 2}, {110, 8}, {74, 2}, {6, 2}, {0, 0}};
    static constexpr PruneGroup Mg_grp[] = {{6, 5},   {18, 2},  {26, 2}, {38, 2}, {50, 2}, {74, 1}, {110, 2},
                                        {146, 4}, {110, 1}, {86, 1}, {38, 2}, {18, 1}, {6, 1},  {0, 0}};
    static constexpr PruneGroup Al_grp[] = {{6, 6},   {18, 2},  {26, 1}, {38, 2}, {50, 2}, {74, 1}, {86, 1}, {146, 2},
                                        {170, 2}, {110, 2}, {86, 1}, {74, 1}, {26, 1}, {18, 1}, {6, 1},  {0, 0}};
    static constexpr PruneGroup Si_grp[] = {{6, 5},   {18, 4},  {38, 4}, {50, 3}, {74, 1}, {110, 2},
                                        {146, 1}, {170, 3}, {86, 1}, {50, 1}, {6, 1},  {0, 0}};
    static constexpr PruneGroup P__grp[] = {{6, 5},   {18, 4},  {38, 4}, {50, 3}, {74, 1}, {110, 2},
                                        {146, 1}, {170, 3}, {86, 1}, {50, 1}, {6, 1},  {0, 0}};
    static constexpr PruneGroup S__grp[] = {{6, 4},   {18, 1},  {26, 8},  {38, 2}, {50, 1}, {74, 2}, {110, 1},
                                        {170, 3}, {146, 1}, {110, 1}, {50, 1}, {6, 1},  {0, 0}};
    static constexpr PruneGroup Cl_grp[] = {{6, 4},   {18, 7},  {26, 2},  {38, 2}, {50, 1}, {74, 1}, {110, 2},
                                        {170, 3}, {146, 1}, {110, 1}, {86, 1}, {6, 1},  {0, 0}};

    // clang-format off
    static constexpr PruneSpec SG0specs[18] = {
        {nullptr,  0,    0,    0},
        {H__grp,  23, 1406, 1.30},
        {nullptr,  0,    0,    0},
        {Li_grp,  23, 1406, 1.95},
        {Be_grp,  23, 1390, 2.20},
        {B__grp,  23, 1426, 1.45},
        {C__grp,  23, 1390, 1.20},
        {N__grp,  23, 1414, 1.10},
        {O__grp,  23, 1154, 1.10},
        {F__grp,  23, 1494, 1.20},
        {nullptr,  0,    0,    0},
        {Na_grp,  26, 1328, 2.30},
        {Mg_grp,  26, 1468, 2.20},  // Warning: The original paper has 1492 points instead of 1468...
        {Al_grp,  26, 1496, 2.10},
        {Si_grp,  26, 1496, 1.30},
        {P__grp,  26, 1496, 1.30},
        {S__grp,  26, 1456, 1.10},
        {Cl_grp,  26, 1480, 1.45}};
    // clang-format on

    // rows[Z] tells you on which row of the periodic table you can find element Z.
    // clang-format off
    static constexpr short rows[] = {
        0,
        1,                      1,
        2, 2,    2, 2, 2, 2, 2, 2,
        3, 3,    3, 3, 3, 3, 3, 3,
    };
    // clang-format on

    //
    // The paper also doesn't provide any guidance on how to handle the elements after Argon.

    // These radii are taken from Table 1 of Gill, Johnson, and Pople.
    // clang-format off
    static constexpr double SG1radii[] = {
        0,
        1.0000,                                                     0.5882,
        3.0769, 2.0513,     1.5385, 1.2308, 1.0256, 0.8791, 0.7692, 0.6838,
        4.0909, 3.1579,     2.5714, 2.1687, 1.8750, 1.6514, 1.4754, 1.3333};
    // clang-format on

    // These angular point counts were selected for compatibility with Q-Chem.
    static constexpr PruneGroup row1spec[] = {{6, 16}, {38, 5}, {86, 4}, {194, 9}, {86, 16}, {0, 0}};  // H and He
    static constexpr PruneGroup row2spec[] = {{6, 14}, {38, 7}, {86, 3}, {194, 9}, {86, 17}, {0, 0}};  // Li to Ne
    static constexpr PruneGroup row3spec[] = {{6, 12}, {38, 7}, {86, 5}, {194, 7}, {86, 19}, {0, 0}};  // Na to Ar
    static constexpr PruneGroup row4spec[] = {{38, 12}, {194, 38}, {0, 0}};  // Everything after Argon

    static constexpr PruneSpec SG1specs[] = {
        {row1spec, 50, 3752, 0}, {row2spec, 50, 3816, 0}, {row3spec, 50, 3760, 0}, {row4spec, 50, 7828, 0}};

    static const MassPoint *SG0_grids_[18];
    static int SG0_sizes_[18];

    static const MassPoint *SG1_grids_[19];
    static int SG1_sizes_[19];

   private:
    static void makeCubatureGridFromPruneSpec(PruneSpec const &spec, int radscheme, MassPoint *grid_out);
    static void Initialize_SG0();
    static void Initialize_SG1();

   public:
    static void Initialize();
    static void ReleaseMemory();
    static int WhichGrid(const char *name);
    static int GetSG0size(int Z);
    static int GetSG1size(int Z);
    static const MassPoint *GetSG0grid(int Z);
    static const MassPoint *GetSG1grid(int Z);
};

class OrientationMgr {
    // "Local" vector, matrix, atom, and molecule definitions.
    // It makes some of the code easier to read.
    struct LVector {
        double x, y, z;
    };

    struct LMatrix {
        double xx, xy, xz;
        double yx, yy, yz;
        double zx, zy, zz;
    };

    ///// These are the only member variables in the whole class! /////
    std::shared_ptr<Molecule> molecule_;
    LMatrix rotation_;
    ///// Everything else is to set these up /////

    struct LAtom {
        LVector pos;  // This position is ALWAYS relative to the center of charge.
        int atomicNumber;
    };

    typedef LAtom LMolecule[];  // An "LMolecule" is an array of LAtoms.

    static const LMatrix lIdentityMatrix;  // = {1,0,0, 0,1,0, 0,0,1};

    enum SymmetryType {
        SINGLE_ATOM,
        LINEAR,
        ICOSAHEDRAL,
        OCTAHEDRAL,
        TETRAHEDRAL,
        SPHERICAL_MYSTERY,
        SYMMETRIC_TOP,
        ASYMMETRIC
    };
    static const char *symmetryNames[];  // symmetryNames must match SymmetryType!

    // Some handy vector operations
    static inline double vdot(LVector v, LVector w) { return v.x * w.x + v.y * w.y + v.z * w.z; }
    static inline double vnorm(LVector v) { return sqrt(vdot(v, v)); }
    static inline LVector vadd(LVector v, LVector w) {
        LVector sum = {v.x + w.x, v.y + w.y, v.z + w.z};
        return sum;
    }
    static inline LVector vsub(LVector v, LVector w) {
        LVector diff = {v.x - w.x, v.y - w.y, v.z - w.z};
        return diff;
    }
    static inline LVector vcross(LVector v, LVector w) {
        LVector cross = {v.y * w.z - v.z * w.y, v.z * w.x - v.x * w.z, v.x * w.y - v.y * w.x};
        return cross;
    }
    static inline LVector vnormalize(LVector v) {
        double norm = vnorm(v);
        LVector u = {v.x / norm, v.y / norm, v.z / norm};
        return u;
    }
    static inline LVector normalToTriangle(LVector a, LVector b, LVector c) { return vcross(vsub(b, a), vsub(c, a)); }

    static inline LVector mvtimes(const LMatrix &Q, const LVector &v) {
        LVector out = {Q.xx * v.x + Q.xy * v.y + Q.xz * v.z, Q.yx * v.x + Q.yy * v.y + Q.yz * v.z,
                       Q.zx * v.x + Q.zy * v.y + Q.zz * v.z};
        return out;
    }

    static inline LMatrix mmtimes(const LMatrix &A, const LMatrix &B) {
        LMatrix C;
        C.xx = A.xx * B.xx + A.xy * B.yx + A.xz * B.zx;
        C.xy = A.xx * B.xy + A.xy * B.yy + A.xz * B.zy;
        C.xz = A.xx * B.xz + A.xy * B.yz + A.xz * B.zz;
        C.yx = A.yx * B.xx + A.yy * B.yx + A.yz * B.zx;
        C.yy = A.yx * B.xy + A.yy * B.yy + A.yz * B.zy;
        C.yz = A.yx * B.xz + A.yy * B.yz + A.yz * B.zz;
        C.zx = A.zx * B.xx + A.zy * B.yx + A.zz * B.zx;
        C.zy = A.zx * B.xy + A.zy * B.yy + A.zz * B.zy;
        C.zz = A.zx * B.xz + A.zy * B.yz + A.zz * B.zz;
        return C;
    }

    static inline LMatrix cycleXtoZ(const LMatrix &in) {
        // We want to swap the x-row with the z-row, but that
        // can interfere with chirality, so we instead rotate
        // x->z, z->x, x->y.
        LMatrix out;
        out.xx = in.yx;
        out.xy = in.yy;
        out.xz = in.yz;
        out.yx = in.zx;
        out.yy = in.zy;
        out.yz = in.zz;
        out.zx = in.xx;
        out.zy = in.xy;
        out.zz = in.xz;
        return out;
    }

// Floating-point comparisons
#define EPSILON 1e-10
    static inline bool fequ(double a, double b) { return std::fabs(a - b) < EPSILON; }
    static inline bool fless(double a, double b) {
        return a < b + EPSILON;
    }  // We are assured that fequ(a, b) implies !fless(a, b)
    static inline bool fequ2(double a, double b) { return std::fabs(a - b) < EPSILON * EPSILON; }
    static inline bool fvIsZero(LVector v) { return fequ(v.x, 0) && fequ(v.y, 0) && fequ(v.z, 0); }
    static inline bool fvequ(LVector v, LVector w) {
        return std::fabs(v.x - w.x) < EPSILON && std::fabs(v.y - w.y) < EPSILON && std::fabs(v.z - w.z) < EPSILON;
    }
    static inline bool fperp(LVector v, LVector w) { return fequ(vdot(v, w), 0); }
#undef EPSILON

    static inline LVector rotate_xy(LVector v, double theta) {
        double s = sin(theta), c = cos(theta);
        LVector ans = {v.x * c - v.y * s, v.y * s + v.y * c, v.z};
        return ans;
    }

    // Reflect a vector across a vertical plane of angle `theta' with the x-y axis...
    static inline LVector reflect_xy(LVector v, double theta) {
        double s = sin(2 * theta), c = cos(2 * theta);
        LVector ans = {v.x * c + v.y * s, v.x * s - v.y * c, v.z};
        return ans;
    }

    static LMatrix RotMatrixFromTwoAxes(LVector ax1, LVector ax2);
    static LVector someUnitVectorPerpendicularTo(LVector v);
    static LMatrix RotMatrixFromOneAxis(LVector axis3);
    static void diagonalize(LMatrix const &M, LMatrix *Q_out, LVector *D_out);
    static LMatrix buildRotationMatrix(LVector const &axis, double angle);

    static bool isAnAtomLocatedAt(LMolecule mol, int natom, LVector const &position, int true_atomic_number);
    static bool bothAnglesResultInEquivalentGrids(LMolecule mol, int natoms, double angle1, double angle2);
    static LMatrix symmetricTopMatrix(std::shared_ptr<Molecule> mol, LMatrix const &Q, LVector const &center);
    static bool TestAxis(LMolecule mol, int natom, LVector const &axis, int n);

    static bool LookForIcosahedralSymmetry(LMolecule mol, int natom, LMatrix *Q_out);
    static bool LookForOctahedralSymmetry(LMolecule mol, int natom, LMatrix *Q_out);
    static bool LookForTetrahedralSymmetry(LMolecule mol, int natom, LMatrix *Q_out);
    static SymmetryType sphericalTopMatrix(std::shared_ptr<Molecule> origmol, LVector const &center, LMatrix *Q_out);

   public:
    // Use the constructor to determine the standard orientation, and use `MoveIntoPosition' to apply it to a
    // mass-point.
    OrientationMgr(std::shared_ptr<Molecule> mol);
    inline MassPoint MoveIntoPosition(MassPoint mp, int A) {
        LVector oldpos = {mp.x, mp.y, mp.z};
        LVector rotated = mvtimes(rotation_, oldpos);
        LVector atompos = {molecule_->x(A), molecule_->y(A), molecule_->z(A)};
        LVector newpos = vadd(rotated, atompos);
        MassPoint newmp = {newpos.x, newpos.y, newpos.z, mp.w};
        return newmp;
    }

    // Or we could use Matrix like smart people
    std::shared_ptr<Matrix> orientation() const {
        auto O = std::make_shared<Matrix>("O", 3, 3);
        double **Op = O->pointer();
        Op[0][0] = rotation_.xx;
        Op[0][1] = rotation_.xy;
        Op[0][2] = rotation_.xz;
        Op[1][0] = rotation_.yx;
        Op[1][1] = rotation_.yy;
        Op[1][2] = rotation_.yz;
        Op[2][0] = rotation_.zx;
        Op[2][1] = rotation_.zy;
        Op[2][2] = rotation_.zz;
        return O;
    }
};

class RadialPruneMgr {
   private:
    int nominal_order_;
    double alpha_;
    double (*pruneFn_)(double, double);

    // These are ad-hoc functions. The idea is that rho = (distance from center of atom)/(Bragg-Slater radius) and alpha
    // = some settable parameter.
    static double flat(double /*rho*/, double /*alpha*/) { return 1; }
    static double p_slater(double rho, double /*alpha*/) { return rho * exp(1 - rho); }
    static double d_slater(double rho, double alpha) {
        double pslater = p_slater(rho, alpha);
        return pslater * pslater;
    }
    static double log_slater(double rho, double alpha) { return exp(-alpha * std::fabs(log(rho))); }
    static double p_gaussian(double rho, double /*alpha*/) {
        return rho * exp((1 - rho * rho) / 2.0);
    }  // Note: The original implementation had (1 - R*rho) instead of (1 - rho*rho)
    static double d_gaussian(double rho, double /*alpha*/) { return rho * rho * exp(1 - rho * rho); }
    static double log_gaussian(double rho, double alpha) { return exp(-alpha * log(rho) * log(rho)); }

    struct PruneFunctionTable {
        const char *name;
        double (*scalFn)(double, double);
    };
    static PruneFunctionTable prunefunctions[];

   public:
    static int WhichPruneFunction(const char *functionname);
    static const char *FunctionName(int which) { return prunefunctions[which].name; }
    RadialPruneMgr(MolecularGrid::MolecularGridOptions const &opt);
    int GetPrunedNumAngPts(double rho);
    int ShellPruning(int ri, int Z, int radial_pts);
    int TreutlerShellPruning(int ri, int Z, int radial_pts);
};

class RadialGridMgr {
   private:
    // See P. Gill and S. Chien,  J. Comput. Chem. 24 (2003) 732-740
    static double multiexp_r(double x) { return -log(x); }
    static double multiexp_dr(double x) {
        double logx = log(x);
        return 1 / (logx * logx * x);
    }

    // See C.W. Murray, N.C. Handy, G.J. Laming,  Mol Phys 78 (1993) 997.
    // Good luck procuring a copy.
    static double em_r(double x) {
        double zzz = x / (1 - x);
        return zzz * zzz;
    }
    static double em_dr(double x) {
        double onemx = 1 - x;
        return 2 * x / (onemx * onemx * onemx);
    }

    // See A. D. Becke, J. Chem. Phys. 88 (1988) 2547
    static double becke_r(double x) { return (1 - x) / (1 + x); }
    static double becke_dr(double x) { return 2 / ((1 + x) * (1 + x)); }

    // See M. E. Mura and P. J. Knowles, J. Chem. Phys. 104 (1996) 9848
    // They have this annoying `alpha' prefactor that is five for some atoms and seven for others.
    static double mura5_r(double x) {
        double alpha = 5;
        return -alpha * log(1 - x * x * x);
    }
    static double mura5_dr(double x) {
        double alpha = 5;
        return 3 * alpha * x * x / (1 - x * x * x);
    }
    static double mura7_r(double x) {
        double alpha = 7;
        return -alpha * log(1 - x * x * x);
    }
    static double mura7_dr(double x) {
        double alpha = 7;
        return 3 * alpha * x * x / (1 - x * x * x);
    }
    // If Z corresponds to an element in group one or two, alpha = 7. Otherwise it's five.
    static bool muraKnowlesAlphaIsReallySeven(int Z) {
        return Z == 3 || Z == 4 || Z == 11 || Z == 12 || Z == 19 || Z == 20;
    }

// See O. Treutler and R. Ahlrichs, J. Chem. Phys. 102 (1995) 346
#define INVLN2 1.4426950408889634074  // = 1/log(2)
    static double ahlrichs_r(double x) {
        double alpha = 0.6;
        return -INVLN2 * pow(1 + x, alpha) * log((1 - x) / 2);
    }
    static double ahlrichs_dr(double x) {
        double alpha = 0.6;
        return INVLN2 * pow(1 + x, alpha) * (-alpha * log((1 - x) / 2) / (1 + x) + 1 / (1 - x));
    }
#undef INVLN2

    static void getTrapezoidalRoots(int n, double r[], double w[]);
    static void getChebychevRootsKind2(int n, double r[], double w[]);
    static void getLegendreRoots(
        int n, double r[],
        double w[]);  // We don't actually need this function, but I don't have the heart to delete it.
    static void getLaguerreRoots(int n, double r[], double w[]);
    static void getMultiExpRoots(int n, double r[], double w[]);

    static double maxRowSumNorm(int n, double a[], double b[]);
    static void GolombWelsch(int n, double a[], double b[], double q[]);

    struct SchemeTable {
        const char *name;
        void (*getRoots)(int n, double r[], double w[]);
        double (*rFn)(double x);
        double (*drdxFn)(double x);
    };
    static SchemeTable radialschemes[];

   public:
    static int WhichScheme(const char *schemename);
    static const char *SchemeName(int which) { return radialschemes[which].name; }
    static int MuraKnowlesHack(int scheme, int Z);
    static void makeRadialGrid(int n, int whichScheme, double r[], double w[], double Rparam);
};

double GetBSRadius(unsigned Z);

// clang-format off
// LebedevGridMgr is a static class---all the members are static.
// These functions don't necessarily *belong* in a class, but
// we put them there for philosophical encapsulation reasons.
class LebedevGridMgr {
   public:
    static void Initialize();
    // If you know the number of points in the grid you want, call this.
    static const MassPoint *findGridByNPoints(int npoints);
    static int findOrderByNPoints(int npoints);

    // If you know the order of the grid you want, call these.
    static bool isUsableOrder(int order);
    static const MassPoint *findGridByOrder(int order);
    static const MassPoint *findGridByOrder_roundUp(int order);
    static int findNPointsByOrder(int order);
    static int findNPointsByOrder_roundUp(int order);

    static void PrintHelp();
    static const int MaxOrder = 131;

   private:
    LebedevGridMgr();
    static inline MassPoint MASSPOINT(double x, double y, double z, double w);
    static int addPoints1(MassPoint point[], double v);
    static int addPoints2(MassPoint point[], double v);
    static int addPoints3(MassPoint point[], double v);
    static int addPoints4(MassPoint point[], double v, double a);
    static int addPoints5(MassPoint point[], double v, double a);
    static int addPoints6(MassPoint point[], double v, double a, double b);

    static const MassPoint *mk1ptGrid();
    static const MassPoint *mk6ptGrid();
    static const MassPoint *mk14ptGrid();
    static const MassPoint *mk18ptGrid_nonstandard();
    static const MassPoint *mk26ptGrid();
    static const MassPoint *mk38ptGrid();
    static const MassPoint *mk50ptGrid();
    static const MassPoint *mk74ptGrid();
    static const MassPoint *mk86ptGrid();
    static const MassPoint *mk110ptGrid();
    static const MassPoint *mk146ptGrid();
    static const MassPoint *mk170ptGrid();
    static const MassPoint *mk194ptGrid();
    static const MassPoint *mk230ptGrid();
    static const MassPoint *mk266ptGrid();
    static const MassPoint *mk302ptGrid();
    static const MassPoint *mk350ptGrid();
    static const MassPoint *mk434ptGrid();
    static const MassPoint *mk590ptGrid();
    static const MassPoint *mk770ptGrid();
    static const MassPoint *mk974ptGrid();
    static const MassPoint *mk1202ptGrid();
    static const MassPoint *mk1454ptGrid();
    static const MassPoint *mk1730ptGrid();
    static const MassPoint *mk2030ptGrid();
    static const MassPoint *mk2354ptGrid();
    static const MassPoint *mk2702ptGrid();
    static const MassPoint *mk3074ptGrid();
    static const MassPoint *mk3470ptGrid();
    static const MassPoint *mk3890ptGrid();
    static const MassPoint *mk4334ptGrid();
    static const MassPoint *mk4802ptGrid();
    static const MassPoint *mk5294ptGrid();
    static const MassPoint *mk5810ptGrid();

    static const MassPoint *nonstandard18PointGrid_;
    struct GridData {
        int order;
        int npoints;
        MassPoint const *(*mkGridFn)();
        MassPoint const *grid;
    };
    static GridData grids_[];
};
// clang-format on

}  // namespace psi
#endif
