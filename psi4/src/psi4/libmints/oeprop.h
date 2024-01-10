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

#ifndef _psi_src_lib_oeprop_h
#define _psi_src_lib_oeprop_h

#include <map>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include "typedefs.h"
#include "psi4/libmints/vector3.h"

namespace psi {

class Wavefunction;
class IntegralFactory;
class MatrixFactory;
class BasisSet;

/**
 * The Prop object, base class of OEProp and GridProp objects
 *
 *  Wavefunction is not finalized, so we have no idea what bases
 *  general density matrices/orbital coefficients will be in.
 *  Additionally, there are questions of natural/local/canonical orbitals,
 *  relaxed/unrelaxed OPDM, and possible special densities such as needed
 *  for transition dipoles.
 *
 *  Therefore, the Prop object explicitly stores a mutable set of D/C
 *  matrices in the SO basis, and knows how to transform between them at will
 *  Additionally, while these are initially ripped from the constructor's
 *  wavefunction, substitutions may be made later, to use different orbitals
 *  or densities.
 */
class PSI_API Prop {
   protected:
    /// The wavefunction object this Prop is built around
    std::shared_ptr<Wavefunction> wfn_;
    /// The basisset for this wavefunction
    std::shared_ptr<BasisSet> basisset_;
    /// Is this wavefunction object spin-restricted? (Actually closed-shell, but this is wavefunction's convention)
    bool same_orbs_;  // This allows pointers to be duplicated/computation skipped
    bool same_dens_;  // This allows pointers to be duplicated/computation skipped
    /// The integral factory for this wavefunction's basisset
    std::shared_ptr<IntegralFactory> integral_;
    /// The matrix factory for this wavefunction's basisset (SO)
    std::shared_ptr<MatrixFactory> factory_;

    /// The AO to USO matrix
    SharedMatrix AO2USO_;

    /**
     * Internally, data is held in the SO basis in Pitzer order
     */

    /// The alpha eigenvalues in the MO basis (used to form Pitzer ordering)
    SharedVector epsilon_a_;
    /// The alpha eigenvalues in the MO basis (used to form Pitzer ordering)
    SharedVector epsilon_b_;
    /// The alpha density matrix in the SO basis
    SharedMatrix Da_so_;
    /// The beta density matrix in the SO basis
    SharedMatrix Db_so_;
    /// The alpha C matrix in the SO basis
    SharedMatrix Ca_so_;
    /// The beta C matrix in the SO basis
    SharedMatrix Cb_so_;

    /// Common initialization
    void common_init();

   public:
    /// Build a Prop object with C, epsilon, and restricted buit from wfn
    Prop(std::shared_ptr<Wavefunction> wfn);
    /// Virtual destructor
    virtual ~Prop();

    // => Wavefunction Modifiers (rarely called, C is usually fixed at HF) <= //

    // Change restricted flag. Resets C/D/epsilon matrices from wfn
    void set_wavefunction(std::shared_ptr<Wavefunction> wfn);
    // Change restricted flag. Resets C/D/epsilon matrices from wfn
    void set_restricted(bool restricted);
    // Set alpha eigenvalues, MO pitzer order basis
    void set_epsilon_a(SharedVector epsilon_a);
    // Set beta eigenvalues, MO pitzer order basis. Throws if restricted
    void set_epsilon_b(SharedVector epsilon_a);
    // Set alpha C matrix, SO/MO pitzer order basis.
    void set_Ca(SharedMatrix Ca);
    // Set beta C matrix, SO/MO pitzer order basis. Throws if restricted
    void set_Cb(SharedMatrix Cb);

    // => Set OPDM/TDM/DDM (often called). These need not be totally symmetric. Note, you are setting Da and/or Db, I do
    // the adding to Dt  <= //

    // TODO Add symmetry is irrep number
    void set_Da_ao(SharedMatrix Da, int symmetry = 0);
    void set_Db_ao(SharedMatrix Db, int symmetry = 0);
    void set_Da_so(SharedMatrix Da);
    void set_Db_so(SharedMatrix Db);
    void set_Da_mo(SharedMatrix Da);
    void set_Db_mo(SharedMatrix Db);

    // => Get routines (useful to quickly change bases) <= //

    /// The alpha eigenvalues in the MO basis (used to form Pitzer ordering)
    SharedVector epsilon_a();
    /// The alpha eigenvalues in the MO basis (used to form Pitzer ordering)
    SharedVector epsilon_b();

    /// The alpha C matrix in the SO basis
    SharedMatrix Ca_so();
    /// The beta C matrix in the SO basis
    SharedMatrix Cb_so();
    /// The alpha C matrix in the AO (Spherical Harmonics, C1) basis. Ordered by eigenvalue
    SharedMatrix Ca_ao();
    /// The beta C matrix in the AO (Spherical Harmonics, C1) basis. Ordered by eigenvalue
    SharedMatrix Cb_ao();

    /// The alpha density matrix in the AO (Spherical Harmonics, C1) basis
    SharedMatrix Da_ao();
    /// The beta density matrix in the AO (Spherical Harmonics, C1) basis
    SharedMatrix Db_ao();
    /// The alpha density matrix in the SO basis
    SharedMatrix Da_so();
    /// The beta density matrix in the SO basis
    SharedMatrix Db_so();
    /// The alpha density matrix in the MO basis
    SharedMatrix Da_mo();
    /// The beta density matrix in the MO basis
    SharedMatrix Db_mo();

    /// The total/spin density matrix in the ao basis, depending on if true or false
    SharedMatrix Dt_so(bool total = true);
    /// The total/spin density matrix in the ao basis, depending on if true or false
    SharedMatrix Dt_mo(bool total = true);

    /// The alpha natural orbital occupations and orbitals in the MO basis
    std::pair<SharedMatrix, SharedVector> Na_mo();
    /// The beta natural orbital occupations and orbitals in the MO basis. Throws if restricted
    std::pair<SharedMatrix, SharedVector> Nb_mo();
    /// The total natural orbital occupations and orbitals in the MO basis
    std::pair<SharedMatrix, SharedVector> Nt_mo();
    /// The alpha natural orbital occupations and orbitals in the SO basis
    std::pair<SharedMatrix, SharedVector> Na_so();
    /// The beta natural orbital occupations and orbitals in the SO basis. Throws if restricted
    std::pair<SharedMatrix, SharedVector> Nb_so();
    /// The total natural orbital occupations and orbitals in the SO basis
    std::pair<SharedMatrix, SharedVector> Nt_so();
    /// The alpha natural orbital occupations and orbitals in the AO basis
    std::pair<SharedMatrix, SharedVector> Na_ao();
    /// The beta natural orbital occupations and orbitals in the AO basis. Throws if restricted
    std::pair<SharedMatrix, SharedVector> Nb_ao();
    /// The total natural orbital occupations and orbitals in the AO basis
    std::pair<SharedMatrix, SharedVector> Nt_ao();

    /// Density Matrix title, used for fallback naming of OEProp compute jobs
    std::string Da_name() const;
    /// Density Matrix title, used for fallback naming of OEProp compute jobs
    std::string Db_name() const;

    // => Some integral helpers <= //
    SharedMatrix overlap_so();
};

/**
 * MultipolePropCalc
 *
 * Class, which calculates multipoles and mo_extents.
 *
 * Historically this class was part of OEProp.
 *
 * It is initialized with a wavefunction and an origin vector. If the origin breaks symmetry,
 * a warning is generated. Apart from this, this class does not have output, it also does
 * not export any values into the environment.
 *
 * If you are looking for previous OEProp functionality (i.e. output and environment exports)
 * OEProp still contains all that.
 *
 */

class MultipolePropCalc : public Prop {
   private:
    MultipolePropCalc();

   protected:
    /// The center about which properties are computed
    Vector3 origin_;
    /// Whether the origin is on a symmetry axis or not
    bool origin_preserves_symmetry_;

   public:
    /// Common initialization
    MultipolePropCalc(std::shared_ptr<Wavefunction> wfn, Vector3 const& origin);
    // Output Type of multipole function: name, elec, nuc, tot, order
    typedef std::vector<std::tuple<std::string, double, double, double, int>> MultipoleOutputTypeBase;
    typedef std::shared_ptr<MultipoleOutputTypeBase> MultipoleOutputType;
    /// Compute arbitrary-order multipoles up to (and including) l=order. returns name, elec, nuc and tot as vector_ptr
    MultipoleOutputType compute_multipoles(int order, bool transition = false, bool print_output = false,
                                           bool verbose = false);
    /// Compute mo extents
    std::vector<SharedVector> compute_mo_extents(bool print_output = false);
};

/**
 * PopulationAnalysisCalc
 *
 * Class, which carries out popular population analysis, such as Mulliken or Loewdin.
 *
 * Historically this class was part of OEProp.
 *
 * It is initialized with a wavefunction. It does not generate any output or populate any
 * environment variables.
 *
 * If you are looking for previous OEProp functionality (i.e. output and environment exports)
 * OEProp still contains all that.
 *
 */

class PopulationAnalysisCalc : public Prop {
   private:
    PopulationAnalysisCalc();

   public:
    typedef std::shared_ptr<std::vector<double>> SharedStdVector;
    PopulationAnalysisCalc(std::shared_ptr<Wavefunction> wfn);
    ~PopulationAnalysisCalc() override;
    /// Compute Mulliken Charges
    std::tuple<SharedStdVector, SharedStdVector, SharedStdVector> compute_mulliken_charges(bool print_output = false);
    /// Compute Lowdin Charges
    std::tuple<SharedStdVector, SharedStdVector, SharedStdVector> compute_lowdin_charges(bool print_output = false);
    /// Compute MBIS Multipoles (doi:10.1021/acs.jctc.6b00456)
    std::tuple<SharedMatrix, SharedMatrix, SharedMatrix, SharedMatrix> compute_mbis_multipoles(bool free_atom_volumes = false, bool print_output = false);
    /// Compute Mayer Bond Indices (non-orthogoal basis)
    std::tuple<SharedMatrix, SharedMatrix, SharedMatrix, SharedVector> compute_mayer_indices(bool print_output = false);
    /// Compute Wiberg Bond Indices using Lowdin Orbitals (symmetrically orthogonal basis)
    std::tuple<SharedMatrix, SharedMatrix, SharedMatrix, SharedVector> compute_wiberg_lowdin_indices(
        bool print_output = false);
    /// Compute/display natural orbital occupations around the bandgap. Displays max_num above and below the bandgap
    std::shared_ptr<std::vector<std::vector<std::tuple<double, int, int>>>> compute_no_occupations(
        int max_noon = 3, bool print_output = false);
};

/**
 * ESPPropCalc
 *
 * Class, which calculates multiple potentials based on a grid.
 *
 * Historically this class was part of OEProp.
 *
 * It is initialized with a wavefunction.
 * environment variables. Most functions currently require a file "grid.dat" to be present.
 * The functions generate their output on disk in files such as grid-esp.dat.
 * No environment variables are populated.
 *
 * If you are looking for previous OEProp functionality (i.e. output and environment exports)
 * OEProp still contains all that.
 *
 */
class ESPPropCalc : public Prop {
   private:
    // Constructing without wavefunction is forbidden:
    ESPPropCalc();

   protected:
    /// The ESP in a.u., computed at each grid point
    std::vector<double> Vvals_;
    /// The field components in a.u. computed at each grid point
    std::vector<double> Exvals_;
    std::vector<double> Eyvals_;
    std::vector<double> Ezvals_;

   public:
    /// Constructor
    ESPPropCalc(std::shared_ptr<Wavefunction> wfn);
    /// Destructor
    ~ESPPropCalc() override;

    std::vector<double> const& Vvals() const { return Vvals_; }
    std::vector<double> const& Exvals() const { return Exvals_; }
    std::vector<double> const& Eyvals() const { return Eyvals_; }
    std::vector<double> const& Ezvals() const { return Ezvals_; }

    /// This function is missing, it should be removed until it is implemented.
    void compute_electric_field_and_gradients();
    /// Compute electrostatic potentials at the nuclei
    std::shared_ptr<std::vector<double>> compute_esp_at_nuclei(bool print_output = false, bool verbose = false);
    /// Compute electrostatic potential at specified grid points
    void compute_esp_over_grid(bool print_output = false);
    /// Compute field at specified grid points
    void compute_field_over_grid(bool print_output = false);
    /// Compute electrostatic potential at grid points based on input grid, OpenMP version. input_grid is Nx3
    SharedVector compute_esp_over_grid_in_memory(SharedMatrix input_grid) const;
    /// Compute field at grid points based on input grid, OpenMP version. input_grid is Nx3
    SharedMatrix compute_field_over_grid_in_memory(SharedMatrix input_grid) const;
};

/**
 * The OEProp object, computes arbitrary expectation values (scalars)
 * analyses (typically vectors)
 **/
class PSI_API OEProp {
   private:
    /// Constructor, uses globals and Process::environment::reference wavefunction, Implementation does not exist.
    OEProp();

   protected:
    /// The wavefunction object this Prop is built around
    std::shared_ptr<Wavefunction> wfn_;
    /// Print flag
    int print_;
    /// OEProp name, for printout purposes.
    /// TODO: Standardize density matrix names across Psi so we can remove `title_`.
    ///   `title_` is redundant. Ideally, we'd print the density matrix names and not bother naming the
    ///   OEProp object itself. But if a user sets an MO density matrix and gives it a name asserting
    ///   it is an MO density matrix, OEProp would need to recognize that and rename the SO basis matrix
    ///   it creates. Without a standard naming scheme for densities, OEProp can't do that. Therefore,
    ///   until density names are standardized, we need to awkwardly keep `title_` around.
    std::string title_;
    /// Variable names, for variable saving purposes. Each should have a single '{}' substring to indicate
    /// where the variable name goes. The full generality of format strings are needed for excited
    /// state purposes.
    std::unordered_set<std::string> names_;
    /// The set of tasks to complete
    std::set<std::string> tasks_;
    /// Common initialization
    void common_init();
    /// Print header and information
    void print_header();

    // Compute routines
    /// Compute arbitrary-order multipoles up to (and including) l=order
    void compute_multipoles(int order, bool transition = false);
    /// Compute mo extents
    void compute_mo_extents();
    /// Compute Mulliken Charges
    void compute_mulliken_charges();
    /// Compute Lowdin Charges
    void compute_lowdin_charges();
    /// Compute MBIS Multipoles (doi:10.1021/acs.jctc.6b00456)
    void compute_mbis_multipoles(bool free_atom_volumes = false);
    /// Compute Mayer Bond Indices (non-orthogoal basis)
    void compute_mayer_indices();
    /// Compute Wiberg Bond Indices using Lowdin Orbitals (symmetrically orthogonal basis)
    void compute_wiberg_lowdin_indices();
    /// Compute/display natural orbital occupations around the bandgap. Displays max_num above and below the bandgap
    void compute_no_occupations();
    /// Compute electric field and electric field gradients, this function is missing, the declaration should be
    /// removed.
    void compute_electric_field_and_gradients();
    /// Compute electrostatic potentials at the nuclei
    void compute_esp_at_nuclei();
    /// Compute electrostatic potential at specified grid points
    void compute_esp_over_grid();
    /// Compute field at specified grid points
    void compute_field_over_grid();

    MultipolePropCalc mpc_;
    PopulationAnalysisCalc pac_;
    ESPPropCalc epc_;

    int max_noon_ = 3;

    // retrieves the Origin vector from the environment.
    Vector3 get_origin_from_environment() const;
    /// Computes the center for a given property, for the current molecule. Weighted center of geometry function
    Vector3 compute_center(const double* property) const;

   public:
    /// Constructor, uses globals
    OEProp(std::shared_ptr<Wavefunction> wfn);
    /// Destructor
    ~OEProp();

    /// Add a single task to the queue
    void add(const std::string& task);
    /// Add a set of tasks to the queue
    void add(std::vector<std::string> tasks);
    /// Clear task queue
    void clear();
    /// Set title for use in printout. Set the same title to be the name for printout.
    /// We do both because in older OEProp, the same member variable was used for both tasks.
    void set_title(const std::string& title) { names_ = {title + (title.empty() ? "" : " ") + "{}"}; title_ = title; }
    /// Set titles for use in saving information
    void set_names(const std::unordered_set<std::string> names) { names_ = names; }
    /// Compute and print/save the properties
    void compute();

    std::vector<double> const& Vvals() const { return epc_.Vvals(); }
    std::vector<double> const& Exvals() const { return epc_.Exvals(); }
    std::vector<double> const& Eyvals() const { return epc_.Eyvals(); }
    std::vector<double> const& Ezvals() const { return epc_.Ezvals(); }

    // These functions need to be overridden to pass on to the feature classes:

    // Change restricted flag. Resets C/D/epsilon matrices from wfn
    void set_wavefunction(std::shared_ptr<Wavefunction> wfn);
    // Change restricted flag. Resets C/D/epsilon matrices from wfn
    void set_restricted(bool restricted);
    // Set alpha eigenvalues, MO pitzer order basis
    void set_epsilon_a(SharedVector epsilon_a);
    // Set beta eigenvalues, MO pitzer order basis. Throws if restricted
    void set_epsilon_b(SharedVector epsilon_a);
    // Set alpha C matrix, SO/MO pitzer order basis.
    void set_Ca(SharedMatrix Ca);
    // Set beta C matrix, SO/MO pitzer order basis. Throws if restricted
    void set_Cb(SharedMatrix Cb);

    // => Set OPDM/TDM/DDM (often called). These need not be totally symmetric. Note, you are setting Da and/or Db, I do
    // the adding to Dt  <= //

    // TODO Add symmetry is irrep number
    void set_Da_ao(SharedMatrix Da, int symmetry = 0);
    void set_Db_ao(SharedMatrix Db, int symmetry = 0);
    void set_Da_so(SharedMatrix Da);
    void set_Db_so(SharedMatrix Db);
    void set_Da_mo(SharedMatrix Da);
    void set_Db_mo(SharedMatrix Db);
};

/**
 * The GridProp object, contains a cartesian grid and
 * associated point properties
 *
 * The grid is built according to the following rules:
 *  - A unit grid with corners (+/-1,0,0) (and all permutations) is built, with its origin at (0,0,0)
 *  - This unit grid is filled with (n_x, n_y, n_z) subintervals, equally spaced in each dimension
 *  - This unit grid is scaled symmetrically so that its edges measure (l_x, l_y, l_z)
 *  - The grid is translated so that the origin is as (o_x,o_y,o_z)
 **/
// class PSI_API GridProp : public Prop {
//
// protected:
//    /// The absolute file path where results from this analysis will be stored
//    std::string filename_;
//    /// The format for the output (defaults to df3)
//    std::string format_;
//
//    /// Is the grid initialized
//    bool initialized_;
//
//    /// The grid (grid_["x"] = <double***> x, for instance)
//    std::map<std::string, double***> grid_;
//    /// The BasisPoints object to evaluate basis functions on the grid
//    //std::shared_ptr<BasisPoints> points_;
//
//
//    /// The number of subintervals
//    int n_[3];
//    /// The dimensions of the grid (bohr)
//    double l_[3];
//    /// The origin of the final grid (bohr)
//    double o_[3];
//    /// The number of points to compute at once (5000 seems reasonable)
//    int block_size_;
//    /// The max/min for d3f writes
//    double caxis_[2];
//
//    /// An nbf x nbf scratch array
//    double** temp_tens_;
//
//    /// AO basis matrices (everything on grids is AO)
//    SharedMatrix Da_ao_;
//    SharedMatrix Db_ao_;
//    SharedMatrix Ca_ao_;
//    SharedMatrix Cb_ao_;
//    /// irrep offsets (for orbitals)
//    int irrep_offsets_[8];
//
//    /// The array of alpha MOs to store
//    std::vector<std::pair<int, int> > alpha_mos_;
//    /// The array of beta MOs to store
//    std::vector<std::pair<int, int> > beta_mos_;
//    /// The array of basis functions to plot (use C1 if you want AO, otherwise SO will be done)
//    std::vector<std::pair<int, int> > basis_funs_;
//
//    /// Common initialization
//    void common_init();
//
//    /// Print header
//    void print_header();
//
//    // Deprecated
//   // // Compute routines (these all work on a block of points)
//   // /// Compute mo values
//   // void compute_mos(std::shared_ptr<GridBlock> g, size_t offset);
//   // /// Compute basis function values
//   // void compute_basis_funs(std::shared_ptr<GridBlock> g, size_t offset);
//   // /// Compute total density
//   // void compute_rho(std::shared_ptr<GridBlock> g, double* results);
//   // /// Compute spin density (rho_a - rho_b)
//   // void compute_rho_s(std::shared_ptr<GridBlock> g, double* results);
//   // /// Compute rho_a (alpha density)
//   // void compute_rho_a(std::shared_ptr<GridBlock> g, double* results);
//   // /// Compute rho_b (beta density)
//   // void compute_rho_b(std::shared_ptr<GridBlock> g, double* results);
//   // /// Compute gamma_aa (\nabla rho_a ^2)
//   // void compute_gamma_aa(std::shared_ptr<GridBlock> g, double* results);
//   // /// Compute gamma_ab (\nabla rho_a \nabla rho_b)
//   // void compute_gamma_ab(std::shared_ptr<GridBlock> g, double* results);
//   // /// Compute gamma_bb (\nabla rho_b ^2)
//   // void compute_gamma_bb(std::shared_ptr<GridBlock> g, double* results);
//   // /// Compute tau_a (KE density)
//   // void compute_tau_a(std::shared_ptr<GridBlock> g, double* results);
//   // /// Compute tau_b (KE density)
//   // void compute_tau_b(std::shared_ptr<GridBlock> g, double* results);
//
//    /// Compute ESP (perhaps more involved, might need a fast Poisson solver)
//    void compute_ESP();
//
//    /// Allocate a grid3 (and zero it out)
//    double*** block_grid(int nx, int ny, int nz);
//    /// Free a grid3
//    void free_grid(double*** grid);
//    /// Actually build the grid
//    void build_grid();
//    /// allocate all the registers
//    void allocate_arrays();
//    /// Write the grid out in data format
//    void write_data_grid();
//    /// Write the grid out in df3 files
//    void write_df3_grid();
//
//
// public:
//    /// Constructor, uses globals
//    GridProp(std::shared_ptr<Wavefunction> wfn);
//    /// Constructor, uses globals and Process::environment::reference wavefunction
//    GridProp();
//    /// Destructor
//    virtual ~GridProp();
//
//    /// Python issue
//    void gridpy_add(const std::string& task) { add(task); }
//    void gridpy_compute() { compute(); }
//
//    /// Set the output filename
//    void set_filename(const std::string& file) { filename_ = file; }
//    /// Set the format
//    void set_format(const std::string& form) { format_ = form; }
//    /// Set a desired MO (use for restricted)
//    void add_alpha_mo(int irrep, int index);
//    /// Set a desired MO
//    void add_beta_mo(int irrep, int index);
//    /// Set a desired basis function
//    void add_basis_fun(int irrep, int index);
//    /// Compute and print/save the properties
//    void compute();
//
//    /**
//    * Grid specification
//    */
//
//    // High-Level
//
//    /// Set the number of subintervals on the grid
//    /// Setting 0 for one index will graph a plane (origin centered)
//    /// Setting 0 for two indices will graph a line (origin centered)
//    /// Defaults to nx,ny,nz = 40
//    void set_n(int nx, int ny, int nz) { n_[0] = nx; n_[1] = ny; n_[2] = nz; }
//
//    /// build the best grid for a given set of overages (bohr) at the edges
//    void build_grid_overages(double overages);
//
//    /// The maximum/minimum data values to plot for the current properties (clamp)
//    void set_caxis(double min, double max) { caxis_[0] = min; caxis_[1] = max; }
//
//    // Low-Level
//
//    /// Set the dimensions of the grid (in bohr)
//    void set_l(double lx, double ly, double lz) { l_[0] = lx; l_[1] = ly; l_[2] = lz; }
//
//    /// Set the origin of the grid (in bohr)
//    void set_o(double ox, double oy, double oz) { o_[0] = ox; o_[1] = oy; o_[2] = oz; }
//
//    /// Get the k-th element of n_
//    int get_n(int k) const { return n_[k]; }
//
//    /// Get the k-th element of l_
//    double get_l(int k) const { return l_[k]; }
//
//    /// Get the k-th element of o_
//    double get_o(int k) const { return o_[k]; }
//
//    /// Free the grid (useful if doing properties sequentially)
//    void reset();
//
//};

}  // namespace psi

#endif
