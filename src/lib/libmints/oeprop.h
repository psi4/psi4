#ifndef _psi_src_lib_oeprop_h
#define _psi_src_lib_oeprop_h

#include <set>
#include <vector>
#include <map>
#include <string>

#include "typedefs.h"

namespace boost {
template<class T> class shared_ptr;
// Forward declarations for boost.python used in the extract_subsets
}

namespace psi {

class Matrix;
class Wavefunction;
class IntegralFactory;
class MatrixFactory;
class BasisSet;

/**
* The Prop object, base class of OEProp and GridProp objects
*
* Word on the street:
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
class Prop {

protected:
    /// Print flag
    int print_;
    /// Debug flag
    int debug_;

    /// The title of this Prop object, for use in saving info
    std::string title_;

    /// The set of tasks to complete
    std::set<std::string> tasks_;

    /// The wavefunction object this Prop is built around
    boost::shared_ptr<Wavefunction> wfn_;
    /// The basisset for this wavefunction
    boost::shared_ptr<BasisSet> basisset_;
    /// Is this wavefunction object spin-restricted? (Actually closed-shell, but this is wavefunction's convention)
    bool restricted_; // This allows pointers to be duplicated/computation skipped
    /// The integral factory for this wavefunction's basisset
    boost::shared_ptr<IntegralFactory> integral_;
    /// The matrix factory for this wavefunction's basisset (SO)
    boost::shared_ptr<MatrixFactory> factory_;

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

    /// Print header
    virtual void print_header() = 0;

public:

    /// Build a Prop object with C, epsilon, and restricted buit from wfn
    Prop(boost::shared_ptr<Wavefunction> wfn);
    /// Virtual destructor
    virtual ~Prop();

    // => Wavefunction Modifiers (rarely called, C is usually fixed at HF) <= //

    // Change restricted flag. Resets C/D/epsilon matrices from wfn
    void set_wavefunction(boost::shared_ptr<Wavefunction> wfn);
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

    // => Set OPDM/TDM/DDM (often called). These need not be totally symmetric. Note, you are setting Da and/or Db, I do the adding to Dt  <= //

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
    SharedMatrix Dt_ao(bool total = true);
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

    // => Queue/Compute Routines <= //

    /// Add a single task to the queue
    void add(const std::string& task);
    /// Add a set of tasks to the queue
    void add(std::vector<std::string> tasks);
    /// Clear task queue
    void clear();

    /// Set title for use in saving information
    void set_title(const std::string& title) { title_ = title; }

    /// Compute properties
    virtual void compute() = 0;

    // => Utility Routines <= //

    void set_print(int print) { print_; }
    void set_debug(int debug) { debug_; }

};


/**
* The OEProp object, computes arbitrary expectation values (scalars)
* analyses (typically vectors)
**/
class OEProp : public Prop {

protected:
    /// Common initialization
    void common_init();
    /// Print header and information
    void print_header();

    // Compute routines
    /// Compute dipole
    void compute_dipole(bool transition = false);
    /// Compute quadrupole
    void compute_quadrupole(bool transition = false);
    /// Compute mo extents
    void compute_mo_extents();
    /// Compute Mulliken Charges
    void compute_mulliken_charges();
    /// Compute Lowdin Charges
    void compute_lowdin_charges();
    /// Compute Mayer Bond Indices (non-orthogoal basis)
    void compute_mayer_indices();
    /// Compute Wiberg Bond Indices using Lowdin Orbitals (symmetrically orthogonal basis)
    void compute_wiberg_lowdin_indices();
    /// Compute/display natural orbital occupations around the bandgap. Displays max_num above and below the bandgap
    void compute_no_occupations(int max_num = 3);
    /// Compute electric field and electric field gradients
    void compute_electric_field_and_gradients();

public:
    /// Constructor, uses globals
    OEProp(boost::shared_ptr<Wavefunction> wfn);
    /// Constructor, uses globals and Process::environment::reference wavefunction
    OEProp();
    /// Destructor
    virtual ~OEProp();

    /// Python issue
    void oepy_add(const std::string& task) { add(task); }
    void oepy_compute() { compute(); }
    void oepy_set_title(const std::string& title) { set_title(title); }

    /// Compute and print/save the properties
    void compute();
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
class GridProp : public Prop {

protected:
    /// The absolute file path where results from this analysis will be stored
    std::string filename_;
    /// The format for the output (defaults to df3)
    std::string format_;

    /// Is the grid initialized
    bool initialized_;

    /// The grid (grid_["x"] = <double***> x, for instance)
    std::map<std::string, double***> grid_;
    /// The BasisPoints object to evaluate basis functions on the grid
    //boost::shared_ptr<BasisPoints> points_;

    /// The number of subintervals
    int n_[3];
    /// The dimensions of the grid (bohr)
    double l_[3];
    /// The origin of the final grid (bohr)
    double o_[3];
    /// The number of points to compute at once (5000 seems reasonable)
    int block_size_;
    /// The max/min for d3f writes
    double caxis_[2];

    /// An nbf x nbf scratch array
    double** temp_tens_;

    /// AO basis matrices (everything on grids is AO)
    SharedMatrix Da_ao_;
    SharedMatrix Db_ao_;
    SharedMatrix Ca_ao_;
    SharedMatrix Cb_ao_;
    /// irrep offsets (for orbitals)
    int irrep_offsets_[8];

    /// The array of alpha MOs to store
    std::vector<std::pair<int, int> > alpha_mos_;
    /// The array of beta MOs to store
    std::vector<std::pair<int, int> > beta_mos_;
    /// The array of basis functions to plot (use C1 if you want AO, otherwise SO will be done)
    std::vector<std::pair<int, int> > basis_funs_;

    /// Common initialization
    void common_init();

    /// Print header
    void print_header();

    // Deprecated
   // // Compute routines (these all work on a block of points)
   // /// Compute mo values
   // void compute_mos(boost::shared_ptr<GridBlock> g, unsigned long int offset);
   // /// Compute basis function values
   // void compute_basis_funs(boost::shared_ptr<GridBlock> g, unsigned long int offset);
   // /// Compute total density
   // void compute_rho(boost::shared_ptr<GridBlock> g, double* results);
   // /// Compute spin density (rho_a - rho_b)
   // void compute_rho_s(boost::shared_ptr<GridBlock> g, double* results);
   // /// Compute rho_a (alpha density)
   // void compute_rho_a(boost::shared_ptr<GridBlock> g, double* results);
   // /// Compute rho_b (beta density)
   // void compute_rho_b(boost::shared_ptr<GridBlock> g, double* results);
   // /// Compute gamma_aa (\nabla rho_a ^2)
   // void compute_gamma_aa(boost::shared_ptr<GridBlock> g, double* results);
   // /// Compute gamma_ab (\nabla rho_a \nabla rho_b)
   // void compute_gamma_ab(boost::shared_ptr<GridBlock> g, double* results);
   // /// Compute gamma_bb (\nabla rho_b ^2)
   // void compute_gamma_bb(boost::shared_ptr<GridBlock> g, double* results);
   // /// Compute tau_a (KE density)
   // void compute_tau_a(boost::shared_ptr<GridBlock> g, double* results);
   // /// Compute tau_b (KE density)
   // void compute_tau_b(boost::shared_ptr<GridBlock> g, double* results);

    /// Compute ESP (perhaps more involved, might need a fast Poisson solver)
    void compute_ESP();

    /// Allocate a grid3 (and zero it out)
    double*** block_grid(int nx, int ny, int nz);
    /// Free a grid3
    void free_grid(double*** grid);
    /// Actually build the grid
    void build_grid();
    /// allocate all the registers
    void allocate_arrays();
    /// Write the grid out in data format
    void write_data_grid();
    /// Write the grid out in df3 files
    void write_df3_grid();


public:
    /// Constructor, uses globals
    GridProp(boost::shared_ptr<Wavefunction> wfn);
    /// Constructor, uses globals and Process::environment::reference wavefunction
    GridProp();
    /// Destructor
    virtual ~GridProp();

    /// Python issue
    void gridpy_add(const std::string& task) { add(task); }
    void gridpy_compute() { compute(); }

    /// Set the output filename
    void set_filename(const std::string& file) { filename_ = file; }
    /// Set the format
    void set_format(const std::string& form) { format_ = form; }
    /// Set a desired MO (use for restricted)
    void add_alpha_mo(int irrep, int index);
    /// Set a desired MO
    void add_beta_mo(int irrep, int index);
    /// Set a desired basis function
    void add_basis_fun(int irrep, int index);
    /// Compute and print/save the properties
    void compute();

    /**
    * Grid specification
    */

    // High-Level

    /// Set the number of subintervals on the grid
    /// Setting 0 for one index will graph a plane (origin centered)
    /// Setting 0 for two indices will graph a line (origin centered)
    /// Defaults to nx,ny,nz = 40
    void set_n(int nx, int ny, int nz) { n_[0] = nx; n_[1] = ny; n_[2] = nz; }

    /// build the best grid for a given set of overages (bohr) at the edges
    void build_grid_overages(double overages);

    /// The maximum/minimum data values to plot for the current properties (clamp)
    void set_caxis(double min, double max) { caxis_[0] = min; caxis_[1] = max; }

    // Low-Level

    /// Set the dimensions of the grid (in bohr)
    void set_l(double lx, double ly, double lz) { l_[0] = lx; l_[1] = ly; l_[2] = lz; }

    /// Set the origin of the grid (in bohr)
    void set_o(double ox, double oy, double oz) { o_[0] = ox; o_[1] = oy; o_[2] = oz; }

    /// Get the k-th element of n_
    int get_n(int k) const { return n_[k]; }

    /// Get the k-th element of l_
    double get_l(int k) const { return l_[k]; }

    /// Get the k-th element of o_
    double get_o(int k) const { return o_[k]; }

    /// Free the grid (useful if doing properties sequentially)
    void reset();

};


}

#endif
