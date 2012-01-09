/*! \file
    \ingroup LMP2
    \brief Enter brief description of file here
*/
#ifndef _psi_src_lib_lmp2_h_
#define _psi_src_lib_lmp2_h_

#include <string>
#include <sstream>

#include <libpsio/psio.hpp>
#include <libmints/mints.h>
#include <libdiis/diismanager.h>
#include <libparallel/parallel.h>
#include <libqt/qt.h>
#include "liblmp2_solver/dist_container.h"
#include "psi4-dec.h"

#ifdef HAVE_MADNESS
    #include <world/worldobj.h>
#endif



namespace boost {
template<class T> class shared_ptr;
}

namespace psi{

class BasisSet;
class Options;
class PSIO;
class Chkpt;

namespace lmp2 {

#ifdef HAVE_MADNESS

#ifdef HAVE_MADNESS
class LMP2 : public Wavefunction, public madness::WorldObject<LMP2> {
#else
class LMP2 : public Wavefunction

#endif

protected:

//    Options options_;
//    boost::shared_ptr<PSIO> psio_;

    // ****  These are for parallel computing ***
    /// The process ID
    int me_;
    /// The total number of processes
    int nproc_;
    /// The number of threads
    int nthread_;
    /// The communicator type
    std::string comm_;
    /// The madness world communicator
#ifdef HAVE_MADNESS
     SharedMadWorld madworld_;
#endif

    // ****  These are global shared pointers ****
    boost::shared_ptr<Wavefunction> wfn_;
    boost::shared_ptr<IntegralFactory> integral_;
    std::vector< boost::shared_ptr<TwoElectronInt> > eri_;

    // ****  These are parameters taken from the user input options ***
    /// The reference that will be used
    std::string reference_;
    /// Wether we are using ri
    bool ri_lmp2_;
    /// The print level
    int print_;
    /// The maximum number of iterations
    int maxiter_;
    /// The energy convergence criteria
    double econv_;
    /// The screening threshold for the integrals (Not sure if this is working)
    double escreen_;
    /// Wether we are screening the integrals or not (Not sure if this is working)
    bool screen_int_;
    /// RMS Convergence criteria
    double rmsconv_;
    /// This will skip parts of the fock matrix if value is below this value
    double fskip_;
    /// Wether we are using diis extrapolation or not
    bool diis_;
    /// The iteration number to start diis extrapolation. Needs to be greater than 2.
    int diis_start_;
    /// The number of matrices to utilize in the diis extrapolation
    int max_diis_vectors_;
    /// The localization cutoff
    double cutoff_;
    /// Whether distant pairs will be neglected
    bool neglect_dp_;
    /// Distant pair cutoff criteria
    double dp_cutoff_;
    /// Only print the domains and then exit
    bool only_print_domains_;

    // ****  These are orbital and molecule information  ****
    /// The number of irreps
    int nirreps_;
    /// The number of doubly occupied orbitals
    int ndocc_;
    /// The number of virtual orbitals
    int nvirt_;
    /// The number of frozed
    int nfocc_;
    /// The number of frozen virtuals
    int nfvir_;
    /// The number of atomic orbitals
    int nso_;
    /// The number of active doubly occupied orbitals
    int nact_docc_;
    /// The number of active virtual orbitals
    int nact_virt_;
    /// The number of atoms in the molecule
    int natom_;
    /// The number of AM shells in the basis set
    int nshell_;

    /// The reference molecular orbital coefficients in the AO basis
    SharedMatrix C_;
    /// The density matrix in the AO basis
    SharedMatrix D_AO_;
    /// The overlap matrix
    SharedMatrix S_;
    /// The Fock matrix in the AO basis
    SharedMatrix F_AO_;
    /// The Fock matrix in the localized basis
    SharedMatrix F_LO_;
    /// The complete virtual space projector
    SharedMatrix Rt_full_;
    /// The scf energy
    double escf_;
    /// The nuc. energy
    double enuc_;

    /// Where the AO starts for atom A
    std::vector<int> ao_start_;
    /// Where the AO stops for atom A
    std::vector<int> ao_stop_;

    //  **** These are for determing the domains ****
    /// Whether the pair domain exists (i.e. if the domains i and j are distant)
    std::vector<int> pairdom_exist_;
    /// Orbital domains
    std::vector< std::vector<int> > domain_;
    /// Size of each domain
    std::vector<int> domain_len_;
    /// Orbital pair domains
    std::vector< std::vector<int> > pair_domain_;
    /// Size of the pair domains
    std::vector<int> pair_domain_len_;
    /// Length of the non-redundant pair domains
    std::vector<int> pair_domain_nr_len_;

    // **** These are some maps for the ij pairs ****
    /// A map of which process owns pair ij
    std::vector<int> ij_owner_;
    /// A map of ij to the local value of ij on a specific process
    std::map<int,int> ij_local_;
    /// map of ij to i after distant pair has been removed
    std::map<int,int> ij_i_map_;
    /// map of ij to j after distant pair has been removed
    std::map<int,int> ij_j_map_;
    /// map of ij to ij neglecting distant pairs
    std::map<int,int> ij_map_neglect_;
    /// map of ij to kj
    std::map<int, std::map<int,int> > ij_kj_map_;
//    std::vector< std::vector<int> > ij_kj_map_;
    /// map of ij to ik
    std::map<int, std::map<int,int> > ij_ik_map_;
//    std::vector< std::vector<int> > ij_ik_map_;

    /// The total number of ij pairs
    int ij_pairs_;
    /// The number of ij pairs each process owns
    int ij_pairs_per_proc_;

    // **** These get global distributed Matrices **** //
    /// The orbital energies
    std::vector< SharedVector > evals_;
    /// The virtual overlap matrix
    std::vector< SharedMatrix > S_virt_;
    /// The virtual overlap matrix
    std::vector< SharedMatrix > F_virt_;
    /// The projection matrices
    std::vector< SharedMatrix > W_;

    // **** These are used to distribute the integrals over MN pairs *** //
    typedef std::pair< int, std::vector<int> > Shell_Pair;
    std::multimap< int, std::vector<int>, std::greater<int> > MN_Pairs_;
    std::multimap< int, std::vector<int>, std::greater<int> >::iterator MN_iter_;
    std::vector<int> MN_Owner_;
    std::map<int,int> MN_local_;
    int mn_pairs_per_proc_;
    int maxshell_;

    /// half transformed integrals, distributed over MN pairs
    std::map<int,SharedMatrix> eri_2_MN_;
//    std::vector< SharedMatrix > eri_2_MN_;
    /// half transformed integrals, distributed over IJ pairs
    std::map<int,SharedMatrix> eri_ij_;
//    std::vector< SharedMatrix > eri_ij_;

    /// Testing a sparse distributed container
    distributed_container Integrals_;

    /// The lmp2 T2 amplitudes
    std::vector<std::map<int, SharedMatrix> > T2_amp_;
//    std::vector< std::vector< SharedMatrix > > T2_amp_;
    /// The DIIS extrapolated T2s
    std::vector<std::map<int, SharedMatrix> > T2_ext_;
//    std::vector< std::vector< SharedMatrix > > T2_ext_;
    /// The DIIS error matrix
    std::vector<std::map<int, SharedMatrix> > error_;
//    std::vector< std::vector< SharedMatrix > > error_;

    /// The lmp2 correlation energy for the current iteration
    double Elmp2_;
    /// The lmp2 correlation energy from the previous iteration
    double Elmp2_old_;
    /// The energy difference between the current and previous iteration
    double Delta_Elmp2_;
    /// The RMS value between the current and previous iteration
    double Drms_T2_;

    // *** This is some info for diis **** //
    int div_;
    int matsize_;
    int dmat1_;
    int dmat2_;
    int nmat_;
    int omat_;
    std::vector< SharedMatrix > Bmat_;

    /// Print the results from the current iteration
    void print_results(const int &iter) const;

    /// Print a summary of the LMP2 computation
    void print_summary() const;


    /// Madness mutex
#ifdef HAVE_MADNESS
    SharedMutex print_mutex;
    SharedMutex mutex_;
    SharedMutex F_mutex_;
#endif

    boost::shared_ptr<MatrixFactory> nso_nso_;
    boost::shared_ptr<MatrixFactory> occ_occ_;

    /// This function gets the MO information
    void moinfo();
    /// This function prints the input parameter information
    void print_params() const;
    /// This function prints the MO/basis set information
    void print_moinfo() const;
    /// Setup the molecule, wfn, basis, and integralfactory
    void common_init();
    /// Prints the LMP2 header to the outfile
    void print_header() const;
    /// Gets the user supplied input parameters
    void params();
    /// Set up all of the required matrix factories
    void setup_factories();
    /// Get information from the wfn (i.e. MOC, Density, and Fock matricies)
    void reference();

    /// Computes the AO overlap matrix (uses libmints)
    void overlap();

    /// Localizes the AO MO coefficients (pipek-mezey)
    void localize_pipek_mezey();

    /// Transform the Fock matrix from the AO to LO basis
    void transform_fock();

    /// Build the orbital domains
    void build_domains();

    /// Print the domains
    void print_domains(const std::vector<double> &);

    /// Builds a map that maps ij to the i/j value after the distant pairs have been removed
    void set_ij_maps();

    /// Builds a map for ij to kj/ik
    void setup_ij_kj_ik_maps();

    /// Sets up ij_owner_ and ij_local_
    void set_ij_owner_local();

    /// Build the pairdomains and determine the length of each pairdomain
    void build_pairdomains();

    /// This initializes some global matrices (e.g. Rt_full_, W_, ...)
    void init_global_mat();

    /// Build the projection matrix W for pair ij
    int build_projection_matrix(const int &);

    /// Build the projection matrices
    void projection();

    /// Set up some maps for the first half integal transformation
    void set_mn_integral_maps();

    /// Perform an integral direct transformation
    int first_half_integral_transformation(const int &M, const int &N,
                                           const int &numm, const int &numn,
                                           const int &MN);

    /// Perform the direct integral transformation
    void integral_direct_transformation();

    /// setup the matrices for the two electron integrals
    void setup_eri_matrices();

    std::vector<double> copy_eri_mn(const int &MN, const int &ij,
                                    const int &numm, const int &numn);

    /// Perform the second half integral transformation
    int second_half_transformation(const int &ij);

    /// Print the two electron integrals
    int print_eri(const int &ij, std::string name);

    /// Set up some info for diis extrapolation
    void setup_diis(const int &iter);


#ifdef HAVE_MADNESS

    /// This copies part of the eri_mn matrix and returns it as a future
    madness::Future< std::vector<double> > copy_eri_mn_mad(const int &MN, const int &ij,
                                          const int &numm, const int &numn);

    /// This sends a task (including the integrals) of redist_integrals to the owner of ij
    madness::Void redist_madness_integrals(const int &MN, const int &ij,
                                           const int &numm, const int &numn,
                                           const int &abs_m, const int &abs_n);

    /// Compute the T2 amplitudes
    SharedMatrix amplitudes_T2(const SharedMatrix T2, const int &ij, const int &iter);

    /// Returns the ij'th T2 amplitudes may or may not be local
    SharedMatrix get_old_T2(const int &ij);

    /// Builds and returns F_sum to compute the new Residules
    SharedMatrix build_F_sum(const int &ij, const int &iter);

    /// Compute the T2 amplitudes
    madness::Void compute_T2(const SharedMatrix & F_sum, const int &ij, const int &iter);

    /// Build Rtilde
    SharedMatrix build_rtilde(const SharedMatrix F_sum,
                                               const int &ij, const int &iter);

    /// Adds the T2_kj and T2_ik to F_sum
    SharedMatrix F_sum_add(const std::vector<madness::Future<SharedMatrix> > &T2_kj,
                                 const std::vector<madness::Future<SharedMatrix> > &T2_ik,
                                 const int &ij);

        //    /// Adds the the T2_kj amplitudes to F_sum
//    madness::Void F_sum_add_kj(SharedMatrix F_sum, const SharedMatrix T2_kj,
//                                     const int &k, const int &ij);
//    /// Adds the the T2_ik amplitudes to F_sum
//    madness::Void F_sum_add_ik(SharedMatrix F_sum, const SharedMatrix T2_ik,
//                                     const int &k, const int &ij);

    madness::Future<int> build_error(const int & ij);

    madness::Void build_Bmat(const int &check, const int &ij);

    madness::Void build_T2_ext(const int &ij, const std::vector<double> &co);

    madness::Void Elmp2_local_sum(const double &val);
    madness::Void Drms_local_sum(const double &val);

#endif

    void perform_diis(const int &iter);

    /// Redistribute a eri_mn, which is a small portion of the full two electron integrals
    int redist_integrals_mad(const SharedMatrix eri_mn,
                         const int &ij, const int &numm, const int &numn,
                         const int &abs_m, const int &abs_n);
    int redist_integrals(const std::vector<double> &eri_mn,
                         const int &ij, const int &numm, const int &numn,
                         const int &abs_m, const int &abs_n);

    /// Redistribute all of the integrals from MN to ij pairs
    void redistribute_integrals();

    /// Allocate memory for T2's
    void allocate_T2_memory();

    /// Build T2 using madness
    int compute_T2_energy(const int &iter);

    /// Compute the energy contribution from the give T2's
    madness::Future<double> energy(const SharedMatrix T2, const int &ij, const int &iter);

    /// Check to see if the energy and T2 amplitudes are converged
    madness::Future<double> T2_rms(const SharedMatrix T2,
                                   const SharedMatrix T2_old,
                                   const int &ij);

    /// Save the current T2 amplitudes
    int store_T2(const SharedMatrix T2, const int &ij);

    /// Prints a shared matrix
    int print_matrix(const SharedMatrix mat, const std::string &str=NULL) const;

public:
    LMP2(Options& options, boost::shared_ptr<Wavefunction> ref_wfn);

    virtual ~LMP2();

    virtual double compute_energy();

    virtual bool same_a_b_orbs() const { return true };
    virtual bool same_a_b_dens() const { return true };

};

#endif // End of if HAVE_MADNESS

}}

#endif /* Header guard */

