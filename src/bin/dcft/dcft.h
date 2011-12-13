#ifndef _PSI4_SRC_BIN_DCFT_DCFT_H_
#define _PSI4_SRC_BIN_DCFT_DCFT_H_

#include <libmints/matrix.h>
#include <libmints/vector.h>
#include <libmints/wavefunction.h>
#include <libdpd/dpd.h>
#include <libciomr/libciomr.h>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi{

class Options;
class PSIO;
class Chkpt;
class Matrix;
class Vector;
class IntegralTransform;

namespace dcft{

class DCFTSolver:public Wavefunction
{
public:
    DCFTSolver(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options);
    ~DCFTSolver();

    double compute_energy();
protected:
    IntegralTransform *_ints;

    void finalize();
    void transform_integrals();
    void build_lambda();
    void init();
    void compute_dcft_energy();
    void update_lambda_from_residual();
    void compute_scf_energy();
    void mp2_guess();
    void build_tau();
    void print_opdm();
    void write_orbitals_to_checkpoint();
    void check_n_representability();
    void print_orbital_energies();
    void find_occupation(SharedVector &evals_a, SharedVector &evals_b, bool forcePrint = false);
    void build_intermediates();
    void process_so_ints();
    void build_G();
    void build_tensors();
    void build_denominators();
    void dump_density();
    void mulliken_charges();
    void dpd_buf4_add(dpdbuf4 *A, dpdbuf4 *B, double alpha);
    void scf_guess();
    void half_transform(dpdbuf4 *A, dpdbuf4 *B, SharedMatrix& C1, SharedMatrix& C2,
                        int *mospi_left, int *mospi_right, int **so_row, int **mo_row,
                        bool backwards, double alpha, double beta);
    void file2_transform(dpdfile2 *A, dpdfile2 *B, SharedMatrix C, bool backwards);
    void AO_contribute(dpdbuf4 *tau1_AO, dpdbuf4 *tau2_AO, int p, int q,
                       int r, int s, double value, dpdfile2* = NULL, dpdfile2* = NULL, dpdfile2* = NULL);
    void compute_tau_squared();
    void compute_energy_tau_squared();
    //void AO_contribute(dpdfile2 *tau1_AO, dpdfile2 *tau2_AO, int p, int q,
    //        int r, int s, double value);
    bool correct_mo_phases(bool dieOnError = true);
    double compute_lambda_residual();
    double compute_scf_error_vector();
    double update_scf_density(bool damp = false);
    /// Whether to force the code to keep the same occupation from SCF
    bool lock_occupation_;
    /// The maximum number of lambda iterations per update
    int lambdamaxiter_;
    /// The maximum number of SCF iterations per update
    int scfmaxiter_;
    /// The amount of information to print
    int print_;
    /// The number of unique pairs of symmetrized atomic orbitals
    int ntriso_;
    /// The number of active alpha electrons
    int nalpha_;
    /// The number of active beta electrons
    int nbeta_;
    /// The number of virtual alpha orbitals
    int navir_;
    /// The number of virtual beta orbitals
    int nbvir_;
    /// The maximum size of the DIIS subspace
    int maxdiis_;
    /// The number of DIIS vectors needed for extrapolation
    int mindiisvecs_;
    /// The maximum number of iterations
    int maxiter_;
    /// The number of occupied alpha orbitals per irrep
    int *naoccpi_;
    /// The number of occupied beta orbitals per irrep
    int *nboccpi_;
    /// The number of virtual alpha orbitals per irrep
    int *navirpi_;
    /// The number of virtual beta orbitals per irrep
    int *nbvirpi_;
    /// The nuclear repulsion energy in Hartree
    double enuc_;
    /// The cutoff below which and integral is assumed to be zero
    double int_tolerance_;
    /// The RMS value of the error vector after the SCF iterations
    double scf_convergence_;
    /// The RMS value of the change in lambda after the lambda iterations
    double lambda_convergence_;
    /// The convergence criterion for the lambda iterations
    double lambda_threshold_;
    /// The convergence criterion for the scf iterations
    double scf_threshold_;
    /// The convergence that must be achieved before DIIS extrapolation starts
    double diis_start_thresh_;
    /// The SCF component of the energy
    double scf_energy_;
    /// The Lambda component of the energy
    double lambda_energy_;
    /// The Tau^2 correction to the SCF component of the energy
    double energy_tau_squared_;
    /// The previous total energy
    double old_total_energy_;
    /// The updated total energy
    double new_total_energy_;
    /// The Tikhonow regularizer used to remove singularities (c.f. Taube and Bartlett, JCP, 2009)
    double regularizer_;
    /// The alpha occupied eigenvectors, per irrep
    SharedMatrix aocc_c_;
    /// The beta occupied eigenvectors, per irrep
    SharedMatrix bocc_c_;
    /// The alpha virtual eigenvectors, per irrep
    SharedMatrix avir_c_;
    /// The beta virtual eigenvectors, per irrep
    SharedMatrix bvir_c_;
    /// The Tau matrix in the AO basis, stored by irrep, to perturb the alpha Fock matrix
    SharedMatrix a_tau_;
    /// The Tau matrix in the AO basis, stored by irrep, to perturb the beta Fock matrix
    SharedMatrix b_tau_;
    /// The Tau^2 correction to the alpha Tau matrix in the AO basis
    SharedMatrix a_tautau_;
    /// The Tau^2 correction to the beta Tau matrix in the AO basis
    SharedMatrix b_tautau_;
    /// The overlap matrix in the AO basis
    SharedMatrix ao_s_;
    /// The one-electron integrals in the SO basis
    SharedMatrix so_h_;
    /// The alpha Fock matrix in the SO basis
    SharedMatrix Fa_;
    /// The beta Fock matrix in the SO basis
    SharedMatrix Fb_;
    /// The copy of the alpha Fock matrix in the SO basis before the Lowdin orthogonalization
    SharedMatrix Fa_copy;
    /// The copy of the alpha Fock matrix in the SO basis before the Lowdin orthogonalization
    SharedMatrix Fb_copy;
    /// The alpha Fock matrix in the MO basis
    SharedMatrix moFa_;
    /// The beta Fock matrix in the MO basis
    SharedMatrix moFb_;
    /// The inverse square root overlap matrix in the SO basis
    SharedMatrix s_half_inv_;
    /// The old full alpha MO coefficients
    SharedMatrix old_ca_;
    /// The old full beta MO coefficients
    SharedMatrix old_cb_;
    /// The alpha kappa matrix in the SO basis
    SharedMatrix kappa_a_;
    /// The beta kappa matrix in the SO basis
    SharedMatrix kappa_b_;
    /// The alpha external potential in the SO basis
    SharedMatrix g_tau_a_;
    /// The beta external potential in the SO basis
    SharedMatrix g_tau_b_;
    /// The alpha SCF error vector
    SharedMatrix scf_error_a_;
    /// The beta SCF error vector
    SharedMatrix scf_error_b_;
    /// Used to align things in the output
    std::string indent;
};

}} // Namespaces

#endif // Header guard
