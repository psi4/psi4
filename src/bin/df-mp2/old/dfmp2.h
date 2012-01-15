/*
 *  dfmp2.h
 *  matrix
 *
 *
 */

//#define _MKL
#ifndef DFMP2_H
#define DFMP2_H

#include <float.h>
#include <libmints/wavefunction.h>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class Options;
class PSIO;
class Chkpt;

namespace  dfmp2 {

class DFMP2 : public Wavefunction {
protected:

public:
    virtual bool same_a_b_orbs() const { return reference_wavefunction_->same_a_b_orbs(); }
    virtual bool same_a_b_dens() const { return reference_wavefunction_->same_a_b_dens(); }
protected:

    int nproc_;
    int rank_;

    unsigned long int df_memory_;

    int print_;
    std::string algorithm_type_;
    std::string fitting_symmetry_;
    std::string fitting_conditioning_;
    std::string fitting_inversion_;

    boost::shared_ptr<BasisSet> ribasis_;
    boost::shared_ptr<BasisSet> zerobasis_;

    int nirrep_;
    int ndocc_;
    int nvirt_;
    int nf_docc_;
    int nf_virt_;
    int nact_docc_;
    int nact_virt_;
    int nso_;
    int nmo_;
    int naux_raw_;
    int naux_fin_;

    int* clsdpi_;
    int* orbspi_;
    int* frzcpi_;
    int* frzvpi_;

    double ss_scale_;
    double os_scale_;

    double E_scf_;
    long double E_ss_;
    long double E_os_;
    double E_;
    double E_tot_;
    double E_sss_;
    double E_sos_;
    double E_scs_;
    double E_tot_scs_;

    double** C_docc_;
    double** C_virt_;
    double* eps_docc_;
    double* eps_virt_;
    int* sym_docc_;
    int* sym_virt_;

    double** Cia_;
    double** Aia_;
    double** Qia_;
    double** W_; //Fitting metric inverse sqrt or cholesky decomposition

    double schwarz_;
    int* schwarz_shell_pairs_;
    unsigned long int sig_shell_pairs_;

    //Slow MPI version
    double compute_E_old();
    //Threaded core algorithm
    double compute_E_core();
    //Threaded disk algorithm
    double compute_E_disk();

    //Form the Schwarz Sieve
    void form_Schwarz();

    // Raw fitting metric
    double** form_W(boost::shared_ptr<BasisSet> b);
    // Raw fitting overlap
    double** form_W_overlap(boost::shared_ptr<BasisSet> b);
    // Form orthogonlization matrix
    double** form_X(double** W, int n, double max_cond, int* clipped);
    // Matrix power, in-place, with clipping
    int matrix_power(double** J, int n, double pow, double max_cond = DBL_MAX);
    // Matrix upper cholesky decomposition (U'U = J, FORTRAN -style), in-place
    void cholesky_decomposition(double** J, int n);
    // Matrix inverse, via cholesky decomposition, in-place
    void cholesky_inverse(double** J, int n);

    //Form the inverse square root of the fitting metric (preconditioned)
    void form_Wm12_fin();
    //Form the inverse square root of the fitting metric (raw)
    void form_Wm12_raw();
    //Form the square root of the fitting metric and it's cholesky decomposition
    void form_Wp12_chol();

    //Build the (A|ia) tensor and stripe it here
    void form_Aia_disk();
    //Embed the fitting  and stripe here
    void form_Qia_disk();
    //Manage the IO and delegate the summation
    void evaluate_contributions_disk();
    //Do the summation
    void find_disk_contributions(double** Qia, double** Qjb, double***I, int* starts, int* sizes, int* stops, int block1, int block2);

    //Build Aia directly on corei (returns double** as transform might be in place)
    double** form_Aia_core();
    //Build Qia directly on core
    void form_Qia_core();
    //Evaluate and sum the energy contributions
    void evaluate_contributions_core_sym();
    //Build Aia directly on corei (returns double** as transform might be in place)
    double** form_Aia_core_parallel();
    //Build Qia directly on core
    void form_Qia_core_parallel();
    //Evaluate and sum the energy contributions NEW PARALLEL ALGORITHM
    void evaluate_contributions_core_parallel();
    //Free Qia_
    void free_Qia_core();

    //Build Cia directly on core
    void form_Cia_core();
    //Evaluate and sum the energy contributions
    void evaluate_contributions_core_asym();
    //Free Cia_
    void free_Cia_core();

    //Figure out all indexing and parameters
    //Only chkpt reads in here for disk-based
    //Called by constructor
    void setup();
    //Print the header
    void print_header();

public:
    DFMP2(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    /// Compute DF-MP2 energy.
    double compute_energy();

    virtual ~DFMP2();
};

}}

#endif
