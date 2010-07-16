/*
 *  dfmp2.h
 *  matrix
 *
 *
 */

#ifndef DFMP2_H
#define DFMP2_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <libdiis/diismanager.h>
#include <psi4-dec.h>
#include <cstring>

using namespace psi;

namespace psi { namespace  dfmp2 {

class DFMP2 : public Wavefunction {
protected:

public:
protected:
    unsigned long int df_memory_;
    
    int print_;
    std::string algorithm_type_;

    shared_ptr<BasisSet> ribasis_;
    shared_ptr<BasisSet> zerobasis_;
    
    int nirreps_;
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
    
    //Build Qia directly on core
    void form_Qia_core();
    //Evaluate and sum the energy contributions
    void evaluate_contributions_core();   
    //Free Qia_ 
    void free_Qia_core();
 
    //Figure out all indexing and parameters
    //Only chkpt reads in here for disk-based
    //Called by constructor
    void setup();
    //Print the header
    void print_header();

public:
    DFMP2(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    /// Compute DF-MP2 energy.
    double compute_E();

    virtual ~DFMP2();
};

}}

#endif
