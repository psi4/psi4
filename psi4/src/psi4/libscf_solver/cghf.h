#ifndef __math_test_cghf_h__
#define __math_test_cghf_h__

#include "psi4/libpsio/psio.hpp"
#include "psi4/libfock/v.h"
#include "hf.h"
#include <Einsums/Config.hpp>
#include "Einsums/Tensor.hpp"

namespace psi {
namespace scf {

class CGHF : public HF {
   public:
    CGHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional);
    CGHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional, Options& options,
         std::shared_ptr<PSIO> psio);
    ~CGHF() override;

    // Declares all BlockTensors, and some useful misc variables needed
    void common_init();

    // Constructs the spin-blocked overlap (EINS_), core Hamiltonian (F0_), and orthogonalization (EINX_) matrices
    void preiterations();

    // Allow SAP and SAD initial guesses
    void sap_guess();
    void compute_SAD_guess(bool natorb) override;
    
    // Empty functions for now -- seems to be resetting and saving matrices, respectively
    void finalize() override;
    void save_density_and_energy() override;

    // Form orbital gradient FDSmSDF_ = [F, D]
    void form_FDSmSDF();

    // Compute the norm from the orbital gradient as a second test of convergence
    double compute_Dnorm();

    // Empty function for now -- needed for DFT later
    void form_V() override;

    // Computes J and K either explicitly (4-index) or with RI (3-index)
    void form_G() override;

    // Forms Fock matrix F = F0_ + J - K
    void form_F() override;

    // Orthogonalizes then diagonalizes the Fock matrix to form the coefficient matrix C_
    void form_C(double shift) override;

    // Constructs 1-particle density matrix using the occupied coefficients Cocc_ and the conjugate (stored in temp1_, not permanently stored)
    void form_D() override;

    // Empty function for now, but for UHF and RHF, scales the density matrix
    void damping_update(double damping_percentage) override;

    // Compute the energy based purely off F0_ = T + V, with no J and K
    double compute_initial_E() override;

    // Compute 1e and 2e energy separately, then combine with nuclear repulsion energy nuclearrep_ to return a total energy 
    double compute_E() override;

    // Empty functions for now
    void setup_potential() override;
    void openorbital_scf() override;

    std::shared_ptr<CGHF> c1_deep_copy(std::shared_ptr<BasisSet> basis);

    // Unsure of what these are, but needed otherwise seg fault
    virtual bool same_a_b_orbs() const { return false; }
    virtual bool same_a_b_dens() const { return false; }

    // Empty functions for now -- sets up external potentials (TODO later)
    std::shared_ptr<UV> potential_;
    std::shared_ptr<VBase> V_potential() const override { return potential_; };

    // If DIIS is enabled (which it always should be), then it will update the orthogonalized Fock matrix Fp_
    std::complex<double> do_diis();

   protected:
    SharedMatrix V_mat;
    SharedMatrix S_mat;
    SharedMatrix T_mat;
    SharedMatrix G_mat;
    SharedMatrix F_mat;

    // DIIS variables
    // All 4 of these containers have a MAX size of DIIS_MAX_VECS
    // TODO ascertain if there's anything different between real and complex DIIS outside of the containers (e.g. error_doubles could be error_complex)
    std::deque<einsums::BlockTensor<std::complex<double>, 2>> Fdiis;     // Holds the grabbed Fock matrices to extrapolate
    std::deque<einsums::BlockTensor<std::complex<double>, 2>> err_vecs;  // Holds FDSmSDF_ at each iteration (orbital gradients) 
    std::vector<std::complex<double>> diis_coeffs;                       // Holds the coefficients for each Fock matrix in Fdiis
    std::vector<std::complex<double>> error_doubles;                     // RMS errors (real)

    double nuclearrep_; // Nuclear repulsion energy

    // NOTE: EINS_ and EINX_ are spin-blocked variants of S_ and X_ from HF, respectively
    // The change of variable names is intentional to avoid confusion
    einsums::BlockTensor<std::complex<double>, 2> F0_;        // Core Hamilton F0 = T + V
    einsums::BlockTensor<std::complex<double>, 2> F_;         // Non-orthogonal Fock matrix with F = T + V + J - K
    einsums::BlockTensor<std::complex<double>, 2> FDSmSDF_;   // Fock gradient FDSmSDF_ = [F, D]
    einsums::BlockTensor<std::complex<double>, 2> Fp_;        // Orthogonalized Fock matrix
    einsums::BlockTensor<double, 2> EINS_;                    // Spin-blocked overlap matrix -- it is never complex
    einsums::BlockTensor<std::complex<double>, 2> EINX_;      // Spin-blocked orthogonalization matrix -- also never complex, but cannot do real-complex matrix-matrix multiplication
    einsums::BlockTensor<std::complex<double>, 2> C_;         // Coefficient matrix built after back-trasnformation C' = XC
    
    // Cocc_ and cCocc_ deprecated since there doesn't appear to be a reason to permanently store these (temp1_ and temp2_ are used instead)
    //einsums::BlockTensor<std::complex<double>, 2> cCocc_; 
    //einsums::BlockTensor<std::complex<double>, 2> Cocc_;      // Occupied coefficient matrix -- needed for forming density matrix
    einsums::BlockTensor<std::complex<double>, 2> D_;         // 1-particle density matrix
    einsums::BlockTensor<std::complex<double>, 2> Fevecs_;    // Eigenvectors of Fock matrix
    einsums::BlockTensor<double, 1> Fevals_;                  // Eigenvalues of Fock matrix
    einsums::BlockTensor<std::complex<double>, 2> J_;         // Coulomb matrix
    einsums::BlockTensor<std::complex<double>, 2> K_;         // Exchange matrix   
    einsums::BlockTensor<std::complex<double>, 2> temp1_;     // temp1_ and temp2_ are temporary storage containers for intermediate steps
    einsums::BlockTensor<std::complex<double>, 2> temp2_;

    std::vector<int> irrep_sizes_;  // Since GHF is spin-blocked, each irrep (h) size will be 2*nsopi_[h]
    std::vector<int> nelecpi_;      // Number of electrons per irrep
};

}  // namespace scf
}  // namespace psi

#endif
