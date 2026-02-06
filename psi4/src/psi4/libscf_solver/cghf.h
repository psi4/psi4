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

    void common_init();

    void form_X();
    void form_S();

    void preiterations();
    void sap_guess();
    // void form_initial_C() override;
    void compute_SAD_guess(bool natorb) override;
    void zero_tensors();
    void finalize() override;
    void save_density_and_energy() override;
    void form_FDSmSDF();
    double compute_Dnorm();
    void form_V() override;
    void form_G() override;
    void form_F() override;
    void form_C(double shift) override;
    void form_D() override;
    void set_init_D();
    void redo_SCF();
    std::tuple<SharedMatrix, SharedMatrix> einsums_to_numpy(std::string mat_str);
    void form_numpy_D();

    void damping_update(double damping_percentage) override;
    double compute_initial_E() override;
    double compute_E() override;
    void setup_potential() override;
    void openorbital_scf() override;

    std::shared_ptr<CGHF> c1_deep_copy(std::shared_ptr<BasisSet> basis);
    virtual bool same_a_b_orbs() const { return false; }
    virtual bool same_a_b_dens() const { return false; }
    std::shared_ptr<UV> potential_;
    std::shared_ptr<VBase> V_potential() const override { return potential_; };

    std::complex<double> do_diis();

   protected:
    SharedMatrix V_mat;
    SharedMatrix S_mat;
    SharedMatrix T_mat;
    SharedMatrix G_mat;
    SharedMatrix F_mat;

    std::deque<einsums::BlockTensor<std::complex<double>, 2>> Fdiis;
    std::deque<einsums::BlockTensor<std::complex<double>, 2>> err_vecs;
    einsums::BlockTensor<std::complex<double>, 2> F_vecs;
    einsums::BlockTensor<std::complex<double>, 2> e_vecs;
    std::vector<std::complex<double>> diis_coeffs;
    std::vector<std::complex<double>> error_doubles;

    double nuclearrep_;

    einsums::BlockTensor<std::complex<double>, 2> F0_;
    einsums::BlockTensor<std::complex<double>, 2> EINT_;
    einsums::BlockTensor<std::complex<double>, 2> F_;
    einsums::BlockTensor<std::complex<double>, 2> FDSmSDF_;
    einsums::BlockTensor<std::complex<double>, 2> Fp_;
    einsums::BlockTensor<double, 2> EINS_;
    einsums::BlockTensor<std::complex<double>, 2> EINX_;
    einsums::BlockTensor<std::complex<double>, 2> C_;
    einsums::BlockTensor<std::complex<double>, 2> cCocc_;
    einsums::BlockTensor<std::complex<double>, 2> Cocc_;
    einsums::BlockTensor<std::complex<double>, 2> D_;
    einsums::BlockTensor<std::complex<double>, 2> Fevecs_;
    einsums::BlockTensor<double, 1> Fevals_;
    // einsums::BlockTensor<std::complex<double>, 1> Fevals_;

    einsums::BlockTensor<double, 1> RealEvals_;
    einsums::BlockTensor<std::complex<double>, 2> J_;
    einsums::BlockTensor<std::complex<double>, 2> K_;
    einsums::BlockTensor<std::complex<double>, 2> temp1_;
    einsums::BlockTensor<std::complex<double>, 2> temp2_;

    // size_t nirrep_;
    std::vector<int> irrep_sizes_;
    std::vector<int> nelecpi_;
    // Dimension nsopi_;
    // Dimension nbetapi_;
    // Dimension nelecpi_;
    // Dimension nvirtpi_;

    // std::shared_ptr<psi::BasisSet> basisset_;

    // einsums::BlockTensor<std::complex<double>, 2> F0_;
};

}  // namespace scf
}  // namespace psi

#endif
