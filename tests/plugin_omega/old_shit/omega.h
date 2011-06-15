#ifndef OMEGA_PROC_H
#define OMEGA_PROC_H

#include <libscf_solver/ks.h>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

namespace scf {

class OmegaWavefunction : public Wavefunction {

protected:

double E_;
double Eold_;
double Drms_;

public:
    OmegaWavefunction();

};

class OmegaKS : public UKS {

protected:
 
    double E_N_;
    double Eold_N_;
    double Drms_N_;

    boost::shared_ptr<Vector> epsilon_a_N_;
    boost::shared_ptr<Vector> epsilon_b_N_;

    boost::shared_ptr<Matrix> Ca_N_;
    boost::shared_ptr<Matrix> Cb_N_;
    boost::shared_ptr<Matrix> Da_N_;
    boost::shared_ptr<Matrix> Db_N_;
    boost::shared_ptr<Matrix> Dt_N_;
    boost::shared_ptr<Matrix> Dtold_N_;
    boost::shared_ptr<Matrix> J_N_;
    boost::shared_ptr<Matrix> Ka_N_;
    boost::shared_ptr<Matrix> Kb_N_;
    boost::shared_ptr<Matrix> wKa_N_;
    boost::shared_ptr<Matrix> wKb_N_;
    boost::shared_ptr<Matrix> Va_N_;
    boost::shared_ptr<Matrix> Vb_N_;
    boost::shared_ptr<Matrix> Ga_N_;
    boost::shared_ptr<Matrix> Gb_N_;
    boost::shared_ptr<Matrix> Fa_N_;
    boost::shared_ptr<Matrix> Fb_N_;

    int nalpha_N_;
    int nbeta_N_;

    int doccpi_N_[8];
    int soccpi_N_[8];
    int nalphapi_N_[8];
    int nbetapi_N_[8];
 
    void common_init(); 

    void initialize_N();
    void find_occupation_N();
    void save_density_and_energy_N();
    bool test_convergency_N();
    double compute_E_N();
    void form_G_N();
    void form_F_N();
    void form_C_N();
    void form_D_N();
    void form_V_N();

    virtual void finalize();

public:
    OmegaKS(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    OmegaKS(Options& options, boost::shared_ptr<PSIO> psio);
    virtual ~OmegaKS();   

    // Override HF::compute_energy to launch custom procedure 
    virtual double compute_energy();
};


}} // Namespaces

#endif

