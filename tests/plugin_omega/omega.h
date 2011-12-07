#ifndef OMEGA_PROC_H
#define OMEGA_PROC_H

#include <libscf_solver/ks.h>
#include <boost/tuple/tuple.hpp>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

namespace scf {

class OmegaDF;
class OmegaWavefunction;
class OmegaKS;
class OmegaIPKS;
class UKSPotential;
class RKSPotential;

class OmegaWavefunction {

protected:

    int print_;
    int debug_;

    Options& options_; 
    boost::shared_ptr<PSIO> psio_;
    boost::shared_ptr<BasisSet> basisset_;

    double Enuc_;
    double Exc_;
    double E_;
    double Eold_;
    double Drms_;
    
    int nalpha_;
    int nbeta_;
    int nso_;
    int nmo_;

    int iteration_;
    int min_diis_vectors_;
    int max_diis_vectors_;
    int diis_start_;
    bool diis_enabled_;
    boost::shared_ptr<DIISManager> diis_manager_;

    boost::shared_ptr<Vector> epsilon_b_;
    boost::shared_ptr<Vector> epsilon_a_;
    SharedMatrix Ca_;
    SharedMatrix Cb_;
    SharedMatrix Da_;
    SharedMatrix Db_;
    SharedMatrix Dt_;
    SharedMatrix Dtold_;
    SharedMatrix J_;
    SharedMatrix Ka_;
    SharedMatrix Kb_;
    SharedMatrix wKa_;
    SharedMatrix wKb_;
    SharedMatrix Va_;
    SharedMatrix Vb_;
    SharedMatrix Ga_;
    SharedMatrix Gb_;
    SharedMatrix Fa_;
    SharedMatrix Fb_;

    // Common hooks
    boost::shared_ptr<OmegaDF> df_;
    boost::shared_ptr<UKSPotential> ks_;

    SharedMatrix S_;
    SharedMatrix X_;
    SharedMatrix H_;

    void common_init();
    void save_info();
    void form_G();
    void form_V();
    void form_F();
    bool diis();
    void compute_E();
    void form_C();
    void form_D();

    SharedMatrix form_FDSmSDF(SharedMatrix F, SharedMatrix D);

public:
    OmegaWavefunction(Options& options,
                      boost::shared_ptr<PSIO> psio,
                      boost::shared_ptr<BasisSet> primary,
                      SharedMatrix Ca, int na,
                      SharedMatrix Cb, int nb,
                      SharedMatrix S,
                      SharedMatrix X, SharedMatrix H,
                      boost::shared_ptr<OmegaDF> df, boost::shared_ptr<UKSPotential> ks);
    virtual ~OmegaWavefunction();
    
    // Take an SCF step
    std::string step();
    // Reset DIIS and whatever
    void reset();
    // Clear, without resetting DIIS
    void clear();

    // Perform a guess using the given KS matrices
    void guess(SharedMatrix Fa, SharedMatrix Fb);

    // Print the orbital energies
    void print_orbitals();

    double deltaE() const { return E_ - Eold_; }
    double deltaD() const { return Drms_; } 

    double E() const { return E_; }
    double koopmansIP();
    double koopmansEA();

    int nso() const { return nso_; }
    int nmo() const { return nmo_; }
    int nalpha() const { return nalpha_; }
    int nbeta() const { return nbeta_; }
    SharedMatrix Fa() const { return Fa_; } 
    SharedMatrix Fb() const { return Fb_; } 
    SharedMatrix Ca() const { return Ca_; } 
    SharedMatrix Cb() const { return Cb_; } 
};

class OmegaKS {

protected:
    Options& options_;
    boost::shared_ptr<PSIO> psio_;
  
    int block_size_; 
    boost::shared_ptr<functional::SuperFunctional> functional_;
 
    double energy_threshold_;
    double density_threshold_;
    
    int print_;
    int debug_;
    long int memory_;
    int nthread_;
    
    std::map<std::string, boost::shared_ptr<OmegaWavefunction> > wfns_;

    double initial_omega_;

    void common_init();
    static SharedMatrix build_X(boost::shared_ptr<BasisSet> primary, double min_S);
    static SharedMatrix build_H(boost::shared_ptr<BasisSet> primary);
    static SharedMatrix build_S(boost::shared_ptr<BasisSet> primary);

public:
    OmegaKS(Options&, boost::shared_ptr<PSIO>);
    virtual ~OmegaKS();

    virtual void run_procedure();

    virtual void print_header() = 0;
    virtual void guess_omega() = 0;
    virtual void form_H() = 0;
    virtual void form_X() = 0;
    virtual void form_DF() = 0;
    virtual void form_KS() = 0;
    virtual void populate() = 0;
    virtual void finalize() = 0;

    virtual void omega_step(int iter) = 0;
    virtual bool is_omega_converged() = 0;
};

class OmegaIPKS : public OmegaKS {

protected:

    boost::shared_ptr<Wavefunction> reference_;
    boost::shared_ptr<MatrixFactory> factory_;
    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<BasisSet> basisset_;
    boost::shared_ptr<BasisSet> auxiliary_;
    
    // Common matrices
    SharedMatrix H_;
    SharedMatrix X_;
    SharedMatrix S_;
    
    // Common interelectronic bits
    boost::shared_ptr<OmegaDF> df_;
    // Common KS bits
    boost::shared_ptr<UKSPotential> ks_;   

    void common_init();
 
    // Have we backeted the target omega yet?
    bool bracketed_;
 
    // History of omegas, with <omega, kIP, IP, status>
    std::vector<boost::tuple<double,double,double,std::string > > omega_trace_;

    // Used in modified Regula Falsi
    int left_;
    int right_;

    // Used in Brent's method
    double omega_a_;
    double omega_b_;
    double omega_c_;
    double omega_d_;
    double delta_a_;
    double delta_b_;
    double delta_c_;

    bool mflag_;

    // Left and right omegas
    double omega_l_;
    double omega_r_;

    // Delta is kIP - IP
    double delta_l_;
    double delta_r_;

    // Left and right F matrices
    SharedMatrix Fa_l_N_; 
    SharedMatrix Fa_l_M_; 
    SharedMatrix Fb_l_N_; 
    SharedMatrix Fb_l_M_; 
    SharedMatrix Fa_r_N_; 
    SharedMatrix Fa_r_M_; 
    SharedMatrix Fb_r_N_; 
    SharedMatrix Fb_r_M_; 

public:
    OmegaIPKS(Options&, boost::shared_ptr<PSIO>);
    virtual ~OmegaIPKS();

    void print_header();
    void guess_omega();
    void form_H();
    void form_X();
    void form_DF();
    void form_KS();
    void populate();
    void step();
    void finalize();

    void omega_step(int iter);
    bool is_omega_converged();
};

class OmegaDF {

protected:
    double omega_;
    // I need the SphericalTransforms from this. 
    boost::shared_ptr<IntegralFactory> factory_;
    std::vector<boost::shared_ptr<ErfERI> > erf_eri_;
    std::vector<const double*> erf_buffer_;
    int nthread_;

    boost::shared_ptr<BasisSet> primary_;
    boost::shared_ptr<BasisSet> auxiliary_;
    boost::shared_ptr<BasisSet> zero_;
    boost::shared_ptr<PSIO> psio_;   
 
    // (mn|A) J_AC ^ -1
    SharedMatrix Cmn_;
    // (mn|A) 
    SharedMatrix Amn_;
    // (mn|erf(wr)/r|A) 
    SharedMatrix Wmn_;
    
    void common_init();
    // Build the static integrals
    void build_static();
    // Build the dynamic integrals
    void build_dynamic();
   
public:
    OmegaDF(boost::shared_ptr<PSIO> psio,
        boost::shared_ptr<BasisSet> primary,
        boost::shared_ptr<BasisSet> auxiliary);
    virtual ~OmegaDF();

    // Set the omega value
    void set_omega(double omega);       

    SharedMatrix J(SharedMatrix D); 
    SharedMatrix K(SharedMatrix C, int nocc); 
    SharedMatrix wK(SharedMatrix C, int nocc); 
};


}} // Namespaces

#endif

