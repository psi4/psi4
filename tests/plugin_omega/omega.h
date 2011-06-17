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
class OmegaV;
class OmegaWavefunction;
class OmegaKS;
class OmegaIPKS;

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
    boost::shared_ptr<Matrix> Ca_;
    boost::shared_ptr<Matrix> Cb_;
    boost::shared_ptr<Matrix> Da_;
    boost::shared_ptr<Matrix> Db_;
    boost::shared_ptr<Matrix> Dt_;
    boost::shared_ptr<Matrix> Dtold_;
    boost::shared_ptr<Matrix> J_;
    boost::shared_ptr<Matrix> Ka_;
    boost::shared_ptr<Matrix> Kb_;
    boost::shared_ptr<Matrix> wKa_;
    boost::shared_ptr<Matrix> wKb_;
    boost::shared_ptr<Matrix> Va_;
    boost::shared_ptr<Matrix> Vb_;
    boost::shared_ptr<Matrix> Ga_;
    boost::shared_ptr<Matrix> Gb_;
    boost::shared_ptr<Matrix> Fa_;
    boost::shared_ptr<Matrix> Fb_;

    // Common hooks
    boost::shared_ptr<OmegaDF> df_;
    boost::shared_ptr<OmegaV> ks_;

    boost::shared_ptr<Matrix> S_;
    boost::shared_ptr<Matrix> X_;
    boost::shared_ptr<Matrix> H_;

    void common_init();
    void save_info();
    void form_G();
    void form_V();
    void form_F();
    bool diis();
    void compute_E();
    void form_C();
    void form_D();

    boost::shared_ptr<Matrix> form_FDSmSDF(boost::shared_ptr<Matrix> F, boost::shared_ptr<Matrix> D);

public:
    OmegaWavefunction(Options& options,
                      boost::shared_ptr<PSIO> psio,
                      boost::shared_ptr<BasisSet> primary,
                      boost::shared_ptr<Matrix> Ca, int na,
                      boost::shared_ptr<Matrix> Cb, int nb,
                      boost::shared_ptr<Matrix> S,
                      boost::shared_ptr<Matrix> X, boost::shared_ptr<Matrix> H,
                      boost::shared_ptr<OmegaDF> df, boost::shared_ptr<OmegaV> ks);
    virtual ~OmegaWavefunction();
    
    // Take an SCF step
    std::string step();
    // Reset DIIS and whatever
    void reset();
    // Clear, without resetting DIIS
    void clear();

    // Perform a guess using the given KS matrices
    void guess(boost::shared_ptr<Matrix> Fa, boost::shared_ptr<Matrix> Fb);

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
    boost::shared_ptr<Matrix> Fa() const { return Fa_; } 
    boost::shared_ptr<Matrix> Fb() const { return Fb_; } 
    boost::shared_ptr<Matrix> Ca() const { return Ca_; } 
    boost::shared_ptr<Matrix> Cb() const { return Cb_; } 
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
    static boost::shared_ptr<Matrix> build_X(boost::shared_ptr<BasisSet> primary, double min_S);
    static boost::shared_ptr<Matrix> build_H(boost::shared_ptr<BasisSet> primary);
    static boost::shared_ptr<Matrix> build_S(boost::shared_ptr<BasisSet> primary);

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
    boost::shared_ptr<Matrix> H_;
    boost::shared_ptr<Matrix> X_;
    boost::shared_ptr<Matrix> S_;
    
    // Common interelectronic bits
    boost::shared_ptr<OmegaDF> df_;
    // Common KS bits
    boost::shared_ptr<OmegaV> ks_;   

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
    boost::shared_ptr<Matrix> Fa_l_N_; 
    boost::shared_ptr<Matrix> Fa_l_M_; 
    boost::shared_ptr<Matrix> Fb_l_N_; 
    boost::shared_ptr<Matrix> Fb_l_M_; 
    boost::shared_ptr<Matrix> Fa_r_N_; 
    boost::shared_ptr<Matrix> Fa_r_M_; 
    boost::shared_ptr<Matrix> Fb_r_N_; 
    boost::shared_ptr<Matrix> Fb_r_M_; 

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
    boost::shared_ptr<Matrix> Cmn_;
    // (mn|A) 
    boost::shared_ptr<Matrix> Amn_;
    // (mn|erf(wr)/r|A) 
    boost::shared_ptr<Matrix> Wmn_;
    
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

    boost::shared_ptr<Matrix> J(boost::shared_ptr<Matrix> D); 
    boost::shared_ptr<Matrix> K(boost::shared_ptr<Matrix> C, int nocc); 
    boost::shared_ptr<Matrix> wK(boost::shared_ptr<Matrix> C, int nocc); 
};

class OmegaV {

protected:
    boost::shared_ptr<PSIO> psio_;
    boost::shared_ptr<BasisSet> primary_;
    boost::shared_ptr<functional::SuperFunctional> functional_;
    boost::shared_ptr<Integrator> integrator_;
    boost::shared_ptr<Properties> properties_;
    
    std::map<std::string, double> quad_values_;
    boost::shared_ptr<Matrix> Va_;
    boost::shared_ptr<Matrix> Vb_;

    void common_init();

public:
    OmegaV(boost::shared_ptr<PSIO> psio,
        boost::shared_ptr<BasisSet> primary,
        boost::shared_ptr<functional::SuperFunctional> functional,
        boost::shared_ptr<Integrator> integrator,
        boost::shared_ptr<Properties> properties);
    virtual ~OmegaV();    

    // Set the omega value
    void set_omega(double omega);       

    void form_V(boost::shared_ptr<Matrix> Da, boost::shared_ptr<Matrix> Ca, int na,
                boost::shared_ptr<Matrix> Db, boost::shared_ptr<Matrix> Cb, int nb); 

    boost::shared_ptr<functional::SuperFunctional> functional() const { return functional_; }

    boost::shared_ptr<Matrix> Va() const { return Va_; }
    boost::shared_ptr<Matrix> Vb() const { return Vb_; }
    double Exc();
    double variable(const std::string& key); 

};

}} // Namespaces

#endif

