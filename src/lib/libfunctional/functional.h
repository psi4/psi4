#ifndef functional_h
#define functional_h

/**********************************************************
* functional.h: declarations for functionals for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/
#include <psi4-dec.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace psi;
using namespace boost;
using namespace std;

namespace psi {

class Properties;

namespace functional {

class Functional {
    public:

        static boost::shared_ptr<Functional> createFunctional(const std::string & name, int npoints = 5000, int deriv = 1);
        static std::string availableFunctionals();
        static std::vector<std::string> availableNames();
        static std::string testFunctionals();

        Functional(int npoints = 5000, int deriv = 1);
        ~Functional();

        std::string getName() const { return name_; }        
        std::string getDescription() const { return description_; }        
        std::string getCitation() const { return citation_; }        
        std::string getParametersString();        
        std::vector<std::pair<std::string,double> > getParameters() const { return params_; }        
        void setName(const std::string & name) { name_ = name; }
        void setDescription(const std::string & description) { description_ = description; }
        void setCitation(const std::string & citation) { citation_ = citation; }
        void setParameters(const std::vector<std::pair<std::string,double> > & params) { params_ = params; }        
        void setParameter(const std::string & key, double value);

        std::string testFunctional(shared_ptr<Properties> props);

        bool isGGA() const { return is_gga_; } 
        bool isMeta() const { return is_meta_; } 

        void setNPoints(int npoints) { reallocate(npoints,deriv_); }
        void setDeriv(int deriv) { reallocate(npoints_,deriv); }
        int getNPoints() const { return npoints_; }
        int getDeriv() const { return deriv_; }

        double getDensityCutoff() const { return cutoff_; } 
        void setDensityCutoff(double cutoff) {cutoff = cutoff; }
        
        virtual void computeRKSFunctional(boost::shared_ptr<Properties> props) {};
        virtual void computeUKSFunctional(boost::shared_ptr<Properties> props) {};

        double* getFunctional() const { return functional_; }       
        double* getV_RhoA() const { return v_rho_a_; }       
        double* getV_RhoB() const { return v_rho_b_; }       
        double* getV_RhoA_RhoA() const { return v_rho_a_rho_a_; }       
        double* getV_RhoA_RhoB() const { return v_rho_a_rho_b_; }       
        double* getV_RhoB_RhoB() const { return v_rho_b_rho_b_; }       
        double* getV_GammaAA() const { return v_gamma_aa_; }       
        double* getV_GammaAB() const { return v_gamma_ab_; }       
        double* getV_GammaBB() const { return v_gamma_bb_; }       
        double* getV_RhoA_GammaAA() const { return v_rho_a_gamma_aa_; }       
        double* getV_RhoA_GammaAB() const { return v_rho_a_gamma_ab_; }       
        double* getV_RhoA_GammaBB() const { return v_rho_a_gamma_bb_; }       
        double* getV_RhoB_GammaAA() const { return v_rho_b_gamma_aa_; }       
        double* getV_RhoB_GammaAB() const { return v_rho_b_gamma_ab_; }       
        double* getV_RhoB_GammaBB() const { return v_rho_b_gamma_bb_; }       
        double* getV_GammaAA_GammaAA() const { return v_gamma_aa_gamma_aa_; }       
        double* getV_GammaAA_GammaAB() const { return v_gamma_aa_gamma_ab_; }       
        double* getV_GammaAA_GammaBB() const { return v_gamma_aa_gamma_bb_; }       
        double* getV_GammaAB_GammaAB() const { return v_gamma_ab_gamma_ab_; }       
        double* getV_GammaAB_GammaBB() const { return v_gamma_ab_gamma_bb_; }       
        double* getV_GammaBB_GammaBB() const { return v_gamma_bb_gamma_bb_; }       
        double* getV_TauA() const { return v_tau_a_; }
        double* getV_TauB() const { return v_tau_b_; }
        double* getV_RhoA_TauA() const { return v_rho_a_tau_a_; }
        double* getV_RhoA_TauB() const { return v_rho_a_tau_b_; }
        double* getV_RhoB_TauA() const { return v_rho_b_tau_a_; }
        double* getV_RhoB_TauB() const { return v_rho_b_tau_b_; }
        double* getV_GammaAA_TauA() const { return v_gamma_aa_tau_a_; }
        double* getV_GammaAA_TauB() const { return v_gamma_aa_tau_b_; }
        double* getV_GammaAB_TauA() const { return v_gamma_ab_tau_a_; }
        double* getV_GammaAB_TauB() const { return v_gamma_ab_tau_b_; }
        double* getV_GammaBB_TauA() const { return v_gamma_bb_tau_a_; }
        double* getV_GammaBB_TauB() const { return v_gamma_bb_tau_b_; }
        double* getV_TauA_TauA() const { return v_tau_a_tau_a_; }
        double* getV_TauA_TauB() const { return v_tau_a_tau_b_; }
        double* getV_TauB_TauB() const { return v_tau_b_tau_b_; }
 
     protected:
        
        void reallocate(int npoints, int deriv);
        void release();
        void allocate();
        
        std::string name_;
        std::string description_;
        std::string citation_;
        std::vector<std::pair<std::string, double> > params_;

        bool is_gga_;
        bool is_meta_;
        
        double cutoff_;

        double* functional_;
        double* v_rho_a_;
        double* v_rho_b_;
        double* v_rho_a_rho_a_;
        double* v_rho_a_rho_b_;
        double* v_rho_b_rho_b_;
        double* v_gamma_aa_;
        double* v_gamma_ab_;
        double* v_gamma_bb_;
        double* v_rho_a_gamma_aa_;
        double* v_rho_a_gamma_ab_;
        double* v_rho_a_gamma_bb_;
        double* v_rho_b_gamma_aa_;
        double* v_rho_b_gamma_ab_;
        double* v_rho_b_gamma_bb_;
        double* v_gamma_aa_gamma_aa_;
        double* v_gamma_aa_gamma_ab_;
        double* v_gamma_aa_gamma_bb_;
        double* v_gamma_ab_gamma_ab_;
        double* v_gamma_ab_gamma_bb_;
        double* v_gamma_bb_gamma_bb_;
        double* v_tau_a_;
        double* v_tau_b_;
        double* v_rho_a_tau_a_;
        double* v_rho_a_tau_b_;
        double* v_rho_b_tau_a_;
        double* v_rho_b_tau_b_;
        double* v_gamma_aa_tau_a_;
        double* v_gamma_aa_tau_b_;
        double* v_gamma_ab_tau_a_;
        double* v_gamma_ab_tau_b_;
        double* v_gamma_bb_tau_a_;
        double* v_gamma_bb_tau_b_;
        double* v_tau_a_tau_a_;
        double* v_tau_a_tau_b_;
        double* v_tau_b_tau_b_;

        int npoints_;
        int deriv_;
};

}}

#endif 
