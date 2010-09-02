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

        static boost::shared_ptr<Functional> createFunctional(const std::string & name, int npoints, int deriv = 1);

        Functional(int npoints, int deriv = 1);
        virtual void init() = 0;
        void allocate();
        ~Functional();

        std::string getName() const { return name_; }        
        std::string getDescription() const { return name_; }        
        std::string getCitation() const { return citation_; }        
        std::string getParametersString();        
        std::vector<std::pair<std::string,double> > getParameters() const { return params_; }        
        void setName(const std::string & name) { name_ = name; }
        void setDescription(const std::string & description) { description_ = description; }
        void setCitation(const std::string & citation) { citation_ = citation; }
        void setParameters( const std::vector<std::pair<std::string,double> > & params) { params_ = params; }        

        bool isGGA() const { return is_gga_; } 
        bool isMeta() const { return is_meta_; } 
        bool isHybrid() const { return is_hybrid_; } 

        double getExactExchange() const { return exact_exchange_; } 
        void setExactExchange(double exch) {exact_exchange_ = exch; }
     
        int getNPoints() const { return npoints_; }
        int getDeriv() const { return deriv_; }

        virtual void computeRKSFunctional(boost::shared_ptr<Properties> props) = 0;
        virtual void computeUKSFunctional(boost::shared_ptr<Properties> props) = 0;

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
 
     protected:
        std::string name_;
        std::string description_;
        std::string citation_;
        std::vector<std::pair<std::string, double> > params_;

        bool is_gga_;
        bool is_meta_;
        bool is_hybrid_;

        double exact_exchange_;

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

        int npoints_;
        int deriv_;
};

}}

#endif 
