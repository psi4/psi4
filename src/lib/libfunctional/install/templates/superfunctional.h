#ifndef super_functional_h
#define super_functional_h

/**********************************************************
* superfunctional.h: declarations for super (aliased) functionals for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/
#include "functional.h"
#include "dispersion.h"
#include <psiconfig.h>
#include <psi4-dec.h>
#include <stdlib.h>
#include <string>
#include <vector>

namespace psi {

class RKSFunctions; 
class UKSFunctions; 

namespace functional {

class SuperFunctional {
    public:

        static boost::shared_ptr<SuperFunctional> createSuperFunctional(const std::string & alias, int npoints = 5000, int deriv = 1);
        static boost::shared_ptr<SuperFunctional> buildSuperFunctional(const std::string & build, int npoints = 5000, int deriv = 1);
        static std::string availableSuperFunctionals();
        static std::vector<std::string> availableNames();
        //static std::string testSuperFunctionals();

        SuperFunctional(int npoints = 5000, int deriv = 1);
        ~SuperFunctional();

        std::string getName() const { return name_; }
        std::string getDescription() const { return description_; }
        std::string getCitation() const { return citation_; }
        void setName(const std::string & name) { name_ = name; }
        void setDescription(const std::string & description) { description_ = description; }
        void setCitation(const std::string & citation) { citation_ = citation; }

        void print(FILE* out = outfile, int level = 2) const {}

        std::string getComposition();

        bool isGGA();
        bool isMeta();

        /// 0 - LSDA, 1 - GGA, 2 - Meta-GGA
        int getLocalAnsatz();

        bool isHybrid() { return  exact_exchange_ > 0.0; }
        bool isDoubleHybrid() { return pt2_ > 0.0; }
        bool isRangeCorrected() { return omega_ > 0.0; }
        bool isDashD() { return dashD_weight_ > 0.0; }

        void setParameter(const std::string & functional, const std::string & param, double val);

        //std::string testSuperFunctional(boost::shared_ptr<Properties> props);

        double getExactExchange() { return exact_exchange_; }
        void setExactExchange(double exch) {exact_exchange_ = exch; }

        double getOmega() { return omega_; }
        void setOmega(double omega);

        double getPT2() { return pt2_; }
        void setPT2(double pt2) {pt2_ = pt2; }

        double getDashDWeight() { return dashD_weight_; }
        boost::shared_ptr<Dispersion> getDashD() { return dashD_; }
        void setDashD(boost::shared_ptr<Dispersion> disp, double weight);

        void setNPoints(int npoints) { reallocate(npoints,deriv_); }
        void setDeriv(int deriv) { reallocate(npoints_,deriv); }
        int getNPoints() const { return npoints_; }
        int getDeriv() const { return deriv_; }

        void computeRKSFunctional(boost::shared_ptr<RKSFunctions> props);
        void computeUKSFunctional(boost::shared_ptr<UKSFunctions> props);
        void collectResults();

        void addFunctional(const boost::shared_ptr<Functional> & f, double weight);
        void addFunctional(const std::string & name, double weight);
        boost::shared_ptr<Functional>& operator[](int index) { return functionals_[index].first; }
        boost::shared_ptr<Functional> getFunctional(int index) { return functionals_[index].first; }
        double getWeight(int index) { return functionals_[index].second; }
        int size() { return functionals_.size(); }

        double* getFunctionalValue() const { return functional_; }
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

        std::vector<std::pair<boost::shared_ptr<Functional>, double> > functionals_;

        std::string name_;
        std::string description_;
        std::string citation_;

        double exact_exchange_;
        double omega_;
        double pt2_;
        double dashD_weight_;
        boost::shared_ptr<Dispersion> dashD_;

        double* restrict functional_;
        double* restrict v_rho_a_;
        double* restrict v_rho_b_;
        double* restrict v_rho_a_rho_a_;
        double* restrict v_rho_a_rho_b_;
        double* restrict v_rho_b_rho_b_;
        double* restrict v_gamma_aa_;
        double* restrict v_gamma_ab_;
        double* restrict v_gamma_bb_;
        double* restrict v_rho_a_gamma_aa_;
        double* restrict v_rho_a_gamma_ab_;
        double* restrict v_rho_a_gamma_bb_;
        double* restrict v_rho_b_gamma_aa_;
        double* restrict v_rho_b_gamma_ab_;
        double* restrict v_rho_b_gamma_bb_;
        double* restrict v_gamma_aa_gamma_aa_;
        double* restrict v_gamma_aa_gamma_ab_;
        double* restrict v_gamma_aa_gamma_bb_;
        double* restrict v_gamma_ab_gamma_ab_;
        double* restrict v_gamma_ab_gamma_bb_;
        double* restrict v_gamma_bb_gamma_bb_;
        double* restrict v_tau_a_;
        double* restrict v_tau_b_;
        double* restrict v_rho_a_tau_a_;
        double* restrict v_rho_a_tau_b_;
        double* restrict v_rho_b_tau_a_;
        double* restrict v_rho_b_tau_b_;
        double* restrict v_gamma_aa_tau_a_;
        double* restrict v_gamma_aa_tau_b_;
        double* restrict v_gamma_ab_tau_a_;
        double* restrict v_gamma_ab_tau_b_;
        double* restrict v_gamma_bb_tau_a_;
        double* restrict v_gamma_bb_tau_b_;
        double* restrict v_tau_a_tau_a_;
        double* restrict v_tau_a_tau_b_;
        double* restrict v_tau_b_tau_b_;

        int npoints_;
        int deriv_;
        bool oldGGA_;
        bool oldMeta_;

};

}}

#endif
