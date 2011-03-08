#ifndef dropbox_h
#define dropbox_h

/**********************************************************
* dropbox.h: declarations for partials storage for KS-DFT
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
namespace functional {

class DropBox {
    public:

        DropBox(int npoints = 5000, int deriv = 1);
        ~DropBox();

        bool isGGA();    
        bool isMeta();    
        
        void setNPoints(int npoints) { reallocate(npoints,deriv_); }
        void setDeriv(int deriv) { reallocate(npoints_,deriv); }
        int getNPoints() const { return npoints_; }
        int getDeriv() const { return deriv_; }

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
        void zero();
 
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
        bool oldGGA_;
        bool oldMeta_;

};

}}

#endif 
