#ifndef _psi_src_lib_libmints_basispoints_h_
#define _psi_src_lib_libmints_basispoints_h_

#include <boost/shared_ptr.hpp>

namespace psi {

class BasisSet;
class GridBlock;

class BasisPoints
{
    protected:
        boost::shared_ptr<BasisSet> basis_;
        double* cutoff_radii_2_;
        double cutoff_epsilon_;
        bool* sig_shells_;
        bool* false_shells_;
        int* rel2abs_shells_;
        int* abs2rel_shells_;
        int nsig_shells_;
        int* rel2abs_functions_;
        int* abs2rel_functions_;
        int nsig_functions_;
        double block_efficiency_;

        double* ao_points_; 
        double* ao_gradX_, *ao_gradY_, *ao_gradZ_;
        double* ao_hessXY_, *ao_hessXZ_, *ao_hessYZ_, *ao_hessXX_, *ao_hessYY_, *ao_hessZZ_;
        double* ao_laplac_; 

        bool spherical_;
        int max_am_;       
        int max_carts_;       
        int max_prim_;       
        int nao_;       
        int nso_;       
        int nshell_;       

        int* nprim_per_shell_;
        int* ncart_per_shell_;
        int* nso_per_shell_;
        int* am_per_shell_;

        int* sotrans_count_;
        int** sotrans_so_;
        int** sotrans_ao_;
        double** sotrans_coef_;
 
        int* a_;       
        int* b_;       
        int* c_;       
        double* xc_;       
        double* yc_;       
        double* zc_;       
        double** coef_;
        double** alpha_;
        double** prims_;
        double* ex_;
        double* dx_pow_;
        double* dy_pow_;
        double* dz_pow_;
 
        double** points_;
        double** gradientsX_;
        double** gradientsY_;
        double** gradientsZ_;
        double** hessiansXY_;
        double** hessiansXZ_;
        double** hessiansYZ_;
        double** hessiansXX_;
        double** hessiansYY_;
        double** hessiansZZ_;
        double** laplacians_;

        bool do_points_;
        bool do_gradients_;
        bool do_hessians_;
        bool do_laplacians_;
        bool have_points_;
        bool have_gradients_;
        bool have_hessians_;
        bool have_laplacians_;

        int block_size_;
        int true_size_;

        void release();
        void allocate();
        void init_basis_data();
    public:

        BasisPoints(boost::shared_ptr<BasisSet> b, int block_size);
        virtual ~BasisPoints();

        void computePoints(boost::shared_ptr<GridBlock> block);

        double** getPoints() const { return points_; }
        double** getGradientsX() const {return gradientsX_; }
        double** getGradientsY() const {return gradientsY_; }
        double** getGradientsZ() const {return gradientsZ_; }
        double** getHessiansXY() const {return hessiansXY_; }
        double** getHessiansXZ() const {return hessiansXZ_; }
        double** getHessiansYZ() const {return hessiansYZ_; }
        double** getHessiansXX() const {return hessiansXX_; }
        double** getHessiansYY() const {return hessiansYY_; }
        double** getHessiansZZ() const {return hessiansZZ_; } 
        double** getLaplacians() const {return laplacians_; }

        double* getCutoffRadii2() const {return cutoff_radii_2_; }
        void setCutoffEpsilon(double epsilon = 0.0) { computeCutoffRadii2(epsilon); }
        double getCutoffEpsilon() const { return cutoff_epsilon_; }
        void computeCutoffRadii2(double epsilon = 0.0);
        void computeSignificantShells(boost::shared_ptr<GridBlock> block);

        bool* getSignificantShells() const {return sig_shells_; }
        int* abs2relShells() const {return abs2rel_shells_;}
        int* rel2absShells() const {return rel2abs_shells_;}
        int* abs2relFunctions() const {return abs2rel_functions_;}
        int* rel2absFunctions() const {return rel2abs_functions_;}
        double getBlockEfficiency() const {return block_efficiency_;}

        void setToComputePoints(bool val) {do_points_ = val; allocate(); release(); }
        void setToComputeGradients(bool val) { do_gradients_ = val; allocate(); release(); }
        void setToComputeHessians(bool val) { do_hessians_ = val; allocate(); release(); }
        void setToComputeLaplacians(bool val) { do_laplacians_ = val; allocate(); release(); }

        int getBlockSize() const { return block_size_;}  
        int getTrueSize() const {return true_size_; }
        int nbf(); 
        int nSignificantShells() const {return nsig_shells_;}
        int nSignificantFunctions() const {return nsig_functions_;}

        void setTrueSize(int size) {true_size_ = size; }
};
}
#endif
