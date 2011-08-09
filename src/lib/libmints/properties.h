#ifndef _psi_src_lib_libmints_properties_h_
#define _psi_src_lib_libmints_properties_h_

#include <boost/shared_ptr.hpp>
#include <libmints/basispoints.h>

namespace psi {

class BasisSet;
class Matrix;
class Vector3;
class GridBlock;

class Properties : public BasisPoints
{
    protected:
        // Local copies of D and Cocc matrices
        // to avoid OTF malloc
        double** Da_;
        double** Db_;
        double** Ca_;
        double** Cb_;

        // temporary scratch tensor for DGEMM
        double** temp_tens_;

        // Fundamental DFT variables
        double* rho_a_;
        double* rho_b_;
        double* gamma_aa_;
        double* gamma_ab_;
        double* gamma_bb_;
        double* tau_a_;
        double* tau_b_;

        // Some other gradients required for DFT
        double* rho_a_x_;
        double* rho_a_y_;
        double* rho_a_z_;
        double* rho_b_x_;
        double* rho_b_y_;
        double* rho_b_z_;

        // Indicator variables
        bool do_density_;
        bool do_density_gradient_;
        bool do_ke_density_;
    public:
        static boost::shared_ptr<Properties> constructProperties(boost::shared_ptr<BasisSet> b, int block_size)
        {
            return (boost::shared_ptr<Properties>)new Properties(b,block_size);
        }
        static boost::shared_ptr<Properties> get_testbed();

        Properties(boost::shared_ptr<BasisSet> b, int block_size);
        virtual ~Properties();

        void computeRKSProperties(boost::shared_ptr<GridBlock> grid, boost::shared_ptr<Matrix> D, boost::shared_ptr<Matrix> C = boost::shared_ptr<Matrix>(), int* docc = NULL );
        void computeUKSProperties(boost::shared_ptr<GridBlock> grid, boost::shared_ptr<Matrix> Da, boost::shared_ptr<Matrix> Db, boost::shared_ptr<Matrix> Ca = boost::shared_ptr<Matrix>(), boost::shared_ptr<Matrix> Cb = boost::shared_ptr<Matrix>(), int* Na = NULL, int* Nb = NULL );

        double* getRhoA() const { return rho_a_; }
        double* getRhoB() const { return rho_b_; }
        double* getGammaAA() const { return gamma_aa_; }
        double* getGammaAB() const { return gamma_ab_; }
        double* getGammaBB() const { return gamma_bb_; }
        double* getTauA() const { return tau_a_; }
        double* getTauB() const { return tau_b_; }

        double* getRhoAX() const { return rho_a_x_; }
        double* getRhoAY() const { return rho_a_y_; }
        double* getRhoAZ() const { return rho_a_z_; }
        double* getRhoBX() const { return rho_b_x_; }
        double* getRhoBY() const { return rho_b_y_; }
        double* getRhoBZ() const { return rho_b_z_; }

        double** getScratch() const { return temp_tens_; }

        void setToComputeDensity(bool v);
        void setToComputeDensityGradient(bool v);
        void setToComputeKEDensity(bool v);
};
typedef boost::shared_ptr<Properties> SharedProperties;
}
#endif
