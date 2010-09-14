#ifndef _psi_src_lib_libmints_properties_h_
#define _psi_src_lib_libmints_properties_h_

#include <boost/shared_ptr.hpp>
#include <libmints/basispoints.h>

namespace psi {

class BasisSet;
class Matrix;
class SimpleMatrix;
class Vector3;
class GridBlock;
class ElectrostaticInt;

class Properties : public BasisPoints
{
    protected:
        double** mos_;
        double* density_;
        double* densityX_;
        double* densityY_;
        double* densityZ_;
        double* densityXY_;
        double* densityXZ_;
        double* densityYZ_;
        double* densityXX_;
        double* densityYY_;
        double* densityZZ_;
        double* density_gradient_2_;
        double* density_laplacian_;
        double* ke_density_;
        double* electrostatic_;
       
        double* rho_a_;
        double* rho_b_;
        double* gamma_aa_;
        double* gamma_ab_;
        double* gamma_bb_;
        double* tau_a_;
        double* tau_b_;
        
        int* mo_inds_;
        int nmo_;
        
        bool do_mos_;
        bool do_density_;
        bool do_density_gradient_;
        bool do_density_hessian_;
        bool do_density_laplacian_;
        bool do_ke_density_;
        bool do_electrostatic_;
        boost::shared_ptr<ElectrostaticInt> e_ints_;
    public:
        static boost::shared_ptr<Properties> constructProperties(boost::shared_ptr<BasisSet> b, int block_size)
        {
            return (boost::shared_ptr<Properties>)new Properties(b,block_size);
        }
        Properties(boost::shared_ptr<BasisSet> b, int block_size);
        virtual ~Properties();
        void computeProperties(boost::shared_ptr<GridBlock> grid, boost::shared_ptr<Matrix> D, boost::shared_ptr<Matrix> C = boost::shared_ptr<Matrix>() );
        static boost::shared_ptr<Properties> get_testbed();
        const double* getDensity() const { return density_; }		
        const double* getDensityX() const { return densityX_; }		
        const double* getDensityY() const { return densityY_; }		
        const double* getDensityZ() const { return densityZ_; }		
        const double* getDensityGradientSquared() const { return density_gradient_2_; }
        const double* getDensityXY() const { return densityXY_; }		
        const double* getDensityXZ() const { return densityXZ_; }		
        const double* getDensityYZ() const { return densityYZ_; }		
        const double* getDensityXX() const { return densityXX_; }		
        const double* getDensityYY() const { return densityYY_; }		
        const double* getDensityZZ() const { return densityZZ_; }		
        const double* getDensityLaplacian() const { return density_laplacian_; }
        const double* getKEDensity() const { return ke_density_; }	
        const double* getElectrostatic() const {return electrostatic_; }
        const double* getDensityA() const { return rho_a_; }	
        const double* getDensityB() const { return rho_b_; }	
        const double* getDensityGradientSquaredAA() const { return gamma_aa_; }
        const double* getDensityGradientSquaredAB() const { return gamma_ab_; }
        const double* getDensityGradientSquaredBB() const { return gamma_bb_; }
        const double* getKEDensityA() const { return tau_a_; }	
        const double* getKEDensityB() const { return tau_b_; }	
        double** getMOs() {return mos_; }
        void setToComputeMOs(bool v, int* indices, int n);
        void setToComputeDensity(bool v);
        void setToComputeDensityGradient(bool v);	
        void setToComputeDensityHessian(bool v);	
        void setToComputeDensityLaplacian(bool v);	
        void setToComputeKEDensity(bool v);	
        void setToComputeElectrostatic(bool v);	
};
typedef boost::shared_ptr<Properties> SharedProperties;
}
#endif
