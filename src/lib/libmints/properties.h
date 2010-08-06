#ifndef _psi_src_lib_libmints_properties_h_
#define _psi_src_lib_libmints_properties_h_

#include <boost/shared_ptr.hpp>

#define RKS_GGA_TESTBED_SIZE_  9
namespace psi {

    class BasisSet;
    class BasisPoints;
    class Matrix;
    class SimpleMatrix;
    class Vector3;
    class GridBlock;
    class ElectrostaticInt;

	class Properties:public BasisPoints
	{
		protected:
			double** mos_;
			int* mo_inds_;
			int nmo_;
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
			bool do_mos_;
			bool do_density_;
			bool do_density_gradient_;
			bool do_density_hessian_;
			bool do_density_laplacian_;
			bool do_ke_density_;
			bool do_electrostatic_;
		    boost::shared_ptr<ElectrostaticInt> e_ints_;
                public:
			static Properties * constructProperties(boost::shared_ptr<BasisSet> b, int block_size)
			{
				return new Properties(b,block_size);
			}
			Properties(boost::shared_ptr<BasisSet> b, int block_size);
			~Properties();
			void computeProperties(boost::shared_ptr<GridBlock> grid, boost::shared_ptr<Matrix> D, boost::shared_ptr<Matrix> C = boost::shared_ptr<Matrix>() );
                        void get_RKS_GGA_Testbed();
                        int get_RKS_GGA_Testbed_Size() const {return RKS_GGA_TESTBED_SIZE_; }
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
