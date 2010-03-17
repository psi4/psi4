#ifndef _psi_src_lib_libmints_basispoints_h_
#define _psi_src_lib_libmints_basispoints_h_
#include <libmints/basisset.h>
#include <libmints/vector3.h>
#include <libutil/ref.h>
namespace psi {

	class BasisPoints
	{
		protected:
			shared_ptr<BasisSet> basis_;
			double* points_;
			double* gradientsX_;
			double* gradientsY_;
			double* gradientsZ_;
			double* jacobiansXY_;
			double* jacobiansXZ_;
			double* jacobiansYZ_;
			double* jacobiansXX_;
			double* jacobiansYY_;
			double* jacobiansZZ_;
			double* laplacians_;
			bool do_points_;
			bool do_gradients_;
			bool do_jacobians_;
			bool do_laplacians_;
			bool have_points_;
			bool have_gradients_;
			bool have_jacobians_;
			bool have_laplacians_;
			void release();
			void allocate();
		public:
			BasisPoints(shared_ptr<BasisSet> b);
			~BasisPoints();
			void computePoints(Vector3 point);
			const double* getPoints() const { return points_; }
			const double* getGradientsX() const {return gradientsZ_; }
			const double* getGradientsY() const {return gradientsY_; }
			const double* getGradientsZ() const {return gradientsX_; }
			const double* getJacobiansXY() const {return jacobiansXY_; }
			const double* getJacobiansXZ() const {return jacobiansXZ_; }
			const double* getJacobiansYZ() const {return jacobiansYZ_; }
			const double* getJacobiansXX() const {return jacobiansXX_; }
			const double* getJacobiansYY() const {return jacobiansYY_; }
			const double* getJacobiansZZ() const {return jacobiansZZ_; } 
			const double* getLaplacians() const {return laplacians_; }
			void setToComputePoints(bool val) {do_points_ = val; allocate(); release(); }
			void setToComputeGradients(bool val) { do_gradients_ = val; allocate(); release(); }
			void setToComputeJacobians(bool val) { do_jacobians_ = val; allocate(); release(); }
			void setToComputeLaplacians(bool val) { do_laplacians_ = val; allocate(); release(); }
			int nbf() { return basis_->nbf();  }
			
	};
}
#endif
