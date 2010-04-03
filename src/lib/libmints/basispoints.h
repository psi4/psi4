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
			double* ao_points_; 
			double* ao_gradX_, *ao_gradY_, *ao_gradZ_;
			double* ao_hessXY_, *ao_hessXZ_, *ao_hessYZ_, *ao_hessXX_, *ao_hessYY_, *ao_hessZZ_;
			double *ao_laplac_; 
			
			double *prims;
			double *ang;
			double *Nam;
			int *a;
			int *b;
			int *c;

			double* points_;
			double* gradientsX_;
			double* gradientsY_;
			double* gradientsZ_;
			double* hessiansXY_;
			double* hessiansXZ_;
			double* hessiansYZ_;
			double* hessiansXX_;
			double* hessiansYY_;
			double* hessiansZZ_;
			double* laplacians_;
			bool do_points_;
			bool do_gradients_;
			bool do_hessians_;
			bool do_laplacians_;
			bool have_points_;
			bool have_gradients_;
			bool have_hessians_;
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
			const double* getHessiansXY() const {return hessiansXY_; }
			const double* getHessiansXZ() const {return hessiansXZ_; }
			const double* getHessiansYZ() const {return hessiansYZ_; }
			const double* getHessiansXX() const {return hessiansXX_; }
			const double* getHessiansYY() const {return hessiansYY_; }
			const double* getHessiansZZ() const {return hessiansZZ_; } 
			const double* getLaplacians() const {return laplacians_; }
			void setToComputePoints(bool val) {do_points_ = val; allocate(); release(); }
			void setToComputeGradients(bool val) { do_gradients_ = val; allocate(); release(); }
			void setToComputeHessians(bool val) { do_hessians_ = val; allocate(); release(); }
			void setToComputeLaplacians(bool val) { do_laplacians_ = val; allocate(); release(); }
			int nbf() { return basis_->nbf();  }
			
	};
}
#endif
