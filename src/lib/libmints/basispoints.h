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
			double* ao_points_; 
			double* ao_gradX_, *ao_gradY_, *ao_gradZ_;
			double* ao_hessXY_, *ao_hessXZ_, *ao_hessYZ_, *ao_hessXX_, *ao_hessYY_, *ao_hessZZ_;
			double* ao_laplac_; 
			
			double *prims;
			double *ang;
			double *Nam;
			int *a;
			int *b;
			int *c;

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
		public:
			BasisPoints(boost::shared_ptr<BasisSet> b, int block_size);
			~BasisPoints();
			void computePoints(boost::shared_ptr<GridBlock> block);
			double** getPoints() { return points_; }
			double** getGradientsX() {return gradientsZ_; }
			double** getGradientsY() {return gradientsY_; }
			double** getGradientsZ() {return gradientsX_; }
			double** getHessiansXY() {return hessiansXY_; }
			double** getHessiansXZ() {return hessiansXZ_; }
			double** getHessiansYZ() {return hessiansYZ_; }
			double** getHessiansXX() {return hessiansXX_; }
			double** getHessiansYY() {return hessiansYY_; }
			double** getHessiansZZ() {return hessiansZZ_; } 
			double** getLaplacians() {return laplacians_; }
			void setToComputePoints(bool val) {do_points_ = val; allocate(); release(); }
			void setToComputeGradients(bool val) { do_gradients_ = val; allocate(); release(); }
			void setToComputeHessians(bool val) { do_hessians_ = val; allocate(); release(); }
			void setToComputeLaplacians(bool val) { do_laplacians_ = val; allocate(); release(); }
			int nbf();
            int getBlockSize() { return block_size_;}  
			int getTrueSize() {return true_size_; }
            void setTrueSize(int size) {true_size_ = size; }
	};
}
#endif
