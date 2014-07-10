#wget http://hpc.pnl.gov/globalarrays/download/ga-5-2.tgz
export ARMCI_OPENIB_DEVICE="mlx4_0"
cd $1
./configure --with-tcgmsg \
            --prefix="/home1/02036/xingliu/GTFock/GAlib" \
            --with-blas="-I/opt/apps/intel/11.1/mkl/include -L/opt/apps/intel/11.1/mkl/lib/em64t/ -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5" \
            --with-mpi="-I/opt/apps/intel11_1/mvapich2/1.6/include -L/opt/apps/intel11_1/mvapich2/1.6/lib/ -lmpich" \
            --with-lapack="-I/opt/apps/intel/11.1/mkl/include -L/opt/apps/intel/11.1/mkl/lib/em64t/ -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5"  \
            --with-scalapck="-I/opt/apps/intel/11.1/mkl/include -L/opt/apps/intel/11.1/mkl/lib/em64t/ -lmkl_scalapack_lp64" \
            --with-openib="-I/opt/ofed/include/ -L/opt/ofed/lib64/ -libumad -libverbs" 
	
make CC= /opt/apps/intel/11.1/bin/intel64/icc FC= /opt/apps/intel/11.1/bin/intel64/ifort MPICC=/opt/apps/intel11_1/mvapich2/1.6/bin/mpicc
make install
