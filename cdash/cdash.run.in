#!/bin/sh -l
#
# Example PBS run file for testing a branch with different 
# compilers and settings. Should be reimplemented in Python.
#
###-------- PBS parameters ----------
#PBS -N cdash.run
#PBS -o cdash.out 
#PBS -e qsub.out
#PBS -j oe
#PBS -lnodes=2:ppn=8
#PBS -lwalltime=01:00:00
#PBS -lpmem=1000MB
###-------- end PBS parameters ----------

dashdir="@DASHBOARD_DIR@"

[ "x$PBS_SERVER" != "x" ] && site=$PBS_SERVER
[ "x$site" = "x" ] && site=`hostname -s`

case $site in 
    stallo*)
    module load cmake
    module load boost
    module load eigen
    module load openmpi
    module load intel-mkl
    module load valgrind
    toolchain_name=Intel
    ;;
    gardar*)
    module load intel-mkl
    module load valgrind
    ;;
esac

export MPIEXEC_PREFLAGS=""
export OMP_NUM_THREADS=1
export OMPI_MCA_btl=self,sm,tcp # disable Infiniband

if [ "x$PBS_JOBID" != "x" ]; then
    $dashdir/cdash.py --site=$site --branch=develop \
        --debug --memcheck --coverage 

    $dashdir/cdash.py --site=$site --branch=develop \
        --release --enable-mpi --enable-openmp 

    OMP_NUM_THREADS=8
    MPIEXEC_PREFLAGS="-pernode"
    $dashdir/cdash.py --site=$site --branch=develop \
        --release --enable-mpi --enable-openmp 
else
    $dashdir/cdash.py --site=$site --branch=develop
fi

# ------ END  ------
# vim:syntax=sh:filetype=sh
