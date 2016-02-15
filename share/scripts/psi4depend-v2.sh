#!/bin/bash

# This is a psi4 dependencies installer based on Robert Harrison's 
# installer from madness. MSM

# You will need ~2GB of free disk to build and
# about 150MB of free disk to install.

#Specify the install and build directories 
#(they will be made if they do not exist)
INSTALLDIR=$HOME/temp-install2

#BUILDDIR should be the fastest possible disk available
#this will be destroyed at the end
BUILDDIR=/scratch/$USER/temp-build2

#If you need to force the build type (apple needs this)
BUILD=""
#BUILD="--build=i686-apple-darwin8"
#BUILD="--build=x86_64-apple-darwin10"

#No. of processes to use when building gcc
NPROC=`cat /proc/cpuinfo | grep processor | wc -l`
#NPROC=4

# Versions of everything
#PSI4 needs the following:
AUTOCONFVER=autoconf-2.68
AUTOMAKEVER=automake-1.11
M4VER=m4-1.4.13
GMPVER=gmp-4.3.1     #mpc depend
MPFRVER=mpfr-2.4.2   #mpc depend
MPCVER=mpc-0.8.2     #gcc depend
GCCVER=gcc-4.1.2
MPICH2VER=mpich2-1.2.1
PYTHONVER=2.6.6
LAPACKVER=lapack-3.2.2
BOOST=1.48.0 

######################################
# Should not need to modify below here
######################################

mkdir -p $INSTALLDIR
mkdir -p $BUILDDIR

#lets do some work
cd $BUILDDIR

export PATH=$INSTALLDIR/bin:$PATH
export LD_LIBRARY_PATH=$INSTALLDIR/lib:$INSTALLDIR/lib64:$LD_LIBRARY_PATH

# m4
echo 'WORKING ON M4'
curl -O http://vergil.chemistry.gatech.edu/download/$M4VER.tar.gz
tar xzf $M4VER.tar.gz
cd $M4VER
./configure --prefix=$INSTALLDIR
make  
make install  
cd ..
/bin/rm -rf $M4VER*

# autoconf
echo 'WORKING ON AUTOCONF'
curl -O http://vergil.chemistry.gatech.edu/download/$AUTOCONFVER.tar.gz
tar xzf $AUTOCONFVER.tar.gz
cd $AUTOCONFVER
./configure --prefix=$INSTALLDIR 
make  
make install  
cd ..
/bin/rm -rf $AUTOCONFVER*

# automake
echo 'WORKING ON AUTOMAKE'
curl -O http://vergil.chemistry.gatech.edu/download/$AUTOMAKEVER.tar.gz
tar xzf $AUTOMAKEVER.tar.gz
cd $AUTOMAKEVER
./configure --prefix=$INSTALLDIR  
make  
make install  
cd ..
/bin/rm -rf $AUTOMAKEVER*

# GMP - needed by gcc
echo 'WORKING ON GMP'
curl -O http://vergil.chemistry.gatech.edu/download/$GMPVER.tar.gz
tar xzf $GMPVER.tar.gz
cd $GMPVER
./configure --prefix=$INSTALLDIR $BUILD 
make -j $NPROC  
make install  
cd ..

# MPFR - needed by gcc
echo 'WORKING ON MPFR'
curl -O http://vergil.chemistry.gatech.edu/download/$MPFRVER.tar.gz
tar xzf $MPFRVER.tar.gz
cd $MPFRVER
./configure --prefix=$INSTALLDIR --with-gmp=$INSTALLDIR $BUILD  
make -j $NPROC  
make install  
cd ..

# MPC - needed by gcc
echo 'WORKING ON MPC'
curl -O http://vergil.chemistry.gatech.edu/download/$MPCVER.tar.gz
tar xzf $MPCVER.tar.gz
cd $MPCVER 
./configure --prefix=$INSTALLDIR --with-gmp=$INSTALLDIR --with-mpfr=$INSTALLDIR $BUILD  
make -j $NPROC  
make install  
cd ..

# gcc ... this will need about 2GB of empty space
echo 'WORKING ON GCC'
curl -O http://vergil.chemistry.gatech.edu/download/$GCCVER.tar.bz2
tar xjf $GCCVER.tar.bz2
mkdir build
cd build
../$GCCVER/configure --prefix=$INSTALLDIR --enable-languages=c,c++,fortran --with-gmp=$INSTALLDIR --with-mpfr=$INSTALLDIR --with-mpc=$INSTALLDIR --disable-checking $BUILD  
make -j $NPROC   
make install  
cd ..
/bin/rm -rf $GCCVER* $GMPVER* $MPFRVER* $MPCVER* build

# mpich2
echo 'WORKING ON MPICH2'
MPICH2NUM=`echo $MPICH2VER | sed -e 's/mpich2-//'`
curl -O http://vergil.chemistry.gatech.edu/download/$MPICH2VER.tar.gz
tar xzf $MPICH2VER.tar.gz
cd $MPICH2VER
./configure --prefix=$INSTALLDIR --with-pm=gforker,mpd --enable-fast --disable-f90 $BUILD  
make  
make install  
cd ..
/bin/rm -rf $MPICH2VER*

# Python
echo 'WORKING ON PYTHON'
curl -O http://vergil.chemistry.gatech.edu/download/Python-$PYTHONVER.tgz
tar xzf Python-$PYTHONVER.tgz
cd Python-$PYTHONVER
./configure --prefix=$INSTALLDIR  
make clean  
make install  
cd ..
/bin/rm -rf Python-$PYTHONVER

# Blas and Lapack
echo 'WORKING ON BLAS AND LAPACK'
curl -O http://vergil.chemistry.gatech.edu/download/$LAPACKVER.tgz
tar xzf $LAPACKVER.tgz
cd $LAPACKVER
cp INSTALL/make.inc.gfortran make.inc
make blaslib  
make lapacklib  
cp lapack_LINUX.a $INSTALLDIR/lib/liblapack.a
cp blas_LINUX.a $INSTALLDIR/lib/libblas.a
cd ..
/bin/rm -rf $LAPACKVER

# BOOST
echo 'WORKING ON BOOST'
curl -O http://vergil.chemistry.gatech.edu/download/boost_1_48_0.tar.gz
tar xzf boost_1_48_0.tar.gz
cd boost_1_48_0
sh bootstrap.sh
./b2 --prefix=$INSTALLDIR
cd ..
/bin/rm -rf build


#destroy build dir
rm -rf $BUILDDIR

echo
echo
echo '======================================'
echo 'Add the following to your *.rc path:' 
echo $INSTALLDIR/bin
echo 'Add the following to your *.rc LD_LIBRARY_PATH:'
echo $INSTALLDIR/lib:$INSTALLDIR/lib64
echo '======================================'



