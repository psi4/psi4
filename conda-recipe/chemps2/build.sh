source /theoryfs2/common/software/intel2015/bin/compilervars.sh intel64
export TLIBC=/theoryfs2/ds/cdsgroup/psi4-compile/psi4cmake/psi4/glibc2.12rpm
REDIST=/theoryfs2/common/software/intel2015/composer_xe_2015.0.090/compiler/lib/intel64/libiomp5.so

mkdir build
cd build
mkdir lib
cp ${REDIST} ${PREFIX}/lib
cmake \
    -DCMAKE_CXX_COMPILER=icpc \
    -DMKL=ON \
    -DBUILD_DOCUMENTATION=OFF \
    -DENABLE_TESTS=OFF \
    -DENABLE_GENERIC=ON \
    -DENABLE_XHOST=OFF \
    -DLIBC_INTERJECT="-L${TLIBC}/usr/lib64 ${TLIBC}/lib64/libpthread.so.0 ${TLIBC}/lib64/libc.so.6" \
    -DCMAKE_INSTALL_PREFIX=${PREFIX} \
    -DCMAKE_INSTALL_LIBDIR=lib \
    ${SRC_DIR}
make -j3  # -j${CPU_COUNT}
make install

