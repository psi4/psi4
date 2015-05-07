source /theoryfs2/common/software/intel2015/bin/compilervars.sh intel64
export TLIBC=/theoryfs2/ds/cdsgroup/psi4-compile/psi4cmake/psi4/glibc2.12rpm

mkdir build
cd build
cmake \
    -DCMAKE_C_COMPILER=icc \
    -DCMAKE_CXX_COMPILER=icpc \
    -DCMAKE_Fortran_COMPILER=ifort \
    -DPYTHON_INTERPRETER=${PYTHON} \
    -DLIBINT_OPT_AM=5 \
    -DENABLE_STATIC_LINKING=ON \
    -DENABLE_XHOST=OFF \
    -DCMAKE_BUILD_TYPE=release \
    -DLIBC_INTERJECT="-L${TLIBC}/usr/lib64 ${TLIBC}/lib64/libpthread.so.0 ${TLIBC}/lib64/libc.so.6" \
    -DCMAKE_INSTALL_PREFIX=${PREFIX} \
    ${SRC_DIR}
make -j3  # -j${CPU_COUNT}
make install

