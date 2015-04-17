source /theoryfs2/common/software/intel2015/bin/compilervars.sh intel64
#export TLIBC=/theoryfs2/ds/cdsgroup/psi4-compile/psi4cmake/psi4/glibc2.12rpm
#    -DLIBC_INTERJECT="-L${TLIBC}/usr/lib64 ${TLIBC}/lib64/libpthread.so.0 ${TLIBC}/lib64/libc.so.6" \

mkdir build
cd build
cmake \
    -DCMAKE_CXX_COMPILER=icpc \
    -DSTATIC_ONLY=ON \
    -DENABLE_XHOST=OFF \
    -DENABLE_TESTS=OFF \
    -DCMAKE_INSTALL_PREFIX=${PREFIX} \
    ${SRC_DIR}
make -j3  # -j${CPU_COUNT}
make install

