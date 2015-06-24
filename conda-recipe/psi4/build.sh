source /theoryfs2/common/software/intel2015/bin/compilervars.sh intel64
export TLIBC=/theoryfs2/ds/cdsgroup/psi4-compile/nightly/glibc2.12

mkdir build
cd build
cmake \
    -DCMAKE_C_COMPILER=icc \
    -DCMAKE_CXX_COMPILER=icpc \
    -DCMAKE_Fortran_COMPILER=ifort \
    -DPYTHON_INTERPRETER=${PYTHON} \
    -DLIBINT_OPT_AM=6 \
    -DENABLE_STATIC_LINKING=ON \
    -DENABLE_XHOST=OFF \
    -DENABLE_PLUGINS=OFF \
    -DENABLE_CHEMPS2=OFF \
    -DSPHINX_ROOT=/theoryfs2/ds/cdsgroup/psi4-install/miniconda/envs/boostenv \
    -DCMAKE_BUILD_TYPE=release \
    -DLIBC_INTERJECT="-L${TLIBC}/usr/lib64 ${TLIBC}/lib64/libpthread.so.0 ${TLIBC}/lib64/libc.so.6" \
    -DBUILDNAME=LAB-intel15.0.3-mkl-release-conda \
    -DSITE=gatech-conda \
    -DCMAKE_INSTALL_PREFIX=${PREFIX} \
    ${SRC_DIR}
make -j3  # -j${CPU_COUNT}  # get incomplete build at full throttle
make sphinxman
make ghfeed
make install

