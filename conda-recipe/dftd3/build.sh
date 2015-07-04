# load Intel compilers and mkl
source /theoryfs2/common/software/intel2015/bin/compilervars.sh intel64

make
# pseudo "make install"
mkdir ${PREFIX}/bin
cp dftd3 ${PREFIX}/bin

