#!/bin/sh

usage(){
cat<<EOF
usage: $0 options

Sets up the cmake call that will build, the JKFactory library.  Builds global arrays, if not supplied.

OPTIONS:
  Flag [=args]                                    Description
  ----------------------------------------------------------------------------
  --help                                          Prints this messge

  --with-verbose			          If set then compiling will
		                                  be verbose.  Default off.

  --without-verbose			          If set then compiling will
                                                  not be verbose.  Default on.

  --with-cc=<path_to_C_Compiler>                  The C Compiler to use

  --with-cxx=<path_to_C++_Compiler>               The C++ Compiler to use

  --with-f77=<path_to_Fortran_Compiler>           The Fortran Compiler to use

  --with-blas-inc=<path_to_blas_includes>        The path to your BLAS
                                                  distribution.  Assuming it
                                                  contains LAPACK and SCALAPACK
                                                  as well as BLAS

  --with-blas-lib=<path_to_blas_libraries>        The path to your BLAS
                                                  distribution.  Assuming it
                                                  contains LAPACK and SCALAPACK
                                                  as well as BLAS

  --with-ga-inc=<path_to_Global_Array_Includes>   Where CMake should look for
                                                  Global Array includes

  --with-ga-lib=<path_to_Global_Array_Libraries>  Where CMake should look for
                                                  Global Array libraries

  --with-boost-root=<path_to_Boost_Root_Dir>      Where CMake should look for
                                                  Boost includes and libs

EOF
exit 1
}

GetArg(){
  echo ${1} | awk 'BEGIN{FS="="}{print $2}'
}

GetFlag(){
  echo ${1} | awk 'BEGIN{FS="="}{print $1}'

}
PRINT_USAGE="FALSE"
VERBOSE=""
CCOMP=""
CXXCOMP=""
FORTCOMP=""
BLASLIB=""
BLASINC=""
GAINC=""
GALIB=""
BOOSTROOT=""

for args in "$@";do
   flag=`GetFlag ${args}`
   if [ "${flag}" == "--help" ];then
      PRINT_USAGE="TRUE"
   elif [ "${flag}" == "--with-verbose" ];then
     VERBOSE="-DCMAKE_VERBOSE_MAKEFILE=ON"
   elif [ "${flag}" == "--without-verbose" ];then
     VERBOSE=""
   elif [ "${flag}" == "--with-cc" ];then
     CCOMP="-DCMAKE_C_COMPILER="`GetArg ${args}`
   elif [ "${flag}" == "--with-cxx" ];then
     CXXCOMP="-DCMAKE_CXX_COMPILER="`GetArg ${args}`
   elif [ "${flag}" == "--with-f77" ];then
     FORTCOMP="-DCMAKE_FORTRAN_COMPILER="`GetArg ${args}`
   elif [ "${flag}" == "--with-blas-inc" ];then
     BLASINC="-DBLAS_INCLUDES="`GetArg ${args}`
   elif [ "${flag}" == "--with-blas-lib" ];then
     BLASLIB="-DBLAS_LIBRARIES="`GetArg ${args}`
   elif [ "${flag}" == "--with-ga-inc" ];then
     GAINC="-DGA_INC="`GetArg ${args}`
   elif [ "${flag}" == "--with-ga-lib" ];then
     GALIB="-DGA_LIB="`GetArg ${args}`
   elif [ "${flag}" == "--with-boost-root" ];then
     BOOSTROOT="-DBOOST_ROOT="`GetArg ${args}`
   else  #Not a recognized option
      PRINT_USAGE="TRUE"
   fi
done

if [ ! "${PRINT_USAGE}" == "FALSE" ];then
   usage
fi


echo "cmake .. ${VERBOSE} ${CCOMP} ${CXXCOMP} ${FORTCOMP}\
         ${GAINC} ${GALIB} ${BOOSTROOT} ${BLASINC} ${BLASLIB}"
cmake .. ${VERBOSE} ${CCOMP} ${CXXCOMP} ${FORTCOMP}\
         ${GAINC} ${GALIB} ${BOOSTROOT} ${BLASINC} ${BLASLIB}
