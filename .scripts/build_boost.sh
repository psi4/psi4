#!/usr/bin/env bash

if [ $# -lt "3" ]
then
  echo You must provide the compiler name, version number and c++ compile command
  exit 1
fi 

echo "using ${1} : ${2} : ${3} : ;" > ${HOME}/user-config.jam

mkdir boost_build

cd boost_build
tar -xjvf ../external/boost/boost-src/boost_1_57_0.tar.bz2 >& /dev/null
cd boost_1_57_0
./bootstrap.sh \
      --with-libraries="filesystem,python,regex,serialization,system,timer,chrono,thread" \
      --with-toolset=${1} \
      --prefix=${TRAVIS_BUILD_DIR}/boost_install
./b2 -j2 -q --toolset=${1} --variant=debug install >& /dev/null
cd ..
cd ..
