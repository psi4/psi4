#!/usr/bin/env bash

if [ $# -lt "1" ]
then
  echo You must provide the GCC version number
  exit 1
fi 

mkdir boost_build

cd boost_build

tar -xjvf ../boost/boost_1_57_0.tar.bz2 >& /dev/null

cd boost_1_57_0

echo "using gcc : ${1} : g++-${1} ;" > ${HOME}/user-config.jam

./bootstrap.sh \
      --with-libraries="filesystem,python,regex,serialization,system,timer,chrono,thread" \
      --with-toolset=gcc \
      --prefix=/Users/andysim/travis/build/boost_comp

./b2 -j2 -q --toolset=gcc

cd ..

cd ..
