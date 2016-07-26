#!/usr/bin/env bash

if [ $# -lt "3" ]
then
  echo You must provide the compiler name, version number and c++ compile command
  exit 1
fi 


echo "using ${1} : ${2} : ${3} : ;" > ${HOME}/user-config.jam
