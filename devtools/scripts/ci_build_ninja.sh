#!/usr/bin/env bash

# The return code will capture an error from ANY of the functions in the pipe
set -o pipefail
cmake --build . 2>&1 | tee build.log | grep "Building"
RESULT=$?

if [ $RESULT -eq 0 ]; then
  echo build succeeded
else
  echo build failed
  cat build.log
  exit 1
fi
