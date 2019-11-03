#!/usr/bin/env bash

# The return code will capture an error from ANY of the functions in the pipe
set -o pipefail
export NINJA_STATUS="[Built edge %f of %t in %e sec] "
cmake --build . -- -j 2 -v -d stats 2>&1 | tee build.log | grep "Built"
RESULT=$?

if [ $RESULT -eq 0 ]; then
  echo build succeeded
else
  echo build failed
  cat build.log
  exit 1
fi
