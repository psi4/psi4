#!/usr/bin/env bash

set -o pipefail
ctest -j2 -L quicktests 2>&1 | tee run.log
RESULT=$?

# Now we scrape the run log, looking for any failures
while read p; do
  set -o pipefail
  echo $p | grep "(Failed)"
  if [ $? -eq 0 ]; then
    badtest=`echo $p | awk '{print $3}'`
    echo $badtest failed.  Here is the output:
    cat tests/${badtest}/output.dat
  fi
done <run.log

exit $RESULT
