#!/bin/sh

OUTPUT_COMPLETED="COMPLETED SUCCESSFULLY"
OUTPUT_MATCH="DOES NOT MATCH"

RED="\033[33;31m"
GREEN="\033[33;32m"
NORMAL="\033[0m"

print_success()
{
	echo -e "${GREEN}SUCCESS: ${TEST}${NORMAL}"
}

print_failure()
{
	echo -e "${RED}FAILURE: ${TEST}${NORMAL}"
}

for TEST in *.in; do
	TEST=`basename ${TEST} .in`
	${EFPMD} ${TEST}.in > ${TEST}.out

	if grep -q "${OUTPUT_COMPLETED}" ${TEST}.out; then
		if grep -q "${OUTPUT_MATCH}" ${TEST}.out; then
			print_failure
		else
			print_success
		fi
	else
		print_failure
	fi
done
