#!/bin/bash
# 
# modeled on https://github.com/pybind/pybind11/blob/master/tools/check-style.sh

# Script to check include/test code for common code style errors.
# 
# This script currently checks for
#
# 1. use of tabs instead of spaces
# 2. trailing spaces
# 
# Invoke as: cmake/check-style.sh
#

errors=0
IFS=$'\n'
found=
# The mt=41 sets a red background for matched tabs:
exec 3< <(GREP_COLORS='mt=41' grep -rin --include \*.h --include \*.cc $'\t' psi4/ --color=always)
while read -u 3 f; do
    if [ -z "$found" ]; then
        echo -e '\e[31m\e[01mError: found tabs instead of spaces in the following files:\e[0m'
        found=1
        errors=1
    fi

    echo "    $f"
done

found=
# The mt=41 sets a red background for matched trailing spaces
exec 3< <(GREP_COLORS='mt=41' grep -rin --include \*.h --include \*.cc '\s\+$' psi4/ --color=always)
while read -u 3 f; do
    if [ -z "$found" ]; then
        echo -e '\e[31m\e[01mError: found trailing spaces in the following files:\e[0m'
        found=1
        errors=1
    fi

    echo "    $f"
done

