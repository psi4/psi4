#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2022 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import qcdb


def true_false_decorator(compare_fn, *args, **kwargs):
    """Turns `compare_fn` that returns `None` on success and raises
    `qcdb.TestComparisonError` on failure into a function that returns
    True/False, suitable for assertions in pytest.

    """

    def true_false_wrapper(*args, **kwargs):
        try:
            compare_fn(*args, **kwargs)
        except qcdb.TestComparisonError as err:
            return False
        else:
            return True

    return true_false_wrapper


compare_values = true_false_decorator(qcdb.compare_values)
compare_strings = true_false_decorator(qcdb.compare_strings)
compare_integers = true_false_decorator(qcdb.compare_integers)
compare_matrices = true_false_decorator(qcdb.compare_matrices)
compare_arrays = true_false_decorator(qcdb.compare_arrays)
compare_dicts = true_false_decorator(qcdb.compare_dicts)
compare_molrecs = true_false_decorator(qcdb.compare_molrecs)
