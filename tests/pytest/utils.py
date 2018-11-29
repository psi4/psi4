import sys
import pprint
pp = pprint.PrettyPrinter(width=120)

import psi4
from psi4.driver import qcdb


a2a = 0.52917721067 / 0.52917720859

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
compare_matrices = true_false_decorator(psi4.compare_matrices)
compare_arrays = true_false_decorator(qcdb.compare_arrays)
compare_dicts = true_false_decorator(qcdb.compare_dicts)
compare_molrecs = true_false_decorator(qcdb.compare_molrecs)


def tnm():
    """Returns the name of the calling function, usually name of test case."""

    return sys._getframe().f_back.f_code.co_name
