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
