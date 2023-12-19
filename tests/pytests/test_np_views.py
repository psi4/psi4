"""
This is a simple script that verifies several ways of accessing numpy arrays
and ensures that their memory is properly cleaned.
"""

import pytest
from addons import uusing

import numpy as np

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]

# If it's too small, something odd happens with the memory manager
mat_size = 10000


def snapshot_memory():
    import memory_profiler as mp

    return mp.memory_usage()[0] * 1048576


def check_leak(func, tol=1.e6):
    start = snapshot_memory()
    func()
    diff = abs(start - snapshot_memory())

    # A megabyte is excusable due to various GC funcs
    if diff > tol:
        raise MemoryError("Function did not correctly clean up: leaked %d bytes of memory!" % diff)
    else:
        print("Function %s: PASSED" % func.__name__)
        return True


def build_mat():
    mat = psi4.core.Matrix(mat_size, mat_size)
    return mat


def build_view_mat():
    mat = psi4.core.Matrix(mat_size, mat_size)
    view = mat.np
    return mat, view


def build_viewh_mat():
    mat = psi4.core.Matrix(mat_size, mat_size)
    view = mat.np
    return mat, view


def build_view_set_mat():
    mat = psi4.core.Matrix(mat_size, mat_size)
    view = mat.np
    view[:] = 5
    return mat, view


def build_arr_mat():
    mat = psi4.core.Matrix(mat_size, mat_size)
    view = np.asarray(mat)
    return mat, view


def build_copy_mat():
    mat = psi4.core.Matrix(mat_size, mat_size)
    view = np.array(mat)
    return mat, view


@uusing("memory_profiler")
def test_build_mat():
    assert check_leak(build_mat)


@uusing("memory_profiler")
def test_build_view_mat():
    assert check_leak(build_view_mat)


@uusing("memory_profiler")
def test_build_viewh_mat():
    assert check_leak(build_viewh_mat)


@uusing("memory_profiler")
def test_build_view_set_mat():
    assert check_leak(build_view_set_mat)


@uusing("memory_profiler")
def test_build_arr_mat():
    assert check_leak(build_arr_mat)


@uusing("memory_profiler")
def test_build_copy_mat():
    assert check_leak(build_copy_mat)


@uusing("memory_profiler")
def test_totals():
    start = snapshot_memory()

    check_leak(build_mat)
    check_leak(build_view_mat)
    check_leak(build_viewh_mat)
    check_leak(build_view_set_mat)
    check_leak(build_arr_mat)
    check_leak(build_copy_mat)

    # Double check totals
    diff = abs(start - snapshot_memory())
    if diff > 1.e6:
        raise MemoryError("\nA function leaked %d bytes of memory!" % diff)
    else:
        print("\nNo leaks detected!")
