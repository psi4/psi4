

import numpy as np
import math

__all__ = [
    "divisors",
    "geometry"
]

def divisors(n):
    """Generate the factors of an integer.

    Examples
    --------
    >>> for i in divisors(12): print i
    ...
    1
    2
    3
    4
    6
    12

    """
    yield 1
    for i in xrange(2, n/2+1):
        if n % i == 0:
            yield i
    yield n

class MathWarning(Warning):
    """ A warning for when something potentially unexpected is going to happen mathematically.
    """

class MatrixMultiplicationWarning(MathWarning):
    """ A warning specifically for when matrix multiplication may not be proceeding as expected.
    If you expect the possibility of element-wise matrix multiplication to be allowed, please catch this warning.
    """


import tensor; from tensor import *
__all__.extend(tensor.__all__)

import vector; from vector import *
__all__.extend(vector.__all__)

import matrix; from matrix import *
__all__.extend(matrix.__all__)

import geometry; from geometry import *
__all__.extend(geometry.__all__)

import einsum; from einsum import *
__all__.extend(einsum.__all__)

import einsum_indices; from einsum_indices import *
__all__.extend(einsum_indices.__all__)

