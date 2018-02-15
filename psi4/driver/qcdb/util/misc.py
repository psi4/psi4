import re
import sys
import math

import numpy as np

from ..physconst import psi_bohr2angstroms

if sys.version_info >= (3, 0):
    basestring = str


def distance_matrix(a, b):
    """Euclidean distance matrix between rows of arrays `a` and `b`. Equivalent to
    `scipy.spatial.distance.cdist(a, b, 'euclidean')`. Returns a.shape[0] x b.shape[0] array.

    """
    assert a.shape[1] == b.shape[1], """Inner dimensions do not match"""
    distm = np.zeros([a.shape[0], b.shape[0]])
    for i in range(a.shape[0]):
        for j in range(b.shape[0]):
            distm[i, j] = np.linalg.norm(a[i] - b[j])
    return distm


def update_with_error(a, b, path=None):
    """Merges `b` into `a` like dict.update; however, raises KeyError if values of a
    key shared by `a` and `b` conflict.

    Adapted from: https://stackoverflow.com/a/7205107

    """
    if path is None:
        path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                update_with_error(a[key], b[key], path + [str(key)])
            elif a[key] == b[key]:
                pass  # same leaf value
            elif a[key] is None:
                a[key] = b[key]
            elif (isinstance(a[key], (list, tuple)) and
                  not isinstance(a[key], basestring) and
                  isinstance(b[key], (list, tuple)) and
                  not isinstance(b[key], basestring) and
                  len(a[key]) == len(b[key]) and
                  all((av is None or av == bv) for av, bv in zip(a[key], b[key]))):  # yapf: disable
                a[key] = b[key]
            else:
                raise KeyError('Conflict at {}: {} vs. {}'.format('.'.join(path + [str(key)]), a[key], b[key]))
        else:
            a[key] = b[key]
    return a


def standardize_efp_angles_units(units, geom_hints):
    """Applies to the pre-validated xyzabc or points hints in `geom_hints`
    the libefp default (1) units of [a0] and (2) radian angle range of
    (-pi, pi]. The latter is handy since this is how libefp returns hints

    """

    def radrge(radang):
        """Adjust `radang` by 2pi into (-pi, pi] range."""
        if radang > math.pi:
            return radang - 2 * math.pi
        elif radang <= -math.pi:
            return radang + 2 * math.pi
        else:
            return radang

    if units == 'Angstrom':
        iutau = 1. / psi_bohr2angstroms
    else:
        iutau = 1.

    hints = []
    for hint in geom_hints:
        if len(hint) == 6:
            x, y, z = [i * iutau for i in hint[:3]]
            a, b, c = [radrge(i) for i in hint[3:]]
            hints.append([x, y, z, a, b, c])
        if len(hint) == 9:
            points = [i * iutau for i in hint]
            hints.append(points)

    return hints


def filter_comments(string):
    """Remove from `string` any Python-style comments ('#' to end of line)."""

    comment = re.compile(r'(^|[^\\])#.*')
    string = re.sub(comment, '', string)
    return string


def unnp(dicary):
    """Return `dicary` with any ndarray values replaced by lists."""

    ndicary = {}
    for k, v in dicary.items():
        try:
            v.shape
        except AttributeError:
            ndicary[k] = v
        else:
            ndicary[k] = v.tolist()
    return ndicary
