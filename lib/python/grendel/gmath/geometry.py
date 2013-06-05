"""
Module containing functions for calculating and transforming various aspects of Cartesian coordinate systems and
Cartesian representations.  May be refactored to be included in another class sometime.
"""
import math
import numpy as np
from grendel.gmath.vector import Vector
from matrix import Matrix

__all__ = [
    "rotate_point_about_axis",
    "angle_between_vectors"
]

def rotate_point_about_axis(point, axis, angle):
    """ Returns the Cartesian vector `point` rotated about `axis` by `angle`

    Source:  http://en.wikipedia.org/wiki/Rotation_matrix
    """
    u = axis.normalized()
    cos = math.cos(angle)
    sin = math.sin(angle)
    R = Matrix(
        [
            [ cos + u.x**2*(1-cos), u.x*u.y*(1-cos) - u.z*sin, u.x*u.z*(1-cos) + u.y*sin ],
            [ u.y*u.x*(1-cos) + u.z*sin, cos + u.y**2*(1-cos), u.y*u.z*(1-cos) - u.x*sin ],
            [ u.z*u.x*(1-cos) - u.y*sin, u.z*u.y*(1-cos) + u.x*sin, cos + u.z**2*(1-cos) ]
        ]
    )
    return R * point

def angle_between_vectors(a, b):
    """ Returns the angle between the two vectors, in radians
    """
    if a.is_zero() or b.is_zero():
        raise ZeroDivisionError("can't get angle between vectors {0} and {1} when one of them is almost 0".format(a, b))
    am = a.magnitude()
    bm = b.magnitude()
    operand = a.dot(b) / (am * bm)
    # Chop values slightly larger than 1 or slightly smaller than -1
    if abs(operand) - 1.0 > 1e-3:
        raise ValueError("math domain error for acos function: " + str(operand) + " is not between -1 and 1.")
    if operand > 1.0: operand = 1.0
    elif operand < -1.0: operand = -1.0
    return math.acos(operand)

def safe_math_func(func, domain_bottoms, domain_tops, cutoff=1e-10):
    # TODO @quick implement this
    raise NotImplementedError

def safe_acos(ang):
    try:
        return math.acos(ang)
    except ValueError as e:
        # Not sure what the actual cutoff is before a math domain error is raised,
        # but this should be sufficient. Probably could go smaller with this cutoff
        # without causing problems, but there isn't really a good reason to do so.
        if ang - 1 < 1e-10:
            return math.pi / 2.0
        elif ang + 1 < 1e-10:
            return 3.0 * math.pi / 2.0
        else:
            raise e

def safe_asin(ang):
    try:
        return math.asin(ang)
    except ValueError as e:
        # Not sure what the actual cutoff is before a math domain error is raised,
        # but this should be sufficient. Probably could go smaller with this cutoff
        # without causing problems, but there isn't really a good reason to do so.
        if ang - 1 < 1e-10:
            return math.pi / 2.0
        elif ang + 1 < 1e-10:
            return 3.0 * math.pi / 2.0
        else:
            raise e

