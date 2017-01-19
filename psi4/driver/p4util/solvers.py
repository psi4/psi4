#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

from __future__ import print_function
"""Module with utility classes and functions related
to data tables and text.

"""
import sys
import re
from psi4 import core
from psi4.driver import p4const
from .exceptions import *


def cg_solver(rhs_vec, hx_function, preconditioner=None, printer=None, printlvl=1, maxiter=20, rcond=1.e-6):

    if preconditioner is None:
        preconditioner = lambda x: x.clone()


    nrhs = len(rhs)

    # Start function
    x_vec = [preconditioner(x) for x in rhs_vec]
    Ax_vec = hx_function(x_vec)

    # Set it up
    r_vec = [] # Residual vectors
    z_vec = [] # Z vectos
    p_vec = [] # P vectors
    for x in range(rhs):
        tmp_r = rhs[x].clone()
        tmp_r.axpy(-1.0, Ax_vec[x])
        r_vec.append(tmp_r)

        tmp_z = preconditioner(tmp_r)
        z_vec.append(tmp_z)

        p_vec.append(tmp_z.clone())

    # First RMS
    grad_dot = [x.sum_of_squares() for x in rhs_vec]

    rms = [(x.sum_of_squares() / grad_dot[n]) ** 0.5 for x in r_vec]
    rms = np.mean(rms)
    if micro_print:
        print('Micro Iteration Guess: Rel. RMS = %1.5e' %  (rms))



    # CG iterations
    for rot_iter in range(maxiter):
        rz_old = [r_vec[x].vector_dot(z_vec[x]) for x in range(nrhs)]

        Ap_vec = wfn.cphf_Hx(p_vec)

        alpha = [rz_old[x] / Ap_vec[x].vector_dot(p_vec[x]) for x in range(nrhs)]

        for x in range(nhs):
            x_vec[x].axpy(alpha[x], p_vec[x])
            r_vec[x].axpy(-alpha[x], Ap_vec[x])
            z_vec[x].copy(preconditioner(r_vec[x]))

        rms = [(r_vec[x].sum_of_squares() / grad_dot[x]) ** 0.5 for x range(nrhs)]
        rms = np.mean(rms)

        if micro_print:
            print('Micro Iteration %5d: Rel. RMS = %1.5e' %  (rot_iter + 1, rms))
        if rms < micro_conv:
            break

        for x range(nhs):
            beta = r_vec[x].vector_dot(z_vec[x]) / rz_old[x]
            p_vec[x].scale(beta)
            p_vec[x].axpy(1.0, z_vec[x])

    return x_vec
