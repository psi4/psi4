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

import numpy as np
import psi4
import os
from MCSCF import *

# Shoody rel imports
import sys
import inspect
path_dir = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe()))[0], '../../')))
sys.path.append(path_dir)

import p4util
from p4util.exceptions import *
np.set_printoptions(precision=5, linewidth=200, threshold=2000, suppress=True)

def OCHx(mcscf_obj, ciwfn, C0, vector, mcscf):

    # Build OPDM, < 0 | \gamma | c >
    c_opdm = ciwfn.opdm(C0, vector, 0, 0)[2]
    nirrep = c_opdm.nirrep()
    for h in range(nirrep):
        c_opdm.nph[h][:] += c_opdm.nph[h].T.copy()

    # Build TPDM and symmetrize, < 0 | \Gamma | c >
    c_tpdm = ciwfn.tpdm(C0, vector, 0, 0)[3]
    tmp_tpdm = np.array(c_tpdm)
    tmp_tpdm += np.einsum('pqrs->qprs', c_tpdm)
    tmp_tpdm += np.einsum('pqrs->rspq', tmp_tpdm).copy()
    tmp_tpdm *= 0.5

    c_tpdm.np[:] = tmp_tpdm
    c_tpdm.np[:] += tmp_tpdm.T

    # Get dim information
    dim_docc = ciwfn.get_dimension('DRC')
    dim_oa = ciwfn.get_dimension('OA')
    dim_av = ciwfn.get_dimension('AV')
    dim_rot = ciwfn.get_dimension('ROT')

    # Build required components
    IF = mcscf_obj.current_IFock()
    rot_AF = mcscf_obj.compute_AFock(c_opdm)
    rot_Q = mcscf_obj.compute_Q(c_tpdm)

    # Setup output matrices
    ci_overlap = C0.vdot(vector, 0, 0)
    matF = psi4.Matrix('OCFx', dim_rot, dim_rot)
    matgrad = psi4.Matrix('OCFx', dim_rot, dim_rot)
    ret = psi4.Matrix('OCFx', dim_oa, dim_av)

    # Build outputs
    for h in range(nirrep):

        # Build generalized Fock
        if dim_docc[h]:
            matF.nph[h][:dim_docc[h], :] += IF.nph[h][:, :dim_docc[0]].T
            matF.nph[h][:dim_docc[h], :] *= ci_overlap
            matF.nph[h][:dim_docc[h], :] += rot_AF.nph[h][:, :dim_docc[0]].T
            matF.nph[h][:dim_docc[h], :] *= 2

        asl = slice(dim_docc[h], dim_oa[h])
        if dim_oa[h] - dim_docc[h]:
            matF.nph[h][asl, :] = np.dot(c_opdm, IF.nph[h][:, asl].T)
            matF.nph[h][asl, :] += rot_Q.nph[h]

        # Build gradient
        matgrad.nph[h][:] = matF.nph[h].T
        matgrad.nph[h][:] -= matF.nph[h]
        matgrad.nph[h][:] *= -1

        # Slice out active rotations
        ret.nph[h][:] = matgrad.nph[h][:dim_oa[h], dim_docc[h]:]
        ret.nph[h][:] *= 2

    mcscf_obj.zero_redundant(ret)
    return ret


def CCHx(ciwfn, vector, output, ci_energy):

    # 2.0 * (< r | H | c > - c * e_ci)
    ciwfn.sigma(vector, output, 0, 0)
    output.axpy(-ci_energy, vector, 0, 0)
    output.scale(2.0, 0)
    return np.array(output)


def COHx(ciwfn, C0, k, output, nact):

    # Build temporaries
    ncitri = nact * (nact + 1) / 2
    ncitri2 = ncitri * (ncitri + 1) / 2
    onel_k = psi4.Vector(ncitri)
    twoel_k = psi4.Vector(ncitri2)

    # Build k as a matrix, not needed? 
    matK = psi4.Matrix(k.shape[0], k.shape[1])
    matK_arr = np.asarray(matK)
    matK_arr.flat[:] = k

    # Build rotated ERI's
    ciwfn.rotate_mcscf_integrals(matK, onel_k, twoel_k)

    # Sigma computation
    ciwfn.sigma(C0, output, 0, 0, onel_k, twoel_k)
    output.scale(2.0, 0)


def compute_vector(proj_x, mcscf_obj, mcscf, ciwfn, orb_grad, ci_grad, C0, dvecs1, dvecs2):

    # Build up a lot of scratch space because im bad
    scratch1 = ciwfn.new_civector(1, 0, False, True)
    scratch2 = ciwfn.new_civector(1, 0, False, True)
    ci_Ap = ciwfn.new_civector(1, 0, False, True)
    C0vec = ciwfn.new_civector(1, 0, False, True)
    ci_grad_vec = ciwfn.new_civector(1, 0, False, True)

    # Set data
    ci_grad_vec.np[:] = ci_grad
    ndet = C0.shape[0]
    ci_p = ciwfn.new_civector(1, 0, False, True)
    ci_p.np[:] = proj_x[:ndet]
    C0vec.np[:] = C0

    # Project P right
    orig_ci_scale_cp = np.vdot(C0, ci_p.np)
    ci_p.axpy(-orig_ci_scale_cp, C0vec, 0, 0)
    ci_p.symnormalize(1.0, 0)
    ci_scale_cp = np.vdot(C0, ci_p.np)

    # Make matrices
    orb_p = psi4.Matrix.from_array(proj_x[ndet:].reshape(orb_grad.shape))
    orb_grad = psi4.Matrix.from_array(orb_grad)

    # Grab orbital info 
    nact = sum(ciwfn.get_dimension('ACT').to_tuple())
    ndocc = sum(ciwfn.get_dimension('DRC').to_tuple())
    nvir = sum(ciwfn.get_dimension('VIR').to_tuple())

    # Build ci_Ap
    CCHx(ciwfn, ci_p, ci_Ap, mcscf_obj.current_ci_energy())
    ci_Ap.axpy(2.0 * ci_scale_cp, C0vec, 0, 0)
    COHx(ciwfn, C0vec, orb_p, scratch1, nact)
    ci_Ap.axpy(1.0, scratch1, 0, 0)
    scale = ci_p.vdot(ci_grad_vec, 0, 0)
    scale += 2 * orb_grad.vector_dot(orb_p)
    ci_Ap.axpy(-scale, C0vec, 0, 0)
    ci_Ap.axpy(-ci_scale_cp, ci_grad_vec, 0, 0)

    # Build orb_Ap
    orb_Ap = OCHx(mcscf_obj, ciwfn, C0vec, ci_p, mcscf)
    orb_Ap.add(mcscf_obj.compute_Hk(orb_p))
    orb_Ap.axpy(-2.0 * ci_scale_cp, orb_grad)
    orb_Ap = np.array(orb_Ap)
    orb_Ap[ndocc:, :nact] = 0

    # Project P left
    proj_scale_cp = ci_Ap.vdot(C0vec, 0, 0)
    ci_Ap.axpy(-proj_scale_cp, C0vec, 0, 0)

    Ap = np.hstack((ci_Ap.np, orb_Ap.ravel()))

    return Ap


def qc_iteration(dvec, ci_grad, ciwfn, mcscf_obj):

    # Build temps
    C0 = np.array(dvec)
    ndet = C0.shape[0]

    ci_precon = ciwfn.new_civector(1, 0, False, True)
    dvecs1 = ciwfn.new_civector(1, 0, False, True)
    dvecs1_buff = np.asarray(dvecs1)
    dvecs2 = ciwfn.new_civector(1, 0, False, True)
    dvecs2_buff = np.asarray(dvecs2)
    dvecs1_buff[:] = C0
    dvecs2_buff[:] = C0

    # Density matrices
    opdm = np.array(ciwfn.get_opdm(-1, -1, 'SUM', False))
    tpdm = np.array(ciwfn.get_tpdm('SUM', True))

    # Orbital Info
    docc = sum(ciwfn.get_dimension('DOCC').to_tuple())
    act = sum(ciwfn.get_dimension('ACT').to_tuple())
    vir = sum(ciwfn.get_dimension('VIR').to_tuple())
    noa = docc + act
    nav = act + vir
    mcscf = None

    # Grab temps
    ci_grad = ci_grad
    CI_elec_energy = mcscf_obj.current_total_energy()
    orb_grad = np.array(mcscf_obj.gradient())
    diag_hess = np.array(mcscf_obj.H_approx_diag())
    norm_ci_grad = ci_grad.norm(0)
    norm_orb_grad = np.linalg.norm(orb_grad)
    dvecs1_buff[:] = C0
    dvecs1.symnormalize(1.0, 0)
    C0[:] = np.asarray(dvecs1_buff)

    # Gradients
    grad = -np.hstack((ci_grad.np, orb_grad.ravel()))
    RHS = grad.ravel()

    # Preconditioners

    # CI preconditioner
    # 2 * (Hd - E + C0 * C0) - 2.0 * C0 * ci_grad
    ci_precon = ciwfn.Hd_vector(0)
    ci_precon.vector_multiply(1.0, dvecs1, dvecs2, 0, 0, 0)
    ci_precon.shift(-mcscf_obj.current_total_energy(), 0)
    ci_precon.scale(2.0, 0)
    ci_precon.vector_multiply(-2.0, dvecs1, ci_grad, 0, 0, 0)

    ci_precon = np.array(ci_precon)
    ci_grad = np.array(ci_grad)
    
    #ci_precon = 2 * (Hd - mcscf_obj.current_total_energy() + C0 * C0)
#    ci_precon -= 2.0 * C0 * ci_grad
    orb_precon = diag_hess.copy()
    precon = np.hstack((ci_precon, orb_precon.ravel()))

    # Initial X
    ci_x = -ci_grad / ci_precon
    orb_x = -orb_grad / orb_precon
    x = RHS / precon
    x[:ndet][np.abs(x[:ndet]) < 1e-15] = 0
    ci_x[:ndet][np.abs(x[:ndet]) < 1e-15] = 0

    # Ax product
    Ax = compute_vector(x, mcscf_obj, mcscf, ciwfn, orb_grad, ci_grad, C0, dvecs1, dvecs2)

    # CG update
    r = RHS - Ax
    z = r / precon
    p = z.copy()
    r[:ndet][np.abs(x[:ndet]) < 1e-14] = 0

    rel_tol = np.vdot(r, r)
    print rel_tol
    raise Exception("")

    # CG loops
    for rot_iter in range(12):
        rz_old = np.vdot(r, z)
        p[:ndet][np.abs(p[:ndet]) < 1e-15] = 0
        Ap = compute_vector(p, mcscf_obj, mcscf, ciwfn, orb_grad, ci_grad, C0, dvecs1, dvecs2)
        alpha = rz_old / np.vdot(Ap, p)
        if alpha > 15:
            print 'Warning! CG Alpha is %f' % alpha
            break
        x += alpha * p
        r -= alpha * Ap
        z = r / precon
        rms = ((np.vdot(r, r) / rel_tol) ** 0.5).real
        print '      Micro Iteration %5d: Rel. RMS = %1.5e' % (rot_iter + 1, rms)
        if rms < 0.0001:
            break
        beta = np.vdot(r, z) / rz_old
        if beta < 0:
            print 'Warning! CG Beta is %f, capping at 0' % beta
            beta = 0
        p = z + beta * p

    step = x

    # CI update
    dvecs1.np[:] = step[:ndet]
    dvec.axpy(1.0, dvecs1, 0, 0)
    norm = dvec.norm(0)
    dvec.symnormalize(1.0 / norm, 0)

    # Orbital step
    orb_step = step[ndet:].reshape(noa, nav)
    orb_step[docc:, :act] = 0
    orb_step = psi4.Matrix.from_array(orb_step)
    orb_step.np[:] *= -1

    return (True,
     rot_iter,
     rot_iter,
     orb_step)
# okay decompiling quadratic_step.pyc
