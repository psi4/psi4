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

# Shoody rel imports
import sys
import inspect
path_dir = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe()))[0], '../../')))
sys.path.append(path_dir)

import p4util
from p4util.exceptions import *
np.set_printoptions(precision=5, linewidth=200, threshold=2000, suppress=True)

def OCHx(mcscf_obj, ciwfn, C0, vector):

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
    dim_docc = ciwfn.get_dimension('DOCC')
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
            matF.nph[h][:dim_docc[h], :] += IF.nph[h][:, :dim_docc[h]].T
            matF.nph[h][:dim_docc[h], :] *= ci_overlap
            matF.nph[h][:dim_docc[h], :] += rot_AF.nph[h][:, :dim_docc[h]].T
            matF.nph[h][:dim_docc[h], :] *= 2

        asl = slice(dim_docc[h], dim_oa[h])
        if dim_oa[h] - dim_docc[h]:
            matF.nph[h][asl, :] = np.dot(c_opdm.nph[h], IF.nph[h][:, asl].T)
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

    # 2.0 * (< ret | H | c > - c * e_ci)
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

    # Build rotated ERI's
    ciwfn.rotate_mcscf_integrals(k, onel_k, twoel_k)

    # Sigma computation
    ciwfn.sigma(C0, output, 0, 0, onel_k, twoel_k)
    output.scale(2.0, 0)


def compute_vector(ci_p, orb_p, mcscf_obj, ciwfn, orb_grad, ci_grad, C0, dvecs1, dvecs2):

    # Build up a lot of scratch space because im bad
    scratch1 = ciwfn.new_civector(1, 0, False, True)
    scratch2 = ciwfn.new_civector(1, 0, False, True)
    ci_Ap = ciwfn.new_civector(1, 0, False, True)
    C0vec = ciwfn.new_civector(1, 0, False, True)

    # Set data
    ndet = C0.shape[0]
    C0vec.np[:] = C0

    # Project P right
    orig_ci_scale_cp = C0vec.vdot(ci_p, 0, 0)
    ci_p.axpy(-orig_ci_scale_cp, C0vec, 0, 0)
    ci_p.symnormalize(1.0, 0)
    ci_scale_cp = np.vdot(C0, ci_p.np)

    # Grab orbital info
    nact = sum(ciwfn.get_dimension('ACT').to_tuple())
    ndocc = sum(ciwfn.get_dimension('DRC').to_tuple())
    nvir = sum(ciwfn.get_dimension('VIR').to_tuple())

    # Build ci_Ap
    CCHx(ciwfn, ci_p, ci_Ap, mcscf_obj.current_ci_energy())
    ci_Ap.axpy(2.0 * ci_scale_cp, C0vec, 0, 0)
    COHx(ciwfn, C0vec, orb_p, scratch1, nact)
    ci_Ap.axpy(1.0, scratch1, 0, 0)
    scale = ci_p.vdot(ci_grad, 0, 0)
    scale += 2 * orb_grad.vector_dot(orb_p)
    ci_Ap.axpy(-scale, C0vec, 0, 0)
    ci_Ap.axpy(-ci_scale_cp, ci_grad, 0, 0)

    # Build orb_Ap
    orb_Ap = OCHx(mcscf_obj, ciwfn, C0vec, ci_p)
    orb_Ap.add(mcscf_obj.compute_Hk(orb_p))
    orb_Ap.axpy(-2.0 * ci_scale_cp, orb_grad)
    mcscf_obj.zero_redundant(orb_Ap)

    # Project P left
    proj_scale_cp = ci_Ap.vdot(C0vec, 0, 0)
    ci_Ap.axpy(-proj_scale_cp, C0vec, 0, 0)

    return ci_Ap, orb_Ap


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
    CI_elec_energy = mcscf_obj.current_total_energy()
    orb_grad = mcscf_obj.gradient()
    norm_ci_grad = ci_grad.norm(0)
    dvecs1_buff[:] = C0
    dvecs1.symnormalize(1.0, 0)
    C0[:] = np.asarray(dvecs1_buff)


    ### Preconditioners

    # CI preconditioner
    # 2 * (Hd - E + C0 * C0) - 2.0 * C0 * ci_grad
    ci_precon = ciwfn.Hd_vector(0)
    ci_precon.vector_multiply(1.0, dvecs1, dvecs2, 0, 0, 0)
    ci_precon.shift(-mcscf_obj.current_total_energy(), 0)
    ci_precon.scale(2.0, 0)
    ci_precon.vector_multiply(-2.0, dvecs1, ci_grad, 0, 0, 0)

    orb_precon = mcscf_obj.H_approx_diag()

    ### Compute initial X
    ci_x = ciwfn.new_civector(1, 0, False, True)
    ci_x.copy(ci_grad, 0, 0)
    ci_x.divide(ci_precon, 1.e-8, 0, 0)
    ci_x.scale(-1.0, 0)

    orb_x = orb_grad.clone()
    orb_x.apply_denominator(orb_precon)
    orb_x.scale(-1.0)

    ### Ax product
    ci_Ax, orb_Ax = compute_vector(ci_x, orb_x, mcscf_obj, ciwfn, orb_grad, ci_grad, C0, dvecs1, dvecs2)

    ### Compute R, Z, P
    ci_r = ciwfn.new_civector(1, 0, False, True)
    ci_r.copy(ci_grad, 0, 0)
    ci_r.scale(-1.0, 0)
    ci_r.axpy(-1.0, ci_Ax, 0, 0)
    ci_r.np[np.abs(ci_r.np) < 1.e-14] = 0

    orb_r = orb_grad.clone()
    orb_r.scale(-1.0)
    orb_r.subtract(orb_Ax)

    ci_z = ciwfn.new_civector(1, 0, False, True)
    ci_z.copy(ci_r, 0, 0)
    ci_z.divide(ci_precon, 1.e-12, 0, 0)

    orb_z = orb_r.clone()
    orb_z.apply_denominator(orb_precon)

    ci_p = ciwfn.new_civector(1, 0, False, True)
    ci_p.copy(ci_z, 0, 0)

    orb_p = orb_z.clone()

    ci_rel_tol = ci_r.norm(0) ** 2
    orb_rel_tol = orb_r.sum_of_squares()
    rel_tol = ci_rel_tol + orb_rel_tol

    psi4.print_out("\n")
    psi4.print_out("                 Starting CI RMS = %1.4e   ORB RMS = %1.4e\n" %
                                        (ci_rel_tol ** 0.5, orb_rel_tol ** 0.5))
    converged = False

    ### CG loops
    for rot_iter in range(20):

        # Compute rz_old
        rz_old = ci_r.vdot(ci_z, 0, 0)
        rz_old += orb_r.vector_dot(orb_z)

        # Compute Hessian vector product
        ci_p.np[:ndet][np.abs(ci_p.np) < 1e-15] = 0
        ci_Ap, orb_Ap = compute_vector(ci_p, orb_p, mcscf_obj, ciwfn, orb_grad, ci_grad, C0, dvecs1, dvecs2)

        # Compute alpha
        alpha_denom = ci_Ap.vdot(ci_p, 0, 0)
        alpha_denom += orb_Ap.vector_dot(orb_p)
        alpha = rz_old / alpha_denom

        if alpha > 15:
            psi4.print_out("Warning! CG Alpha is %f\n" % alpha)
            break

        # Update x, r, z
        ci_x.axpy(alpha, ci_p, 0, 0)
        ci_r.axpy(-alpha, ci_Ap, 0, 0)
        ci_z.copy(ci_r, 0, 0)
        ci_z.divide(ci_precon, 1.e-15, 0, 0)

        orb_x.axpy(alpha, orb_p)
        orb_r.axpy(-alpha, orb_Ap)
        orb_z = orb_r.clone()
        orb_z.apply_denominator(orb_precon)

        # Compute relative RMS and print
        ci_numer = ci_r.norm(0) ** 2
        orb_numer = orb_r.sum_of_squares()
        rms_numer = ci_numer + orb_numer

        ci_rms = (ci_numer / ci_rel_tol) ** 0.5
        orb_rms = (orb_numer / orb_rel_tol) ** 0.5
        rms = (rms_numer / rel_tol) ** 0.5


        psi4.print_out("      Micro Iter %2d: Rel. CI RMS = %1.4e   ORB RMS = %1.4e\n" %
                                            (rot_iter + 1, ci_rms, orb_rms))
        if rms < 0.0001:
            psi4.print_out("                 MCSCF microiteration have converged!\n\n")
            converged = True
            break

        # Compute beta
        beta_numer = ci_r.vdot(ci_z, 0, 0)
        beta_numer += orb_r.vector_dot(orb_z)
        beta = beta_numer / rz_old

        if beta < 0:
            psi4.print_out("Warning! CG Beta is %f, capping at 0\n" % beta)
            beta = 0

        # Update p
        ci_p.scale(beta, 0)
        ci_p.axpy(1.0, ci_z, 0, 0)

        orb_p.scale(beta)
        orb_p.axpy(1.0, orb_z)

    if not converged:
        psi4.print_out("                 Warning! MCSCF OS microiterations did not converge.\n\n")

    # CI update
    dvec.axpy(1.0, ci_x, 0, 0)
    norm = dvec.norm(0)
    dvec.symnormalize(1.0 / norm, 0)

    # Orbital step
    orb_x.scale(-1.0)
    mcscf_obj.zero_redundant(orb_x)

    rot_iter *= 2
    rot_iter += 1

    return (True, rot_iter, rot_iter, orb_x)
