import numpy as np
import psi4
import os
from MCSCF import *

# Relative hack for now
import sys, inspect
path_dir = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../../")))
sys.path.append(path_dir)
import p4util
from p4util.exceptions import *
np.set_printoptions(precision=5, linewidth=200, threshold=2000, suppress=True)

def OCHx(mcscf_obj, ciwfn, C0, vector, mcscf):

    # Build the density matrices we need
    c_opdm = ciwfn.opdm(C0, vector, 0, 0)[2]
    c_tpdm = ciwfn.tpdm(C0, vector, 0, 0)[3]
    nirrep = c_opdm.nirrep()

    # Transpose OPDM
    for h in range(nirrep):
        c_opdm.nph[h][:] += c_opdm.nph[h].T.copy()

    # Symmetrize and tranpose TPDM
    tmp_tpdm = np.array(c_tpdm) 
    tmp_tpdm += np.einsum('pqrs->qprs', c_tpdm)
    tmp_tpdm += np.einsum('pqrs->rspq', tmp_tpdm).copy()
    tmp_tpdm *= 0.5

    c_tpdm.np[:] = tmp_tpdm
    c_tpdm.np[:] += tmp_tpdm.T

    dim_docc = ciwfn.get_dimension("DRC")
    dim_oa = ciwfn.get_dimension("OA")
    dim_av = ciwfn.get_dimension("AV")
    dim_rot = ciwfn.get_dimension("ROT")

    IF = mcscf_obj.current_IFock()
    rot_AF = mcscf_obj.compute_AFock(c_opdm)
    rot_Q = mcscf_obj.compute_Q(c_tpdm)

    ci_overlap = C0.vdot(vector, 0, 0)

    matF = psi4.Matrix("OCFx", dim_rot, dim_rot)
    matgrad = psi4.Matrix("OCFx", dim_rot, dim_rot)
    ret = psi4.Matrix("OCFx", dim_oa, dim_av)

    # Construct rotated generalized Fock matrix and result
    for h in range(nirrep):
        # F_in = 2.0 * (< i | 0 > IF + AF^k)

        if dim_docc[h]:
            matF.nph[h][:dim_docc[h], :] += IF.nph[h][:, :dim_docc[0]].T
            matF.nph[h][:dim_docc[h], :] *= ci_overlap 
            matF.nph[h][:dim_docc[h], :] += rot_AF.nph[h][:, :dim_docc[0]].T
            matF.nph[h][:dim_docc[h], :] *= 2

        # F_vn = IF_nw \gamma_wv + Qk
        asl = slice(dim_docc[h], dim_oa[h])
        if dim_oa[h] - dim_docc[h]:
            matF.nph[h][asl, :] = np.dot(c_opdm, IF.nph[h][:, asl].T)
            matF.nph[h][asl, :] += rot_Q.nph[h]

        matgrad.nph[h][:] = matF.nph[h].T
        matgrad.nph[h][:] -= matF.nph[h]
        matgrad.nph[h][:] *= -1
        ret.nph[h][:] = matgrad.nph[h][:dim_oa[h], dim_docc[h]:]
        ret.nph[h][:] *= 2

    mcscf_obj.zero_redundant(ret)
    return ret

def CCHx(ciwfn, vector, output, ci_energy):

    ciwfn.sigma(vector, output, 0, 0)
    output.axpy(-ci_energy, vector, 0, 0)
    output.scale(2.0, 0)
    return np.array(output)

def COHx(ciwfn, C0, k, output, nact):
    ncitri = nact * (nact + 1) / 2
    ncitri2 = (ncitri * (ncitri + 1) ) / 2

    onel_k = psi4.Vector(ncitri)
    twoel_k = psi4.Vector(ncitri2)

    matK = psi4.Matrix(k.shape[0], k.shape[1])
    matK_arr = np.asarray(matK)
#    k.scale(-1.0)
    matK_arr.flat[:] = k

    ciwfn.rotate_mcscf_integrals(matK, onel_k, twoel_k)

    ciwfn.sigma(C0, output, 0, 0, onel_k, twoel_k)
    output.scale(2.0, 0)
#    return output

def compute_vector(proj_x_nr, mcscf_obj, mcscf, ciwfn, eta, orb_grad, ci_grad, C0, rd_bool, dvecs1, dvecs2):
    scratch1 = ciwfn.new_civector(1, 0, False, True)        
    scratch2 = ciwfn.new_civector(1, 0, False, True)        
    ci_Ap = ciwfn.new_civector(1, 0, False, True)        
    C0vec = ciwfn.new_civector(1, 0, False, True)        
    ci_grad_vec = ciwfn.new_civector(1, 0, False, True)        

    ci_grad_vec.np[:] = ci_grad

    # Project out reference state
    proj_x = np.zeros_like(eta)
    proj_x[~rd_bool] = proj_x_nr
    ndet = C0.shape[0] 

    #eta = np.hstack((C0, np.zeros((orb_grad.size))))
    ci_p = ciwfn.new_civector(1, 0, False, True)        
    ci_p.np[:] = proj_x[:ndet]
    C0vec.np[:] = C0

    orig_ci_scale_cp = np.vdot(C0, ci_p.np) 

    ci_p.axpy(-orig_ci_scale_cp, C0vec, 0, 0) 
    ci_p.symnormalize(1.0, 0)
    ci_scale_cp = np.vdot(C0, ci_p.np) 

    orb_p = psi4.Matrix.from_array(proj_x[ndet:].reshape(orb_grad.shape))
    orb_grad = psi4.Matrix.from_array(orb_grad)

    nact = sum(ciwfn.get_dimension("ACT").to_tuple())
    ndocc = sum(ciwfn.get_dimension("DRC").to_tuple())
    nvir = sum(ciwfn.get_dimension("VIR").to_tuple())

    # Build sigma_c = Hc_c
    CCHx(ciwfn, ci_p, ci_Ap, mcscf_obj.current_ci_energy())
    ci_Ap.axpy(2.0 * ci_scale_cp, C0vec, 0, 0)
    COHx(ciwfn, C0vec, orb_p, scratch1, nact)
    ci_Ap.axpy(1.0, scratch1, 0, 0)

    # Subtract out current state
    scale = ci_p.vdot(ci_grad_vec, 0, 0)
    scale += 2 * orb_grad.vector_dot(orb_p) 
    ci_Ap.axpy( -scale, C0vec, 0, 0)
    ci_Ap.axpy(-ci_scale_cp, ci_grad_vec, 0, 0)

    # Build sigm_o = Hc_o
    orb_Ap  = OCHx(mcscf_obj, ciwfn, C0vec, ci_p, mcscf)
    orb_Ap.add(mcscf_obj.compute_Hk(orb_p))
    orb_Ap.axpy(-2.0 * ci_scale_cp, orb_grad)
    orb_Ap = np.array(orb_Ap)
    orb_Ap[ndocc:, :nact] = 0

    # Project out ref state so its symmetric
    proj_scale_cp = ci_Ap.vdot(C0vec, 0, 0)
    ci_Ap.axpy(-proj_scale_cp, C0vec, 0, 0)
   
     
    Ap = np.hstack((ci_Ap.np, orb_Ap.ravel()))


    return Ap[~rd_bool]


def qc_iteration(dvec, ci_grad, ciwfn, mcscf_obj):

    C0 = np.array(dvec)
    ndet = C0.shape[0]

    dvecs1 = ciwfn.new_civector(1, 0, False, True)
    dvecs1_buff = np.asarray(dvecs1)
    dvecs2 = ciwfn.new_civector(1, 0, False, True)
    dvecs2_buff = np.asarray(dvecs2)

    # Update MCSCF object
    dvecs1_buff[:] = C0
    dvecs2_buff[:] = C0

    # Build H, CCgradient, and CC energy
    Hd = np.array(ciwfn.Hd_vector(0))

    opdm = np.array(ciwfn.get_opdm(-1, -1, "SUM", False))
    tpdm = np.array(ciwfn.get_tpdm("SUM", True))

    docc = sum(ciwfn.get_dimension("DOCC").to_tuple())
    act = sum(ciwfn.get_dimension("ACT").to_tuple())
    vir = sum(ciwfn.get_dimension("VIR").to_tuple())

    noa = docc + act
    nav = act + vir

    mcscf = None
    #mcscf_energy = mcscf.update(np.array(ciwfn.get_orbitals("ALL")), opdm, tpdm)

    ci_grad = np.array(ci_grad)

    CI_elec_energy = mcscf_obj.current_total_energy()

    orb_grad = np.array(mcscf_obj.gradient())
    diag_hess = np.array(mcscf_obj.H_approx_diag())

    ci_tol = 1.e-14
    orb_tol = 1.e-14
   # rd_bool = np.hstack((np.abs(ci_grad) < ci_tol, np.zeros_like(orb_grad, dtype=np.bool).ravel()))
    rd_bool = np.hstack((np.abs(ci_grad) < ci_tol, np.abs(orb_grad).ravel() < orb_tol))
    rd_bool[ndet:] = False
    rd_bool[:] = False

    # Check for convergence and print
    norm_ci_grad = np.linalg.norm(ci_grad)
    norm_orb_grad = np.linalg.norm(orb_grad)

    dvecs1_buff[:] = C0
    dvecs1.symnormalize(1.0, 0)
    C0[:] = np.asarray(dvecs1_buff)

    eta = np.hstack((C0, np.zeros((orb_grad.size))))
    grad = -np.hstack((ci_grad, orb_grad.ravel()))

    RHS = grad.ravel()[~rd_bool]

    ci_precon  = 2 * (Hd - mcscf_obj.current_total_energy() + C0 * C0)
    ci_precon -= 2.0 * C0 * ci_grad
    #ci_precon[~crdm] = 1
    #print Hd
    #print mcscf_obj.current_ci_energy()
    #print 'ci_precon'
    #print ci_precon

    orb_precon = diag_hess.copy()

    precon = np.hstack((ci_precon, orb_precon.ravel()))[~rd_bool]
    #xCI = -ci_grad / ci_precon
    #xORB = - orb_grad / orb_precon

    x = RHS / precon
    x[:ndet][np.abs(x[:ndet]) < 1.e-15] = 0
    Ax = compute_vector(x, mcscf_obj, mcscf, ciwfn, eta, orb_grad, ci_grad, C0, rd_bool, dvecs1, dvecs2)
    r = RHS - Ax
    z = r / precon
    p = z.copy()
    r[:ndet][np.abs(x[:ndet]) < 1.e-14] = 0

    rel_tol = np.vdot(r, r)
    fx = np.zeros_like(rd_bool, dtype=np.double)
    fx[~rd_bool] = Ax
    pprint = False
    if pprint:
        print '----'
        print fx[:ndet] * 1.e4
        print fx[ndet:].reshape(orb_grad.shape) * 1.e2
        print orb_grad * 1.e2
        #raise Exception("")


    for rot_iter in range(12):
        rz_old = np.vdot(r, z)

        p[:ndet][np.abs(p[:ndet]) < 1.e-15] = 0
        Ap = compute_vector(p, mcscf_obj, mcscf, ciwfn, eta, orb_grad, ci_grad, C0, rd_bool, dvecs1, dvecs2)


        alpha = rz_old / np.vdot(Ap, p)
        if alpha > 15:
            print('Warning! CG Alpha is %f' % alpha)
            break

        x += alpha * p
        r -= alpha * Ap
        z = r / precon

        rms = ((np.vdot(r, r) / rel_tol) ** 0.5).real

        print('      Micro Iteration %5d: Rel. RMS = %1.5e' %  (rot_iter + 1, rms))
        if rms < 1.e-8:
            break

        beta = np.vdot(r, z) / rz_old
        if beta < 0:
            print('Warning! CG Beta is %f, capping at 0' % beta)
            beta = 0

        p = z + beta * p


    step = np.zeros_like(eta)
    step[~rd_bool] = x

    dvecs1.np[:] = step[:ndet]

    dvec.axpy(1.0, dvecs1, 0, 0)
    norm = dvec.norm(0)
    dvec.symnormalize(1.0 / norm, 0)

    orb_step = step[ndet:].reshape(noa, nav)
    orb_step[docc:, :act] = 0    
    
    orb_step = psi4.Matrix.from_array(orb_step)
    if pprint:
        print orb_grad
        Ax = compute_vector(x, mcscf_obj, mcscf, ciwfn, eta, orb_grad, ci_grad, C0, rd_bool, dvecs1, dvecs2)
        ret = np.zeros_like(eta)
        ret[~rd_bool] = Ax
        print ret[ndet:].reshape(noa, nav)

        print np.array(step[:ndet]) * 1.e4
        print np.array(orb_step) * 1.e2
        raise Exception("")

    orb_step.np[:] *= -1

    return True, rot_iter, rot_iter, orb_step

