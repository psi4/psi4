"""Validate the meta-GGA Vx (fxc) contraction against finite differences
of the XC potential.

The reference is Psi4's own compute_V — whose meta-GGA terms are
long-standing — differentiated numerically along a random density
displacement on the same grid, with Richardson extrapolation:

  RKS singlet:  d/dh Vxc[D + h X]            = 2 compute_Vx(X)
  RKS triplet:  d/dh Vxc_a[(D + hX, D - hX)] = 2 compute_Vx_full(X, singlet=False)
  UKS:          d/dh Vxc_s[(Da + hXa, Db + hXb)] = compute_Vx([Xa, Xb])_s

The factors of two follow the convention that the restricted Vx routines
return the alpha block divided by two. The GGA case is included as a
control of the harness itself.
"""
import numpy as np
import pytest

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, pytest.mark.dft]


def _mat(a):
    return psi4.core.Matrix.from_array(a)


def _fd_v(vpot, dens_plus, dens_minus, h, nbf, comp):
    out = []
    for dens in (dens_plus, dens_minus):
        vpot.set_D([_mat(d) for d in dens])
        vs = [psi4.core.Matrix(nbf, nbf) for _ in dens]
        vpot.compute_V(vs)
        out.append(vs[comp].to_array())
    return (out[0] - out[1]) / (2 * h)


def _richardson(f, h):
    return (4.0 * f(h / 2) - f(h)) / 3.0


def _scf(func, ref, geom):
    psi4.core.clean()
    psi4.set_options({"basis": "6-31G", "reference": ref, "scf_type": "pk",
                      "dft_spherical_points": 194, "dft_radial_points": 40,
                      "e_convergence": 1e-9, "d_convergence": 1e-8,
                      "maxiter": 200})
    mol = psi4.geometry(geom)
    _, wfn = psi4.energy(func, return_wfn=True, molecule=mol)
    return wfn


@pytest.mark.parametrize("func", ["PBE", "SCAN", "TPSS", "M06-L"])
def test_rks_vx_fd(func):
    wfn = _scf(func, "rks", """
        O  0.000  0.000  0.117
        H  0.000  0.755 -0.469
        H  0.000 -0.755 -0.469
        symmetry c1
        noreorient""")
    vpot = wfn.V_potential()
    nbf = wfn.nmo()
    Da = wfn.Da().to_array()
    rng = np.random.default_rng(7)
    X = rng.standard_normal((nbf, nbf)) * 0.01
    X = 0.5 * (X + X.T)

    # singlet
    vpot.set_D([wfn.Da()])
    ret = [psi4.core.Matrix(nbf, nbf)]
    vpot.compute_Vx([_mat(X)], ret)
    fd = _richardson(lambda h: _fd_v(vpot, [Da + h * X], [Da - h * X], h, nbf, 0), 2e-4)
    assert np.allclose(2.0 * ret[0].to_array(), fd, atol=5e-7)

    # triplet, referenced through an unrestricted twin potential
    vpot.set_D([wfn.Da()])
    ret = [psi4.core.Matrix(nbf, nbf)]
    vpot.compute_Vx_full([_mat(X)], ret, False)
    sup_u = psi4.driver.dft.build_superfunctional(func, False)[0]
    vpot_u = psi4.core.VBase.build(wfn.basisset(), sup_u, "UV")
    vpot_u.initialize()
    fd = _richardson(lambda h: _fd_v(vpot_u, [Da + h * X, Da - h * X],
                                     [Da - h * X, Da + h * X], h, nbf, 0), 2e-4)
    assert np.allclose(2.0 * ret[0].to_array(), fd, atol=5e-7)


@pytest.mark.parametrize("func", ["PBE", "SCAN", "TPSS", "M06-L"])
def test_uks_vx_fd(func):
    wfn = _scf(func, "uks", """
        0 2
        N  0.000  0.000  0.000
        H  0.000  0.800  0.600
        H  0.000 -0.800  0.600
        symmetry c1
        noreorient""")
    vpot = wfn.V_potential()
    nbf = wfn.nmo()
    Da, Db = wfn.Da().to_array(), wfn.Db().to_array()
    rng = np.random.default_rng(11)
    Xa = rng.standard_normal((nbf, nbf)) * 0.01
    Xb = rng.standard_normal((nbf, nbf)) * 0.01
    Xa, Xb = 0.5 * (Xa + Xa.T), 0.5 * (Xb + Xb.T)

    vpot.set_D([wfn.Da(), wfn.Db()])
    ret = [psi4.core.Matrix(nbf, nbf), psi4.core.Matrix(nbf, nbf)]
    vpot.compute_Vx([_mat(Xa), _mat(Xb)], ret)
    for comp in (0, 1):
        fd = _richardson(lambda h: _fd_v(vpot, [Da + h * Xa, Db + h * Xb],
                                         [Da - h * Xa, Db - h * Xb], h, nbf, comp), 2e-4)
        assert np.allclose(ret[comp].to_array(), fd, atol=5e-7)
