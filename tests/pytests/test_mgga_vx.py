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


@pytest.mark.parametrize("func", ["PBE", "TPSS"])
def test_uks_rks_hessian_consistency(func):
    """Closed shell through the UV path must reproduce the RV explicit
    XC Hessian to machine precision."""
    wfn = _scf(func, "rks", """
        O  0.000  0.000  0.117
        H  0.000  0.755 -0.469
        H  0.000 -0.755 -0.469
        symmetry c1
        noreorient""")
    vpot = wfn.V_potential()
    vpot.set_D([wfn.Da()])
    H_rv = vpot.compute_hessian().to_array()
    sup_u = psi4.driver.dft.build_superfunctional(func, False)[0]
    vpot_u = psi4.core.VBase.build(wfn.basisset(), sup_u, "UV")
    vpot_u.initialize()
    vpot_u.set_D([wfn.Da(), wfn.Da()])
    nbf = wfn.nmo()
    vpot_u.compute_V([psi4.core.Matrix(nbf, nbf), psi4.core.Matrix(nbf, nbf)])
    H_uv = vpot_u.compute_hessian().to_array()
    assert np.allclose(H_rv, H_uv, atol=1e-10)


def test_uks_hessian_explicit_fd():
    """The explicit UKS XC Hessian against mixed finite differences of
    the quadrature energy with both spin densities frozen. The residual
    is the fixed-grid convention (grid-motion + weight classes), a few
    1e-4 at this grid; term bugs show up orders of magnitude larger."""
    geom = """
        0 2
        N  0.000  0.000  0.142
        H  0.000  0.802 -0.496
        H  0.000 -0.802 -0.496
        symmetry c1
        noreorient"""
    psi4.core.clean()
    psi4.set_options({"basis": "6-31G", "reference": "uks", "scf_type": "pk",
                      "dft_spherical_points": 590, "dft_radial_points": 99,
                      "e_convergence": 1e-9, "d_convergence": 1e-8,
                      "maxiter": 200})
    mol = psi4.geometry(geom)
    _, wfn = psi4.energy("TPSS", return_wfn=True, molecule=mol)
    Da, Db = wfn.Da().to_array(), wfn.Db().to_array()
    vpot = wfn.V_potential()
    vpot.set_D([wfn.Da(), wfn.Db()])
    Hx = vpot.compute_hessian().to_array()
    coords0 = mol.geometry().to_array()

    def exc_at(coords):
        m = psi4.geometry(geom)
        m.set_geometry(psi4.core.Matrix.from_array(coords))
        m.update_geometry()
        basis = psi4.core.BasisSet.build(m, "ORBITAL", "6-31G")
        sup = psi4.driver.dft.build_superfunctional("TPSS", False)[0]
        vp = psi4.core.VBase.build(basis, sup, "UV")
        vp.initialize()
        vp.set_D([_mat(Da), _mat(Db)])
        nbf = basis.nbf()
        vp.compute_V([psi4.core.Matrix(nbf, nbf), psi4.core.Matrix(nbf, nbf)])
        return vp.quadrature_values()["FUNCTIONAL"]

    h = 2e-3
    for (A, x, B, y) in [(0, 2, 0, 2), (0, 2, 1, 1)]:
        def E(sa, sb):
            c = coords0.copy()
            c[A, x] += sa * h
            c[B, y] += sb * h
            return exc_at(c)
        if (A, x) == (B, y):
            fd = (E(1, 0) - 2 * E(0, 0) + E(-1, 0)) / h**2
        else:
            fd = (E(1, 1) - E(1, -1) - E(-1, 1) + E(-1, -1)) / (4 * h**2)
        assert abs(Hx[3 * A + x][3 * B + y] - fd) < 2e-3


# ==========================================================================
# End-to-end DFT_GRID_RESPONSE coverage: with grid response on, the analytic
# XC gradient/Hessian become the exact derivatives of the computed energy.
# These exercise the full driver path (CPKS + nuclear + JK + XC), the grid-
# response term classes, and the analytic first/second quadrature-weight
# derivatives -- none of which the fixed-grid unit tests above touch.
# ==========================================================================

# C2v water / NH2, kept symmetric so the driver's finite-difference builds only
# the symmetry-unique displacements. The frame is fixed (noreorient/no_com), as
# grid response requires (the grid rides its parent atoms rigidly).
_WATER = """
    O  0.000  0.000  0.117
    H  0.000  0.755 -0.469
    H  0.000 -0.755 -0.469
    noreorient
    no_com"""

_NH2 = """
    0 2
    N  0.000  0.000  0.142
    H  0.000  0.802 -0.496
    H  0.000 -0.802 -0.496
    noreorient
    no_com"""


def _gr_opts(ref, response):
    # Atomic blocking + orientation off are required for the analytic weight
    # Jacobian; toggling `response` on the same grid isolates the grid-response
    # terms. `points 5` gives the driver a 5-point findif stencil, accurate
    # enough to resolve the ~1e-10 agreement (the default 3-point is not). The
    # grid is deliberately coarse: this is an analytic-vs-findif consistency
    # check on one and the same grid, so density costs rigor nothing.
    return {"basis": "6-31G", "reference": ref, "scf_type": "pk",
            "dft_spherical_points": 110, "dft_radial_points": 30,
            "e_convergence": 1e-10, "d_convergence": 1e-9, "maxiter": 200,
            "dft_block_scheme": "atomic", "dft_grid_orientation": False,
            "dft_grid_response": response, "points": 5}


@pytest.mark.parametrize("ref,geom", [("rks", _WATER), ("uks", _NH2)])
def test_gradient_grid_response(ref, geom):
    """With DFT_GRID_RESPONSE on, the analytic XC gradient is the exact
    derivative of the energy: it matches the driver's 5-point finite-difference
    gradient (dertype=0) to the findif floor, whereas the fixed-grid gradient on
    the same grid is off by the grid-motion/weight convention (~1e-4)."""
    func = "TPSS"
    psi4.core.clean()
    psi4.set_options(_gr_opts(ref, True))
    mol = psi4.geometry(geom)
    g_analytic = np.array(psi4.gradient(func, dertype=1, molecule=mol))

    psi4.core.clean()
    psi4.set_options(_gr_opts(ref, True))
    g_findif = np.array(psi4.gradient(func, dertype=0, molecule=mol))
    assert np.max(np.abs(g_analytic - g_findif)) < 1e-8      # grid-response gradient is exact

    psi4.core.clean()
    psi4.set_options(_gr_opts(ref, False))
    g_fixed = np.array(psi4.gradient(func, dertype=1, molecule=mol))
    assert np.max(np.abs(g_fixed - g_findif)) > 1e-5         # fixed-grid gradient is not


@pytest.mark.parametrize("ref,geom", [("rks", _WATER), ("uks", _NH2)])
def test_hessian_grid_response(ref, geom):
    """End-to-end analytic mGGA Hessian (dertype=2) through the full driver
    (CPKS + nuclear + JK + XC) with grid response on, against the driver's
    5-point finite difference of the analytic gradients (dertype=1). Exercises
    the Hessian grid-response classes, including the second weight derivatives,
    over the whole matrix."""
    func = "TPSS"
    psi4.core.clean()
    psi4.set_options(_gr_opts(ref, True))
    mol = psi4.geometry(geom)
    H_analytic = np.array(psi4.hessian(func, dertype=2, molecule=mol))
    assert np.allclose(H_analytic, H_analytic.T, atol=1e-8)

    psi4.core.clean()
    psi4.set_options(_gr_opts(ref, True))
    H_findif = np.array(psi4.hessian(func, dertype=1, molecule=mol))
    assert np.max(np.abs(H_analytic - H_findif)) < 1e-6
