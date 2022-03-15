"""
This file tests the ao_multipole_potential integrals
using finite differences
"""
import numpy as np
import itertools
import pytest
import psi4


pytestmark = pytest.mark.quick


def multipole_offset(k):
    o = 0
    for i in range(k):
        o += (i + 1) * (i + 2) / 2
    return int(o)


def compare_charge_integrals(mints, point):
    pot = psi4.core.ExternalPotential()
    pot.addCharge(-1, *point)
    ref = pot.computePotentialMatrix(mints.basisset()).np
    test = mints.ao_multipole_potential(point, max_k=0)[0].np
    np.testing.assert_allclose(ref, test)


def compare_field_integrals_fd(mints, point, step=1e-5):
    # get the electric field integrals
    o = multipole_offset(1)
    test = mints.ao_multipole_potential(point, max_k=1)[o:]
    assert len(test) == 3
    # loop over components and generate reference with FD
    for c in [0, 1, 2]:
        offset = np.zeros_like(point)
        offset[c] = step
        step_plus = mints.ao_multipole_potential(point + offset, max_k=0)[0]
        step_minus = mints.ao_multipole_potential(point - offset, max_k=0)[0]
        field_finite = (step_plus.np - step_minus.np) / (2 * step)
        np.testing.assert_allclose(field_finite, test[c].np, atol=1e-9)


def compare_field_gradient_fd(mints, point, step=1e-5):
    # get the electric field gradient integrals
    o = multipole_offset(2)
    test = mints.ao_multipole_potential(point, max_k=2)[o:]
    assert len(test) == 6
    components = list(itertools.combinations_with_replacement((0, 1, 2), 2))
    comp_idx = 0
    # loop over components and generate reference
    for c1, c2 in components:
        offset2 = np.zeros_like(point)
        offset2[c2] = step
        step_plus = mints.ao_multipole_potential(
            point + offset2, max_k=1
        )[1 + c1].np
        step_minus = mints.ao_multipole_potential(
            point - offset2, max_k=1
        )[1 + c1].np
        field_grad_finite = (step_plus - step_minus) / (2 * step)
        np.testing.assert_allclose(
            field_grad_finite, test[comp_idx].np, atol=1e-9
        )
        comp_idx += 1


def compare_field_hessian_fd(mints, point, step=1e-5):
    # get the electric field gradient integrals
    o = multipole_offset(3)
    test = mints.ao_multipole_potential(point, max_k=3)[o:]
    assert len(test) == 10
    components = list(itertools.combinations_with_replacement((0, 1, 2), 3))
    comp_grad = list(itertools.combinations_with_replacement((0, 1, 2), 2))
    comp_idx = 0
    # loop over components and generate reference
    for c1, c2, c3 in components:
        offset3 = np.zeros_like(point)
        offset3[c3] = step
        step_plus = mints.ao_multipole_potential(
            point + offset3, max_k=2
        )[4 + comp_grad.index((c1, c2))].np
        step_minus = mints.ao_multipole_potential(
            point - offset3, max_k=2
        )[4 + comp_grad.index((c1, c2))].np
        field_grad_finite = (step_plus - step_minus) / (2 * step)
        np.testing.assert_allclose(
            field_grad_finite, test[comp_idx].np, atol=1e-9
        )
        comp_idx += 1


def test_multipole_potential_integrals():
    mol = psi4.geometry("""
    units bohr
    0 1
    O1     0.000000000000     0.000000000000     0.224348285559
    H2    -1.423528800232     0.000000000000    -0.897393142237
    H3     1.423528800232     0.000000000000    -0.897393142237
    symmetry c1
    no_com
    no_reorient
    """)
    basis_obj = psi4.core.BasisSet.build(mol, 'ORBITAL', "cc-pvdz")
    mints_helper = psi4.core.MintsHelper(basis_obj)
    compare_charge_integrals(mints_helper, point=[1.0, 2.0, 3.0])
    compare_field_integrals_fd(mints_helper, point=[1.0, 2.0, 3.03])
    compare_field_gradient_fd(mints_helper, point=[1.0, 2.0, 3.0])
    compare_field_hessian_fd(mints_helper, point=[1.0, 2.0, 3.0])
