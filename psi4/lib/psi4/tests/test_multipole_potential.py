"""
This file tests the ao_multipole_potential integrals
using finite differences
"""
import numpy as np
import itertools
import pytest
import psi4


pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def multipole_offset(k):
    o = 0
    for i in range(k):
        o += (i + 1) * (i + 2) / 2
    return int(o)


def matlist_to_ndarray(mats):
    """Converts a list of psi4.core.Matrix to np.ndarray"""
    return np.array([x.np for x in mats])


def compare_charge_integrals(mints, point):
    pot = psi4.core.ExternalPotential()
    pot.addCharge(-1, *point)
    ref = pot.computePotentialMatrix(mints.basisset()).np
    test = mints.ao_multipole_potential(order=0, origin=point)[0].np
    np.testing.assert_allclose(ref, test)


def compare_field_integrals(mints, point):
    # electric_field is evaluated using Libint2
    ref = matlist_to_ndarray(mints.electric_field(point))
    test = matlist_to_ndarray(
        mints.ao_multipole_potential(order=1, origin=point)[1:]
    )
    np.testing.assert_allclose(ref, test, atol=1e-14)


def compare_arbitrary_order_fd(mints, point, step=1e-5):
    """Compares the integral of order m_order with finite differences
    of m_order - 1."""
    for m_order in range(1, 6):
        o = multipole_offset(m_order)
        o_prev = multipole_offset(m_order - 1)
        test = matlist_to_ndarray(
            mints.ao_multipole_potential(order=m_order, origin=point)[o:]
        )

        components = list(itertools.combinations_with_replacement((0, 1, 2), m_order))
        comp_grad = list(itertools.combinations_with_replacement((0, 1, 2), m_order - 1))

        prefactors_5p = np.array([1.0, -8.0, 8.0, -1.0]) / 12.0
        multipliers_5p = [-2, -1, 1, 2]
        grad_fd = np.zeros_like(test)
        for icomp, cc in enumerate(components):
            for f, p in zip(multipliers_5p, prefactors_5p):
                pert = point.copy()
                pert[cc[-1]] += f * step
                idxp = o_prev + comp_grad.index(tuple(cc[:-1]))
                der = mints.ao_multipole_potential(m_order - 1, pert)[idxp].np
                grad_fd[icomp] += p * der / step
        np.testing.assert_allclose(grad_fd, test, atol=1e-9)


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
    point = [1.0, 2.0, 3.0]
    compare_charge_integrals(mints_helper, point=point)
    compare_field_integrals(mints_helper, point=point)
    compare_arbitrary_order_fd(mints_helper, point=point)
