"""
This file tests electric_field_value and induction_operator
agains the canonical electric_field integral evaluation
"""
import numpy as np
import pytest
import psi4


pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def test_elec_fields():
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
    mints = psi4.core.MintsHelper(basis_obj)
    # generate random points and dipole moments
    coords = 5 * np.random.rand(50, 3)
    moments = 0.1 * np.random.rand(50, 3)

    # run the new implementation
    coordinates = psi4.core.Matrix.from_array(coords)
    dips = psi4.core.Matrix.from_array(moments)
    ret = mints.induction_operator(coordinates, dips).np

    # old implementation (used for EFP, for example)
    nbf = mints.basisset().nbf()
    V2 = np.zeros((nbf, nbf))

    field_ints = np.zeros((3, nbf, nbf))
    for c, m in zip(coords, moments):
        # get electric field integrals from Psi4
        p4_field_ints = mints.electric_field(origin=c)
        for pole in range(3):
            field_ints[pole] = np.asarray(p4_field_ints[pole])

        # scale field integrals by induced dipole magnitudes.
        for pole in range(3):
            field_ints[pole] *= -m[pole]
            V2 += field_ints[pole]
    np.testing.assert_allclose(V2, ret)

    # electric field expectation values
    nbf = basis_obj.nbf()
    mock_dmat = psi4.core.Matrix.from_array(np.random.rand(nbf, nbf))
    field_val = mints.electric_field_value(coordinates, mock_dmat).np

    # Electric field at points
    points = coords
    npt = len(points)
    field_ref = np.zeros((npt, 3))
    for ipt in range(npt):
        p4_field_ints = mints.electric_field(origin=points[ipt])
        field_ref[ipt] = [
            np.vdot(mock_dmat.np, np.asarray(p4_field_ints[0])),  # Ex
            np.vdot(mock_dmat.np, np.asarray(p4_field_ints[1])),  # Ey
            np.vdot(mock_dmat.np, np.asarray(p4_field_ints[2]))   # Ez
        ]
    np.testing.assert_allclose(field_ref, field_val)
