import numpy as np
import pytest

import psi4


pytestmark = [pytest.mark.psi, pytest.mark.api]


def test_petitelist_i1763_shapes_and_inverse_behavior():
    mol = psi4.geometry(
        """
0 3
symmetry c1
C  0.0000000000  0.0000000000 -0.5928430915
H -0.0000000000  0.9469373770 -1.1509808737
H  0.0000000000 -0.9469373770 -1.1509808737
"""
    )

    basis = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pvdz", quiet=True)
    mints = psi4.core.MintsHelper(basis)

    petite_true = mints.petite_list1(True)
    so2ao_true = petite_true.sotoao().nph[0]
    ao2so_true = petite_true.aotoso().nph[0]

    assert so2ao_true.shape == (basis.nbf(), basis.nao())
    assert ao2so_true.shape == (basis.nao(), basis.nbf())
    assert so2ao_true.shape[1] > so2ao_true.shape[0]

    id_like_true = np.dot(so2ao_true, ao2so_true)
    assert not np.allclose(id_like_true, np.eye(basis.nbf()), atol=1.0e-12)

    petite_false = mints.petite_list1(False)
    so2ao_false = petite_false.sotoao().nph[0]
    ao2so_false = petite_false.aotoso().nph[0]

    assert so2ao_false.shape == (basis.nbf(), basis.nbf())
    assert ao2so_false.shape == (basis.nbf(), basis.nbf())

    id_like_false = np.dot(so2ao_false, ao2so_false)
    assert np.allclose(id_like_false, np.eye(basis.nbf()), atol=1.0e-12)


def test_cartao_to_ao_transform_matches_petitelist_true():
    mol = psi4.geometry(
        """
0 3
symmetry c1
C  0.0000000000  0.0000000000 -0.5928430915
H -0.0000000000  0.9469373770 -1.1509808737
H  0.0000000000 -0.9469373770 -1.1509808737
"""
    )

    basis = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pvdz", quiet=True)
    mints = psi4.core.MintsHelper(basis)

    cart_to_ao = mints.cartao_to_ao_transform().nph[0]
    aotoso_true = mints.petite_list1(True).aotoso().nph[0]
    aotoso_false = mints.petite_list1(False).aotoso().nph[0]

    assert np.allclose(aotoso_true, cart_to_ao.T @ aotoso_false, atol=1.0e-12)


def test_sphere_puream_transform_sentinel():
    water = psi4.geometry(
        """
O  0.000000000000  0.000000000000 -0.075791843589
H  0.000000000000 -0.866811828967  0.601435779270
H  0.000000000000  0.866811828967  0.601435779270
symmetry c1
no_reorient
no_com
"""
    )

    psi4.set_options(
        {
            "basis": "cc-pvdz",
            "reference": "rhf",
            "scf_type": "df",
            "damping_percentage": 0,
            "perturb_h": True,
            "perturb_with": "sphere",
            "theta_points": 10,
            "phi_points": 10,
        }
    )

    energy = psi4.energy("scf", molecule=water)
    assert psi4.compare_values(-76.04926605031216, energy, 6, "puream sphere perturbation sentinel")
