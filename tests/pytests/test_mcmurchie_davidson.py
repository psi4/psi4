"""
This file tests one-electron integrals from libmints computed with
the McMurchie-Davidson scheme
"""
import numpy as np
import pytest
import psi4
from pathlib import Path
import json


pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


@pytest.fixture
def reference_data():
    # Reference data generated using Psi4 1.5, still using OS86
    with open(Path(__file__).parent / "oei_reference_data.json") as f:
        reference_data = json.load(f)
    return reference_data


def test_mcmurchie_davidson_consistency_angmom(reference_data):
    assert reference_data['version'] == '1.5'
    refdata = reference_data['data']
    for molname in refdata:
        moldata = refdata[molname]
        mol = psi4.core.Molecule.from_string(moldata.pop('psi4string'))
        for basis_name in moldata:
            ref = moldata[basis_name]
            basis = psi4.core.BasisSet.build(mol, 'orbital', basis_name)
            mints = psi4.core.MintsHelper(basis)

            shape = ref[f'L_shape']
            integral_ref = np.array(ref["L"]).reshape(shape)
            ret = mints.ao_angular_momentum()
            ret_np = np.array([x.np for x in ret])
            np.testing.assert_allclose(ret_np, integral_ref, atol=1e-14)


def test_mcmurchie_davidson_multipoles():
    mol = psi4.geometry("""
        units bohr
        O  0.000000000000  0.000000000000 -0.149544924691
        H  0.000000000000 -1.336238163149  1.186693238458
        H  0.000000000000  1.336238163149  1.186693238458
        no_com
        no_reorient
        symmetry c1
    """)
    basis = psi4.core.BasisSet.build(mol, 'orbital', 'cc-pvtz')
    mints = psi4.core.MintsHelper(basis)
    order = 6
    M = mints.ao_multipoles(origin=[0.0, 0.0, 0.0], order=order)
    assert len(M) == (order + 1) * (order + 2) * (order + 3) / 6 - 1

    to_ndarray = lambda mats: np.array([x.np for x in mats])

    Mnp = to_ndarray(M)
    # reference from l2
    dips = to_ndarray(mints.ao_dipole())
    quads = to_ndarray(mints.ao_quadrupole())
    
    np.testing.assert_allclose(dips, Mnp[:3], atol=1e-14)
    np.testing.assert_allclose(quads, Mnp[3:9], atol=1e-14)


def test_mcmurchie_davidson_multipoles_order10():
    mol = psi4.geometry("""
        units bohr
        O  0.000000000000  0.000000000000 -0.149544924691
        H  0.000000000000 -1.336238163149  1.186693238458
        H  0.000000000000  1.336238163149  1.186693238458
        no_com
        no_reorient
        symmetry c1
    """)
    basis = psi4.core.BasisSet.build(mol, 'orbital', 'sto-3g')
    psi4.set_options({'puream' : False})
    mints = psi4.core.MintsHelper(basis)
    order = 10
    M = mints.ao_multipoles(origin=[1.0, 2.0, 3.0], order=order)
    assert len(M) == (order + 1) * (order + 2) * (order + 3) / 6 - 1
    to_ndarray = lambda mats: np.array([x.np for x in mats])
    Mnp = to_ndarray(M)
    for i, m in enumerate(M):
        print(m.name, np.max(Mnp[i]), Mnp[i].shape)