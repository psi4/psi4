"""
This file tests the F12 two-electron integrals from libmints
"""
import numpy as np
import pytest
import psi4
from pathlib import Path
import json


pytestmark = [pytest.mark.api, pytest.mark.quick]


@pytest.fixture
def reference_data():
    with open(Path(__file__).parent / "f12_libint1.json") as f:
        reference_data = json.load(f)
    return reference_data


def matlist_to_ndarray(mats):
    """Converts a list of psi4.core.Matrix to np.ndarray"""
    return psi4.driver.p4util.numpy_helper._to_array(mats) 


def test_f12_integrals(reference_data):
    assert reference_data['version'] == '1.5'
    refdata = reference_data['data']
    moldata = refdata['h2o']
    mol = psi4.core.Molecule.from_string(moldata.pop('psi4string'))

    basis = psi4.core.BasisSet.build(mol, 'orbital', '6-31G*')
    mints = psi4.core.MintsHelper(basis)
    f12 = mints.f12_cgtg(1.0)

    for int_type in ['F', 'F2', 'FG', 'Uf']:
        ref = moldata[int_type]
        shape = ref[f'shape']
        integral_ref = np.array(ref["I"]).reshape(shape)

        mapping = {
            "F": mints.ao_f12,
            "F2": mints.ao_f12_squared,
            "FG": mints.ao_f12g12,
            "Uf": mints.ao_f12_double_commutator,
        }
        libint2 = mapping[int_type](f12)

        libint2_np = matlist_to_ndarray(libint2)
        np.testing.assert_allclose(libint2_np, integral_ref, atol=1.e-14)
