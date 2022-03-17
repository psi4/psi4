"""
This file tests the integrals using the McMurchie-Davidson scheme
against hard-coded values from previous Psi4 versions for consistency check
"""
import numpy as np
import pytest
import psi4
from pathlib import Path
import json
from dataclasses import dataclass


pytestmark = pytest.mark.quick


@pytest.fixture
def reference_data():
    # Reference data generated using Psi4 1.5, still using OS86
    with open(Path(__file__).parent / "oei_reference_data.json") as f:
        reference_data = json.load(f)
    return reference_data


@dataclass
class IntegralInfo:
    dictkey: str
    function_name: str
    extra_args: list = None


integrals = [
    IntegralInfo("L", "ao_angular_momentum"),
]


def test_mcmurchie_davidson_consistency(reference_data):
    assert reference_data['version'] == '1.5'
    refdata = reference_data['data']
    for molname in refdata:
        moldata = refdata[molname]
        mol = psi4.core.Molecule.from_string(moldata.pop('psi4string'))
        for basis_name in moldata:
            ref = moldata[basis_name]
            basis = psi4.core.BasisSet.build(mol, 'orbital', basis_name)
            mints = psi4.core.MintsHelper(basis)
            for integral in integrals:
                shape = ref[f'{integral.dictkey}_shape']
                integral_ref = np.array(ref[integral.dictkey]).reshape(shape)
                psi4_fun = getattr(mints, integral.function_name)
                if integral.extra_args is not None:
                    ret = psi4_fun(*integral.extra_args)
                else:
                    ret = psi4_fun()
                ret_np = np.array([x.np for x in ret])
                np.testing.assert_allclose(
                    ret_np, integral_ref, atol=1e-13,
                    err_msg=f'Consistency error in {integral.function_name}.'
                )