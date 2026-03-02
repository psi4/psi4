"""
This file tests the ZORA spin-orbit coupling integrals from libmints
"""
import numpy as np
np.set_printoptions(suppress=True, linewidth=200, precision = 6)
import pytest
import psi4
from pathlib import Path
import json

pytestmark = [pytest.mark.api, pytest.mark.quick]


def test_zora_integrals():
    to_array = psi4.driver.p4util.numpy_helper._to_array
    mol = psi4.core.Molecule.from_string("""
    0 1
        O 0 0 0
    nocom
    noreorient
    """)

    basis = psi4.core.BasisSet.build(mol, 'orbital', '3-21G')
    mints = psi4.core.MintsHelper(basis)
    ints = mints.ao_zora_spin_orbit()

    Hx, Hy, Hz = tuple([to_array(m) for m in ints])
    Hso = np.block([[Hz, Hx + 1j*Hy],[Hx-1j*Hy, -Hz]])
    np.testing.assert_allclose(Hso, -Hso.conj().T, atol=1e-10)
