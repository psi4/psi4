"""
This file tests the ao_multipole_potential integrals
using finite differences
"""
import numpy as np
import pytest
import psi4


pytestmark = pytest.mark.quick


def test_potential_erf_integrals():
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
    basis_obj = psi4.core.BasisSet.build(mol, 'ORBITAL', "cc-pvtz")
    mints = psi4.core.MintsHelper(basis_obj)
    
    # corner cases
    # erf(inf) = 1
    x = mints.ao_potential_erf([1, 2, 3], omega=1e20).np
    y = mints.ao_multipole_potential([1, 2, 3], max_k=0)[0].np
    np.testing.assert_allclose(x, y, atol=1e-14)

    # erfc(0) = 1
    x = mints.ao_potential_erf_complement([1, 2, 3], omega=0.0).np
    y = mints.ao_multipole_potential([1, 2, 3], max_k=0)[0].np
    np.testing.assert_allclose(x, y, atol=1e-14)

    # 1/R - erf(R)/R - erfc(R)/R = 0
    erf = mints.ao_potential_erf([1, 2, 3], omega=1.5).np
    erfc = mints.ao_potential_erf_complement([1, 2, 3], omega=1.5).np
    diff = y - erf - erfc
    np.testing.assert_allclose(diff, np.zeros_like(diff), atol=1e-14)