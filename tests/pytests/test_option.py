import numpy as np
import pytest

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


@pytest.mark.smoke
def test_spacious_option():
    """tu1-h2o-energy"""
    #! Sample HF/cc-pVDZ H2O computation

    h2o = psi4.geometry("""
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """)

    psi4.set_options({' basis ': "cc-pVDZ"})
    psi4.energy('scf')

    assert psi4.compare_values(-76.0266327341067125, psi4.variable('SCF TOTAL ENERGY'), 6, 'SCF energy')
