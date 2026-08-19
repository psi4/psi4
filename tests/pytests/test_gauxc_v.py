import pytest

import psi4
from addons import uusing

pytestmark = [pytest.mark.psi, pytest.mark.api]


@pytest.fixture
def water():
    return psi4.geometry("""
O 0.0  0.000 0.000
H 0.0  0.757 0.587
H 0.0 -0.757 0.587
units angstrom
symmetry c1
no_reorient
no_com
""")


def test_dft_v_algorithm_default_is_internal():
    """Der Default darf sich nicht aendern."""
    psi4.core.clean_options()
    assert psi4.core.get_option("SCF", "DFT_V_ALGORITHM") == "INTERNAL"


@uusing("gauxc")
def test_gauxc_v_selectable(water):
    """DFT_V_ALGORITHM=GAUXC muss akzeptiert werden und den GauXC-Pfad waehlen."""
    psi4.core.clean_options()
    psi4.set_options({"basis": "cc-pVDZ", "scf_type": "df",
                      "dft_v_algorithm": "gauxc", "reference": "rks"})
    with pytest.raises(Exception, match="not implemented"):
        psi4.energy("svwn", molecule=water)
