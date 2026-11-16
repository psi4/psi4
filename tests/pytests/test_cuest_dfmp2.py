"""End-to-end check of the cuEST DF-MP2 path (cuEST SCF + GPU MP2 energy loop)
against plain CPU Psi4 (CPU SCF + CPU DF-MP2).
"""
import psi4
import pytest

from utils import *
from addons import uusing

pytestmark = [pytest.mark.psi, pytest.mark.api]

_GEOM = """
0 1
O  -1.551007  -0.114520   0.000000
H  -1.934259   0.762503   0.000000
H  -0.599677   0.040712   0.000000
O   1.350625   0.111469   0.000000
H   1.680398  -0.373741  -0.758561
H   1.680398  -0.373741   0.758561
units angstrom
symmetry c1
no_reorient
no_com
"""

_ATOL = 1e-7


def _run_mp2(use_cuest, freeze_core):
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.core.clean_variables()
    psi4.core.set_num_threads(4)
    psi4.geometry(_GEOM)
    psi4.set_options({
        "basis": "def2-svp",
        "df_basis_scf": "def2-universal-jkfit",
        "df_basis_mp2": "def2-svp-ri",
        "scf_type": "df",
        "mp2_type": "df",
        "qc_module": "dfmp2",
        "reference": "rhf",
        "freeze_core": freeze_core,
        "e_convergence": 1e-10,
        "d_convergence": 1e-10,
        "use_cuest": use_cuest,
        "cuest_mixed_precision": False,
    })
    total = psi4.energy("mp2")
    return {
        "scf": psi4.variable("SCF TOTAL ENERGY"),
        "ss": psi4.variable("MP2 SAME-SPIN CORRELATION ENERGY"),
        "os": psi4.variable("MP2 OPPOSITE-SPIN CORRELATION ENERGY"),
        "total": total,
    }


@uusing("cuest")
@uusing("cuda_cc8")
@pytest.mark.quick
@pytest.mark.parametrize("freeze_core", [False, True], ids=["all_electron", "frozen_core"])
def test_cuest_dfmp2_matches_cpu(freeze_core):
    cpu = _run_mp2(use_cuest=False, freeze_core=freeze_core)
    gpu = _run_mp2(use_cuest=True, freeze_core=freeze_core)

    # Sanity: a real MP2 calculation happened.
    assert cpu["os"] < -1e-2 and cpu["ss"] < -1e-3

    assert compare_values(cpu["scf"], gpu["scf"], _ATOL, "SCF energy: cuEST vs CPU")
    assert compare_values(cpu["ss"], gpu["ss"], _ATOL, "MP2 SS: cuEST vs CPU")
    assert compare_values(cpu["os"], gpu["os"], _ATOL, "MP2 OS: cuEST vs CPU")
    assert compare_values(cpu["total"], gpu["total"], _ATOL, "MP2 total: cuEST vs CPU")
