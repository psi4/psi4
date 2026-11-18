import pytest

import numpy as np
import psi4

from utils import compare_arrays

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]

def test_fock_subset_mo():
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
    """)

    rhf_e, wfn = psi4.energy('SCF/cc-pVDZ', molecule=h2o, return_wfn=True)

    F_diagonals = []
    for h in wfn.epsilon_a().nph:
        F_diagonals.append(np.diag(h))
    F_expected = psi4.core.Matrix.from_array(F_diagonals)
    assert psi4.compare_matrices(F_expected, wfn.Fa_subset("MO"), 8, "Alpha Fock Matrix")
    assert psi4.compare_matrices(F_expected, wfn.Fb_subset("MO"), 8, "Beta Fock Matrix")


def test_module_roles():
    """set_module/module are the "qc_module" entry of the role map, and other roles ride along."""
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
    """)

    _, wfn = psi4.energy('SCF/cc-pVDZ', molecule=h2o, return_wfn=True)

    # set_module/module are a view on one role, so the two agree and neither loses the other
    assert wfn.module() == "scf"
    assert wfn.module_roles()["qc_module"] == "scf"
    assert wfn.module_role("qc_module") == "scf"

    # an unfilled role is empty rather than a lookup error
    assert wfn.module_role("no_such_role") == ""

    wfn.set_module_role("widget_maker", "acme")
    assert wfn.module_role("widget_maker") == "acme"
    assert wfn.module() == "scf"

    wfn.set_module("dfmp2")
    assert wfn.module() == "dfmp2"
    assert wfn.module_role("widget_maker") == "acme"

    # roles survive a file round trip, which carries them flattened under a prefix
    reloaded = psi4.core.Wavefunction.from_file(wfn.to_file())
    assert reloaded.module() == "dfmp2"
    assert reloaded.module_role("widget_maker") == "acme"


def test_merge_module_roles():
    """Roles from several results (findif displacements, CBS stages) consolidate per role:
    agreement keeps the module, disagreement is "(mixed)", and a result without a role doesn't vote."""
    from qcelemental.models.v2 import Provenance
    from psi4.driver.p4util import merge_module_roles, provenance_stamp, set_wfn_module_roles

    provs = [
        {"module_roles": {"qc_module": "scf", "orbital_optimizer": "internal", "coulomb_eri_integrals": "libint2"}},
        Provenance(
            creator="Psi4",
            version="1.0",
            routine="x",
            module="scf",
            module_roles={
                "qc_module": "scf",
                "orbital_optimizer": "openorbitaloptimizer",
                "second_order_orbital_optimizer": "opentrustregion",
                "coulomb_eri_integrals": "libint2",
            },
        ),
        {"module": "scf"},  # no module_roles dict; contributes only qc_module via module
    ]
    merged = merge_module_roles(provs)
    assert merged == {
        "coulomb_eri_integrals": "libint2",
        "orbital_optimizer": "(mixed)",
        "qc_module": "scf",
        "second_order_orbital_optimizer": "opentrustregion",
    }
    assert merge_module_roles([]) == {}

    # provenance_stamp accepts a module string or a module_roles dict (copying qc_module to module)
    assert "module_roles" not in provenance_stamp("x", module="scf")
    assert provenance_stamp("x", module="scf")["module"] == "scf"
    stamped = provenance_stamp("x", module=merged)
    assert stamped["module"] == "scf"
    assert stamped["module_roles"] == merged

    # onto a wfn: set_wfn_module_roles copies all roles, including qc_module (wfn.module())
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
    """)
    wfn = psi4.core.Wavefunction.build(h2o, "sto-3g")
    set_wfn_module_roles(wfn, Provenance(**stamped))
    assert wfn.module() == "scf"
    assert wfn.module_role("coulomb_eri_integrals") == "libint2"
    assert wfn.module_role("orbital_optimizer") == "(mixed)"
    assert wfn.module_role("second_order_orbital_optimizer") == "opentrustregion"

    # round-trip through QCPortal's Provenance and QCSchema v1 -> v2 AtomicResult conversion
    try:
        from qcportal.record_models import Provenance as QCPProvenance
        from qcportal.qcschema_v1 import AtomicResult as QCPAtomicResult, AtomicResultProperties
    except ImportError:
        pass
    else:
        qcp_res = QCPAtomicResult(
            driver="energy",
            model={"method": "hf", "basis": "sto-3g"},
            molecule={"symbols": ["H"], "geometry": [0.0, 0.0, 0.0], "molecular_multiplicity": 2},
            properties=AtomicResultProperties(),
            return_result=-0.5,
            provenance=QCPProvenance(**stamped).model_dump(),
            success=True,
        )
        v2_res = qcp_res.convert_v(2)
        assert merge_module_roles([v2_res.provenance]) == merged

