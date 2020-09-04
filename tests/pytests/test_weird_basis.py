import pytest

import psi4

from .utils import compare_values


pytestmark = [pytest.mark.quick, pytest.mark.scf]


@pytest.fixture
def mols():
    smols = {
        "H": """
0 2
H
""",
        "H2": """
1 2
H
H 1 1.0
""",
        "HeH": """
1 1
He
H 1 1.0
""",
    }

    return {k: psi4.core.Molecule.from_string(v) for k, v in smols.items()}


@pytest.mark.parametrize(
    "subject,bas,ans",
    [
        pytest.param(
            "H",
            """
spherical
****
H 0
S 1 1.0
 1.0 1.0
****
""",
            -0.09576912160573,
            id="01: H energy, S basis",
        ),
        pytest.param(
            "H",
            """
spherical
****
H 0
P 1 1.0
 1.0 1.0
****
""",
            1.43615391892951,
            id="02: H energy, P basis",
        ),
        pytest.param(
            "H",
            """
spherical
****
H 0
D 1 1.0
 1.0 1.0
****
""",
            2.64892313514361,
            id="03: H energy, D basis",
        ),
        pytest.param(
            "H",
            """
spherical
****
H 0
F 1 1.0
 1.0 1.0
****
""",
            3.77050554440881,
            id="04: H energy, F basis",
        ),
        pytest.param(
            "H2",
            """
spherical
****
H 0
S 1 1.0
 1.0 1.0
****
""",
            -0.33308154102469,
            id="05: H2+ energy, S basis",
        ),
        pytest.param(
            "H2",
            """
spherical
****
H 0
P 1 1.0
 1.0 1.0
****
""",
            1.17859086038869,
            id="06: H2+ energy, P basis",
        ),
        pytest.param(
            "H2",
            """
spherical
****
H 0
D 1 1.0
 1.0 1.0
****
""",
            2.13041388298912,
            id="07: H2+ energy, D basis",
        ),
        pytest.param(
            "H2",
            """
spherical
****
H 0
F 1 1.0
 1.0 1.0
****
""",
            3.15390046335277,
            id="08: H2+ energy, F basis",
        ),
        pytest.param(
            "HeH",
            """
spherical
****
H 0
S 1 1.0
 1.0 1.0
****
He 0
S 1 1.0
 1.0 1.0
****
""",
            -2.38847095431346,
            id="09: HeH+ energy, S basis",
        ),
        pytest.param(
            "HeH",
            """
spherical
****
H 0
P 1 1.0
 1.0 1.0
****
He 0
P 1 1.0
 1.0 1.0
****
""",
            1.36088614069560,
            id="10: HeH+ energy, P basis",
        ),
        pytest.param(
            "HeH",
            """
spherical
****
H 0
D 1 1.0
 1.0 1.0
****
He 0
D 1 1.0
 1.0 1.0
****
""",
            3.45306620376334,
            id="11: HeH+ energy, D basis",
        ),
        pytest.param(
            "HeH",
            """
spherical
****
H 0
F 1 1.0
 1.0 1.0
****
He 0
F 1 1.0
 1.0 1.0
****
""",
            5.65071627649024,
            id="12: HeH+ energy, F basis",
        ),
    ],
)
def test_weird_basis(subject, bas, ans, mols, request):
    """non-consecutive angular momentum basis set patterns"""

    mol = mols[subject]
    ref = "rhf" if mol.multiplicity() == 1 else "uhf"
    psi4.set_options(
        {"reference": ref, "guess": "core", "scf_type": "pk", "df_scf_guess": "false", "basis": "anonymous1234",}
    )

    def basisspec_psi4_yo__anonymous1234(mol, role):
        mol.set_basis_all_atoms("test", role=role)
        return {"test": bas}

    psi4.driver.qcdb.libmintsbasisset.basishorde["ANONYMOUS1234"] = basisspec_psi4_yo__anonymous1234

    ene = psi4.energy("scf", molecule=mol)
    assert compare_values(ans, ene, 6, request.node.name)

