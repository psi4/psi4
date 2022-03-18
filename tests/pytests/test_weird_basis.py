from collections import Counter

import pytest

import psi4

from utils import compare_values


pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, pytest.mark.scf]


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


_av10z_feller_H = """
spherical
****
H 0
S 19 1.0
 4.2749598000e+06  2.0850302000e-09
 2.8489434000e+05  1.2800000000e-07
 1.3950281000e+04  2.8880000000e-06
 3.4719458000e+03  5.7980000000e-06
 1.9025267000e+03  1.7364000000e-05
 6.2356338000e+02  7.2599000000e-05
 2.9036861000e+02  1.1927900000e-04
 1.2650425000e+02  5.2892400000e-04
 4.4545793000e+01  2.0203280000e-03
 1.5983776000e+01  6.8656000000e-03
 6.0889730000e+00  2.0853603000e-02
 2.4835050000e+00  5.5989638000e-02
 1.0775260000e+00  1.3009211500e-01
 4.9412500000e-01  2.4488996600e-01
 2.3815700000e-01  3.3268970900e-01
 1.1919700000e-01  2.6640337200e-01
 6.0125000000e-02  8.2526348000e-02
 2.2060000000e-02  2.2063960000e-03
 1.1706000000e-02 -2.8636900000e-04
S 1 1.0
 6.0889730000e+00  1.0000000000e+00
S 1 1.0
 2.4835050000e+00  1.0000000000e+00
S 1 1.0
 1.0775260000e+00  1.0000000000e+00
S 1 1.0
 4.9412500000e-01  1.0000000000e+00
S 1 1.0
 2.3815700000e-01  1.0000000000e+00
S 1 1.0
 1.1919700000e-01  1.0000000000e+00
S 1 1.0
 6.0125000000e-02  1.0000000000e+00
S 1 1.0
 2.2060000000e-02  1.0000000000e+00
S 1 1.0
 1.1706000000e-02  1.0000000000e+00
S 1 1.0
 7.7880000000e-03  1.0000000000e+00
P 1 1.0
 3.0150122000e+01  1.0000000000e+00
P 1 1.0
 1.5390349000e+01  1.0000000000e+00
P 1 1.0
 7.8561158000e+00  1.0000000000e+00
P 1 1.0
 4.0102115000e+00  1.0000000000e+00
P 1 1.0
 2.0470417000e+00  1.0000000000e+00
P 1 1.0
 1.0449273000e+00  1.0000000000e+00
P 1 1.0
 5.3339080000e-01  1.0000000000e+00
P 1 1.0
 2.7227320000e-01  1.0000000000e+00
P 1 1.0
 1.3898380000e-01  1.0000000000e+00
P 1 1.0
 4.9224000000e-02  1.0000000000e+00
D 1 1.0
 3.4729594000e+01  1.0000000000e+00
D 1 1.0
 1.7575838000e+01  1.0000000000e+00
D 1 1.0
 8.8947218000e+00  1.0000000000e+00
D 1 1.0
 4.5014113000e+00  1.0000000000e+00
D 1 1.0
 2.2780593000e+00  1.0000000000e+00
D 1 1.0
 1.1528727000e+00  1.0000000000e+00
D 1 1.0
 5.8344200000e-01  1.0000000000e+00
D 1 1.0
 2.9526640000e-01  1.0000000000e+00
D 1 1.0
 1.0643900000e-01  1.0000000000e+00
F 1 1.0
 1.6883840000e+01  1.0000000000e+00
F 1 1.0
 9.0325589000e+00  1.0000000000e+00
F 1 1.0
 4.8322610000e+00  1.0000000000e+00
F 1 1.0
 2.5851751000e+00  1.0000000000e+00
F 1 1.0
 1.3830235000e+00  1.0000000000e+00
F 1 1.0
 7.3989340000e-01  1.0000000000e+00
F 1 1.0
 3.9583000000e-01  1.0000000000e+00
F 1 1.0
 1.5436000000e-01  1.0000000000e+00
G 1 1.0
 9.1827775000e+00  1.0000000000e+00
G 1 1.0
 5.3126124000e+00  1.0000000000e+00
G 1 1.0
 3.0735636000e+00  1.0000000000e+00
G 1 1.0
 1.7781823000e+00  1.0000000000e+00
G 1 1.0
 1.0287512000e+00  1.0000000000e+00
G 1 1.0
 5.9517470000e-01  1.0000000000e+00
G 1 1.0
 2.3193300000e-01  1.0000000000e+00
H 1 1.0
 9.6155582000e+00  1.0000000000e+00
H 1 1.0
 5.6307069000e+00  1.0000000000e+00
H 1 1.0
 3.2972460000e+00  1.0000000000e+00
H 1 1.0
 1.9308110000e+00  1.0000000000e+00
H 1 1.0
 1.1306500000e+00  1.0000000000e+00
H 1 1.0
 4.0090300000e-01  1.0000000000e+00
I 1 1.0
 6.8540732000e+00  1.0000000000e+00
I 1 1.0
 3.6114556000e+00  1.0000000000e+00
I 1 1.0
 1.9028994000e+00  1.0000000000e+00
I 1 1.0
 1.0026500000e+00  1.0000000000e+00
I 1 1.0
 4.3501600000e-01  1.0000000000e+00
{L7} 1 1.0
 7.3870330000e+00  1.0000000000e+00
{L7} 1 1.0
 3.7343330000e+00  1.0000000000e+00
{L7} 1 1.0
 1.8878000000e+00  1.0000000000e+00
{L7} 1 1.0
 7.1571800000e-01  1.0000000000e+00
{L8} 1 1.0
 4.5331500000e+00  1.0000000000e+00
{L8} 1 1.0
 2.2155200000e+00  1.0000000000e+00
{L8} 1 1.0
 1.0367000000e+00  1.0000000000e+00
{L9} 1 1.0
 4.7209680000e+00  1.0000000000e+00
{L9} 1 1.0
 1.3433090000e+00  1.0000000000e+00
****
"""


_hijhik = pytest.mark.xfail(reason="HIJ/HIK convention uncertain", raises=psi4.driver.qcdb.exceptions.ValidationError, strict=True)


@pytest.mark.parametrize(
    "subject,bas",
    [
        pytest.param("H", _av10z_feller_H.format(L7="L=7", L8="L=8", L9="L=9"), id="HI7_num"),
        pytest.param("H", _av10z_feller_H.format(L7="K",   L8="L",   L9="L=9"), id="HIK_num"),
        pytest.param("H", _av10z_feller_H.format(L7="J",   L8="K",   L9="L=9"), id="HIJ_num"),
        pytest.param("H", _av10z_feller_H.format(L7="L=7", L8="L=8", L9="M"  ), id="HI7_charM", marks=_hijhik),
        pytest.param("H", _av10z_feller_H.format(L7="L=7", L8="L=8", L9="L"  ), id="HI7_charL", marks=_hijhik),
        pytest.param("H", _av10z_feller_H.format(L7="K",   L8="L",   L9="M"  ), id="HIK_char"),
        pytest.param("H", _av10z_feller_H.format(L7="J",   L8="K",   L9="L"  ), id="HIJ_char"),
    ],
)
def test_high_angmom_basis(subject, bas, mols, request):
    """accommodate different high-angmom conventions"""

    # 11s / 10p / 9d / 8f / 7g / 6h / 5i / 4k / 3l / 2m
    #  11 +  30 + 45 + 56 + 63 + 66 + 65 + 60 + 51 + 38 = 485

    mol = mols[subject]

    def basisspec_psi4_yo__anonymous1234(mol, role):
        mol.set_basis_all_atoms("test", role=role)
        return {"test": bas}

    psi4.driver.qcdb.libmintsbasisset.basishorde["FELLER"] = basisspec_psi4_yo__anonymous1234

    wert = psi4.core.BasisSet.build(mol, "BASIS", "FELLER")

    assert wert.nbf() == 485
    assert wert.max_am() == 9

    # Psi uses HI_K for export
    count = Counter(wert.shell(ish).AMCHAR for ish in range(wert.nshell()))
    assert dict(count) == {'S': 11, 'P': 10, 'D': 9, 'F': 8, 'G': 7, 'H': 6, 'I': 5, 'K': 4, 'L': 3, 'M': 2}


@pytest.mark.parametrize(
    "subject,bas",
    [
        pytest.param("H", "PSDG"),
        pytest.param("H", "PSDGIF"),
        pytest.param("H", "PSDGKIF", marks=_hijhik),
        pytest.param("H", "PSDGLIF", marks=_hijhik),
    ],
)
def test_skipped_angmom(subject, bas, mols, request):

    def fake_basis_string(seq):
        amchar = "SPDFGHIKLM"

        nbf = 0
        basis_string = ["spherical", "****", "H 0"]

        for am in seq:
            basis_string.append(f"{am} 1 1.0")
            basis_string.append(" 2.2 1.0")
            l = amchar.index(am)
            nbf += 2 * l + 1

        basis_string.append("****")

        return nbf, "\n".join(basis_string)

    nbf, basis_string = fake_basis_string(bas)

    mol = mols[subject]

    def basisspec_psi4_yo__anonymous1234(mol, role):
        mol.set_basis_all_atoms("test", role=role)
        return {"test": basis_string}

    psi4.driver.qcdb.libmintsbasisset.basishorde["FAKE"] = basisspec_psi4_yo__anonymous1234

    wert = psi4.core.BasisSet.build(mol, "BASIS", "FAKE")

    assert wert.nbf() == nbf
