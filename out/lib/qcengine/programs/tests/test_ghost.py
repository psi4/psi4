import re
import pprint
import pytest
import qcelemental as qcel
from qcelemental.testing import compare, compare_values
from qcengine.programs.tests.test_dftd3_mp2d import eneyne_ne_qcschemamols

import qcengine as qcng
from qcengine.testing import using


@pytest.fixture
def hene():
    smol = """
 0 1
 He 0 0 0
 @Ne 2.5 0 0
"""
    return qcel.models.Molecule.from_data(smol)


@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param("cfour", "aug-pvdz", {}, marks=using("cfour")),
        pytest.param("gamess", "ccd", {"contrl__ispher": 1}, marks=using("gamess")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True}, marks=using("nwchem")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True, "qc_module": "tce"}, marks=using("nwchem")),
        pytest.param("psi4", "aug-cc-pvdz", {"scf_type": "direct"}, marks=using("psi4")),
    ],
)
def test_simple_ghost(program, basis, keywords, hene):
    resi = {"molecule": hene, "driver": "energy", "model": {"method": "hf", "basis": basis}, "keywords": keywords}

    if program == "gamess":
        with pytest.raises(qcng.exceptions.InputError) as e:
            res = qcng.compute(resi, program, raise_error=True, return_dict=True)
        pytest.xfail("no ghosts with gamess")

    res = qcng.compute(resi, program, raise_error=True, return_dict=True)

    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    pprint.pprint(res, width=200)

    atol = 1.0e-6
    assert compare_values(0.0, res["properties"]["nuclear_repulsion_energy"], atol=atol, label="nre")
    assert compare(32, res["properties"]["calcinfo_nbasis"], label="nbas")
    assert compare(32, res["properties"]["calcinfo_nmo"], label="nmo")
    assert compare_values(-2.8557143339397539, res["return_result"], atol=atol, label="ene")


bimol_ref = {}
dmm = ["dimer", "mA", "mB", "mAgB", "gAmB"]
bimol_ref["eneyne"] = {}
bimol_ref["eneyne"]["natom"] = dict(zip(dmm, [10, 6, 4, 10, 10]))
bimol_ref["eneyne"]["nreal"] = dict(zip(dmm, [10, 6, 4, 6, 4]))
bimol_ref["eneyne"]["nre"] = dict(zip(dmm, [85.1890645313, 33.3580722134, 24.6979461998, 33.3580722134, 24.6979461998]))
bimol_ref["eneyne"]["pg"] = dict(zip(dmm, ["C2v", ["D2h", "C2v"], ["D4h", "C4v", "C2v"], "C2v", "C2v"]))
# 6-31G*
bimol_ref["eneyne"]["nbasis"] = dict(zip(dmm, [72, 38, 34, 72, 72]))
bimol_ref["eneyne"]["nmo"] = dict(zip(dmm, [72, 38, 34, 72, 72]))
bimol_ref["eneyne"]["mp2"] = dict(
    zip(dmm, [-155.3738614073, -78.2942587466, -77.0760547257, -78.2957592762, -77.0762583352])
)


@pytest.mark.parametrize("subject", dmm)
@pytest.mark.parametrize(
    "qcprog, basis, keywords",
    [
        pytest.param("cfour", "6-31g*", {"spherical": 0, "scf_conv": 12}, id="cfour", marks=using("cfour")),
        pytest.param(
            "gamess",
            "n31",
            {"basis__ngauss": 6, "basis__ndfunc": 1, "mp2__nacore": 0},
            id="gamess",
            marks=using("gamess"),
        ),
        pytest.param("nwchem", "6-31g*", {"scf__thresh": 1.0e-8}, id="nwchem", marks=using("nwchem")),
        pytest.param(
            "psi4",
            "6-31g*",
            {
                "scf_type": "pk",
                "mp2_type": "conv",
            },
            id="psi4",
            marks=using("psi4"),
        ),
    ],
)
def test_tricky_ghost(qcprog, subject, basis, keywords):
    kmol = qcel.models.Molecule(**eneyne_ne_qcschemamols()["eneyne"][subject])
    ref = bimol_ref["eneyne"]

    assert len(kmol.symbols) == ref["natom"][subject]
    assert sum([int(at) for at in kmol.real]) == ref["nreal"][subject]

    atin = qcel.models.AtomicInput(
        **{"molecule": kmol, "model": {"method": "mp2", "basis": basis}, "driver": "energy", "keywords": keywords}
    )

    if qcprog == "gamess" and subject in ["mAgB", "gAmB"]:
        with pytest.raises(qcng.exceptions.InputError) as e:
            res = qcng.compute(atin, qcprog, raise_error=True)
        pytest.xfail("no ghosts with gamess")

    atres = qcng.compute(atin, qcprog)
    pprint.pprint(atres.dict(), width=200)

    assert compare_values(
        ref["nre"][subject], atres.properties.nuclear_repulsion_energy, atol=1.0e-4, label="nre"
    ), f'nre: {atres.properties.nuclear_repulsion_energy} != {ref["nre"][subject]}'
    assert compare(
        ref["nbasis"][subject], atres.properties.calcinfo_nbasis, label="nbasis"
    ), f'nbasis: {atres.properties.calcinfo_nbasis} != {ref["nbasis"][subject]}'
    assert compare(
        ref["nmo"][subject], atres.properties.calcinfo_nmo, label="nmo"
    ), f'nmo: {atres.properties.calcinfo_nmo} != {ref["nmo"][subject]}'
    assert compare_values(
        ref["mp2"][subject], atres.return_result, atol=3.0e-6, label="ene"
    ), f'ene: {atres.return_result} != {ref["mp2"][subject]}'

    pgline = {
        "cfour": r"Computational point group: (?P<pg>\w+)",
        "gamess": r"THE POINT GROUP IS (?P<pg>[\w\s,=]+)",
        "nwchem": r"Group name\s+(?P<pg>\w+)",
        "psi4": r"Running in (?P<pg>\w+) symmetry.",
    }
    mobj = re.search(pgline[qcprog], atres.stdout)
    if mobj:
        pg = mobj.group("pg").strip()
        if pg == "CNV, NAXIS= 2, ORDER= 4":
            pg = "C2v"
        elif pg == "C1 , NAXIS= 0, ORDER= 1":
            pg = "C1"
        pg = pg.capitalize()

    if qcprog == "gamess":  # and subject in ["mAgB", "gAmB"]:
        # don't know how to get master frame w/ghosts in gamess, so C1 forced
        assert pg == "C1", f"pg: {pg} != C1"
    else:
        assert pg in ref["pg"][subject], f'pg: {pg} != {ref["pg"][subject]}'


@pytest.mark.parametrize(
    "qcprog, basis, keywords",
    [
        pytest.param("cfour", "aug-pvdz", {"scf_conv": 12}, id="cfour", marks=using("cfour")),
        pytest.param(
            "gamess",
            "accd",
            {"mp2__nacore": 0, "contrl__ispher": 1},
            id="gamess",
            marks=using("gamess"),
        ),
        pytest.param(
            "nwchem", "aug-cc-pvdz", {"scf__nr": 1.0, "scf__thresh": 1.0e-8}, id="nwchem", marks=using("nwchem")
        ),
        pytest.param(
            "psi4",
            "aug-cc-pvdz",
            {
                "scf_type": "pk",
                "mp2_type": "conv",
            },
            id="psi4",
            marks=using("psi4"),
        ),
    ],
)
def test_atom_labels(qcprog, basis, keywords):
    kmol = qcel.models.Molecule.from_data(
        """
      H       0 0 0
      H5      5 0 0
      H_other 0 5 0
      H_4sq   5 5 0
      units au
    """
    )

    assert compare(["H", "H", "H", "H"], kmol.symbols, "elem")
    assert compare(["", "5", "_other", "_4sq"], kmol.atom_labels, "elbl")

    atin = qcel.models.AtomicInput(
        **{"molecule": kmol, "model": {"method": "mp2", "basis": basis}, "driver": "energy", "keywords": keywords}
    )

    atres = qcng.compute(atin, qcprog)
    pprint.pprint(atres.dict(), width=200)

    nre = 1.0828427
    assert compare_values(
        nre, atres.properties.nuclear_repulsion_energy, atol=1.0e-4, label="nre"
    ), f"nre: {atres.properties.nuclear_repulsion_energy} != {nre}"

    nmo = 36
    assert compare(nmo, atres.properties.calcinfo_nmo, label="nmo"), f"nmo: {atres.properties.calcinfo_nmo} != {nmo}"

    scf = -1.656138508
    assert compare_values(
        scf, atres.properties.scf_total_energy, atol=3.0e-6, label="scf ene"
    ), f"scf ene: {atres.properties.scf_total_energy} != {scf}"

    mp2 = -1.7926264513
    assert compare_values(mp2, atres.return_result, atol=3.0e-6, label="ene"), f"ene: {atres.return_result} != {mp2}"
