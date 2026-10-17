"""Basis construction parses each (source, entry) once and keeps per-atom shells independent."""

from collections import Counter

import pytest

import psi4
from psi4.driver import qcdb
from psi4.driver.qcdb.libmintsbasisset import basishorde
from psi4.driver.qcdb.libmintsbasissetparser import Gaussian94BasisSetParser

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def _count_parses(monkeypatch):
    calls = Counter()
    original = Gaussian94BasisSetParser.parse

    def parse(self, entry, lines):
        calls[(entry, lines if isinstance(lines, str) else tuple(lines))] += 1
        return original(self, entry, lines)

    monkeypatch.setattr(Gaussian94BasisSetParser, "parse", parse)
    return calls


def _center_shells(basis, center):
    return [(sh.am(), list(sh.PYexp), list(sh.PYoriginal_coef), sh.rpowers)
            for sh in (basis.shell(i, center) for i in range(basis.nshell_on_center(center)))]


def _construct(geom, assignments, **kwargs):
    mol = qcdb.Molecule(geom + "\nsymmetry c1\nno_reorient\nno_com")
    for atom, name in enumerate(assignments):
        mol.set_basis_by_number(atom, name, role="BASIS")
    return qcdb.BasisSet.construct(Gaussian94BasisSetParser(), mol, "BASIS", **kwargs)


@pytest.mark.parametrize("atomlist", [False, True])
def test_parse_once_per_source_entry(monkeypatch, atomlist):
    calls = _count_parses(monkeypatch)
    mol = psi4.geometry("H 0 0 0\nH 0 0 1\nH 0 0 2\nH 0 0 3\nsymmetry c1")
    basis = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pvdz", puream=True, return_atomlist=atomlist)
    nbf = sum(b.nbf() for b in basis) if atomlist else basis.nbf()
    assert nbf == 20
    assert calls
    assert max(calls.values()) == 1


def test_label_ghost_and_redefined_source(monkeypatch):
    template = """spherical
****
H 0
S 1 1.0
{ordinary} 1.0
****
H_SPECIAL 0
S 1 1.0
2.0 1.0
****
"""
    source = {"basis": template.format(ordinary=1.0)}

    def spec(mol, role):
        mol.set_basis_all_atoms("basis", role=role)
        return source

    monkeypatch.setitem(basishorde, "PARSE_REUSE", spec)
    mol = psi4.geometry("0 2\nH_SPECIAL 0 0 0\nH 0 0 1\nH 0 0 2\n@H 0 0 3\nsymmetry c1")
    for ordinary in (1.0, 3.0):
        source["basis"] = template.format(ordinary=ordinary)
        basis = psi4.core.BasisSet.build(mol, "ORBITAL", "PARSE_REUSE", puream=True)
        assert basis.nbf() == 4
        assert [basis.shell(i).exp(0) for i in range(4)] == [2.0, ordinary, ordinary, ordinary]


def test_parse_reuse_respects_puream(monkeypatch):

    def spec(mol, role):
        mol.set_basis_all_atoms("basis", role=role)
        return {"basis": "spherical\n****\nH 0\nD 1 1.0\n1.0 1.0\n****"}

    monkeypatch.setitem(basishorde, "PARSE_PUREAM", spec)
    mol = psi4.geometry("H 0 0 0\nH 0 0 1\nsymmetry c1")
    assert psi4.core.BasisSet.build(mol, "ORBITAL", "PARSE_PUREAM", puream=True).nbf() == 10
    assert psi4.core.BasisSet.build(mol, "ORBITAL", "PARSE_PUREAM", puream=False).nbf() == 12


@pytest.mark.parametrize("forced_puream,expected", [(True, 10), (False, 12)])
def test_parser_configuration(forced_puream, expected):
    mol = qcdb.Molecule("H 0 0 0\nH 0 0 1\nsymmetry c1")
    mol.set_basis_all_atoms("inline", role="BASIS")
    basis, _, _ = qcdb.BasisSet.construct(
        Gaussian94BasisSetParser(forced_puream=forced_puream),
        mol,
        "BASIS",
        basstrings={"inline": "spherical\n****\nH 0\nD 1 1.0\n1.0 1.0\n****"},
    )
    assert basis.nbf() == expected


def test_two_sources_for_same_element(monkeypatch):

    def spec(mol, role):
        for atom, name in enumerate(("first", "second", "first", "second")):
            mol.set_basis_by_number(atom, name, role=role)
        return {
            name: f"spherical\n****\nH 0\nS 1 1.0\n{exponent} 1.0\n****"
            for name, exponent in (("first", 1.0), ("second", 3.0))
        }

    monkeypatch.setitem(basishorde, "TWO_SOURCES", spec)
    calls = _count_parses(monkeypatch)
    mol = psi4.geometry("H1 0 0 0\nH2 0 0 1\nH3 0 0 2\nH4 0 0 3\nsymmetry c1")
    basis = psi4.core.BasisSet.build(mol, "ORBITAL", "TWO_SOURCES", puream=True)
    assert [basis.shell(i).exp(0) for i in range(4)] == [1.0, 3.0, 1.0, 3.0]
    assert max(calls.values()) == 1


def test_decontracted_and_contracted_share_source(monkeypatch):
    geom = "O1 0 0 0\nO2 0 0 2\nO3 0 0 4"
    ref_con, _, _ = _construct("O 0 0 0", ["cc-pvdz"])
    ref_dec, _, _ = _construct("O 0 0 0", ["cc-pvdz-decon"])
    assert ref_dec.nbf() > ref_con.nbf()

    calls = _count_parses(monkeypatch)
    for order in (["cc-pvdz-decon", "cc-pvdz", "cc-pvdz-decon"], ["cc-pvdz", "cc-pvdz-decon", "cc-pvdz"]):
        calls.clear()
        basis, _, _ = _construct(geom, order)
        assert max(calls.values()) == 1
        for center, name in enumerate(order):
            ref = ref_dec if name.endswith("-decon") else ref_con
            assert _center_shells(basis, center) == _center_shells(ref, 0)


def test_mixed_library_bases_for_one_element(monkeypatch):
    names = ["cc-pvtz", "cc-pvdz", "6-31g*", "cc-pvdz", "cc-pvdz-decon", "cc-pvtz"]
    geom = "\n".join(f"C{i + 1} 0 0 {1.5 * i}" for i in range(len(names)))
    refs = {name: _construct("C 0 0 0", [name])[0] for name in set(names)}

    calls = _count_parses(monkeypatch)
    basis, _, _ = _construct(geom, names)
    assert max(calls.values()) == 1
    # cc-pvdz and cc-pvdz-decon read one file, so "C" is parsed once per distinct file
    assert sum(n for (entry, _), n in calls.items() if entry == "C") == 3
    for center, name in enumerate(names):
        assert _center_shells(basis, center) == _center_shells(refs[name], 0)


def test_repeated_ecp_atoms():
    _, _, ref = _construct("I 0 0 0", ["def2-svp"])
    _, _, ecp = _construct("I1 0 0 0\nI2 0 0 3\nI3 0 0 6", ["def2-svp"] * 3)
    assert ecp.ecp_coreinfo == {"I1": 28, "I2": 28, "I3": 28}
    expected = _center_shells(ref, 0)
    assert expected
    for center in range(3):
        assert _center_shells(ecp, center) == expected



def test_atoms_do_not_share_parsed_shells():
    """Each label keeps its own ShellInfo objects in atom_basis_shell, as before parse reuse."""
    basis, _, ecp = _construct("I1 0 0 0\nI2 0 0 3\nI3 0 0 6", ["def2-svp"] * 3)
    for bs in (basis, ecp):
        per_label = [bs.atom_basis_shell[label]["def2-svp"] for label in ("I1", "I2", "I3")]
        shells = [sh for label_shells in per_label for sh in label_shells]
        assert shells
        assert len({id(lst) for lst in per_label}) == 3
        assert len({id(sh) for sh in shells}) == len(shells)
        assert len({id(sh.PYexp) for sh in shells}) == len(shells)
