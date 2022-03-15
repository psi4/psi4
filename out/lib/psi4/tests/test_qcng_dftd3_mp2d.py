import copy

import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.testing import compare, compare_recursive, compare_values, tnm

import qcengine as qcng
from qcengine.programs import empirical_dispersion_resources
from qcengine.testing import is_program_new_enough, using
from .addons import using

from qcengine.programs.tests import test_dftd3_mp2d
ref, gref = test_dftd3_mp2d.ref, test_dftd3_mp2d.gref


pytestmark = [pytest.mark.quick]

@using("dftd3")
@pytest.mark.parametrize("method", [
    "b3lyp-d3",
    "b3lyp-d3m",
    "b3lyp-d3bj",
    "b3lyp-d3mbj",
])
def test_dftd3_task(method):
    json_data = {"molecule": qcng.get_molecule("eneyne"), "driver": "energy", "model": {"method": method}}

    ret = qcng.compute(json_data, "dftd3", raise_error=True, return_dict=True)

    assert ret["driver"] == "energy"
    assert "provenance" in ret
    assert "normal termination of dftd3" in ret["stdout"]

    for key in ["cpu", "hostname", "username", "wall_time"]:
        assert key in ret["provenance"]

    assert ret["success"] is True



seneyne = """
C   0.000000  -0.667578  -2.124659
C   0.000000   0.667578  -2.124659
H   0.923621  -1.232253  -2.126185
H  -0.923621  -1.232253  -2.126185
H  -0.923621   1.232253  -2.126185
H   0.923621   1.232253  -2.126185
--
C   0.000000   0.000000   2.900503
C   0.000000   0.000000   1.693240
H   0.000000   0.000000   0.627352
H   0.000000   0.000000   3.963929
"""

sne = """
Ne 0 0 0
"""


def eneyne_ne_qcdbmols():
    if not is_program_new_enough("psi4", "1.4a1.dev55"):
        pytest.skip("Psi4 requires at least Psi4 v1.3rc2")
    from psi4.driver import qcdb

    eneyne = qcdb.Molecule(seneyne)
    ne = qcdb.Molecule(sne)
    mols = {
        'eneyne': {
            'dimer': eneyne,
            'mA': eneyne.extract_subsets(1),
            'mB': eneyne.extract_subsets(2),
            'mAgB': eneyne.extract_subsets(1, 2),
            'gAmB': eneyne.extract_subsets(2, 1),
        },
        'ne': {
            'atom': ne,
        }
    }
    return mols


def eneyne_ne_psi4mols():
    if not is_program_new_enough("psi4", "1.4a1.dev55"):
        pytest.skip("Psi4 requires at least Psi4 v1.3rc2")
    import psi4

    eneyne = psi4.core.Molecule.from_string(seneyne)
    ne = psi4.core.Molecule.from_string(sne)
    mols = {
        'eneyne': {
            'dimer': eneyne,
            'mA': eneyne.extract_subsets(1),
            'mB': eneyne.extract_subsets(2),
            'mAgB': eneyne.extract_subsets(1, 2),
            'gAmB': eneyne.extract_subsets(2, 1),
        },
        'ne': {
            'atom': ne,
        }
    }
    return mols


def eneyne_ne_qcschemamols():

    eneyne = qcel.molparse.to_schema(qcel.molparse.from_string(seneyne)['qm'], dtype=2)
    mA = qcel.molparse.to_schema(qcel.molparse.from_string('\n'.join(seneyne.splitlines()[:7]))['qm'], dtype=2)
    mB = qcel.molparse.to_schema(qcel.molparse.from_string('\n'.join(seneyne.splitlines()[-4:]))['qm'], dtype=2)
    ne = qcel.molparse.to_schema(qcel.molparse.from_string(sne)['qm'], dtype=2)

    mAgB = qcel.molparse.from_string(seneyne)['qm']
    mAgB['real'] = [(iat < mAgB['fragment_separators'][0])
                    for iat in range(len(mAgB['elem']))]  # works b/c chgmult doesn't need refiguring
    mAgB = qcel.molparse.to_schema(mAgB, dtype=2)

    gAmB = qcel.molparse.from_string(seneyne)['qm']
    gAmB['real'] = [(iat >= gAmB['fragment_separators'][0]) for iat in range(len(gAmB['elem']))]
    gAmB = qcel.molparse.to_schema(gAmB, dtype=2)

    mols = {
        'eneyne': {
            'dimer': eneyne,
            'mA': mA,
            'mB': mB,
            'mAgB': mAgB,
            'gAmB': gAmB,
        },
        'ne': {
            'atom': ne,
        }
    }
    return mols


db3lypd3bj = {
    'dashlevel': 'd3bj',
    'dashparams': {
        's8': 1.9889,
        's6': 1.0,
        'a2': 4.4211,
        'a1': 0.3981
    },
    'dashparams_citation': '',
    'fctldash': 'b3lyp-d3(bj)'
}
db3lypd3bjcustom = copy.deepcopy(db3lypd3bj)
db3lypd3bjcustom['fctldash'] = ''
db3lypd3bjcustom['dashparams']['a2'] = 5.4211

dpbed3zero = {
    'dashlevel': 'd3zero',
    'dashparams': {
        's6': 1.0,
        's8': 0.722,
        'sr6': 1.217,
        'sr8': 1.0,
        'alpha6': 14.0
    },
    'dashparams_citation': '',
    'fctldash': 'pbe-d3'
}

atmgr = {
    'dashlevel': 'atmgr',
    'dashparams': {
        'alpha6': 14.0,
    },
    'dashparams_citation': '',
    'fctldash': 'atm(gr)',
}

chg = {
    'dashlevel': 'chg',
    'dashparams': {
        's6': 1.0,
    },
    'dashparams_citation': '',
    'fctldash': 'chg',
}

dmp2dmp2 = {
    'dashlevel': 'dmp2',
    'dashparams': {
        's8': 1.187,
        'a1': 0.944,
        'a2': 0.480,
        'rcut': 0.72,
        'w': 0.20,
    },
    'dashparams_citation': '',
    'fctldash': 'mp2-dmp2'
}


def _compute_key(pjrec):
    return pjrec['fctldash'].upper()


## Tests


@pytest.mark.parametrize("inp,expected", [
    (({'name_hint': 'b3lyp', 'level_hint': 'd3bj'}, 'B3LYP-D3(BJ)'), db3lypd3bj),
    (({'name_hint': 'b3LYP', 'level_hint': 'D3bj'}, 'B3LYP-D3(BJ)'), db3lypd3bj),
    (({'param_tweaks': {'s8': 1.9889, 's6': 1.0, 'a2': 4.4211, 'a1': 0.3981}, 'level_hint': 'd3bj'}, 'B3LYP-D3(BJ)'), db3lypd3bj),
    (({'name_hint': 'b3lyp', 'level_hint': 'd3bJ', 'param_tweaks': {'a2': 4.4211}}, 'B3LYP-D3(BJ)'), db3lypd3bj),
    (({'verbose': 3, 'name_hint': 'b3lyp', 'level_hint': 'd3bJ', 'param_tweaks': {'a2': 5.4211}}, ''), db3lypd3bjcustom),
    (({'name_hint': 'b3lyp-d3bj', 'param_tweaks': {'a2': 4.4211}}, 'B3LYP-D3(BJ)'), db3lypd3bj),
    (({'name_hint': 'pbe', 'level_hint': 'd3zero'}, 'PBE-D3'), dpbed3zero),
    (({'name_hint': 'pbe', 'level_hint': 'd3'}, 'PBE-D3'), dpbed3zero),
    (({'name_hint': 'pbe-d3'}, 'PBE-D3'), dpbed3zero),
    (({'name_hint': 'atm(gr)', 'level_hint': 'atmgr'}, 'ATM(GR)'), atmgr),
    (({'name_hint': 'atmgr'}, 'ATM(GR)'), atmgr),
    (({'name_hint': 'bp86-atmgr'}, 'ATM(GR)'), atmgr),
    (({'name_hint': 'asdf-chg'}, 'CHG'), chg),
    (({'name_hint': 'mp2-dmp2'}, 'MP2-DMP2'), dmp2dmp2),
    (({'name_hint': 'MP2', 'level_hint': 'dmp2'}, 'MP2-DMP2'), dmp2dmp2),
])  # yapf: disable
def test_dftd3__from_arrays(inp, expected):
    res = empirical_dispersion_resources.from_arrays(**inp[0])
    assert compare_recursive(expected, res, atol=1.e-4)
    assert compare(inp[1], _compute_key(res), 'key')
    res = empirical_dispersion_resources.from_arrays(name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    assert compare_recursive(expected, res, tnm() + ' idempotent', atol=1.e-4)


@pytest.mark.parametrize("inp", [
    ({'name_hint': 'b3lyp', 'level_hint': 'd3bJ', 'param_tweaks': {'a3': 5.4211}}),
    ({'name_hint': 'fakeb3lyp', 'level_hint': 'd3bJ', 'param_tweaks': {'s6': 5.4211}}),
    ({'level_hint': 'd3bJ', 'param_tweaks': {'s6': 5.4211}}),
    ({'name_hint': 'b3lyp-d3bj', 'param_tweaks': {'a2': 4.4211, 'zzz': 0.0}}),
    ({'name_hint': 'asdf-d4'}),
    ({'name_hint': 'atm(gr)', 'level_hint': 'chg'}),
])  # yapf:disable
def test_dftd3__from_arrays__error(inp):
    with pytest.raises(qcng.exceptions.InputError):
        empirical_dispersion_resources.from_arrays(**inp)


def test_dftd3__from_arrays__supplement():
    ans = {
        'dashlevel': 'chg',
        'dashparams': {
            's6': 4.05
        },
        'fctldash': 'asdf-d4',
        'dashparams_citation': '    mypaper\n'
    }
    supp = {'chg': {'definitions': {'asdf-d4': {'params': {'s6': 4.05}, 'citation': '    mypaper\n'}}}}

    res = empirical_dispersion_resources.from_arrays(name_hint='asdf-d4', level_hint='chg', dashcoeff_supplement=supp)
    assert compare_recursive(ans, res, atol=1.e-4)
    with pytest.raises(qcng.exceptions.InputError) as e:
        empirical_dispersion_resources.from_arrays(name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    assert "Can't guess -D correction level" in str(e.value)
    res = empirical_dispersion_resources.from_arrays(
        name_hint=res['fctldash'],
        level_hint=res['dashlevel'],
        param_tweaks=res['dashparams'],
        dashcoeff_supplement=supp)
    assert compare_recursive(ans, res, tnm() + ' idempotent', atol=1.e-4)


@using("dftd3")
def test_3():
    sys = qcel.molparse.from_string(seneyne)['qm']

    resinp = {
        'schema_name': 'qcschema_input',
        'schema_version': 1,
        'molecule': qcel.molparse.to_schema(sys, dtype=2),
        'driver': 'energy',
        'model': {
            'method': 'b3lyp',
        },
        'keywords': {
            'level_hint': 'd3bj'
        },
    }
    res = qcng.compute(resinp, 'dftd3', raise_error=True)
    res = res.dict()

    #res = dftd3.run_dftd3_from_arrays(molrec=sys, name_hint='b3lyp', level_hint='d3bj')
    assert compare('B3LYP-D3(BJ)', _compute_key(res['extras']['local_keywords']), 'key')


@using("dftd3")
@pytest.mark.parametrize(
    "subjects",
    [
        pytest.param(eneyne_ne_psi4mols, marks=using("psi4")),
        pytest.param(eneyne_ne_qcdbmols,
                     marks=using("psi4")),  # needs qcdb.Molecule, presently more common in psi4 than in qcdb
    ],
    ids=['qmol', 'pmol'])
@pytest.mark.parametrize(
    "inp", [
        ({'first': 'b3lyp', 'second': 'd', 'parent': 'eneyne', 'subject': 'dimer', 'lbl': 'B3LYP-D2'}),
        ({'first': 'b3lyp', 'second': 'd3bj', 'parent': 'eneyne', 'subject': 'mA', 'lbl': 'B3LYP-D3(BJ)'}),
        ({'first': 'pbe', 'second': 'd3zero', 'parent': 'eneyne', 'subject': 'mB', 'lbl': 'PBE-D3'}),
        ({'first': 'pbe', 'second': 'd3zero', 'parent': 'eneyne', 'subject': 'gAmB', 'lbl': 'PBE-D3'}),
        ({'first': 'pbe', 'second': 'd2', 'parent': 'eneyne', 'subject': 'mAgB', 'lbl': 'PBE-D2'}),
        ({'first': 'b3lyp', 'second': 'd3bj', 'parent': 'ne', 'subject': 'atom', 'lbl': 'B3LYP-D3(BJ)'}),
        #({'first': '', 'second': 'atmgr', 'parent': 'eneyne', 'subject': 'dimer', 'lbl': 'ATM'}),
        #({'first': 'b3lyp', 'second': 'atmgr', 'parent': 'eneyne', 'subject': 'mA', 'lbl': 'ATM'}),
        #({'first': 'pbe', 'second': 'atm(gr)', 'parent': 'eneyne', 'subject': 'mB', 'lbl': 'ATM'}),
        #({'first': '', 'second': 'ATMgr', 'parent': 'eneyne', 'subject': 'mAgB', 'lbl': 'ATM'}),
        # below two xfail until dftd3 that's only 2-body is out of psi4 proper
        pytest.param({'first': 'atmgr', 'second': 'atmgr', 'parent': 'eneyne', 'subject': 'gAmB', 'lbl': 'ATM'}, marks=[using("dftd3_321"), pytest.mark.xfail]),
        pytest.param({'first': 'pbe-atmgr', 'second': None, 'parent': 'ne', 'subject': 'atom', 'lbl': 'ATM'}, marks=[using("dftd3_321"), pytest.mark.xfail]),
    ])  # yapf: disable
def test_molecule__run_dftd3__23body(inp, subjects):
    subject = subjects()[inp['parent']][inp['subject']]
    expected = ref[inp['parent']][inp['lbl']][inp['subject']]
    gexpected = gref[inp['parent']][inp['lbl']][inp['subject']]

    E, G = subject.run_dftd3(inp['first'], inp['second'])
    assert compare_values(expected, E, atol=1.e-7)
    assert compare_values(gexpected, G, atol=1.e-7)


@using("qcdb")
def test_qcdb__energy_d3():
    eneyne = qcdb.set_molecule(seneyne)
    eneyne.update_geometry()

    E, jrec = qcdb.energy('d3-b3lyp-d2', return_wfn=True)
    assert compare_values(ref['eneyne']['B3LYP-D2']['dimer'], E, 7, 'P: Ethene-Ethyne -D2')
    assert compare_values(ref['eneyne']['B3LYP-D2']['dimer'], jrec['qcvars']['DISPERSION CORRECTION ENERGY'].data, 7,
                          tnm())
    assert compare_values(ref['eneyne']['B3LYP-D2']['dimer'],
                          jrec['qcvars']['B3LYP-D2 DISPERSION CORRECTION ENERGY'].data, 7, tnm())

    mA = eneyne.extract_subsets(1)

    E, jrec = qcdb.energy('d3-b3lyp-d3bj', return_wfn=True, molecule=mA)
    assert compare_values(ref['eneyne']['B3LYP-D3(BJ)']['mA'], E, 7, tnm())
    assert compare_values(ref['eneyne']['B3LYP-D3(BJ)']['mA'], jrec['qcvars']['DISPERSION CORRECTION ENERGY'].data, 7,
                          tnm())
    assert compare_values(ref['eneyne']['B3LYP-D3(BJ)']['mA'],
                          jrec['qcvars']['B3LYP-D3(BJ) DISPERSION CORRECTION ENERGY'].data, 7, tnm())


@using("mp2d")
@pytest.mark.parametrize(
    "subjects",
    [
        pytest.param(eneyne_ne_psi4mols, marks=using("psi4")),
        pytest.param(eneyne_ne_qcdbmols,
                     marks=using("psi4")),  # needs qcdb.Molecule, presently more common in psi4 than in qcdb
        pytest.param(eneyne_ne_qcschemamols),
    ],
    ids=['qmol', 'pmol', 'qcmol'])
@pytest.mark.parametrize("inp", [
    ({'parent': 'eneyne', 'name': 'mp2d-mp2-dmp2', 'subject': 'dimer', 'lbl': 'MP2-DMP2'}),
    ({'parent': 'eneyne', 'name': 'mp2d-mp2-dmp2', 'subject': 'mA', 'lbl': 'MP2-DMP2'}),
    ({'parent': 'eneyne', 'name': 'mp2d-mp2-dmp2', 'subject': 'mB', 'lbl': 'MP2-DMP2'}),
    ({'parent': 'eneyne', 'name': 'mp2d-mp2-dmp2', 'subject': 'gAmB', 'lbl': 'MP2-DMP2'}),
    ({'parent': 'eneyne', 'name': 'mp2d-mp2-dmp2', 'subject': 'mAgB', 'lbl': 'MP2-DMP2'}),
    ({'parent': 'ne', 'name': 'mp2d-mp2-dmp2', 'subject': 'atom', 'lbl': 'MP2-DMP2'}),
])  # yapf: disable
def test_mp2d__run_mp2d__2body(inp, subjects, request):
    subject = subjects()[inp['parent']][inp['subject']]
    expected = ref[inp['parent']][inp['lbl']][inp['subject']]
    gexpected = gref[inp['parent']][inp['lbl']][inp['subject']]

    if 'qcmol' in request.node.name:
        mol = subject
    else:
        mol = subject.to_schema(dtype=2)

    resinp = {
        'schema_name': 'qcschema_input',
        'schema_version': 1,
        'molecule': mol,
        'driver': 'gradient',
        'model': {
            'method': inp['name']
        },
        'keywords': {},
    }
    jrec = qcng.compute(resinp, 'mp2d', raise_error=True)
    jrec = jrec.dict()

    #assert len(jrec['extras']['qcvars']) == 8

    assert compare_values(expected, jrec['extras']['qcvars']['CURRENT ENERGY'], atol=1.e-7)
    assert compare_values(expected, jrec['extras']['qcvars']['DISPERSION CORRECTION ENERGY'], atol=1.e-7)
    assert compare_values(expected, jrec['extras']['qcvars'][inp['lbl'] + ' DISPERSION CORRECTION ENERGY'], atol=1.e-7)

    assert compare_values(gexpected, jrec['extras']['qcvars']['CURRENT GRADIENT'], atol=1.e-7)
    assert compare_values(gexpected, jrec['extras']['qcvars']['DISPERSION CORRECTION GRADIENT'], atol=1.e-7)
    assert compare_values(
    gexpected, jrec['extras']['qcvars'][inp['lbl'] + ' DISPERSION CORRECTION GRADIENT'], atol=1.e-7)


@using("dftd3")
@pytest.mark.parametrize(
    "subjects",
    [
        pytest.param(eneyne_ne_psi4mols, marks=using("psi4")),
        pytest.param(eneyne_ne_qcdbmols,
                     marks=using("psi4")),  # needs qcdb.Molecule, presently more common in psi4 than in qcdb
        pytest.param(eneyne_ne_qcschemamols),
    ],
    ids=['qmol', 'pmol', 'qcmol'])
@pytest.mark.parametrize("inp", [
    ({'parent': 'eneyne', 'name': 'd3-b3lyp-d', 'subject': 'dimer', 'lbl': 'B3LYP-D2'}),
    ({'parent': 'eneyne', 'name': 'd3-b3lyp-d3bj', 'subject': 'mA', 'lbl': 'B3LYP-D3(BJ)'}),
    ({'parent': 'eneyne', 'name': 'd3-PBE-D3zero', 'subject': 'mB', 'lbl': 'PBE-D3'}),
    ({'parent': 'eneyne', 'name': 'd3-PBE-D3zero', 'subject': 'gAmB', 'lbl': 'PBE-D3'}),
    ({'parent': 'eneyne', 'name': 'd3-PBE-D2', 'subject': 'mAgB', 'lbl': 'PBE-D2'}),
    ({'parent': 'ne', 'name': 'd3-b3lyp-d3bj', 'subject': 'atom', 'lbl': 'B3LYP-D3(BJ)'}),
])  # yapf: disable
def test_dftd3__run_dftd3__2body(inp, subjects, request):
    subject = subjects()[inp['parent']][inp['subject']]
    expected = ref[inp['parent']][inp['lbl']][inp['subject']]
    gexpected = gref[inp['parent']][inp['lbl']][inp['subject']]

    if 'qcmol' in request.node.name:
        mol = subject
    else:
        mol = subject.to_schema(dtype=2)

    resinp = {
        'schema_name': 'qcschema_input',
        'schema_version': 1,
        'molecule': mol,
        'driver': 'gradient',
        'model': {
            'method': inp['name']
        },
        'keywords': {},
    }
    jrec = qcng.compute(resinp, 'dftd3', raise_error=True)
    jrec = jrec.dict()

    assert len(jrec['extras']['qcvars']) == 8

    assert compare_values(expected, jrec['extras']['qcvars']['CURRENT ENERGY'], atol=1.e-7)
    assert compare_values(expected, jrec['extras']['qcvars']['DISPERSION CORRECTION ENERGY'], atol=1.e-7)
    assert compare_values(expected, jrec['extras']['qcvars']['2-BODY DISPERSION CORRECTION ENERGY'], atol=1.e-7)
    assert compare_values(expected, jrec['extras']['qcvars'][inp['lbl'] + ' DISPERSION CORRECTION ENERGY'], atol=1.e-7)

    assert compare_values(gexpected, jrec['extras']['qcvars']['CURRENT GRADIENT'], atol=1.e-7)
    assert compare_values(gexpected, jrec['extras']['qcvars']['DISPERSION CORRECTION GRADIENT'], atol=1.e-7)
    assert compare_values(gexpected, jrec['extras']['qcvars']['2-BODY DISPERSION CORRECTION GRADIENT'], atol=1.e-7)
    assert compare_values(
        gexpected, jrec['extras']['qcvars'][inp['lbl'] + ' DISPERSION CORRECTION GRADIENT'], atol=1.e-7)


@using("dftd3_321")
@pytest.mark.parametrize(
    "subjects",
    [
        pytest.param(eneyne_ne_psi4mols, marks=using("psi4")),
        pytest.param(eneyne_ne_qcdbmols,
                     marks=using("psi4")),  # needs qcdb.Molecule, presently more common in psi4 than in qcdb
        pytest.param(eneyne_ne_qcschemamols),
    ],
    ids=['qmol', 'pmol', 'qcmol'])
@pytest.mark.parametrize("inp", [
    ({'parent': 'eneyne', 'name': 'd3-atmgr', 'subject': 'dimer', 'lbl': 'ATM'}),
    ({'parent': 'eneyne', 'name': 'd3-b3lyp-atmgr', 'subject': 'mA', 'lbl': 'ATM'}),
    ({'parent': 'eneyne', 'name': 'd3-pbe-atm(gr)', 'subject': 'mB', 'lbl': 'ATM'}),
    ({'parent': 'eneyne', 'name': 'd3-ATMgr', 'subject': 'mAgB', 'lbl': 'ATM'}),
    ({'parent': 'eneyne', 'name': 'd3-atmgr', 'subject': 'gAmB', 'lbl': 'ATM'}),
    ({'parent': 'ne', 'name': 'd3-atmgr', 'subject': 'atom', 'lbl': 'ATM'}),
])  # yapf: disable
def test_dftd3__run_dftd3__3body(inp, subjects, request):
    subject = subjects()[inp['parent']][inp['subject']]
    expected = ref[inp['parent']][inp['lbl']][inp['subject']]
    gexpected = gref[inp['parent']][inp['lbl']][inp['subject']]

    if 'qcmol' in request.node.name:
        mol = subject
    else:
        mol = subject.to_schema(dtype=2)

    resinp = {
        'schema_name': 'qcschema_input',
        'schema_version': 1,
        'molecule': mol,
        'driver': 'gradient',
        'model': {
            'method': inp['name']
        },
        'keywords': {},
    }
    jrec = qcng.compute(resinp, 'dftd3', raise_error=True)
    jrec = jrec.dict()

    assert len(jrec['extras']['qcvars']) == 8

    assert compare_values(expected, jrec['extras']['qcvars']['CURRENT ENERGY'], atol=1.e-7)
    assert compare_values(expected, jrec['extras']['qcvars']['DISPERSION CORRECTION ENERGY'], atol=1.e-7)
    assert compare_values(expected, jrec['extras']['qcvars']['3-BODY DISPERSION CORRECTION ENERGY'], atol=1.e-7)
    assert compare_values(
        expected, jrec['extras']['qcvars']['AXILROD-TELLER-MUTO 3-BODY DISPERSION CORRECTION ENERGY'], atol=1.e-7)

    assert compare_values(gexpected, jrec['extras']['qcvars']['CURRENT GRADIENT'], atol=1.e-7)
    assert compare_values(gexpected, jrec['extras']['qcvars']['DISPERSION CORRECTION GRADIENT'], atol=1.e-7)
    assert compare_values(gexpected, jrec['extras']['qcvars']['3-BODY DISPERSION CORRECTION GRADIENT'], atol=1.e-7)
    assert compare_values(
        gexpected, jrec['extras']['qcvars']['AXILROD-TELLER-MUTO 3-BODY DISPERSION CORRECTION GRADIENT'], atol=1.e-7)
